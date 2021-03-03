# In this file, I build a function to generate exact draws from the integrated variance (IV) given Vs and Vt, see Broadie/Kaya p. 220 - 221

# The following package is required since the base R besselI function does not evaluate at complex arguments. 
# With this package, the modified Bessel function of the first kind is calculated with exactly the same algorithm as used by Broadie/Kaya, see p . 221
library(Bessel)

# Required for drawing starting values when numerically inverting the approximated cdf of the IV
library(actuar) 

# Required for numerical inversion of the IV's approximated c.d.f. via Netwon's method
library(nleqslv) 

# Functions to evaluate C-type code from Mathematica (not needed anymore)
# Power = function(a,b) a^b
# E = exp(1)
# Sqrt = function(a) sqrt(a)


# Note: The argument of the modified Bessel functions (A and B in the following) are too large to handle if Delta gets too small

# The potential solution to this problem is outlined below:

# kappa = 3; theta = 0.19; sigma = 0.4; Delta = seq(1 / 130, 1, len = 100) ; Vs = 0.2; Vt = 0.202; u = 100 # Here, Delta is not "small enough" to become a problem
# d = 4 * kappa * theta / sigma^2
# gamma = sqrt(complex(real = kappa^2, imaginary = - 2 * sigma^2 * u))
#  
# A = (4 * gamma * exp(- 0.5 * gamma * Delta)) / (sigma^2 * (1 - exp(- gamma * Delta))) * sqrt(Vs * Vt)
# B = (4 * kappa * exp(- 0.5 * kappa * Delta)) / (sigma^2 * (1 - exp(- kappa * Delta))) * sqrt(Vs * Vt)
#  
# test.unscaled = BesselI(A, nu = d/2 - 1) / BesselI(B, nu = d/2 - 1)
# test.scaled = BesselI(A, nu = d/2 - 1, expon.scaled = T) / BesselI(B, nu = d/2 - 1, expon.scaled = T)

# plot(Re(test.unscaled) - Re(test.scaled) ~ Delta, type = "l")
# plot(Im(test.unscaled) - Im(test.scaled) ~ Delta, type = "l") 

# Ergo:  I(A) / I(B) = I_scaled(A) / I_scaled(B), if Delta <= 0.2 = 1 / 5 and thus (more important) B \approx Re(A) >= 23.75
# In the function, I will switch btw. the two calculation approaches if B > 50 since the overflow will only occure for much larger values of B -> approximation is more save


# Finally, define the characteristic function of the IV

Phi.IV = function(u, kappa, theta, sigma, Delta, Vs, Vt){
  
  d = 4 * kappa * theta / sigma^2
  gamma = sqrt(complex(real = kappa^2, imaginary = - 2 * sigma^2 * u))
  
  # A is the complex argument of the modified Bessel function -> complex logarithm will be calculated. Therefore, we will keep track of Arg(A) to not have discontinuities, see p. 221 on the right side
  # If -pi < Arg(A) < pi does NOT hold, evaluate exp(complex(real = 0, imaginary = m * (d/2 - 1) * pi)) * BesselI(A, nu = d/2 - 1, expon.scaled = T) instead of BesselI(A, nu = d/2 - 1, expon.scaled = T)
  A = (4 * gamma * exp(- 0.5 * gamma * Delta)) / (sigma^2 * (1 - exp(- gamma * Delta))) * sqrt(Vs * Vt)
  B = (4 * kappa * exp(- 0.5 * kappa * Delta)) / (sigma^2 * (1 - exp(- kappa * Delta))) * sqrt(Vs * Vt)
  C = ((Vs + Vt) / sigma^2) * ( (kappa * (1 + exp(- kappa * Delta))) / (1 - exp(- kappa * Delta)) - (gamma * (1 + exp(- gamma * Delta))) /  (1 - exp(- gamma * Delta)) )
  
  
  if(B >= 50){
    
    # if(Arg(A) < - pi || Arg(A) > pi){
      
      
      # Psi = (gamma * exp(- 0.5 * (gamma - kappa) * Delta) * (1 - exp(- kappa * Delta))) / (kappa * (1 - exp( - gamma * Delta))) *  exp (complex(real = 0, imaginary = (d/2 - 1) * pi)) * BesselI(A, nu = d/2 - 1, expon.scaled = T) / BesselI(B, nu = d/2 - 1, expon.scaled = T) * exp(C)
      
      
    # }else{
      
        Psi = (gamma * exp(- 0.5 * (gamma - kappa) * Delta) * (1 - exp(- kappa * Delta))) / (kappa * (1 - exp( - gamma * Delta))) * BesselI(A, nu = d/2 - 1, expon.scaled = T) / BesselI(B, nu = d/2 - 1, expon.scaled = T) * exp(C)
      
    # }  
      
    }else{
    
    # if(Arg(A) < - pi || Arg(A) > pi)
  
    Psi = (gamma * exp(- 0.5 * (gamma - kappa) * Delta) * (1 - exp(- kappa * Delta))) / (kappa * (1 - exp( - gamma * Delta))) * BesselI(A, nu = d/2 - 1) / BesselI(B, nu = d/2 - 1) * exp(C)
  
  }  
  
}
Phi.IV = Vectorize(Phi.IV)

# old stuff, where I compared scaling everytime/no scaling/scaling when it is necessary

# test  = Psi.IV(u = 100, kappa = 3, theta = 0.19, sigma = 0.4, Delta = 1 / 12, Vs = 0.24, Vt = 0.25)
# test1 = Psi.IV(u = 100, kappa = 3, theta = 0.19, sigma = 0.4, Delta = 1 / 12, Vs = 0.24, Vt = 0.25, scaled = F)
# test2 = Psi.IV.V2(u = 100, kappa = 3, theta = 0.19, sigma = 0.4, Delta = 1 / 12, Vs = 0.24, Vt = 0.25)
# 
# plot(Re(test) ~ Delta, type = "l", col = 2)
# lines(Re(test1) ~ Delta, col = 1)
# lines(Re(test2) ~ Delta, col = 3)
 
# Define the first two moments of the IV to be able to calculate an optimal grid size in the next steps

# Old, complete copy pasted from Mathematica code

# m1IV = function(kappa, theta, sigma, Delta, Vs, Vt){
#   
#   (0.5*(1 - Power(E,-(Delta*kappa)))*((2.*(Power(E,2.5*Delta*kappa)*(-2. + 1.*Delta*kappa) + Power(E,1.5*Delta*kappa)*(2. + 1.*Delta*kappa))*
#                                          Sqrt(Vs*Vt)*(BesselI( (4*Power(E,0.5*Delta*kappa)*kappa*Sqrt(Vs*Vt))/((-1 + Power(E,Delta*kappa))*Power(sigma,2)) , (2.*kappa*theta)/Power(sigma,2)) + 
#                                                         BesselI((4*Power(E,0.5*Delta*kappa)*kappa*Sqrt(Vs*Vt))/((-1 + Power(E,Delta*kappa))*Power(sigma,2)), -2 + (2.*kappa*theta)/Power(sigma,2)))
#   )/(Power(E,2.*Delta*kappa)*(-1. + Power(E,Delta*kappa))) + 
#     (2.*Delta*Power(sigma,2)*BesselI((4*Power(E,0.5*Delta*kappa)*kappa*Sqrt(Vs*Vt))/((-1 + Power(E,Delta*kappa))*Power(sigma,2)), -1 + (2.*kappa*theta)/Power(sigma,2)))/Power(E,Delta*kappa) - 
#     (2.*(-1. + Power(E,Delta*kappa))*Power(sigma,2)*BesselI((4*Power(E,0.5*Delta*kappa)*kappa*Sqrt(Vs*Vt))/((-1 + Power(E,Delta*kappa))*Power(sigma,2)), -1 + (2.*kappa*theta)/Power(sigma,2)))/(Power(E,Delta*kappa)*kappa) + 
#     (1.*(1.*Delta*Power(1. - 1.*Power(E,Delta*kappa),2)*Power(sigma,2) + (-2. + 2.*Power(E,2*Delta*kappa))*Vs + 
#            Delta*Power(E,Delta*kappa)*kappa*(-4.*Vs - 4.*Vt) + (-2. + 2.*Power(E,2*Delta*kappa))*Vt)*
#        BesselI((4*Power(E,0.5*Delta*kappa)*kappa*Sqrt(Vs*Vt))/((-1 + Power(E,Delta*kappa))*Power(sigma,2)), -1 + (2.*kappa*theta)/Power(sigma,2)))/
#     (Power(E,Delta*kappa)*(-1. + Power(E,Delta*kappa)))))/
#     (Power(-1 + Power(E,-(Delta*kappa)),2)*kappa*BesselI((4*Power(E,0.5*Delta*kappa)*kappa*Sqrt(Vs*Vt))/((-1 + Power(E,Delta*kappa))*Power(sigma,2)), -1 + (2.*kappa*theta)/Power(sigma,2)))
#   
# }

# Bit optimized version w.o. C-functions


m1IV = function(kappa, theta, sigma, Delta, Vs, Vt, scaled = T){
  
  ekt = exp(Delta * kappa)
  Arg = (4 * sqrt(ekt) * kappa * sqrt(Vs*Vt)) / ((-1 + ekt) * sigma^2)
  Feller = (2 * kappa * theta)/ sigma^2
  BM1 = BesselI(Arg, Feller - 1, expon.scaled = scaled) 
  
  
  (0.5 * (1 - ekt^-1) * ((2 * (ekt^2.5 * (- 2 + Delta * kappa) + ekt^1.5 * (2 + Delta * kappa)) *
  sqrt(Vs * Vt) * (BesselI(Arg , Feller, expon.scaled = scaled) + BesselI(Arg, Feller - 2, expon.scaled = scaled))) /(ekt^2 * (- 1 + ekt)) + 
  (2 * Delta * sigma^2 * BM1) / ekt - (2 * (-1 + ekt) * sigma^2 * BM1)/( ekt * kappa) + 
  ((Delta * (1 - ekt)^2 * sigma^2 + (-2 + 2 * ekt^2) * Vs + 
   Delta * ekt * kappa * (-4 * Vs - 4 * Vt) + (-2 + 2 * ekt^2) * Vt) * BM1)/
  (ekt * (-1 + ekt)))) / ((-1 + ekt^-1)^2 * kappa * BM1)
  
}



# m1IV(3, 0.19, 0.4, 1/12, 0.2, 0.202)

 
# Here, I tested whether the exponential scaling of the Bessel function (to avoid overflow) has negative impact on the results
# Turns out, it does not!
# Ergo: use scaled version all the time
# See also the moment computation approach of Glasserman/Kim (2011) also used in Tse/Wan (2013). From their formulas it is clear, why scaling does not change the result

# Delta = seq(1 / (252 * 6.5), 3, len = 10000)
# 
# test = m1IV(kappa, theta, sigma, Delta, Vs, Vt, scaled = F)
# test2 = m1IV(kappa, theta, sigma, Delta, Vs, Vt)
# 
# plot(test ~ Delta, type = "l", ylab = "m1IV")
# lines(test2 ~ Delta, col = "red3")
# 
# head(test)
# head(test2)
# 
# tail(test)
# tail(test2)

# Again old, copy pasted from Mathematica code

# m2IV = function(kappa, theta, sigma, Delta, Vs, Vt, scaled = T){

# ((4*Power(-1 + Power(E,Delta*kappa),2)*(1 + Power(E,Delta*kappa))*(Vs + Vt)*
#     ((-1 + Power(E,Delta*kappa))*Power(sigma,2) + (1 + Power(E,Delta*kappa))*kappa*(Vs + Vt)) +
#     Delta^2*kappa*(Power(-1 + Power(E,Delta*kappa),2)*(1 + 6*Power(E,Delta*kappa) + Power(E,2*Delta*kappa))*Power(sigma,4) -
#                             16*Power(E,Delta*kappa)*(-1 + Power(E,2*Delta*kappa))*kappa*Power(sigma,2)*(Vs + Vt) +
#                             16*Power(E,2*Delta*kappa)*Power(kappa,2)*Power(Vs + Vt,2)) +
#     2*Delta*(-1 + Power(E,Delta*kappa))*(Power(-1 + Power(E,Delta*kappa),2)*(1 + Power(E,Delta*kappa))*Power(sigma,4) +
#                                            2*(-1 - 3*Power(E,Delta*kappa) + 3*Power(E,2*Delta*kappa) + Power(E,3*Delta*kappa))*kappa*Power(sigma,2)*(Vs + Vt) -
#                                            8*Power(E,Delta*kappa)*(1 + Power(E,Delta*kappa))*Power(kappa,2)*Power(Vs + Vt,2)))/Power(-1 + Power(E,Delta*kappa),4) +
#    (4*(2*(-1 + Power(E,2*Delta*kappa))*(Vs + Vt) + Delta*
#          ((-1 + Power(E,2*Delta*kappa))*Power(sigma,2) - 4*Power(E,Delta*kappa)*kappa*(Vs + Vt)))*
# (Power(E,(Delta*kappa)/2.)*kappa*(2 + Delta*kappa + Power(E,Delta*kappa)*(-2 + Delta*kappa))*Sqrt(Vs*Vt)*
# BesselI((4*Power(E,(Delta*kappa)/2.)*kappa*Sqrt(Vs*Vt))/((-1 + Power(E,Delta*kappa))*Power(sigma,2)), (2*kappa*theta)/Power(sigma,2), expon.scaled = scaled)  #\
# +  Power(E,(Delta*kappa)/2.)*kappa*(2 + Delta*kappa + Power(E,Delta*kappa)*(-2 + Delta*kappa))*Sqrt(Vs*Vt)*
# BesselI((4*Power(E,(Delta*kappa)/2.)*kappa*Sqrt(Vs*Vt))/
# ((-1 + Power(E,Delta*kappa))*Power(sigma,2)), -2 + (2*kappa*theta)/Power(sigma,2), expon.scaled = scaled) -
# Power(-1 + Power(E,Delta*kappa),2)*Power(sigma,2)*
# BesselI((4*Power(E,(Delta*kappa)/2.)*kappa*Sqrt(Vs*Vt))/
# ((-1 + Power(E,Delta*kappa))*Power(sigma,2)), -1 + (2*kappa*theta)/Power(sigma,2), expon.scaled = scaled)))/
# (Power(-1 + Power(E,Delta*kappa),4)*BesselI(
# (4*Power(E,(Delta*kappa)/2.)*kappa*Sqrt(Vs*Vt))/((-1 + Power(E,Delta*kappa))*Power(sigma,2)), -1 + (2*kappa*theta)/Power(sigma,2), expon.scaled = scaled)) -
# (4*(1 - Power(E,-(Delta*kappa)))*((2*Power(E,(Delta*kappa)/2.)*kappa*(2 + Delta*kappa + Power(E,Delta*kappa)*(-2 + Delta*kappa))*
# Power(sigma,2)*Sqrt(Vs*Vt)*(BesselI(
# (4*Power(E,(Delta*kappa)/2.)*kappa*Sqrt(Vs*Vt))/((-1 + Power(E,Delta*kappa))*Power(sigma,2)), (2*kappa*theta)/Power(sigma,2), expon.scaled = scaled) +
# BesselI(
# (4*Power(E,(Delta*kappa)/2.)*kappa*Sqrt(Vs*Vt))/((-1 + Power(E,Delta*kappa))*Power(sigma,2)), -2 + (2*kappa*theta)/Power(sigma,2), expon.scaled = scaled)))/
# Power(-1 + Power(E,Delta*kappa),2) + Power(sigma,4)*
# BesselI((4*Power(E,(Delta*kappa)/2.)*kappa*Sqrt(Vs*Vt))/
# ((-1 + Power(E,Delta*kappa))*Power(sigma,2)), -1 + (2*kappa*theta)/Power(sigma,2), expon.scaled = scaled) +
# (Power(E,(Delta*kappa)/2.)*kappa*((1 - Power(E,Delta*kappa))*
# (-4 + 2*Delta*kappa + Power(Delta,2)*Power(kappa,2) +
# Power(E,2*Delta*kappa)*(-4 - 2*Delta*kappa + Power(Delta,2)*Power(kappa,2)) +
# Power(E,Delta*kappa)*(8 + 6*Power(Delta,2)*Power(kappa,2)))*Power(sigma,2)*Sqrt(Vs*Vt)*
# (BesselI((4*Power(E,(Delta*kappa)/2.)*kappa*Sqrt(Vs*Vt))/
# ((-1 + Power(E,Delta*kappa))*Power(sigma,2)), (2*kappa*theta)/Power(sigma,2), expon.scaled = scaled) +
# BesselI(
# (4*Power(E,(Delta*kappa)/2.)*kappa*Sqrt(Vs*Vt))/((-1 + Power(E,Delta*kappa))*Power(sigma,2)), -2 + (2*kappa*theta)/Power(sigma,2), expon.scaled = scaled)) -
# 2*Power(E,(Delta*kappa)/2.)*kappa*Power(2 + Delta*kappa + Power(E,Delta*kappa)*(-2 + Delta*kappa),2)*Vs*Vt*
# (BesselI(
# (4*Power(E,(Delta*kappa)/2.)*kappa*Sqrt(Vs*Vt))/((-1 + Power(E,Delta*kappa))*Power(sigma,2)), -3 + (2*kappa*theta)/Power(sigma,2), expon.scaled = scaled) +
# 2*BesselI(
# (4*Power(E,(Delta*kappa)/2.)*kappa*Sqrt(Vs*Vt))/((-1 + Power(E,Delta*kappa))*Power(sigma,2)), -1 + (2*kappa*theta)/Power(sigma,2), expon.scaled = scaled) +
# BesselI(
# (4*Power(E,(Delta*kappa)/2.)*kappa*Sqrt(Vs*Vt))/((-1 + Power(E,Delta*kappa))*Power(sigma,2)), 1 + (2*kappa*theta)/Power(sigma,2), expon.scaled = scaled))))/
# (2.*Power(-1 + Power(E,Delta*kappa),4))))/
# ((kappa - kappa/Power(E,Delta*kappa))*BesselI(
# (4*Power(E,(Delta*kappa)/2.)*kappa*Sqrt(Vs*Vt))/((-1 + Power(E,Delta*kappa))*Power(sigma,2)), -1 + (2*kappa*theta)/Power(sigma,2), expon.scaled = scaled)))/(4.*Power(kappa,3))

#}

m2IV = function(kappa, theta, sigma, Delta, Vs, Vt, scaled = T){
  
  ekt = exp(kappa * Delta)
  Feller = 2 * kappa * theta / sigma^2
  Arg = (4*sqrt(ekt)*kappa*sqrt(Vs*Vt))/((-1 + ekt)*sigma^2)
  
  B = BesselI(Arg, Feller, expon.scaled = scaled) 
  BM1 = BesselI(Arg, -1 + Feller, expon.scaled = scaled)
  BM2 = BesselI(Arg, -2 + Feller, expon.scaled = scaled)
  
  ((4*(-1 + ekt)^2*(1 + ekt)*(Vs + Vt)* ((-1 + ekt)*sigma^2 + (1 + ekt)*kappa*(Vs + Vt)) +Delta^2*kappa*((-1 + ekt)^2*(1 + 6*ekt + exp(2 *kappa * Delta))*sigma^4 - 
  16 * ekt * (-1 + ekt^2) * kappa*sigma^2*(Vs + Vt) + 16*ekt^2*kappa^2*(Vs + Vt)^2) + 2*Delta*(-1 + ekt)*((-1 + ekt)^2*(1 + ekt)*sigma^4 + 
  2*(-1 - 3*ekt + 3*exp(2 * kappa * Delta) + ekt^3)*kappa*sigma^2*(Vs + Vt) - 8*ekt*(1 + ekt)*kappa^2*(Vs + Vt)^2))/ (-1 + ekt)^4 +     
  (4*(2*(-1 + ekt^2)*(Vs + Vt) + Delta*((-1 + ekt^2)*sigma^2 - 4*ekt*kappa*(Vs + Vt)))*(sqrt(ekt)*kappa*(2 + Delta*kappa + ekt*(-2 + Delta*kappa))*sqrt(Vs*Vt)*B+ 
  sqrt(ekt)*kappa*(2 + Delta*kappa + ekt*(-2 + Delta*kappa))*sqrt(Vs*Vt)*BM2 -(-1 + ekt)^2*sigma^2*BM1))/( (-1 + ekt)^4*BM1) - 
  (4*(1 - exp(- kappa * Delta))*((2*sqrt(ekt)*kappa*(2 + Delta*kappa + ekt*(-2 + Delta*kappa))*sigma^2*sqrt(Vs*Vt)*(B + BM2))/
  (-1 + ekt)^2 + sigma^4*BM1 +(sqrt(ekt)*kappa*((1 - ekt)*(-4 + 2*Delta*kappa + Delta^2*kappa^2 +     
  exp(2 * kappa * Delta)*(-4 - 2*Delta*kappa + Delta^2*kappa^2) + ekt*(8 + 6*Delta^2*kappa^2))*sigma^2*sqrt(Vs*Vt)*(B + BM2) - 
  2*sqrt(ekt)*kappa* (2 + Delta*kappa + ekt*(-2 + Delta*kappa))^2*Vs*Vt*(BesselI(Arg, -3 + Feller, expon.scaled = scaled) + 
  2*BM1 + BesselI(Arg, 1 + Feller, expon.scaled = scaled))))/(2 *(-1 + ekt)^4)))/((kappa - kappa/ekt)*BM1))/(4*kappa^3)

}

# m2IV(3, 0.19, 0.4, 1/12, 0.2, 0.2)

# Here, I tested whether the exponential scaling of the Bessel function (to avoid overflow) has negative impact on the results (same as for m1IV)
# Turns out, it does not either!
# Ergo: use scaled version all the time

# test = m2IV(kappa, theta, sigma, Delta, Vs, Vt, scaled = F)
# test2 = m2IV(kappa, theta, sigma, Delta, Vs, Vt)
# 
# plot(test ~ Delta, type = "l", ylab = "m2IV")
# lines(test2 ~ Delta, col = "red3")
# 
# head(test)
# head(test2)
# 
# tail(test)
# tail(test2)


# Define the (numerical approximation of the) cdf of the IV
# N is the number of evaluations and h the step size
# h and N have to be calculated with the help of the first two moments of the IV for EVERY argument x!
# After that, for any chosen discretization error \epsilon, N can be determined. 
# The approximation follows the POISSON algorithm in Abate and Whitt (1992)

Cdf.IV.approx = function(x, Vs, Vt, Delta, kappa, theta, sigma, u_epsilon){
  
  h = 2 * pi / (x + u_epsilon)
  
  epsilon = 10^-5
  
  N = 0
  
  while(abs(Phi.IV(h * N, kappa, theta, sigma, Delta, Vs, Vt)) / N > pi * epsilon / 2){
     
   N = N + 1
     
  }
  
  
  A = h * x / pi
  vec = (1:N) * h
  B = (2 / pi) * sum( (sin(vec * x) / (1:N)) * Re(Phi.IV(vec, kappa, theta, sigma, Delta, Vs, Vt))) - 2 * epsilon # I substract 2 * epsilon because we have an discretization and a truncation error of epsilon
  
  return(A + B)
  
}


# Define the function for the inverse transform method, i. e. the difference between the cdf at x and some uniform rv U
# The root of this function has to be found numerically and is at the same time a draw of the IV

Cdf.IV.approx.diff = function(x, Vs, Vt, Delta, kappa, theta, sigma, u_epsilon, U){
  
  Cdf.IV.approx(x, Vs, Vt, Delta, kappa, theta, sigma, u_epsilon) - U
  
}

# Note: If numerical inversion with Newton's algorithm fails to converge, I will use a bisectional method (see function uniroot).  
# Note: uniroot does not use bisectional search but Brent's derivative-free optimizer..
# Also: I did not implement this in Matlab. If IV < 0 I do the Newton - way again with another U
# This is according to B/K, p. 211 right side on the bottom

GenerateIV = function(Vs, Vt, Delta, kappa, theta, sigma, U){
  
  m1 = m1IV(kappa, theta, sigma, Delta, Vs, Vt)
  m2 = m2IV(kappa, theta, sigma, Delta, Vs, Vt) 
  
  u_epsilon = m1 + 5 * sqrt(m2 - m1^2)
    
  StartVal = rinvgauss(1, mean = m1, dispersion = sqrt(m2 - m1^2))
  
  if(StartVal < 0) StartVal = 0.01 * m1
  
  IV = nleqslv(x = StartVal,  fn = Cdf.IV.approx.diff, Vs = Vs, Vt = Vt, Delta = Delta, kappa = kappa, theta = theta, sigma = sigma, u_epsilon = u_epsilon, U = U, method = "Newton")$x
  
  if(IV < 0) IV = uniroot(f = Cdf.IV.approx.diff, interval = c(0, 2), Vs = Vs, Vt = Vt, Delta = Delta, kappa = kappa, theta = theta, sigma = sigma, u_epsilon = u_epsilon, U = U)$root
  
  return(IV)
}
GenerateIV = Vectorize(GenerateIV, vectorize.args = c("Vs", "Vt", "U"))

# test = GenerateIV(Vs = 0.1, Vt = 0.11, Delta = 1  , kappa = 0.5, theta = 0.19, sigma = 1, U = seq(0, 1, len = 100))
