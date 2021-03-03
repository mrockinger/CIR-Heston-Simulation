library(Bessel) # Is needed because the base R besselI function does not evaluate at complex arguments
library(actuar) # For drawing starting values when numerically inverting the approximated cdf of the IV
library(nleqslv) # For numerical inversion of the approximated cdf of the IV

# Here, I do basically the same as in the B-K-Scheme, but with a different characeristic function. For details of the original algorithm implementation see. BroadieKayaAddendum.r
# Define z = w * (Vs + Vt) / 2 + (1 - w) * sqrt(Vs * Vt), where w \in [0,1]


Psi.IV.Smith = function(u,  kappa, theta, sigma, Delta, z){
  
  d = 4 * kappa * theta / sigma^2
  gamma = sqrt(complex(real = kappa^2, imaginary = - 2 * sigma^2 * u))
  
  A = (4 * gamma * exp(- 0.5 * gamma * Delta)) / (sigma^2 * (1 - exp(- gamma * Delta))) * z
  B = (4 * kappa * exp(- 0.5 * kappa * Delta)) / (sigma^2 * (1 - exp(- kappa * Delta))) * z
  C = (2 * z / sigma^2) * ( (kappa * (1 + exp(- kappa * Delta))) / (1 - exp(- kappa * Delta)) - (gamma * (1 + exp(- gamma * Delta))) /  (1 - exp(- gamma * Delta)) )
  
  if(B >= 50){
    
    Psi = (gamma * exp(- 0.5 * (gamma - kappa) * Delta) * (1 - exp(- kappa * Delta))) / (kappa * (1 - exp( - gamma * Delta))) * BesselI(A, nu = d/2 - 1, expon.scaled = T) / BesselI(B, nu = d/2 - 1, expon.scaled = T) * exp(C)
    
  }else{
    
    Psi = (gamma * exp(- 0.5 * (gamma - kappa) * Delta) * (1 - exp(- kappa * Delta))) / (kappa * (1 - exp( - gamma * Delta))) * BesselI(A, nu = d/2 - 1) / BesselI(B, nu = d/2 - 1) * exp(C)
    
  }  
  
}
Psi.IV.Smith = Vectorize(Psi.IV.Smith, vectorize.args = "u")

m1IV.Smith = function(kappa, theta, sigma, Delta, z, scaled = F) {
  
  
  ekt = exp(kappa * Delta)
  sekt = sqrt(ekt)
  
  Feller = 2 * kappa * theta / sigma^2
  Arg = (4*sekt*kappa*z)/((-1 + ekt)*sigma^2)
  
  BM1 = BesselI(Arg, Feller - 1, expon.scaled = scaled)
  
  
  (2*sekt*kappa*(2 + Delta*kappa + ekt*(-2 + Delta*kappa))*z*BesselI(Arg, Feller) + 
      2*sekt*kappa*(2 + Delta*kappa + ekt*(-2 + Delta*kappa))*z*
      BesselI(Arg, -2 + Feller) + 
      ((-1 + ekt)*(2 + Delta*kappa + ekt*(-2 + Delta*kappa))*sigma^2 - 4*kappa*(1 - ekt^2 + 2*Delta*ekt*kappa)*z)*
      BM1)/(2.* (-1 + ekt)^2*kappa^2*BM1)
  
  
  
  
}

# m1IV.Smith(3, 0.19, 0.4, 1 / 12, 0.2)

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




m2IV.Smith= function(kappa, theta, sigma, Delta, z, scaled = F){
  
  z2 = z^2
  
  D2 = Delta^2
  
  s2 = sigma^2
  s4 = sigma^4
  
  k2 = kappa^2
  k3 = kappa^3
  k4 = kappa^4
  
  ekt = exp(kappa * Delta)
  sekt = sqrt(ekt)
  ekt2 = ekt^2
  ekt3 = ekt^3
  ekt4 = ekt^4
  thekt = ekt^(3/2)
  fhekt = ekt^(5/2)
  shekt = ekt^(7/2)
  
  Feller = 2 * kappa * theta / s2
  Arg = (4*sekt*kappa*z)/((-1 + ekt)*s2)
  
  B1 = BesselI(Arg, Feller + 1, expon.scaled = scaled)
  BM1 = BesselI(Arg, Feller - 1, expon.scaled = scaled)
  BM2 = BesselI(Arg, Feller - 2, expon.scaled = scaled)
  
  
  
  (2*sekt*kappa*z*((-1 + ekt)*(4 + 10*Delta*kappa + 3*D2*k2 + ekt2*(4 - 10*Delta*kappa + 3*D2*k2) + 
                                 2*ekt*(-4 + 5*D2*k2))*s2 - 
                     8*kappa*(2 + Delta*kappa + ekt3*(2 - Delta*kappa) + ekt2*(-2 - 5*Delta*kappa + 2*D2*k2) + 
                                ekt*(-2 + 5*Delta*kappa + 2*D2*k2))*z)*BesselI(Arg, Feller) + 
      4*ekt*k2* (2 + Delta*kappa + ekt*(-2 + Delta*kappa))^2*z2*
      BesselI(Arg, -3 + Feller) - 
      8*sekt*kappa*s2*z*BM2 + 
      24*thekt*kappa*s2*z*BM2 - 
      24*fhekt*kappa*s2*z*BM2 + 
      8*shekt*kappa*s2*z*BM2 - 
      20*Delta*sekt*k2*s2*z*BM2 + 
      20*Delta*thekt*k2*s2*z*BM2 + 
      20*Delta*fhekt*k2*s2*z*BM2 - 
      20*Delta*shekt*k2*s2*z*BM2 - 
      6*D2*sekt*k3*s2*z*BM2 - 
      14*D2*thekt*k3*s2*z*BM2 + 
      14*D2*fhekt*k3*s2*z*BM2 + 
      6*D2*shekt*k3*s2*z*BM2 - 
      32*sekt*k2*z2*BM2 + 
      32*thekt*k2*z2*BM2 + 
      32*fhekt*k2*z2*BM2 - 
      32*shekt*k2*z2*BM2 - 
      16*Delta*sekt*k3*z2*BM2 - 
      80*Delta*thekt*k3*z2*BM2 + 
      80*Delta*fhekt*k3*z2*BM2 + 
      16*Delta*shekt*k3*z2*BM2 - 
      32*D2*thekt*k4*z2*BM2 - 
      32*D2*fhekt*k4*z2*BM2 - 
      4*s4*BM1 + 
      16*ekt*s4*BM1 - 
      24*ekt2*s4*BM1 + 
      16*ekt3*s4*BM1 - 
      4*ekt4*s4*BM1 + 
      2*Delta*kappa*s4*BM1 - 
      4*Delta*ekt*kappa*s4*BM1 + 
      4*Delta*ekt3*kappa*s4*BM1 - 
      2*Delta*ekt4*kappa*s4*BM1 + 
      D2*k2*s4*BM1 + 
      4*D2*ekt*k2*s4*BM1 - 
      10*D2*ekt2*k2*s4*BM1 + 
      4*D2*ekt3*k2*s4*BM1 + 
      D2*ekt4*k2*s4*BM1 + 
      8*kappa*s2*z*BM1 - 
      16*ekt*kappa*s2*z*BM1 + 
      16*ekt3*kappa*s2*z*BM1 - 
      8*ekt4*kappa*s2*z*BM1 + 
      8*Delta*k2*s2*z*BM1 + 
      48*Delta*ekt*k2*s2*z*BM1 - 
      112*Delta*ekt2*k2*s2*z*BM1 + 
      48*Delta*ekt3*k2*s2*z*BM1 + 
      8*Delta*ekt4*k2*s2*z*BM1 + 
      32*D2*ekt*k3*s2*z*BM1 - 
      32*D2*ekt3*k3*s2*z*BM1 + 
      16*k2*z2*BM1 + 
      32*ekt*k2*z2*BM1 - 
      96*ekt2*k2*z2*BM1 + 
      32*ekt3*k2*z2*BM1 + 
      16*ekt4*k2*z2*BM1 + 
      96*Delta*ekt*k3*z2*BM1 - 
      96*Delta*ekt3*k3*z2*BM1 + 
      8*D2*ekt*k4*z2*BM1 + 
      80*D2*ekt2*k4*z2*BM1 + 
      8*D2*ekt3*k4*z2*BM1 + 
      16*ekt*k2*z2*B1 - 
      32*ekt2*k2*z2*B1 + 
      16*ekt3*k2*z2*B1 + 
      16*Delta*ekt*k3*z2*B1 - 
      16*Delta*ekt3*k3*z2*B1 + 
      4*D2*ekt*k4*z2*B1 + 
      8*D2*ekt2*k4*z2*B1 + 
      4*D2*ekt3*k4*z2*B1)/
    (4.* (-1 + ekt)^4*k4*BM1)
  
}
# m2IV.Smith(3, 0.19, 0.4, 1 / 12, 0.2, scaled = F)


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


# Cdf.IV.approx.Smith = function(x, z, Delta, kappa, theta, sigma, u_epsilon){
#   
#   h_z = 2 * pi / (x + u_epsilon)
#   
#   epsilon = 10^-5
#   
#   N_z = 0
#   
#   while(abs(Psi.IV.Smith(h_z * N_z, kappa, theta, sigma, Delta, z)) / N_z > pi * epsilon / 2){
#     
#     N_z = N_z + 1
#     
#   }
#   
#   A = h_z* x / pi
#   vec = (1:N_z) * h_z
#   B = (2 / pi) * sum( (sin(vec * x) / (1:N_z)) * Re(Psi.IV.Smith(vec, kappa, theta, sigma, Delta, z)))
#   
#   Cdf.appprox = A + B
#   
# }
# 
# 
# Cdf.IV.approx.diff.Smith = function(x, z, Delta, kappa, theta, sigma, u_epsilon, U){
#   
#   Cdf.IV.approx.Smith(x, z, Delta, kappa, theta, sigma, u_epsilon) - U
#   
# }
# 
# GenerateIV.Smith = function(z, Delta, kappa, theta, sigma, U){
#   
#   m1 = m1IV.Smith(kappa, theta, sigma, Delta, z)
#   m2 = m2IV.Smith(kappa, theta, sigma, Delta, z)
#   
#   u_epsilon = m1 + 5 * sqrt(m2 - m1^2)
#   
#   StartVal = rinvgauss(1, mean = m1, dispersion = sqrt(m2 - m1^2))
#   
#   if(StartVal < 0) StartVal = 0.01 * m1
#   
#   IV = nleqslv(x = StartVal,  fn = Cdf.IV.approx.diff.Smith, z = z, Delta = Delta, kappa = kappa, theta = theta, sigma = sigma, u_epsilon = u_epsilon, U = U, method = "Newton")$x
#   
#   if(IV < 0) IV = uniroot(f = Cdf.IV.approx.diff.Smith, interval = c(0, 2), z = z, Vt = Vt, Delta = Delta, kappa = kappa, theta = theta, sigma = sigma, u_epsilon = u_epsilon, U = U)$root
#   
#   return(IV)
# }
# GenerateIV.Smith = Vectorize(GenerateIV.Smith)

# test = GenerateIV.Smith(z = 0.2, Delta = 1/12, kappa = 3, theta = 0.19, sigma = 0.3, U = seq(0, 1, len = 100))

