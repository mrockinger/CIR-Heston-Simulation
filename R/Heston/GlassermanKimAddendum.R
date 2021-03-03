# Addendum to the Gamma expansion (GE) scheme of Glasserman / Kim (2011)

library(Bessel)
library(sfsmisc) # to determine the quadrant in which phi2 is

# Define the Laplace transforms of X2 and Z

L = function(b, kappa, sigma) sqrt(2 * sigma^2 * b + kappa^2)
Phi2 = function(b, kappa, sigma, Delta, d) (L(b, kappa, sigma) * sinh(kappa * Delta /2) / (kappa * sinh(L(b, kappa, sigma) * Delta / 2)))^(d / 2)
PhiZ = function(b, kappa, sigma, Delta) (L(b, kappa, sigma) * sinh(kappa * Delta /2) / (kappa * sinh(L(b, kappa, sigma) * Delta / 2))) * (L(b, kappa, sigma) * sinh(kappa * Delta /2) / (kappa * sinh(L(b, kappa, sigma) * Delta / 2)))

# Define approximations of the cdf of X2 and Z via the POISSON algorithm of Abate and Whitt (see also Broadie / Kaya), respectively
# I included the fix from p. 285
FX2 =  function(x, kappa, sigma, Delta, d, u_epsilon_X2){ 
  
  h = 2 * pi / (x  + u_epsilon_X2)
  
  epsilon = 10^-5
  
  N = 0
  
  while(abs(Phi2(b = -1i * h * N, kappa, sigma, Delta, d)) / N > pi * epsilon / 2){
    
    N = N + 1
    
  }
  
  
  k = 1:N
  
  phi2 = Phi2(b = -1i * h * k, kappa, sigma, Delta, d)
  
  # Find entries, where we move from second to third quadrant
  # Then replace argument of Phi2 there to get new values.
  for(j in 1:(N-1)){
    
    if(quadrant(phi2[j + 1]) == 3 && quadrant(phi2[j]) == 2){
      
      phi2[j + 1] = Phi2(b = -1i * (h * (j + 1) + 2 * pi), kappa, sigma, Delta, d)
      
    }
    
  }
  
  (h * x   + 2 * sum( sin(h * k * x) / k * Re(phi2) ) ) / pi
  
  
} 
FX2 = Vectorize(FX2)

FZ =  function(x, kappa, sigma, Delta, u_epsilon_Z){ 
  
  h = 2 * pi / (x  + u_epsilon_Z)
  
  epsilon = 10^-5
  
  N = 0
  
  while(abs(PhiZ(b = -1i * h * N, kappa, sigma, Delta)) / N > pi * epsilon / 2){
    
    N = N + 1
    
  }
  
  
  k = 1:N
  
  h * x / pi  + 2 / pi * sum( sin(h * k * x) / k * Re(PhiZ(b = -1i * h * k, kappa, sigma, Delta)) )
  
} 
FZ = Vectorize(FZ)

# Define the pdf of eta for n \in N, a  r.v. that follows the Bessel distribution with paramters nu and z
# pdf.eta.recursive = function(n, nu, z){
# 
#   if(n == 0){
# 
#     ((z / 2)^nu) / (BesselI(z = z, nu = nu) * gamma(nu + 1))
# 
#   }else if(n > 0){
# 
#     z^2 / (4 * n * (n + nu)) *  pdf.eta.recursive(n - 1, nu, z)
# 
#   }
# 
# }
# pdf.eta.recursive = Vectorize(pdf.eta.recursive)

# Use alternative pdf definition by Yuan / Kalbfleisch (2000), which is evaluated much faster
# Note: This function will yield NaNs for extreme values of d and z, say z = 200 and d = 0.3, since the first term is 0 and the second term is Inf

# pdf.eta = function(n, nu, z) (BesselI(z, nu) * factorial(n) * gamma(n + nu + 1))^-1 * (z / 2)^(2 * n + nu)

# plot(pdf.eta(n = 0, z = seq(0, 5, len = 100), nu = - 0.96) ~ seq(0, 5, len = 100), type = "l", xlab = "z", ylab = "Pr(eta = 0)")


# Solution: Combine both approaches
# pdf.eta.proper = function(n, nu , z){
# 
#  if(n < 0) stop("n must be an integer >= 0")
#  if(n %% 1 != 0) stop("n must be an integer")
# 
#  if(n < 50){ # Use definition of the pdf from Yuan / Kalbfleisch which is way faster to calculate but unfortunately gets instable for big values of n
# 
#    (BesselI(z, nu) * factorial(n) * gamma(n + nu + 1))^-1 * (z / 2)^(2 * n + nu)
# 
#  }else{
# 
#    pdf.eta.recursive(n, nu, z)
# 
#  }
# 
# }
# pdf.eta.proper = Vectorize(pdf.eta.proper)

# Even better: Approach suggested by Michael. It is extremely fast and stable!
# I also activate the exponential scaling in the modified Bessel function all the time to avoid overflow
# I cancel out the scaling factor by multiplying by exp(-abs(Re(z)))^-1

pdf.eta.proper.2 = function(n, nu, z){
  
  l = -lfactorial(n) + log(z) * (2 * n + nu) - log(exp(-abs(Re(z)))^-1 * BesselI(z, nu, expon.scaled = T)) - lgamma(n + nu + 1) + log(2) * (-2 * n - nu )
  exp(l);
  
}

# test = pdf.eta.proper(0:300, nu = -0.81, z = 170.809)
# test1 = pdf.eta.recursive(0:300, nu = -0.81, z = 170.809)
# test2 = pdf.eta.proper.2(0:300, nu = -0.81, z = 170.809)

# all.equal(test, test1)
# all.equal(test, test2)