# Broadie/Kaya; Alfonsi; Glasserman; Johnson et al's approach of simulating the non-central chi-square distribution directly
# I will use the implemented algorithm to draw a chi-squared RV in case of d > 1 (this will use the Gamma distribution algorithm)
# I will also use the implemented algorithm to draw a non-central chi-squared RV if d<=1
# Thie latter one is the same algorithm as described in Glasserman. So this is really just a "lazy shortcut"


Exact = function(v0, kappa, theta, sigma, T, n, M){

  Delta = T / n
  emkt = exp(- kappa * Delta)
  d = 4 * kappa * theta / sigma^2
  h = 4 * kappa / (sigma^2 * (1 - emkt))
  v = matrix(, M, n + 1)
  v[,1] = v0
  
  
  if(d > 1){
    
    Z_v = matrix(rnorm(M * n), M, n)
    X = matrix(rchisq(M * n, df = d - 1), M, n)
    
    for(i in 1:n){
      
      lambda = h * v[,i] * emkt
      # Here, it is important to omit the ncp parameter. rchisq will use the algorithm from rgamma. 
      v[,i + 1] = h^-1 * ((Z_v[,i] + sqrt(lambda))^2 + X[,i])
        
    }
    
    
  }else{
    
    for(i in 1:n){
      
      lambda = h * v[,i] * emkt
      N = sapply(lambda / 2, rpois, n = 1)
      v[,i + 1] = h^-1 * sapply(d + 2 * N, rchisq, n = 1) 
      
      #Used algorithm corresponds to drawing a Poisson mixture of central chi-squared RVs
      # v[,i + 1] = h^-1 * sapply(h * v[,i] * emkt, rchisq, n = 1, df = d)
      
    }
  
  } 
  
  return(v)
    
  
}

# test = Exact(0.1, 0.5, 0.19, 1, 1, 1, 1000000)
# MomentsCIR(p = 1:4, kappa = 0.5, theta = 0.19, sigma = 1, v0 = 0.1, t = 1)
# Evt^1|v0   Evt^2|v0   Evt^3|v0   Evt^4|v0 
# 0.13541224 0.09548216 0.11331171 0.18573172 