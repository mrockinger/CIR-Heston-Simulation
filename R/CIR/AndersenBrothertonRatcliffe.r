ABRscheme = function(v0, kappa, theta, sigma, T, n, M, Z = F){
  
  Delta = T / n
  
  v = matrix(, M, n + 1)
  
  v[,1]   = v0
  
  sqrt.dt = sqrt(Delta)
  Z_v = matrix(rnorm(M * n), M, n)
  
  emkt = exp(- kappa * Delta)
  
    
  for(i in 1:n){
      
    Gamma2 = log(1 + 0.5 * sigma^2 * kappa^-1 * v[, i] * ( 1- emkt^2) / (emkt * v[, i] + (1 - emkt) *theta)^2)
      
    v[, i + 1] = (emkt * v[, i] + (1 - emkt) * theta) * exp(-0.5 * Gamma2 + sqrt(Gamma2) * Z_v[, i])
    
  }
    
  if(Z){
    
    return(list(v, Z_v))
    
  }else{
    
    return(v)
    
  }
  
}

# test = ABRscheme(0.1, 3, 0.19, 0.4, 1, 1, 1000000)

# > MomentsCIR(p = 1:4, kappa = 3, theta = 0.19, sigma = 0.4, v0 = 0.1, t = 1)
# Evt^1|v0    Evt^2|v0    Evt^3|v0    Evt^4|v0 
# 0.185519164 0.039244388 0.009322616 0.002457098

