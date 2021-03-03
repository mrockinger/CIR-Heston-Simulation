# Zhu: Pure CIR simulation -> only innovation when compared to Euler.

TVSScheme = function(v0, kappa, theta, sigma, T, n, M, Z = F){
  
  Delta = T / n
  
  v = matrix(, M, n + 1) #Note: v is now the VOLATILITY, not the variance. This is the "trick" of the scheme.
  v[,1]   = sqrt(v0)
  
  sqrt.dt = sqrt(Delta)
  Z_v = matrix(rnorm(M * n), M, n)
  
  emkt = exp(- kappa * Delta)
  
  m2 = sigma^2 / (4 * kappa) * (1 - emkt) # variance of VOLATILITY
  
  for(i in 1:n){
      
    m1 = theta + (v[, i]^2 - theta) * emkt #expected VARIANCE
      
    beta = sqrt(pmax(m1 - m2, 0))
    theta_v_star = (beta - v[, i] * emkt^(1/2)) / (1 - emkt^(1/2))
      
    v[, i + 1] = v[, i] + 1 / 2 * kappa * (theta_v_star - v[, i]) * Delta + 1 / 2 * sigma * sqrt.dt * Z_v[, i]
 
  }
    
  if(Z){
    
    return(list(v^2, Z_v))
    
  }else{
    
    return(v^2)
    
  }
  
}

# test = ZhuScheme(v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, T = 1 / 12, n = 1, M = 1000000)
# MomentsCIR(p = 1:4, kappa = 3, theta = 0.19, sigma = 0.4, v0 = 0.1, t = 1 / 12)
# Evt^1|v0     Evt^2|v0     Evt^3|v0     Evt^4|v0 
# 0.1199079295 0.0155445930 0.0021628919 0.0003210907 