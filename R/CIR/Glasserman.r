G = function(v0, kappa, theta, sigma, T, n, M, Z = F){
  
  v = matrix(, M, n + 1)
  v[, 1] = v0
  
  Delta = T / n
  sqrt.dt = sqrt(Delta)
  
  Z_v = matrix(rnorm(M * n), M, n)
  
  for(i in 1:n){
    
    v[,i + 1] = 
      
      abs(
        v[,i] + kappa * (theta - v[,i]) * Delta + sigma * sqrt(v[,i]) * sqrt.dt * Z_v[,i] + sigma^2 / 4 * Delta * (Z_v[,i]^2 - 1) - 0.5 * kappa^2 * (theta - v[,i]) * Delta^2 + ((kappa * theta / 4 - sigma^2 / 16) / sqrt(v[,i]) - 3 * kappa / 2 * sqrt(v[,i])) * Z_v[,i] * sigma * Delta^(3 / 2)
      )
    
  }  
  
  if(Z){
    
    return(list(v, Z_v))
    
  }else{
    
    return(v)
    
  }
  
}

# test = G(0.1, 3, 0.19, 0.4, 1 / 12, 1, 1000000)
# MomentsCIR(p = 1:4, kappa = 3, theta = 0.19, sigma = 0.4, v0 = 0.1, t = 1 / 12)
# Evt^1|v0     Evt^2|v0     Evt^3|v0     Evt^4|v0 
# 0.1199079295 0.0155445930 0.0021628919 0.0003210907 