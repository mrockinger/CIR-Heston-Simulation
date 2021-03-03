# Implicit Milstein for CIR process discretization proposed in Kahl/Jäckel (2006)

IM = function(v0, kappa, theta, sigma, T, n, M, Z = F){
  
  Delta = T / n
  sqrt.dt = sqrt(Delta)
  v = matrix(, M, n + 1)
  v[,1] = v0
  Z_v = matrix(rnorm(M * n), M, n)
  
  for(i in 1:n){
    
    v[,i + 1] = v[,i] + kappa * (theta - v[,i]) * Delta + sigma * sqrt(v[,i]) * Z_v[,i] * sqrt.dt + 0.25 * sigma^2 * ((Z_v[,i] * sqrt.dt)^2 - Delta)
    
  }
  
  if(Z){
    
    return(list(v, Z_v))
    
  }else{
    
    return(v)
    
  }
  
}

# test = IM(0.1, 3, 0.19, 0.4, 1 / 12, 1, 1000000) 
# MomentsCIR(p = 1:4, kappa = 3, theta = 0.19, sigma = 0.4, v0 = 0.1, t = 1 / 12)
# Evt^1|v0     Evt^2|v0     Evt^3|v0     Evt^4|v0 
# 0.1199079295 0.0155445930 0.0021628919 0.0003210907 