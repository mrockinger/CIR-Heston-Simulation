# CIR discretization scheme of Deelstra and Delbaen (1998)
# Defintion follows Lord et al (2010), Table 1

DD = function(v0, kappa, theta, sigma, T, n, M, Z = F){
  
  Delta = T / n 
  sqrt.dt = sqrt(Delta)
  v = matrix(, M, n + 1) # This will store the auxiliary variance process
  v[,1] = v0
  Z_v = matrix(rnorm(M * n), M, n)
  
  for(i in 1:n){  
    
    v[,i + 1] = v[,i] + kappa * (theta - v[,i]) * Delta + sigma * sqrt(pmax(v[,i], 0 )) * Z_v[,i] * sqrt.dt
  
  }
  
  if(Z){
    
    return(list(pmax(v, 0), Z_v))
    
  }else{
    
    return(pmax(v, 0))
    
  }
  
}  

# test = DD(0.1, 3, 0.19, 0.4, T = 1/12, n = 1, M = 1000000)
# MomentsCIR(p = 1:4, kappa = 3, theta = 0.19, sigma = 0.4, v0 = 0.1, t = 1 / 12)
# Evt^1|v0     Evt^2|v0     Evt^3|v0     Evt^4|v0 
# 0.1199079295 0.0155445930 0.0021628919 0.0003210907 