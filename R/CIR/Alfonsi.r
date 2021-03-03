# Simulate CIR process via Alfonsi's (2005) explicit scheme E, see formula (5)
# The explicit scheme is nonnegativ for 0 <= lambda <= kappa * theta - sigma^2 / 4
# To have nonnegativity for all parameter combinations, I also implemented the "trick" mentioned at the beginning of section 5
# Note: I do the fix as stated in this paper NOT as in the Lord et al. I. e. I do not simulate an an independent auxiliary process which is pushed into  the right domain after the recursion but I do fix immediately
# See p. 373 at the bottom 

AE0 = function(v0, kappa, theta, sigma, T, n, M, Z = F){
  
  #if(4 * kappa * theta < sigma^2) print("Warning. 4 * kappa * theta >= sigma^2 not fullfilled")
  
  v = matrix(, M, n + 1)
  v[, 1] = v0
  
  
  Delta = T / n
  sqrt.dt = sqrt(Delta)
  f1mkT2n = 1 - 0.5 * kappa * Delta
  a = kappa * theta
  
  Z_v = matrix(rnorm(M * n), M, n)
  
  for(i in 1:n){
    
     v[, i + 1] = 
       
       pmax((f1mkT2n * sqrt(v[, i]) + sigma * Z_v[,i] * sqrt.dt / (2 * f1mkT2n))^2 + (a - sigma^2 / 4) * Delta, 0)
    
 }
  
  if(Z){
    
    return(list(v, Z_v))
    
  }else{
    
    return(v)
    
  }
  

}





# test = AE0(0.1, 3, 0.19, 0.4, 1 / 12, 1, 1000000)
# MomentsCIR(p = 1:4, kappa = 3, theta = 0.19, sigma = 0.4, v0 = 0.1, t = 1 / 12)
# Evt^1|v0     Evt^2|v0     Evt^3|v0     Evt^4|v0 
# 0.1199079295 0.0155445930 0.0021628919 0.0003210907 

# To be sure, compute it like Michael does in matlab (I think some steps are redundant) 

# AE02 = function(v0, kappa, theta, sigma, T, n, M, Z = F){
#   
#   #if(4 * kappa * theta < sigma^2) print("Warning. 4 * kappa * theta >= sigma^2 not fullfilled")
#   
#   vtD = matrix(, M, n + 1)
#   vt = v0
#   vtD[, 1] = vt
#   
#   
#   
#   Delta = T / n
#   sqrt.dt = sqrt(Delta)
#   f1mkT2n = 1 - 0.5 * kappa * Delta
#   a = kappa * theta
#   
#   Z_v = matrix(rnorm(M * n), M, n)
#   
#   for(i in 1:n){
#     
#     vt = (f1mkT2n * sqrt(vtD[, i]) + sigma * Z_v[,i] * sqrt.dt / (2 * f1mkT2n))^2 + (a - sigma^2 / 4) * Delta
#     
#     vtD[, i + 1] = pmax(vt, 0)
#       
#   }
#   
#   if(Z){
#     
#     return(list(vtD, Z_v))
#     
#   }else{
#     
#     return(vtD)
#     
#   }
#   
#   
# }


# > set.seed(123)
# > test = AE0(0.19, 3, 0.19, 0.4, 1/12, 10, 50000)
# > set.seed(123)
# > test2 = AE02(0.19, 3, 0.19, 0.4, 1/12, 10, 50000)
# > all.equal(test, test2)
# [1] TRUE