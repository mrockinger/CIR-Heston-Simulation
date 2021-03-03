# Approximation of the non-central chi-squared via Sankaran (1963)

Sankaran = function(v0, kappa, theta, sigma, T, n, M, mode){
  
  Delta = T / n
  emkt = exp(- kappa * Delta)
  d = 4 * kappa * theta / sigma^2
  h = 4 * kappa / (sigma^2 * (1 - emkt))
  v = matrix(, M, n + 1)
  v[,1] = v0
  
  
  Z = matrix(rnorm(M * n, 0, 1), M, n)
  
  if(mode == "simple"){ #corresponds to equation (12)
    
    for(i in 1:n){
      
      lambda = h * v[,i] * emkt
      v[,i + 1] = h^-1 * (0.5 * (d - 1) + (sqrt(lambda + 0.5 * (d - 1)) + Z[, i])^2)
      
    }
    
  }else if(mode == "normal"){ #corresponds to equation (7)
    
    
    for(i in 1:n){
      
      lambda = h * v[,i] * emkt
      r = d + lambda
      k1 = 1 - (d + 2) / (6 * r) - 1 / (72 * r^2) * (d^2 - 2 * d + 10) - 1 / (432 * r^3) * (d^3 - 12 * d^2 - 6 * d + 44) - 5 / (10368 * r^4) * (d^4 - 28 * d^3 + 24 * d^2 + 1112 * d - 1028) 
      k2 = 1 / r * ( 1 - (d - 1) / (6 * r) - 1 / (18 * r^2) * (d^2 + d - 2) - 1 / (216 * r^3) * (4 * d^3 - 9 * d^2 - 228 * d + 233) )
      
      
      # ifelse(k2 < 0, stop("k2 < 0"))
      
      v[,i + 1] = h^-1 * (1 / 3 * (d - 1) + (d + lambda) * (k1 + Z[, i] * sqrt(k2))^2)
      

    }
    
  }
  
  return(v)
  
}

# test.SANK1 = Sankaran(0.1, 3, 0.19, 0.4, 1 / 12, 1, 1000000, mode = "simple")
# test.SANK2 = Sankaran(0.1, 3, 0.19, 0.4, 1 / 12, 1, 1000000, mode = "normal")
# MomentsCIR(p = 1:4, kappa = 3, theta = 0.19, sigma = 0.4, v0 = 0.1, t = 1 / 12)
# Evt^1|v0     Evt^2|v0     Evt^3|v0     Evt^4|v0 
# 0.1199079295 0.0155445930 0.0021628919 0.0003210907 