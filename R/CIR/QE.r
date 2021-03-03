# use quadratic exponential scheme from Andersen (2008) to discretize the CIR process

QE = function(v0, kappa, theta, sigma, T, n, M){
    
    Delta = T / n
    
    U = matrix(runif(M * n), M, n)
    
    psic = 1.5
    
    emkt = exp(- kappa * Delta)
    c1  = sigma * sigma * emkt * (1 - emkt) / kappa
    c2  = theta * sigma * sigma * ((1 - emkt)^2) / (2 * kappa)
    
    v   = matrix(, M, n + 1)
    v[, 1]= v0

    for(i in 1:n){
      
      s2  = v[, i] * c1 + c2
      m   = v[, i] * emkt + theta * (1 - emkt)
      psi = s2 / (m * m)
      
      v[,i + 1] = 
      
        ifelse(
          
          psi <= psic,
          
          # a * (b + Z_v)^2, where a = m / (1 + b^2)
          m / (1 + 2 * psi^-1 - 1 + sqrt(2 * psi^-1) * sqrt(2 * psi^-1 - 1))  * (sqrt(2 * psi^-1 - 1 + sqrt(2 * psi^-1) * sqrt(2 * psi^-1 - 1)) + qnorm(U[, i]))^2,
          
          # p = (psi - 1) / (psi + 1) 
          ifelse(U[, i] < (psi - 1) / (psi + 1), 0, ((1 - (psi - 1) / (psi + 1)) / m)^-1 * log( (1 - (psi - 1) / (psi + 1)) / (1 - U[,i]) ))
          
        )
      
    } 
       
    # Old, not vectorized code:
    
    #   if(psi < psic){
    #     
    #     psii = 1 / psi
    #     b2   = 2 * psii - 1 + sqrt(2 * psii * (2 * psii - 1))
    #     a    = m / (1 + b2)
    #     v[i] = a * (sqrt(b2) + qnorm(U[i - 1]))^2
    #     
    #   }else{
    #     
    #     p    = (psi - 1) / (psi + 1)
    #     beta = (1 - p) / m
    #     
    
    
    #     if(U[i - 1] < p){
    #       
    #       v[i] = 0
    #       
    #     }else{
    #       
    #       v[i] = (1 / beta) * log( (1 - p) / (1 - U[i - 1]) )
    #       
    #     }
    #   }
    # }
    
    return(v)
    
}

# test = QE(0.04, 0.5, 0.04, 1, 1, 1, 1000000)
# > MomentsCIR(p = 1:4, kappa = 0.5, theta = 0.04, sigma = 1, v0 = 0.04, t = 1)
# Evt^1|v0   Evt^2|v0   Evt^3|v0   Evt^4|v0 
# 0.04000000 0.02688482 0.03050794 0.04777093 
