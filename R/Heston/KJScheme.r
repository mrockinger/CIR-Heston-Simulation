IKJIMMScheme = function(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, v = F, Z_v = F, method = "Absorption"){
  
  Delta = T / n
  sqrt.dt = sqrt(Delta)
  
  # Define compulsory placeholders
  
  X = matrix(, M, n + 1)
  X[,1] = log(S0)
  
  # If no CIR values are fed into the function, simulate them now
  if(!is.matrix(v)){
    
    v =  matrix(, M, n + 1)
    v[,1] = v0
    
    if(!is.matrix(Z_v)){
      Z_v = matrix(rnorm(M * n), M, n)
    }
    
    if(method == "Standard"){
      
      for(i in 1:n){
        
        v[, i+1] = (v[, i] + kappa*theta*Delta + sigma*sqrt(v[, i])*sqrt.dt*Z_v[, i] + 0.25*sigma^2*Delta*(Z_v[, i]^2-1)) / (1+kappa*Delta)
        
      }
      
    }else if(method == "Absorption"){
      
      for(i in 1:n){
        
        v[, i + 1] = pmax((v[, i] + kappa * theta * Delta + sigma * sqrt(v[, i]) * sqrt.dt * Z_v[, i] + 0.25 * sigma^2 * Delta * (Z_v[, i]^2 - 1)) / (1 + kappa * Delta), 0)
        
      }
      
    }
    
  }
  
  # Simulate log-prices
  
  # If not fed into the function, simulate variance noises (should never be the case)
  if(!is.matrix(Z_v)){
    
    Z_v = matrix(rnorm(M * n), M, n)
    
  }
  
  Z_x = rho * Z_v + sqrt(1-rho^2) * matrix(rnorm(M * n), M, n) 
  
  for(i in 1:n){
    
    X[,i + 1] = X[,i] + (mu - (v[,i] + v[,i + 1]) / 4) * Delta + rho * sqrt(v[,i]) * sqrt.dt * Z_v[,i] + (sqrt(v[,i]) + sqrt(v[,i + 1])) / 2 * sqrt.dt * (Z_x[,i] - rho * Z_v[,i]) + rho * sigma / 4 * Delta * (Z_v[,i]^2 - 1)
    
  }
  
  return(X)
  
}

# set.seed(123)
# system.time({test = IKJIMMScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1, n = 12, M = 10000)})
# 
# set.seed(123)
# V = KJ(0.1, 3, 0.19, 0.4, 1, n = 12, Z = T, M = 10000)
# system.time({test2 = IKJIMMScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1, n = 12, M = 10000, v = V[[1]], Z_v = V[[2]])})
# 
# all.equal(test, test2)



# set.seed(123)
# system.time({test = IKJIMMScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1/12, n = 21, M = 10000)})
# test.r = test[,ncol(test)] - log(100)

# mean(test.r)
# mean(test.r^2)
# mean(test.r^3)
# mean(test.r^4)
# var(test.r)
# EnvStats::skewness(test.r)
# EnvStats::kurtosis(test.r, excess = F)

#  MomentsBates(0, 3, 0.19, 0.4, -0.7, 0, 0, 0, 1/12, 0.1)
# Ert.v0       Ert2.v0       Ert3.v0       Ert4.v0      Varrt.v0     Skewrt.v0     Kurtrt.v0 
# -0.0045986784  0.0093152301 -0.0004226896  0.0002827454  0.0092940822 -0.2811156802  3.1969308454 

# > mean(test.r)
# [1] -0.00455419
# > mean(test.r^2)
# [1] 0.004731807
# > mean(test.r^3)
# [1] -0.0002364582
# > mean(test.r^4)
# [1] 7.977798e-05
# > var(test.r)
# [1] 0.004711538
# > EnvStats::skewness(test.r)
# [1] -0.531999
# > EnvStats::kurtosis(test.r, excess = F)
# [1] 3.427755
