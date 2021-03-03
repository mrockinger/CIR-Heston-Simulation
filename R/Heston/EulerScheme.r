EulerScheme = function(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, v = F, Z_v = F, mode = "FT"){
  
      Delta = T / n
      sqrt.dt = sqrt(Delta)
      
      X = matrix(, M, n + 1)
      X[,1] = log(S0) 
      
      # If no CIR values are fed into the function, simulate them now
      if(!is.matrix(v)){
        
        v =  matrix(, M, n + 1)
        v[,1] = v0
        
        if(!is.matrix(Z_v)){
          Z_v = matrix(rnorm(M * n), M, n)
        }

        if(mode == "naive"){
          
          for(i in 1:n){  
            
            v[,i + 1] = v[,i] + kappa * (theta - v[,i]) * Delta + sigma * sqrt(v[,i]) * Z_v[,i] * sqrt.dt
            
          }
      
          
        }else if(mode == "absorp"){
          
          for(i in 1:n){  
            
            v[,i + 1] =  pmax(v[,i], 0) + kappa * (theta - pmax(v[,i], 0)) * Delta + sigma * sqrt(pmax(v[,i], 0)) * Z_v[,i] * sqrt.dt
            
          }
          
          v = pmax(v,0)  
          
          
        }else if(mode == "reflect"){
          
          for(i in 1:n){  
            
            v[,i + 1] = abs(v[,i]) + kappa * (theta - abs(v[,i])) * Delta + sigma * sqrt(abs(v[,i])) * Z_v[,i] * sqrt.dt
            
          }
          

          v = abs(v, 0)
          
          
        }else if(mode == "FT"){
          
          for(i in 1:n){  
            v[,i + 1] = v[,i] + kappa * (theta - pmax(v[,i],0)) * Delta + sigma * sqrt(pmax(v[,i],0)) * Z_v[,i] * sqrt.dt
            
          }
          

          v = pmax(v, 0)
          
        } 
        
        
      }
        
      

      # Generate log-price noise  
      
      # If not fed into the function, simulate variance noises (should never be the case)
      if(!is.matrix(Z_v)){
        
        Z_v = matrix(rnorm(M * n), M, n)
        
      }
      
      Z_x = rho * Z_v + sqrt(1-rho^2) * matrix(rnorm(M * n), M, n) 
        
      for(i in 1:n){
          
        X[, i+1] <- X[, i] + (mu - v[, i]/2) * Delta + sqrt(v[, i]) * sqrt.dt * Z_x[, i]
          
      }
              
    return(X)  
      
}

# test = EulerScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1/12, n = 5, M = 10000, mode = "reflect")
# 
# test.r = test[,ncol(test)] - log(100)
# 
# mean(test.r)
# mean(test.r^2)
# mean(test.r^3)
# mean(test.r^4)
# var(test.r)
# skewness(test.r)
# kurtosis(test.r, excess = F)


#  MomentsBates(0, 3, 0.19, 0.4, -0.7, 0, 0, 0, 1/12, 0.1)
# Ert.v0       Ert2.v0       Ert3.v0       Ert4.v0      Varrt.v0     Skewrt.v0     Kurtrt.v0 
# -0.0045986784  0.0093152301 -0.0004226896  0.0002827454  0.0092940822 -0.2811156802  3.1969308454 

# > mean(test.r)
# [1] -0.005172206
# > mean(test.r^2)
# [1] 0.00836624
# > mean(test.r^3)
# [1] -0.000148341
# > mean(test.r^4)
# [1] 0.0002151312
# > var(test.r)
# [1] 0.008340323
# > skewness(test.r)
# [1] -0.02469217
# > kurtosis(test.r, excess = F)
# [1] 3.069102


# set.seed(123)
# test = EulerScheme(S0 = 100, v0 = 0.1, kappa = 0.5, theta = 0.19, sigma = 1, mu = 0, rho = -0.7, T = 1/12, n = 5, M = 10000, mode = "FT")
# 
# test.r = test[,ncol(test)] - log(100)
# 
# mean(test.r)
# mean(test.r^2)
# mean(test.r^3)
# mean(test.r^4)
# var(test.r)
# EnvStats::skewness(test.r)
# EnvStats::kurtosis(test.r, excess = F)

# The following is equal to the above which is what we want

# set.seed(123)
# V = EulerCIR(0.1, 0.5, 0.19, 1, 1 / 12, 5, 10000, Z = T, mode = "FT")
# X = EulerScheme(100, 0.1, 0.5, 0.19, 1, 0, -0.7, 1/12, 5, 10000, v = V[[1]], Z_v = V[[2]], mode = "FT")

# all.equal(X, test)

# MomentsBates(0, 0.5, 0.19, 1, -0.7, 0, 0, 0, 1 / 12, 0.1)
# Ert.v0       Ert2.v0       Ert3.v0       Ert4.v0      Varrt.v0     Skewrt.v0     Kurtrt.v0 
# -0.0042437178  0.0087528535 -0.0008825076  0.0003632022  0.0087348444 -0.9004037725  4.5763823279

# > mean(test.r)
# [1] -0.003436663
# > mean(test.r^2)
# [1] 0.008522391
# > mean(test.r^3)
# [1] -0.0006517223
# > mean(test.r^4)
# [1] 0.0002843367
# > var(test.r)
# [1] 0.008511431
# > EnvStats::skewness(test.r)
# [1] -0.7183855
# > EnvStats::kurtosis(test.r, excess = F)
# [1] 3.811324