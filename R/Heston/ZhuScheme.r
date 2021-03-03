# Scheme from Zhu (2011): A Simple and Accurate Simulation Approach to the Heston model in J. of Derivatives
# "Transformed volatility scheme" 

ZhuScheme = function(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, v = F, Z_v = F){
  
  Delta = T / n
  sqrt.dt = sqrt(Delta)
  
  lnSt = matrix(, M, n + 1) 
  lnSt[,1] = log(S0) 
  
  # Simulate variance if necessary. Thereby, use the TVS scheme, i. e. simulate the volatility instead of the variance (therfore the v = v^2 at the end)
  if(!is.matrix(v)){
    
    v = matrix(, M, n + 1) # plays role of volatility for now
    v[,1] = sqrt(v0)
    
    if(!is.matrix(Z_v)){
      Z_v = matrix(rnorm(M * n), M, n)
    }

    
    emkt = exp(- kappa * Delta)
    m2 = sigma^2 / (4 * kappa) * (1 - emkt) # variance of VOLATILITY
  
    for(i in 1:n){
    
      m1 = theta + (v[, i]^2 - theta) * emkt #expected VARIANCE
      
      beta = sqrt(pmax(m1 - m2, 0))
      theta_v_star = (beta - v[, i] * emkt^(1/2)) / (1 - emkt^(1/2))
      
      v[, i + 1] = v[, i] + 1 / 2 * kappa * (theta_v_star - v[, i]) * Delta + 1 / 2 * sigma * sqrt.dt * Z_v[, i]
    
    
    }
    
    v = v^2  
  
  }
    
  # Do log-price simulation
  
  # If not fed into the function, simulate variance noises (should never be the case)
  if(!is.matrix(Z_v)){
    Z_v = matrix(rnorm(M * n), M, n)
  }
    
  Z_x = rho * Z_v + sqrt(1-rho^2) * matrix(rnorm(M * n), M, n) 

  for(i in 1:n){

    lnSt[, i+1] = lnSt[, i] + (mu - 1 / 2 * v[, i]) * Delta + sqrt(v[, i]) * sqrt.dt * Z_x[, i]
      
      
  }

  return(lnSt)
  
}

# set.seed(123)
# test = ZhuScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1/12, n = 1, M = 10000)
# 
# set.seed(123)
# V = TVSScheme(v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, T = 1/12, n = 1, M = 10000, Z = T)
# test2 = ZhuScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1/12, n = 1, M = 10000, v = V[[1]], Z_v = V[[2]])

# all.equal(test, test2)

# test = ZhuScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1/12, n = 1, M = 10000)
# 
# test.r = test[,2] - log(100)
# 
# mean(test.r)
# mean(test.r^2)
# mean(test.r^3)
# mean(test.r^4)
# var(test.r)
# skewness(test.r)
# kurtosis(test.r, excess = F)


#  MomentsBates(0, 3, 0.19, 0.4, -0.7, 0, 0, 0, 1/2, 0.1)
# Ert.v0       Ert2.v0       Ert3.v0       Ert4.v0      Varrt.v0     Skewrt.v0     Kurtrt.v0 
# -0.0045986784  0.0093152301 -0.0004226896  0.0002827454  0.0092940822 -0.2811156802  3.1969308454 

# > mean(test.r)
# [1] -0.002506914
# > mean(test.r^2)
# [1] 0.008338496
# > mean(test.r^3)
# [1] -5.422176e-05
# > mean(test.r^4)
# [1] 0.0002101927
# > var(test.r)
# [1] 0.008333045
# > skewness(test.r)
# [1] 0.01112277
# > kurtosis(test.r, excess = F)
# [1] 3.024898

# test = ZhuScheme(S0 = 100, v0 = 0.1, kappa = 0.5, theta = 0.19, sigma = 1, mu = 0, rho = -0.7, T = 1/12, n = 1, M = 10000)
# 
# test.r = test[,2] - log(100)
# 
# mean(test.r)
# mean(test.r^2)
# mean(test.r^3)
# mean(test.r^4)
# var(test.r)
# skewness(test.r)
# kurtosis(test.r, excess = F)

# MomentsBates(0, 0.5, 0.19, 1, -0.7, 0, 0, 0, 1 / 12, 0.1)
# Ert.v0       Ert2.v0       Ert3.v0       Ert4.v0      Varrt.v0     Skewrt.v0     Kurtrt.v0 
# -0.0042437178  0.0087528535 -0.0008825076  0.0003632022  0.0087348444 -0.9004037725  4.5763823279

# > mean(test.r)
# [1] -0.005213979
# > mean(test.r^2)
# [1] 0.00818759
# > mean(test.r^3)
# [1] -0.0001268406
# > mean(test.r^4)
# [1] 0.00019632
# > var(test.r)
# [1] 0.00816122
# > skewness(test.r)
# [1] 0.001283015
# > kurtosis(test.r, excess = F)
# [1] 2.928955


