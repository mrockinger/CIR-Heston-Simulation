# Exact scheme from Broadie/Kaya (2006)

source("BroadieKayaSchemeAddendum.r")


BKScheme = function(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, v = F, Z_v = F){
  
  
  # Define compulsory placeholders
  
  lnSt = matrix(, M, n + 1)
  lnSt[,1] = log(S0) 
  
  # Define some constants
  
  Delta = T / n
  emkt   = exp(- kappa * Delta)
  nu     = 4 * theta * kappa / sigma^2
  h = 4 * kappa / (sigma^2 * (1 - emkt))
  
  
  # If necessary, simulate CIR process. Do so by exact method
  
  if(!is.matrix(v)){
    
    v = matrix(, M, n + 1)
    v[,1] = v0
    
    if(nu > 1){
      
      if(!is.matrix(Z_v)){
        
        Z_v = matrix(rnorm(M * n), M, n)
        
      }
      
      # Here, important to omit parameter ncp. rchisq will use the algorithm from rgamma. 
      X2.rv = matrix(rchisq(n * M, df = nu - 1), M, n)
      
      
      for(i in 1:n){
        
        lambda = h * v[, i] * emkt
        v[,i + 1] = h^-1 * ((Z_v[,i] + sqrt(lambda))^2 + X2.rv[,i])
        
        
      }
      
    }else{
      
      for(i in 1:n){
        # Used algorithm corresponds to drawing a Poisson mixture of central chi-squared RVs
        v[,i + 1] = h^-1 * sapply(h * v[,i] * emkt, rchisq, n = 1, df = nu)
        
      }
      
    }
    
  }
  
  # Step 2 of Heston simulation: Sample X_t given X_s, v_s, v_t
  
  Z_x = matrix(rnorm(n * M), M, n)
  U = matrix(runif(M * n), M, n)
  
  for(i in 1:n){
    
    # Firt, generate a draw of \int_s^t v_s ds given v_s and v_t
    IV = GenerateIV(Vs = v[, i], Vt = v[, i+1], Delta = Delta, kappa = kappa, theta = theta, sigma = sigma, U = U[, i])
    
    # This could be also written like LnSt = rnorm(mu, sigma^2) where mu and sigma^2 are defined as in section 3.3 of the original paper
    lnSt[, i + 1] = lnSt[, i] + mu * Delta - 0.5 * IV + rho / sigma * (v[, i + 1] - v[, i] - kappa * theta * Delta + kappa * IV) + sqrt((1 - rho^2) * IV) * Z_x[, i] 
    
  }
  
  
  return(lnSt)
  
}

# set.seed(123)
# system.time({test = BKScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1, n = 12, M = 100)})
# 
# 
# set.seed(123)
# V = Exact(0.1, 3, 0.19, 0.4, T = 1, n = 12, M = 100)
# system.time({test2 = BKScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1, n = 12, M = 100, v = V)})
# 
# all.equal(test, test2)

# "Fair" parameters

# set.seed(123)
# system.time({test = BKScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1/12, n = 1, M = 10000)})

# test.r = test[,2] - log(100)
# plot(density(test.r))

# mean(test.r)
# mean(test.r^2)
# mean(test.r^3)
# mean(test.r^4)
# var(test.r)
# EnvStats::skewness(test.r)
# EnvStats::kurtosis(test.r, excess = F)


#  MomentsBates(0, 3, 0.19, 0.4, -0.7, 0, 0, 0, 1 / 12, 0.1)
# Ert.v0       Ert2.v0       Ert3.v0       Ert4.v0      Varrt.v0     Skewrt.v0     Kurtrt.v0 
# -0.0045986784  0.0093152301 -0.0004226896  0.0002827454  0.0092940822 -0.2811156802  3.1969308454 

# > mean(test.r)
# [1] -0.00459868
# > mean(test.r^2)
# [1] 0.009197773
# > mean(test.r^3)
# [1] -0.000441345
# > mean(test.r^4)
# [1] 0.0002774504
# > var(test.r)
# [1] 0.009177543
# > EnvStats::skewness(test.r)
# [1] -0.3579842
# > EnvStats::kurtosis(test.r, excess = F)
# [1] 3.21287



# Tough parameters

# system.time({test = BKScheme(S0 = 100, v0 = 0.1, kappa = 0.5, theta = 0.19, sigma = 1, mu = 0, rho = -0.7, T = 1/12, n = 1, M = 10000)})
# User      System verstrichen 
# 363.65        0.36      364.42 

# test.r = test[,2] - log(100)
# plot(density(test.r))
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
# [1] -0.004378141
# > mean(test.r^2)
# [1] 0.008943675
# > mean(test.r^3)
# [1] -0.0008899693
# > mean(test.r^4)
# [1] 0.0003638285
# > skewness(test.r)
# [1] -0.916603
# > kurtosis(test.r, excess = F)
# [1] 4.386533

