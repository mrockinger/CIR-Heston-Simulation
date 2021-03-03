# Scheme of Tse/Wan (2013)

library(Bessel)

# S0 = 100; v0 = 0.1; kappa = 3; theta = 0.19; sigma = 0.4; mu = 0; rho = -0.7; T = 1/12; n = 1; M = 10000

IGScheme = function(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, v = F, Z_v = F){
  
  # Define compulsory placeholders
  
  lnSt = matrix(, M, n + 1)
  lnSt[,1] = log(S0) 
  
  # IV = matrix(, M, n)
  
  # Define constants I will use in the remainder
  
  Delta = T / n
  
  emkt = exp(- kappa * Delta)
  h = 4 * kappa * emkt / (sigma^2 * (1 - emkt)) # is n in original paper
  d = 4 * kappa * theta / sigma^2
  nu = d / 2 - 1 
  
  N_u = 2^(15 + ceiling(log2(n))) + 1
  u = seq(0, 1, len = N_u)
  
  C1 = 1 + 2 / (exp(kappa * Delta) - 1)
  C2 = 4 / (exp(kappa * Delta / 2) - exp(- kappa * Delta / 2))^2
  Cz = 2 * kappa / (sigma^2 * sinh(kappa * Delta / 2))
  
  
  FEX1 = (C1 / kappa - Delta * C2 / 2)
  FVarX1 = (sigma^2 * C1 / kappa^3 + sigma^2 * Delta *C2 /(2 * kappa^2) - sigma^2 * Delta^2 * C1 * C2 / (2 * kappa)) 
  
  EX2 = d * sigma^2 * (- 2 + kappa * Delta * C1) /(4 * kappa^2)
  VarX2 = d * sigma^4 * (- 8 + 2 * kappa * Delta * C1 + kappa^2 * Delta^2 * C2) / (8 * kappa^4)
  
  EZ = 4 * EX2 / d
  VarZ = 4 * VarX2 / d

  
  # Define grid for Algorithm 2 and Algorithm 4
  
  v_min = 0.0001
  v_max = 8 * sigma
  
  v.v = seq(v_min, v_max, len = N_u)
  
  # Precompute quantiles of v_t|v_s via nearest-neighbor interpolation (Algorithm 2) if necessary. 
  # These quantiles will only be used, if the Poisson variate m_p calculated in the Alorithm 2 (below) is 0
  # This will take a while because of the for loop! p is generated quite fast even for the original N_u value. Can this be vectorized?

  if(!is.matrix(v)){
    
    # Here, I use shape = d / 2 since we will only use this approximation if the Possion variate in eq. (21) is 0
    
    p = pgamma(h / (2 * emkt) * v.v, shape = d / 2) 
    q = rep(NA, N_u)
    
    for(i in 1:N_u){
      
      #print(i)
      
      if(u[i] < p[1]){
        
        q[i] = 0
        
      }else{
        
        q[i] = v.v[which.min(abs(p-u[i]))] # which.min is not a vectorized function
        
      }
      
    }
    
  }
  
  # Do precomputation of EIVVsVt and VarIVVsVt on v.v for algorithm 4, that is the "fast moment calculation"
  # Note: To avoid computational cost of evaluating the modified Bessel function, calculate EIVVsVt/ VarIVVsVt only at 1 / 4 of the nodes and calculate the other points by linear interpolation.
  
  v.v.sparse = v.v[seq(1, length(v.v), by = 4)]
  
  Eeta     = Cz * v.v.sparse * BesselI(z = Cz * v.v.sparse, nu = nu + 1, expon.scaled = T) / (2 * BesselI(z = Cz * v.v.sparse, nu = nu, expon.scaled = T))
  Eeta2    = Cz^2 * v.v.sparse^2 * BesselI(z = Cz * v.v.sparse, nu = nu + 2, expon.scaled = T) / (4 * BesselI(z = Cz * v.v.sparse, nu = nu, expon.scaled = T))  + Eeta 
  
  EIVVsVt = approx(x = v.v.sparse, y = EX2 + EZ * Eeta, n = length(v.v))$y
  VarIVVsVt = approx(x = v.v.sparse, y = VarX2 + Eeta * VarZ + (Eeta2 - Eeta^2) * EZ^2, n = length(v.v))$y
  
  ################ Exit precomputations and start simulation ################
  
  # If necessary, simulate CIR process via IPZ scheme
  
  if(!is.matrix(v)){
    
    v = matrix(, M, n + 1)
    v[,1] = v0
  
    ### Simulate v_t | v_s with Algorithm 1 and 2
    
    for(j in 1:M){
      
      for(i in 1:n){
        
        if(v[j, i] == 0){
          
           # print("IPZ")
          
          U = runif(1)
          v[j, i + 1] = q[which.min(abs(U-u))]
          
        }else if(v[j, i] > 0){
          
          m_p = rpois(1, v[j, i] * h / 2)
          
          if(m_p == 0){
            
             # print("IPZ")
            
            U = runif(1)
            v[j, i + 1] = q[which.min(abs(U-u))]
            
          }else if(m_p > 0){
            
            v[j, i + 1] = 2 * emkt / h * rgamma(1, shape = m_p + d  / 2)
            
          } 
          
        }
      
      }
      
    }
  
  }  
    
  ### Simulate the IV (and thus the log-prices) via algorithms 4 and 3
  
  for(j in 1:M){
    
    ### Do the "fast moment calculation" via Algorithm 4 
    
    # Compute EX1 and Sigma2X1  
    
    SumV = v[j,][-ncol(v)] + v[j,][-1]
    
    EX1   = SumV * FEX1
    VarX1 = SumV * FVarX1

  
    # Calculate EIV and VarIV 
    # Thereby, if none of v_s or v_t are 0 use nearest neighbour interpolation on the precomputed values of EIVVsVt and VarIVVsVt
    
    EIV = VarIV = c()
    
    for(i in 1:n){
      
      if(v[j, i] == 0 || v[j, i + 1] == 0){
        
        EIV[i] = EX1[i] + EX2 
        VarIV[i] = VarX1[i] + VarX2
        
        
      }else{
        
        index = which.min(abs(sqrt(v[j, i] * v[j, i + 1]) - v.v))
        
        EIV[i] = EX1[i] + EIVVsVt[index]
        VarIV[i] = VarX1[i] + VarIVVsVt[index]
          
        
      }
      
    }
    
    ### Calculate the IV via algorithm 3
  
    N.rd = rnorm(n)
    U = runif(n)
    
    m = EIV
    s = EIV^3 / VarIV
    
    x = 1 + N.rd^2 / (2 * s / m) - sqrt((2 * s / m + 2 * s / m) * N.rd^2 + (N.rd^2)^2) / (2 * s / m)
    
    
    IV = ifelse(U * (1 + x) > 1, m / x, m * x)  
    
    
    ### Do log-price simulation
    
    for(i in 1:n){
      
      lnSt[j, i + 1] = lnSt[j, i] + mu * Delta - 0.5 * IV[i] + rho / sigma * (v[j, i + 1] - v[j, i] - kappa * theta * Delta + kappa * IV[i]) + sqrt((1 - rho^2) * IV[i]) * rnorm(1) 
      
    }
    
  }  
  
  return(lnSt)
    
}

# set.seed(123)
# system.time({test = IGScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1 / 12, n = 1, M = 10000)})
# 
# set.seed(123)
# v = IPZScheme(v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, T = 1 / 12, n = 1, M = 10000)
# system.time({test2 = IGScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1 / 12, n = 1, M = 10000, v = v)})
# 
# all.equal(test, test2)



# "Fair" parameters

# set.seed(123)
# system.time({test = IGScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1 / 12, n = 1, M = 10000)})



# test.r = test[,2] - log(100)
# plot(density(test.r))
# 
# mean(test.r)
# mean(test.r^2)
# mean(test.r^3)
# mean(test.r^4)
# EnvStats::skewness(test.r)
# EnvStats::kurtosis(test.r, excess = F)
# 
# MomentsBates(0, 3, 0.19, 0.4, -0.7, 0, 0, 0, 1 / 12, 0.1)
# Ert.v0       Ert2.v0       Ert3.v0       Ert4.v0      Varrt.v0     Skewrt.v0     Kurtrt.v0 
# -0.0045986784  0.0093152301 -0.0004226896  0.0002827454  0.0092940822 -0.3285375760  3.1969308454 

# > mean(test.r)
# [1] -0.00454238
# > mean(test.r^2)
# [1] 0.009200512
# > mean(test.r^3)
# [1] -0.0003923663
# > mean(test.r^4)
# [1] 0.0002762398
# > EnvStats::skewness(test.r)
# [1] -0.3038148
# > EnvStats::kurtosis(test.r, excess = F)
# [1] 3.207632

# "Tough" parameters

# system.time({test = IGScheme(S0 = 100, v0 = 0.1, kappa = 0.5, theta = 0.19, sigma = 1, mu = 0, rho = -0.7, T = 1/12, n = 1, M = 10000)})



# test.r = test[,2] - log(100)
# plot(density(test.r))
# 
# mean(test.r)
# mean(test.r^2)
# mean(test.r^3)
# mean(test.r^4)
# EnvStats::skewness(test.r)
# EnvStats::kurtosis(test.r, excess = F)

# MomentsBates(0, 0.5, 0.19, 1, -0.7, 0, 0, 0, 1 / 12, 0.1)
# Ert.v0       Ert2.v0       Ert3.v0       Ert4.v0      Varrt.v0     Skewrt.v0     Kurtrt.v0 
# -0.0042437178  0.0087528535 -0.0008825076  0.0003632022  0.0087348444 -0.9447114876  4.5763823279 

# > mean(test.r)
# [1] -0.005482966
# > mean(test.r^2)
# [1] 0.008856072
# > mean(test.r^3)
# [1] -0.0009193924
# > mean(test.r^4)
# [1] 0.0003653529
# > EnvStats::skewness(test.r)
# [1] -0.933657
# > EnvStats::kurtosis(test.r, excess = F)
# [1] 4.453067
