source("GlassermanKimAddendum.r")


# Define the Gamma expansion scheme

GEScheme = function(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, v = F, Z_v = F){
  
  # S0 = 100; v0 = 0.1; kappa = 3; theta = 0.19; sigma = 0.4; mu = 0; rho = -0.7; T = 1/252; n = 1; M = 1000
  # S0 = 100; v0 = 0.1; kappa = 0.5; theta = 0.19; sigma = 1; mu = 0; rho = -0.7; T = 30; n = 30 * 12; M = 500
  
  # Define placeholder and noise for returns
  
  lnSt = matrix(, M, n + 1)
  lnSt[,1] = log(S0) 
  
  
  ### Define constants for variance simulation
  
  Delta = T / n
  
  emkt = exp(- kappa * Delta)
  d    = 4 * theta * kappa / sigma^2
  h    = 4 * kappa / (sigma^2 * (1 - emkt))
  
  
  ### Define some constants 
  
  nu = d / 2 - 1
  
  # C1 = 1 + 2 / (exp(kappa * Delta) - 1) # coth(kappa * Delta / 2)
  # try this to get closer to matlab results:
  C1 = cosh(kappa * Delta / 2) / sinh(kappa * Delta / 2)
  C2 = 4 / (exp(kappa * Delta / 2) - exp(- kappa * Delta / 2))^2 # csch^2(kappa * Delta / 2)
  # 2 / (cosh(kappa * Delta) - 1)
  Cz = (2 * kappa / sigma^2) / sinh(kappa * Delta / 2)
  
  ### Define stuff for simulating X1 via truncated Gamma expansion
  
  K = 10 # Truncation level
  n.X1 = 1:K
  
  # Note: The intensity lambda as well as the variate gamma are dependend on the summation index n given above
  
  lambda_n =  16 * pi^2 * n.X1^2  / (sigma^2 * Delta * (kappa^2 * Delta^2 + 4 * pi^2 * n.X1^2)) 
  gamma_n =  (kappa^2 * Delta^2 + 4 * pi^2 * n.X1^2) / (2 * sigma^2 * Delta^2)
  
  
  ### Precalculate a table of cdf values of X2 and Z 
  # Thereby, approximate the cdf via Abate/Whitts POISSON algorithm
  
  # muX1_star = C1 / kappa - Delta * C2 / 2
  # sigma2X1_star = sigma^2 * C1 / kappa^3 + sigma^2 * Delta * C2 /(2 * kappa^2) - sigma^2 * Delta^2 * C1 * C2 / (2 * kappa) 
  
  muX2_star = sigma^2 / (4 * kappa^2) * (-2 + kappa * Delta * C1)
  sigma2X2_star = sigma^4 / (8 * kappa^4) * (-8 + 2 * kappa * Delta * C1 + kappa^2 * Delta^2 * C2)
  
  muX2 = d * muX2_star
  sigma2X2 = d *sigma2X2_star
  
  muZ = 4 * muX2_star
  sigma2Z =  4 * sigma2X2_star
  
  # Define some terms for the POISSON algorithm 
  
  u_epsilon_X2 = muX2 + 12 * sqrt(sigma2X2) 
  u_epsilon_Z = muZ + 12 * sqrt(sigma2Z) 
  
  M.grid = 200
  w = 0.01
  index_CDF_eval = 1:(M.grid+1)
  X2_CDF_eval_grid = w * muX2 + (index_CDF_eval - 1) / M.grid * (u_epsilon_X2 - w * muX2)
  Z_CDF_eval_grid = w * muZ + (index_CDF_eval - 1) / M.grid * (u_epsilon_Z - w * muZ)
  
  # Do the CDF tabulation for X2 and Z on the given grids
  
  CDF.X2 = FX2(X2_CDF_eval_grid, kappa, sigma, Delta, d, u_epsilon_X2) 
  CDF.Z = FZ(Z_CDF_eval_grid, kappa, sigma, Delta, u_epsilon_Z) 
  
  # Define the linear interpolation rule
  
  CDF.X2.li = approxfun(CDF.X2 ~ X2_CDF_eval_grid)
  CDF.Z.li = approxfun(CDF.Z ~ Z_CDF_eval_grid)
  
  
  
  ### Calculate m cutpoints/pointers for the cutpoint method used to generate draws of X2 and Z
  
  m = 100 # no. of cutpoints for numerically inverting the CDFs of X2 and Z
  
  # X2
  
  # I have to define some domain a<= i <=b where a and b are the starting and end points of the sparse grid respectively
  # In Fishman, this would be a discrete space, where in our case it should be continuous (CM is for sampling from DISCRETE distributions!)
  # Note: The length / precision of the artificial "domains" is crucial to computing time when X2/Z is drawn!
  # Note: There is a merge of the CPM and an acceptance - rejection method (Ahrens (1993)) in sec. 3.6
  
  # X2
  
  i.X2 = seq(X2_CDF_eval_grid[1], X2_CDF_eval_grid[M.grid+1], len = 1e4) 
  q.i.X2 = CDF.X2.li(i.X2)
  
  I.X2 = c()
  for(j in 1:m){
    
    I.X2[j] = min(i.X2[which(q.i.X2 > (j - 1) / m)])
    
  }
  I.X2[m + 1] = X2_CDF_eval_grid[M.grid + 1]
  
  # Z
  
  i.Z = seq(Z_CDF_eval_grid[1], Z_CDF_eval_grid[M.grid+1], len = 1e4)
  q.i.Z = CDF.Z.li(i.Z)
  
  I.Z = c()
  for(j in 1:m){
    
    I.Z[j] = min(i.Z[which(q.i.Z > (j - 1) / m)])
    
  }
  I.Z[m + 1] = Z_CDF_eval_grid[M.grid + 1]
  
  
  ############################# Start simulation #####################################

  ############ If necessary, simulate CIR process. Do so by exact method
  
  if(!is.matrix(v)){
    
    v = matrix(, M, n + 1)
    v[,1] = v0
    
    if(d > 1){
      
      if(!is.matrix(Z_v)){
        
        Z_v = matrix(rnorm(M * n), M, n)
        
      }
      
      X2.rv = matrix(rchisq(n * M, df = d - 1), M, n)
      
      
      for(i in 1:n){
        
        lambda = h * v[, i] * emkt
        # Here, important to omit parameter ncp. rchisq will use the algorithm from rgamma. 
        v[,i + 1] = h^-1 * ((Z_v[,i] + sqrt(lambda))^2 + X2.rv[,i])
        
        
      }
      
    }else{
      
      for(i in 1:n){
        # Used algorithm corresponds to drawing a Poisson mixture of central chi-squared RVs
        v[,i + 1] = h^-1 * sapply(h * v[,i] * emkt, rchisq, n = 1, df = d)
        
      }
      
    }
    
  }
  
  
  ############# Start the Gamma expansion (GE) scheme and thus log-price simulation
  
  for(j in 1:M){
    
    for(i in 1:n){
      
      # Prevent simulation from generating NaNs because of division through 0
      
      if(v[j, i] <= 0) v[j,i] = 1e-10
      if(v[j, i+1] <= 0) v[j, i+1] = 1e-10
      
      ##### Simulate X_1 through truncated Gamma expansion
      
      # Mean and variance of remainder term (truncation correction). See Lemma 3.1 on p. 276
      
      EX1K = 2 * (v[j, i] + v[j, i+1]) * Delta / (pi^2 * K)
      VarX1K = 2 * (v[j, i] + v[j, i+1]) * sigma^2 * Delta^3 / (3 * pi^4 * K^3)
      
      # Draw the number of exponential r.v.s N_n for each summation index n = 1,...,K from a Poisson distribution 
      # Note further: Intensity (mean) of the Poisson r.v. is a function of the summation index n
      
      N_n = sapply((v[j, i] + v[j, i+1]) * lambda_n, rpois, n = 1) 
      
      # Calculate the SUM of N_n unit mean exponential distributed variates for each n = 1,...,K -> so SumExp will have n elements.
      
      SumExp = unlist(lapply(sapply(N_n, rexp, rate = 1), sum))
      
      # Approximating the remainder term with a gamma random variable which has the same mean and variance as the true remainder term (EX1K, VarX1K)
      
      X1 = sum(gamma_n^-1 * SumExp) + rgamma(1, shape = EX1K^2 / VarX1K, scale = VarX1K / EX1K)
      
      ##### Simulate X_2

      U.X2 = runif(1)
      L.cp.X2 = ceiling(U.X2 * m)
      X2 = I.X2[L.cp.X2]
      
      # N = 0
      while(U.X2 > CDF.X2.li(X2)){
        
        # N = N + 1
        # print(N)
        X2 = i.X2[which(i.X2 == X2) + 1]
        
      }

      ##### Simulate X_3
      
      # First, calculate eta, the number of r.v.'s Z, we have to generate
      
      z = Cz * sqrt(v[j, i] * v[j, i + 1])
      
      U.eta = runif(1)
      
      eta = 1

      while(U.eta -  sum( pdf.eta.proper.2(n = 0:(eta-1), nu = nu, z = z) ) > pdf.eta.proper.2(n = eta, nu = nu, z = z)){ 
        
        eta = eta + 1
        
      }

      # If we know the number of Z we have to draw, proceed as in X2, i. e. draw from the pre-computed distribution table
      
      Z = c()
      
      for(k in 1:eta){
        
        U.Z = runif(1)
        L.cp.Z = ceiling(U.Z * m)
        Z[k] = I.Z[L.cp.Z]
        
        while(U.Z > CDF.Z.li(Z[k])){
          
          Z[k] = i.Z[which(i.Z == Z[k]) + 1]
          
        }
        
      }
      
      X3 = sum(Z)
      
      
      #### Calculate log-price. Thereby I set IV = X1 + X2 + X3 for each step
      
      Z_x = matrix(rnorm(n * M), M, n) # Define noise for log-prices
      
      lnSt[j, i + 1] = lnSt[j, i] + mu * Delta - 0.5 * (X1 + X2 + X3) + kappa * rho / sigma * (X1 + X2 + X3) + rho / sigma * (v[j, i+1] - v[j, i] - kappa * theta * Delta) + sqrt((1 - rho^2) * (X1 + X2 + X3)) * Z_x[j, i] 
      
      
    }
    
  }  

  return(lnSt)
  
}


# set.seed(123)
# system.time({test = GEScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1, n = 12, M = 100)})
# 
# set.seed(123)
# V = Exact(0.1, 3, 0.19, 0.4, T = 1, n = 12, M = 100)
# system.time({test2 = GEScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1, n = 12, M = 100, v = V)})
# 
# all.equal(test, test2)

# system.time({test = GEScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1/12, n = 1, M = 10000)})

# test.r = test[,ncol(test)] - log(100)
# plot(density(test.r))

# mean(test.r)
# mean(test.r^2)
# mean(test.r^3)
# mean(test.r^4)
# var(test.r)
# EnvStats::skewness(test.r)
# EnvStats::kurtosis(test.r, excess = F)

# MomentsBates(0, 3, 0.19, 0.4, -0.7, 0, 0, 0, 1 / 12, 0.1)
# Ert.v0       Ert2.v0       Ert3.v0       Ert4.v0      Varrt.v0     Skewrt.v0     Kurtrt.v0
# -0.0045986784  0.0093152301 -0.0004226896  0.0002827454  0.0092940822 -0.3285375760   3.1969308454

# > mean(test.r)
# [1] -0.004931132
# > mean(test.r^2)
# [1] 0.009356741
# > mean(test.r^3)
# [1] -0.0004085293
# > mean(test.r^4)
# [1] 0.0002826897
# > var(test.r)
# [1] 0.009333359
# > EnvStats::skewness(test.r)
# [1] -0.2999173
# > EnvStats::kurtosis(test.r, excess = F)
# [1] 3.169612

# system.time({test = GEScheme(S0 = 100, v0 = 0.1, kappa = 0.5, theta = 0.19, sigma = 1, mu = 0, rho = -0.7, T = 1/12, n = 1, M = 10000)})

# test.r = test[,ncol(test)] - log(100)
# plot(density(test.r))

# mean(test.r)
# mean(test.r^2)
# mean(test.r^3)
# mean(test.r^4)
# var(test.r)
# EnvStats::skewness(test.r)
# EnvStats::kurtosis(test.r, excess = F)

# MomentsBates(0, 0.5, 0.19, 1, -0.7, 0, 0, 0, 1 / 12, 0.1)
# Ert.v0       Ert2.v0       Ert3.v0       Ert4.v0      Varrt.v0     Skewrt.v0     Kurtrt.v0 
# -0.0042437178  0.0087528535 -0.0008825076  0.0003632022  0.0087348444 -0.9447114876  4.5763823279 

# > mean(test.r)
# [1] -0.003720025
# > mean(test.r^2)
# [1] 0.008574424
# > mean(test.r^3)
# [1] -0.0007912285
# > mean(test.r^4)
# [1] 0.0003361253
# > var(test.r)
# [1] 0.008561442
# > EnvStats::skewness(test.r)
# [1] -0.8784047
# > EnvStats::kurtosis(test.r, excess = F)
# [1] 4.437004
