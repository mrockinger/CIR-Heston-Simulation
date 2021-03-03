# Bégin et al. (2015): Gamma approximation of the IV

# library(Bessel)

# S0 = 100; v0 = 0.1; kappa = 3; theta = 0.19; sigma = 0.4; mu = 0; rho = -0.7; T = 1 / 12; n = 1; M = 10000

GAScheme = function(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, zeta = 52, v = F, Z_v = F){
  
  # Define compulsory placeholders 
  
  lnSt = matrix(, M, n + 1)
  lnSt[,1] = log(S0) 
  
  # Define some constants needed in the remainder
  
  Delta = T / n
  emkt = exp(- kappa * Delta)
  d = 4 * kappa * theta / sigma^2
  h = 4 * kappa / (sigma^2 * (1 - emkt))
  
  # Define constants for calculating the moments of IV via the Gamma expansion of Glasserman/Kim 
  
  C1 = 1 + 2 / (exp(kappa * Delta) - 1)
  C2 = 4 / (exp(kappa * Delta / 2) - exp(- kappa * Delta / 2))^2
  Cz = (2 * kappa / sigma^2) / sinh(kappa * Delta / 2)
  
  FEX1 = (C1 / kappa - Delta * C2 / 2)
  FVarX1 = sigma^2 * (C1 / kappa^3 + Delta *C2 /(2 * kappa^2) - Delta^2 * C1 * C2 / (2 * kappa)) 
  
  EX2 = d * sigma^2 * (- 2 + kappa * Delta * C1) /(4 * kappa^2)
  VarX2 = d * sigma^4 * (- 8 + 2 * kappa * Delta * C1 + kappa^2 * Delta^2 * C2) / (8 * kappa^4)
  
  EZ = 4 * EX2 / d
  VarZ = 4 * VarX2 / d
  
  # Define cache parameters
  # Note: delta in the paper is defined as being dependent on some parameter N, which indicates the number of simulated price paths, NOT time steps taken! For us: N = 1?!
  # Note: I think the resulting grid reaches way too high values in the right, since it will always end at 4
  
  delta = exp(log(10 * sigma) / ((75 - zeta) + log2(M)))
  xi = log(10 * sigma) / log(delta)
  
  # Precompute EIVVsVt and VarIVVsVt (that is EIV* and VarIV*) on the grid defined in Definition 4.2
  # Note: Grid is for values of vs * vt (Caution: Tse/wan (2013) made a grid for sqrt(vs * vt))
  
  grid =  c(0, delta^(-zeta + (0:zeta)), delta^(xi - ((xi-1):0))) 
  
  # Here, the first value of Eeta and Eeta2 will give NaNs since we divide through 0
  
  Eeta     = Cz * sqrt(grid) * BesselI(z = Cz * sqrt(grid), nu = d / 2 , expon.scaled = T) / (2 * BesselI(z = Cz * sqrt(grid), nu = d / 2 - 1, expon.scaled = T))
  Eeta2    = Cz^2 * grid * BesselI(z = Cz * sqrt(grid), nu = d / 2 + 1, expon.scaled = T) / (4 * BesselI(z = Cz * sqrt(grid), nu = d / 2 - 1, expon.scaled = T))  + Eeta 
  
  # Nevertheless, we know from Tse/wan (2013), that these moments are 0 anyway if vs * vt = 0
  
  Eeta[1] = Eeta2[1] = 0

  EIV_star = EX2  + Eeta * EZ 
  VarIV_star =  VarX2 + Eeta * VarZ + (Eeta2 - Eeta^2) * EZ^2
  
  # Define a linear interpolation rule which will take values of vs * vt
  
  EIV_star_li = approxfun(EIV_star ~ grid)
  VarIV_star_li = approxfun(VarIV_star ~ grid)
  
  
  ################ Exit precomputations and start simulation ################
  
  # If necessary, simulate CIR process. Do so by exact method
  
  if(!is.matrix(v)){
    
    v = matrix(, M, n + 1)
    v[,1] = v0
    
    if(d > 1){
      
      if(!is.matrix(Z_v)){zeta 
        
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
    
  ### Simulation of IV and log-prices
  
  
  rgamma.v = Vectorize(rgamma) # An auxiliary function so that I can give multiple shape and scale parameters at the same time. I checked the functionality with some toy examples. It works just the way we want it.
  
  Z_x = matrix(rnorm(n * M), M, n) # Generate log - price noise
  
  IV = matrix(, M, n)
  
  for(i in 1:n){
    
    ### Do the "fast moment calculation" (see Proposition 9)
    
    SumV  = v[,i] + v[,i + 1]
    ProdV = v[,i] * v[,i + 1]
    
    EX1      = SumV * FEX1
    VarX1    = SumV * FVarX1
    
    EIV = EIV_star_li(ProdV) + EX1
    VarIV = VarIV_star_li(ProdV) + VarX1
    
    ### Sample IV from moment matched Gamma distribution
    IV[,i] = rgamma.v(n = 1, shape = EIV^2 / VarIV, scale = VarIV / EIV)
    
    ### Do log-price simulation ###################
    lnSt[, i + 1] = lnSt[, i] + mu * Delta - 0.5 * IV[,i] + rho / sigma * (v[, i + 1] - v[, i] - kappa * theta * Delta + kappa * IV[,i]) + sqrt((1 - rho^2) * IV[,i]) * Z_x[, i] 
    
  }
  
  
  return(lnSt)
  
}

# grid und ProdV gegenüberstellen! grid ist weit weg von ProdV, was schlecht ist!

# grid
# 
# summary(ProdV)
# summary(sqrt(ProdV))
# 
# 
# Evv0 = MomentsCIR(p = 1, kappa = 3, theta = 0.19, sigma = 0.4, v0 = 0.1, t = 1 / 12)
# mean(v[,n+1])
# 
# EIVv0 = m1IV(3, 0.19, 0.4, 1 / 12, 0.1, Evv0)
# mean(IV)
# 
# 
# IV.cdf = Cdf.IV.approx(seq(0.001, 0.02, length.out = 100), 0.1, v[2,2], 1/12, 3, 0.19, 0.4, m1IV(3, 0.19, 0.4, 1 / 12, 0.1, v[2,2]) + 5 * sqrt(m2IV(3, 0.19, 0.4, 1 / 12, 0.1, v[2,2]) - m1IV(3, 0.19, 0.4, 1 / 12, 0.1, v[2,2])^2))
# plot(IV.cdf ~ (seq(0.001, 0.02, length.out = 100)), type = "l")
# 
# 
# # Gamma approximation mit den Momenten vom Grid, also via linearer Interpolation
# IV.GA = rgamma(n = 1000, shape = EIV[2]^2 / VarIV[2], scale = VarIV[2] / EIV[2])
# lines(ecdf(IV.GA), col = "red")
# 
# # mit exakten Momenten
# IV.GA.proper = rgamma(n = 1000, shape = m1IV(3, 0.19, 0.4, 1 / 12, 0.1, v[2,2])^2 / (m2IV(3, 0.19, 0.4, 1 / 12, 0.1, v[2,2]) - m1IV(3, 0.19, 0.4, 1 / 12, 0.1, v[2,2])^2), scale = (m2IV(3, 0.19, 0.4, 1 / 12, 0.1, v[2,2]) - m1IV(3, 0.19, 0.4, 1 / 12, 0.1, v[2,2])^2) / m1IV(3, 0.19, 0.4, 1 / 12, 0.1, v[2,2]))
# lines(ecdf(IV.GA.proper), col = "green") #passt ziemlich gut! D. h. die Approximation an sich ist ok!
# 
# 
# EIVv0 = m1IV(3, 0.19, 0.4, 1 / 12, 0.1, v[1,2]) # vgl. EIV[1] -> ziemlich große Abweichung!
# varIVv0 = m2IV(3, 0.19, 0.4, 1 / 12, 0.1, v[1,2]) - m1IV(3, 0.19, 0.4, 1 / 12, 0.1, v[2,2])^2 # vgl. VarIV[1] -> ziemlich große Abweichung!


# Test influence of zeta on goodness of fit in terms of moments

# set.seed(123)
# system.time({test = GAScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1, n = 400, M = 10000)})

# set.seed(123)
# system.time({test = GAScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1, n = 400, M = 10000, zeta = 85)})


# set.seed(123)
# V = Exact(0.1, 3, 0.19, 0.4, T = 1, n = 400, M = 10000)
# system.time({test2 = GAScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1, n = 400, M = 10000, v = V)})
#  
# all.equal(test, test2)

# test.r = test[,ncol(test)] - log(100)
# plot(density(test.r))
# mean(test.r)
# mean(test.r^2)
# mean(test.r^3)
# mean(test.r^4)
# var(test.r)
# EnvStats::skewness(test.r)
# EnvStats::kurtosis(test.r, excess = F)

# MomentsBates(0, 3, 0.19, 0.4, -0.7, 0, 0, 0, 1, 0.1)
# Ert.v0     Ert2.v0     Ert3.v0     Ert4.v0    Varrt.v0   Skewrt.v0   Kurtrt.v0 
# -0.08074681  0.17824452 -0.07568704  0.11869089  0.17172447 -0.47162731  3.42803659 

# > mean(test.r)
# [1] 0.02085847
# > mean(test.r^2)
# [1] 0.1603797
# > mean(test.r^3)
# [1] -0.0249363
# > mean(test.r^4)
# [1] 0.08955008
# > var(test.r)
# [1] 0.1599606
# > EnvStats::skewness(test.r)
# [1] -0.5465216
# > EnvStats::kurtosis(test.r, excess = F)
# [1] 3.599043



# set.seed(123)
# system.time({test = GAScheme(S0 = 100, v0 = 0.1, kappa = 3, theta = 0.19, sigma = 0.4, mu = 0, rho = -0.7, T = 1 / 12, n = 1, M = 10000)})


# test.r = test[,1 + 1] - log(100)
# plot(density(test.r))

# mean(test.r)
# mean(test.r^2)
# mean(test.r^3)
# mean(test.r^4)
# EnvStats::skewness(test.r)
# EnvStats::kurtosis(test.r, excess = F)

# MomentsBates(0, 3, 0.19, 0.4, -0.7, 0, 0, 0, 1 / 12, 0.1)
# Ert.v0       Ert2.v0       Ert3.v0       Ert4.v0      Varrt.v0     Skewrt.v0     Kurtrt.v0 
# -0.0045986784  0.0093152301 -0.0004226896  0.0002827454  0.0092940822 -0.3285375760   3.1969308454 

# > mean(test.r)
# [1] 0.001545376
# > mean(test.r^2)
# [1] 0.008401503
# > mean(test.r^3)
# [1] -0.0002266442
# > mean(test.r^4)
# [1] 0.0002228258
# > EnvStats::skewness(test.r)
# [1] -0.3450819
# > EnvStats::kurtosis(test.r, excess = F)
# [1] 3.180885


# set.seed(123)
# system.time({test = GAScheme(S0 = 100, v0 = 0.1, kappa = 0.5, theta = 0.19, sigma = 1, mu = 0, rho = -0.7, T = 1/12, n = 1, M = 10000)})

# set.seed(123)
# system.time({test = GAScheme(S0 = 100, v0 = 0.1, kappa = 0.5, theta = 0.19, sigma = 1, mu = 0, rho = -0.7, T = 1/12, n = 1, M = 10000, zeta = 85)})

# test.r = test[,ncol(test)] - log(100)
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
# [1] -0.003891651
# > mean(test.r^2)
# [1] 0.008391665
# > mean(test.r^3)
# [1] -0.0009040571
# > mean(test.r^4)
# [1] 0.0003662229
# > EnvStats::skewness(test.r)
# [1] -1.051753
# > EnvStats::kurtosis(test.r, excess = F)
# [1] 5.031278
