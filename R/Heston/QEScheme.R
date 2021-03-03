# In this file, I try to find out which programming approach deliviers the fastest way to simulate from the full Heston model while optimizing 
# gamma1 and gamma2 at each step by matching the true conditional skewness with the skewness of Andersens discretization.


source("getGammaOpt.r")

QEScheme = function(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, gamma1, gamma2, optimized = F, target = F, mode = F, V = F, v = F, Z_v = F){
  
  if(optimized && (mode == F || target == F)) stop("You have to choose a target and mode if gamma optimization is desired")
  
  psic = 1.5
  
  Delta = T / n
  
  emkt = exp(-kappa*Delta)
  c1  = sigma*sigma*emkt*(1-emkt)/kappa
  c2  = theta*sigma*sigma*((1-emkt)^2)/(2*kappa)
  
  lnSt = matrix(, M, n + 1)
  lnSt[,1] = log(S0) 
  
  
  # If necessary, simulate CIR process
  
  if(!is.matrix(v)){
    
    v = matrix(, M, n + 1)
    v[,1] = v0
    
    U = matrix(runif(M * n), M, n)
  
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
    
  }
  
  # Do the log-price simulation
  
  Z_x = matrix(rnorm(M * n), M, n) # generate log-price noise
  
  if(optimized && mode == "conditional"){
      
    K0  = -rho * kappa * theta * Delta / sigma;
    
    for(i in 1:n){
      
      if(i == 1){
        
        gamma1 = get.opt.gamma1(Delta = Delta, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, vt = v0, target = target, mode = mode)
        
      }else{
        
        gamma1 = get.opt.gamma1(Delta = Delta, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, vt = v[, i], target = target, mode = mode)
        
        
      }
        
      gamma2 = 1 - gamma1
      
      K1  = gamma1 * Delta * ( kappa * rho / sigma-0.5 ) - rho / sigma;
      K2  = gamma2 * Delta * ( kappa * rho / sigma-0.5 ) + rho / sigma;
      K3  = gamma1 * Delta * (1 - rho * rho);
      K4  = gamma2 * Delta * (1 - rho * rho);
        
        
      lnSt[, i+1] = lnSt[, i] +  mu * Delta + K0 + K1 * v[, i] + K2*v[, i+1] + sqrt(K3 * v[, i] + K4 * v[, i+1]) * Z_x[, i]      
        
    }
      
  }else if(optimized && mode == "unconditional"){
      
      
    gamma1 = get.opt.gamma1(Delta = Delta, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, vt = 0, target = target, mode = mode) #vt does not matter here
    gamma2 = 1 - gamma1
    
    K0  = -rho * kappa * theta * Delta / sigma;
    K1  = gamma1 * Delta * ( kappa * rho / sigma-0.5 ) - rho / sigma;
    K2  = gamma2 * Delta * ( kappa * rho / sigma-0.5 ) + rho / sigma;
    K3  = gamma1 * Delta * (1 - rho * rho);
    K4  = gamma2 * Delta * (1 - rho * rho);  
    
    for(i in 1:n){
      
      lnSt[, i+1] = lnSt[, i] +  mu * Delta + K0 + K1 * v[, i] + K2*v[, i+1] + sqrt(K3 * v[, i] + K4 * v[, i+1]) * Z_x[, i]     
      
    }  
      
  }else if(!optimized){
    
    K0  = -rho * kappa * theta * Delta / sigma;
    K1  = gamma1 * Delta * ( kappa * rho / sigma-0.5 ) - rho / sigma;
    K2  = gamma2 * Delta * ( kappa * rho / sigma-0.5 ) + rho / sigma;
    K3  = gamma1 * Delta * (1 - rho * rho);
    K4  = gamma2 * Delta * (1 - rho * rho);

    for(i in 1:n){
        
      lnSt[, i+1] = lnSt[, i] +  mu * Delta + K0 + K1*v[, i] + K2*v[, i+1] + sqrt(K3*v[, i] + K4*v[, i+1]) * Z_x[, i]        
        
    }
      
  }
  
  if(V == TRUE){
    
    return(list(lnSt,v))
    
  }else if(V == FALSE){
    
    return(lnSt)
    
  }
  
}


# Alternative, unvectorized approach to conduct the QE scheme since first ifelse() gives a Warning message without actually generating NaN values (see discussion of this issue below function)


# QEScheme2 = function(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, gamma1, gamma2, optimized = F, target = F, mode = F, V = F, v = F, Z_v = F){
# 
#   if(optimized && (mode == F || target == F)) stop("You have to choose a target and mode if gamma optimization is desired")
# 
#   psic = 1.5
# 
#   Delta = T / n
# 
#   emkt = exp(-kappa*Delta)
#   c1  = sigma*sigma*emkt*(1-emkt)/kappa
#   c2  = theta*sigma*sigma*((1-emkt)^2)/(2*kappa)
# 
#   lnSt = matrix(, M, n + 1)
#   lnSt[,1] = log(S0)
# 
# 
#   # If necessary, simulate CIR process
#   
#   if(!is.matrix(v)){
#     
#     v = matrix(, M, n + 1)
#     v[,1] = v0
#     
#     U = matrix(runif(M * n), M, n)
# 
#     for(j in 1:M){
#   
#       for(i in 1:n){
#   
#         s2  = v[j,i] * c1 + c2
#         m   = v[j,i] * emkt + theta*(1-emkt)
#         psi = s2/(m*m);
#   
#         if(psi <= psic){
#   
#           psii = 1/psi
#           b2   = 2*psii - 1 + sqrt(2 * psii * (2*psii-1) )
#           a    = m/(1+b2)
#           Z_v    = qnorm(U[j,i])
#           v[j, i+1] = a * (sqrt(b2) + Z_v)^2
#   
#         }else{
#   
#           p    = (psi-1)/(psi+1)
#           beta = (1-p)/m
#   
#           if(U[j,i]<p){
#   
#             v[j, i+1] = 0
#   
#           }else{
#   
#             v[j, i+1] = (1/beta) * log( (1-p)/(1-U[j,i]) )
#   
#           }
#   
#         }
#   
#       }
#   
#     }
#   
#   }
# 
#   
#   Z_x = matrix(rnorm(M * n), M, n) # generate log-price noise   
# 
#   if(optimized && mode == "conditional"){
# 
#     for(i in 1:n){
# 
#       gamma1 = get.opt.gamma1(Delta = Delta, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, vt = v[, i], target = target, mode = mode)
#       gamma2 = 1 - gamma1
# 
#       K0  = -rho * kappa * theta * Delta / sigma;
#       K1  = gamma1 * Delta * ( kappa * rho / sigma-0.5 ) - rho / sigma;
#       K2  = gamma2 * Delta * ( kappa * rho / sigma-0.5 ) + rho / sigma;
#       K3  = gamma1 * Delta * (1 - rho * rho);
#       K4  = gamma2 * Delta * (1 - rho * rho);
# 
# 
#       lnSt[, i+1] = lnSt[, i] +  mu * Delta + K0 + K1 * v[, i] + K2*v[, i+1] + sqrt(K3 * v[, i] + K4 * v[, i+1]) * Z_x[, i]
# 
#     }
# 
#   }else if(optimized && mode == "unconditional"){
# 
# 
#     gamma1 = get.opt.gamma1(Delta = Delta, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, vt = 0, target = target, mode = mode) #vt does not matter here
#     gamma2 = 1 - gamma1
# 
#     K0  = -rho * kappa * theta * Delta / sigma;
#     K1  = gamma1 * Delta * ( kappa * rho / sigma-0.5 ) - rho / sigma;
#     K2  = gamma2 * Delta * ( kappa * rho / sigma-0.5 ) + rho / sigma;
#     K3  = gamma1 * Delta * (1 - rho * rho);
#     K4  = gamma2 * Delta * (1 - rho * rho);
# 
#     for(i in 1:n){
# 
#       lnSt[, i+1] = lnSt[, i] +  mu * Delta + K0 + K1 * v[, i] + K2*v[, i+1] + sqrt(K3 * v[, i] + K4 * v[, i+1]) * Z_x[, i]
# 
#     }
# 
#   }else if(!optimized){
# 
#     K0  = -rho * kappa * theta * Delta / sigma;
#     K1  = gamma1 * Delta * ( kappa * rho / sigma-0.5 ) - rho / sigma;
#     K2  = gamma2 * Delta * ( kappa * rho / sigma-0.5 ) + rho / sigma;
#     K3  = gamma1 * Delta * (1 - rho * rho);
#     K4  = gamma2 * Delta * (1 - rho * rho);
# 
#     for(i in 1:n){
# 
#       lnSt[, i+1] = lnSt[, i] +  mu * Delta + K0 + K1*v[, i] + K2*v[, i+1] + sqrt(K3*v[, i] + K4*v[, i+1]) * Z_x[, i]
# 
#     }
# 
#   }
# 
#   if(V == TRUE){
# 
#     return(list(lnSt,v))
# 
#   }else if(V == FALSE){
# 
#     return(lnSt)
# 
#   }
# 
# }

# QEScheme give an akward Warning messeage.
# The following will show, that even if the message is occuring, the result is identical to the non-vectorized approach:
# Take as parameters: 
# S0 = 100; kappa = 0.5; theta = 0.19; sigma = 1; mu = 0; rho = -0.7; T = 1; n = 1; M = 1e5; v0 = rgamma(M, shape = 2 * kappa * theta / sigma^2, scale = theta * (2 * kappa * theta / sigma^2)^-1)

# Run:
# set.seed(123)
# v.V = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, gamma1 = 0.5, gamma2 = 0.5, V = T) # will(most likely) produce 2 or more Warning messages
# set.seed(123)
# v.NV = QEScheme2(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, gamma1 = 0.5, gamma2 = 0.5, V = T)
 
# all.equal(v.V, v.NV) # will give TRUE
# sum(is.na(v.V)) # gives 0
# sum(is.na(v.NV)) # gives also 0


# set.seed(123)
# test = QEScheme(100, 0.1, 3, 0.19, 0.4, 0, -0.7, 1/12, 1, 10000, gamma1 = 0.5, gamma2 = 0.5)
# 
# 
# set.seed(123)
# v = QE(0.1, 3, 0.19, 0.4, T = 1/12, n = 1, M = 10000)
# test2 = QEScheme(100, 0.1, 3, 0.19, 0.4, 0, -0.7, 1/12, 1, 10000, gamma1 = 0.5, gamma2 = 0.5, v = v)
# 
# all.equal(test, test2)
