# In this file, I try to find out which programming approach deliviers the fastest way to simulate from the full Heston model while optimizing 
# gamma1 and gamma2 at each step by matching the true conditional skewness with the skewness of Andersens discretization.


source("getGammaOpt.r")

QESchemeV2 = function(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, gamma1, gamma2, optimized = F, target = F, mode = F, V = F, v = F, Z_v = F){
  
  if(optimized && (mode == F || target == F)) stop("You have to choose a target and mode if gamma optimization is desired")
  
  psic = 1.5
  
  t = Delta = T / n
  
  
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
      
      
      lnSt_tilde = lnSt[, i] +  mu * Delta + K0 + K1 * v[, i] + K2*v[, i+1] + sqrt(K3 * v[, i] + K4 * v[, i+1]) * Z_x[, i]     
      
      # das kann alles mit gamma1 oben berechnet werden
      mu1 =  (2 * t * kappa * mu + theta - t * kappa * theta - v[,i] + (-theta + v[,i]) / exp(t * kappa))/(2 * kappa) + lnSt[,i]
      mu1.disc =  CondMomQE("mu1", Xt = lnSt[,i], vt = v[,i], kappa, theta, sigma, rho, mu, Delta, gamma1, gamma2, R = F)
      zeta1 = (sigma^2 * (theta - 2*v[,i]) + 4 * exp(kappa * t) * (sigma^2 * theta - 2 * kappa^2 * (-1 + t * rho * sigma) * (theta - v[,i]) + kappa * sigma * (t * sigma * (theta - v[,i]) + 2 * rho*(-2 * theta + v[,i]))) + 
                 exp(2 * kappa * t) * (2 * t * kappa^3 * (4 * theta) - 8 * kappa^2 * (theta + t * rho * sigma * theta - v[,i]) + sigma^2 * (-5 * theta + 2 * v[,i]) +
                                         2 * kappa * sigma * (8 * rho * theta + t * sigma * theta - 4 * rho * v[,i]))) / (8 * exp( 2 * kappa * t) * kappa^3)                                                                
      zeta1.disc = CondMomQE("var", Xt = 0, vt = v[,i], kappa, theta, sigma, rho, mu, Delta, gamma1, gamma2)
      b = sqrt(zeta1 / zeta1.disc)
      a = mu1 - mu1.disc * b
      
      lnSt[,i + 1] = a + b * lnSt_tilde
      
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

