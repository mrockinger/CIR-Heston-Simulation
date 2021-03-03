source("getGammaOpt.r")
source("DiscretizationConditionalReturnMomentsV2.r")

ORS = function(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, target = F, mode = F, V = F, v = F, Z_v = F){
  
  t = Delta = T / n
  
  emkt = exp(- kappa * Delta)
  nu   = 4 * theta * kappa / sigma^2
  h    = 4 * kappa / (sigma^2 * (1 - emkt))
  
  lnSt = matrix(, M, n + 1)
  lnSt[,1] = log(S0) 
  
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
  
  Z_x = matrix(rnorm(M * n), M, n)
    
  if(mode == "conditional"){
      
    K0  = -rho * kappa * theta * Delta / sigma;
      
    for(i in 1:n){
        
      if(i == 1){
          
        gamma1 = get.opt.gamma1(Delta = Delta, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, vt = v0, target = target, mode = mode)
        gamma2 = 1 - gamma1
        mu1 =  (2 * Delta * kappa * mu + theta - Delta * kappa * theta - v0 + (-theta + v0) / exp(Delta * kappa))/(2 * kappa) + log(S0)
        mu1.disc =  CondMomQE("mu1", Xt = log(S0), vt = v0, kappa, theta, sigma, rho, mu, Delta, gamma1, gamma2, R = F)
        zeta1 = (sigma^2 * (theta - 2*v[,i]) + 4 * exp(kappa * Delta) * (sigma^2 * theta - 2 * kappa^2 * (-1 + t * rho * sigma) * (theta - v0) + kappa * sigma * (t * sigma * (theta - v0) + 2 * rho*(-2 * theta + v0))) + 
                   exp(2 * kappa * Delta) * (2 * Delta * kappa^3 * (4 * theta) - 8 * kappa^2 * (theta + Delta * rho * sigma * theta - v0) + sigma^2 * (-5 * theta + 2 * v0) +
                                           2 * kappa * sigma * (8 * rho * theta + Delta * sigma * theta - 4 * rho * v0))) / (8 * exp( 2 * kappa * Delta) * kappa^3)  
        # Here, Xt does not matter since the variance of the conditional return X_t+Delta - X_t | X_t is equal to the variance of the conditional log-price at time t+Delta
        zeta1.disc = CondMomQE("var", Xt = log(S0), vt = v0, kappa, theta, sigma, rho, mu, Delta, gamma1, gamma2)
        b = sqrt(zeta1 / zeta1.disc)
        a = mu1 - mu1.disc * b
          
      }else{
          
        gamma1 = get.opt.gamma1(Delta = Delta, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, vt = v[, i], target = target, mode = mode)
        gamma2 = 1 - gamma1
        mu1 =  (2 * Delta * kappa * mu + theta - Delta * kappa * theta - v[,i] + (-theta + v[,i]) / exp(Delta * kappa))/(2 * kappa) + lnSt[,i]
        mu1.disc =  CondMomQE("mu1", Xt = lnSt[,i], vt = v[,i], kappa, theta, sigma, rho, mu, Delta, gamma1, gamma2, R = F)
        zeta1 = (sigma^2 * (theta - 2*v[,i]) + 4 * exp(kappa * Delta) * (sigma^2 * theta - 2 * kappa^2 * (-1 + Delta * rho * sigma) * (theta - v[,i]) + kappa * sigma * (Delta * sigma * (theta - v[,i]) + 2 * rho*(-2 * theta + v[,i]))) + 
                   exp(2 * kappa * Delta) * (2 * Delta * kappa^3 * (4 * theta) - 8 * kappa^2 * (theta + Delta * rho * sigma * theta - v[,i]) + sigma^2 * (-5 * theta + 2 * v[,i]) +
                                           2 * kappa * sigma * (8 * rho * theta + Delta * sigma * theta - 4 * rho * v[,i]))) / (8 * exp( 2 * kappa * Delta) * kappa^3)      
        # Here, Xt does not matter since the variance of the conditional return X_t+Delta - X_t | X_t is equal to the variance of the conditional log-price at time t+Delta
        zeta1.disc = CondMomQE("var", Xt = 0, vt = v[,i], kappa, theta, sigma, rho, mu, Delta, gamma1, gamma2)
        b = sqrt(zeta1 / zeta1.disc)
        a = mu1 - mu1.disc * b
          
      }


        
      K1  = gamma1 * Delta * ( kappa * rho / sigma-0.5 ) - rho / sigma;
      K2  = gamma2 * Delta * ( kappa * rho / sigma-0.5 ) + rho / sigma;
      K3  = gamma1 * Delta * (1 - rho * rho);
      K4  = gamma2 * Delta * (1 - rho * rho);
        
        
      lnSt_tilde = lnSt[, i] +  mu * Delta + K0 + K1 * v[, i] + K2*v[, i+1] + sqrt(K3 * v[, i] + K4 * v[, i+1]) * Z_x[, i]     
        
      lnSt[,i + 1] = a + b * lnSt_tilde
        
    }
      
  }else if(mode == "unconditional"){
      
      
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
    
  }
    
  return(lnSt)  
    
}