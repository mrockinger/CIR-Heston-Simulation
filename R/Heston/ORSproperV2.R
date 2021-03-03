source("getGammaOpt.r")
source("DiscretizationConditionalReturnMomentsV2.r")

# kappa = 0.5; sigma = 1; theta = 0.19; T = 4; n = 4; v0 = 0.05; M = 1e6; rho = -0.7; mu = 0; S0 = 100; gamma1 = 0.5; gamma2 = 0.5

ORS = function(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, target = F, mode = F, V = F, v = F, Z_v = F){
  
  t = Delta = T / n
  
  emkt = exp(- kappa * Delta)
  nu   = 4 * theta * kappa / sigma^2
  h    = 4 * kappa / (sigma^2 * (1 - emkt))
  gi = sigma^2 / (2 * kappa * theta)
  
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
    
    K0  = -rho * kappa * theta * Delta / sigma; # Include mu here!
    
    for(i in 1:n){

      # calculate optimized gammas and corresponding Ks
      
      gamma1 = get.opt.gamma1(Delta = i * Delta, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, vt = v0, target = target, mode = mode)
      gamma2 = 1 - gamma1
      
      K1  = gamma1 * Delta * ( kappa * rho / sigma-0.5 ) - rho / sigma;
      K2  = gamma2 * Delta * ( kappa * rho / sigma-0.5 ) + rho / sigma;
      K3  = gamma1 * Delta * (1 - rho * rho);
      K4  = gamma2 * Delta * (1 - rho * rho);
      
      # calculate theoretical process mean and variance
      
      mu1   = log(S0) + (2 * i * Delta * kappa * mu + theta - i * Delta * kappa * theta - v0 + (-theta + v0) / exp(i * Delta * kappa))/(2 * kappa)
      zeta1 = (sigma^2 * (theta - 2*v0) + 4 * exp(kappa * i * Delta) * (sigma^2 * theta - 2 * kappa^2 * (-1 + i * Delta * rho * sigma) * (theta - v0) + kappa * sigma * (i * Delta * sigma * (theta - v0) + 2 * rho*(-2 * theta + v0))) + 
                 exp(2 * kappa * i * Delta) * (2 * i * Delta * kappa^3 * (4 * theta) - 8 * kappa^2 * (theta + i * Delta * rho * sigma * theta - v0) + sigma^2 * (-5 * theta + 2 * v0) +
                                                 2 * kappa * sigma * (8 * rho * theta + i * Delta * sigma * theta - 4 * rho * v0))) / (8 * exp( 2 * kappa * i * Delta) * kappa^3)      
      
      
      #### Calculate theoretical discretization mean and variance
      

      # j will be the time index. I. e. if we calculate the value for 1 * Delta, we have on the LHS an R-index of 1 but time index of 0!
      j = i - 1
        
         
      emjkt = exp(- kappa * j * Delta)
      emjpokt = exp(- kappa * (j + 1) * Delta) 
      emjptkt = exp(- kappa * (j + 2) * Delta) 
         
      # Therefore, calculate all necessary CIR moments first
         
      # E[v_{(j) * Delta} | v0]
      EvjDv0 = v0 * emjkt + theta * (1 - emjkt) 
      # E[v_{(j + 1) * Delta} | v0]
      EvjpoDv0 = v0 * emjpokt + theta * (1 - emjpokt) 
      # E[v^2_{j*Delta} | v0]
      Ev2jDv0 = v0^2 * emjkt^2 + (1 + gi) * (theta^2 * (1 - emjkt)^2 + 2 * v0 * theta *(emjkt - emjkt^2))
      # E[v^2_{(j + 1) * Delta} | v0]
      Ev2jpoDv0 = v0^2 * emjpokt^2 + (1 + gi) * (theta^2 * (1 - emjpokt)^2 + 2 * v0 * theta * (emjpokt - emjpokt^2))
      # E[v_{j* Delta} * v_{(j+1) * Delta}|v0]
      EvjDvjpoDv0 = Ev2jDv0 * emkt + EvjDv0 * theta * (1- emkt)
      
      if(i > 1){
        # E[v_{(j + 1) * Delta} | v0]
        EvjptDv0 = v0 * emjptkt + theta * (1 - emjptkt) 
        # E[v^2_{j*Delta} | v0]
        Ev2jDv0 = v0^2 * emjkt^2 + (1 + gi) * (theta^2 * (1 - emjkt)^2 + 2 * v0 * theta *(emjkt - emjkt^2))
        # E[v^2_{(j+1)*Delta} | v0]
        Ev2jpoDv0 = v0^2 * emjpokt^2 + (1 + gi) * (theta^2 * (1 - emjpokt)^2 + 2 * v0 * theta *(emjpokt - emjpokt^2))
        # E[v_{j* Delta} * v_{(j+2) * Delta}|v0]
        EvjDvjptDv0 = Ev2jDv0 * emjpokt + EvjDv0 * theta * (1- emjpokt)
        # E[v_{(j+1) * Delta} * v_{(j+2) * Delta}|v0]
        EvjpoDvjptDv0 = Ev2jpoDv0 * emkt + EvjpoDv0 * theta * (1- emkt)
      }   
         
      # Now, calculate all moments of the discretization (for increments and prices)
         
      mu1.disc.incr = K0 + K1 * EvjDv0 + K2 * EvjpoDv0
      if(i == 1){
        mu1.disc =  log(S0) + mu1.disc.incr
      }else{
        mu1.disc = mu1.disc +  mu1.disc.incr
      }

      # Here, it suffices to calculate m2 for the increment only, since X0 is just a shifting parameter which does not influence the variance
      # mu2.disc.incr =  K0^2 + (2 * K0 * K1 + K3) * EvjDv0 + K1^2 * EvjDv0^2 + (2 * K0 * K2 + K4 + 2 * K1 * K2 * EvjDv0) * EvjpoDv0 + K2^2 * Ev2jpoDv0
      
      mu2.disc.incr =  K0^2 + (2 * K0 * K1 + K3) * EvjDv0 + (2 * K0 * K2 + K4) * EvjpoDv0 + K1^2 * Ev2jDv0 + K2^2 * Ev2jpoDv0 + 2 * K1 * K2 * EvjDvjpoDv0
      
      if(i == 1){
        # Therefore, we must use m1 of the increments only
        zeta1.disc = mu2.disc.incr - mu1.disc.incr^2
      }else{
        # Calculate the covariance of all available returns first (matrix-like object).
        # Therefore, compute the necessary comoments E[R_i*R_j|v_0] \forall i \neq j
        ERjDRjpoDv0 = K0^2 + K0 * K1 * EvjpoDv0 + K0 * K2 * EvjptDv0 + K0 * K1 * v0 +
          K1^2 * v0 * EvjpoDv0  + K1 * K2 * v0 * EvjpoDv0  + K0 * K2 * EvjDv0 +
          K1 * K2 * Ev2jDv0 + K2^2 * EvjDvjpoDv0 
        
        # Calculate Cov[R_i*R_j|v_0] = E[R_i * R_j|v_0] - E[R_i|v_0] * E[R_j|v_0]  \forall i \neq j] 
        CovD2Dv0 = ERjDRjpoDv0 - mu1.disc.incr * mu1.disc.incr.jm1
        
        zeta1.disc.incr = mu2.disc.incr - mu1.disc.incr^2 
        zeta1.disc = zeta1.disc + zeta1.disc.incr + 2 * CovD2Dv0
      }
      
      # Calculate shift and scale parameter
      
      b = sqrt(zeta1 / zeta1.disc)
      a = mu1 - mu1.disc * b
      
      # Calculate "auxiliary" log-price at each time i * Delta
      
      lnSt_tilde[, i + 1] = lnSt_tilde[, i] +  mu * Delta + K0 + K1 * v[, i] + K2*v[, i+1] + sqrt(K3 * v[, i] + K4 * v[, i+1]) * Z_x[, i]     
      
      # do the scaling afterwards. Since the discretization moments are all those of the standard QE scheme. No shifting etc involved
      # lnSt[,i + 1] = a + b * lnSt_tilde
      
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