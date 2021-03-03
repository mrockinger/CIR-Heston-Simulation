# names must be a factor of strings
# momth is just needed if setup is C (conditonal) or UC AND approach == 3
# Need function MomError for UC settings

getCIRMErrors = function(v0, kappa,theta, sigma, T, n, M, NRuns, names, nE, nI, nO, setup, approach){
  
  Results = list()
  no.cores = detectCores() - 1
  
  if(setup == 'C'){
    
    Evp    = MomentsCIR(p = 1:4, q = 0, kappa = kappa, theta = theta, sigma = sigma, v0 = v0, t = T, conditional = TRUE, comoments = FALSE)
    varT   = Evp[2] - Evp[1]^2
    skewT  = getSkFromMoments(Evp[1], Evp[2], Evp[3])
    kurtT  = getKuFromMoments(Evp[1], Evp[2], Evp[3], Evp[4])
    
    momth = cbind(Evp[1], varT, skewT, kurtT)
    
    Errors = list()
    
    if(no.cores > length(names)){
      
      cl = makeCluster(length(names))
      
    }else{
      
      cl = makeCluster(no.cores)
      
    }
    
    registerDoParallel(cl)
    clusterExport(cl,ls(globalenv()))
    
    for(k in 1:NRuns){
      
      for(j in 1:length(n)){
        
        Errors[[j]] = foreach(i = 1:length(names), .combine = 'cbind') %dopar% TreatVariance(Model = names[i], v0 = v0, kappa = kappa, theta = theta, sigma = sigma, T = T, n = n[j], M = M, momth = momth, setup = setup)
        
      }
      
      Results[[k]] = Errors
      
      print(k)
      
    }
    
    stopCluster(cl)
    
  }else if(setup == 'UC'){
    
    if(approach == 1){ # fix n to 1, change T/Delta. v0 is a Gamma random variate here 
      
      # Plug in dummy for v0 and T since the unconditional moments are not functions of these two inputs
      Evp    = MomentsCIR(p = 1:4, q = 0, kappa = kappa, theta = theta, sigma = sigma, v0 = 1, t = 1, conditional = FALSE, comoments = FALSE)
      varT   = Evp[2] - Evp[1]^2
      skewT  = getSkFromMoments(Evp[1], Evp[2], Evp[3])
      kurtT  = getKuFromMoments(Evp[1], Evp[2], Evp[3], Evp[4])
      
      momth = cbind(Evp[1], varT, skewT, kurtT)
      
      Errors = list()
      
      
      if(no.cores > length(T)){
        
        cl = makeCluster(length(T))
        
      }else{
        
        cl = makeCluster(no.cores)
        
      }
      
      registerDoParallel(cl)
      clusterExport(cl,ls(globalenv()))
      
      for(k in 1:NRuns){
        
        for(j in 1:length(T)){
          
          # To set setup = 'C' here is not a mistake, since we are interested in the distribution of v_T again
          Errors[[j]] = foreach(i = 1:length(names), .combine = 'cbind') %dopar% TreatVariance(Model = names[i], v0 = v0, kappa = kappa, theta = theta, sigma = sigma, T = T[j], n = n, M = M, momth = momth, setup = 'C')
          
        }
        
        Results[[k]] = Errors
        
        print(k)
        
      }
      
      stopCluster(cl)
    
    }else if(approach == 2){ # Fix n to a high number, change T/Delta
      
      # plug in some dummy for T here, since the unconditional moments are no function of time
      Evp    = MomentsCIR(p = 1:4, q = 0, kappa = kappa, theta = theta, sigma = sigma, v0 = v0, t = 1, conditional = FALSE, comoments = FALSE)
      varT   = Evp[2] - Evp[1]^2
      skewT  = getSkFromMoments(Evp[1], Evp[2], Evp[3])
      kurtT  = getKuFromMoments(Evp[1], Evp[2], Evp[3], Evp[4])
      
      momth = cbind(Evp[1], varT, skewT, kurtT)
      
      Errors = list()
      
      if(no.cores > length(names)){
        
        cl = makeCluster(length(names))
        
      }else{
        
        cl = makeCluster(no.cores)
        
      }
      
      registerDoParallel(cl)
      clusterExport(cl,ls(globalenv()))
      
      for(k in 1:NRuns){
        
        for(j in 1:length(T)){
          
          Errors[[j]] = foreach(i = 1:length(names), .combine = 'cbind', .packages = "EnvStats") %dopar% TreatVariance(Model = names[i], v0 = v0, kappa = kappa, theta = theta, sigma = sigma, T = T[j], n = n, M = M, momth = momth, setup = setup)
          
        }
        
        Results[[k]] = Errors
        print(k)
        
      }
      
      stopCluster(cl)
        
    }else if(approach == 3){ # Fix T, change n/Delta 
      
      Evp    = MomentsCIR(p = 1:4, q = 0, kappa = kappa, theta = theta, sigma = sigma, v0 = v0, t = T, conditional = FALSE, comoments = FALSE)
      varT   = Evp[2] - Evp[1]^2
      skewT  = getSkFromMoments(Evp[1], Evp[2], Evp[3])
      kurtT  = getKuFromMoments(Evp[1], Evp[2], Evp[3], Evp[4])
      
      momth = cbind(Evp[1], varT, skewT, kurtT)
      
      Errors = list()
      
      if(no.cores > length(names)){
        
        cl = makeCluster(length(names))
        
      }else{
        
        cl = makeCluster(no.cores)
        
      }
      
      registerDoParallel(cl)
      clusterExport(cl,ls(globalenv()))
      
      for(k in 1:NRuns){
        
        for(j in 1:length(n)){
          
          Errors[[j]] = foreach(i = 1:length(names), .combine = 'cbind') %dopar% TreatVariance(Model = names[i], v0 = v0, kappa = kappa, theta = theta, sigma = sigma, T = T, n = n[j], M = M, momth = momth, setup = setup)
          
        }
        
        Results[[k]] = Errors
        print(k)
        
      }
      
      stopCluster(cl)
      
    }
    
  }
  
  Bsmean = list() 
  Bsvar  = list()
  Bsskew = list()
  Bskurt = list()
  
  SBsmean = list() 
  SBsvar  = list()
  SBsskew = list()
  SBskurt = list()
  
  APEsmean = list() 
  APEsvar  = list()
  APEsskew = list()
  APEskurt = list()
  
  etimes = list()
  
  
  for(i in 1:NRuns){
    
    Bsmean[[i]] = do.call(rbind, lapply(Results[[i]], `[`, 1,))
    Bsvar[[i]]  = do.call(rbind, lapply(Results[[i]], `[`, 2,))
    Bsskew[[i]] = do.call(rbind, lapply(Results[[i]], `[`, 3,))
    Bskurt[[i]] = do.call(rbind, lapply(Results[[i]], `[`, 4,))
    
    SBsmean[[i]] = Bsmean[[i]]^2
    SBsvar[[i]]  = Bsvar[[i]]^2
    SBsskew[[i]] = Bsskew[[i]]^2
    SBskurt[[i]] = Bskurt[[i]]^2
    
    APEsmean[[i]] = do.call(rbind, lapply(Results[[i]], `[`, 5,))
    APEsvar[[i]]  = do.call(rbind, lapply(Results[[i]], `[`, 6,))
    APEsskew[[i]] = do.call(rbind, lapply(Results[[i]], `[`, 7,))
    APEskurt[[i]] = do.call(rbind, lapply(Results[[i]], `[`, 8,))
    
    etimes[[i]]   = do.call(rbind, lapply(Results[[i]], `[`, 9,))

  } 

  
  
  # Calculate several error measures
  
  
  MAEmean =  Reduce('+', lapply(Bsmean, abs)) / NRuns
  MAEvar =  Reduce('+', lapply(Bsvar, abs)) / NRuns
  MAEskew =  Reduce('+', lapply(Bsskew, abs)) / NRuns
  MAEkurt =  Reduce('+', lapply(Bskurt, abs)) / NRuns
  
  
  RMSEmean = sqrt(Reduce('+', SBsmean) / NRuns) 
  RMSEvar  = sqrt(Reduce('+', SBsvar) / NRuns)
  RMSEskew = sqrt(Reduce('+', SBsskew) / NRuns)
  RMSEkurt = sqrt(Reduce('+', SBskurt) / NRuns)
  
  
  MAPEmean = Reduce('+', APEsmean) / NRuns # Alternatively: MAPEmean = apply(simplify2array(APEsmean), 1:2, mean)
  MAPEvar  = Reduce('+', APEsvar) / NRuns
  MAPEskew = Reduce('+', APEsskew) / NRuns
  MAPEkurt = Reduce('+', APEskurt) / NRuns
  
  # Calculate mean elapsed time
  
  Metime = Reduce('+', etimes) / NRuns
  
  # varAPEsmean = apply(simplify2array(APEsmean), 1:2, var)
  
  # Calculate respective confidence bounds
  
  q05AEsmean = apply(simplify2array(lapply(Bsmean, abs)), 1:2, quantile, 0.05)
  q95AEsmean = apply(simplify2array(lapply(Bsmean, abs)), 1:2, quantile, 0.95)
  
  q05AEsvar = apply(simplify2array(lapply(Bsvar, abs)), 1:2, quantile, 0.05)
  q95AEsvar = apply(simplify2array(lapply(Bsvar, abs)), 1:2, quantile, 0.95)
  
  q05AEsskew = apply(simplify2array(lapply(Bsskew, abs)), 1:2, quantile, 0.05)
  q95AEsskew = apply(simplify2array(lapply(Bsskew, abs)), 1:2, quantile, 0.95)
  
  q05AEskurt = apply(simplify2array(lapply(Bskurt, abs)), 1:2, quantile, 0.05)
  q95AEskurt = apply(simplify2array(lapply(Bskurt, abs)), 1:2, quantile, 0.95)
  
  q05Bsmean = sqrt(apply(simplify2array(SBsmean), 1:2, quantile, 0.05))
  q95Bsmean = sqrt(apply(simplify2array(SBsmean), 1:2, quantile, 0.95))
  
  q05Bsvar = sqrt(apply(simplify2array(SBsvar), 1:2, quantile, 0.05))
  q95Bsvar = sqrt(apply(simplify2array(SBsvar), 1:2, quantile, 0.95))
  
  q05Bsskew = sqrt(apply(simplify2array(SBsskew), 1:2, quantile, 0.05))
  q95Bsskew = sqrt(apply(simplify2array(APEsskew), 1:2, quantile, 0.95))
  
  q05Bskurt = sqrt(apply(simplify2array(SBskurt), 1:2, quantile, 0.05))
  q95Bskurt = sqrt(apply(simplify2array(SBskurt), 1:2, quantile, 0.95))
  
  q05APEsmean = apply(simplify2array(APEsmean), 1:2, quantile, 0.05)
  q95APEsmean = apply(simplify2array(APEsmean), 1:2, quantile, 0.95)
  
  q05APEsvar = apply(simplify2array(APEsvar), 1:2, quantile, 0.05)
  q95APEsvar = apply(simplify2array(APEsvar), 1:2, quantile, 0.95)
  
  q05APEsskew = apply(simplify2array(APEsskew), 1:2, quantile, 0.05)
  q95APEsskew = apply(simplify2array(APEsskew), 1:2, quantile, 0.95)
  
  q05APEskurt = apply(simplify2array(APEskurt), 1:2, quantile, 0.05)
  q95APEskurt = apply(simplify2array(APEskurt), 1:2, quantile, 0.95)
  
  Moments = factor(c("Mean", "Variance", "Skewness", "Kurtosis"), levels = c("Mean", "Variance", "Skewness", "Kurtosis"))
  groups =  factor(c("Explicit", "Implicit", "Others"), levels = c("Explicit", "Implicit", "Others")) 

  Delta = T / n
  
  df.MErrors = data.frame("Scheme" = rep(rep(names, each = length(Delta)), length(Moments)), "Group" = rep(rep(groups, c(nE * length(Delta), nI * length(Delta), nO * length(Delta))), length(Moments)), 
                          "Delta" = rep(rep(Delta, length(names)), length(Moments)), "Moment" = rep(Moments, each = length(names) * length(Delta)), 
                          "MAE" = c(c(MAEmean), c(MAEvar), c(MAEskew), c(MAEkurt)), "q05AE" = c(c(q05AEsmean), c(q05AEsvar), c(q05AEsskew), c(q05AEskurt)),
                          "q95AE" = c(c(q95AEsmean), c(q95AEsvar), c(q95AEsskew), c(q95AEskurt)),
                          "RMSE" = c(c(RMSEmean), c(RMSEvar), c(RMSEskew), c(RMSEkurt)), "q05RSE" = c(c(q05Bsmean), c(q05Bsvar), c(q05Bsskew), c(q05Bskurt)),
                          "q95RSE" = c(c(q95Bsmean), c(q95Bsvar), c(q95Bsskew), c(q95Bskurt)),
                          "MAPE" = c(c(MAPEmean), c(MAPEvar), c(MAPEskew), c(MAPEkurt)), "q05APE" = c(c(q05APEsmean), c(q05APEsvar), c(q05APEsskew), c(q05APEskurt)),
                          "q95APE" = c(c(q95APEsmean), c(q95APEsvar), c(q95APEsskew), c(q95APEskurt)),
                          "Mean Time" = rep(c(Metime), length(Moments)))
  
  return(df.MErrors)
  
}
