TreatReturns = function(Model, S0, v0, kappa, theta, sigma, mu, rho, T, n, M, momth, setup){
  
  if(setup == 'C'){
    
    if(Model == 'EulerFT'){ 
      
      tic = Sys.time()
      rT = EulerScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, mode = "FT")[,n+1] - log(S0)
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'KJ'){ 
      
      tic = Sys.time()
      rT = IKJIMMScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M)[,n+1] - log(S0)
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
     }else if(Model == 'QE'){  
      
      tic = Sys.time()
      rT = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, gamma1  = 0.5, gamma2 = 0.5)[,n+1] - log(S0)
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'BK'){ 
      
      tic = Sys.time()
      rT = BKScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M)[,n+1] - log(S0)
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'Smith' && T / n >= 1 / 12){ 
    
      tic = Sys.time()
      rT = Smith(S0, v0, kappa, theta, sigma, mu, rho, T, n, M)[,n+1] - log(S0)
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'GK'){ 
      
      tic = Sys.time()
      rT = GEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M)[,n+1] - log(S0)
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'TW'){ 
      
      tic = Sys.time()
      v  = Exact(v0, kappa, theta, sigma, T, n, M) #Use the exact method instead of the IPZ - scheme since this approach will take forever
      rT = IGScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, v = v)[,n+1] - log(S0)
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'BBG'){ 
      
      tic = Sys.time()
      rT = GAScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, zeta = 80)[,n+1] - log(S0)
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'ORS'){ 
      
      tic = Sys.time()
      rT = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, optimized = TRUE, target = 'skew', mode = 'unconditional')[,n+1] - log(S0)  
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'ORS + SUB'){ 
      
      tic = Sys.time()
      v  = Exact(v0, kappa, theta, sigma, T, n, M)
      rT = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, optimized = TRUE, target = 'skew', mode = 'conditional', v = v)[,n+1] - log(S0)  
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'Zhu'){
      
      tic = Sys.time()
      rT = ZhuScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M)[,n+1] - log(S0)
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else{
      
      rT = rep(-99, M)
      
    }
  
  }else if(setup == 'UC'){
    
    if(Model == 'EulerFT'){ 
      
      tic = Sys.time()
      rT = diff(EulerScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, mode = "FT")[1,])
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'KJ'){ 
      
      tic = Sys.time()
      rT = diff(IKJIMMScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M)[1,])
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'QE'){  
      
      tic = Sys.time()
      rT = diff(QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, gamma1  = 0.5, gamma2 = 0.5)[1,])
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'BK'){ 
      
      tic = Sys.time()
      rT = diff(BKScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M)[1,])
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'Smith' && T / n >= 1 / 12){ 
      
      tic = Sys.time()
      rT = diff(Smith(S0, v0, kappa, theta, sigma, mu, rho, T, n, M)[1,])
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'GK'){ 
      
      tic = Sys.time()
      rT = diff(GEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M)[1,])
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'TW'){ 
      
      tic = Sys.time()
      v  = Exact(v0, kappa, theta, sigma, T, n, M) #Use the exact method instead of the IPZ - scheme since this approach will take forever
      rT = diff(IGScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, v = v)[1,])
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'BBG'){ 
      
      tic = Sys.time()
      rT = diff(GAScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, zeta = 80)[1,])
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'ORS'){ 
      
      tic = Sys.time()
      rT = diff(QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, optimized = TRUE, target = 'skew', mode = 'conditional')[1,])  
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'ORS + SUB'){ 
      
      tic = Sys.time()
      v  = Exact(v0, kappa, theta, sigma, T, n, M)
      rT = diff(QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, optimized = TRUE, target = 'skew', mode = 'conditional', v = v)[1,]) 
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else if(Model == 'Zhu'){
      
      tic = Sys.time()
      rT = diff(ZhuScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M)[1,])
      toc = Sys.time()
      etime = as.numeric(difftime(toc, tic, units = "secs"))
      
    }else{
      
      rT = rep(-99, T)
      
    }
   
  
 }
  
  Bm  = mean(rT) - momth[1]
  Bv  = var(rT) - momth[2]
  Bsk = EnvStats::skewness(rT) - momth[3]
  Bku = EnvStats::kurtosis(rT, excess = F) - momth[4]
  
  APEm  = abs( Bm / momth[1] )
  APEv  = abs( Bv / momth[2] )
  APEsk = abs( Bsk / momth[3] )
  APEku = abs( Bku / momth[4] )
  
  return(c(Bm, Bv, Bsk, Bku, APEm, APEv, APEsk, APEku, etime))


}                  
