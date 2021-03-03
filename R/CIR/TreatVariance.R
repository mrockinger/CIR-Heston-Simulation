TreatVariance = function(Model, v0, kappa, theta, sigma, T, n, M, momth, setup){
  
  if(setup == "C"){
    
    if(Model == 'AE0'){
      
      tic = Sys.time()
      vT = AE0(v0, kappa, theta, sigma, T, n, M)[,n+1]
      etime = difftime(Sys.time(), tic, units = "secs")  
      
    }else if(Model == 'DD'){
      
      tic = Sys.time()
      vT = DD(v0, kappa, theta, sigma, T, n, M)[,n+1]
      etime = difftime(Sys.time(), tic, units = "secs")  
      
    }else if(Model == 'Euler'){ 
      
      tic = Sys.time()
      vT = EulerCIR(v0, kappa, theta, sigma, T, n, M, mode = "naive")[,n+1]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'EulerAF'){ 
      
      tic = Sys.time()
      vT = EulerCIR(v0, kappa, theta, sigma, T, n, M, mode = "absorp")[,n+1]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'EulerRF'){ 
      
      tic = Sys.time()
      vT = EulerCIR(v0, kappa, theta, sigma, T, n, M, mode = "reflect")[,n+1]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'EulerFT'){ 
      
      tic = Sys.time()
      vT = EulerCIR(v0, kappa, theta, sigma, T, n, M, mode = "FT")[,n+1]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'TVS'){ 
      
      tic = Sys.time()
      vT = TVSScheme(v0, kappa, theta, sigma, T, n, M)[,n+1]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'HM'){ 
      
      tic = Sys.time()
      vT = HM(v0, kappa, theta, sigma, T, n, M)[,n+1]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'M'){ 
      
      tic = Sys.time()
      vT = IM(v0, kappa, theta, sigma, T, n, M)[,n+1]
      etime = difftime(Sys.time(), tic, units = "secs")  
      
    }else if(Model == 'KJ'){ 
      
      tic = Sys.time()
      vT = KJ(v0, kappa, theta, sigma, T, n, M)[,n+1]
      etime = difftime(Sys.time(), tic, units = "secs")  
      
    }else if(Model == 'G'){ 
      
      tic = Sys.time()
      vT = G(v0, kappa, theta, sigma, T, n, M)[,n+1]
      etime = difftime(Sys.time(), tic, units = "secs")   
      
    }else if(Model == 'ABR'){  
      
      tic = Sys.time()
      vT = ABRscheme(v0, kappa, theta, sigma, T, n, M)[,n+1]
      etime = difftime(Sys.time(), tic, units = "secs")  
      
    }else if(Model == 'QE'){  
      
      tic = Sys.time()
      vT = QE(v0, kappa, theta, sigma, T, n, M)[,n+1]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'SSA'){  
      
      tic = Sys.time()
      vT = Sankaran(v0, kappa, theta, sigma, T, n, M, mode = "simple")[,n+1]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'SNA'){  
      
      tic = Sys.time()
      vT = Sankaran(v0, kappa, theta, sigma, T, n, M, mode = "normal")[,n+1]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'Exact'){ 
      
      tic = Sys.time()
      vT = Exact(v0, kappa, theta, sigma, T, n, M)[,n+1]
      etime = difftime(Sys.time(), tic, units = "secs")  
      
    }  
    
  }else if(setup == "UC"){
    
    if(Model == 'AE0'){
      
      tic = Sys.time()
      vT = AE0(v0, kappa, theta, sigma, T, n, M)[1,]
      etime = difftime(Sys.time(), tic, units = "secs")  
      
    }else if(Model == 'DD'){
      
      tic = Sys.time()
      vT = DD(v0, kappa, theta, sigma, T, n, M)[1,]
      etime = difftime(Sys.time(), tic, units = "secs")  
      
    }else if(Model == 'Euler'){ 
      
      tic = Sys.time()
      vT = EulerCIR(v0, kappa, theta, sigma, T, n, M, mode = "naive")[1,]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'EulerAF'){ 
      
      tic = Sys.time()
      vT = EulerCIR(v0, kappa, theta, sigma, T, n, M, mode = "absorp")[1,]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'EulerRF'){ 
      
      tic = Sys.time()
      vT = EulerCIR(v0, kappa, theta, sigma, T, n, M, mode = "reflect")[1,]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'EulerFT'){ 
      
      tic = Sys.time()
      vT = EulerCIR(v0, kappa, theta, sigma, T, n, M, mode = "FT")[1,]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'TVS'){ 
      
      tic = Sys.time()
      vT = TVSScheme(v0, kappa, theta, sigma, T, n, M)[1,]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'HM'){ 
      
      tic = Sys.time()
      vT = HM(v0, kappa, theta, sigma, T, n, M)[1,]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'M'){ 
      
      tic = Sys.time()
      vT = IM(v0, kappa, theta, sigma, T, n, M)[1,]
      etime = difftime(Sys.time(), tic, units = "secs")  
      
    }else if(Model == 'KJ'){ 
      
      tic = Sys.time()
      vT = KJ(v0, kappa, theta, sigma, T, n, M)[1,]
      etime = difftime(Sys.time(), tic, units = "secs")  
      
    }else if(Model == 'G'){ 
      
      tic = Sys.time()
      vT = G(v0, kappa, theta, sigma, T, n, M)[1,]
      etime = difftime(Sys.time(), tic, units = "secs")   
      
    }else if(Model == 'ABR'){  
      
      tic = Sys.time()
      vT = ABRscheme(v0, kappa, theta, sigma, T, n, M)[1,]
      etime = difftime(Sys.time(), tic, units = "secs")  
      
    }else if(Model == 'QE'){  
      
      tic = Sys.time()
      vT = QE(v0, kappa, theta, sigma, T, n, M)[1,]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'SSA'){  
      
      tic = Sys.time()
      vT = Sankaran(v0, kappa, theta, sigma, T, n, M, mode = "simple")[1,]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'SNA'){  
      
      tic = Sys.time()
      vT = Sankaran(v0, kappa, theta, sigma, T, n, M, mode = "normal")[1,]
      etime = difftime(Sys.time(), tic, units = "secs") 
      
    }else if(Model == 'Exact'){ 
      
      tic = Sys.time()
      vT = Exact(v0, kappa, theta, sigma, T, n, M)[1,]
      etime = difftime(Sys.time(), tic, units = "secs")  
      
    }  
    
  }
  
  
  Bm  = mean(vT) - momth[1]
  Bv  = var(vT) - momth[2]
  Bsk = EnvStats::skewness(vT) - momth[3]
  Bku = EnvStats::kurtosis(vT, excess = F) - momth[4]
  
  
  APEm  = abs( Bm / momth[1] )
  APEv  = abs( Bv / momth[2] )
  APEsk = abs( Bsk / momth[3] )
  APEku = abs( Bku / momth[4] )
  
  return(c(Bm, Bv, Bsk, Bku, APEm, APEv, APEsk, APEku, etime))
  
}
