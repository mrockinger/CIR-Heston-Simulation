MomError = function(v0, kappa, theta, sigma, T, n, M, p, q = 0, conditional = FALSE, comoments = FALSE, data = FALSE){
  
  # Plug in some dummy for v0 here, it will not make any difference since the unconditional moments are not functions of v0
  Evp    = MomentsCIR(p = p, q = q, kappa = kappa, theta = theta, sigma = sigma, v0 = 0.1, t = T / n, conditional = conditional, comoments = comoments)
  varT   = Evp[2] - Evp[1]^2
  skewT  = getSkFromMoments(Evp[1], Evp[2], Evp[3])
  kurtT  = getKuFromMoments(Evp[1], Evp[2], Evp[3], Evp[4])
  
  Alfonsi    = AE0(v0, kappa, theta, sigma, T, n, M)[,n+1]
  Deelstra   = DD(v0, kappa, theta, sigma, T, n, M)[,n+1]
  # Euler      = EulerCIR(v0, kappa, theta, sigma, T, n, M, mode = "naive")[,n+1]
  EulerAF    = EulerCIR(v0, kappa, theta, sigma, T, n, M, mode = "absorp")[,n+1]
  EulerRF    = EulerCIR(v0, kappa, theta, sigma, T, n, M, mode = "reflect")[,n+1]
  EulerFT    = EulerCIR(v0, kappa, theta, sigma, T, n, M, mode = "FT")[,n+1]
  Zhu        = TVSScheme(v0, kappa, theta, sigma, T, n, M)[,n+1]
  Higham     = HM(v0, kappa, theta, sigma, T, n, M)[,n+1]
  KJcir      = KJ(v0, kappa, theta, sigma, T, n, M)[,n+1]
  # Milstein   = IM(v0, kappa, theta, sigma, T, n, M)[,n+1]
  Glasserman = G(v0, kappa, theta, sigma, T, n, M)[,n+1]
  ABR        = ABRscheme(v0, kappa, theta, sigma, T, n, M)[,n+1]
  Andersen   = QE(v0, kappa, theta, sigma, T, n, M)[,n+1]
  Sanks      = Sankaran(v0, kappa, theta, sigma, T, n, M, mode = "simple")[,n+1]
  Sankn      = Sankaran(v0, kappa, theta, sigma, T, n, M, mode = "normal")[,n+1]
  # IPZ        = IPZScheme(v0, kappa, theta, sigma, T, n, M)[,2] will not use this, because it is very slow!
  ExAppr     = Exact(v0, kappa, theta, sigma, T, n, M)[,n+1]
  
  # res = cbind(Alfonsi, Deelstra, Euler, EulerAF, EulerRF, EulerFT, Zhu, Higham, KJcir, Milstein, Glasserman, ABR, Andersen, Sanks, Sankn, ExAppr)
  res = cbind(Alfonsi, Deelstra, EulerAF, EulerRF, EulerFT, Zhu, Higham, KJcir, Glasserman, ABR, Andersen, Sanks, Sankn, ExAppr)
  
  RelAbsDiff = matrix(NA, nrow = length(p), ncol = ncol(res))
  
  RelAbsDiff[1,] = abs(apply(res, 2, mean) - Evp[1]) / Evp[1]
  RelAbsDiff[2,] = abs(apply(res, 2, var) - varT) / varT
  RelAbsDiff[3,] = abs(apply(res, 2, EnvStats::skewness) - skewT) / skewT
  RelAbsDiff[4,] = abs(apply(res, 2, EnvStats::kurtosis, excess = F) - kurtT) / kurtT
  
  if(data == TRUE){
    return(list(RelAbsDiff, res))
  }else{
    return(RelAbsDiff)
  } 
  
}  