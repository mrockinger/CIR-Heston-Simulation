MomErrorReturns = function(S0, v0, kappa, theta, sigma, mu, rho, T, n, M){
  
  # Plug in some dummy for v0 here
  ErT = MomentsBates(mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, lambda = 0, muj = 0, vj = 0, t = T, v0 = 1, conditional = FALSE)[c(1, 5:7)]
  varT = ErT[2] - ErT[1]^2

  EulerFT = EulerScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M)[,n+1] - log(100)
  Zhu     = ZhuScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M)[,n+1] - log(100)
  KJ      = IKJIMMScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M)[,n+1] - log(100)
  QE      = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, gamma1 = 0.5, gamma2 = 0.5, optimized = F)[,n+1] - log(100)
  BK      = BKScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M)[,n+1] - log(100)
  S       = Smith(S0, v0, kappa, theta, sigma, mu, rho, T, n, M)[,n+1] - log(100)
  GK      = GEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M)[,n+1] - log(100)
  BBG     = GAScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, zeta = 85)[,n+1] - log(100)
  ORS     = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n, M, optimized = TRUE, target = 'skew', mode = 'unconditional')[,n+1] - log(100)
    
  res = cbind(EulerFT, Zhu, KJ, QE, BK, S, GK, BBG, ORS)
  
  RelAbsDiff = matrix(NA, nrow = length(ErT), ncol = ncol(res))
  
  RelAbsDiff[1,] = abs( (apply(res, 2, mean) - ErT[1]) / ErT[1])
  RelAbsDiff[2,] = abs( (apply(res, 2, var) - varT) / varT)
  RelAbsDiff[3,] = abs( (apply(res, 2, EnvStats::skewness) - ErT[5]) / ErT[5])
  RelAbsDiff[4,] = abs( (apply(res, 2, EnvStats::kurtosis, excess = F) - ErT[6]) / ErT[6])

  return(RelAbsDiff)

}  