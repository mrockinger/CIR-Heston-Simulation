function x =  skewUCdist(gamma1, kappa, theta, sigma, rho, t)
  
  % zeta2 is defined as the exact skewness of our process
  zeta2 = skewUnBates(t, 0, 0,0, kappa, theta, sigma, rho);
  
  h = 2 * kappa * rho - sigma;
  as = (sigma*t*(-12*(-1 + exp(kappa * t))*(-1 + rho^2)*sigma - 6*h^2*rho*t + exp(kappa * t)*h^3*t^2 + 6*h*(-2 + exp(kappa * t)*(2 + 2*kappa*t - rho*sigma*t)))*theta)/...
    (16 *exp(kappa * t)*kappa^2);
  bs = (3*(-1 + exp(- kappa *t))*sigma*t*(-8*(-1 + rho^2)*(h + sigma) + 4*h*(2*kappa - rho*sigma)*t +h^3*t^2)*theta)/(16.*kappa^2);
  cs = (3*(1 - exp(- kappa *t))*h*sigma*t^2*(8*kappa - 4*rho*(h + sigma) + h^2*t)*theta)/(16.*kappa^2);
  a = (theta*(t^2*h^2+8 *kappa *t+8 *rho^2-4 *rho *sigma *t)-4* rho *theta *exp(- kappa *t)*(t*h+2 *rho))/(8 *kappa);
  b = (t *theta *h*(exp(- kappa* t)-1)*(t*h+4 *rho))/(4 *kappa);
  c = (t^2 *theta *h^2*(1-exp(- kappa *t)))/(4 *kappa);
  
  
 x  = (as + bs * gamma1 + cs * gamma1^2) / (a + b * gamma1 + c * gamma1^2)^(3/2) - zeta2;