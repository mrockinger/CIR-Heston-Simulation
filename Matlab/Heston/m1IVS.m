function X = m1IVS(kappa, theta, sigma, Delta, z)
% modification required for Smith 2007.
            
ekt = exp(kappa * Delta);
sekt = sqrt(ekt);
  
Feller = 2 * kappa * theta / sigma^2;
Arg = (4*sekt*kappa*z)/((-1 + ekt)*sigma^2);
  
BM1 = besseli(Feller - 1, Arg);
  
  
X = (2*sekt*kappa*(2 + Delta*kappa + ekt*(-2 + Delta*kappa))*z*besseli(Feller, Arg) + ...
      2*sekt*kappa*(2 + Delta*kappa + ekt*(-2 + Delta*kappa))*z* ...
      besseli(-2 + Feller, Arg) + ...
      ((-1 + ekt)*(2 + Delta*kappa + ekt*(-2 + Delta*kappa))*sigma^2 - 4*kappa*(1 - ekt^2 + 2*Delta*ekt*kappa)*z)* ...
      BM1)/(2.* (-1 + ekt)^2*kappa^2*BM1);
  

end

