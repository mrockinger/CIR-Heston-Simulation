function CF = CFGK2(a, kappa, theta, sigma, Delta)
% this is the Phi^2 CF of Glasserman-Kim

L  = sqrt( kappa^2 + 2 * sigma^2 * a ); 

nu = 4*kappa*theta/(sigma^2);

CF = ( L.*sinh( kappa.*Delta/2 ) ./ (kappa.*sinh( L.*Delta/2 )) ).^(nu / 2);