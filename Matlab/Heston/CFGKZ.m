function CF = CFGKZ( a, kappa, sigma, Delta)
% this is the Phi^2 CF of Glasserman-Kim

L  = sqrt( kappa^2 + 2 * sigma^2 * a ); 

CF = ( L.*sinh( kappa*Delta/2 ) ./ (kappa.*sinh( L.*Delta/2 )) );

CF = CF.*CF; % the tric to avoir exponential and log.