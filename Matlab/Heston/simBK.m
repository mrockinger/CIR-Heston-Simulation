function IV = simBK( kappa, theta, sigma, D, v1, v2 )
% simulates M draws from Broadie-Kaya at once

% we know that m1 and m2 are the correct moments since we have verified it also
% via numerical computation
% We calculate these moments not with the moments from the "CF" way but by
% the approach of Tse/Wan
m1   = m1IVq(kappa, theta, sigma, D, v1, v2);
varI = VarIVq(kappa, theta, sigma, D, v1, v2);
stdI = sqrt(varI);

% determine u_epsilon to calculate the optimal bandwith h and number of steps N
ueps = m1 + 5.*stdI;

% get starting value for optimization from "moment-matched" inverse
% Gaussian distribution
StartVal = random('InverseGaussian', m1, stdI.^-1);

ii = (StartVal < 0);
StartVal(ii) = 0.01 .* m1(ii);

IV = -1*ones(length(v1),1);

% since we have convergence problems if U is near 1 (IV is getting
% negative) just do it until we have a reasonable value of IV (i.e. positiv), since bisection search is unstable as well
for j=1:length(v1)
    
    while IV(j) < 0 
        IV(j) = fzero(@(x) FBK(x, v1(j), v2(j), kappa, theta, sigma,D, ueps(j)) - unifrnd(0,1), StartVal(j));
    end
    
end

