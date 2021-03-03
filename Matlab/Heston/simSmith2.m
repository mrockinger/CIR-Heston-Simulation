function IV = simSmith2( kappa, theta, sigma, D, v1, v2, zgrid, hzv, Nzv, fv, m1Sv, stdISv)
% simulates M draws of the IV with the approach of Smith at once

omega = 0.5;

% get z and verify for this z the moments
z = omega.*(v1+v2)/2 + (1-omega)*sqrt(v1.*v2);


m1BK = m1IVq(kappa, theta, sigma, D, v1, v2); % do not use CF stuff, use GK approach
varIBK = VarIVq(kappa, theta, sigma, D, v1,v2); % do not use CF stuff, use GK approach.
stdIBK = sqrt(varIBK);

% Here, index must be of length M not a single number.
% I. e. find min(abs(z - zgrid)) for each z!
[~,index] = min(abs(bsxfun(@minus, z', zgrid)), [], 2);

m = m1Sv(index);
sd = stdISv(index);

StartVal = random('InverseGaussian', m, sd.^-1);

ii = (StartVal < 0);
StartVal(ii) = 0.01 * m(ii);

IV = -1*ones(length(v1),1);
% U = unifrnd(0, 1, length(v1), 1); % this will not work, since we may need
% more unif r.v. than that

% since we have convergence problems if U is near 1 (IV is getting
% negative) just do it again, since bisection search is unstable as well
for j=1:length(v1)
    
    while IV(j) < 0 
        IV(j) = fzero(@(x) FBKS(x, hzv(index(j)), Nzv(index(j)), fv{index(j)}) - unifrnd(0,1), StartVal(j));
    end
    
    IV(j) = stdIBK(j) * (IV(j) - m1Sv(index(j))) / stdISv(index(j)) + m1BK(j);
    
end


