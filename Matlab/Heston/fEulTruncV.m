function [vtD,Zv] = fEulTruncV(M,n,T,kappa,theta,sigma,v0)
% Performs the Euler truncation 
% returns volatilities

D = T/n;

Zv = randn(n,M);
vtD= zeros(n+1,M);
vt = v0.*ones(1,M); % will play role of vtilde.
vtD(1,:) = vt;

% in the following, since vtD(t,:) is positive, there is no need to truncate it
% again
for t=1:n
    
    vtp = max([vt; zeros(1,M)]);
    vt = vt + kappa*(theta-vtp)*D + sigma*sqrt(vtp.*D).*Zv(t,:); % vt_Delta
    vtD(t+1,:) = max([vt; zeros(1,M)]);
    
end


