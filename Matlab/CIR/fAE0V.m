function vtD = fAE0V(M,n,T,kappa,theta,sigma,v0)
% Alfonsi

D = T/n;

Zv = randn(n,M);
vtD= zeros(n+1,M);
vt = v0.*ones(1,M); % will play role of vtilde.
vtD(1,:) = vt;

% in the following, since vtD(t,:) is positive, there is no need to truncate it
% again

% compute the constants
c1 = 1-kappa*D/2;
c2 = sigma*sqrt(D)/(2*c1);
c3 = (kappa*theta - 0.25*sigma^2)*D;

% and run through the samples
for t=1:n
    
    vt = ( c1*sqrt(vtD(t,:)) + c2.*Zv(t,:) ).^2 + c3;
    vtD(t+1,:) = max([vt; zeros(1,M)]);
    
end


