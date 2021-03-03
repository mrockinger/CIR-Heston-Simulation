function vtD = fKJV(M,n,T,kappa,theta,sigma,v0)
% Kahl-Jäckel

D = T/n;

Zv = randn(n,M);
vtD= zeros(n + 1,M);
vt = v0.*ones(1,M); % will play role of vtilde.
vtD(1,:) = vt;

% in the following, since vtD(t,:) is positive, there is no need to truncate it
% again


% compute the constants
c1 = 1+kappa*D;

% and run through the samples
for t=1:n
    
    vt = (vtD(t,:)+kappa*theta*D+sigma*sqrt(vtD(t,:)).*Zv(t,:)*sqrt(D)+0.25*D*sigma^2*(Zv(t,:).^2-1))/c1;
    vtD(t+1,:) = max([vt; zeros(1,M)]);
    
end

