function vtD=fG(M,n,T,kappa,theta,sigma,v0)
% Glasserman book, P 357.

D = T/n;

Zv = randn(n,M);
vtD= zeros(n+1,M);
vt = v0.*ones(1,M); % will play role of vtilde.
vtD(1,:) = vt;

% in the following, since vtD(t,:) is positive, there is no need to truncate it
% again

% compute the constants
c1 = kappa*theta/4-(sigma^2)/16;
c2 = 1.5*kappa;
% and run through the samples
for t=1:n
    
    sqvt=sqrt(vtD(t,:));
    
    vt = vtD(t,:) + kappa*( theta -  vtD(t,:)).*D + sigma.*sqvt.*sqrt(D).*Zv(t,:) ...
        +0.25.*D.*(sigma^2).*( Zv(t,:).^2 - 1 ) ...
        -0.5.*kappa.^2.*(theta-vtD(t,:))*(D^2) + (c1./sqvt-c2.*sqvt).*Zv(t,:)*sigma*(D^1.5);
    
    vtD(t+1,:) = abs(vt);
end
