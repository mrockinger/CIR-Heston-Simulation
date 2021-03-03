function vtD = fABRV(M,n,T,kappa,theta,sigma,v0)
% Anderson Brotherton-Ratcliffe

D = T/n;

Zv = randn(n,M);
vtD= zeros(n+1,M);
vt = v0.*ones(1,M); % will play role of vtilde.
vtD(1,:) = vt;

% in the following, since vtD(t,:) is positive, there is no need to truncate it
% again


% compute the constants
c1  = 0.5*sigma^2*(1-exp(-2*kappa*D))/kappa;
emk = exp(-kappa*D);

% and run through the samples
for t=1:n
    
    Num = c1*vtD(t,:);
    Den = emk*vtD(t,:)+(1-emk)*theta;
    G2  = log( 1  + Num./(Den.^2) );
    
    vt = Den.*exp( -0.5*G2 + sqrt(G2).*Zv(t,:));
    
    vtD(t+1,:) = max([vt; zeros(1,M)]);
    
end

