function vtD=fSank1(M,n,T,kappa,theta,sigma,v0)
% Sankaran approximation 1
D = T/n;

nu=4*kappa*theta/(sigma^2);

emk = exp(-kappa*D);
h   = 4*kappa/(sigma^2*(1-emk));
    
    
Zv = randn(n,M);
vtD= zeros(n+1,M);
vt = v0.*ones(1,M); % will play role of vtilde.
vtD(1,:) = vt;
    

for t=1:n
        
     Lt = h*emk.*vtD(t,:); %Lambda 
     vtD(t+1,:) = ( 0.5*(nu-1) + ( sqrt(Lt+0.5*(nu-1)) + Zv(t,:)).^2 )/h;
        
end
