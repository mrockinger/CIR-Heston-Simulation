function vtD = fSUB(M,n,T,kappa,theta,sigma,v0)
% Johnson's Subordinated model

D = T/n;

vtD= zeros(n+1,M);
vt = v0.*ones(1,M); % will play role of vtilde.
vtD(1,:) = vt;

emk = exp(-kappa*D);
nu = 4*kappa*theta/(sigma^2); %degree of freedom
c  = sigma^2*(1-emk)/(4*kappa); % Glasserman notation, Book p 124.

if nu>1
   
    Zv = randn(n,M);
    X  = chi2rnd(nu-1, n, M); 
    
    for t=1:n
        L = vtD(t,:).*emk/c; %Lambda
        vtD(t+1,:) = c.*( (Zv(t,:) + sqrt(L)).^2 + X(t,:));
    end
    
else
    
    for t=1:n
        
        L = vtD(t,:).*emk/c; %Lambda
        N = poissrnd(L/2);
        vtD(t+1,:) = c.*chi2rnd(nu+2*N);
        
    end    
end


