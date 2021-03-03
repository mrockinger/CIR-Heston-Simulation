function vtD=fSank2(M,n,T,kappa,theta,sigma,v0)
% Sankaran approximation 2

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
    rt = nu+Lt;
    k1 = 1-(nu+2)./(6*rt)-(nu^2-2*nu+10)./(72*rt.^2)-(nu^3-12*nu^2-6*nu+44)./(432*rt.^3) ...
                 -5*(nu^4-28*nu^3+24*nu^2+1112*nu-1028)./(10368*rt.^4);
    k2 = (1-(nu-1)./(6*rt)-(nu^2+nu-2)./(18*rt.^2)-(4*nu^3-9*nu^2-228*nu+223)./(216*rt.^3))./rt;

    vtD(t+1,:) = ( (nu-1)/3 + rt .* ( k1 + sqrt(k2).*Zv(t,:) ).^2 )/h;
    
end

