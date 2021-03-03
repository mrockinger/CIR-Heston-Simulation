function vt = fEulTrunc(M,n,T,kappa,theta,sigma,v0)
% Euler discretization with full truncation. 

D = T/n;

Zv = randn(n,M);
vt= zeros(n+1,M);
vt(1,:) = v0; % will play role of vtilde.

for t=1:n
    
    vtp = max([vt(t,:); zeros(1,M)]);
    vt(t+1,:) = vt(t,:) + kappa.*(theta-vtp).*D + sigma.*sqrt(vtp.*D).*Zv(t,:); 
   
end
