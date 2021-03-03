function vtD=fEulRef(M,n,T,kappa,theta,sigma,v0)
% insert the Euler discretization with reflection

D = T/n;

Zv = randn(n,M);
vtD= zeros(n+1,M);
vtD(1,:) = v0;


for t=1:n
    
    vtp = abs(vtD(t,:));
    vtD(t+1,:) = vtp + kappa.*(theta-vtp).*D + sigma.*sqrt(vtp.*D).*Zv(t,:); 

end
