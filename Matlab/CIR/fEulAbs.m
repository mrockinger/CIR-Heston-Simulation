function vtD=fEulAbs(M,n,T,kappa,theta,sigma,v0)
% insert the Euler discretization with absorption

D = T/n;

Zv = randn(n,M);
vtD= zeros(n+1,M);
vtD(1,:) = v0;

for t=1:n
    vtp = max([vtD(t,:); zeros(1,M)]);
    vtD(t+1,:) = vtp + kappa.*(theta-vtp).*D + sigma.*sqrt(vtp.*D).*Zv(t,:); 
end

