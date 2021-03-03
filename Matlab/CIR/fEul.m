function vtD=fEul(M,n,T,kappa,theta,sigma,v0)
% insert the Euler discretization
D = T/n;

Zv = randn(n,M);
vtD= zeros(n+1,M);
vtD(1,:) = v0;

for t=1:n
    vtD = vtD(t,:) + kappa*(theta-vtD(t,:))*D + sigma*sqrt(vtD(t,:).*D).*Zv(t,:);
end
