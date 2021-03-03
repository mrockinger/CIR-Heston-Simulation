function vtD=fM(M,n,T,kappa,theta,sigma,v0)
% insert the Milstein discretization

D = T/n;

Zv = randn(n,M);
vtD= zeros(n+1,M);
vtD(1,:) = v0;

for t=1:n
    vtD(t+1,:) = vtD(t,:) + kappa * (theta - vtD(t,:)) * D + sigma * sqrt(vtD(t,:)) .* Zv(t,:) * sqrt(D) + 0.25 * sigma^2 * ((Zv(t,:) .* sqrt(D)).^2 - D);
end
