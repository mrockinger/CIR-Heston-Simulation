function vt = fDD(M,n,T,kappa,theta,sigma,v0)
% insert the discretization of Deelstra/Delbaen

D = T/n;

Zv = randn(n,M);
vt= zeros(n+1,M);
vt(1,:) = v0;


for t=1:n
    vt(t+1,:) = vt(t,:) + kappa.*(theta-vt(t,:)).*D + sigma.*sqrt(max([vt(t,:); zeros(1,M)]).*D).*Zv(t,:); % vt_Delta
end
