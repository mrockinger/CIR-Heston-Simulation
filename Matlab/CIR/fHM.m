function vt = fHM(M,n,T,kappa,theta,sigma,v0)
% insert Higham/Mao

D = T/n;

Zv = randn(n,M);
vt= zeros(n+1,M);
vt(1,:) = v0;


for t=1:n
    vt(t+1,:) = vt(t,:) + kappa.*(theta-vt(t,:)).*D + sigma.*sqrt(abs(vt(t,:)).*D).*Zv(t,:); % vt_Delta
end
