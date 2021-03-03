function X =fEulTruncRet(M,n,T,V,Z,mu,rho)
fprintf('EulerFT\n')

D = T / n;

Z_x = rho * Z + sqrt(1 - rho^2) * randn(n,M);
X=zeros(n + 1,M);

Vh=sqrt(V);
Dsq=sqrt(D);

for t=1:n

        X(t+1,:) = X(t,:) + (mu - 0.5 .* V(t,:)) .* D + Vh(t,:) .* Dsq .* Z_x(t,:);
end