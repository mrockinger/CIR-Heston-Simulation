function X=fBKRet(M, n, T, V,kappa,theta,sigma,mu,rho)
% containts the implementation of the Broadie-Kaya simulation
fprintf('BK\n')

D = T/n;

IV = zeros( n, M );

for i=1:n
        IV(i,:) = simBK( kappa, theta, sigma, D, V(i,:), V(i+1,:) ); % this is the draw from the famous integral
end

S0=100;

X=zeros(n+1,M);

X(1,:) = log(S0)*ones(1,M); % in the end, we compute moments for log-returns so the actual level does not matter.

for i=1:n
    
    X(i+1,:)=  X(i,:) + mu * D - 0.5 * IV(i,:) + rho / sigma * (V(i+1,:) - V(i,:) - kappa *theta * D + kappa * IV(i,:)) + sqrt(1 - rho^2) * sqrt(IV(i,:)).* randn(1,M); 
     
end

