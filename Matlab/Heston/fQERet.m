function X = fQERet(M,n,T,V,kappa,theta,sigma,mu,rho)        
% containts the implementation of the QE scheme
fprintf('QE\n')

D = T/n;

gam1 = 0.5;
gam2 = 0.5;

K0  = -rho * kappa * theta * D / sigma;
K1a = D * ( kappa * rho / sigma - 0.5 );
K1b = rho / sigma;
K1  = gam1 * K1a - K1b;
K2  = gam2 * K1a + K1b;
Kc  = D * (1 - rho * rho);
K3  = gam1 * Kc;
K4  = gam2 * Kc;

S0=100;

X=zeros(n+1,M);

X(1,:) = log(S0)*ones(1,M); % in the end, we compute moments for log-returns so the actual level does not matter.

% do the simulations
for t = 1:n
    
    X(t+1,:) = X(t,:) + mu*D + K0 + K1*V(t,:) + K2*V(t+1,:) + sqrt( K3*V(t,:) + K4*V(t+1,:) ).*randn(1,M);
    
end


