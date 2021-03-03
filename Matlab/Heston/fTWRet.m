function X = fTWRet(M, n, T, V, kappa, theta, sigma, mu, rho, grid, C1, C2, EIV_star, VarIV_star, EX2, sig2X2)

fprintf('TW\n')

D = T/n;

IV = zeros( n, M );

for j=1:M
    for i=1:n
        
        % this is the draw from the famous integral
        IV(i,j) = simTW(V(i,j), V(i+1,j), D, kappa, sigma, C1 ,C2, grid, EIV_star, VarIV_star, EX2, sig2X2);
              
    end
end


S0=100;

X=zeros(n+1,M);

X(1,:) = log(S0)*ones(1,M); % in the end, we compute moments for log-returns so the actual level does not matter.

for i=1:n
                  
  X(i+1,:) = X(i,:) + mu * D - 0.5 * IV(i,:) + rho / sigma * (V(i+1,:) - V(i,:) - kappa * theta * D + kappa * IV(i,:)) + sqrt(1 - rho^2) * sqrt(IV(i,:)).* randn(1,M); 
     
end

