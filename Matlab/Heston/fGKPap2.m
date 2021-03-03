function X = fGKPap2(M, n, T, V, kappa, theta, sigma, mu, rho, nu, Cz, K, lambda_n, gamma_n, i_X2, i_Z, q_i_X2, q_i_Z, mcp, I_X2, I_Z)

fprintf('GKPap\n')

D = T/n;


IV = zeros(n, M);

for i=1:n

   IV(i,:) = simGKPap2(V(i,:), V(i+1,:), nu, Cz, K, lambda_n, gamma_n, i_X2, i_Z, q_i_X2, q_i_Z, mcp, I_X2, I_Z, D, sigma);

end

S0=100;

X=zeros(n+1,M);

X(1,:) = log(S0)*ones(1,M); % in the end, we compute moments for log-returns so the actual level does not matter.

for i=1:n
    
    X(i+1,:) = X(i,:) + mu * D - 0.5 .* IV(i,:) + rho / sigma .* (V(i+1,:) - V(i,:) - kappa *theta * D + kappa .* IV(i,:)) + sqrt(1 - rho^2) * sqrt(IV(i,:)).* randn(1,M); 
     
end

