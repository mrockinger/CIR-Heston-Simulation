function X = fSmithRet(M,n,T,V,kappa,theta,sigma,mu,rho, zgrid, hzv, Nzv, fv, m1Sv, stdISv)
% containts the implementation following the method of Smith
fprintf('Smith\n')

D = T/n;

IV = zeros( n, M ); % container for the integral I_v_t^{v_t+D} v_u du

% for j=1:M
%     for t=1:n
%         % this is the draw from the famous integral
%         IV(t,j) = simSmith( kappa, theta, sigma, D, V(t,j), V(t + 1,j), zgrid, hzv, Nzv, fv, m1Sv, stdISv); 
%     end
% end


for i=1:n
    
    % this is the draw from the famous integral
    IV(i,:) = simSmith2( kappa, theta, sigma, D, V(i,:), V(i + 1,:), zgrid, hzv, Nzv, fv, m1Sv, stdISv); 
  
 end


S0=100;

X=zeros(n + 1,M);

X(1,:) = log(S0)*ones(1,M); % in the end, we compute moments for log-returns so the actual level does not matter.

for i=1:n
    
    X(i+1,:) = X(i,:) + mu * D - 0.5 * IV(i,:) + rho / sigma * (V(i+1,:) - V(i,:) - kappa *theta * D + kappa * IV(i,:)) + sqrt(1 - rho^2) * sqrt(IV(i,:)).* randn(1,M); 
    
end
