function [nu, Cz, K, lambda_n, gamma_n, i_X2, i_Z, q_i_X2, q_i_Z, mcp, I_X2, I_Z] = simGKInit(T,n,kappa,theta,sigma)

D = T/n;

d  = 4*kappa*theta/(sigma^2);
nu = d / 2 - 1;

kd2 = kappa * D / 2;
C1  = coth( kd2 );
C2  = csch( kd2 )^2;
Cz = (2 * kappa / sigma^2) / sinh(kd2);

% Define stuff for simulating X1 via (truncated) Gamma expansion
  
K = 10; % Truncation level
nX1 = 1:K;

  
lambda_n =  16 * pi^2 .* nX1.^2  ./ (sigma^2 * D * (kappa^2 * D^2 + 4 * pi^2 .* nX1.^2));
gamma_n =  (kappa^2 * D^2 + 4 * pi^2 .* nX1.^2) / (2 * sigma^2 * D^2);
  

% Precalculate a table of cdf values of X2 and Z 
% Thereby, approximate the cdf via Abate/Whitts POISSON algorithm
  
mSX2    = sigma^2*(-2 + kappa*D*C1)/(4*kappa^2);
sig2SX2 = sigma^4*(-8+2*kappa*D*C1+kappa^2*D^2*C2)/(8*kappa^4);

mX2 = d * mSX2;
sig2X2 = d * sig2SX2;

mZ = 4 * mSX2;
sig2Z = 4 * sig2SX2;

u_epsilon_X2 = mX2 + 12 * sqrt(sig2X2); 
u_epsilon_Z = mZ + 12 * sqrt(sig2Z); 

% Define evaluation grid for the cdf approximations of X2 and Z (see p.
% 284)
M = 200;
w = 0.01;

index = 1:(M+1);

X2_CDF_eval_grid = w * mX2 + (index - 1) / M * (u_epsilon_X2 - w * mX2);
Z_CDF_eval_grid = w * mZ + (index - 1) / M * (u_epsilon_Z - w * mZ);

% Do the CDF tabulation for X2 and Z on the given grids
CDF.X2 = zeros(M+1,1); 
CDF.Z = zeros(M+1,1); 

% Not sure how to vectorize the cdf approximations properly, so use for loop
% For really small arguments in CDF.X2, we get not the same values as in R. 
% See Function - file for discussion
for j=1:(M+1) 
    CDF.X2(j) =  FGKX2(X2_CDF_eval_grid(j), kappa, theta, sigma, D, u_epsilon_X2); 
    CDF.Z(j) =  FGKZ(Z_CDF_eval_grid(j), kappa, sigma, D, u_epsilon_Z);
end

  
% Define everything for the "cutpoint method"

% i_X2 and i_Z are arguments for the continous cdfs (not sure if this is
% right since this thishshould be)
% q_i_X2 / q_i_Z are cdf function values on these arguments
% I_X2 and I_Z are quantiles of the respective r.v.s, i. e. we have a
% quantile function approximation

mcp = 100; %number of cutpoints

% For X2
  
i_X2 = linspace(X2_CDF_eval_grid(1), X2_CDF_eval_grid(M+1), 1e4); % I am not sure if this is right because of the choosen n (which is now arbitrarly)
q_i_X2 = interp1(X2_CDF_eval_grid, CDF.X2, i_X2);
 
I_X2 = zeros(mcp, 1);
for j=1:mcp
    
   index = find(q_i_X2 > (j - 1) / mcp);
   I_X2(j) = min(i_X2(index));
    
end
I_X2(mcp + 1) = X2_CDF_eval_grid(M + 1);
  
% For Z
  
i_Z = linspace(Z_CDF_eval_grid(1), Z_CDF_eval_grid(M+1), 1e4);
q_i_Z = interp1(Z_CDF_eval_grid, CDF.Z, i_Z);
  
I_Z = zeros(mcp, 1);
for j=1:mcp
    
    index = find(q_i_Z > (j - 1) / mcp);
    I_Z(j) = min(i_Z(index));
    
end
I_Z(mcp + 1) = Z_CDF_eval_grid(M + 1);
