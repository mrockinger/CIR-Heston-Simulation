function   [C1, C2, grid, EX2,sig2X2, EIV_star, VarIV_star] = simTWInit(n,T,kappa,theta,sigma)
% Initialization for Tse-Wan

D = T/n;

d = 4 * kappa * theta / sigma^2;

kd2 = kappa * D / 2;
C1  = coth( kd2 );
C2  = csch( kd2 )^2;
Cz  = 2 * kappa / ( sigma^2 * sinh(kd2) );

EX2 = d * sigma^2 * (- 2 + kappa * D * C1) /(4 * kappa^2);
sig2X2 = d * sigma^4 * (- 8 + 2 * kappa * D * C1 + kappa^2 * D^2 * C2) / (8 * kappa^4);
  
EZ = 4 * EX2 / d;
VarZ = 4 * sig2X2 / d;

N_u = 2^(15 + ceil(log2(n))) + 1;

% Define the grid for sqrt(Vs * Vt) (where the grid is the same as for the
% IPZ scheme)
v_min = 0.0001;
v_max = 8 * sigma;
grid = linspace(v_min, v_max, N_u);

% Calculate the first two moments of eta for a sparse grid of sqrt(Vs * Vt)
% I do not use getEIVI function to do so
grid_sparse = grid(1,1:4:N_u);
  
Eeta     = Cz.* grid_sparse.* besseli(d / 2, Cz.* grid_sparse , 1) ./ (2 * besseli(d / 2 - 1, Cz.* grid_sparse, 1));
Eeta2    = Cz.^2 .* grid_sparse.^2 .* besseli(d / 2 + 1, Cz.* grid_sparse, 1) ./ (4 * besseli(d / 2 - 1, Cz.* grid_sparse, 1))  + Eeta; 


% Afterwards, to a linear interpolation btw. calculated points to have the
% same length as the original grid
EIV_star = interp1(grid_sparse, EX2  + Eeta .* EZ, grid);
VarIV_star = interp1(grid_sparse, sig2X2 + Eeta .* VarZ + (Eeta2 - Eeta.^2) .* EZ.^2, grid);








