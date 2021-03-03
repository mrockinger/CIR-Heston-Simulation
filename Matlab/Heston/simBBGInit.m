function [C1, C2, grid, EIV_star, VarIV_star] = simBBGInit(M, n, T, kappa, theta, sigma, zeta)
% Initialization for BBG 
% Output will be the grid for vs * vt, EIV* and VarIV* for linear
% interpolation in the simBBG algo


D = T/n;

d = 4 * kappa * theta / sigma^2;

kd2 = kappa * D / 2;
C1  = coth( kd2 );
C2  = csch( kd2 )^2;
Cz  = 2 * kappa / ( sigma^2 * sinh(kd2) );

EX2 = d * sigma^2 * (- 2 + kappa * D * C1) /(4 * kappa^2);
VarX2 = d * sigma^4 * (- 8 + 2 * kappa * D * C1 + kappa^2 * D^2 * C2) / (8 * kappa^4);
  
EZ = 4 * EX2 / d;
VarZ = 4 * VarX2 / d;

delta =  exp(log(10 * sigma) / ((75 - zeta) + log2(M)));

xi = log(10 * sigma) / log(delta);

% grid is for vs * vt not for sqrt(vs * vt) as in TW
grid =  [0, delta.^(-zeta + (0:zeta)), delta.^(xi - ((xi-1):-1:0))]; 

% The first value of these two vectors will be NaN since dividing through
% 0, but we know from TW that the moments are 0 if vs*vt = 0
Eeta     = Cz.* sqrt(grid).* besseli(d / 2, Cz.* sqrt(grid) , 1) ./ (2 * besseli(d / 2 - 1, Cz.* sqrt(grid), 1));
Eeta2    = Cz.^2 .* grid .* besseli(d / 2 + 1, Cz.* sqrt(grid), 1) ./ (4 * besseli(d / 2 - 1, Cz.* sqrt(grid), 1))  + Eeta; 
  
Eeta(1,1) = 0;
Eeta2(1,1) = 0;

EIV_star = EX2  + Eeta .* EZ; 
VarIV_star =  VarX2 + Eeta .* VarZ + (Eeta2 - Eeta.^2) .* EZ.^2;

%linear interpolation rule. I have to use this within the function, not
%here. Input will be only the "real" value of vs * vt
%vq = interp1(grid, EIV_star, grid(end,end)-0.1); 




