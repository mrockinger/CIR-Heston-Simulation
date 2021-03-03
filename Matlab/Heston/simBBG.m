function x = simBBG(v1, v2, D, kappa ,sigma, C1, C2, grid, EIV_star, VarIV_star)
% here we suppose to be in a simulation from BBG model

EX1    = ( v1 + v2 ) .* ( C1/kappa - D * C2 / 2 );
sig2X1 = ( v1 + v2 ) .* sigma^2 * (C1 / (kappa^3) + 0.5 * D * C2 / (kappa^2) ...
                        -0.5 * D^2 * C1 * C2 / kappa);
            
EIc = interp1(grid, EIV_star, v1 .* v2) + EX1;
ii = isnan(EIc);
EIc(ii) = interp1(grid, EIV_star, v1(ii) .* v2(ii), 'spline','extrap') + EX1(ii);
VIc = interp1(grid, VarIV_star, v1 .* v2) + sig2X1;
VIc(ii) = interp1(grid, VarIV_star, v1(ii) .* v2(ii), 'spline','extrap') + sig2X1(ii);

scale = VIc./EIc;
shape = EIc./scale;

% perform the draw of the Gamma Distribution
x = gamrnd(shape, scale);

