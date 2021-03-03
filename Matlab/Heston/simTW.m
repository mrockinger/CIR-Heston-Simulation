function x = simTW( v1, v2, D, kappa, sigma, C1 ,C2, grid, EIV_star, VarIV_star, EX2, sig2X2)
% here we suppose to be in a simulation from Tan and Wan model

EX1    = ( v1 + v2 ) .* ( C1/kappa - D * C2 / 2 );
sig2X1 = ( v1 + v2 ) .* ( sigma^2 * C1 / (kappa^3) + 0.5 * sigma^2 * D * C2 / (kappa^2) ...
                -0.5 * sigma^2 * D^2 * C1 * C2 / kappa);
 
if v1==0 || v2==0
    
    EIc = EX1 + EX2;
    VIc = sig2X1 + sig2X2;
    
else
    
    [~,miidx] = min(abs(sqrt(v1 * v2)-grid));
    EIc = EX1 + EIV_star(miidx);
    VIc = sig2X1 + VarIV_star(miidx);
    
end    


m = EIc;
s = m^3/VIc;

% perform the draw of the Inverse Gaussian Distribution
n = randn(1,1);
u = rand(1,1);
y = 2*s/m;
x = 1 + (n^2)/y - sqrt(2*y*n^2+n^4)/y;
if u>1/(1+x)
    x=m/x;
else
    x=m*x;
end

%fprintf('this is the x from the Tan-Wan simulation %20.16f\n', x)