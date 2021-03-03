function Phi = CFBK(a, v1, v2, kappa, theta, sigma, D)
% Function is only vectorized in a

d = 4 * kappa * theta / sigma^2;
nu = d / 2 - 1;
gamma = sqrt(kappa^2 - 2 * sigma^2 * a .* 1i);
  
% A is the complex argument of the modified Bessel function -> complex logarithm will be calculated. Therefore, we will keep track of Arg(A) to not have discontinuities, see p. 221 on the right side
% If -pi < Arg(A) < pi does NOT hold, evaluate exp(complex(real = 0, imaginary = m * (d/2 - 1) * pi)) * BesselI(A, nu = d/2 - 1, expon.scaled = T) instead of BesselI(A, nu = d/2 - 1, expon.scaled = T)
A = (4 .* gamma .* exp(- 0.5 .* gamma .* D)) ./ (sigma^2 * (1 - exp(- gamma .* D))) * sqrt(v1 * v2);
B = (4 * kappa * exp(- 0.5 * kappa * D)) ./ (sigma^2 * (1 - exp(- kappa * D))) * sqrt(v1 * v2);
C = ((v1 + v2) / sigma^2) .* ( (kappa * (1 + exp(- kappa * D))) ./ (1 - exp(- kappa * D)) - (gamma .* (1 + exp(- gamma .* D))) ./  (1 - exp(- gamma .* D)) );
  
  
if(B >= 50)

    Phi = (gamma .* exp(- 0.5 .* (gamma - kappa) * D) .* (1 - exp(- kappa * D))) ./ (kappa * (1 - exp( - gamma .* D))) .* besseli(nu, A, 1) ./ besseli(nu, B, 1) .* exp(C);
    
else
    
    Phi = (gamma .* exp(- 0.5 .* (gamma - kappa) * D) .* (1 - exp(- kappa * D))) ./ (kappa * (1 - exp( - gamma .* D))) .* besseli(nu, A) ./ besseli(nu, B) .* exp(C);
  
end  

