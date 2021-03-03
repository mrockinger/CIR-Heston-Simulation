
function y = pdf_eta(n, nu, z)
% pdf of Bessel distribution
% Thereby, I use the log-transform of the pdf stated in Yuan/Kalbfleisch
% I also use exponential scaling within the besseli call to make this more
% stable. I manually devide by the exponential factor afterwards, see
% "exp(-abs(real(z)))^-1"
l = - gammaln(n + 1) + log(z) .* (2 .* n + nu) - log(exp(-abs(real(z))).^-1 .* besseli(nu, z, 1)) -...
    gammaln(n + nu + 1) + log(2) .* (-2 .* n - nu);
y = exp(l);