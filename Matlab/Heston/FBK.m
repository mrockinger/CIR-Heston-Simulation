function F = FBK(x, v1, v2, kappa, theta, sigma,D, u_epsilon)
% this is the approximation of the IVs cdf as it is done in Broadie/Kaya
% function is NOT VECTORIZED at all

h = 2 * pi ./ (x + u_epsilon);

epsilon = 10^-5;

 N = 0;
  
  while abs(CFBK(h .* N, v1, v2, kappa, theta, sigma, D)) / N > pi * epsilon / 2
    
    N = N + 1;
    
  end 

jgrid = (1:N)';

X1 = sin(h.*jgrid*x)./jgrid;

X2 = CFBK(h.*jgrid, v1, v2, kappa, theta, sigma, D);

X2  = real(X2);
X   = sum( X1 .* X2 );

F = ( h.*x + 2*X )./pi - 2 * epsilon;








