function FZ = FGKZ(x,kappa,sigma, D, u_epsilon) 

h = 2 * pi ./ (x + u_epsilon);

epsilon = 10^-5;

 N = 0;
  
  while abs(CFGKZ(-1i * h * N, kappa, sigma, D)) / N > pi * epsilon / 2
    
    N = N + 1;
    
  end 

jgrid = (1:N)';

X1 = sin(h*jgrid*x)./jgrid;

X2 = CFGKZ( -1i*h*jgrid, kappa, sigma, D); % 2, 7, 14 sind anderse als in R: Complex Logarithm?

X2  = real(X2);
X   = sum( X1 .* X2 );

FZ = ( h*x + 2*X )/pi;