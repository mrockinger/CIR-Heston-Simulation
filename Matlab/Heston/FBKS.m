function F = FBKS(x, h, N, fv)
% this is the F function of Smith
jgrid = (1:N);

X1 = sin(h.*jgrid.*x)./jgrid;

X = sum( X1 .* fv );

F = ( h*x + 2*X )/pi;








