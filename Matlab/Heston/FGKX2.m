function FX2 = FGKX2(x, kappa, theta, sigma,D, u_epsilon)
% this is the F_X2 function of Glasserman-Kim
% take care of winding around origin, crossing negative axis

h = 2 * pi ./ (x + u_epsilon);

epsilon = 10^-5;

 N = 0;
  
  while abs(CFGK2(-1i * h * N, kappa, theta, sigma, D)) / N > pi * epsilon / 2
    
    N = N + 1;
    
  end 

jgrid = (1:N)';

X1 = sin(h*jgrid*x)./jgrid;

X2 = CFGK2( -1i*h*jgrid, kappa, theta, sigma, D); % 2, 7, 14 sind anderse als in R: Complex Logarithm?

%X2initial=X2;
%nowinding=true;

% gives same values as R function
for j=1:N-1
    if inquadrant(X2(j+1))==3 && inquadrant(X2(j))==2
        X2(j+1)= CFGK2( -1i*h*jgrid(j+1) + 2*pi, kappa, theta, sigma, D);
        %nowinding=false;        
    end
end

% if ~nowinding
%         figure()
%         plot(X2initial)
%         title('X2 initial')
%         figure()
%         plot(X2)
%          title('X2 after')
%          error('kurr')
% end



X2  = real(X2);
X   = sum( X1 .* X2 );

FX2 = ( h*x + 2*X )/pi;

