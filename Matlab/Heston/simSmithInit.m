function [zgrid, m1Sv, stdISv, hzv, Nzv, fv] = simSmithInit(kappa, theta, sigma, T, n)

D = T / n;
v0 = theta;

% Determine values for the calculation of Vmax
epsilon = 10^-5;
  
emkt = exp(- kappa * D);
d = 4 * kappa * theta / sigma^2;

a = d / 2;
b = 2 * kappa * emkt * v0 / (sigma^2 * (1 - emkt));
h = 1 - 2 / 3 * (a + b) * (a + 3 * b) / (a + 2 * b)^2;
p = 0.5 * (a + 2 * b) / (a + b)^2;
m = (h - 1) * (1 - 3 * h);
z = - norminv(epsilon);
  
Vmax = (a + b) * (z * h * sqrt(2 * p * (1 + m * p)) + 1 - h * p * (1 - h + (2 - h) * m * p / 2))^(1 / h) * sigma^2 * (1 - emkt) / (2 * kappa);
  
% Generate grid for z. "...we found 100 to be a sufficient number of slices in the z dimension."
zgrid = linspace(0.0000001, Vmax, 100); 

m1Sv = zeros(100,1);
% m2Sv = zeros(100,1);
stdISv = zeros(100,1);
Nzv= zeros(100,1);
hzv= zeros(100,1);
fv = {};

% for all the z on zgrid determine the optimal h and Nz and calculate the
% c.f. on the corresponding grid {h * j}, j = 1,..,Nz
for k=1:100
    
    zk = zgrid(k);
    
    % this could be accelerated with the use of Glasserman-Kim formula
    m1S = m1IVS(kappa, theta, sigma, D, zk);
    m2S = m2IVS(kappa, theta, sigma, D, zk);
    
    varIS = m2S - m1S^2; %this happens to get negativ for small Delta? But why will BK work then? Checl in R!
    stdIS = sqrt(varIS);
    
    % determine the optimal h and associated N
    ueps = m1S + 5*stdIS;
    hz   =  pi / ueps; % this is the lower boundary of h since otherwise it would be a function of the argument of the cdf
    
    N = 0;
    
   while abs(CFBKS(hz * N, zk, kappa, theta, sigma, D)) / N > pi * epsilon / 2
    
    N = N + 1;
    
   end 
   
       
   f = CFBKS(hz .* (1:N), zk, kappa, theta, sigma, D); 
   
   m1Sv(k) = m1S;
   % m2Sv(k) = m2S;
   stdISv(k) = stdIS;
   hzv(k) = hz;
   Nzv(k) = N;
   fv{k} = real(f);
    
end
% at this stage for each z in zgrid, hzv and Nzv contain the optimal h and
% N and fv contains the real parts of the corresponding cf evaluations



