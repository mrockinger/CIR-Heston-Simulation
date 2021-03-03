function X = m1IVq(kappa, theta, sigma, D, Vs, Vt)
% this is the quick implementation based on the Glasserman-Kim decomposition

d  = 4*kappa*theta/(sigma^2);
nu = d/2-1;

kd2 = kappa * D / 2;
C1  = coth( kd2 );
C2  = csch( kd2 )^2;
CZ  = 2 * kappa / ( sigma^2 * sinh(kd2) );

z  = CZ*sqrt(Vs.*Vt);

I1 = find(z < 50); % Calculate indizes to avoid if else statement
I2 = find(z > 50);

Eeta(I1) = 0.5 * z(I1) .* besseli(nu+1, z(I1)) ./ besseli(nu, z(I1));
Eeta(I2)  = 0.5 * z(I2) .* besseli(nu+1, z(I2), 1) ./ besseli(nu, z(I2), 1);    

EX2    = 0.25 * d * sigma^2 * ( -2 + kappa * D * C1 ) / ( kappa^2 );

EZ     = 4 * EX2 / d;

EX1    = ( Vs + Vt ) .* ( C1/kappa - D*C2/2 );

% now this is the first moment
X = EX1 + EX2 + Eeta * EZ;

