function VGK = VarIVq(kappa, theta, sigma, D, Vs, Vt)
% this is the quick implementation of the Var[Ic] based on the Glasserman-Kim decomposition


d  = 4*kappa*theta/(sigma^2);
nu = d/2-1;

kd2 = kappa * D / 2;
C1  = coth( kd2 );
C2  = csch( kd2 )^2;
CZ  = 2 * kappa / ( sigma^2 * sinh(kd2) );

z  = CZ * sqrt( Vs .* Vt );

I1 = find(z < 50); % Calculate indizes to avoid if else statement
I2 = find(z > 50);

Eeta(I1)  = 0.5 * z(I1) .* besseli(nu+1, z(I1)) ./ besseli(nu,z(I1));
Eeta(I2)  = 0.5 * z(I2) .* besseli(nu+1, z(I2), 1) ./ besseli(nu,z(I2),1);    

Eeta2(I1) = 0.25 * (z(I1).^2) .* besseli(nu+2, z(I1)) ./ besseli(nu,z(I1))  + Eeta(I1);   
Eeta2(I2) = 0.25 * (z(I2).^2) .* besseli(nu+2, z(I2), 1) ./ besseli(nu,z(I2),1)  + Eeta(I2);   
   
sig2X1 = ( Vs + Vt ) .* ( sigma^2 * C1 / (kappa^3) + 0.5 * sigma^2 * D * C2 / (kappa^2) ...
                          -0.5 * sigma^2 * D^2 * C1 * C2 / kappa);

sig2X2 = d * sigma^4 * (-8 + 2 * kappa * D * C1 + kappa^2 * D^2 * C2) / ( 8 * kappa^4 );

sig2Z  = 4 * sig2X2 / d;

EX2    = 0.25 * d * sigma^2 * ( -2 + kappa * D * C1 ) / ( kappa^2 );

EZ     = 4 * EX2 / d;

% now this is the variance made by Glasserman-Kim.  
VGK = sig2X1 + sig2X2 + Eeta.*sig2Z + (Eeta2-Eeta.^2).*EZ.^2;


