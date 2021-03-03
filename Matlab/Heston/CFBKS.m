function Phi = CFBKS(a, z, kappa, theta, sigma, Delta)
% returns the characteristic function of Smith


nu  = 4*kappa*theta/(sigma^2);
ga  = sqrt( kappa^2 - 2 * sigma^2 * 1i * a ); % gamma(a)
emk = exp( -kappa * Delta);
emg = exp( -ga * Delta);
opemk = 1 + emk;
omemk = 1 - emk;
opemg = 1 + emg;
omemg = 1 - emg;

X1 = (ga.*exp(-0.5.*(ga-kappa).*Delta).*omemk)./(kappa.*omemg);
X2 = exp( 2 * z./(sigma^2).*(kappa.*opemk./omemk - ga.*opemg./omemg ));
AN = z.*4.*ga.*exp(-0.5.*ga.*Delta)./((sigma^2).*omemg); 
AD = z.*4.*kappa.*exp(-0.5.*kappa.*Delta)./((sigma^2).*omemk); 

if AN>50
    BN=besseli(0.5*nu-1,AN,1);
else
    BN=besseli(0.5*nu-1,AN);
end

if AD>50
    BD=besseli(0.5*nu-1,AD,1);
else
    BD=besseli(0.5*nu-1,AD);
end

X3 = BN/BD; % Besseli=modified Bessel first kinds

Phi = X1 .* X2 .* X3;

% fprintf('%20.16f+%20.16fi\n', real(ga), imag(ga))
% fprintf('%20.16f+%20.16fi\n', real(X1), imag(X1))
% fprintf('%20.16f+%20.16fi\n', real(X2), imag(X2))
% fprintf('%20.16f+%20.16fi\n', real(X3), imag(X3))
% 
% fprintf('%20.16f+%20.16fi\n', real(AN), imag(AN))
