function X = m1IV(kappa, theta, sigma, Delta, Vs, Vt)
            
Arg = (4 * exp(0.5 * Delta * kappa) * kappa * sqrt(Vs*Vt)) / ((-1 + exp(Delta * kappa)) * sigma^2);
Feller = (2 * kappa * theta)/ sigma^2;
if Arg>50
   B3 = besseli(Feller - 1, Arg, 1);
else
   B3 = besseli(Feller - 1, Arg, 0);    
end

toScale=1;
X = (0.5 * (1 - exp(- kappa * Delta)) * ((2 * (exp(2.5 * Delta * kappa) * (- 2 + Delta * kappa) + exp(1.5 * kappa * Delta) * (2 + Delta * kappa)) * ...
sqrt(Vs * Vt) * (besseli(Feller, Arg , toScale) + besseli(Feller - 2, Arg,  toScale))) /(exp(2 * kappa * Delta) * (- 1 + exp(kappa * Delta))) + ...
(2 * Delta * sigma^2 * B3) / exp(kappa * Delta) - (2 * (-1 + exp(kappa * Delta)) * sigma^2 * B3)/( exp(kappa * Delta) * kappa) + ...
((Delta * (1 - exp(kappa * Delta))^2 * sigma^2 + (-2 + 2 * exp(2 * kappa * Delta)) * Vs + ...
Delta * exp(kappa * Delta) * kappa * (-4 * Vs - 4 * Vt) + (-2 + 2 * exp(2 * kappa * Delta)) * Vt) * B3)/ ...
(exp(kappa * Delta) * (-1 + exp(kappa * Delta))))) / ((-1 + exp(- kappa * Delta))^2 * kappa * B3);

% fprintf('%20.16f\n',Arg)
% fprintf('%20.16f\n',Feller)
% fprintf('%20.16f\n',B3)
% fprintf('%20.16f\n',X)

end

