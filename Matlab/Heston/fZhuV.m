function [vtD,Zv] = fZhuV(M,n,T,kappa,theta,sigma, v0)
% Zhu's simulation of volatility instead of variance

D = T/n;

Zv = randn(n,M);
xtD= zeros(n+1,M); % placeholder for volat this time
xt = sqrt(v0).*ones(1,M); % now as volatility and not variance.
xtD(1,:) = xt;


% simplifying notations
emk  = exp(-kappa*D);
emhk = exp(-0.5*kappa*D); % half a kappa

% do the simulations
for t=1:n
    
    m1 = theta + (xtD(t,:).^2-theta)*emk;
    m2 = sigma^2*(1-emk)/(4*kappa);
    beta = sqrt( max([m1-m2;zeros(1,M)]) );
    thetax = ( beta-xtD(t,:).*emhk )./( 1-emhk );
    xtD(t+1,:) = xtD(t,:) + 0.5*kappa*( thetax-xtD(t,:) )*D + 0.5*sigma*sqrt(D)*Zv(t,:);
    
end

vtD = xtD.^2; % go from volat to variance
