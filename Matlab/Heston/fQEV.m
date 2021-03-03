function vtD = fQEV(M,n,T,kappa,theta,sigma,v0)
% Andersen's QE scheme. Instead of simulating U simulate Z~N(0,1). To get U use
% cdfnormal(Z)

D = T/n;

U = unifrnd(0,1,[n M]);
vtD = NaN(n+1,M);
vtD(1,:) = v0.*ones(1,M); 


psic=1.5;

% simplifying notations
emt = exp(-kappa*D);
c1  = sigma*sigma*emt*(1-emt)/kappa;
c2  = theta*sigma*sigma*((1-emt)^2)/(2*kappa);

% do the simulations
for t=1:n
   
    s2 = vtD(t,:) * c1 + c2;
    m  = vtD(t,:) * emt + theta * (1 - emt);
    psi = s2 ./ (m .* m);
    
    % Use logical indexing instead of if else statements in the subsequent
    I = find(psi < psic); 

    vtD(t+1,I) = (m(I)./(1+2.*(1./psi(I)) - 1 + sqrt(2 .* (1./psi(I)) .* (2.*(1./psi(I))-1)))) .* (sqrt(2.*(1./psi(I)) - 1 + sqrt(2 .* (1./psi(I)) .* (2*(1./psi(I))-1))) + norminv(U(t,I))).^2;
    
    I = find(psi > psic);

        I2 = I(U(t,I) < (psi(I) - 1) ./ (psi(I) + 1));

        vtD(t+1,I2) = 0;

        I2 = I(U(t,I) > (psi(I) - 1) ./ (psi(I) + 1));

         vtD(t+1,I2) = (1 ./ ((1 - (psi(I2) - 1) ./ (psi(I2) + 1))  ./ m(I2))) .* log( (1 - (psi(I2) - 1) ./ (psi(I2) + 1)) ./ (1 - U(t, I2)));

        
end

