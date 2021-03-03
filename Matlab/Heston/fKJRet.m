% Caution: I modified this function (26.08)
function X = fKJRet(M,n,T,V,Z,sigma,mu,rho)
           
% implements Kahl-Jäkel simulation of log-prices. P. 527, formula 131.
% see Lord, Koekkoek and Van Dijk formula (8)
fprintf('KJ\n')

D = T/n;

Vh=sqrt(V);
Dsq=sqrt(D);

S0=100;

X=zeros(n + 1,M);

X(1,:) = log(S0)*ones(1,M); % in the end, we compute moments for log-returns so the actual level does not matter.
Z_x = rho * Z + sqrt(1 - rho^2) * randn(n,M) ;

for t=1:n
    
    X(t+1,:) = X(t,:) + (mu - 0.25*(V(t+1,:) + V(t,:))) * D ...
                      + rho*Vh(t,:).*Dsq.*Z(t,:) ...                   
                      + 0.5*(Vh(t+1,:) + Vh(t,:)).*Dsq.*(Z_x(t,:)-rho*Z(t,:))  ...
                      + 0.25*sigma*rho.*D.*( Z(t,:).^2 - 1 );
end

