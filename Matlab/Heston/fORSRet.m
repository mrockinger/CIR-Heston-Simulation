function X = fORSRet(M,n,T,V,kappa,theta,sigma,mu,rho,target,mode)
                   
D = T/n;

% containts the implementation of the ORS simulation 
fprintf('ORS\n')

K0  = -rho * kappa * theta * D / sigma;
K1a = D * ( kappa * rho / sigma - 0.5 );
K1b = rho / sigma;
Kc  = D * (1 - rho * rho);

S0=100;

X=zeros(n+1,M);

X(1,:) = log(S0)*ones(1,M); % in the end, we compute moments for log-returns so the actual level does not matter.

if strcmp(mode, 'conditional')
    
    % do the simulations
    for t = 1:n

        vt = V(t,:);
        
        % This must be vectorized!
        
        % tic
        for k=1:size(vt,2)
            gam1(k) = simORSInit2(D, mu, kappa, theta, sigma, rho, vt(k) , target, mode);
        end
        % toc
        
        % This will take longer for big sizes of vt then the for - loop
        % Both approaches are way slower than in R 
        
        % tic
        % gam11 = simORSInit2(D, mu, kappa, theta, sigma, rho, vt, target, mode);
        % toc
        
        gam2 = 1 - gam1;

        K1  = gam1 .* K1a - K1b;
        K2  = gam2 .* K1a + K1b;
        K3  = gam1 .* Kc;
        K4  = gam2 .* Kc;

        X(t+1,:) = X(t,:) + mu*D + K0 + K1.*V(t,:) + K2.* V(t+1,:) + sqrt( K3.*V(t,:) + K4.*V(t+1,:) ).*randn(1,M);

    end

elseif strcmp(mode, 'unconditional')   
    
    vt = 1; % This is a dummy, since vt does not matter when matching unconditional moments
    gam1 = simORSInit2(D, mu, kappa, theta, sigma, rho, vt , target, mode);
    gam2 = 1 - gam1;
    
    K1  = gam1 * K1a - K1b;
    K2  = gam2 * K1a + K1b;
    K3  = gam1 * Kc;
    K4  = gam2 * Kc;
    
    for t = 1:n
        
        X(t+1,:) = X(t,:) + mu*D + K0 + K1*V(t,:) + K2*V(t+1,:) + sqrt( K3*V(t,:) + K4*V(t+1,:) ).*randn(1,M);

    end   
    
end

