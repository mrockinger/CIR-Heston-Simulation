function Res = TreatReturns(Model,M,n,T,kappa,theta,sigma,mu,rho,v0,momth,setup)

    if strcmp(setup, 'C')

        if     strcmp(Model, 'EulerFT')
            
            tic
            [V, Z] = fEulTruncV(M,n,T,kappa,theta,sigma,v0);
            X = fEulTruncRet(M,n,T,V,Z,mu,rho);
            R = X(end,:) - X(1,:);
            etime = toc;
            
        elseif strcmp(Model, 'Zhu') 

            tic
            [V, Z] = fZhuV(M,n,T,kappa,theta,sigma,v0);
            X = fZhuRet(M,n,T,V,Z,mu,rho);
            R = X(end,:) - X(1,:);
            etime = toc;

        elseif strcmp(Model, 'KJ') 

            tic
            [V, Z] = fKJV(M,n,T,kappa,theta,sigma,v0);
            X = fKJRet(M,n,T,V,Z,sigma,mu,rho);
            R = X(end,:) - X(1,:);
            etime = toc;

        elseif strcmp(Model, 'QE') 

            tic
            V = fQEV(M,n,T,kappa,theta,sigma,v0);
            X = fQERet(M,n,T,V,kappa,theta,sigma,mu,rho);  
            R = X(end,:) - X(1,:);
            etime = toc;

        elseif strcmp(Model, 'BK') 

            tic
            V = fSUBV(M,n,T,kappa,theta,sigma, v0);
            X = fBKRet(M,n,T,V,kappa,theta,sigma,mu,rho);
            R = X(end,:) - X(1,:);
            etime = toc;

        elseif strcmp(Model, 'Smith') && T / n >= 1 / 12 
            
            tic
            V = fSUBV(M,n,T,kappa,theta,sigma, v0);
            [zgrid, m1Sv, stdISv, hzv, Nzv, fv] = simSmithInit(kappa, theta, sigma, T, n);
            X = fSmithRet(M,n,T,V,kappa,theta,sigma,mu,rho, zgrid, hzv, Nzv, fv, m1Sv, stdISv);
            R = X(end,:) - X(1,:);
            etime = toc;

        elseif strcmp(Model, 'GKPap') 
     
            tic
            V = fSUBV(M,n,T,kappa,theta,sigma, v0);
            [nu, Cz, K, lambda_n, gamma_n, i_X2, i_Z, q_i_X2, q_i_Z, mcp, I_X2, I_Z] = simGKInit(T,n,kappa,theta,sigma);
            X = fGKPap2(M, n, T, V, kappa, theta, sigma, mu, rho, nu, Cz, K, lambda_n, gamma_n, i_X2, i_Z, q_i_X2, q_i_Z, mcp, I_X2, I_Z);
            R = X(end,:) - X(1,:);
            etime = toc;

        elseif strcmp(Model, 'TW') 
            
            tic
            V = fSUBV(M,n,T,kappa,theta,sigma, v0); % Use exact scheme instead of IPZ scheme to draw variance, since IPZ is unreasonable slow
            [C1, C2, grid, EX2,sig2X2, EIV_star, VarIV_star] = simTWInit(n,T,kappa,theta,sigma);
            X = fTWRet(M, n, T, V, kappa, theta, sigma, mu, rho, grid, C1, C2, EIV_star, VarIV_star, EX2, sig2X2);
            R = X(end,:) - X(1,:);
            etime = toc;

        elseif strcmp(Model, 'BBG')  

            tic
            V = fSUBV(M,n,T,kappa,theta,sigma, v0);
            % Choose zeta "wisely". If e. g. M is small (e.g. 1e3), zeta =
            % 85 will result in a non available grid -> choose zeta  = 80
            if 2 * kappa * theta / sigma^2 < 1 
                zeta = 60;
            else
                zeta = 80;
            end
            [C1, C2, grid, EIV_star, VarIV_star] = simBBGInit(M, n, T, kappa, theta, sigma, zeta);
            X = fBBGRet(M,n,T,V,kappa,theta,sigma,mu,rho, C1, C2, grid, EIV_star, VarIV_star);
            R = X(end,:) - X(1,:);
            etime = toc;

        elseif strcmp(Model, 'ORS')
            
            mode = 'unconditional'; % Here, write unconditional since even for the UC setup, we use the 'C' setup mentioned here...
            target = 'skew';

            tic
            V = fQEV(M,n,T,kappa,theta,sigma,v0);
            X = fORSRet(M,n,T,V,kappa,theta,sigma,mu,rho,target,mode);
            R = X(end,:) - X(1,:);
            etime = toc;
            
        elseif strcmp(Model, 'ORS + SUB')

            mode = 'unconditional';  % Here, write unconditional since even for the UC setup, we use the 'C' setup mentioned here...
            target = 'skew';

            tic
            V = fSUBV(M,n,T,kappa,theta,sigma,v0);
            X = fORSRet(M,n,T,V,kappa,theta,sigma,mu,rho,target,mode);
            R = X(end,:) - X(1,:);
            etime = toc;
            
        end   

    elseif strcmp(setup, 'UC')
        
        if     strcmp(Model, 'EulerFT')
            
            tic
            [V, Z] = fEulTruncV(M,n,T,kappa,theta,sigma,v0);
            X = fEulTruncRet(M,n,T,V,Z,sigma,mu,rho);
            R = X(: ,2:end) - X(: , 1:end-1);
            etime = toc;
            
        elseif strcmp(Model, 'Zhu') 

            tic
            [V, Z] = fZhuV(M,n,T,kappa,theta,sigma,v0);
            X = fZhuRet(M,n,T,V,Z,sigma,mu,rho);
            R = X(: ,2:end) - X(: , 1:end-1);
            etime = toc;

        elseif strcmp(Model, 'KJ') 

            tic
            [V, Z] = fKJV(M,n,T,kappa,theta,sigma,v0);
            X = fKJRet(M,n,T,V,Z,sigma,mu,rho);
            R = X(: ,2:end) - X(: , 1:end-1);
            etime = toc;

        elseif strcmp(Model, 'QE') 

            tic
            V = fQEV(M,n,T,kappa,theta,sigma,v0);
            X = fQERet(M,n,T,V,kappa,theta,sigma,mu,rho);  
            R = X(: ,2:end) - X(: , 1:end-1);
            etime = toc;

        elseif strcmp(Model, 'BK') 

            tic
            V = fSUBV(M,n,T,kappa,theta,sigma, v0);
            X = fBKRet(M,n,T,V,kappa,theta,sigma,mu,rho);
            R = X(: ,2:end) - X(: , 1:end-1);
            etime = toc;

        elseif strcmp(Model, 'Smith') && T/n >= 1 / 12 
            
            tic
            V = fSUBV(M,n,T,kappa,theta,sigma, v0);
            [zgrid, m1Sv, stdISv, hzv, Nzv, fv] = simSmithInit(kappa, theta, sigma, T, n);
            X = fSmithRet(M,n,T,V,kappa,theta,sigma,mu,rho, zgrid, hzv, Nzv, fv, m1Sv, stdISv);
            R = X(: ,2:end) - X(: , 1:end-1);
            etime = toc;

        elseif strcmp(Model, 'GKPap') 

            tic
            V = fSUBV(M,n,T,kappa,theta,sigma, v0);
            [nu, Cz, K, lambda_n, gamma_n, i_X2, i_Z, q_i_X2, q_i_Z, mcp, I_X2, I_Z] = simGKInit(T,n,kappa,theta,sigma);
            X = fGKPap2(M, n, T, V, kappa, theta, sigma, mu, rho, nu, Cz, K, lambda_n, gamma_n, i_X2, i_Z, q_i_X2, q_i_Z, mcp, I_X2, I_Z);
            R = X(: ,2:end) - X(: , 1:end-1);
            etime = toc;

        elseif strcmp(Model, 'TW') 

            tic
            V = fSUBV(M,n,T,kappa,theta,sigma, v0); % Use exact scheme instead of IPZ scheme to draw variance, since IPZ is unreasonable slow
            [C1, C2, grid, EX2,sig2X2, EIV_star, VarIV_star] = simTWInit(n,T,kappa,theta,sigma);
            X = fTWRet(M, n, T, V, kappa, theta, sigma, mu, rho, grid, C1, C2, EIV_star, VarIV_star, EX2, sig2X2);
            R = X(: ,2:end) - X(: , 1:end-1);
            etime = toc;
            

        elseif strcmp(Model, 'BBG')  

            tic
            V = fSUBV(M,n,T,kappa,theta,sigma, v0);
            [C1, C2, grid, EIV_star, VarIV_star] = simBBGInit(M, n, T, kappa, theta, sigma, 80);
            X = fBBGRet(M,n,T,V,kappa,theta,sigma,mu,rho, C1, C2, grid, EIV_star, VarIV_star);
            R = X(: ,2:end) - X(: , 1:end-1);
            etime = toc;

        elseif strcmp(Model, 'ORS')
            
            mode = 'unconditional';
            target = 'skew';

            tic
            V = fQEV(M,n,T,kappa,theta,sigma,v0);
            X = fORSRet(M,n,T,V,kappa,theta,sigma,mu,rho,v0,target,mode);
            R = X(: ,2:end) - X(: , 1:end-1);
            etime = toc;
            
        elseif strcmp(Model, 'ORS + SUB')

            mode = 'unconditional';
            target = 'skew';

            tic
            V = fSUBV(M,n,T,kappa,theta,sigma,v0);
            X = fORSRet(M,n,T,V,kappa,theta,sigma,mu,rho,v0,target,mode);
            R = X(: ,2:end) - X(: , 1:end-1);
            etime = toc;
            
        end   
         
    end
    
   Errors = getStats(R, momth);
   Res = [Errors etime];
   
end