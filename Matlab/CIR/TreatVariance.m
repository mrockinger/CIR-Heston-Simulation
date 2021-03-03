function Res = TreatVariance(Model,M,n,T,kappa,theta,sigma,v0,momth,setup)

    if strcmp(setup, 'C')

        if     strcmp(Model, 'AE0')
            
            tic
            vtD = fAE0V(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(end,:);
            etime = toc;
            
        elseif strcmp(Model, 'DD') 

            tic
            vtD  = fDD(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(end,:);
            etime = toc;

        elseif strcmp(Model, 'E') 

            tic
            vtD  = fEul(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(end,:);
            etime = toc;

        elseif strcmp(Model, 'EA') 

            tic
            vtD  = fEulAbs(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(end,:);
            etime = toc;

        elseif strcmp(Model, 'ER') 

            tic
            vtD  = fEulRef(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(end,:);
            etime = toc;

        elseif strcmp(Model, 'EFT') 
            
            tic
            vtD  = fEulTrunc(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(end,:);
            etime = toc;

        elseif strcmp(Model, 'TVS') 

            tic
            vtD  = fZhu(M,n,T,kappa,theta,sigma, v0);
            vtD = vtD(end,:);
            etime = toc;

        elseif strcmp(Model, 'HM') 

            tic
            vtD  = fHM(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(end,:);
            etime = toc;

        elseif strcmp(Model, 'KJ')  

            tic
            vtD  = fKJV(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(end,:);
            etime = toc;

        elseif strcmp(Model, 'M')

            tic
            vtD  = fM(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(end,:);
            etime = toc;
            
        elseif strcmp(Model, 'G')  

            tic
            vtD = fG(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(end,:);
            etime = toc;

        elseif strcmp(Model, 'ABR')     

            tic
            vtD  = fABRV(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(end,:);
            etime = toc;

        elseif strcmp(Model, 'QE') 

            tic
            vtD  = fQEV(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(end,:);
            etime = toc;
            

        elseif strcmp(Model, 'SAS') 

            tic
            vtD  = fSank1(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(end,:);
            etime = toc;

        elseif strcmp(Model, 'SAN')  

            tic
            vtD  = fSank2(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(end,:);
            etime = toc;

        elseif strcmp(Model, 'BK')
            
            tic
            vtD  = fSUB(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(end,:);
            etime = toc;

        end   

    elseif strcmp(setup, 'UC')
        
         if     strcmp(Model, 'AE0')
            
            tic
            vtD = fAE0V(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(:,1);
            etime = toc;
            
        elseif strcmp(Model, 'DD') 

            tic
            vtD  = fDD(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(:,1);
            etime = toc;

        elseif strcmp(Model, 'E') 

            tic
            vtD  = fEul(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(:,1);
            etime = toc;

        elseif strcmp(Model, 'EA') 

            tic
            vtD  = fEulAbs(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(:,1);
            etime = toc;

        elseif strcmp(Model, 'ER') 

            tic
            vtD  = fEulRef(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(:,1);
            etime = toc;

        elseif strcmp(Model, 'EFT') 
            
            tic
            vtD  = fEulTrunc(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(:,1);
            etime = toc;

        elseif strcmp(Model, 'TVS') 

            tic
            vtD  = fZhu(M,n,T,kappa,theta,sigma, v0);
            vtD = vtD(:,1);
            etime = toc;

        elseif strcmp(Model, 'HM') 

            tic
            vtD  = fHM(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(:,1);
            etime = toc;

        elseif strcmp(Model, 'KJ')  

            tic
            vtD  = fKJV(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(:,1);
            etime = toc;

        elseif strcmp(Model, 'M')

            tic
            vtD  = fM(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(:,1);
            etime = toc;
            
        elseif strcmp(Model, 'G')  

            tic
            vtD  = fG(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(:,1);
            etime = toc;

        elseif strcmp(Model, 'ABR')     

            tic
            vtD  = fABRV(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(:,1);
            etime = toc;

        elseif strcmp(Model, 'QE') 

            tic
            vtD  = fQEV(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(:,1);
            etime = toc;

        elseif strcmp(Model, 'SAS') 

            tic
            vtD  = fSank1(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(:,1);
            etime = toc;

        elseif strcmp(Model, 'SAN')  

            tic
            vtD  = fSank2(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(:,1);
            etime = toc;

        elseif strcmp(Model, 'BK')
            
            tic
            vtD  = fSUB(M,n,T,kappa,theta,sigma,v0);
            vtD = vtD(:,1);
            etime = toc;


        end  
               
    end
    
   Errors = getStats(vtD, momth);
   Res = [Errors etime];
   
end