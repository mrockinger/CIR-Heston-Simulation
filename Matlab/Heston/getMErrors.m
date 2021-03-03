function [MAEmean, MAEvar, MAEskew, MAEkurt, RMSEmean, RMSEvar, RMSEskew, RMSEkurt, MAPEmean, MAPEvar, MAPEskew, MAPEkurt, Metime] = getMErrors(v0, mu, kappa, theta, sigma, rho, n, T, M, NRuns, RetM, setup) 

if strcmp(setup, 'C')

    m1 = Mom1Bates(T, 0, 0, 0, kappa, theta, sigma, rho, mu, v0);
    m2 = Mom2Bates(T, 0, 0, 0, kappa, theta, sigma, rho, mu, v0);
    m3 = Mom3Bates(T, 0, 0, 0, kappa, theta, sigma, rho, mu, v0);
    m4 = Mom4Bates(T, 0, 0, 0, kappa, theta, sigma, rho, mu, v0);

    varR = m2-m1^2;

    skthR = getSkFromMoments(m1, m2, m3);
    kuthR = getKuFromMoments(m1, m2, m3, m4);

    momth.m = m1; momth.v = varR; momth.sk = skthR; momth.ku = kuthR;

    % feature('numcores'); % Edit the number of used cores (NumWorkers) in Cluster profile 'Local' if it is not the number of all logical cores matlab is using
    if feature('numcores') < length(RetM)
        parpool(feature('numcores'));
    else
        parpool(length(RetM));
    end


    % since we are making many steps in each algo, I think it would be better
    % to parallelize the nested for-loop instead of the bigger one (here, it is
    % therefore exactly the other way around compared to the "unconditional"
    % setting)

    for k=1:NRuns

        for j=1:size(n,2)

            parfor i = 1:size(RetM,1)

                Stats(:,i) = TreatReturns(RetM{i}, M, n(j), T, kappa, theta, sigma, mu, rho, v0, momth, setup);

            end

            Errors{j} = Stats;

        end

        Res{k} = Errors;    
        disp(k)

    end

    delete(gcp('nocreate'));
    
    for k=1:NRuns 

       for j = 1:length(n)

          auxBsmean(j,:) = Res{k}{j}(1,:);
          auxBsvar(j,:) = Res{k}{j}(2,:);
          auxBsskew(j,:) = Res{k}{j}(3,:);
          auxBskurt(j,:) = Res{k}{j}(4,:);

          auxSBsmean = auxBsmean.^2;
          auxSBsvar  =  auxBsvar.^2;
          auxSBsskew = auxBsskew.^2;
          auxSBskurt = auxBskurt.^2;

          auxAPEsmean(j,:) = Res{k}{j}(5,:);
          auxAPEsvar(j,:) = Res{k}{j}(6,:);
          auxAPEsskew(j,:) = Res{k}{j}(7,:);
          auxAPEskurt(j,:) = Res{k}{j}(8,:);

          auxTimes(j,:) = Res{k}{j}(9,:);

        end    

        Bsmean{k} = auxBsmean;
        Bsvar{k} = auxBsvar;
        Bsskew{k} = auxBsskew;
        Bskurt{k} = auxBskurt;

        SBsmean{k} = auxSBsmean;
        SBsvar{k}  = auxSBsvar;
        SBsskew{k} = auxSBsskew;
        SBskurt{k} = auxSBskurt;

        APEsmean{k} = auxAPEsmean;
        APEsvar{k} = auxAPEsvar;
        APEsskew{k} = auxAPEsskew;
        APEskurt{k} = auxAPEskurt;

        times{k} = auxTimes;

    end

elseif strcmp(setup, 'UC')
    
    parpool(size(RetM,1));
    
     for k=1:NRuns

        for j=1:size(T,2)
            
            m1 = m1UnBates(T(j), 0, 0, 0, kappa, theta, sigma, rho, mu);
            m2 = m2UnBates(T(j), 0, 0, 0, kappa, theta, sigma, rho, mu);
            m3 = m3UnBates(T(j), 0, 0, 0, kappa, theta, sigma, rho, mu);
            m4 = m4UnBates(T(j), 0, 0, 0, kappa, theta, sigma, rho, mu);

            varR = m2-m1^2;

            skthR = getSkFromMoments(m1, m2, m3);
            kuthR = getKuFromMoments(m1, m2, m3, m4);

            momth.m = m1; momth.v = varR; momth.sk = skthR; momth.ku = kuthR;

            parfor i = 1:size(RetM,1)

                    % To set setup = 'C' here is not a mistake, since we are interested in the distribution of v_T again
                    Stats(:,i) = TreatReturns(RetM{i}, M, n, T(j), kappa, theta, sigma, mu, rho, v0, momth, 'C');

            end

            Errors{j} = Stats;

        end

        Res{k} = Errors;    
        disp(k)

    end

    delete(gcp('nocreate'));
    
    for k=1:NRuns 

       for j = 1:length(T)

          auxBsmean(j,:) = Res{k}{j}(1,:);
          auxBsvar(j,:) = Res{k}{j}(2,:);
          auxBsskew(j,:) = Res{k}{j}(3,:);
          auxBskurt(j,:) = Res{k}{j}(4,:);

          auxSBsmean = auxBsmean.^2;
          auxSBsvar  =  auxBsvar.^2;
          auxSBsskew = auxBsskew.^2;
          auxSBskurt = auxBskurt.^2;

          auxAPEsmean(j,:) = Res{k}{j}(5,:);
          auxAPEsvar(j,:) = Res{k}{j}(6,:);
          auxAPEsskew(j,:) = Res{k}{j}(7,:);
          auxAPEskurt(j,:) = Res{k}{j}(8,:);

          auxTimes(j,:) = Res{k}{j}(9,:);

        end    

        Bsmean{k} = auxBsmean;
        Bsvar{k} = auxBsvar;
        Bsskew{k} = auxBsskew;
        Bskurt{k} = auxBskurt;

        SBsmean{k} = auxSBsmean;
        SBsvar{k}  = auxSBsvar;
        SBsskew{k} = auxSBsskew;
        SBskurt{k} = auxSBskurt;

        APEsmean{k} = auxAPEsmean;
        APEsvar{k} = auxAPEsvar;
        APEsskew{k} = auxAPEsskew;
        APEskurt{k} = auxAPEskurt;

        times{k} = auxTimes;

    end
    
    
end 


MAEmean = sum(abs(cat(3,Bsmean{:})),3) / NRuns;
MAEvar = sum(abs(cat(3,Bsvar{:})),3) / NRuns;
MAEskew = sum(abs(cat(3,Bsskew{:})),3) / NRuns;
MAEkurt = sum(abs(cat(3,Bskurt{:})),3) / NRuns;

RMSEmean = sqrt(sum(cat(3,SBsmean{:}),3) / NRuns);
RMSEvar = sqrt(sum(cat(3,SBsvar{:}),3) / NRuns);
RMSEskew = sqrt(sum(cat(3,SBsskew{:}),3) / NRuns);
RMSEkurt = sqrt(sum(cat(3,SBskurt{:}),3) / NRuns);

MAPEmean = sum(cat(3,APEsmean{:}),3) / NRuns;
MAPEvar = sum(cat(3,APEsvar{:}),3) / NRuns;
MAPEskew = sum(cat(3,APEsskew{:}),3) / NRuns;
MAPEkurt = sum(cat(3,APEskurt{:}),3) / NRuns;

Metime =  sum(cat(3,times{:}),3) / NRuns;

