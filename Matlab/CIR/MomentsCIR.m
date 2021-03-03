function res = MomentsCIR(p, q, kappa, theta, sigma, v0, t, conditional, comoments)

if conditional && ~comoments
    
    [res, namesres] = Evp_v0_A(p, kappa, theta, sigma, v0, t);
    for k=1:length(res)
        fprintf('%20s %16.10f\n',namesres{k},res(k))
    end
    
elseif conditional && comoments
    
    [res,colnames,rownames] = Evpvq_v0(p, q, kappa, theta, sigma, v0, t);
    for i=1:size(res,1)
        for j=1:size(res,2)
            fprintf('%16.10f  ',res(i,j))
        end
        fprintf('\n')
    end
    
    
elseif ~conditional && ~comoments
    
    [res, namesres] = Evp_A(p, kappa, theta, sigma, v0);
    for k=1:length(res)
        fprintf('%20s %16.10f\n',namesres{k},res(k))
    end
    
    
elseif ~conditional && comoments
    
    [res, rownames, colnames] = Evpvq(p, q, kappa, theta, sigma, v0, t);
    for i=1:size(res,1)
        for j=1:size(res,2)
            fprintf('%16.10f  ',res(i,j))
        end
        fprintf('\n')
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function x = aFun(p, kappa, v0)
x = v0.^(p - 1)./kappa;
end

function x = bFun(kappa, t)
x = exp(kappa .* t) - 1;
end

function x = cFun(p, kappa, theta, sigma)
x = (p - 1) .* (kappa .* theta + 0.5 .* (p - 2) .* sigma.^2);
end

function x = IFun(k, kappa, t)
x = bFun(kappa, t).^k ./ (factorial(k) .* kappa.^(k-1));
% transform this into gama for large k
end

function res = mFun(k, p, kappa, theta, sigma, v0)
res = zeros(length(k),1);
for i = 1:length(k)
    if (k(i) == 0)
        res(i) = 0;
    else
        res(i) = aFun(p - k(i), kappa, v0) * prod( cFun(((p - k(i) + 1):p)', kappa, theta, sigma));
    end
end
end

function x = nuFun (p, kappa, theta, sigma, v0, t)
x = aFun(p, kappa, v0) .* bFun(kappa, t) + sum(mFun( (0:(p - 1))', p, kappa, theta, sigma, v0) .* IFun((1:p)', kappa, t));
end

function [res, namesres] = Evp_v0_A(p, kappa, theta, sigma, v0, t)
res = zeros(length(p),1);
namesres = cell(length(p),1);
for i = 1:length(p)
    res(i) = exp(-p(i) * kappa * t) * (v0^p(i) + cFun(p(i) + 1, kappa, theta, sigma) * nuFun(p(i), kappa, theta, sigma, v0, t));
    namesres{i} = strcat("Evt^",num2str(p(i)),"|v0");
end
end

function res = Evp_v0_B(p, kappa, theta, sigma, v0, t)
res = zeros(length(p),1);
for i =1:length(p)
    res(i) = exp(-p(i) * kappa * t) * (v0^p(i) + cFun(p(i) + 1, kappa, theta, sigma) * nuFun(p(i), kappa, theta, sigma, v0, t));
end
end


function [res,colnames,rownames]=Evpvq_v0(p, q, kappa, theta, sigma, v0, t)
res = zeros(length(p), length(q));
colnames=cell(length(q),1);
rownames=cell(length(p),1);
for i = 1:length(p)
    for j = 1:length(q)
        res(i,j) = v0^p(i) * Evp_v0_B(q(j), kappa, theta, sigma, v0, t);
        colnames{j} = strcat("q = ", num2str(q(j)));
        rownames{i} = strcat("p = ", num2str(p(i)));
    end
end
end

function [res, namesres] = Evp_A(p, kappa, theta, sigma, v0)
res = zeros(length(p),1);
namesres=cell(length(p),1);
for i = 1:length(p)
    if p(i) == 1
        res(i) = cFun(2, kappa, theta, sigma) * aFun(1, kappa, v0);
    else
        res(i) = cFun(p(i) + 1, kappa, theta, sigma) * mFun(p(i) - 1, p(i), kappa, theta, sigma, v0) / (factorial(p(i)) * kappa^(p(i) - 1));
    end
    namesres{i}= strcat("Evt^",num2str(p(i)));
end
end

function res = Evp_B(p, kappa, theta, sigma, v0)
res = zeros(length(p),1);
for i = 1:length(p)
    if p(i) == 1
        res(i) = cFun(2, kappa, theta, sigma) * aFun(1, kappa, v0);
    else
        res(i) = cFun(p(i) + 1, kappa, theta, sigma) * mFun(p(i) - 1, p(i), kappa, theta, sigma, v0) / (factorial(p(i)) * kappa^(p(i) - 1));
    end
end
end

function res = Evp_m(k, p, q, kappa, theta, sigma, v0, t)
res = zeros(length(k), 1);
for i = 1:length(k)
    if k(i) == 0
        res(i) = 0;
    else
        res(i) = Evp_B(p + q - k(i) - 1, kappa, theta, sigma, v0) / kappa * prod(cFun(((q - k(i) + 1):q)', kappa, theta, sigma));
    end
end
end

function x = Evp_nuq(p, q, kappa, theta, sigma, v0, t)
x = Evp_B(p+q-1, kappa, theta, sigma, v0) .* bFun(kappa, t) ./ kappa + sum(Evp_m((0:(q-1))', p , q, kappa, theta, sigma, v0).*IFun((1:q)', kappa, t));
end

function [res, rownames, colnames] = Evpvq(p, q, kappa, theta, sigma, v0, t)
res = zeros( length(p), length(q));
rownames = cell(length(p),1);
colnames = cell(length(q),1);
for i = 1:length(p)
    for j = 1:length(q)
        res(i,j) = exp(-q(j) * kappa * t) * (Evp_B(p(i) + q(j), kappa, theta, sigma, v0) + cFun(q(j) + 1, kappa, theta, sigma) * Evp_nuq(p(i), q(j), kappa, theta, sigma, v0, t));
        colnames{j} = strcat("q = ", num2str(q(j)));
        rownames{i} = strcat("p = ", num2str(p(i)));
    end
end
end



