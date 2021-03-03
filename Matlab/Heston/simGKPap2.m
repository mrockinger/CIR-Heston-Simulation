function IV = simGKPap2(v1, v2, nu, Cz, K, lambda_n, gamma_n, i_X2, i_Z, q_i_X2, q_i_Z, mcp, I_X2, I_Z, D, sigma)

% Note: pdf_eta is only definied for z > 0 (see Yuan/Kalbfleisch). To
% ensure this, fix both v1 /v2 values ad hoc if necessary

v1(v1 <= 0) = 1e-5;
v2(v2 <= 0) = 1e-5;
 
%%% Simulate X1

% Mean and variance of remainder term of X1(truncation correction). See Lemma 3.1 on p. 276

EX1K = 2.*(v1+v2).*D./(pi^2.*K);
VX1K = 2.*(v1+v2).*sigma^2.*D^3./(3*pi^4*K^3);

% Draw the NUMBER of exponential r.v.s N_n for each summation index n = 1,...,K from a Poisson distribution 
% Note further: Intensity (mean) of the Poisson r.v. is a function of the summation index n

N_n =  arrayfun(@poissrnd, bsxfun(@times,  (v1 + v2)',  lambda_n));


% Draw N_K exponential variables with unit mean for each K. Take sum of them afterwards 

% mu = 1;
% sz1 = 1;
Exp = arrayfun(@(sz2) exprnd(1, 1, sz2), N_n, 'UniformOutput', false);
SumExp = cellfun(@sum, Exp);

% Finally, calculate X1

X1 = sum(bsxfun(@times,  gamma_n.^-1,  SumExp), 2)' + gamrnd(EX1K.^2 ./ VX1K, VX1K ./ EX1K);

%%% Simulate X2

UX2 = unifrnd(0, 1, 1, length(v1));
LcpX2 = ceil(UX2 .* mcp);
X2 = I_X2(LcpX2);

% Logical indexing
LIX2 = (UX2' > interp1(i_X2, q_i_X2, X2));
 
while sum(LIX2) > 0
    
    X2(LIX2) = i_X2(arrayfun(@(x) find(i_X2==x), X2(LIX2)) + 1);
    LIX2 = (UX2' > interp1(i_X2, q_i_X2, X2));
    
end


%%% Simulate X3

% On the other side, if z is getting too large, pdf_eta gives values of 0 all the time. 
% This is because of the Bessel-function term. 
% Therefore, make sure, to keep z "reasonable" small.
% I tested the function. It crashed for values of z > 705


z = Cz * sqrt(v1 .* v2);
z(z>705) = 705;

% First, simulate the Bessel r.v. eta

Ueta = unifrnd(0,1, 1, length(v1));
eta = ones(1, length(v1));
% eta = poissrnd(5, 1, 100);

% calculate pdf values for each eta and each z
% This will result in a cell array since for each z, we will have a (possibly) different number of pdf
% values/evaluations, depending on the value of eta!
pdfs = arrayfun(@(n, x) pdf_eta(0:(n-1), nu, x), eta, z, 'UniformOutput', false);
CDF = cellfun(@sum, pdfs);


LIZ = (Ueta - CDF > pdf_eta(eta, nu, z));

while sum(LIZ) > 0
    
    eta(LIZ) = eta(LIZ) + 1;
    pdfs = arrayfun(@(n, x) pdf_eta(0:(n-1), nu, x), eta, z, 'UniformOutput', false);
    CDF = cellfun(@sum, pdfs);
    LIZ = (Ueta - CDF > pdf_eta(eta, nu, z));
    
end    

% If we know the number of Z we have to draw (which is eta) for each z, proceed as in X2, i. e. draw from the pre-computed distribution table

UZ = arrayfun(@(x) unifrnd(0, 1, 1, x), eta, 'UniformOutput', false);
LcpZ = cellfun(@(x) ceil(x .* mcp), UZ, 'UniformOutput', false);
Z = cellfun(@(x) I_Z(x), LcpZ, 'UniformOutput', false);

% Store Logical Indizes for each eta and each z
LIX3 = cellfun(@(x, y) (x > interp1(i_Z, q_i_Z, y)'), UZ, Z, 'UniformOutput', false);
 
% Calculate the sum as it is the driver of the while loop
S = cellfun(@sum, LIX3);

for i=length(v1)
   
    while S(i) > 0
    
        Z{i}(LIX3{i}) = i_Z(arrayfun(@(x) find(i_Z==x), Z{i}(LIX3{i})) + 1);
        LIX3{i} = (UZ{i} > interp1(i_Z, q_i_Z, Z{i})');
        S(i) = sum(LIX3{i});
        
    end
end


 X3 = cellfun(@sum, Z);
 
 IV = X1 + X2' + X3;



