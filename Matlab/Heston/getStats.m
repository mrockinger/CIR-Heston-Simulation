function Errors = getStats(X,momth)

m1 = mean(X);    
v = var(X);       
sk = skewness(X); 
ku = kurtosis(X); 

mth  = momth.m;
vth  = momth.v;
skth = momth.sk;
kuth = momth.ku;

% calculate the bias (error) and absolute percentage error

Bm  = (m1-mth);
Bv  = (v-vth);
Bsk = (sk-skth);
Bku  = (ku-kuth);

APEm   = abs(Bm / mth) ;
APEv   = abs(Bv / vth) ;
APEsk  = abs( Bsk / skth);
APEku  = abs( Bku / kuth);

Errors = [Bm Bv Bsk Bku APEm APEv APEsk APEku];


