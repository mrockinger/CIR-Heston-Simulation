% Comparison QE/ORS with huge samples

%HNC

n = 1;
T = [12 4 2 1].^-1;

M = 1e6; 
NRuns = 100;

kappaHNC = 3;
theta = 0.19;
sigmaHNC = 0.4;
rho = -0.7;
mu = 0;

shapeHNC = 2 * kappaHNC * theta / sigmaHNC^2;
scaleHNC = theta / shapeHNC; 

v0HNCUC = gamrnd(shapeHNC, scaleHNC, 1, M);

RetMUC={'QE'; 'ORS + SUB'};

[MAEmeanHNCUC, MAEvarHNCUC, MAEskewHNCUC, MAEkurtHNCUC, RMSEmeanHNCUC, RMSEvarHNCUC, RMSEskewHNCUC, RMSEkurtHNCUC, MAPEmeanHNCUC, MAPEvarHNCUC, MAPEskewHNCUC, MAPEkurtHNCUC, MetimeHNCUC] = getMErrors(v0HNCUC, mu, kappaHNC, theta, sigmaHNC, rho, n, T, M, NRuns, RetMUC, 'UC'); 

save('ResultsCompQEORSHNC.mat', 'MAEmeanHNCUC', 'MAEvarHNCUC', 'MAEskewHNCUC', 'MAEkurtHNCUC', 'RMSEmeanHNCUC', 'RMSEvarHNCUC', 'RMSEskewHNCUC', 'RMSEkurtHNCUC', 'MAPEmeanHNCUC', 'MAPEvarHNCUC', 'MAPEskewHNCUC', 'MAPEkurtHNCUC', 'MetimeHNCUC')


% LNC

kappaLNC = 0.5;
sigmaLNC = 1;

shapeLNC = 2 * kappaLNC * theta / sigmaLNC^2;
scaleLNC = theta / shapeLNC; 

v0LNCUC = gamrnd(shapeLNC, scaleLNC, 1, M);

[MAEmeanLNCUC, MAEvarLNCUC, MAEskewLNCUC, MAEkurtLNCUC, RMSEmeanLNCUC, RMSEvarLNCUC, RMSEskewLNCUC, RMSEkurtLNCUC, MAPEmeanLNCUC, MAPEvarLNCUC, MAPEskewLNCUC, MAPEkurtLNCUC, MetimeLNCUC] = getMErrors(v0LNCUC, mu, kappaLNC, theta, sigmaLNC, rho, n, T, M, NRuns, RetMUC, 'UC'); 

save('ResultsCompQEORSLNC.mat', 'MAEmeanLNCUC', 'MAEvarLNCUC', 'MAEskewLNCUC', 'MAEkurtLNCUC', 'RMSEmeanLNCUC', 'RMSEvarLNCUC', 'RMSEskewLNCUC', 'RMSEkurtLNCUC', 'MAPEmeanLNCUC', 'MAPEvarLNCUC', 'MAPEskewLNCUC', 'MAPEkurtLNCUC', 'MetimeLNCUC')
