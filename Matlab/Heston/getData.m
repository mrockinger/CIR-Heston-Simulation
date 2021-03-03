% Simulation for return processes

n = [12 4 2 1];
T = 1;

M = 1e4; 
NRuns = 100;

kappaHNC = 3;
theta = 0.19;
sigmaHNC = 0.4;
rho = -0.7;
mu = 0;
v0HNC = 0.1;

% RetM={'EulerFT'; 'Zhu'; 'KJ'; 'QE';  'BK'; 'Smith'; 'GKPap'; 'TW'; 'BBG'; 'ORS'; 'ORS + SUB'};
% Exclude ORS in this conditional setting, since it does not make any sense
RetM={'EulerFT'; 'Zhu'; 'KJ'; 'QE';  'BK'; 'Smith'; 'GKPap'; 'TW'; 'BBG'};

[MAEmeanHNC, MAEvarHNC, MAEskewHNC, MAEkurtHNC, RMSEmeanHNC, RMSEvarHNC, RMSEskewHNC, RMSEkurtHNC, MAPEmeanHNC, MAPEvarHNC, MAPEskewHNC, MAPEkurtHNC, MetimeHNC] = getMErrors(v0HNC, mu, kappaHNC, theta, sigmaHNC, rho, n, T, M, NRuns, RetM, 'C'); 

save('ResultsReturnsHNCC.mat', 'MAEmeanHNC', 'MAEvarHNC', 'MAEskewHNC', 'MAEkurtHNC', 'RMSEmeanHNC', 'RMSEvarHNC', 'RMSEskewHNC', 'RMSEkurtHNC', 'MAPEmeanHNC', 'MAPEvarHNC', 'MAPEskewHNC', 'MAPEkurtHNC', 'MetimeHNC')

% Same for LNC parameters

kappaLNC = 0.5;
sigmaLNC = 1;
v0LNC = 0.05;

[MAEmeanLNC, MAEvarLNC, MAEskewLNC, MAEkurtLNC, RMSEmeanLNC, RMSEvarLNC, RMSEskewLNC, RMSEkurtLNC, MAPEmeanLNC, MAPEvarLNC, MAPEskewLNC, MAPEkurtLNC, MetimeLNC] = getMErrors(v0LNC, mu, kappaLNC, theta, sigmaLNC, rho, n, T, M, NRuns, RetM, 'C'); 

save('ResultsReturnsLNCC.mat', 'MAEmeanLNC', 'MAEvarLNC', 'MAEskewLNC', 'MAEkurtLNC', 'RMSEmeanLNC', 'RMSEvarLNC', 'RMSEskewLNC', 'RMSEkurtLNC', 'MAPEmeanLNC', 'MAPEvarLNC', 'MAPEskewLNC', 'MAPEkurtLNC', 'MetimeLNC')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Now, do the same for the UC setting

RetMUC={'EulerFT'; 'Zhu'; 'KJ'; 'QE';  'BK'; 'Smith'; 'GKPap'; 'TW'; 'BBG'; 'ORS + SUB'};

% HNC

n = 1;
T = [12 4 2 1].^-1;

shapeHNC = 2 * kappaHNC * theta / sigmaHNC^2;
scaleHNC = theta / shapeHNC; 

v0HNCUC = gamrnd(shapeHNC, scaleHNC, 1, M);

[MAEmeanHNCUC, MAEvarHNCUC, MAEskewHNCUC, MAEkurtHNCUC, RMSEmeanHNCUC, RMSEvarHNCUC, RMSEskewHNCUC, RMSEkurtHNCUC, MAPEmeanHNCUC, MAPEvarHNCUC, MAPEskewHNCUC, MAPEkurtHNCUC, MetimeHNCUC] = getMErrors(v0HNCUC, mu, kappaHNC, theta, sigmaHNC, rho, n, T, M, NRuns, RetMUC, 'UC'); 

save('ResultsReturnsHNCUC.mat', 'MAEmeanHNCUC', 'MAEvarHNCUC', 'MAEskewHNCUC', 'MAEkurtHNCUC', 'RMSEmeanHNCUC', 'RMSEvarHNCUC', 'RMSEskewHNCUC', 'RMSEkurtHNCUC', 'MAPEmeanHNCUC', 'MAPEvarHNCUC', 'MAPEskewHNCUC', 'MAPEkurtHNCUC', 'MetimeHNCUC')

% LNC

shapeLNC = 2 * kappaLNC * theta / sigmaLNC^2;
scaleLNC = theta / shapeLNC; 

v0LNCUC = gamrnd(shapeLNC, scaleLNC, 1, M);

[MAEmeanLNCUC, MAEvarLNCUC, MAEskewLNCUC, MAEkurtLNCUC, RMSEmeanLNCUC, RMSEvarLNCUC, RMSEskewLNCUC, RMSEkurtLNCUC, MAPEmeanLNCUC, MAPEvarLNCUC, MAPEskewLNCUC, MAPEkurtLNCUC, MetimeLNCUC] = getMErrors(v0LNCUC, mu, kappaLNC, theta, sigmaLNC, rho, n, T, M, NRuns, RetMUC, 'UC'); 

save('ResultsReturnsLNCUC.mat', 'MAEmeanLNCUC', 'MAEvarLNCUC', 'MAEskewLNCUC', 'MAEkurtLNCUC', 'RMSEmeanLNCUC', 'RMSEvarLNCUC', 'RMSEskewLNCUC', 'RMSEkurtLNCUC', 'MAPEmeanLNCUC', 'MAPEvarLNCUC', 'MAPEskewLNCUC', 'MAPEkurtLNCUC', 'MetimeLNCUC')
