% Simulation for the CIR process

n = [12 4 2 1];
T = 1;

M = 1e5; 
NRuns = 100;

kappaHNC = 3;
theta = 0.19;
sigmaHNC = 0.4;
v0HNCC = 0.1;

VarMHNC={'AE0'; 'DD'; 'EA'; 'ER'; 'EFT'; 'TVS'; 'HM'; 'KJ'; 'G'; 'ABR'; 'QE'; 'SAS'; 'SAN'; 'BK'};

[MAEmeanHNCC, MAEvarHNCC, MAEskewHNCC, MAEkurtHNCC, RMSEmeanHNCC, RMSEvarHNCC, RMSEskewHNCC, RMSEkurtHNCC, MAPEmeanHNCC, MAPEvarHNCC, MAPEskewHNCC, MAPEkurtHNCC, MetimeHNCC] = getMErrors(v0HNCC, kappaHNC, theta, sigmaHNC, n, T, M, NRuns, VarMHNC, 'C'); 

save('ResultsVarianceHNCC.mat', 'MAEmeanHNCC', 'MAEvarHNCC', 'MAEskewHNCC', 'MAEkurtHNCC', 'RMSEmeanHNCC', 'RMSEvarHNCC', 'RMSEskewHNCC', 'RMSEkurtHNCC', 'MAPEmeanHNCC', 'MAPEvarHNCC', 'MAPEskewHNCC', 'MAPEkurtHNCC', 'MetimeHNCC')

% Same for LNC parameters

kappaLNC = 0.5;
sigmaLNC = 1;
v0LNCC = 0.05;

VarMLNC={'AE0'; 'DD'; 'EA'; 'ER'; 'EFT'; 'TVS'; 'HM'; 'KJ'; 'G'; 'ABR'; 'QE'; 'BK'};

[MAEmeanLNCC, MAEvarLNCC, MAEskewLNCC, MAEkurtLNCC, RMSEmeanLNCC, RMSEvarLNCC, RMSEskewLNCC, RMSEkurtLNCC, MAPEmeanLNCC, MAPEvarLNCC, MAPEskewLNCC, MAPEkurtLNCC, MetimeLNCC] = getMErrors(v0LNCC, kappaLNC, theta, sigmaLNC, n, T, M, NRuns, VarMLNC, 'C'); 

save('ResultsVarianceLNCC.mat', 'MAEmeanLNCC', 'MAEvarLNCC', 'MAEskewLNCC', 'MAEkurtLNCC', 'RMSEmeanLNCC', 'RMSEvarLNCC', 'RMSEskewLNCC', 'RMSEkurtLNCC', 'MAPEmeanLNCC', 'MAPEvarLNCC', 'MAPEskewLNCC', 'MAPEkurtLNCC', 'MetimeLNCC')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Now, do the same for the UC setting

% HNC

n = 1;
T = [12 4 2 1].^-1;

shapeHNC = 2 * kappaHNC * theta / sigmaHNC^2;
scaleHNC = theta / shapeHNC; 

v0HNCUC = gamrnd(shapeHNC, scaleHNC, 1, M);

[MAEmeanHNCUC, MAEvarHNCUC, MAEskewHNCUC, MAEkurtHNCUC, RMSEmeanHNCUC, RMSEvarHNCUC, RMSEskewHNCUC, RMSEkurtHNCUC, MAPEmeanHNCUC, MAPEvarHNCUC, MAPEskewHNCUC, MAPEkurtHNCUC, MetimeHNCUC] = getMErrors(v0HNCUC, kappaHNC, theta, sigmaHNC, n, T, M, NRuns, VarMHNC, 'UC'); 

save('ResultsVarianceHNCUC.mat', 'MAEmeanHNCUC', 'MAEvarHNCUC', 'MAEskewHNCUC', 'MAEkurtHNCUC', 'RMSEmeanHNCUC', 'RMSEvarHNCUC', 'RMSEskewHNCUC', 'RMSEkurtHNCUC', 'MAPEmeanHNCUC', 'MAPEvarHNCUC', 'MAPEskewHNCUC', 'MAPEkurtHNCUC', 'MetimeHNCUC')

% LNC

% Glasserman suffers from erratic behaviour here, so I leave this approach
% out of further considerations

VarMLNC={'AE0'; 'DD'; 'EA'; 'ER'; 'EFT'; 'TVS'; 'HM'; 'KJ'; 'ABR'; 'QE'; 'BK'};

shapeLNC = 2 * kappaLNC * theta / sigmaLNC^2;
scaleLNC = theta / shapeLNC; 

v0LNCUC = gamrnd(shapeLNC, scaleLNC, 1, M);

[MAEmeanLNCUC, MAEvarLNCUC, MAEskewLNCUC, MAEkurtLNCUC, RMSEmeanLNCUC, RMSEvarLNCUC, RMSEskewLNCUC, RMSEkurtLNCUC, MAPEmeanLNCUC, MAPEvarLNCUC, MAPEskewLNCUC, MAPEkurtLNCUC, MetimeLNCUC] = getMErrors(v0LNCUC, kappaLNC, theta, sigmaLNC, n, T, M, NRuns, VarMLNC, 'UC'); 

save('ResultsVarianceLNCUC.mat', 'MAEmeanLNCUC', 'MAEvarLNCUC', 'MAEskewLNCUC', 'MAEkurtLNCUC', 'RMSEmeanLNCUC', 'RMSEvarLNCUC', 'RMSEskewLNCUC', 'RMSEkurtLNCUC', 'MAPEmeanLNCUC', 'MAPEvarLNCUC', 'MAPEskewLNCUC', 'MAPEkurtLNCUC', 'MetimeLNCUC')


