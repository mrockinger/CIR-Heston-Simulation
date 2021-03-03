library(doParallel)
library(ggplot2)
library(EnvStats)

source("SupplementCodeV2.r")
source("DiscretizationConditionalReturnMomentsV2.r")
source("QEScheme.r")
source('../CIR/Exact.r')
source("getSkFromMoments.r")4
source("getKuFromMoments.r")

source("TreatReturns.r")


set.seed(123)
test = TreatReturns('QE', 100, 0.1, 3, 0.19, 0.4, 0, -0.7, 1, 1, 1e6, momth)
set.seed(123)
test1 = TreatReturns('ORS', 100, 0.1, 3, 0.19, 0.4, 0, -0.7, 1, 1, 1e6, momth)
set.seed(123)
test2 = TreatReturns(Model = 'ORS + SUB', 100, 0.1, 3, 0.19, 0.4, 0, -0.7, 1, 1, 1e6, momth)

# > test
# Ert.v0   Varrt.v0  Skewrt.v0  Kurtrt.v0            
# 1.33388129 0.02388774 0.03758334 0.01033055 0.80469298 
# > test1
# Ert.v0     Varrt.v0    Skewrt.v0    Kurtrt.v0              
# 2.21604977   0.21012828   0.13333494   0.03457302 283.02202702 
# > test2
# Ert.v0     Varrt.v0    Skewrt.v0    Kurtrt.v0              
# 2.221786e+00 2.096618e-01 1.974465e-03 5.876428e-03 2.685515e+02 

# Fehler in Varianz und Mean bei optimiertem QE können durch rescaling kompensiert werden und treten wahrscheinlich nur bei hohen Deltas auf ( siehe unten )
# Es bringt hier sicherlich sehr viel über die SUB scheme V zu generieren, da die QE scheme bei hohem Delta doch einen recht hohen Schiefe - Bias hat ( siehe Grafiken zu CIR-Simulation )
# Das könnte sich mit kleinerem Delta ändern (siehe auch CIR-Grafiken)


set.seed(123)
test = TreatReturns('QE', 100, 0.1, 3, 0.19, 0.4, 0, -0.7, 1, 2, 1e6, momth)
set.seed(123)
test1 = TreatReturns('ORS', 100, 0.1, 3, 0.19, 0.4, 0, -0.7, 1, 2, 1e6, momth)
set.seed(123)
test2 = TreatReturns(Model = 'ORS + SUB', 100, 0.1, 3, 0.19, 0.4, 0, -0.7, 1, 2, 1e6, momth)

# > test
# Ert.v0   Varrt.v0  Skewrt.v0  Kurtrt.v0            
# 0.36679736 0.03487483 0.04059892 0.02563560 1.75327706 
# > test1
# Ert.v0     Varrt.v0    Skewrt.v0    Kurtrt.v0              
# 6.520257e-01 4.800361e-03 3.894437e-02 2.454327e-02 4.767690e+02 
# > test2
# Ert.v0     Varrt.v0    Skewrt.v0    Kurtrt.v0              
# 6.466142e-01 6.707003e-03 6.430275e-02 1.204946e-02 5.207602e+02 

set.seed(123)
test = TreatReturns('QE', 100, 0.1, 3, 0.19, 0.4, 0, -0.7, 1, 4, 1e6, momth)
set.seed(123)
test1 = TreatReturns('ORS', 100, 0.1, 3, 0.19, 0.4, 0, -0.7, 1, 4, 1e6, momth)
set.seed(123)
test2 = TreatReturns(Model = 'ORS + SUB', 100, 0.1, 3, 0.19, 0.4, 0, -0.7, 1, 4, 1e6, momth)

# > test
# Ert.v0   Varrt.v0  Skewrt.v0  Kurtrt.v0            
# 0.09268304 0.01214243 0.03820370 0.01807258 4.98224401 
# > test1
# Ert.v0     Varrt.v0    Skewrt.v0    Kurtrt.v0              
# 1.862090e-01 7.790213e-03 2.814202e-02 1.607309e-02 1.448701e+03 
# > test2
# Ert.v0     Varrt.v0    Skewrt.v0    Kurtrt.v0              
# 1.905326e-01 4.107601e-03 2.074294e-02 1.984074e-03 1.121811e+03 

# Fazit:
# Der Fehler im mean ist nahezu immer doppelt so groß wie bei der Standard QE scheme. 

# Für Delta = 1:

# Für beide modifizierten QEs ist Fehler in der Varianz ist knapp 10x größer als für Standard QE (wieso?!)
# Das ist an sich nicht schlimm, da beides skaliert werden könnte 
# Für die optimierte QE ist allerdings auch Fehler in Schiefe größer als bei normaler QE
# Die opt QE + SUB hat hingegen einen viel kleineren Fehler in der Schiefe als die QE
# Alles für die Schiefe gilt acuh für die Kurtosis

# Für Delta = 0.5:

# opt QE hat kleineren Skewness-Bias als normale, wohingegen die opt QE + SUB einen größeren Skewness-Bias hat
# Varianzen sind hier bei beiden deutlich besser als bei der normalen QE

# Für Delta  = 1 / 4

# Var / Schiefe / Kurtosis sind bei beiden besser als bei QE
# QE + SUB macht besten Job 


names = c('QE', 'ORS', 'ORS + SUB')
n = c(4, 2, 1)
M = 1e6
T = 1
NRuns = 100
S0 = 100
v0 = 0.1
kappa = 3
theta = 0.19
sigma = 0.4
mu = 0
rho  = -0.7

momth = MomentsBates(mu = 0, kappa = 3, theta = 0.19, sigma = 0.4, rho = -0.7, lambda = 0, muj = 0, vj = 0, t = 1, v0 = 0.1, conditional = TRUE)[c(1, 5:7)]


cl = makeCluster(length(names))
registerDoParallel(cl)

APEs = list()
Results = list()

clusterExport(cl,ls(globalenv()))

for(k in 1:NRuns){
  
  for(j in 1:length(n)){
    
          APEs[[j]] = foreach(i = 1:length(names), .combine = 'cbind') %dopar% TreatReturns(Model = names[i], S0 = S0, v0 = v0, kappa = kappa, theta = theta, sigma = sigma, mu = mu, rho = rho, T = T, n = n[j], M = M, momth = momth)
  }
  
  
  Results[[k]] = APEs
  
  print(k)
  
}

stopCluster(cl)

APEsmean = list() 
APEsvar  = list()
APEsskew = list()
APEskurt = list()
etimes   = list()


for(i in 1:NRuns){
  
  APEsmean[[i]] = do.call(rbind, lapply(Results[[i]], `[`, 1,))
  APEsvar[[i]]  = do.call(rbind, lapply(Results[[i]], `[`, 2,))
  APEsskew[[i]] = do.call(rbind, lapply(Results[[i]], `[`, 3,))
  APEskurt[[i]] = do.call(rbind, lapply(Results[[i]], `[`, 4,))
  etimes[[i]]   = do.call(rbind, lapply(Results[[i]], `[`, 5,))
  
} 


# rows = n; column = names
MAPEmean = Reduce('+', APEsmean) / NRuns # Alternatively: MAPEmean = apply(simplify2array(APEsmean), 1:2, mean)
MAPEvar  = Reduce('+', APEsvar) / NRuns
MAPEskew = Reduce('+', APEsskew) / NRuns
MAPEkurt = Reduce('+', APEskurt) / NRuns

Metime = Reduce('+', etimes) / NRuns

q05APEsmean = apply(simplify2array(APEsmean), 1:2, quantile, 0.05)
q95APEsmean = apply(simplify2array(APEsmean), 1:2, quantile, 0.95)

q05APEsvar = apply(simplify2array(APEsvar), 1:2, quantile, 0.05)
q95APEsvar = apply(simplify2array(APEsvar), 1:2, quantile, 0.95)

q05APEsskew = apply(simplify2array(APEsskew), 1:2, quantile, 0.05)
q95APEsskew = apply(simplify2array(APEsskew), 1:2, quantile, 0.95)

q05APEskurt = apply(simplify2array(APEskurt), 1:2, quantile, 0.05)
q95APEskurt = apply(simplify2array(APEskurt), 1:2, quantile, 0.95)

Moments = factor(c("Mean", "Variance", "Skewness", "Kurtosis"), levels = c("Mean", "Variance", "Skewness", "Kurtosis"))
groups =  factor(c("Explicit", "Implicit", "Others"), levels = c("Explicit", "Implicit", "Others")) 

Delta = T / n

df.MAPE = data.frame("Scheme" = rep(rep(names, each = length(Delta)), length(Moments)), "Group" = rep(rep(groups, c(0 * length(Delta), 0 * length(Delta), 3 * length(Delta))), length(Moments)), "Delta" = rep(rep(Delta, length(names)), length(Moments)), "Moment" = rep(Moments, each = length(names) * length(Delta)) ,"MAPE" = c(c(MAPEmean), c(MAPEvar), c(MAPEskew), c(MAPEkurt)), "q05" = c(c(q05APEsmean), c(q05APEsvar), c(q05APEsskew), c(q05APEskurt)), "q95" = c(c(q95APEsmean), c(q95APEsvar), c(q95APEsskew), c(q95APEskurt)), "etime" = rep(c(Metime), length(Moments)))

gg.MAPEHNC = ggplot(df.MAPE, aes(x = Delta, y = MAPE, colour = Scheme, shape = Scheme)) + geom_point() + geom_line() + facet_grid(.~Moment, scales = "free") + theme_bw() +  scale_x_reverse() + xlab(expression(Delta))
gg.MAPEHNC.cb = gg.MAPEHNC + geom_ribbon(aes(ymin = q05, ymax = q95,  fill = Scheme), alpha = 0.1, colour = NA)


gg.MAPE.var =  ggplot(subset(df.MAPE, Moment == "Variance"), aes(Delta, y = MAPE, colour = Scheme, shape = Scheme)) + geom_point() + geom_line() + theme_bw() +  scale_x_reverse() + xlab(expression(Delta)) + geom_ribbon(aes(ymin = q05, ymax = q95,  fill = Scheme), alpha = 0.1, colour = NA)
gg.MAPE.skew =  ggplot(subset(df.MAPE, Moment == "Skewness"), aes(Delta, y = MAPE, colour = Scheme, shape = Scheme)) + geom_point() + geom_line() + theme_bw() +  scale_x_reverse() + xlab(expression(Delta)) + geom_ribbon(aes(ymin = q05, ymax = q95,  fill = Scheme), alpha = 0.1, colour = NA)

pdf("QEComparison.pdf")
plot(gg.MAPEHNC.cb)
plot(gg.MAPE.var)
plot(gg.MAPE.skew)
dev.off()

# Test whether my theoretical QE moments are flawed.
# Therefore, simulate really huge sample of conditional returns and compare
# To this for many gammas / deltas

T = 1

RetQEemp = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n = 1, M = 1e7, .5, .5, optimized = F)[,2] - log(100)

# Do the return simulation with exact simulation of the variance
v = Exact(v0, kappa, theta, sigma, T, n  = 1, M = 1e7)
mean(v[,2])
mean(v[,2]^2)
mean(v[,2]^3)
mean(v[,2]^4)
MomentsCIR(p = 1:4, kappa = kappa, theta = theta, sigma = sigma, v0 = v0, t = T)
# > mean(v[,2])
# [1] 0.1855178
# > mean(v[,2]^2)
# [1] 0.03924366
# > mean(v[,2]^3)
# [1] 0.009322297
# > mean(v[,2]^4)
# [1] 0.002456967
# > MomentsCIR(p = 1:4, kappa = kappa, theta = theta, sigma = sigma, v0 = v0, t = T)
# Evt^1|v0    Evt^2|v0    Evt^3|v0    Evt^4|v0 
# 0.185519164 0.039244388 0.009322616 0.002457098 

RetQESUBemp = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n = 1, M = 1e7, .5, .5, optimized = F, v = v)[,2] - log(100)

MomQETheory = c(
                  CondMomQE("mu1", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
                  CondMomQE("mu2", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
                  CondMomQE("mu3", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
                  CondMomQE("mu4", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
                  CondMomQE("var", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
                  CondMomQE("skew", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
                  CondMomQE("kurt", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5)
)

ResQE = cbind(MomQETheory, 
              c(mean(RetQEemp), mean(RetQEemp^2), mean(RetQEemp^3), mean(RetQEemp^4), var(RetQEemp), skewness(RetQEemp), kurtosis(RetQEemp, excess = F)),
              c(mean(RetQESUBemp), mean(RetQESUBemp^2), mean(RetQESUBemp^3), mean(RetQESUBemp^4), var(RetQESUBemp), skewness(RetQESUBemp), kurtosis(RetQESUBemp, excess = F)),
              MomentsBates(mu, kappa, theta, sigma, rho, lambda = 0, muj = 0, vj = 0, t = T, v0 = v0)
              )
# > ResQE
# MomQETheory                                    
# Ert.v0     0.02697387  0.02693350  0.02714965 -0.08074681
# Ert2.v0    0.17678812  0.17688730  0.17678194  0.17824452
# Ert3.v0   -0.02764548 -0.02177093 -0.02758692 -0.07568704
# Ert4.v0    0.10734910  0.10185684  0.10741421  0.11869089
# Varrt.v0   0.17606053  0.17616191  0.17604485  0.17172447
# Skewrt.v0 -0.56734536 -0.48722476 -0.56787345 -0.47162731
# Kurtrt.v0  3.58425009  3.38254831  3.58773629  3.42803659

# Ergo: The matching will only work, if the variance simulation is conducted exactly. 
# Only then it is assured, that the simulated variance has the same moments as the theoretical diffusion 

get.opt.gamma1(Delta = T, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, vt = v0, target = "skew", mode = "conditional")
# [1] 0.645008

print(CondMomQE("skew", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.645008, 1-0.645008))
# [1] -0.471619

RetQESUBoptemp = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n = 1, M = 1e7, optimized = T, target = "skew", mode = "conditional", v = v)[,2] - log(100)

c(mean(RetQESUBoptemp), mean(RetQESUBoptemp^2), mean(RetQESUBoptemp^3), mean(RetQESUBoptemp^4), var(RetQESUBoptemp), skewness(RetQESUBoptemp), kurtosis(RetQESUBoptemp, excess = F))
# [1]  0.09835956  0.14561252  0.01746968  0.06256066  0.13593793 -0.47075619  3.45574452

# Do scaling

mu1 =  MomentsBates(mu, kappa, theta, sigma, rho, lambda = 0, muj = 0, vj = 0, t = T, v0 = v0)[1]
mu1.disc =  CondMomQE("mu1", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.645008, 1-0.645008)
zeta1 = MomentsBates(mu, kappa, theta, sigma, rho, lambda = 0, muj = 0, vj = 0, t = T, v0 = v0)[5]
zeta1.disc = CondMomQE("var", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.645008, 1-0.645008)
b = sqrt(zeta1 / zeta1.disc)
a = mu1 * (1 - b)

RetQESUBoptemp.scaled = a + b * RetQESUBoptemp
c(mean(RetQESUBoptemp.scaled), var(RetQESUBoptemp.scaled), skewness(RetQESUBoptemp.scaled), kurtosis(RetQESUBoptemp.scaled, excess = F))
# 0.1206148  0.1718194 -0.4707562  3.4557445 # Variance is working, mean is not

a.new = mu1 - mu1.disc * b
RetQESUBoptemp.scaled.new = a.new + b * RetQESUBoptemp
c(mean(RetQESUBoptemp.scaled.new), var(RetQESUBoptemp.scaled.new), skewness(RetQESUBoptemp.scaled.new), kurtosis(RetQESUBoptemp.scaled.new, excess = F))
# [1] -0.08065669  0.17181935 -0.47075619  3.45574452
# Which is absolutely perfect!


# Do the same for T = 0.5

T = 0.5

RetQEemp = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n = 1, M = 1e7, .5, .5, optimized = F)[,2] - log(100)

# Do the return simulation with exact simulation of the variance
v = Exact(v0, kappa, theta, sigma, T, n  = 1, M = 1e7)

RetQESUBemp = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n = 1, M = 1e7, .5, .5, optimized = F, v = v)[,2] - log(100)

MomQETheory = c(
  CondMomQE("mu1", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("mu2", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("mu3", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("mu4", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("var", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("skew", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("kurt", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5)
)

ResQE = cbind(MomQETheory, 
              c(mean(RetQEemp), mean(RetQEemp^2), mean(RetQEemp^3), mean(RetQEemp^4), var(RetQEemp), skewness(RetQEemp), kurtosis(RetQEemp, excess = F)),
              c(mean(RetQESUBemp), mean(RetQESUBemp^2), mean(RetQESUBemp^3), mean(RetQESUBemp^4), var(RetQESUBemp), skewness(RetQESUBemp), kurtosis(RetQESUBemp, excess = F)),
              MomentsBates(mu, kappa, theta, sigma, rho, lambda = 0, muj = 0, vj = 0, t = T, v0 = v0)
)
# > ResQE
#           MomQETheory    
# Ert.v0    -0.01161454 -0.01157245 -0.01159412 -0.03584695
# Ert2.v0    0.07501093  0.07497602  0.07503301  0.07609649
# Ert3.v0   -0.01342981 -0.01201863 -0.01341489 -0.01802916
# Ert4.v0    0.02031242  0.01928015  0.02030751  0.02128456
# Varrt.v0   0.07487604  0.07484211  0.07489860  0.07481148
# Skewrt.v0 -0.52806156 -0.46001799 -0.52728125 -0.48566637
# Kurtrt.v0  3.52259583  3.35348458  3.51988513  3.44505631


get.opt.gamma1(Delta = T, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, vt = v0, target = "skew", mode = "conditional")
# [1] 0.5771123

print(CondMomQE("skew", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5771123, 1-0.5771123))
# [1] -0.4856671

RetQESUBoptemp = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n = 1, M = 1e7, optimized = T, target = "skew", mode = "conditional", v = v)[,2] - log(100)

c(mean(RetQESUBoptemp), mean(RetQESUBoptemp^2), mean(RetQESUBoptemp^3), mean(RetQESUBoptemp^4), var(RetQESUBoptemp), skewness(RetQESUBoptemp), kurtosis(RetQESUBoptemp, excess = F))
# [1]  0.003906993  0.068042268 -0.007807213  0.015907225  0.068027010 -0.484964037  3.465126041

# Do scaling

mu1 =  MomentsBates(mu, kappa, theta, sigma, rho, lambda = 0, muj = 0, vj = 0, t = T, v0 = v0)[1]
mu1.disc =  CondMomQE("mu1", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5771123, 1-0.5771123)
zeta1 = MomentsBates(mu, kappa, theta, sigma, rho, lambda = 0, muj = 0, vj = 0, t = T, v0 = v0)[5]
zeta1.disc = CondMomQE("var", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5771123, 1-0.5771123)
b = sqrt(zeta1 / zeta1.disc)
a = mu1 * (1 - b)

RetQESUBoptemp.scaled = a + b * RetQESUBoptemp
c(mean(RetQESUBoptemp.scaled), var(RetQESUBoptemp.scaled), skewness(RetQESUBoptemp.scaled), kurtosis(RetQESUBoptemp.scaled, excess = F))
# 0.005829536  0.074765829 -0.484964037  3.465126041 # Variance is working, mean is not

a.new = mu1 - mu1.disc * b
RetQESUBoptemp.scaled.new = a.new + b * RetQESUBoptemp
c(mean(RetQESUBoptemp.scaled.new), var(RetQESUBoptemp.scaled.new), skewness(RetQESUBoptemp.scaled.new), kurtosis(RetQESUBoptemp.scaled.new, excess = F))
# [1] -0.03582515  0.07476583 -0.48496404  3.46512604
# Which is absolutely perfect!



# Do the same for T = 4

T = 4

RetQEemp = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n = 1, M = 1e7, .5, .5, optimized = F)[,2] - log(100)

# Do the return simulation with exact simulation of the variance
v = Exact(v0, kappa, theta, sigma, T, n  = 1, M = 1e7)

RetQESUBemp = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n = 1, M = 1e7, .5, .5, optimized = F, v = v)[,2] - log(100)

MomQETheory = c(
  CondMomQE("mu1", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("mu2", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("mu3", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("mu4", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("var", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("skew", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("kurt", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5)
)

ResQE = cbind(MomQETheory, 
              c(mean(RetQEemp), mean(RetQEemp^2), mean(RetQEemp^3), mean(RetQEemp^4), var(RetQEemp), skewness(RetQEemp), kurtosis(RetQEemp, excess = F)),
              c(mean(RetQESUBemp), mean(RetQESUBemp^2), mean(RetQESUBemp^3), mean(RetQESUBemp^4), var(RetQESUBemp), skewness(RetQESUBemp), kurtosis(RetQESUBemp, excess = F)),
              MomentsBates(mu, kappa, theta, sigma, rho, lambda = 0, muj = 0, vj = 0, t = T, v0 = v0)
)
# > ResQE
# MomQETheory                                 
# Ert.v0      0.4975073  0.4970245  0.4970216 -0.3650001
# Ert2.v0     1.4328245  1.4328049  1.4338332  0.9282700
# Ert3.v0     1.0582293  1.2071888  1.0554060 -1.1408104
# Ert4.v0     5.3491707  5.2593204  5.3619016  2.9950944
# Varrt.v0    1.1853109  1.1857716  1.1868029  0.7950449
# Skewrt.v0  -0.6462855 -0.5294710 -0.6473569 -0.3126076
# Kurtrt.v0   3.6921534  3.4137656  3.6959833  3.1929892


get.opt.gamma1(Delta = T, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, vt = v0, target = "skew", mode = "conditional")
# [1] 0.8310602

print(CondMomQE("skew", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.8310602, 1-0.8310602))
# [1] -0.3126338

RetQESUBoptemp = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n = 1, M = 1e7, optimized = T, target = "skew", mode = "conditional", v = v)[,2] - log(100)

c(mean(RetQESUBoptemp), mean(RetQESUBoptemp^2), mean(RetQESUBoptemp^3), mean(RetQESUBoptemp^4), var(RetQESUBoptemp), skewness(RetQESUBoptemp), kurtosis(RetQESUBoptemp, excess = F))
# [1]   1.1827406  1.7951524  2.9826491  5.4267152  0.3962771 -0.3124229  3.2637484

# Do scaling

mu1 =  MomentsBates(mu, kappa, theta, sigma, rho, lambda = 0, muj = 0, vj = 0, t = T, v0 = v0)[1]
mu1.disc =  CondMomQE("mu1", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.8310602, 1-0.8310602)
zeta1 = MomentsBates(mu, kappa, theta, sigma, rho, lambda = 0, muj = 0, vj = 0, t = T, v0 = v0)[5]
zeta1.disc = CondMomQE("var", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.8310602, 1-0.8310602)
b = sqrt(zeta1 / zeta1.disc)
a = mu1 * (1 - b)

RetQESUBoptemp.scaled = a + b * RetQESUBoptemp
c(mean(RetQESUBoptemp.scaled), var(RetQESUBoptemp.scaled), skewness(RetQESUBoptemp.scaled), kurtosis(RetQESUBoptemp.scaled, excess = F))
#  1.8282211  0.7957332 -0.3124229  3.2637484 # Variance is working, mean is not

a.new = mu1 - mu1.disc * b
RetQESUBoptemp.scaled.new = a.new + b * RetQESUBoptemp
c(mean(RetQESUBoptemp.scaled.new), var(RetQESUBoptemp.scaled.new), skewness(RetQESUBoptemp.scaled.new), kurtosis(RetQESUBoptemp.scaled.new, excess = F))
# [1] -0.3650811  0.7957332 -0.3124229  3.2637484
# Which is (still) absolutely perfect!

# Do the same for T  = 4 and ugly parameters

T = 4
kappa = 0.5
sigma = 1

RetQEemp = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n = 1, M = 1e7, .5, .5, optimized = F)[,2] - log(100)

# Do the return simulation with exact simulation of the variance
v = Exact(v0, kappa, theta, sigma, T, n, M)

RetQESUBemp = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n , M, .5, .5, optimized = F, v = v)[,2] - log(100)

MomQETheory = c(
  CondMomQE("mu1", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("mu2", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("mu3", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("mu4", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("var", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("skew", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("kurt", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5)
)

ResQE = cbind(MomQETheory, 
              c(mean(RetQEemp), mean(RetQEemp^2), mean(RetQEemp^3), mean(RetQEemp^4), var(RetQEemp), skewness(RetQEemp), kurtosis(RetQEemp, excess = F)),
              c(mean(RetQESUBemp), mean(RetQESUBemp^2), mean(RetQESUBemp^3), mean(RetQESUBemp^4), var(RetQESUBemp), skewness(RetQESUBemp), kurtosis(RetQESUBemp, excess = F)),
              MomentsBates(mu, kappa, theta, sigma, rho, lambda = 0, muj = 0, vj = 0, t = T, v0 = v0)
)
# > ResQE
# MomQETheory                                 
# Ert.v0     -0.2607676 -0.2604465 -0.2610872 -0.3021802
# Ert2.v0     1.3044053  1.3049130  1.3057614  1.3600981
# Ert3.v0    -6.4352825 -5.4379597 -6.4426656 -6.2856300
# Ert4.v0    50.2045687 33.3353347 50.2870413 45.3815161
# Varrt.v0    1.2364056  1.2370807  1.2375950  1.2687853
# Skewrt.v0  -3.9644186 -3.2368703 -3.9624835 -3.5739988
# Kurtrt.v0  28.7894766 18.4187318 28.7787960 23.9183216


get.opt.gamma1(Delta = T, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, vt = v0, target = "skew", mode = "conditional")
# [1] 0.629

print(CondMomQE("skew", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.629, 1-0.629))
# [1] -3.573974

RetQESUBoptemp = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n = 1, M = 1e7, optimized = T, target = "skew", mode = "conditional", v = v)[,2] - log(100)

c(mean(RetQESUBoptemp), mean(RetQESUBoptemp^2), mean(RetQESUBoptemp^3), mean(RetQESUBoptemp^4), var(RetQESUBoptemp), skewness(RetQESUBoptemp), kurtosis(RetQESUBoptemp, excess = F))
# [1]   -0.2266990  0.9510105 -3.6674934 23.5176818  0.8996182 -3.5674642 25.3021466

# Do scaling

mu1 =  MomentsBates(mu, kappa, theta, sigma, rho, lambda = 0, muj = 0, vj = 0, t = T, v0 = v0)[1]
mu1.disc =  CondMomQE("mu1", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.629, 1-0.629)
zeta1 = MomentsBates(mu, kappa, theta, sigma, rho, lambda = 0, muj = 0, vj = 0, t = T, v0 = v0)[5]
zeta1.disc = CondMomQE("var", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.629, 1-0.629)
b = sqrt(zeta1 / zeta1.disc)
a = mu1 * (1 - b)

RetQESUBoptemp.scaled = a + b * RetQESUBoptemp
c(mean(RetQESUBoptemp.scaled), var(RetQESUBoptemp.scaled), skewness(RetQESUBoptemp.scaled), kurtosis(RetQESUBoptemp.scaled, excess = F))
#  -0.2125301  1.2690603 -3.5674642 25.3021466

a.new = mu1 - mu1.disc * b
RetQESUBoptemp.scaled.new = a.new + b * RetQESUBoptemp
c(mean(RetQESUBoptemp.scaled.new), var(RetQESUBoptemp.scaled.new), skewness(RetQESUBoptemp.scaled.new), kurtosis(RetQESUBoptemp.scaled.new, excess = F))
# [1] -0.3022553  1.2690603 -3.5674642 25.3021466
# Here, the Kurtosis is not so well fitting anymore but everything else still works great


# Try for PRICES instead as returns:

XQEemp = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n = 1, M = 1e7, .5, .5, optimized = F)[,2] 

# Do the return simulation with exact simulation of the variance
v = Exact(v0, kappa, theta, sigma, T, n  = 1, M = 1e7)

XQESUBemp = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n = 1, M = 1e7, .5, .5, optimized = F, v = v)[,2]

MomQETheory = c(
  CondMomQE("mu1", Xt = log(S0), vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5, R = FALSE),
  CondMomQE("var", Xt = log(S0), vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("skew", Xt = log(S0), vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5),
  CondMomQE("kurt", Xt = log(S0), vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.5, 0.5)
)

ResQE = cbind(MomQETheory, 
              c(mean(XQEemp), var(XQEemp), skewness(XQEemp), kurtosis(XQEemp, excess = F)),
              c(mean(XQESUBemp), var(XQESUBemp), skewness(XQESUBemp), kurtosis(XQESUBemp, excess = F)),
              c(MomentsBates(mu, kappa, theta, sigma, rho, lambda = 0, muj = 0, vj = 0, t = T, v0 = v0)[1] + log(S0), MomentsBates(mu, kappa, theta, sigma, rho, lambda = 0, muj = 0, vj = 0, t = T, v0 = v0)[5:7])
)
# > ResQE
# MomQETheory                              
# Ert.v0       4.344403  4.344150  4.344681  4.302990
# Varrt.v0     1.236406  1.237051  1.234728  1.268785
# Skewrt.v0   -3.964419 -3.230789 -3.962354 -3.573999
# Kurtrt.v0   28.789477 18.324891 28.826002 23.918322


get.opt.gamma1(Delta = T, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, vt = v0, target = "skew", mode = "conditional")
# [1] 0.629

print(CondMomQE("skew", vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.629, 1-0.629))
# [1] -3.573974

XQESUBoptemp = QEScheme(S0, v0, kappa, theta, sigma, mu, rho, T, n = 1, M = 1e7, optimized = T, target = "skew", mode = "conditional", v = v)[,2]

c(mean(XQESUBoptemp), var(XQESUBoptemp), skewness(XQESUBoptemp), kurtosis(XQESUBoptemp, excess = F))
# [1]  4.3788149  0.8967231 -3.5625953 25.3605558  

# Do scaling

mu1 =  MomentsBates(mu, kappa, theta, sigma, rho, lambda = 0, muj = 0, vj = 0, t = T, v0 = v0)[1] + log(100)
mu1.disc =  CondMomQE("mu1", Xt = log(S0), vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.629, 1-0.629, R = F)
zeta1 = MomentsBates(mu, kappa, theta, sigma, rho, lambda = 0, muj = 0, vj = 0, t = T, v0 = v0)[5]
zeta1.disc = CondMomQE("var", Xt = log(100), vt = v0, kappa, theta, sigma, rho, mu, Delta = T, 0.629, 1-0.629)
b = sqrt(zeta1 / zeta1.disc)
a = mu1 * (1 - b)

XQESUBoptemp.scaled = a + b * XQESUBoptemp
c(mean(XQESUBoptemp.scaled), var(XQESUBoptemp.scaled), skewness(XQESUBoptemp.scaled), kurtosis(XQESUBoptemp.scaled, excess = F))
#  4.393048  1.264976 -3.562595 25.360556

a.new = mu1 - mu1.disc * b
XQESUBoptemp.scaled.new = a.new + b * XQESUBoptemp
c(mean(XQESUBoptemp.scaled.new), var(XQESUBoptemp.scaled.new), skewness(XQESUBoptemp.scaled.new), kurtosis(XQESUBoptemp.scaled.new, excess = F))
# [1] 4.303323  1.264976 -3.562595 25.360556
# Here, the Kurtosis is not so well fitting anymore but everything else still works great


###################### Directly with novel QE #########

v = Exact(v0, kappa, theta, sigma, T, n, M = 1e7)
XQESUBoptemp.scaled.new.direct = QESchemeV2(S0, v0, kappa, theta, sigma, mu, rho, T, n, M = 1e7, optimized = T, target = "skew", mode = "conditional", v = v)[,n+1]

c(mean(XQESUBoptemp.scaled.new.direct),  var(XQESUBoptemp.scaled.new.direct), skewness(XQESUBoptemp.scaled.new.direct), kurtosis(XQESUBoptemp.scaled.new.direct, excess = F))
# 4.302920  1.268995 -3.574287 25.459549
RQESUBoptemp.scaled.new.direct = XQESUBoptemp.scaled.new.direct - log(S0)
c(mean(RQESUBoptemp.scaled.new.direct),  var(RQESUBoptemp.scaled.new.direct), skewness(RQESUBoptemp.scaled.new.direct), kurtosis(RQESUBoptemp.scaled.new.direct, excess = F))
# -0.3022498  1.2689950 -3.5742874 25.4595491
MomentsBates(mu, kappa, theta, sigma, rho, lambda = 0, muj= 0, vj = 0, t = T, v0)
# Ert.v0    Ert2.v0    Ert3.v0    Ert4.v0   Varrt.v0  Skewrt.v0  Kurtrt.v0 
# -0.3021802  1.3600981 -6.2856300 45.3815161  1.2687853 -3.5739988 23.9183216 

# do more than just 1 step
v = Exact(v0, kappa, theta, sigma, T, n = 2, M = 1e6)
XQESUBoptemp.scaled.new.direct = QESchemeV2(S0, v0, kappa, theta, sigma, mu, rho, T, n = 2, M = 1e6, optimized = T, target = "skew", mode = "conditional", v = v)[,2+1]

c(mean(XQESUBoptemp.scaled.new.direct),  var(XQESUBoptemp.scaled.new.direct), skewness(XQESUBoptemp.scaled.new.direct), kurtosis(XQESUBoptemp.scaled.new.direct, excess = F))
# 4.303214  1.282292 -3.667577 25.287383

RQESUBoptemp.scaled.new.direct = XQESUBoptemp.scaled.new.direct - log(S0)
c(mean(RQESUBoptemp.scaled.new.direct),  var(RQESUBoptemp.scaled.new.direct), skewness(RQESUBoptemp.scaled.new.direct), kurtosis(RQESUBoptemp.scaled.new.direct, excess = F))
# -0.3019557  1.2822917 -3.6675772 25.2873829

XQEemp = QESchemeV2(S0, v0, kappa, theta, sigma, mu, rho, T, n = 2, M = 1e6, optimized = F, gamma1 = 0.5, gamma2 = 0.5)[,2+1]
RQEemp = XQEemp - log(S0)
c(mean(RQEemp),  var(RQEemp), skewness(RQEemp), kurtosis(RQEemp, excess = F))
# -0.2908899  1.2861815 -3.2468834 19.6147666

############################### Proper?!? 08.10.20 ######

test = ORS(100, 0.1, 0.5, 0.19, 1, 0, -0.7, T = 4, n = 2, M = 5e6, target = "skew", mode = "conditional")
MomentsBates(0, 0.5, 0.19, 1, -0.7, 0, 0, 0, 2, 0.1)
# -0.1331091  0.4537311 -1.1251123  4.4629514  0.4360130 -3.2949825 20.5736072
ORSemp = test[,2] - log(100)
c(mean(ORSemp),  var(ORSemp), skewness(ORSemp), kurtosis(ORSemp, excess = F)) # so this "works"
# -0.1334196  0.4368516 -3.3038563 22.3763097

MomentsBates(0, 0.5, 0.19, 1, -0.7, 0, 0, 0, 4, 0.1)
# -0.3021802  1.3600981 -6.2856300 45.3815161  1.2687853 -3.5739988 23.9183216 
ORSemp = test[,3] - log(100)
c(mean(ORSemp),  var(ORSemp), skewness(ORSemp), kurtosis(ORSemp, excess = F)) # so this "works"
# -0.3827038  1.8088464 -3.5643183 23.9550698

# mean and variance does not go very well. But this is just a matter of scaling and shaping!!! The optimization part is working great here.
# After the second iteration, the mean(lnSt_tilde) and var(lnSt_tilde) does not match with mu1.disc and zeta1.disc anymore. Therefore, the scaling/shaping fails!
# I think this is because we use lnSt[,2], which is already perfect and not biased?

test = Exact(0.1, 0.5, 0.19, 1, 4, 2, 1e6)
MomentsCIR(p = 1:4, kappa  = 0.5, theta  = 0.19, sigma = 1, v0 = 0.1, t = 2)
# Evt^1|v0  Evt^2|v0  Evt^3|v0  Evt^4|v0 
# 0.1568909 0.1470431 0.2456636 0.5842566 

# > mean(test[,2])
# [1] 0.1568575
# > mean(test[,2]^2)
# [1] 0.1471393
# > mean(test[,2]^3)
# [1] 0.245693
# > mean(test[,2]^4)
# [1] 0.5804882

MomentsCIR(p = 1:4, kappa  = 0.5, theta  = 0.19, sigma = 1, v0 = 0.1, t = 4)
# Evt^1|v0  Evt^2|v0  Evt^3|v0  Evt^4|v0 
# 0.1778198 0.1970764 0.4002523 1.1796324 

# > mean(test[,3])
# [1] 0.1782068
# > mean(test[,3]^2)
# [1] 0.1985944
# > mean(test[,3]^3)
# [1] 0.4085799
# > mean(test[,3]^4)
# [1] 1.22892