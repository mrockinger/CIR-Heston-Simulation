# Test skewness estimator

library(EnvStats)
library(ggplot2)
library(yarrr)

source("SupplementCodeV2.R")
source("ORSproper.R")


M = 10000
T = 100
Delta = 1 / 252
n = T / Delta

######################################################################## HNC

S0 = 100
mu = 0
rho = -0.7
kappaHNC = 3
theta = 0.19
sigmaHNC = 0.4
v0HNC = 0.19

XHNC = ORS(S0, v0HNC, kappaHNC, theta, sigmaHNC, mu, rho, T, n, M, target = "skew", mode = "unconditional")
RHNC = apply(XHNC, 1, diff)

momTHNC = MomentsBates(mu, kappaHNC, theta, sigmaHNC, rho, 0, 0, 0, Delta, 1, F)
skewTHNC = momTHNC[5]

 # Now start to shrink the sample size 

Sub = c(1, 0.75, 0.5, 0.25, 0.1, 0.05, 0.01)
n.sub = n * Sub
T.sub = T * Sub 

skewEstHNC = matrix(NA, nrow = M, ncol = length(Sub))
skewEstHNC.r.m.s.e = skewEstHNC.m.a.e = skewEstHNC.m.a.p.e = c()

for(i in 1:length(Sub)){
  
  R.sub = RHNC[1:(Sub[i] * n),]
  skewEstHNC[,i] = apply(R.sub, 2, skewness)
  skewEstHNC.r.m.s.e[i] = sqrt(mean((skewEstHNC[,i] - skewTHNC)^2)) 
  skewEstHNC.m.a.e[i] = mean(abs(skewEstHNC[,i] - skewTHNC)) #MAE
  skewEstHNC.m.a.p.e[i] = mean( abs( (skewEstHNC[,i] - skewTHNC) / skewTHNC ) ) #MAPE
  
}

# save(RHNC, skewEstHNC, file = "ResHNCskewEst.RData")
load("C:/Users/mschmid/Desktop/ResHNCskewEst.RData")

# Generate pirate pltots of estimators density

dfHNCEst = data.frame("Observations" = c(skewEstHNC), T =  rep(T.sub, each = M))

pdf("SkewEstHNCPirate.pdf")
pirateplot(formula = Observations ~ T,
           data = dfHNCEst,
           theme = 1,
           main = "Distribution of Skewness Estimators",
           avg.line.o = 1,
           bean.b.o = 0.4,
           point.o = 0,
           quant = c(0.1,0.9),
           quant.col = "black",
           quant.lwd = 1,
           inf.f.o = 0,
           inf.b.o = 0,
           inf.method = 'ci',
           avg.line.col = "black",
           avg.line.fun = mean
           )
abline(h = skewTHNC, col = "red3", lty = 2, lwd = 2)
dev.off()

# Plot mean absolute (percentage) bias (Y) against sample size in T (X)

dfHNC = data.frame("MAE" = skewEstHNC.m.a.e, "RMSE" = skewEstHNC.r.m.s.e, "MAPE" = skewEstHNC.m.a.p.e, T = T.sub, n = n.sub)

ggHNC1 = ggplot(dfHNC, aes(y = MAE, x = T)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggHNC1.2 = ggplot(dfHNC, aes(y = MAE, x = n)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggHNC2 = ggplot(dfHNC, aes(y = log10(MAE), x = log10(T))) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggHNC2.2 = ggplot(dfHNC, aes(y = log10(MAE), x = log10(n))) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggHNC3 = ggplot(dfHNC, aes(y = MAPE, x = T)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggHNC3.2 = ggplot(dfHNC, aes(y = MAPE, x = n)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggHNC4 = ggplot(dfHNC, aes(y = log10(MAPE), x = log10(T))) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggHNC4.2 = ggplot(dfHNC, aes(y = log10(MAPE), x = log10(n))) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()

pdf("skewEstHNC.pdf")
plot(ggHNC1)
plot(ggHNC1.2)
plot(ggHNC2)
plot(ggHNC2.2)
plot(ggHNC3)
plot(ggHNC3.2)
plot(ggHNC4)
plot(ggHNC4.2)
dev.off()


###################################################################### LNC

kappaLNC = 0.5
sigmaLNC = 1
v0LNC = 0.05

XLNC = ORS(S0, v0LNC, kappaLNC, theta, sigmaLNC, mu, rho, T, n, M, target = "skew", mode = "unconditional")
RLNC = apply(XLNC, 1, diff)

rm(XLNC)

momTLNC = MomentsBates(mu, kappaLNC, theta, sigmaLNC, rho, 0, 0, 0, Delta, 1, F)
skewTLNC = momTLNC[5]


# Now start to shrink the sample size 
skewEstLNC = matrix(NA, nrow = M, ncol = length(Sub))
skewEstLNC.r.m.s.e = skewEstLNC.m.a.e = skewEstLNC.m.a.p.e = c()

for(i in 1:length(Sub)){
  
  R.sub = RLNC[1:(Sub[i] * n),]
  skewEstLNC[,i] = apply(R.sub, 2, skewness)
  skewEstLNC.r.m.s.e[i] = sqrt(mean((skewEstLNC[,i] - skewTLNC)^2)) 
  skewEstLNC.m.a.e[i] = mean(abs(skewEstLNC[,i] - skewTLNC)) #MAE
  skewEstLNC.m.a.p.e[i] = mean( abs( (skewEstLNC[,i] - skewTLNC) / skewTLNC ) ) #MAPE
  
}


# save(RLNC, skewEstLNC, file = "ResLNCskewEst.RData")
load("C:/Users/mschmid/Desktop/ResLNCskewEst.RData")

# Generate pirate pltots of estimators density

dfLNCEst = data.frame("Observations" = c(skewEstLNC), T =  rep(T.sub, each = M))

pdf("SkewEstLNCPirate.pdf")
pirateplot(formula = Observations ~ T,
           data = dfLNCEst,
           theme = 1,
           main = "Distribution of Skewness Estimators",
           avg.line.o = 1,
           bean.b.o = 0.4,
           point.o = 0,
           quant = c(0.1,0.9),
           quant.col = "black",
           quant.lwd = 1,
           inf.f.o = 0,
           inf.b.o = 0,
           inf.method = 'ci',
           avg.line.col = "black",
           avg.line.fun = mean
)
abline(h = skewTLNC, col = "red3", lty = 2, lwd = 2)
dev.off()


# Plot mean absolute bias (Y) against sample size in T (X)

dfLNC = data.frame("MAE" = skewEstLNC.m.a.e, "RMSE" = skewEstLNC.r.m.s.e, "MAPE" = skewEstLNC.m.a.p.e, T = T.sub, n = n.sub)

ggLNC1 = ggplot(dfLNC, aes(y = MAE, x = T)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggLNC1.2 = ggplot(dfLNC, aes(y = MAE, x = n)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggLNC2 = ggplot(dfLNC, aes(y = log10(MAE), x = log10(T))) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggLNC2.2 = ggplot(dfLNC, aes(y = log10(MAE), x = log10(n))) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggLNC3 = ggplot(dfLNC, aes(y = MAPE, x = T)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggLNC3.2 = ggplot(dfLNC, aes(y = MAPE, x = n)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggLNC4 = ggplot(dfLNC, aes(y = log10(MAPE), x = log10(T))) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggLNC4.2 = ggplot(dfLNC, aes(y = log10(MAPE), x = log10(n))) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()

pdf("skewEstLNC.pdf")
plot(ggLNC1)
plot(ggLNC1.2)
plot(ggLNC2)
plot(ggLNC2.2)
plot(ggLNC3)
plot(ggLNC3.2)
plot(ggLNC4)
plot(ggLNC4.2)
dev.off()

dfComplete = merge(dfHNC, dfLNC, all = T)
dfComplete$Parameter = rep(c("HNC", "LNC"), each = length(T.sub))

ggCom1 = ggplot(dfComplete, aes(y = MAE, x = T, color = Parameter)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggCom1.2 = ggplot(dfComplete, aes(y = MAE, x = n, color = Parameter)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggCom2 = ggplot(dfComplete, aes(y = log10(MAE), x = log10(T), color = Parameter)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggCom2.2 = ggplot(dfComplete, aes(y = log10(MAE), x = log10(n), color = Parameter)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggCom3 = ggplot(dfComplete, aes(y = MAPE, x = T, color = Parameter)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggCom3.2 = ggplot(dfComplete, aes(y = MAPE, x = n, color = Parameter)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggCom4 = ggplot(dfComplete, aes(y = log10(MAPE), x = log10(T), color = Parameter)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggCom4.2 = ggplot(dfComplete, aes(y = log10(MAPE), x = log10(n), color = Parameter)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggCom5 = ggplot(dfComplete, aes(y = RMSE, x = T, color = Parameter)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggCom5.2 = ggplot(dfComplete, aes(y = RMSE, x = n, color = Parameter)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggCom6 = ggplot(dfComplete, aes(y = log10(RMSE), x = log10(T), color = Parameter)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()
ggCom6.2 = ggplot(dfComplete, aes(y = log10(RMSE), x = log10(n), color = Parameter)) + geom_point(size = 2) + geom_line(lty = 2) + theme_bw()

pdf("skewEstComplete.pdf")
plot(ggCom1)
plot(ggCom1.2)
plot(ggCom2)
plot(ggCom2.2)
plot(ggCom3)
plot(ggCom3.2)
plot(ggCom4)
plot(ggCom4.2)
plot(ggCom5)
plot(ggCom5.2)
plot(ggCom6)
plot(ggCom6.2)
dev.off()

