library(doParallel)
library(ggplot2)

source("SupplementCodeV2.r")

source("Alfonsi.r")
source("DeelstraDelbaen.r")
source("EulerCIR.r")
source("Zhu.r")
source("HighamMao.r")
source("Milstein.r")
source("KJ.r")
source("Glasserman.r")
source("AndersenBrothertonRatcliffe.r")
source("QE.r")
source("Sankaran.r")
source("TseWanIPZ.r")
source("Exact.r")

source("getSkFromMoments.r")
source("getKuFromMoments.r")

source("TreatVariance.r")
# source("MomError.r")

source("getMErrors.r")

###### Generate MAPE plots for a high noncentrality (HNC) parameter setting and the conditional (option pricing way of thinking) simulation setup

# n = c(252, 50, 12, 4, 2, 1)
n = c(12, 4, 2, 1)
T = 1
M = 1e5
NRuns = 100

kappaHNC = 3
theta = 0.19
sigmaHNC = 0.4
v0HNC = 0.1

names = factor(c("AE0", "DD", "EAF", "ERF", "EFT", "TVS", "HM", "KJ", "G", "ABR", "QE", "SSA", "SNA", "BK"),  levels =c("AE0", "DD", "EAF", "ERF", "EFT", "TVS", "HM", "KJ", "G", "ABR", "QE", "SSA", "SNA", "BK"))

nE = 7
nI = 2
nO = 5

df.MErrorsHNC.C = getCIRMErrors(v0 = v0HNC, kappa = kappaHNC, theta = theta, sigma = sigmaHNC, T = T, n = n, M = M, NRuns = NRuns, names = names, nE = nE, nI = nI, nO = nO, setup = 'C', approach = 0)

save(df.MErrorsHNC.C, file = "df.MErrorsHNC.C.RData")


# Define some stuff to generate nice looking graphics

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colors = gg_color_hue(3) # Mimic the colors, ggplot would choose, if we mapped 3 variables to colors
shapes =  c(seq(0, length(names) -1))


# Generate the MAPE plots wo and w confidence bounds of 90 %

gg.MErrorsHNC.C = ggplot(df.MErrorsHNC.C, aes(x = Delta, y = MAPE, color = Group, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") +  guides(color = "none", shape = guide_legend(override.aes = list(color = rep(colors, c(7,2,5))))) + theme_bw() +  scale_x_reverse() + xlab(expression(Delta))
gg.MErrorsHNC.C.f = ggplot(subset(df.MErrorsHNC.C, Scheme != 'AE0'), aes(x = Delta, y = MAPE, color = Group, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes[-1]) + facet_grid(.~Moment, scales = "free") +  guides(color = "none", shape = guide_legend(override.aes = list(color = rep(colors, c(6,2,5))))) + theme_bw() +  scale_x_reverse() + xlab(expression(Delta))
gg.MErrorsHNC.C.ff = ggplot(subset(df.MErrorsHNC.C, Group == 'Others'), aes(x = Delta, y = MAPE, color = Scheme, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes[c(10:14)]) + facet_grid(.~Moment, scales = "free") + theme_bw() +  scale_x_reverse() +  xlab(expression(Delta))
gg.MErrorsHNC.C.fff = ggplot(subset(df.MErrorsHNC.C, Scheme %in% c('QE', 'SNA', 'Exact')), aes(x = Delta, y = MAPE, color = Scheme, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes[c(11, 13:14)]) + facet_grid(.~Moment, scales = "free") + theme_bw() +  scale_x_reverse() +  xlab(expression(Delta))

gg.MErrorsHNC.C.cb = gg.MErrorsHNC.C + geom_ribbon(aes(ymin = q05, ymax = q95,  fill = Scheme), alpha = 0.1, colour = NA)
gg.MErrorsHNC.C.f.cb = gg.MErrorsHNC.C.f + geom_ribbon(aes(ymin = q05, ymax = q95, fill = Scheme), alpha = 0.1, colour = NA)
gg.MErrorsHNC.C.ff.cb = gg.MErrorsHNC.C.ff + geom_ribbon(aes(ymin = q05, ymax = q95, fill = Scheme), alpha = 0.1, colour = NA)
gg.MErrorsHNC.C.fff.cb = gg.MErrorsHNC.C.fff + geom_ribbon(aes(ymin = q05, ymax = q95, fill = Scheme), alpha = 0.1, colour = NA)

# Plots which perfectly blend in with other paper: time against bias, exactly what we want
pap1HNC.C = ggplot(df.MErrorsHNC.C, aes(x = Mean.Time, y = MAE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() 

pap1HNC.C.cb = pap1 + geom_ribbon(aes(ymin = q05AE, ymax = q95AE, fill = Scheme), alpha = 0.1, colour = NA)

pap1.2HNC.C = ggplot(df.MErrorsHNC.C, aes(x = log10(Mean.Time), y = log10(MAE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab("log10(Mean Time)") + ylab("log10(MAE)")

pap1.2HNC.C.cb = pap1.2HNC.C + geom_ribbon(aes(ymin = log10(q05AE), ymax = log10(q95AE), fill = Scheme), alpha = 0.1, colour = NA)

pap2HNC.C = ggplot(df.MErrorsHNC.C, aes(x = Mean.Time, y = RMSE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() 

pap2.2HNC.C = ggplot(df.MErrorsHNC.C, aes(x = log10(Mean.Time), y = log10(RMSE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab("log10(Mean Time)") + ylab("log10(RMSE)")

pap2.2HNC.C.cb = pap2.2HNC.C + geom_ribbon(aes(ymin = log10(q05RSE), ymax = log10(q95RSE), fill = Scheme), alpha = 0.1, colour = NA)

pap3HNC.C = ggplot(df.MErrorsHNC.C, aes(x = Mean.Time, y = MAPE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() 

pap3.2HNC.C = ggplot(df.MErrorsHNC.C, aes(x = log10(Mean.Time), y = log10(MAPE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab("log10(Mean Time)") + ylab("log10(MAPE)")


pdf("MErrorsHNC.C.pdf", width = 16, height = 9)
plot(pap1.2HNC.C)
plot(pap2.2HNC.C)
plot(pap3.2HNC.C)
dev.off()


###### Generate MAPE plots for a high noncentrality (HNC) parameter setting and the unconditional simulation setup

names = factor(c("AE0", "DD", "EulerAF", "EulerRF", "EulerFT", "TVS", "HM", "KJ", "G", "ABR", "QE", "SSA", "SNA", "Exact"),  levels = c("AE0", "DD", "Euler", "EulerAF", "EulerRF", "EulerFT", "TVS", "HM", "M", "KJ", "G", "ABR", "QE", "SSA", "SNA", "Exact"))

nE = 7
nI = 2
nO = 5

# Calculate shape and scale of steady state distribution of v

shape = 2 * kappa * theta / sigma^2
scale = theta / shape 

n = 1
T = c(252, 50, 12, 4, 2, 1)^-1 * n

v0 = rgamma(M, shape = shape, scale = scale)
df.MErrorsHNC.UC = getCIRMErrors(v0 = v0, kappa = kappa, theta = theta, sigma = sigma, T = T, n = n, M = M, NRuns = NRuns, names = names, nE = nE, nI = nI, nO = nO, setup = 'UC', approach = 1)

gg.MAEHNC.UC = ggplot(df.ErrorsHNC.UC, aes(x = Delta, y = MAE, color = Group, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") +  guides(color = "none", shape = guide_legend(override.aes = list(color = rep(colors, c(7,2,5))))) + theme_bw() +  scale_x_reverse() + xlab(expression(Delta))
gg.MAEHNC.UC.f = ggplot(subset(df.ErrorsHNC.UC, Scheme != 'AE0'), aes(x = Delta, y = MAE, color = Group, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes[-1]) + facet_grid(.~Moment, scales = "free") +  guides(color = "none", shape = guide_legend(override.aes = list(color = rep(colors, c(6,2,5))))) + theme_bw() +  scale_x_reverse() + xlab(expression(Delta))
gg.MAEHNC.UC.ff = ggplot(subset(df.ErrorsHNC.UC, Group == 'Others'), aes(x = Delta, y = MAE, color = Scheme, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes[c(10:14)]) + facet_grid(.~Moment, scales = "free") + theme_bw() +  scale_x_reverse() +  xlab(expression(Delta))
gg.MAEHNC.UC.fff = ggplot(subset(df.ErrorsHNC.UC, Scheme %in% c('QE', 'SNA', 'Exact')), aes(x = Delta, y = MAE, color = Scheme, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes[c(11, 13:14)]) + facet_grid(.~Moment, scales = "free") + theme_bw() +  scale_x_reverse() +  xlab(expression(Delta))

gg.MAEHNC.UC.cb = gg.MAEHNC.UC + geom_ribbon(aes(ymin = q05AE, ymax = q95AE,  fill = Scheme), alpha = 0.1, colour = NA)
gg.MAEHNC.UC.f.cb = gg.MAEHNC.UC.f + geom_ribbon(aes(ymin = q05AE, ymax = q95AE, fill = Scheme), alpha = 0.1, colour = NA)
gg.MAEHNC.UC.ff.cb = gg.MAEHNC.UC.ff + geom_ribbon(aes(ymin = q05AE, ymax = q95AE, fill = Scheme), alpha = 0.1, colour = NA)
gg.MAEHNC.UC.fff.cb = gg.MAEHNC.UC.fff + geom_ribbon(aes(ymin = q05AE, ymax = q95AE, fill = Scheme), alpha = 0.1, colour = NA)


# such a plot does not make much sense here, because the time is always quite the same since n is equal
# pap = ggplot(subset(df.ErrorsHNC.UC, Moment == "Skewness"), aes(x = Mean.Time, y = MAE, color = Group, shape = Scheme)) + 
#  geom_point() + geom_line() + scale_shape_manual(values = shapes) +  
#  guides(color = "none", shape = guide_legend(override.aes = list(color = rep(colors, c(7,2,5))))) + theme_bw() + xlab("Mean Time (in sec)")

pdf("MAEHNCsetup1UC.pdf", width = 16, height = 9)
plot(gg.MAEHNC.UC)
plot(gg.MAEHNC.UC.cb)
plot(gg.MAEHNC.UC.f)
plot(gg.MAEHNC.UC.f.cb)
plot(gg.MAEHNC.UC.ff)
plot(gg.MAEHNC.UC.ff.cb)
plot(gg.MAEHNC.UC.fff)
plot(gg.MAEHNC.UC.fff.cb)
dev.off()


###### Generate MAPE plots for a low noncentrality (LNC) parameter setting and the conditional (option pricing way of thinking) simulation setup

n = c(12, 4, 2, 1)
T = 1
M = 1e5
NRuns = 100


kappaLNC = 0.5
theta = 0.19
sigmaLNC = 1
v0LNC = 0.05

names = factor(c("AE0", "DD", "EAF", "ERF", "EFT", "TVS", "HM", "KJ", "G", "ABR", "QE", "BK"),  levels =c("AE0", "DD", "EAF", "ERF", "EFT", "TVS", "HM", "KJ", "G", "ABR", "QE", "BK"))

nE = 7
nI = 2
nO = 3


df.MErrorsLNC.C = getCIRMErrors(v0 = v0LNC, kappa = kappaLNC, theta = theta, sigma = sigmaLNC, T = T, n = n, M = M, NRuns = NRuns, names = names, nE = nE, nI = nI, nO = nO, setup = 'C', approach = 0)

save(df.MErrorsLNC.C, file = "df.MErrorsLNC.C.RData")

# Define some stuff to generate nice looking graphics

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colors = gg_color_hue(3) # Mimic the colors, ggplot would choose, if we mapped 3 variables to colors
shapes =  c(seq(0, length(names) -1))


# Generate the MAPE plots w and wo confidence bounds of 90 %

gg.MErrorsLNC.C = ggplot(df.MErrorsLNC.C, aes(x = Delta, y = MAPE, color = Group, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") +  guides(color = "none", shape = guide_legend(override.aes = list(color = rep(colors, c(7,2,5))))) + theme_bw() +  scale_x_reverse() + xlab(expression(Delta))
gg.MErrorsLNC.C.f = ggplot(subset(df.MErrorsLNC.C, Scheme != 'AE0'), aes(x = Delta, y = MAPE, color = Group, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes[-1]) + facet_grid(.~Moment, scales = "free") +  guides(color = "none", shape = guide_legend(override.aes = list(color = rep(colors, c(6,2,5))))) + theme_bw() +  scale_x_reverse() + xlab(expression(Delta))
gg.MErrorsLNC.C.ff = ggplot(subset(df.MErrorsLNC.C, Group == 'Others'), aes(x = Delta, y = MAPE, color = Scheme, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes[c(10:14)]) + facet_grid(.~Moment, scales = "free") + theme_bw() +  scale_x_reverse() +  xlab(expression(Delta))
gg.MErrorsLNC.C.fff = ggplot(subset(df.MErrorsLNC.C, Scheme %in% c('QE', 'SNA', 'Exact')), aes(x = Delta, y = MAPE, color = Scheme, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes[c(11, 13:14)]) + facet_grid(.~Moment, scales = "free") + theme_bw() +  scale_x_reverse() +  xlab(expression(Delta))

gg.MErrorsLNC.C.cb = gg.MErrorsLNC.C + geom_ribbon(aes(ymin = q05, ymax = q95,  fill = Scheme), alpha = 0.1, colour = NA)
gg.MErrorsLNC.C.f.cb = gg.MErrorsLNC.C.f + geom_ribbon(aes(ymin = q05, ymax = q95, fill = Scheme), alpha = 0.1, colour = NA)
gg.MErrorsLNC.C.ff.cb = gg.MErrorsLNC.C.ff + geom_ribbon(aes(ymin = q05, ymax = q95, fill = Scheme), alpha = 0.1, colour = NA)
gg.MErrorsLNC.C.fff.cb = gg.MErrorsLNC.C.fff + geom_ribbon(aes(ymin = q05, ymax = q95, fill = Scheme), alpha = 0.1, colour = NA)


# Plots which perfectly blend in with other paper: time against bias, exactly what we want

pap1LNC.C = ggplot(df.MErrorsLNC.C, aes(x = Mean.Time, y = MAE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw()

pap1LNC.C.cb = pap1LNC.C + geom_ribbon(aes(ymin = q05AE, ymax = q95AE, fill = Scheme), alpha = 0.1, colour = NA)

pap1.2LNC.C = ggplot(df.MErrorsLNC.C, aes(x = log10(Mean.Time), y = log10(MAE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab("log10(Mean Time)") + ylab("log10(MAE)")

pap1.2LNC.C.cb = pap1.2 + geom_ribbon(aes(ymin = log10(q05AE), ymax = log10(q95AE), fill = Scheme), alpha = 0.1, colour = NA)

pap2LNC.C = ggplot(df.MErrorsLNC.C, aes(x = Mean.Time, y = RMSE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() 

pap2.2LNC.C = ggplot(df.MErrorsLNC.C, aes(x = log10(Mean.Time), y = log10(RMSE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab("log10(Mean Time)") + ylab("log10(RMSE)")

pap2.2LNC.C.cb = pap2.2 + geom_ribbon(aes(ymin = log10(q05RSE), ymax = log10(q95RSE), fill = Scheme), alpha = 0.1, colour = NA)

pap3LNC.C = ggplot(df.MErrorsLNC.C, aes(x = Mean.Time, y = MAPE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() 

pap3.2LNC.C = ggplot(df.MErrorsLNC.C, aes(x = log10(Mean.Time), y = log10(MAPE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab("log10(Mean Time)") + ylab("log10(MAPE)")

pdf("MErrorsLNC.C.pdf", width = 16, height = 9)
plot(pap1.2LNC.C)
plot(pap2.2LNC.C)
plot(pap3.2LNC.C)
dev.off()



###### Generate MAPE plots for a low noncentrality (LNC) parameter setting and the unconditional simulation setup

# I left out Euler, Milstein, Glasserman, SSA and SNA here, for G we have huge errors and for SSA and SNA NAs
# Leave in Euler and Milstein!

names = factor(c("AE0", "DD", "EulerAF", "EulerRF", "EulerFT", "TVS", "HM", "KJ", "ABR", "QE", "Exact"),  levels = c("AE0", "DD", "EulerAF", "EulerRF", "EulerFT", "TVS", "HM", "KJ", "ABR", "QE", "Exact"))

shapes =  c(seq(0, length(names) - 1))

# Calculate shape and scale of steady state distribution of v

shape = 2 * kappa * theta / sigma^2
scale = theta / shape 

n = 1
T = c(252, 50, 12, 4, 2, 1)^-1 * n

v0 = rgamma(M, shape = shape, scale = scale)
df.MAPELNC.UC = getCIRMAPE(v0 = v0, kappa = kappa, theta = theta, sigma = sigma, T = T, n = n, M = M, NRuns = NRuns, names = names, nE = 7, nI = 1, nO = 3, setup = 'UC', approach = 1)

gg.MAPELNC.UC = ggplot(df.MAPELNC.UC, aes(x = Delta, y = MAPE, color = Group, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") +  guides(color = "none", shape = guide_legend(override.aes = list(color = rep(colors, c(7,1,3))))) + theme_bw() +  scale_x_reverse() + xlab(expression(Delta))
gg.MAPELNC.UC.f = ggplot(subset(df.MAPELNC.UC, Group == 'Others'), aes(x = Delta, y = MAPE, color = Scheme, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes[c(8:10)]) + facet_grid(.~Moment, scales = "free") + theme_bw() +  scale_x_reverse() +  xlab(expression(Delta))
gg.MAPELNC.UC.ff = ggplot(subset(df.MAPELNC.UC, Scheme %in% c('QE', 'Exact')), aes(x = Delta, y = MAPE, color = Scheme, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes[c(9:10)]) + facet_grid(.~Moment, scales = "free") + theme_bw() +  scale_x_reverse() +  xlab(expression(Delta))


gg.MAPELNC.UC.cb = gg.MAPELNC.UC + geom_ribbon(aes(ymin = q05, ymax = q95,  fill = Scheme), alpha = 0.1, colour = NA)
gg.MAPELNC.UC.f.cb = gg.MAPELNC.UC.f + geom_ribbon(aes(ymin = q05, ymax = q95, fill = Scheme), alpha = 0.1, colour = NA)
gg.MAPELNC.UC.ff.cb = gg.MAPELNC.UC.ff + geom_ribbon(aes(ymin = q05, ymax = q95, fill = Scheme), alpha = 0.1, colour = NA)


pdf("MAPELNCsetup1UC.pdf", width = 16, height = 9)
plot(gg.MAPELNC.UC)
plot(gg.MAPELNC.UC.cb)
plot(gg.MAPELNC.UC.f)
plot(gg.MAPELNC.UC.f.cb)
plot(gg.MAPELNC.UC.ff)
plot(gg.MAPELNC.UC.ff.cb)
dev.off()