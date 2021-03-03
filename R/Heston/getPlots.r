library(doParallel)
library(ggplot2)

source("SupplementCodeV2.r")

source("EulerScheme.r")
source("ZhuScheme.r")
source("KJScheme.r")
source("QEScheme.r")
source("BroadieKayaScheme.r")
source("SmithScheme.r")
source("GKScheme.r")
source("BBGScheme.r")
source("TWScheme.r")

source("../CIR/Exact.r") # for QE + SUB

source("getSkFromMoments.r")
source("getKuFromMoments.r")

source("TreatReturns.r")
source("MomErrorReturns.r")

source("getMErrors.r")

###### Generate MAPE plots for a high noncentrality (HNC) parameter setting and the conditional (option pricing way of thinking) simulation setup

n = c(12, 4, 2, 1)
T = 1
M = 1e3
NRuns = 5

mu = 0
kappa = 3
theta = 0.19
sigma = 0.4
rho = -0.7

v0 = 0.1
S0 = 100

names = factor(c("EulerFT", "Zhu", "KJ", "QE", 'BK', 'Smith', 'GK', 'TW', 'BBG', 'ORS', 'ORS + SUB'),  levels = c("EulerFT", "Zhu", "KJ", "QE", 'BK', 'Smith', 'GK', 'TW', 'BBG', 'ORS', 'ORS + SUB'))

nE = 2
nI = 1
nO = 7

df.MErrorsHNC.C = getHestonMErrors(S0 = S0, v0 = v0, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, T = T, n = n, M = M, NRuns = NRuns, names = names, nE = nE, nI = nI, nO = nO, setup = 'C', approach = 0)

# Define some stuff to generate nice looking graphics

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colors = gg_color_hue(3) # Mimic the colors, ggplot would choose, if we mapped 3 variables to colors
shapes =  c(seq(0, length(names) -1))


# Generate the MAPE plots w and wo confidence bounds of 90 %

gg.MAPEHNC = ggplot(df.MAPEHNC, aes(x = Delta, y = MAPE, colour = Scheme, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() +  scale_x_reverse() + xlab(expression(Delta))
# gg.MAPEHNC.t = ggplot(df.MAPEHNC, aes(x = etime, y = MAPE, colour = Scheme, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw()  + xlab("Time in sec")
# gg.MAPEHNC.t.f = ggplot(subset(df.MAPEHNC, !(Scheme %in% c('GK', 'BK'))), aes(x = etime, y = MAPE, colour = Scheme, shape = Scheme)) + geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw()  + xlab("Time in sec")

gg.MAPEHNC.cb = gg.MAPEHNC + geom_ribbon(aes(ymin = q05, ymax = q95,  fill = Scheme), alpha = 0.1, colour = NA)

# Plots which perfectly blend in with other paper: time against bias, exactly what we want

pap1 = ggplot(df.MErrorsHNC.C, aes(x = Mean.Time, y = MAE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw()

pap1.cb = pap1 + geom_ribbon(aes(ymin = q05AE, ymax = q95AE, fill = Scheme), alpha = 0.1, colour = NA)

pap1.2 = ggplot(df.MErrorsHNC.C, aes(x = log10(Mean.Time), y = log10(MAE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() 

pap1.2.cb = pap1.2 + geom_ribbon(aes(ymin = log10(q05AE), ymax = log10(q95AE), fill = Scheme), alpha = 0.1, colour = NA)

pap2 = ggplot(df.MErrorsHNC.C, aes(x = Mean.Time, y = RMSE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() 

pap2.2 = ggplot(df.MErrorsHNC.C, aes(x = log10(Mean.Time), y = log10(RMSE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw()

pap2.2.cb = pap2.2 + geom_ribbon(aes(ymin = log10(q05RSE), ymax = log10(q95RSE), fill = Scheme), alpha = 0.1, colour = NA)

pap3 = ggplot(df.MErrorsHNC.C, aes(x = Mean.Time, y = MAPE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() 

pap3.2 = ggplot(df.MErrorsHNC.C, aes(x = log10(Mean.Time), y = log10(MAPE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() 



# pdf("MAPEHNCsetup1C.pdf", width = 16, height = 9)
# plot(gg.MAPEHNC)
# plot(gg.MAPEHNC.cb)
# plot(gg.MAPEHNC.t)
# plot(gg.MAPEHNC.t.f)
# dev.off()


###### Generate MAPE plots for a low noncentrality (HNC) parameter setting and the conditional (option pricing way of thinking) simulation setup

