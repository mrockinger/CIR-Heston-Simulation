# Handle matlab results

library(R.matlab)
library(ggplot2)

DataHNCC = readMat("ResultsReturnsHNCC.mat")

MAEmeanHNC = DataHNCC$MAEmeanHNC
MAEvarHNC  = DataHNCC$MAEvarHNC
MAEskewHNC = DataHNCC$MAEskewHNC
MAEkurtHNC = DataHNCC$MAEkurtHNC
RMSEmeanHNC = DataHNCC$RMSEmeanHNC
RMSEvarHNC = DataHNCC$RMSEvarHNC
RMSEskewHNC = DataHNCC$RMSEskewHNC
RMSEkurtHNC = DataHNCC$RMSEkurtHNC
MAPEmeanHNC = DataHNCC$MAPEmeanHNC
MAPEvarHNC = DataHNCC$MAPEvarHNC
MAPEskewHNC = DataHNCC$MAPEskewHNC
MAPEkurtHNC = DataHNCC$MAPEkurtHNC
MetimeHNC  = DataHNCC$MetimeHNC

n = c(12, 4, 2, 1)
T = 1
Delta = T / n

nEHNC = 2
nIHNC = 1
nOHNC = 6


namesC = factor(c('EFT', 'TVS', 'KJ', 'QE',  'BK', 'S', 'GK', 'TW', 'BBG'),  levels =c('EFT', 'TVS', 'KJ', 'QE',  'BK', 'S', 'GK', 'TW', 'BBG'))
Moments = factor(c("Mean", "Variance", "Skewness", "Kurtosis"), levels = c("Mean", "Variance", "Skewness", "Kurtosis"))
groups =  factor(c("Explicit", "Implicit", "Others"), levels = c("Explicit", "Implicit", "Others")) 

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colors = gg_color_hue(3) # Mimic the colors, ggplot would choose, if we mapped 3 variables to colors
shapes =  c(seq(0, length(namesC) -1))

df.MErrorsHNC.C.matlab = data.frame("Scheme" = rep(rep(namesC, each = length(Delta)), length(Moments)), "Group" = rep(rep(groups, c(nEHNC * length(Delta), nIHNC * length(Delta), nOHNC * length(Delta))), length(Moments)), 
                        "Delta" = rep(rep(Delta, length(namesC)), length(Moments)), "Moment" = rep(Moments, each = length(namesC) * length(Delta)), 
                        "MAE" = c(c(MAEmeanHNC), c(MAEvarHNC), c(MAEskewHNC), c(MAEkurtHNC)), 
                        "RMSE" = c(c(RMSEmeanHNC), c(RMSEvarHNC), c(RMSEskewHNC), c(RMSEkurtHNC)), 
                        "MAPE" = c(c(MAPEmeanHNC), c(MAPEvarHNC), c(MAPEskewHNC), c(MAPEkurtHNC)), 
                        "Mean Time" = rep(c(MetimeHNC), length(Moments)))

save(df.MErrorsHNC.C.matlab, file = "df.MErrorsHNC.C.matlab.RData")

# Plots which perfectly blend in with other paper: time against bias, exactly what we want

pap1HNC.C = ggplot(df.MErrorsHNC.C.matlab, aes(x = Mean.Time, y = MAE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() 

pap1.2HNC.C = ggplot(df.MErrorsHNC.C.matlab, aes(x = log10(Mean.Time), y = log10(MAE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab("log10(Mean Time)") + ylab("log10(MAE)")

pap2HNC.C = ggplot(df.MErrorsHNC.C.matlab, aes(x = Mean.Time, y = RMSE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() 

pap2.2HNC.C = ggplot(df.MErrorsHNC.C.matlab, aes(x = log10(Mean.Time), y = log10(RMSE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab("log10(Mean Time)") + ylab("log10(RMSE)")

pap3HNC.C = ggplot(df.MErrorsHNC.C.matlab, aes(x = Mean.Time, y = MAPE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() 

pap3.2HNC.C = ggplot(df.MErrorsHNC.C.matlab, aes(x = log10(Mean.Time), y = log10(MAPE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab("log10(Mean Time)") + ylab("log10(MAPE)")


pdf("MErrorsHNC.C.matlab.pdf", width = 16, height = 9)
plot(pap1.2HNC.C)
plot(pap2.2HNC.C)
plot(pap3.2HNC.C)
dev.off()




############################################################################ LNC

DataLNCC = readMat("ResultsReturnsLNCC.mat")

MAEmeanLNC = DataLNCC$MAEmeanLNC
MAEvarLNC  = DataLNCC$MAEvarLNC
MAEskewLNC = DataLNCC$MAEskewLNC
MAEkurtLNC = DataLNCC$MAEkurtLNC
RMSEmeanLNC = DataLNCC$RMSEmeanLNC
RMSEvarLNC = DataLNCC$RMSEvarLNC
RMSEskewLNC = DataLNCC$RMSEskewLNC
RMSEkurtLNC = DataLNCC$RMSEkurtLNC
MAPEmeanLNC = DataLNCC$MAPEmeanLNC
MAPEvarLNC = DataLNCC$MAPEvarLNC
MAPEskewLNC = DataLNCC$MAPEskewLNC
MAPEkurtLNC = DataLNCC$MAPEkurtLNC
MetimeLNC  = DataLNCC$MetimeLNC


nELNC = 7
nILNC = 2
nOLNC = 3


namesLNC = factor(c("AE0", "DD", "EAF", "ERF", "EFT", "TVS", "HM", "KJ", "G", "ABR", "QE", "BK"),  levels =c("AE0", "DD", "EAF", "ERF", "EFT", "TVS", "HM", "KJ", "G", "ABR", "QE", "BK"))


df.MErrorsLNC.C.matlab = data.frame("Scheme" = rep(rep(namesLNC, each = length(Delta)), length(Moments)), "Group" = rep(rep(groups, c(nELNC * length(Delta), nILNC * length(Delta), nOLNC * length(Delta))), length(Moments)), 
                                    "Delta" = rep(rep(Delta, length(namesLNC)), length(Moments)), "Moment" = rep(Moments, each = length(namesLNC) * length(Delta)), 
                                    "MAE" = c(c(MAEmeanLNC), c(MAEvarLNC), c(MAEskewLNC), c(MAEkurtLNC)), 
                                    "RMSE" = c(c(RMSEmeanLNC), c(RMSEvarLNC), c(RMSEskewLNC), c(RMSEkurtLNC)), 
                                    "MAPE" = c(c(MAPEmeanLNC), c(MAPEvarLNC), c(MAPEskewLNC), c(MAPEkurtLNC)), 
                                    "Mean Time" = rep(c(MetimeLNC), length(Moments)))

save(df.MErrorsLNC.C.matlab, file = "df.MErrorsLNC.C.matlab.RData")

# Plots which perfectly blend in with other paper: time against bias, exactly what we want

pap1LNC.C = ggplot(df.MErrorsLNC.C.matlab, aes(x = Mean.Time, y = MAE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() 

pap1.2LNC.C = ggplot(df.MErrorsLNC.C.matlab, aes(x = log10(Mean.Time), y = log10(MAE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab("log10(Mean Time)") + ylab("log10(MAE)")

pap2LNC.C = ggplot(df.MErrorsLNC.C.matlab, aes(x = Mean.Time, y = RMSE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() 

pap2.2LNC.C = ggplot(df.MErrorsLNC.C.matlab, aes(x = log10(Mean.Time), y = log10(RMSE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab("log10(Mean Time)") + ylab("log10(RMSE)")

pap3LNC.C = ggplot(df.MErrorsLNC.C.matlab, aes(x = Mean.Time, y = MAPE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() 

pap3.2LNC.C = ggplot(df.MErrorsLNC.C.matlab, aes(x = log10(Mean.Time), y = log10(MAPE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes) + facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab("log10(Mean Time)") + ylab("log10(MAPE)")


pdf("MErrorsLNC.C.matlab.pdf", width = 16, height = 9)
plot(pap1.2LNC.C)
plot(pap2.2LNC.C)
plot(pap3.2LNC.C)
dev.off()

########################### Do the same for the UC setup

### HNC

DataHNCUC = readMat("ResultsReturnsHNCUC.mat")

MAEmeanHNCUC = DataHNCUC$MAEmeanHNCUC
MAEvarHNCUC  = DataHNCUC$MAEvarHNCUC
MAEskewHNCUC = DataHNCUC$MAEskewHNCUC
MAEkurtHNCUC = DataHNCUC$MAEkurtHNCUC
RMSEmeanHNCUC = DataHNCUC$RMSEmeanHNCUC
RMSEvarHNCUC = DataHNCUC$RMSEvarHNCUC
RMSEskewHNCUC = DataHNCUC$RMSEskewHNCUC
RMSEkurtHNCUC = DataHNCUC$RMSEkurtHNCUC
MAPEmeanHNCUC = DataHNCUC$MAPEmeanHNCUC
MAPEvarHNCUC = DataHNCUC$MAPEvarHNCUC
MAPEskewHNCUC = DataHNCUC$MAPEskewHNCUC
MAPEkurtHNCUC = DataHNCUC$MAPEkurtHNCUC
MetimeHNCUC  = DataHNCUC$MetimeHNCUC

n = 1
T = c(1/12, 1/4, 1/2, 1)
Delta = T / n

nEUC = 2
nIUC = 1
nOUC = 7


namesUC = factor(c('EFT', 'TVS', 'KJ', 'QE',  'BK', 'S', 'GK', 'TW', 'BBG', 'ORS'),  levels =c('EFT', 'TVS', 'KJ', 'QE',  'BK', 'S', 'GK', 'TW', 'BBG', 'ORS'))
Moments = factor(c("Mean", "Variance", "Skewness", "Kurtosis"), levels = c("Mean", "Variance", "Skewness", "Kurtosis"))
groups =  factor(c("Explicit", "Implicit", "Others"), levels = c("Explicit", "Implicit", "Others")) 

shapes =  c(seq(0, length(namesUC) -1))

df.MErrorsHNC.UC.matlab = data.frame("Scheme" = rep(rep(namesUC, each = length(Delta)), length(Moments)), "Group" = rep(rep(groups, c(nEUC * length(Delta), nIUC * length(Delta), nOUC * length(Delta))), length(Moments)), 
                                    "Delta" = rep(rep(Delta, length(namesUC)), length(Moments)), "Moment" = rep(Moments, each = length(namesUC) * length(Delta)), 
                                    "MAE" = c(c(MAEmeanHNCUC), c(MAEvarHNCUC), c(MAEskewHNCUC), c(MAEkurtHNCUC)), 
                                    "RMSE" = c(c(RMSEmeanHNCUC), c(RMSEvarHNCUC), c(RMSEskewHNCUC), c(RMSEkurtHNCUC)), 
                                    "MAPE" = c(c(MAPEmeanHNCUC), c(MAPEvarHNCUC), c(MAPEskewHNCUC), c(MAPEkurtHNCUC)), 
                                    "Mean Time" = rep(c(MetimeHNCUC), length(Moments)))

save(df.MErrorsHNC.UC.matlab, file = "df.MErrorsHNC.UC.matlab.RData")

# Plots which perfectly blend in with other paper: time against bias, exactly what we want

pap1HNC.UC = ggplot(df.MErrorsHNC.UC.matlab, aes(x = Delta, y = MAE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(Delta)) + scale_x_reverse()

pap1.2HNC.UC = ggplot(df.MErrorsHNC.UC.matlab, aes(x = log10(Delta), y = log10(MAE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(log10(Delta))) + ylab("log10(MAE)") + scale_x_reverse()

pap2HNC.UC = ggplot(df.MErrorsHNC.UC.matlab, aes(x = Delta, y = RMSE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(Delta)) + ylab("RMSE") + scale_x_reverse()

pap2.2HNC.UC = ggplot(df.MErrorsHNC.UC.matlab, aes(x = log10(Delta), y = log10(RMSE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(log10(Delta))) + ylab("log10(RMSE)") + scale_x_reverse()

pap3HNC.UC = ggplot(df.MErrorsHNC.UC.matlab, aes(x = Delta, y = MAPE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(Delta)) + ylab("MAPE") + scale_x_reverse()

pap3.2HNC.UC = ggplot(df.MErrorsHNC.UC.matlab, aes(x = log10(Delta), y = log10(MAPE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(log10(Delta))) + ylab("log10(MAPE)") + scale_x_reverse()

pdf("MErrorsHNC.UC.matlab.pdf", width = 16, height = 9)
plot(pap1.2HNC.UC)
plot(pap2.2HNC.UC)
plot(pap3.2HNC.UC)
dev.off()

### LNC

DataLNCUC = readMat("ResultsReturnsLNCUC.mat")

MAEmeanLNCUC = DataLNCUC$MAEmeanLNCUC
MAEvarLNCUC  = DataLNCUC$MAEvarLNCUC
MAEskewLNCUC = DataLNCUC$MAEskewLNCUC
MAEkurtLNCUC = DataLNCUC$MAEkurtLNCUC
RMSEmeanLNCUC = DataLNCUC$RMSEmeanLNCUC
RMSEvarLNCUC = DataLNCUC$RMSEvarLNCUC
RMSEskewLNCUC = DataLNCUC$RMSEskewLNCUC
RMSEkurtLNCUC = DataLNCUC$RMSEkurtLNCUC
MAPEmeanLNCUC = DataLNCUC$MAPEmeanLNCUC
MAPEvarLNCUC = DataLNCUC$MAPEvarLNCUC
MAPEskewLNCUC = DataLNCUC$MAPEskewLNCUC
MAPEkurtLNCUC = DataLNCUC$MAPEkurtLNCUC
MetimeLNCUC  = DataLNCUC$MetimeLNCUC

df.MErrorsLNC.UC.matlab = data.frame("Scheme" = rep(rep(namesUC, each = length(Delta)), length(Moments)), "Group" = rep(rep(groups, c(nEUC * length(Delta), nIUC * length(Delta), nOUC * length(Delta))), length(Moments)), 
                                     "Delta" = rep(rep(Delta, length(namesUC)), length(Moments)), "Moment" = rep(Moments, each = length(namesUC) * length(Delta)), 
                                     "MAE" = c(c(MAEmeanLNCUC), c(MAEvarLNCUC), c(MAEskewLNCUC), c(MAEkurtLNCUC)), 
                                     "RMSE" = c(c(RMSEmeanLNCUC), c(RMSEvarLNCUC), c(RMSEskewLNCUC), c(RMSEkurtLNCUC)), 
                                     "MAPE" = c(c(MAPEmeanLNCUC), c(MAPEvarLNCUC), c(MAPEskewLNCUC), c(MAPEkurtLNCUC)), 
                                     "Mean Time" = rep(c(MetimeLNCUC), length(Moments)))

save(df.MErrorsLNC.UC.matlab, file = "df.MErrorsLNC.UC.matlab.RData")

pap1LNC.UC = ggplot(df.MErrorsLNC.UC.matlab, aes(x = Delta, y = MAE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(Delta)) + scale_x_reverse()

pap1.2LNC.UC = ggplot(df.MErrorsLNC.UC.matlab, aes(x = log10(Delta), y = log10(MAE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(log10(Delta))) + ylab("log10(MAE)") + scale_x_reverse()

pap2LNC.UC = ggplot(df.MErrorsLNC.UC.matlab, aes(x = Delta, y = RMSE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(Delta)) + ylab("RMSE") + scale_x_reverse()

pap2.2LNC.UC = ggplot(df.MErrorsLNC.UC.matlab, aes(x = log10(Delta), y = log10(RMSE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(log10(Delta))) + ylab("log10(RMSE)") + scale_x_reverse()

pap3LNC.UC = ggplot(df.MErrorsLNC.UC.matlab, aes(x = Delta, y = MAPE, color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(Delta)) + ylab("MAPE") + scale_x_reverse()

pap3.2LNC.UC = ggplot(df.MErrorsLNC.UC.matlab, aes(x = log10(Delta), y = log10(MAPE), color = Scheme, shape = Scheme)) + 
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(log10(Delta))) + ylab("log10(MAPE)") + scale_x_reverse()

pdf("MErrorsLNC.UC.matlab.pdf", width = 16, height = 9)
plot(pap1.2LNC.UC)
plot(pap2.2LNC.UC)
plot(pap3.2LNC.UC)
dev.off()



###################################### Plot the comparison with a higher sample size (M = 1e6)


namesUCadd = factor(c('EFT', 'TVS', 'KJ', 'QE', 'BBG', 'ORS'),  levels =c('EFT', 'TVS', 'KJ', 'QE', 'BBG', 'ORS'))

nEUCadd = 2
nIUCadd = 1
nOUCadd = 3


DataHNCUCadd = readMat("ResultsReturnsHNCUCadd.mat")

MAEmeanHNCUCadd = DataHNCUCadd$MAEmeanHNCUC
MAEvarHNCUCadd  = DataHNCUCadd$MAEvarHNCUC
MAEskewHNCUCadd = DataHNCUCadd$MAEskewHNCUC
MAEkurtHNCUCadd = DataHNCUCadd$MAEkurtHNCUC
RMSEmeanHNCUCadd = DataHNCUCadd$RMSEmeanHNCUC
RMSEvarHNCUCadd = DataHNCUCadd$RMSEvarHNCUC
RMSEskewHNCUCadd = DataHNCUCadd$RMSEskewHNCUC
RMSEkurtHNCUCadd = DataHNCUCadd$RMSEkurtHNCUC
MAPEmeanHNCUCadd = DataHNCUCadd$MAPEmeanHNCUC
MAPEvarHNCUCadd = DataHNCUCadd$MAPEvarHNCUC
MAPEskewHNCUCadd = DataHNCUCadd$MAPEskewHNCUC
MAPEkurtHNCUCadd = DataHNCUCadd$MAPEkurtHNCUC
MetimeHNCUCadd  = DataHNCUCadd$MetimeHNCUC

df.MErrorsHNCUCadd.matlab = data.frame("Scheme" = rep(rep(namesUCadd, each = length(Delta)), length(Moments)), "Group" = rep(rep(groups, c(nEUCadd * length(Delta), nIUCadd * length(Delta), nOUCadd * length(Delta))), length(Moments)),
                                     "Delta" = rep(rep(Delta, length(namesUCadd)), length(Moments)), "Moment" = rep(Moments, each = length(namesUCadd) * length(Delta)),
                                     "MAE" = c(c(MAEmeanHNCUCadd), c(MAEvarHNCUCadd), c(MAEskewHNCUCadd), c(MAEkurtHNCUCadd)),
                                     "RMSE" = c(c(RMSEmeanHNCUCadd), c(RMSEvarHNCUCadd), c(RMSEskewHNCUCadd), c(RMSEkurtHNCUCadd)),
                                     "MAPE" = c(c(MAPEmeanHNCUCadd), c(MAPEvarHNCUCadd), c(MAPEskewHNCUCadd), c(MAPEkurtHNCUCadd)),
                                     "Mean Time" = rep(c(MetimeHNCUCadd), length(Moments)))

shapes =  c(seq(0, length(namesUCadd) -1))


pap1HNCUCadd = ggplot(df.MErrorsHNCUCadd.matlab, aes(x = Delta, y = MAE, color = Scheme, shape = Scheme)) +
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(Delta)) + scale_x_reverse()

pap1.2HNCUCadd = ggplot(df.MErrorsHNCUCadd.matlab, aes(x = log10(Delta), y = log10(MAE), color = Scheme, shape = Scheme)) +
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(log10(Delta))) + ylab("log10(MAE)") + scale_x_reverse()

pap2HNCUCadd = ggplot(df.MErrorsHNCUCadd.matlab, aes(x = Delta, y = RMSE, color = Scheme, shape = Scheme)) +
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(Delta)) + ylab("RMSE") + scale_x_reverse()

pap2.2HNCUCadd = ggplot(df.MErrorsHNCUCadd.matlab, aes(x = log10(Delta), y = log10(RMSE), color = Scheme, shape = Scheme)) +
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(log10(Delta))) + ylab("log10(RMSE)") + scale_x_reverse()

pap3HNCUCadd = ggplot(df.MErrorsHNCUCadd.matlab, aes(x = Delta, y = MAPE, color = Scheme, shape = Scheme)) +
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(Delta)) + ylab("MAPE") + scale_x_reverse()

pap3.2HNCUCadd = ggplot(df.MErrorsHNCUCadd.matlab, aes(x = log10(Delta), y = log10(MAPE), color = Scheme, shape = Scheme)) +
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(log10(Delta))) + ylab("log10(MAPE)") + scale_x_reverse()

pdf("MErrorsHNC.UC.matlab.add.pdf", width = 16, height = 9)
plot(pap1.2HNCUCadd)
plot(pap2.2HNCUCadd)
plot(pap3.2HNCUCadd)
dev.off()


#### LNC 

DataLNCUCadd = readMat("ResultsReturnsLNCUCadd.mat")

MAEmeanLNCUCadd = DataLNCUCadd$MAEmeanLNCUC
MAEvarLNCUCadd  = DataLNCUCadd$MAEvarLNCUC
MAEskewLNCUCadd = DataLNCUCadd$MAEskewLNCUC
MAEkurtLNCUCadd = DataLNCUCadd$MAEkurtLNCUC
RMSEmeanLNCUCadd = DataLNCUCadd$RMSEmeanLNCUC
RMSEvarLNCUCadd = DataLNCUCadd$RMSEvarLNCUC
RMSEskewLNCUCadd = DataLNCUCadd$RMSEskewLNCUC
RMSEkurtLNCUCadd = DataLNCUCadd$RMSEkurtLNCUC
MAPEmeanLNCUCadd = DataLNCUCadd$MAPEmeanLNCUC
MAPEvarLNCUCadd = DataLNCUCadd$MAPEvarLNCUC
MAPEskewLNCUCadd = DataLNCUCadd$MAPEskewLNCUC
MAPEkurtLNCUCadd = DataLNCUCadd$MAPEkurtLNCUC
MetimeLNCUCadd  = DataLNCUCadd$MetimeLNCUC

df.MErrorsLNCUCadd.matlab = data.frame("Scheme" = rep(rep(namesUCadd, each = length(Delta)), length(Moments)), "Group" = rep(rep(groups, c(nEUCadd * length(Delta), nIUCadd * length(Delta), nOUCadd * length(Delta))), length(Moments)),
                                       "Delta" = rep(rep(Delta, length(namesUCadd)), length(Moments)), "Moment" = rep(Moments, each = length(namesUCadd) * length(Delta)),
                                       "MAE" = c(c(MAEmeanLNCUCadd), c(MAEvarLNCUCadd), c(MAEskewLNCUCadd), c(MAEkurtLNCUCadd)),
                                       "RMSE" = c(c(RMSEmeanLNCUCadd), c(RMSEvarLNCUCadd), c(RMSEskewLNCUCadd), c(RMSEkurtLNCUCadd)),
                                       "MAPE" = c(c(MAPEmeanLNCUCadd), c(MAPEvarLNCUCadd), c(MAPEskewLNCUCadd), c(MAPEkurtLNCUCadd)),
                                       "Mean Time" = rep(c(MetimeLNCUCadd), length(Moments)))



pap1LNCUCadd = ggplot(df.MErrorsLNCUCadd.matlab, aes(x = Delta, y = MAE, color = Scheme, shape = Scheme)) +
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(Delta)) + scale_x_reverse()

pap1.2LNCUCadd = ggplot(df.MErrorsLNCUCadd.matlab, aes(x = log10(Delta), y = log10(MAE), color = Scheme, shape = Scheme)) +
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(log10(Delta))) + ylab("log10(MAE)") + scale_x_reverse()

pap2LNCUCadd = ggplot(df.MErrorsLNCUCadd.matlab, aes(x = Delta, y = RMSE, color = Scheme, shape = Scheme)) +
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(Delta)) + ylab("RMSE") + scale_x_reverse()

pap2.2LNCUCadd = ggplot(df.MErrorsLNCUCadd.matlab, aes(x = log10(Delta), y = log10(RMSE), color = Scheme, shape = Scheme)) +
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(log10(Delta))) + ylab("log10(RMSE)") + scale_x_reverse()

pap3LNCUCadd = ggplot(df.MErrorsLNCUCadd.matlab, aes(x = Delta, y = MAPE, color = Scheme, shape = Scheme)) +
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(Delta)) + ylab("MAPE") + scale_x_reverse()

pap3.2LNCUCadd = ggplot(df.MErrorsLNCUCadd.matlab, aes(x = log10(Delta), y = log10(MAPE), color = Scheme, shape = Scheme)) +
  geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
  xlab(expression(log10(Delta))) + ylab("log10(MAPE)") + scale_x_reverse()

pdf("MErrorsLNC.UC.matlab.add.pdf", width = 16, height = 9)
plot(pap1.2LNCUCadd)
plot(pap2.2LNCUCadd)
plot(pap3.2LNCUCadd)
dev.off()

# namesQEORSUC =  factor(c('QE', 'ORS'),  levels = c('QE', 'ORS'))
# 
# nEQEORSUC = 0
# nIQEORSUC = 0
# nOQEORSUC = 2
# 
# DataComQEORSHNC = readMat("ResultsCompQEORSHNC.mat")
# 
# MAEmeanQEORSHNC = DataComQEORSHNC$MAEmeanHNCUC
# MAEvarQEORSHNC  = DataComQEORSHNC$MAEvarHNCUC
# MAEskewQEORSHNC = DataComQEORSHNC$MAEskewHNCUC
# MAEkurtQEORSHNC = DataComQEORSHNC$MAEkurtHNCUC
# RMSEmeanQEORSHNC = DataComQEORSHNC$RMSEmeanHNCUC
# RMSEvarQEORSHNC = DataComQEORSHNC$RMSEvarHNCUC
# RMSEskewQEORSHNC = DataComQEORSHNC$RMSEskewHNCUC
# RMSEkurtQEORSHNC = DataComQEORSHNC$RMSEkurtHNCUC
# MAPEmeanQEORSHNC = DataComQEORSHNC$MAPEmeanHNCUC
# MAPEvarQEORSHNC = DataComQEORSHNC$MAPEvarHNCUC
# MAPEskewQEORSHNC = DataComQEORSHNC$MAPEskewHNCUC
# MAPEkurtQEORSHNC = DataComQEORSHNC$MAPEkurtHNCUC
# MetimeQEORSHNC  = DataComQEORSHNC$MetimeHNCUC
# 
# df.MErrorsQEORSHNC.matlab = data.frame("Scheme" = rep(rep(namesQEORSUC, each = length(Delta)), length(Moments)), "Group" = rep(rep(groups, c(nEQEORSUC * length(Delta), nIQEORSUC * length(Delta), nOQEORSUC * length(Delta))), length(Moments)), 
#                                      "Delta" = rep(rep(Delta, length(namesQEORSUC)), length(Moments)), "Moment" = rep(Moments, each = length(namesQEORSUC) * length(Delta)), 
#                                      "MAE" = c(c(MAEmeanQEORSHNC), c(MAEvarQEORSHNC), c(MAEskewQEORSHNC), c(MAEkurtQEORSHNC)), 
#                                      "RMSE" = c(c(RMSEmeanQEORSHNC), c(RMSEvarQEORSHNC), c(RMSEskewQEORSHNC), c(RMSEkurtQEORSHNC)), 
#                                      "MAPE" = c(c(MAPEmeanQEORSHNC), c(MAPEvarQEORSHNC), c(MAPEskewQEORSHNC), c(MAPEkurtQEORSHNC)), 
#                                      "Mean Time" = rep(c(MetimeQEORSHNC), length(Moments)))
# 
# 
# shapes =  c(seq(0, length(namesQEORSUC) -1))
# 
# 
# pap1QEORSHNC = ggplot(df.MErrorsQEORSHNC.matlab, aes(x = Delta, y = MAE, color = Scheme, shape = Scheme)) + 
#   geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
#   xlab(expression(Delta)) + scale_x_reverse()
# 
# pap1.2QEORSHNC = ggplot(df.MErrorsQEORSHNC.matlab, aes(x = log10(Delta), y = log10(MAE), color = Scheme, shape = Scheme)) + 
#   geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
#   xlab(expression(log10(Delta))) + ylab("log10(MAE)") + scale_x_reverse()
# 
# pap2QEORSHNC = ggplot(df.MErrorsQEORSHNC.matlab, aes(x = Delta, y = RMSE, color = Scheme, shape = Scheme)) + 
#   geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
#   xlab(expression(Delta)) + ylab("RMSE") + scale_x_reverse()
# 
# pap2.2QEORSHNC = ggplot(df.MErrorsQEORSHNC.matlab, aes(x = log10(Delta), y = log10(RMSE), color = Scheme, shape = Scheme)) + 
#   geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
#   xlab(expression(log10(Delta))) + ylab("log10(RMSE)") + scale_x_reverse()
# 
# pap3QEORSHNC = ggplot(df.MErrorsQEORSHNC.matlab, aes(x = Delta, y = MAPE, color = Scheme, shape = Scheme)) + 
#   geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
#   xlab(expression(Delta)) + ylab("MAPE") + scale_x_reverse()
# 
# pap3.2QEORSHNC = ggplot(df.MErrorsQEORSHNC.matlab, aes(x = log10(Delta), y = log10(MAPE), color = Scheme, shape = Scheme)) + 
#   geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
#   xlab(expression(log10(Delta))) + ylab("log10(MAPE)") + scale_x_reverse()
# 
# pdf("CompQEORSHNC.pdf")
# plot(pap1.2QEORSHNC)
# plot(pap2.2QEORSHNC)
# plot(pap3.2QEORSHNC)
# dev.off()
# 
# DataComQEORSLNC = readMat("ResultsCompQEORSLNC.mat")
# 
# MAEmeanQEORSLNC = DataComQEORSLNC$MAEmeanLNCUC
# MAEvarQEORSLNC  = DataComQEORSLNC$MAEvarLNCUC
# MAEskewQEORSLNC = DataComQEORSLNC$MAEskewLNCUC
# MAEkurtQEORSLNC = DataComQEORSLNC$MAEkurtLNCUC
# RMSEmeanQEORSLNC = DataComQEORSLNC$RMSEmeanLNCUC
# RMSEvarQEORSLNC = DataComQEORSLNC$RMSEvarLNCUC
# RMSEskewQEORSLNC = DataComQEORSLNC$RMSEskewLNCUC
# RMSEkurtQEORSLNC = DataComQEORSLNC$RMSEkurtLNCUC
# MAPEmeanQEORSLNC = DataComQEORSLNC$MAPEmeanLNCUC
# MAPEvarQEORSLNC = DataComQEORSLNC$MAPEvarLNCUC
# MAPEskewQEORSLNC = DataComQEORSLNC$MAPEskewLNCUC
# MAPEkurtQEORSLNC = DataComQEORSLNC$MAPEkurtLNCUC
# MetimeQEORSLNC  = DataComQEORSLNC$MetimeLNCUC
# 
# df.MErrorsQEORSLNC.matlab = data.frame("Scheme" = rep(rep(namesQEORSUC, each = length(Delta)), length(Moments)), "Group" = rep(rep(groups, c(nEQEORSUC * length(Delta), nIQEORSUC * length(Delta), nOQEORSUC * length(Delta))), length(Moments)), 
#                                        "Delta" = rep(rep(Delta, length(namesQEORSUC)), length(Moments)), "Moment" = rep(Moments, each = length(namesQEORSUC) * length(Delta)), 
#                                        "MAE" = c(c(MAEmeanQEORSLNC), c(MAEvarQEORSLNC), c(MAEskewQEORSLNC), c(MAEkurtQEORSLNC)), 
#                                        "RMSE" = c(c(RMSEmeanQEORSLNC), c(RMSEvarQEORSLNC), c(RMSEskewQEORSLNC), c(RMSEkurtQEORSLNC)), 
#                                        "MAPE" = c(c(MAPEmeanQEORSLNC), c(MAPEvarQEORSLNC), c(MAPEskewQEORSLNC), c(MAPEkurtQEORSLNC)), 
#                                        "Mean Time" = rep(c(MetimeQEORSLNC), length(Moments)))
# 
# 
# 
# 
# pap1QEORSLNC = ggplot(df.MErrorsQEORSLNC.matlab, aes(x = Delta, y = MAE, color = Scheme, shape = Scheme)) + 
#   geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
#   xlab(expression(Delta)) + scale_x_reverse()
# 
# pap1.2QEORSLNC = ggplot(df.MErrorsQEORSLNC.matlab, aes(x = log10(Delta), y = log10(MAE), color = Scheme, shape = Scheme)) + 
#   geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
#   xlab(expression(log10(Delta))) + ylab("log10(MAE)") + scale_x_reverse()
# 
# pap2QEORSLNC = ggplot(df.MErrorsQEORSLNC.matlab, aes(x = Delta, y = RMSE, color = Scheme, shape = Scheme)) + 
#   geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
#   xlab(expression(Delta)) + ylab("RMSE") + scale_x_reverse()
# 
# pap2.2QEORSLNC = ggplot(df.MErrorsQEORSLNC.matlab, aes(x = log10(Delta), y = log10(RMSE), color = Scheme, shape = Scheme)) + 
#   geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
#   xlab(expression(log10(Delta))) + ylab("log10(RMSE)") + scale_x_reverse()
# 
# pap3QEORSLNC = ggplot(df.MErrorsQEORSLNC.matlab, aes(x = Delta, y = MAPE, color = Scheme, shape = Scheme)) + 
#   geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
#   xlab(expression(Delta)) + ylab("MAPE") + scale_x_reverse()
# 
# pap3.2QEORSLNC = ggplot(df.MErrorsQEORSLNC.matlab, aes(x = log10(Delta), y = log10(MAPE), color = Scheme, shape = Scheme)) + 
#   geom_point() + geom_line() + scale_shape_manual(values = shapes)+ facet_grid(.~Moment, scales = "free") + theme_bw() +
#   xlab(expression(log10(Delta))) + ylab("log10(MAPE)") + scale_x_reverse()
# 
# pdf("CompQEORSLNC.pdf")
# plot(pap1.2QEORSLNC)
# plot(pap2.2QEORSLNC)
# plot(pap3.2QEORSLNC)
# dev.off()
