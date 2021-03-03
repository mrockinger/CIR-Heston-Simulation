# Handle matlab results

library(R.matlab)
library(ggplot2)

DataHNCC = readMat("ResultsVarianceHNCC.mat")

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

nEHNC = 7
nIHNC = 2
nOHNC = 5


namesHNC = factor(c("AE0", "DD", "EAF", "ERF", "EFT", "TVS", "HM", "KJ", "G", "ABR", "QE", "SSA", "SNA", "BK"),  levels =c("AE0", "DD", "EAF", "ERF", "EFT", "TVS", "HM", "KJ", "G", "ABR", "QE", "SSA", "SNA", "BK"))
Moments = factor(c("Mean", "Variance", "Skewness", "Kurtosis"), levels = c("Mean", "Variance", "Skewness", "Kurtosis"))
groups =  factor(c("Explicit", "Implicit", "Others"), levels = c("Explicit", "Implicit", "Others")) 

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colors = gg_color_hue(3) # Mimic the colors, ggplot would choose, if we mapped 3 variables to colors
shapes =  c(seq(0, length(namesHNC) -1))

df.MErrorsHNC.C.matlab = data.frame("Scheme" = rep(rep(namesHNC, each = length(Delta)), length(Moments)), "Group" = rep(rep(groups, c(nEHNC * length(Delta), nIHNC * length(Delta), nOHNC * length(Delta))), length(Moments)), 
                        "Delta" = rep(rep(Delta, length(namesHNC)), length(Moments)), "Moment" = rep(Moments, each = length(namesHNC) * length(Delta)), 
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

DataLNCC = readMat("ResultsVarianceLNCC.mat")

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

DataHNCUC = readMat("ResultsVarianceHNCUC.mat")

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

n = c(12, 4, 2, 1)
T = 1
Delta = T / n

nEHNC = 7
nIHNC = 2
nOHNC = 5


namesHNC = factor(c("AE0", "DD", "EAF", "ERF", "EFT", "TVS", "HM", "KJ", "G", "ABR", "QE", "SSA", "SNA", "BK"),  levels =c("AE0", "DD", "EAF", "ERF", "EFT", "TVS", "HM", "KJ", "G", "ABR", "QE", "SSA", "SNA", "BK"))
Moments = factor(c("Mean", "Variance", "Skewness", "Kurtosis"), levels = c("Mean", "Variance", "Skewness", "Kurtosis"))
groups =  factor(c("Explicit", "Implicit", "Others"), levels = c("Explicit", "Implicit", "Others")) 


df.MErrorsHNC.UC.matlab = data.frame("Scheme" = rep(rep(namesHNC, each = length(Delta)), length(Moments)), "Group" = rep(rep(groups, c(nEHNC * length(Delta), nIHNC * length(Delta), nOHNC * length(Delta))), length(Moments)), 
                                    "Delta" = rep(rep(Delta, length(namesHNC)), length(Moments)), "Moment" = rep(Moments, each = length(namesHNC) * length(Delta)), 
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

DataLNCUC = readMat("ResultsVarianceLNCUC.mat")

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

nELNC = 7
nILNC = 1
nOLNC = 3


namesLNC = factor(c("AE0", "DD", "EAF", "ERF", "EFT", "TVS", "HM", "KJ", "ABR", "QE", "BK"),  levels =c("AE0", "DD", "EAF", "ERF", "EFT", "TVS", "HM", "KJ", "ABR", "QE", "BK"))

df.MErrorsLNC.UC.matlab = data.frame("Scheme" = rep(rep(namesLNC, each = length(Delta)), length(Moments)), "Group" = rep(rep(groups, c(nELNC * length(Delta), nILNC * length(Delta), nOLNC * length(Delta))), length(Moments)), 
                                     "Delta" = rep(rep(Delta, length(namesLNC)), length(Moments)), "Moment" = rep(Moments, each = length(namesLNC) * length(Delta)), 
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
