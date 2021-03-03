# Get intuition on how the grid for the precomputation of moments in the BBG Algo is changing with zeta

BBGgrid = function(zeta, sigma, M){
  
  delta = exp(log(10 * sigma) / ((75 - zeta) + log2(M)))
  xi = log(10 * sigma) / log(delta)
  grid =  c(0, delta^(-zeta + (0:zeta)), delta^(xi - ((xi-1):0))) 
  
}

test =  lapply(52:94, BBGgrid, sigma = 0.4, M = 1e6)


pdf("BBGgrid.pdf")
plot(test[[1]], type = "l", ylab = expression(v[t] * v[t + Delta]), xlab = "Gridindex")
for(i in 2:length(test)){
  
  lines(test[[i]])
  
}
dev.off()
# Be careful with values in the grid becoming too large! 
# If this is the case, the moments cannot be computed, because the argument of the Bessel function is getting too large
# I. e. a "moderate" choice of zeta will be the most stable -> e. g. zeta = 85 
# Note: the grid will not exist if M is too small, i. e. M = 1e3
