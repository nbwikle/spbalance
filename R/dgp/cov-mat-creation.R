### cov-mat-creation.R
### Nathan Wikle

### Create and save Gaussian process covariance matrices used
###   in Tec et al. (2023) simulations.

### Important: You only have to perform this once!


#-----------------------------------#
#--- 1. Packages and source code ---#
#-----------------------------------#

# R package
library(here)
# source code
source(here::here("R", "dgp", "data-gen-1.R"))
# check that directories exist for matrix decompositions
dir.create(here::here("data"), showWarnings = FALSE)
dir.create(here::here("data", "cov-mats"), showWarnings = FALSE)


#-------------------------------------------------------------------------------#
#--- 2. Visual comparisons of cov. functions with different range parameters ---#
#-------------------------------------------------------------------------------#

# Matern Covariance:
par(mfrow = c(1,1))
x.vals <- seq(from = 0, to = 20, by = 0.01)

png(here::here("data", "cov-mats", "Matern-plot.png"),
    width = 7, height = 5, units = "in", res = 300)
plot(x.vals, covFun(x.vals, theta = c(3, 1.5), type = "Matern"), type = "l",
     bty = "l", xlab = "distance", ylab = "correlation")
lines(x.vals, covFun(x.vals, theta = c(5, 1.5), type = "Matern"), type = "l", lty = 2, col = "red")
lines(x.vals, covFun(x.vals, theta = c(10, 1.5), type = "Matern"), type = "l", col = "blue", lty = 3)
title(main = "Matern (nu = 1.5)", line = 0, adj = 0.05)
legend(15, 1, legend=c("rho = 3", "rho = 5", "rho = 10"),
       col=c("black", "red", "blue"), lty=1:3, cex=0.8, box.lty=0)
dev.off()

# Squared Exponential Covariance:
png(here::here("data", "cov-mats", "SE-plot.png"),
    width = 7, height = 5, units = "in", res = 300)
plot(x.vals, covFun(x.vals, theta = c(3 / 2), type = "SE"), type = "l",
     bty = "l", xlab = "distance", ylab = "correlation")
lines(x.vals, covFun(x.vals, theta = c(5 / 2), type = "SE"), type = "l", col = "red", lty = 2)
lines(x.vals, covFun(x.vals, theta = c(10 / 2), type = "SE"), type = "l", col = "blue", lty = 3)
title(main = "Squared Exponential", line = 0, adj = 0.05)
legend(15, 1, legend=c("rho = 1.5", "rho = 2.5", "rho = 5.5"),
       col=c("black", "red", "blue"), lty=1:3, cex=0.8, box.lty=0)
dev.off()

# Exponential Covariance:
png(here::here("data", "cov-mats", "Expo-plot.png"),
    width = 7, height = 5, units = "in", res = 300)
plot(x.vals, covFun(x.vals, theta = c(3 / 2), type = "Expo"), type = "l",
     bty = "l", xlab = "distance", ylab = "correlation")
lines(x.vals, covFun(x.vals, theta = c(5 / 2), type = "Expo"), type = "l", col = "red", lty = 2)
lines(x.vals, covFun(x.vals, theta = c(10 / 2), type = "Expo"), type = "l", col = "blue", lty = 3)
title(main = "Exponential", line = 0, adj = 0.05)
legend(15, 1, legend=c("rho = 1.5", "rho = 2.5", "rho = 5.5"),
       col=c("black", "red", "blue"), lty=1:3, cex=0.8, box.lty=0)
dev.off()

#--------------------------------------------------------#
#--- 3. Create GP covariances and save decompositions ---#
#--------------------------------------------------------#

# Important: This only needs to be performed once!
#   -> requires ~ 5.5 GB of space
#   -> time on my laptop: ~11 hours

# GP1 (squared exponential, rho = 5.1, on bigger lattice)
simDecomps(type = "GP1")
# Squared exponential (rho = 3/2, 5/2, 10/2)
simDecomps(type = "SE")
# Exponential (rho = 3/2, 5/2, 10/2)
simDecomps(type = "Expo")
# Matern (rho = 3, 5, 10)
simDecomps(type = "Matern")
