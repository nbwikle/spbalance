### spbal-krr-demo.R
### Nathan Wikle

### A short demo of 'spbalance' using a kernel ridge regression term.

#--------------------------------------------#
#--- 1. Load packages and R scripts       ---#
#--------------------------------------------#

library(here)           # 1.0.1
library(terra)          # 1.7.29
library(geoR)           # 1.9.2
library(mgcv)           # 1.8.39
library(MASS)           # 7.3.55
library(matrixStats)    # 0.61.0
library(RColorBrewer)   # 1.1.2
library(raster)         # 3.5.15
library(reticulate)     # 1.28
library(sf)
library(Matrix)

source(here::here("R", "dgp", "data-gen-1.R"))
source(here::here("R", "spbalance.R"))
source(here::here("R", "tuning-params.R"))
source(here::here("R", "ate-est.R"))

# indicate that we want to use a specific Python virtualenv
use_virtualenv("r-reticulate")

# load python libraries
scipy <- import("scipy")
np <- import("numpy")
cv <- import("cv2")
# source Python functions
source_python(here::here("src", "tec-simdata.py"))

#--------------------------------------------#
#--- 2. Simulate data from Zhao (2019)    ---#
#--------------------------------------------#

### Simulate data from Zhao et al.
ZhaoSim <- function(n.sim, covariance){

  # 1. Simulate covariates
  X <- matrix(rnorm(n.sim * 5), nrow = n.sim, ncol = 5)
  colnames(X) <- c("x1", "x2", "x3", "x4", "x5")

  # 2. Obtain covariance matrix
  d.mat <- as.matrix(dist(X, diag = TRUE, upper = TRUE))
  Sigma <- kernelEval(d.mat, k.params = covariance$params, type = covariance$type, centered = FALSE)

  # 3. Propensity score
  gp.sims <- mvrnorm(n = 2, mu = rep(0, n.sim), Sigma)
  pi.0 <- expit(gp.sims[1,])

  # 4. Treatment assignment
  trt <- rbinom(n = n.sim, size = 1, prob = pi.0)

  # 5. Generate outcomes
  outcome <- gp.sims[2,] + rnorm(n.sim)

  # 6. Combine results
  list(
    X = X, prop.score = pi.0, trt = trt, out = outcome
  )
}



# set seed
set.seed(39)

# simulate data
zhao.sim <- ZhaoSim(
  n.sim = 500,
  covariance = list(
    type = "SE",
    params = c(5)
  )
)
sim.k <- cbind(zhao.sim$trt, zhao.sim$X)
colnames(sim.k)[1] <- "trt"



#----------------------------------------------------#
#--- 3. Fit balancing weights model with KRR term ---#
#----------------------------------------------------#


# Fit model
fit.krr <- spBalance(
  formula = trt ~ krr(x1, x2, x3, x4, x5, kern = "SE", kp = 5, centered = FALSE),
  data = sim.k,
  lambda = c(0.1, 0.25, 0.5, 0.75, 1, 2, 5, 10),
  tuning = "coefvar", coefvar.r = 0.5,
  hide.details = TRUE,
  opt.params = list(tol = 1e-7, max.iter = 1000, alpha = 0.5, beta = 0.5)
)
fit.krr$lambda
ateEst(out = zhao.sim$out, trt = zhao.sim$trt, pr.trt = fit.krr$pi.hat)$ate


#----------------------------------------------------#
#--- 4. Repeat, but using RFF approximation.      ---#
#----------------------------------------------------#


# Using random Fourier features (RFF); number of features = 250

fit.rff <- spBalance(
  formula = trt ~ krr(x1, x2, x3, x4, x5, kern = "SE", kp = 5, centered = FALSE, rff = TRUE, nf = 250),
  data = sim.k,
  lambda = c(0.1, 0.25, 0.5, 0.75, 1, 2, 5, 10),
  tuning = "coefvar", coefvar.r = 0.5,
  hide.details = TRUE,
  opt.params = list(tol = 1e-7, max.iter = 1000, alpha = 0.5, beta = 0.5)
)

fit.rff$lambda
ateEst(out = zhao.sim$out, trt = zhao.sim$trt, pr.trt = fit.rff$pi.hat)$ate



#-----------------------------------------------------------#
#--- 5. Simulate data from Tec et al. generating process ---#
#-----------------------------------------------------------#

# set seed
set.seed(27)

# simulate data
#   -> 25x25 raster
#   -> error = GP with Matern covariance; r.type = 3 (rho = 5)
#   -> linear version of Tec et al. sim
#   -> ATE = 0.1

sim.data <- dataGen1(
  sim.args = list(grid = 25, cov = "Matern", r.type = 3, tau = 0.1),
  K.mat = "default", nonlin = FALSE
)

# convert to SpatRast object
sim.rast <- dataGen1.SpatRaster(sim.data)

# # plot of covariates, y, propensity score, and treatment
# plot(sim.rast)

# create data frame
sim.df <- data.frame(
  treat = values(sim.rast$z),
  outcome = values(sim.rast$y),
  pi = values(sim.rast$pi),
  v1 = values(sim.rast$x1),
  v2 = values(sim.rast$x2),
  crds(sim.rast$y)
)
colnames(sim.df) <- c("treat", "outcome", "pi", "v1", "v2", "x", "y")


#-----------------------------------------------------------------#
#--- 6. Estimate ATE using spline + krr term (with RFF approx) ---#
#-----------------------------------------------------------------#

# select tuning parameter using coefficient of variation method
spb.rff <- spBalance(
  formula = treat ~ 1 + s(x, y, k = 500) +
    krr(v1, v2, kern = "SE", kp = 5, centered = TRUE, rff = TRUE, nf = 250),
  data = sim.df,
  lambda = c(0.05, 0.1, 0.25, 0.5, 0.75, 1, 5, 10),
  tuning = "coefvar",
  hide.details = TRUE,
  opt.params = list(tol = 1e-7, max.iter = 1000, alpha = 0.5, beta = 0.5)
)

# chosen lambda:
spb.rff$lambda

# ATE estimate (true value = 0.1)
ateEst(out = sim.df$outcome, trt = sim.df$treat, pr.trt = spb.rff$pi.hat)$ate


