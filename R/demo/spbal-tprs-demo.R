### spbal-tprs-demo.R
### Nathan Wikle

### A short demo of 'spbalance', with a thin-plate regression spline (TPRS)
###   to model the unmeasured spatial confounder.

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

# source(here::here("R", "aIPTW.R"))

# indicate that we want to use a specific Python virtualenv
use_virtualenv("r-reticulate")

# load python libraries
scipy <- import("scipy")
np <- import("numpy")
cv <- import("cv2")
# source Python functions
source_python(here::here("src", "tec-simdata.py"))

#--------------------------------------------#
#--- 2. Simulate data from Tec et al. DGP ---#
#--------------------------------------------#

# set seed
set.seed(33)

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


#--------------------------------------------#
#--- 3. ATE estimation using GAM          ---#
#--------------------------------------------#

# estimate propensity score via mgcv::gam
gam.fit <- mgcv::gam(form = treat ~ s(x, y, k = 500), family = binomial, data = sim.df)
# gam.check(gam.fit)
pi.hat <- gam.fit$fitted.values

# ATE estimate (true value = 0.1)
ateEst(out = sim.df$outcome, trt = sim.df$treat, pi.hat)$ate

# # Note: the above is all performed by the 'ateGAM' function:
# iptw.gam <- ateGAM(
#   form = treat ~ s(x, y, k = 300),
#   data = sim.df, out = sim.df$outcome, est.type = "all"
# )
# iptw.gam$ate

#--------------------------------------------#
#--- 4. ATE estimation using spbalance    ---#
#--------------------------------------------#

# select tuning parameter using coefficient of variation method
spb.cvr <- spBalance(
  formula = treat ~ s(x, y, k = 500),
  data = sim.df,
  lambda = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 5, 10),
  tuning = "coefvar",
  hide.details = TRUE,
  opt.params = list(tol = 1e-7, max.iter = 1000, alpha = 0.5, beta = 0.5)
)

# chosen lambda:
spb.cvr$lambda

# ATE estimate (true value = 0.1)
ateEst(out = sim.df$outcome, trt = sim.df$treat, pr.trt = spb.cvr$pi.hat)$ate

# # Note: the above can be performed using the 'ateBAL' function:
# iptw.bal <- ateBAL(
#   form = treat ~ s(x, y, k = 500),
#   data = sim.df,
#   out = sim.df$outcome,
#   est.type = "all",
#   opt.params = list(tol = 1e-7, max.iter = 1000, alpha = 0.5, beta = 0.5),
#   lambda.vals = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 5, 10),
#   tuning = "coefvar",
#   coefvar.r = 0.9
# )
# iptw.bal$ate
