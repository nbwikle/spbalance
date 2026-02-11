### dgp1-sim-study.R
### Nathan Wikle

### Recreate simulation studies from Tec et al.

#------------------------------------------------------------------------#
#--- 0. Process command line arguments                                ---#
#------------------------------------------------------------------------#

# extract command line args
args <- commandArgs(TRUE)
# save command line arguments
cov.type = args[1] # one of SE, Matern, or Expo
r.type = as.numeric(args[2]) # 1, 2, or 3
grid.size = as.numeric(args[3]) # specify n for an n x n grid
random.seed = as.numeric(args[4]) # where to start random seed
n.sims = as.numeric(args[5]) # number of simulations to perform
nonlin = args[6] # whether simulation is nonlinear or not

n.args <- length(args)
if (n.args >= 7){
  sparse = as.numeric(args[7])
} else if (n.args < 7){
  sparse = 0
}

if (nonlin == "T"){
  nonlinear = TRUE
} else {
  nonlinear = FALSE
}


#------------------------------------------------------------------------#
#--- 1. Load packages, Python, and source functions                   ---#
#------------------------------------------------------------------------#

# required packages
library(here)           # 1.0.1
library(terra)          # 1.8.54
library(geoR)           # 1.9.4
library(mgcv)           # 1.9.1
library(MASS)           # 7.3.61
library(matrixStats)    # 1.5.0
library(RColorBrewer)   # 1.1.3
library(raster)         # 3.6.30
library(reticulate)     # 1.37.0
library(sf)             # 1.0.16
library(Matrix)         # 1.7.0
library(spmodel)        # 0.8.0

# source code
source(here::here("R", "dgp", "data-gen-1.R"))
source(here::here("R", "spbalance.R"))
source(here::here("R", "tuning-params.R"))
source(here::here("R", "ate-est.R"))
source(here::here("R", "aIPTW.R"))

# indicate that we want to use a specific Python virtualenv
use_virtualenv("r-reticulate")

# load python libraries
scipy <- import("scipy")
np <- import("numpy")
cv <- import("cv2")
# source Python functions
source_python(here::here("src", "tec-simdata.py"))

# check that directories exist for saving results
dir.create(here::here("output"), showWarnings = FALSE)
dir.create(here::here("output", "dgp1"), showWarnings = FALSE)

#------------------------------------------------------------------------#
#--- 2. Simulation study in Tec et al. (2023)                         ---#
#------------------------------------------------------------------------#

tecIPTW <- function(
  n.sims, sim.start = 0, nonlin = FALSE,
  data.args = list(grid = 100, cov = "Matern", r.type = 3, tau = 0.1),
  car.fit = FALSE, sparse = 0
){

  # collect results
  results.iptw <- array(NA_real_, dim = c(n.sims, 3, 8))
  dimnames(results.iptw)[[1]] <- paste("sim-", sim.start + 1:n.sims, sep = "")
  dimnames(results.iptw)[[2]] <- c("HT", "Hajek", "OW")
  dimnames(results.iptw)[[3]] <- c("naive", "oracle", "glm", "gam", "gam0",
                                   "bal.cvr", "bal.mb", "bal.min")

  if (car.fit){
    results.aiptw <- matrix(NA_real_, nrow = n.sims, ncol = 9)
    rownames(results.aiptw) <- paste("sim-", sim.start + 1:n.sims, sep = "")
    colnames(results.aiptw) <- c("out.gam", "out.gam0", "aug.gbcvr", "aug.gbm",
      "aug.gg", "out.car", "aug.cbcvr", "aug.cbm", "aug.cg")
  } else {
    results.aiptw <- matrix(NA_real_, nrow = n.sims, ncol = 5)
    rownames(results.aiptw) <- paste("sim-", sim.start + 1:n.sims, sep = "")
    colnames(results.aiptw) <- c("out.gam", "out.gam0", "aug.gbcvr", "aug.gbm", "aug.gg")
  }

  for (k in 1:n.sims){

    ### 1. simulate data

    # set seed
    set.seed(sim.start + k)

    # simulate data
    sim.data <- dataGen1(sim.args = data.args, K.mat = "default", nonlin)

    # convert to SpatRast object
    sim.rast <- dataGen1.SpatRaster(sim.data)

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

    if (sparse > 0){
      keep.k <- sample(nrow(sim.df), sparse)
      sim.df <- sim.df[keep.k, ]
    }

    ### 2. naive differences
    results.iptw[k,,1] <- rep(mean(sim.df$outcome[which(sim.df$treat == 1)]) -
      mean(sim.df$outcome[which(sim.df$treat == 0)]), 3)

    ### 3. oracle estimate
    results.iptw[k,,2] <- ateEst(out = sim.df$outcome, trt = sim.df$treat, pr.trt = sim.df$pi)[[2]]

    ### 4. naive GLM
    ate.glm <- ateGLM(
      form = treat ~ v1 + v2 + x + y,
      data = sim.df, out = sim.df$outcome, est.type = "all"
    )
    results.iptw[k,,3] <- ate.glm

    ### 5. estimate using MGCV

    if (sparse > 0){
      n.knots = 300
    } else {
      n.knots = 500
    }
    weight.gam <- ateGAM(
      form = treat ~ s(x, y, k = n.knots),
      data = sim.df, out = sim.df$outcome, est.type = "all"
    )
    results.iptw[k,,4] <- weight.gam$ate



    ### 5b. estimate using MGCV with lambda = 0...

    weight.gam0 <- ateGAM(
      form = treat ~ s(x, y, k = n.knots, fx = TRUE),
      data = sim.df, out = sim.df$outcome, est.type = "all"
    )
    results.iptw[k,,5] <- weight.gam0$ate

    ### 6. estimate using spatial balance

    # coefvar method
    spb.cvr <- spBalance(
      formula = treat ~ s(x, y, k = 500),
      data = sim.df,
      lambda = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 5, 10),
      tuning = "coefvar",
      hide.details = TRUE,
      opt.params = list(tol = 1e-7, max.iter = 1000, alpha = 0.5, beta = 0.5)
    )
    results.iptw[k,,6] <- ateEst(
      out = sim.df$outcome, trt = sim.df$treat, pr.trt = spb.cvr$pi.hat
    )[[2]]

    # max balance method
    spb.mb <- spBalance(
      formula = treat ~ s(x, y, k = 500), #s(x,y,k=500),
      data = sim.df,
      lambda = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 5, 10),
      tuning = "max.bal",
      hide.details = TRUE,
      opt.params = list(tol = 1e-7, max.iter = 1000, alpha = 0.5, beta = 0.5)
    )
    results.iptw[k,,7] <- ateEst(
      out = sim.df$outcome, trt = sim.df$treat, pr.trt = spb.mb$pi.hat
    )[[2]]

    # minimum method
    spb.min <- spBalance(
      formula = treat ~ s(x, y, k = 500), #s(x,y,k=500),
      data = sim.df,
      lambda = c(0.001),
      tuning = "none",
      hide.details = TRUE,
      opt.params = list(tol = 1e-7, max.iter = 1000, alpha = 0.5, beta = 0.5)
    )
    results.iptw[k,,8] <- ateEst(
      out = sim.df$outcome, trt = sim.df$treat, pr.trt = spb.min$pi.hat
    )[[2]]


    ### 7. Estimates using outcome models

    # with penalty
    out.gam <- outcomeATE(
      formula = outcome ~ treat + s(x, y, k = 500),
      data = sim.df,
      family = gaussian(),
      trt.var = "treat",
      type = "gam"
    )
    results.aiptw[k,1] <- out.gam$ate


    # no penalty (lambda = 0)
    out.gam0 <- outcomeATE(
      formula = outcome ~ treat + s(x, y, k = 500, fx = TRUE),
      data = sim.df,
      family = gaussian(),
      trt.var = "treat",
      type = "gam"
    )

    results.aiptw[k,2] <- out.gam0$ate


    if (car.fit){
      out.car <- outcomeATE(
        formula = outcome ~ treat + car(x, y, grid = TRUE, icar = TRUE),
        data = sim.df,
        trt.var = "treat",
        type = "car",
        n.post = 20000,
        save.iter = 100,
        priors = list(
          beta = 1000, # variance of beta prior
          s2 = c(0.001, 0.001), # shape and rate parameter of IG prior on sigma^2
          tau2 = 0.1 * sqrt(pi/2), # half-normal prior on tau^2 with variance = 100
          rho = c(18, 2) # shape parameters of beta prior on rho
        ),
        s2.prop = 0.05,
        tau2.prop = 0.01,
        rho.prop = 0.01,
        keep.fit = TRUE
      )
      results.aiptw[k,6] <- out.car$ate
    }

    ### 8. Estimates using AIPTW

    # (i) gam, balancing (cvr)
    results.aiptw[k,3] <- aiptwATE(
      out = sim.df$outcome, trt = sim.df$treat, mu0 = out.gam$mu0,
      mu1 = out.gam$mu1, pi.hat = spb.cvr$pi.hat
    )

    # (ii) gam, balancing (min)
    results.aiptw[k,4] <- aiptwATE(
      out = sim.df$outcome, trt = sim.df$treat, mu0 = out.gam$mu0,
      mu1 = out.gam$mu1, pi.hat = as.vector(spb.min$pi.hat)
    )

    # (iii) gam, gam
    results.aiptw[k,5] <- aiptwATE(
      out = sim.df$outcome, trt = sim.df$treat, mu0 = out.gam$mu0,
      mu1 = out.gam$mu1, pi.hat = weight.gam$pi.hat
    )

    if (car.fit){

      # (iv) car, balancing (cvr)
      results.aiptw[k,7] <- aiptwATE(
        out = sim.df$outcome, trt = sim.df$treat, mu0 = out.car$mu0,
        mu1 = out.car$mu1, pi.hat = spb.cvr$pi.hat
      )

      # (v) car, balancing (min)
      results.aiptw[k,8] <- aiptwATE(
        out = sim.df$outcome, trt = sim.df$treat, mu0 = out.car$mu0,
        mu1 = out.car$mu1, pi.hat = as.vector(spb.min$pi.hat)
      )

      # (vi) car, gam
      results.aiptw[k,9] <- aiptwATE(
        out = sim.df$outcome, trt = sim.df$treat, mu0 = out.car$mu0,
        mu1 = out.car$mu1, pi.hat = weight.gam$pi.hat
      )
    }

  }

  # combine results
  results <- list(
    iptw = results.iptw,
    aiptw = results.aiptw
  )

  # return results
  results
}

#------------------------------------------------------------------------#
#--- 3. Perform simulation study                                      ---#
#------------------------------------------------------------------------#

res <- tecIPTW(
  n.sims = n.sims,
  sim.start = random.seed,
  nonlin = nonlinear,
  data.args = list(grid = grid.size, cov = cov.type, r.type = r.type, tau = 0.1),
  car.fit = FALSE,
  sparse = sparse
)

#------------------------------------------------------------------------#
#--- 4. Save results                                                  ---#
#------------------------------------------------------------------------#

if (sparse > 0){
  file.name <- paste("TS_", cov.type, "-r", r.type, "-n", grid.size, "-nl_", nonlin,
    "sp", sparse, "-s", random.seed + 1,"-f", random.seed + n.sims, ".RDS", sep = "")
} else {
  file.name <- paste("TS_", cov.type, "-r", r.type, "-n", grid.size, "-nl_", nonlin,
          "-s", random.seed + 1,"-f", random.seed + n.sims, ".RDS", sep = "")
}

# save results
saveRDS(res, here::here("output", "dgp1", file.name))

