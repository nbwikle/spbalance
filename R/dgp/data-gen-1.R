### data-gen-1.R
### Nathan Wikle

### Simulate from Tec et al.'s data generating process.

simDecomps <- function(type, method = "svd"){
  # Creates covariance matrices for fast simulation of GPs. Matrices are generated
  #   at grid sizes 10x10, 25x25, 50x50, and 100x100. Corresponding spectral
  #   decomps are saved in "./data/cov-mats/"
  # Input
  #   type: specifies the type of covariance matrix to generate; options include
  #     GP1 = squared exponential kernel with k(d) = exp(-(d / 5.1)^2);
  #       this is used in Tec et al. to generate "velocity vector fields"
  #     SE = squared exponential kernel with k(d) = exp(-(d / theta)^2);
  #       range param -> theta = 3/2, 5/2, and 10/2
  #     Expo = exponential kernel with k(d) = exp(- |d| / theta);
  #       range param -> theta = 3/2, 5/2, 10/2
  #     Matern = Matern kernel with range rho = 3, 5, 10; smoothness nu = 1.5
  #   method: matrix decomposition method; options include
  #     'eigen', 'svd', 'chol'
  # Output
  #   Matrix decompositions are saved in "./data/cov-mats/".

  if (type == "GP1"){

    # 14x14 grid
    my.pts <- expand.grid(
      x = seq(1, 14, by = 1), y = seq(1, 14, by = 1)
    )
    D <- as.matrix(dist(my.pts, diag = TRUE, upper = TRUE))
    K.mat <- covFun(D, theta = 5.1, type = "SE")
    R.mat <- decomp(K.mat, method)
    saveRDS(R.mat, file = here::here("data", "cov-mats", "GP1-n10.RDS"))

    # 29x29 grid
    my.pts <- expand.grid(
      x = seq(1, 29, by = 1), y = seq(1, 29, by = 1)
    )
    D <- as.matrix(dist(my.pts, diag = TRUE, upper = TRUE))
    K.mat <- covFun(D, theta = 5.1, type = "SE")
    R.mat <- decomp(K.mat, method)
    saveRDS(R.mat, file = here::here("data", "cov-mats", "GP1-n25.RDS"))

    # 54x54 grid
    my.pts <- expand.grid(
      x = seq(1, 54, by = 1), y = seq(1, 54, by = 1)
    )
    D <- as.matrix(dist(my.pts, diag = TRUE, upper = TRUE))
    K.mat <- covFun(D, theta = 5.1, type = "SE")
    R.mat <- decomp(K.mat, method)
    saveRDS(R.mat, file = here::here("data", "cov-mats", "GP1-n50.RDS"))

    # 104x104 grid
    my.pts <- expand.grid(
      x = seq(1, 104, by = 1), y = seq(1, 104, by = 1)
    )
    D <- as.matrix(dist(my.pts, diag = TRUE, upper = TRUE))
    K.mat <- covFun(D, theta = 5.1, type = "SE")
    R.mat <- decomp(K.mat, method)
    saveRDS(R.mat, file = here::here("data", "cov-mats", "GP1-n100.RDS"))

  } else if (type == "SE"){

    rho.vals <- c(3, 5, 10)

    for (k in 1:3){
      # 10x10 grid
      my.pts <- expand.grid(
        x = seq(1, 10, by = 1), y = seq(1, 10, by = 1)
      )
      D <- as.matrix(dist(my.pts, diag = TRUE, upper = TRUE))
      K.mat <- covFun(D, theta = rho.vals[k] / 2, type = "SE")
      R.mat <- decomp(K.mat, method)
      saveRDS(R.mat, file = here::here("data", "cov-mats",
        paste("SE-n10-r", k, ".RDS", sep = "")))
    }

    for (k in 1:3){
      # 25x25 grid
      my.pts <- expand.grid(
        x = seq(1, 25, by = 1), y = seq(1, 25, by = 1)
      )
      D <- as.matrix(dist(my.pts, diag = TRUE, upper = TRUE))
      K.mat <- covFun(D, theta = rho.vals[k] / 2, type = "SE")
      R.mat <- decomp(K.mat, method)
      saveRDS(R.mat, file = here::here("data", "cov-mats",
        paste("SE-n25-r", k, ".RDS", sep = "")))
    }

    # 50x50 grid
    for (k in 1:3){
      my.pts <- expand.grid(
        x = seq(1, 50, by = 1), y = seq(1, 50, by = 1)
      )
      D <- as.matrix(dist(my.pts, diag = TRUE, upper = TRUE))
      K.mat <- covFun(D, theta = rho.vals[k] / 2, type = "SE")
      R.mat <- decomp(K.mat, method)
      saveRDS(R.mat, file = here::here("data", "cov-mats",
        paste("SE-n50-r", k, ".RDS", sep = "")))
    }

    # 100x100 grid
    for (k in 1:3){
      my.pts <- expand.grid(
        x = seq(1, 100, by = 1), y = seq(1, 100, by = 1)
      )
      D <- as.matrix(dist(my.pts, diag = TRUE, upper = TRUE))
      K.mat <- covFun(D, theta = rho.vals[k] / 2, type = "SE")
      R.mat <- decomp(K.mat, method)
      saveRDS(R.mat, file = here::here("data", "cov-mats",
        paste("SE-n100-r", k, ".RDS", sep = "")))
    }

  } else if (type == "Expo"){

    rho.vals <- c(3, 5, 10)

    # 10x10 grid
    for (k in 1:3){
      my.pts <- expand.grid(
        x = seq(1, 10, by = 1), y = seq(1, 10, by = 1)
      )
      D <- as.matrix(dist(my.pts, diag = TRUE, upper = TRUE))
      K.mat <- covFun(D, theta = rho.vals[k] / 2, type = "Expo")
      R.mat <- decomp(K.mat, method)
      saveRDS(R.mat, file = here::here("data", "cov-mats",
        paste("Expo-n10-r", k, ".RDS", sep = "")))
    }

    # 25x25 grid
    for (k in 1:3){
      my.pts <- expand.grid(
        x = seq(1, 25, by = 1), y = seq(1, 25, by = 1)
      )
      D <- as.matrix(dist(my.pts, diag = TRUE, upper = TRUE))
      K.mat <- covFun(D, theta = rho.vals[k] / 2, type = "Expo")
      R.mat <- decomp(K.mat, method)
      saveRDS(R.mat, file = here::here("data", "cov-mats",
        paste("Expo-n25-r", k, ".RDS", sep = "")))
    }

    # 50x50 grid
    for (k in 1:3){
      my.pts <- expand.grid(
        x = seq(1, 50, by = 1), y = seq(1, 50, by = 1)
      )
      D <- as.matrix(dist(my.pts, diag = TRUE, upper = TRUE))
      K.mat <- covFun(D, theta = rho.vals[k] / 2, type = "Expo")
      R.mat <- decomp(K.mat, method)
      saveRDS(R.mat, file = here::here("data", "cov-mats",
        paste("Expo-n50-r", k, ".RDS", sep = "")))
    }

    # 100x100 grid
    for (k in 1:3){
      my.pts <- expand.grid(
        x = seq(1, 100, by = 1), y = seq(1, 100, by = 1)
      )
      D <- as.matrix(dist(my.pts, diag = TRUE, upper = TRUE))
      K.mat <- covFun(D, theta = rho.vals[k] / 2, type = "Expo")
      R.mat <- decomp(K.mat, method)
      saveRDS(R.mat, file = here::here("data", "cov-mats",
        paste("Expo-n100-r", k, ".RDS", sep = "")))
    }

  } else if (type == "Matern"){

    rho.vals <- c(3, 5, 10)

    # 10x10 grid
    for (k in 1:3){
      my.pts <- expand.grid(
        x = seq(1, 10, by = 1), y = seq(1, 10, by = 1)
      )
      D <- as.matrix(dist(my.pts, diag = TRUE, upper = TRUE))
      K.mat <- covFun(D, theta = c(rho.vals[k], 1.5), type = "Matern")
      R.mat <- decomp(K.mat, method)
      saveRDS(R.mat, file = here::here("data", "cov-mats",
        paste("Matern-n10-r", k, ".RDS", sep = "")))
    }

    # 25x25 grid
    for (k in 1:3){
      my.pts <- expand.grid(
        x = seq(1, 25, by = 1), y = seq(1, 25, by = 1)
      )
      D <- as.matrix(dist(my.pts, diag = TRUE, upper = TRUE))
      K.mat <- covFun(D, theta = c(rho.vals[k], 1.5), type = "Matern")
      R.mat <- decomp(K.mat, method)
      saveRDS(R.mat, file = here::here("data", "cov-mats",
        paste("Matern-n25-r", k, ".RDS", sep = "")))
    }

    # 50x50 grid
    for (k in 1:3){
      my.pts <- expand.grid(
        x = seq(1, 50, by = 1), y = seq(1, 50, by = 1)
      )
      D <- as.matrix(dist(my.pts, diag = TRUE, upper = TRUE))
      K.mat <- covFun(D, theta = c(rho.vals[k], 1.5), type = "Matern")
      R.mat <- decomp(K.mat, method)
      saveRDS(R.mat, file = here::here("data", "cov-mats",
        paste("Matern-n50-r", k, ".RDS", sep = "")))
    }

    # 100x100 grid
    for (k in 1:3){
      my.pts <- expand.grid(
        x = seq(1, 100, by = 1), y = seq(1, 100, by = 1)
      )
      D <- as.matrix(dist(my.pts, diag = TRUE, upper = TRUE))
      K.mat <- covFun(D, theta = c(rho.vals[k], 1.5), type = "Matern")
      R.mat <- decomp(K.mat, method)
      saveRDS(R.mat, file = here::here("data", "cov-mats",
        paste("Matern-n100-r", k, ".RDS", sep = "")))
    }
  }
}

covFun <- function(d, theta, type = "Matern"){
  # Evaluate covariance function at distance d.
  # Input:
  #   d: distance (scalar, vector, or matrix)
  #   theta: vector of covariance function parameters
  #   type: type of covariance function (one of 'Matern', 'SE', or 'Expo')
  # Output:
  #   Covariance (scalar, vector, or matrix)

  if (type == "Matern"){
    # Matern (using definition from Fuglstad et al.)
    nu <- theta[2] # smoothness parameter
    rho <- theta[1] # range parameter (correlation is ~ 0.1 at distance = rho)
    log.ratio <- (log(8) + log(nu)) / 2 + log(d) - log(rho)
    log.cov <- log(2^(1 - nu)) - lgamma(nu) + nu * log.ratio +
      log(besselK(exp(log.ratio), nu))
    log.cov[is.na(log.cov)] <- 0 # Corr(0) = 1 in the limit
  } else if (type == "SE"){
    # squared exponential
    log.cov <- -(d / theta[1])^2 / 2
  } else if (type == "Expo"){
    # exponential
    log.cov <- -abs(d) / theta[1]
  }
  # return covariance
  exp(log.cov)
}

decomp <- function(Sigma, method = "eigen"){
  # Performs a matrix decomposition using the available methods
  #   from mvtnorm::mvnorm.
  # Input
  #   Sigma: covariance matrix
  #   method: decomposition method; options include
  #     'eigen', 'svd', 'chol'
  # Output
  #   Covariance matrix decomposition.

  # create matrix R such that Sigma = R'R
  if (method == "eigen"){
    # eigendecomposition
    eS <- eigen(Sigma, symmetric = TRUE)
    # check that it's positive definite
    if (!all(eS$values >= -sqrt(.Machine$double.eps) * abs(eS$values[1]))){
      warning("Sigma is numerically not positive semidefinite")
    }
    t(eS$vectors %*% (t(eS$vectors) * sqrt(pmax(eS$values, 0))))
  } else if (method == "svd"){
    s. <- svd(Sigma)
    if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))){
      warning("Sigma is numerically not positive semidefinite")
    }
    t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
  } else if (method == "chol"){
    R <- chol(Sigma, pivot = TRUE)
    R[, order(attr(R, "pivot"))]
  }
}

sampleGP <- function(n = 1, R){
  # Generates samples from GP with Sigma = R'R
  # Input
  #   n: number of samples to generate
  #   R: decomposition of covariance matrix (Sigma = R'R)
  # Output
  #   Samples from a mean-zero GP with given covariance matrix.

  mvn.samps <- matrix(rnorm(n * ncol(R)), nrow = n, byrow = TRUE) %*%  R
  return(drop(mvn.samps))
}

dataGen1 <- function(sim.args, K.mat = "default", nonlin = FALSE){
  # Simulates data from the data generating process considered in Tec et al. (2023).
  # Input
  #   sim.args: list containing the following:
  #     (i) grid: specifies the dimension of a square grid (10, 25, 50, or 100)
  #     (ii) cov: type of covariance function ('SE', 'Expo', or 'Matern')
  #     (iii) r.type: specifies which spatial range parameter to use;
  #             must be one of 3, 5, 10.
  #     (iv) tau: magnitude of treatment effect
  #   K.mat: convolution matrix ("default" = matrix used in Tec et al.)
  #   nonlin: whether the "non-local" confounding is a linear or nonlinear
  #       function of the underlying gradient field (default = FALSE)
  # Output
  #   A list with the following simulated data:
  #     v1: covariate 1 (x-gradient of potential surface)
  #     v2: covariate 2 (y-gradient of potential surface)
  #     pot: potential surface
  #     y: outcome
  #     pi: propensity for treatment
  #     z: treatment status
  #     effect: true causal effect (tau)


  # 1. simulation set-up

  n.val <- sim.args$grid
  r.val <- sim.args$r.type
  sigma.type <- sim.args$cov

  # Sigma1 = t(R1) %*% R1
  R1 <- readRDS(here::here("data", "cov-mats",
          paste("GP1-n", n.val, ".RDS", sep = "")))

  # Sigma2 = t(R2) %*% R2
  R2 <- readRDS(here::here("data", "cov-mats",
          paste(sigma.type, "-n", n.val, "-r", r.val, ".RDS", sep = "")))


  # 2. sample from GP(0, Cov = Squared Exp, rho = 5.1)
  f.sim <- matrix(sampleGP(n = 1, R1), nrow = n.val + 4, ncol = n.val + 4)

  # 2. calculate gradient fields
  grad.x <- -rowDiffs(cbind(rep(0, nrow(f.sim)), f.sim)) # gradient in x direction
  grad.y <- colDiffs(rbind(rep(0, ncol(f.sim)), f.sim)) # gradient in y direction

  # 3. restrict to original dimensions
  orig.x <- 3:(n.val+2); orig.y <- 3:(n.val+2)

  gradx.f <- grad.x[orig.x, orig.y]
  grady.f <- grad.y[orig.x, orig.y]

  f.surf <- f.sim[orig.x, orig.y]

  # 4. convolve vector fields

  # convolution matrix
  if (K.mat == "default"){
    phi.upper <- 6

    if (phi.upper %% 2 == 0){
      phi.upper <- phi.upper + 1
    }

    K1 <- rbind(
      matrix(-1, nrow = (phi.upper - 1) / 2, ncol = phi.upper),
      rep(0, phi.upper),
      matrix(1, nrow = (phi.upper - 1) / 2, ncol = phi.upper)
    )

    K2 <- t(K1)
  }

  # perform convolution
  if (nonlin){
    x.conv <- scipy$signal$correlate2d(sign(gradx.f), K2, mode = 'same')
    y.conv <- scipy$signal$correlate2d(sign(grady.f), -K1, mode = 'same')
  } else {
    x.conv <- scipy$signal$correlate2d(gradx.f, K2, mode = 'same')
    y.conv <- scipy$signal$correlate2d(grady.f, -K1, mode = 'same')
  }

  conv.surf <- x.conv + y.conv

  # conv.surf <- rnorm(n = nrow(R2))

  # 5. treatment assignment

  # define mean surface
  logit <- (conv.surf - mean(conv.surf)) / sd(conv.surf)
  p.vals <-  exp(-log1p(exp(-logit)))
  # sample treatment assignment
  treat <- matrix(rbinom(n = length(p.vals), size = 1, prob = as.vector(p.vals)),
                  nrow = n.val, ncol = n.val)

  # 6. outcome

  # spatial random effect
  eta <- matrix(sampleGP(n = 1, R2), nrow = n.val, ncol = n.val)
  # local (measurement) error
  eps <- matrix(rnorm(n = length(p.vals)), nrow = n.val, ncol = n.val)
  # combined error terms
  error <-  0.5 * eta + 0.5 * eps

  # outcome
  y <- sqrt(0.5) * (-logit + error) + sim.args$tau * treat
  # y <- sqrt(0.5) * (-logit) + sim.args$tau * treat + 0.5 * error

  # return simulated data
  list(
    x1 = gradx.f, # covariate 1 (x gradient)
    x2 = grady.f, # covariate 2 (y gradient)
    pot = f.surf, # potential surface
    y = y, # simulated outcome
    pi = p.vals, # true probability of treatment
    z = treat, # simulated treatment assignment
    effect = sim.args$tau # effect size
  )
}

dataGen1.SpatRaster <- function(sim){
  # Reformat data from 'dataGen1' into raster layers.
  # Input
  #   sim: output (a list) from the simulation function, 'simData.Tec'
  # Output
  #   A SpatRaster object with the following layers:
  #     x1: covariate 1 (x gradient)
  #     x2: covariate 2 (y gradient)
  #     pot: potential surface
  #     y: simulated outcome
  #     pi: true probability of treatment
  #     z: treatment assignment

  ### create SpatRaster with the listed elements

  # dimensions
  n.rows <- nrow(sim[[1]])
  n.cols <- ncol(sim[[1]])

  raster.obj <- terra::rast(
    ncol = n.cols, nrow = n.rows,
    xmin = 0, xmax = n.cols,
    ymin = 0, ymax = n.rows
  )

  # save results in a raster
  values(raster.obj) <- sim[[1]]
  names(raster.obj)[1] <- "x1"
  raster.obj$x2 <- sim[[2]]
  raster.obj$pot <- sim[[3]]
  raster.obj$y <- sim[[4]]
  raster.obj$pi <- sim[[5]]
  raster.obj$z <- sim[[6]]

  # return raster
  return(raster.obj)
}


simDG1 <- function(
  data.args = list(grid = 50, cov = "Matern", r.type = 3, tau = 0.1),
  se.method = "all",
  sim.seed = 1,
  save.res = FALSE,
  p.known = TRUE
){

  ### 1. set seed for reproducibilty
  set.seed(sim.seed)

  ### 2. simulate data
  sim.data <- dataGen1(sim.args = data.args, K.mat = "default", nonlin = FALSE)

  # SpatRast object
  sim.rast <- dataGen1.SpatRaster(sim.data)

  # data frame
  sim.df <- data.frame(
    treat = values(sim.rast$z),
    outcome = values(sim.rast$y),
    pi = values(sim.rast$pi),
    v1 = values(sim.rast$x1),
    v2 = values(sim.rast$x2),
    crds(sim.rast$y)
  )
  colnames(sim.df) <- c("treat", "outcome", "pi", "v1", "v2", "x", "y")

  ### 3. Estimate treatment effect
  ate.hat <- ateEst(out = sim.df$outcome, trt = sim.df$treat, pr.trt = sim.df$pi)

  ### 4. Estimate standard error
  if (se.method == "true" | se.method == "all"){
    if (save.res){
      saveRDS(ate.hat$ate, here::here("output", "var-est",
        paste("GP1_ATE-", data.args[[2]], "-r", data.args[[3]],
              "-n", data.args[[1]], "-sim", sim.seed, ".RDS", sep = "")))
    }
  }
  if (se.method == "sandwich" | se.method == "all"){

    # sandwich estimate of standard error
    sand.est <- c(
      sandwichEst(ate.hat[[2]][1], z = sim.df$treat,
                  y = sim.df$outcome, p = sim.df$pi),
      sandwichEst(ate.hat[[2]][2], z = sim.df$treat,
                  y = sim.df$outcome, p = sim.df$pi, type = 2)
    )

    if (save.res){
      saveRDS(sand.est, here::here("output", "var-est",
        paste("GP1_sandwich-", data.args[[2]], "-r", data.args[[3]],
              "-n", data.args[[1]], "-sim", sim.seed, ".RDS", sep = "")))
    }

  }
  if (se.method == "window" | se.method == "all"){

    # window subsampling estimate of standard error
    window.est <- windowSE(
      theta.hat = ate.hat[[1]][1:2,],
      out = sim.df$outcome,
      trt = sim.df$treat,
      p.trt = sim.df$pi,
      data = sim.df,
      coord.names = c("x","y"),
      window = "square",
      ipw = c("ht", "haj")
    )

    if (save.res){
      saveRDS(sand.est, here::here("output", "var-est",
        paste("GP1_window-", data.args[[2]], "-r", data.args[[3]],
              "-n", data.args[[1]], "-sim", sim.seed, ".RDS", sep = "")))
    }
  }
  if (se.method == "bootstrap" | se.method == "all"){

    # generate bootstraps
    bts.prelim <- spbbGrid(
      data = sim.df,
      coord.names = c("x", "y"),
      n.boot = 100
    )

    # determine optimal block size
    bts.se1 <- bootSE(
      data = sim.df,
      bootstraps = bts.prelim$bootstraps[[1]],
      coords = bts.prelim$coords,
      coord.names = c('x', 'y')
    )

    bts.se2 <- bootSE(
      data = sim.df,
      bootstraps = bts.prelim$bootstraps[[2]],
      coords = bts.prelim$coords,
      coord.names = c('x', 'y')
    )

    bts.se3 <- bootSE(
      data = sim.df,
      bootstraps = bts.prelim$bootstraps[[3]],
      coords = bts.prelim$coords,
      coord.names = c('x', 'y')
    )

    # determine the optimal block size
    bts.bias <-  2 * ceiling(0.5 * nrow(sim.df)^(1 / 6)) * (bts.se2 - bts.se3)
    bsize.nol <- (bts.bias^2 * nrow(sim.df) / 2 / bts.se1^2)^(1 / 4)
    bsize.ol <- bsize.nol * sqrt(3/2)

    # redo bootstrap with correct block size
    bts.final <- spbbGrid(
      data = sim.df,
      coord.names = c("x", "y"),
      n.boot = 100,
      block.l = unique(ceiling(bsize.ol))
    )

    bsize.ceiling <- ceiling(bsize.ol)
    bsize.unique <- unique(bsize.ceiling)
    se.bootstrap <- rep(NA_real_, length = length(bts.se3))

    for (opt in 1:length(bsize.unique)){

      bt.opt <- bootSE(
        data = sim.df,
        bootstraps = bts.final$bootstraps[[opt]],
        coords = bts.final$coords,
        coord.names = c('x', 'y')
      )

      # determine which elements to save
      save.opt <- which(bsize.ceiling == bsize.unique[opt])
      se.bootstrap[save.opt] <- bt.opt[save.opt]
    }

    if (save.res){
      saveRDS(se.bootstrap, here::here("output", "var-est",
        paste("GP1_bootstrap-", data.args[[2]], "-r", data.args[[3]],
              "-n", data.args[[1]], "-sim", sim.seed, ".RDS", sep = "")))
    }
  }

  ### 5. Return results
  if (se.method == "all"){
    final.results <- list(
      ate = ate.hat,
      sandwich = sand.est,
      window = window.est,
      bootstrap = se.bootstrap
    )
  } else if (se.method == "true"){
    final.results <- list(
      ate = ate.hat
    )
  } else if (se.method == "sandwich"){
    final.results <- list(
      ate = ate.hat,
      sandwich = sand.est
    )
  } else if (se.method == "window"){
    final.results <- list(
      ate = ate.hat,
      window = window.est
    )
  } else if (se.method == "bootstrap"){
    final.results <- list(
      ate = ate.hat,
      sandwich = sand.est
    )
  }

  return(final.results)
}


# ### Example
#
# # Simulates data with given grid size and covariance type.
# set.seed(75)
# sim.data <- dataGen1(
#   sim.args = list(grid = 50, cov = "Matern", r.type = 1, tau = 0.1),
#   K.mat = "default", nonlin = FALSE
# )
#
# # Convert to raster
# sim.rast <- dataGen1.SpatRaster(sim.data)
# plot(sim.rast)


