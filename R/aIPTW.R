### aIPTW.R
### Nathan Wikle

### Develop augmented IPTW estimators of the treatment effect.

#---------------------------------------------------------#
#--- 1. Estimate ATE using outcome regression model(s) ---#
#---------------------------------------------------------#

outcomeATE <- function(formula, data = list(), trt.var, type = "gam", keep.fit = TRUE, ...){
  # Estimate the ATE using outcome regression adjustment.
  # Input
  #   formula: formula object passed to 'gam', 'bayesCAR', or 'spmodel'
  #   data: data frame used in the analysis
  #   trt.var: column name of treatment variable
  #   type: one of 'gam', 'car', or 'spmodel'
  #   keep.fit: Boolean indicating if the fitted model should be returned
  # Output
  #   A list containing the ATE, E(Y(1)), and E(Y(0)) estimates, as well as
  # the type of model fitted and the fitted model object (if keep.fit = TRUE).

  # set model fit to NULL (will change if keep.fit == TRUE)
  model.fit <- NULL

  if (type == "gam"){

    # fitted gam model
    gam.fit <- gam(formula, data = data, ...)

    # fitted model with treatment
    trt.data <- data
    trt.data[,trt.var] <- 1
    mu1 <- predict(gam.fit, newdata = trt.data)
    # fitted model without treatment
    ctrl.data <- data
    ctrl.data[,trt.var] <- 0
    mu0 <- predict(gam.fit, newdata = ctrl.data)

    # estimated ATE
    ate.est <- mean(mu1 - mu0)

    # fitted model
    if (keep.fit){
      model.fit <- gam.fit
    }

  } else if (type == "car"){

    car.fit <- bayesCAR(formula = formula, data = data, ...)

    # fitted model with treatment
    trt.data <- data
    trt.data[,trt.var] <- 1
    mu1 <- predictCAR(car.fit, newdata = trt.data)
    # fitted model without treatment
    ctrl.data <- data
    ctrl.data[,trt.var] <- 0
    mu0 <- predictCAR(car.fit, newdata = ctrl.data)

    # estimated ATE
    ate.est <- mean(mu1 - mu0)

    # fitted model
    if (keep.fit){
      model.fit <- car.fit
    }
  } else if (type == "spmodel"){

    # fitted 'splm' object from spmodel package
    splm.fit <- splm(formula = formula, data = data, ...)
    # fitted model with treatment (ignore additive RE)
    trt.data <- data
    trt.data[,trt.var] <- 1
    X.m1 <- model.matrix(formula, data = trt.data)
    mu1 <- X.m1 %*% splm.fit$coefficients$fixed
    # fitted model without treatment (ignore additive RE)
    ctrl.data <- data
    ctrl.data[,trt.var] <- 0
    X.m0 <- model.matrix(formula, data = ctrl.data)
    mu0 <- X.m0 %*% splm.fit$coefficients$fixed

    # estimated ATE
    ate.est <- mean(mu1 - mu0)
  }

  # return results
  list(
    ate = ate.est,
    mu1 = mu1,
    mu0 = mu0,
    type = type,
    model.fit = model.fit
  )
}

#------------------------------------------------------#
#--- 2. Estimate ATE using augmented IPTW estimator ---#
#------------------------------------------------------#

aiptw <- function(
    form.out, form.trt, data = list(), trt.var,
    type.out = "gam", type.trt = "gam", bal.input = NULL,
    ...
){
  # Estimate the ATE using an augmented IPTW (aIPTW) estimator.
  # Input
  #   form.out: formula object for outcome regression model
  #   form.trt: formula object for propensity score model
  #   data: data frame used in the analysis
  #   trt.var: column name of treatment variable
  #   type.out: type of OR to fit ('gam', 'car', or 'spmodel')
  #   type.trt: type of PS model to fit ('gam' or 'bal')
  #   bal.input: a list with input given to spBalance.
  # Output
  #   A list containing the ATE, as well as the fitted outcome and PS
  # regression models.


  ### 1. Fitted outcome model
  if (type.out == "gam"){
    # fitted gam model
    out.fit <- gam(form.out, data = data, ...)
    # estimate y(1)
    trt.data <- data
    trt.data[,trt.var] <- 1
    mu1 <- predict(out.fit, newdata = trt.data)
    # estimate y(0)
    ctrl.data <- data
    ctrl.data[,trt.var] <- 0
    mu0 <- predict(out.fit, newdata = ctrl.data)

  } else if (type.out == "car"){
    # fitted CAR model
    out.fit <- bayesCAR(formula = form.out, data = data, ...)
    # estimate y(1)
    trt.data <- data
    trt.data[,trt.var] <- 1
    mu1 <- predictCAR(out.fit, newdata = trt.data)
    # estimate y(0)
    ctrl.data <- data
    ctrl.data[,trt.var] <- 0
    mu0 <- predictCAR(out.fit, newdata = ctrl.data)

  } else if (type.out == "spmodel") {
    # fitted splm object
    out.fit <- splm(formula = formula, data = data, ...)

    # determine if local fit was used
    if (length(unique(out.fit$local_index)) > 1){
      local.fit <- TRUE
    } else {
      local.fit <- FALSE
    }

    # estimate y(1)
    trt.data <- data
    trt.data[,trt.var] <- 1
    mu1 <- predict(object = out.fit, newdata = trt.data, local = local.fit)

    # check that predict function returned real values
    if (sum(is.na(mu1)) > 0){
      stop("predict.spmodel returned NaN values.")
    }

    # estimate y(0)
    ctrl.data <- data
    ctrl.data[,trt.var] <- 0
    mu0 <- predict(object = out.fit, newdata = ctrl.data, local = local.fit)

    # check that predict function returned real values
    if (sum(is.na(mu0)) > 0){
      stop("predict.spmodel returned NaN values.")
    }

  } else {
    stop("type.out must be one of 'gam', 'car', or 'spmodel'.")
  }

  ### 2. Fitted propensity score model

  if (type.trt == "gam"){
    # estimate propensity score via mgcv::gam
    trt.fit <- mgcv::gam(form.trt, family = binomial, data = data)
    # gam.check(gam.fit)
    pi.hat <- trt.fit$fitted.values
  } else if (type.trt == "bal"){
    # balancing score model
    trt.fit <- spBalance(
      formula = form.trt,
      data = data,
      lambda = bal.input$lambda,
      tuning = bal.input$tuning,
      hide.details = TRUE,
      opt.params = bal.input$opt.params
    )
    # fitted propensity score
    pi.hat <- trt.fit$pi.hat
  } else {
    stop("type.trt must be one of 'gam' or 'bal'.")
  }

  ### 3. ATE estimate

  outcome <- data[, all.vars(form.out)[1]]
  treat <- data[, trt.var]
  ate.est <- aiptwATE(out = outcome, trt = treat, mu0 = mu0, mu1 = mu1, pi.hat = pi.hat)

  ### 4. return results
  results <- list(
    ate = ate.est,
    out.fit = out.fit,
    out.type = type.out,
    out.form = form.out,
    trt.fit = trt.fit,
    trt.type = type.trt,
    trt.form = form.trt
  )
  results
}

aiptwATE <- function(out, trt, mu0, mu1, pi.hat){
  # Calculates the aIPTW estimate using specified nuisance function
  #   estimates.
  # Input
  #   out: outcome vector
  #   trt: treatment vector
  #   mu0: estimates of E(Y | A = 0, X = x)
  #   mu1: estimates of E(Y | A = 1, X = x)
  #   pi.hat: PS estimate
  # Output
  #   Estimated ATE

  # calculate aiptw estimate
  ate.est <- mean((trt * out - (trt - pi.hat) * mu1) / pi.hat) -
    mean(((1 - trt) * out + (trt - pi.hat) * mu0) / (1 - pi.hat))
  ate.est
}

#-------------------------------------------------------#
#--- 3. CAR Bayesian LM code                         ---#
#-------------------------------------------------------#

sampleMVN <- function(b, Q, sum.to.zero = FALSE){
    # Generate a sample from MVN(mean = Q^{-1} * b, Var = Q^{-1}).
    # Input
    # b: transformed mean vector, i.e., mu = Q^{-1} %*% b
    # Q: precision matrix
    # sum.to.zero: Boolean indicating if a sum-to-zero constraint has
    #   been imposed on the random vector.
    # Output
    #   Samples from MVN using Algorithm 2.5 or 2.6 (sum-to-zero) of
    #   Rue and Held (2005).

    if (sum.to.zero){
      # Algorithm 2.6 of Rue and Held
      L.t <- chol(Q); L <- t(L.t); one.vec <- rep(1, nrow(Q))
      w <- solve(L, b)
      mu <- solve(L.t, w)
      z <- rnorm(n = nrow(Q))
      v <- solve(L.t, z)
      x <- mu[,1] + v

      V.vec <- solve(L.t, solve(L, one.vec))
      W.val <- as.numeric(crossprod(one.vec, V.vec))
      U.vec <- t(V.vec) / W.val

      c.val <- as.numeric(crossprod(one.vec, x))
      x.star <- x - (U.vec[1,] * c.val)
    } else {
      # Algorithm 2.5 of Rue and Held
      L.t <- chol(Q); L <- t(L.t)
      w <- solve(L, b)
      mu <- solve(L.t, w)
      z <- rnorm(n = nrow(Q))
      v <- solve(L.t, z)
      x.star <- mu + v
    }
  return(as.vector(x.star))
}

dMVN <- function(x, mu, Q, log = TRUE){
  # MVN likelihood evaluation.
  # Input
  #   x: random vector
  #   mu: mean vector
  #   Q: precision matrix
  #   log: Boolean indicating if log-likelihood should be returned.
  # Output
  #   MVN likelihood.

  n <- nrow(Q)
  L.t <- chol(Q)
  x.m.mu <- x - mu
  log.dens <- as.numeric(-n * log(2 * pi) / 2 +
      sum(log(diag(L.t))) - crossprod(x.m.mu, (Q %*% x.m.mu)) / 2)

  if (log){
    return(log.dens)
  } else {
    return(exp(log.dens))
  }
}

dMVN.tau2 <- function(x, mu, Q, tau2, log = TRUE){
  # MVN likelihood evaluation to be used when updating tau^2.
  #   Note that this avoids the evaluation of determ(Q), which is useful
  #   when using an intrinsic CAR prior for eta.
  # Input
  #   x: random vector
  #   mu: mean vector
  #   Q: precision matrix
  #   tau2: value of tau^2 (note: Q = Q.tilde / tau^2)
  #   log: Boolean indicating if log-likelihood should be returned.
  # Output
  #   MVN likelihood.

  n <- nrow(Q)
  x.m.mu <- x - mu

  log.dens <- as.numeric(-n * log(2 * pi) - (n * log(tau2)) -
      crossprod(x.m.mu, (Q %*% x.m.mu))) / 2

  if (log){
    return(log.dens)
  } else{
    return(exp(log.dens))
  }
}

# Parse formula with CAR term.
parseCAR <- function(form){
  # Parses formula object with a CAR term and returns CAR ojbect in a list
  # Input
  #   form: formual object containing 'car' term
  # Output
  #   A list with parsed formula terms

  # specials attribute indicates which terms are smooth / krr
  tf <- terms.formula(form, specials=c("car"))
  # labels of the model terms
  terms <- attr(tf, "term.labels")

  nt <- length(terms) # how many terms?
  if (attr(tf,"response") > 0){  # start the parametric model formula
    pf <- paste(as.character(attr(tf, "variables")[2]), "~", sep = "")
  } else {
    pf <- "~"
  }

  ### tells you which formula terms match "car"
  carp <- attr(tf, "specials")$car  # array of indices of car terms

  ### deals with offset (not applicable to balancing weights)
  off <- attr(tf, "offset") # location of offset in formula
  if (!is.null(off)){
    # have to remove the offset from this index list
    carp[carp > off] <- carp[carp > off] - 1
  }

  nc <- length(carp) # number of smooths
  kc <- 1; kp <- 1 # counters for terms in the formulae

  if (nt){
    for (i in 1:nt){
      if (kc <= nc && carp[kc] == i + 1){
        # it's a smooth
        car.obj <- eval(parse(text = terms[i]))
        kc <- kc + 1
      } else {
        # parametric term
        if (kp > 1) {
          # add to parametric formula
          pf <- paste(pf, "+", terms[i], sep = "")
        } else {
          pf <- paste(pf, terms[i], sep = "")
        }
        kp <- kp + 1
      }
    }
  }

  if (nc == nt){
    pf <- paste(pf, "1", sep = "")
  } else if (nc == 0){
    car.obj <- NULL
  }

  # return parsed formula
  ret <- list(
    pftext = pf,
    pf = as.formula(pf), # parametric forumla
    car.obj = car.obj # car object
  )
  class(ret) <- "split.formula"
  ret
}

car <- function(..., grid = "FALSE", A = NULL, rook = TRUE, icar = FALSE){
  # Format the CAR term in a formula into a car.spec object.
  # Input
  #   ... : car terms returned by 'parseCAR'
  #   grid: whether the areal data is defined on a grid
  #   A: adjacency matrix
  #   rook: type of CAR object to create (define connections via rook movements)
  #   icar: whether to use intrinsic CAR
  # Output
  #   A processed CAR object ('car.spec').

  # Format CAR term into a 'car.spec' object.
  if (grid){
    # create adjacency matrix
    vars <- as.list(substitute(list(...)))[-1]
    d <- length(vars)
    if (d != 2){
      stop("Input must include two named coordinate variables OR an adjacency matrix.")
    } else {
      coords <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
      if (coords[1] == "."){
        stop("Input must include two named coordinate variables OR an adjacency matrix.")
      }
      coords[2] <- deparse(vars[[2]], backtick = TRUE, width.cutoff = 500)
      if (coords[2] == "."){
        stop("Input must include two named coordinate variables OR an adjacency matrix.")
      }
    }
  } else {
    if (is.null(A)){
      stop("Input must include two named coordinate variables OR an adjacency matrix.")
    }
    coords <- NULL
  }
  # return processed CAR object
  res <- list(
    coords = coords,
    A = A,
    rook = rook,
    icar = icar
  )
  class(res) <- paste("car.spec", sep = "")
  res
}

adjMat <- function(coords, adj.type = "rook", sparse = TRUE){
  # Create an adjacency matrix for a given set of gridded coordinates.
  # Input
  #   coords: matrix with x,y coordinates
  #   adj.type: specifies how connections are defined; default is 'rook'
  #   sparse: whether to use sparse matrices for adjacency matrix
  # Output
  #   An adjacency matrix, A.

  # number of x coordinates
  n.x <- length(unique(coords[,1]))
  # number of y coordinates
  n.y <- length(unique(coords[,2]))
  # adjacency matrix
  if (sparse){
    A <- Matrix(0, nrow = n.x*n.y, ncol = n.x*n.y, sparse = TRUE)
  } else {
    A <- matrix(0, n.x*n.y, n.x*n.y)
  }

  for (j in 1:n.y){
    for (i in 1:n.x){
      # current position
      k <- (j - 1) * n.x + i

      # adjaceny positions
      if (adj.type == "rook"){
        if (i == 1){
          x.neigh <- c(k + 1)
        } else if (i == n.x){
          x.neigh <- c(k-1)
        } else {
          x.neigh <- c(k - 1, k + 1)
        }
        if (j == 1){
          y.neigh <- c(k + n.x)
        } else if (j == n.y){
          y.neigh <- c(k - n.x)
        } else {
          y.neigh <- c(k - n.x, k + n.x)
        }
        k.neigh <- c(x.neigh, y.neigh)
        # modify A matrix
        A[k, k.neigh] <- 1
        A[k.neigh, k] <- 1
      }
    }
  }
  A
}

carMatrices <- function(car.obj, data){
  # Create sparse matrices used to construct CAR precision matrix.
  # Input
  #   car.obj: a car.spec object returned by 'car'
  #   data: data frame
  # Output
  #   Matrices used to define CAR precision.

  # adjacency matrix
  if (is.null(car.obj$A)){
    # coordinates
    coord.mat <- data[,car.obj$coords]
    # create  adjacency matrix
    if (car.obj$rook){
      A.m <- adjMat(coord.mat, adj.type = "rook", sparse = TRUE)
    }
  } else {
    # adjacency matrix already included
    A.m <- car.obj$A
  }

  # diagonal matrix
  M <- Diagonal(n = nrow(A.m), x = rowSums(A.m))

  # return CAR matrices
  list(
    M = M,
    W = A.m
  )
}

### Priors for CAR model

# sigma^2 prior: half-normal prior distribution
#   -> corresponding to zero-mean normal with var = 100)
prior.s2 <- function(s2, theta = 0.1 * sqrt(pi/2)) {
  dhalfnorm(s2, theta, log = TRUE)
}

# tau^2 prior: half-normal prior distribution
#   -> corresponding to zero-mean normal with var = 100)
prior.tau2 <- function(tau2, theta = 0.1 * sqrt(pi/2)) {
  dhalfnorm(tau2, theta, log = TRUE)
}

# rho prior: beta prior distribution
prior.rho <- function(rho, shape1 = 18, shape2 = 2) {
  dbeta(rho, shape1, shape2, log = TRUE)
}

# beta prior: normal prior distribution
prior.beta <- function(beta, s2, kappa){
  sum(dnorm(beta, mean = 0, sd = sqrt(s2 / kappa), log = TRUE))
}

acceptance <- function(samps){
  # Calculate acceptance rate for given set of MCMC samples
  # Input
  #   samps: a vector or matrix of MCMC samples
  # Output
  #   MCMC acceptance rate.

  if(is.matrix(samps)){
    s.vec <- samps[,1]
  } else{
    s.vec <- samps
  }
  1 - sum(s.vec[-1] == s.vec[-length(s.vec)]) / (length(s.vec) - 1)
}

betaOrig <- function(beta, X.m){
  # Convert beta estimates to scale of unstandardized X matrix.
  # Input
  #   beta: coefficient estimates for standardized X matrix
  #   X.m: standardized X matrix
  # Output
  #   Coefficient estimates on the original scale of X

  # transform beta back to original scale of X
  if (is.matrix(X.m)){
    beta.new <- beta[-1] / apply(X.m, 2, sd)
    beta0.new <- beta[1] - beta.new %*% colMeans(X.m)
  } else {
    beta.new <- beta[-1] / sd(X.m)
    beta0.new <- beta[1] - beta.new * mean(X.m)
  }
  c(beta0.new, beta.new)
}

bayesCAR <- function(
  formula, data = list(), n.post = 1000, save.iter = 100, priors,
  s2.prop, tau2.prop, rho.prop, var.init = NULL, A = NULL
){
  # Metropolis-Hastings sampler for Bayesian LMM with CAR prior.
  # Input
  #   formula: a formula with parametric and 'car' term(s)
  #   n.post: number of posterior samples
  #   save.iter: how often to save spatial random effect and residual terms
  #     - default = every 100 iterations
  #   s2.prop: standard deviation of s2's proposal distribution
  #   tau2.prop: standard deviation of tau2's proposal distribution
  #   rho.prop: standard deviation of rho's proposal distribution
  #   var.init: list with initial values for
  #     (i) beta, (ii) s2, (iii) tau2, and (iv) rho.
  #     Notes:
  #       - var.init = NULL -> initial values are assigned to
  #           pre-specified defaults.
  #   A: adjacency matrix for CAR prior
  #     - default = NULL; use coordinates in car formula to create A
  # Output
  #   Return posterior samples as a list, containing:
  #       (i) beta, (ii) s2, (iii) tau2, (iv) rho, (v) beta.cov,
  #       (vi) eta and (vii) resids.

  ### 1. Parse formula
  parseform <- parseCAR(formula)

  ### 2. Create y vector; X matrix; CAR adjacency matrices

  # outcome data
  y.vec <- data[,all.vars(parseform$pf)[1]]
  # design matrix
  X <- model.matrix(parseform$pf, data)
  # adjacency matrices
  car.mats <- carMatrices(parseform$car.obj, data)

  ### 3. Initialize parameters

  # dimension of beta vector
  n.p <- ncol(X)
  # number of observations
  n.y <- length(y.vec)

  # standardize columns of predictor matrix and add on column of ones
  pred.scale <- scale(X[,-1])
  X.full <- cbind(X[,1], pred.scale)
  colnames(X.full) <- colnames(X)

  if (is.null(var.init)){

    # create initial parameter estimates (from OLS)
    lm.fit <- lm(parseform$pf, data)
   # beta.init <- coef(lm.fit)
    s2.k <- summary(lm.fit)$sigma^2
    # re-adjust beta
   # beta.k <- beta.init
    # initial tau and rho parameterization
    tau2.k <- 1
    if (parseform$car.obj$icar){
      rho.k <- 1
    } else {
      rho.k <- 0.95
    }

  } else {
    # initial parameter values supplied by var.init
    beta.k <- var.init$beta
    s2.k <- var.init$s2
    tau2.k <- var.init$tau2
    rho.k <- var.init$rho
  }

  # car precision matrix
  Q.k <- (car.mats$M - rho.k * car.mats$W) / tau2.k

  # initialize eta to 0 and then initialize mean
  eta.k <- rep(0, n.y)
  # mu.k <- as.vector(X.full %*% beta.k + eta.k)

  ### 4. Set up prior/posterior data structures
  post.samps <- matrix(NA_real_, nrow = n.post, ncol = n.p + 3)
  eta.samps <- matrix(NA_real_, nrow = n.post / save.iter, ncol = n.y)
  resid.samps <- matrix(NA_real_, nrow = n.post / save.iter, ncol = n.y)

  ### 5. Generate posterior samples
  for (k in 1:n.post){

    # (a) Sample from beta full conditional

    # sample beta.k
    # Q.beta <- (crossprod(X.full) + priors$beta * diag(ncol(X.full))) / s2.k
    Q.beta <- (crossprod(X.full) + s2.k / priors$beta * diag(ncol(X.full))) / s2.k
    b.beta <- crossprod(X.full, y.vec - eta.k) / s2.k
    beta.k <- sampleMVN(b = b.beta, Q = Q.beta)

    # adjust mean for new beta.k
    mu.k <- as.vector(X.full %*% beta.k + eta.k)

    # (b) Posterior sample of s^2
    alpha.s2 <- n.y + priors$s2[1]
    beta.s2 <- as.numeric(crossprod(y.vec - mu.k) / 2) + priors$s2[2]
    s2.k <- rgamma(n = 1, shape = alpha.s2, rate = beta.s2)

    # (c) Joint sample of spatial random effect & tau^2 from full conditional
    # Based on algorithm 2.5 on pg. 46 of Rue & Held

    #------------------------------------------------------------#
    #--- Block update of tau2, rho, and eta.                  ---#
    #------------------------------------------------------------#

    # (i) proposal from normal for tau^2
    tau2.star <- rnorm(n = 1, mean = tau2.k, sd = tau2.prop)

    # check that tau2 is in region of appropriate support
    if (tau2.star > 0){

      # (ii) Propose rho (if rho < 1; otherwise, hold constant)
      if (rho.k < 1){
        # Proposal from normal
        rho.star <- rtruncnorm(n = 1, a = 0, b = 1, mean = rho.k, sd = rho.prop)
        # Update Q
        Q.star <- (car.mats$M - rho.star * car.mats$W) / tau2.star
        # Q.star <- solve(M.plus, (I.mat - rho.star * W.plus)) / tau2.star
      } else {
        # Update Q
        Q.star <- (Q.k * tau2.k) / tau2.star
      }

      # Used to sample eta.star and evaluate full conditional density
      Q.eta.star <- Q.star + (Diagonal(nrow(Q.star)) / s2.k)
      b.eta.star <- (y.vec - (X.full %*% beta.k)) / s2.k
      mu.full.star <- solve(Q.eta.star, b.eta.star) # mean of eta full conditional

      if (rho.k < 1){
        # CAR prior
        eta.star <- sampleMVN(b = b.eta.star, Q = Q.eta.star, sum.to.zero = FALSE)
      } else {
        # ICAR prior with sum-to-zero constraint
        eta.star <- sampleMVN(b = b.eta.star, Q = Q.eta.star, sum.to.zero = TRUE)
      }

      # Update model mean
      mu.star <- as.vector(X.full %*% beta.k + eta.star)

      # used in evaluation of etafull conditional
      Q.eta.k <- Q.k + (Diagonal(nrow(Q.k)) / s2.k)
      b.eta.k <- (y.vec - (X.full %*% beta.k)) / s2.k
      mu.full.k <- solve(Q.eta.k, b.eta.k)

      # acceptance rate
      if (rho.k < 1){
        # includes prior for rho
        accept.num <- sum(dnorm(y.vec, mean = mu.star, sd = sqrt(s2.k), log = TRUE)) +
          dMVN(eta.star, mu = rep(0, n.y), Q.star) +
          prior.rho(rho.star, shape1 = priors$rho[1], shape2 = priors$rho[1]) +
          prior.tau2(tau2.star, theta = priors$tau2) +
          dMVN(eta.k, mu.full.k, Q.eta.k, log = TRUE) +
          log(dtruncnorm(rho.k, a = 0, b = 1, mean = rho.star, sd = rho.prop))

        accept.denom <- sum(dnorm(y.vec, mean = mu.k, sd = sqrt(s2.k), log = TRUE)) +
          dMVN(eta.k, mu = rep(0, n.y), Q.k) +
          prior.rho(rho.k, shape1 = priors$rho[1], shape2 = priors$rho[1]) +
          prior.tau2(tau2.k, theta = priors$tau2) +
          dMVN(eta.star, mu.full.star, Q.eta.star, log = TRUE) +
          log(dtruncnorm(rho.star, a = 0, b = 1, mean = rho.k, sd = rho.prop))

      } else {

        # no need to update rho
        accept.num <- sum(dnorm(y.vec, mean = mu.star, sd = sqrt(s2.k), log = TRUE)) +
          dMVN.tau2(eta.star, mu = rep(0, n.y), Q = Q.star, tau2 = tau2.star) +
          prior.tau2(tau2.star, theta = priors$tau2) +
          dMVN(eta.k, mu.full.k, Q.eta.k, log = TRUE)

        accept.denom <- sum(dnorm(y.vec, mean = mu.k, sd = sqrt(s2.k), log = TRUE)) +
          dMVN.tau2(eta.k, mu = rep(0, n.y), Q = Q.k, tau2 = tau2.k) +
          prior.tau2(tau2.k, theta = priors$tau2) +
          dMVN(eta.star, mu.full.star, Q.eta.star, log = TRUE)
      }

      # log(acceptance ratio)
      accept.tau2 <- accept.num - accept.denom

      # Decide whether to accept proposals or not
      if (log(runif(1)) < accept.tau2){

        tau2.k <- tau2.star
        Q.k <- Q.star
        mu.k <- mu.star
        eta.k <- eta.star

        if (rho.k < 1){
          rho.k <- rho.star
        }
      }
    }

    # (d) Save results- all beta and s^2 values and every 100th y vector
    post.samps[k,] <- c(beta.k, s2.k, tau2.k, rho.k)
    if (k %% save.iter == 0) {
      eta.samps[k/save.iter,] <- as.vector(eta.k)
      resid.samps[k/save.iter,] <- as.vector(y.vec - (X.full %*% beta.k) - eta.k)
    }
  }

  ### 6. Prepare function output

  # hyperparameters & sample covariance
  s2.post <- post.samps[, n.p + 1]
  tau2.post <- post.samps[, n.p + 2]
  rho.post <- post.samps[, n.p + 3]
  beta.cov <- cov(post.samps[, c(1:n.p)])
  # transform beta back to original scale
  beta.post <- t(apply(post.samps[, c(1:n.p)], 1, betaOrig, X.m = X[,-1]))

  ### 7. Return posterior samples
  list(
    beta = beta.post,
    s2 = s2.post,
    tau2 = tau2.post,
    rho = rho.post,
    beta.cov = beta.cov,
    y = y.vec,
    eta = eta.samps,
    resids = resid.samps,
    save.iter = save.iter,
    fm = parseform
  )
}

predictCAR <- function(car.mcmc, newdata, burnin = 0, post.mean = TRUE, sp.effect = TRUE){
  # Make predictions using CAR output
  # Input
  #   car.mcmc: output from 'bayesCAR'
  #   newdata: data on which to make predictions
  #   burnin: number of samples to discard as burnin
  #   post.mean: Boolean indicating if only the mean prediction should be returned
  #   sp.effect: Boolean indicating if spatial random effect should be included
  #     in the prediction
  # Output
  #   Estimates of E(Y | newdata).


  if (is.matrix(car.mcmc$beta)){
    # design matrix
    X.new <- model.matrix(car.mcmc$fm$pf, newdata)
    # determine which beta values to keep
    keep.beta <- which(1:nrow(car.mcmc$beta) %% car.mcmc$save.iter == 0)

    if (sp.effect){
      # include spatial random effect in prediction
      mu.post <- tcrossprod(X.new, car.mcmc$beta[keep.beta,]) + t(car.mcmc$eta)
    } else {
      # fixed effects only
      mu.post <- tcrossprod(X.new, car.mcmc$beta)
    }
  } else {
    # design matrix
    X.new <- model.matrix(car.mcmc$fm$pf, newdata)
    # determine which beta values to keep
    keep.beta <- which(1:length(car.mcmc$beta) %% car.mcmc$save.iter == 0)

    if (sp.effect){
      # include spatial random effect in prediction
      mu.post <- tcrossprod(X.new, car.mcmc$beta[keep.beta]) + t(car.mcmc$eta)
    } else {
      # fixed effects only
      mu.post <- tcrossprod(X.new, car.mcmc$beta)
    }
  }

  if (post.mean){
    rowMeans(mu.post)
  } else {
    mu.post
  }
}




