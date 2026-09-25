### hac-var-est.R
### Nathan Wikle

# ----------------------------------------------#
# --- required packages                      ---#
# ----------------------------------------------#

# library(fields)
# library(gstat)
# library(conleyreg)
# library(mgcv)

# ----------------------------------------------#
# --- functions                              ---#
# ----------------------------------------------#

score_bal_loss <- function(trt, ps, X, P, lambda, theta){
  # Returns the score contributions for a propensity score model fitted
  #   using 'spBalance'.
  # Input:
  #   fit: fitted 'bal' model
  # Output:
  #   A score matrix with score contributions for ach data point.

  n.data <- length(trt)
  score.mat <- -t((trt - ps) / ps / (1 - ps) * X) +
    matrix(lambda * (P %*% theta) / n.data, nrow = length(theta), ncol = n.data)

  return(score.mat)
}


make_S_lambda <- function(fit, average_scale = TRUE) {
  # Creates the penalty matrix for a fitted gam model.
  # Input:
  #   fit: fitted 'gam' model
  #   average_scale: Boolean indicating if empirical average is used
  #     when linearizing ATE estimator; defaults to TRUE
  # Output:
  #   The penalty matrix used by 'mgcv' when fitting the propensity score model.

  p <- length(coef(fit))
  S_lam <- matrix(0, p, p)

  # Use full.sp when available; it is the full set multiplying penalties
  sp <- if (!is.null(fit$full.sp)) fit$full.sp else fit$sp

  # If min.sp was used, mgcv documentation says it must be added to full.sp
  # to get the smoothing parameters actually multiplying penalties.
  if (!is.null(fit$min.sp) && length(fit$min.sp) == length(sp)) {
    sp <- sp + fit$min.sp
  }

  for (sm in fit$smooth) {
    if (length(sm$S) == 0L) next

    coef_ind <- sm$first.para:sm$last.para

    sp_ind <- if (!is.null(sm$first.sp) && !is.null(sm$last.sp)) {
      sm$first.sp:sm$last.sp
    } else {
      seq_along(sm$S)
    }

    for (k in seq_along(sm$S)) {
      S_lam[coef_ind, coef_ind] <-
        S_lam[coef_ind, coef_ind] + sp[sp_ind[k]] * sm$S[[k]]
    }
  }

  # If you supplied a fixed quadratic penalty H to gam(), add it if stored.
  if (!is.null(fit$H)) {
    S_lam <- S_lam + fit$H
  }

  # mgcv penalties are on the summed likelihood/deviance scale.
  # For an estimating equation written as P_n s_i = 0, divide by n.
  if (average_scale) {
    S_lam <- S_lam / nobs(fit)
  }

  S_lam
}


score_gam_likelihood <- function(fit, X) {
  # Returns the score contributions for a propensity score model fitted
  #   using 'gam'.
  # Input:
  #   fit: fitted 'gam' model
  # Output:
  #   A list containing (i) score -- a matrix with score contributions for
  #     each data point, and (ii) B, the estimated Hessian matrix.

  beta  <- coef(fit)
  A     <- fit$y
  e     <- fitted(fit)
  n     <- nrow(X)

  Sbar <- make_S_lambda(fit, average_scale = TRUE)
  pen_grad <- drop(Sbar %*% beta)

  # per-observation score, rows are s_i^T
  S_i <- sweep(X * as.numeric(e - A), 2, pen_grad, "+")

  # Hessian/Jacobian B
  W <- as.numeric(e * (1 - e))
  B <- crossprod(X, X * W) / n + Sbar

  list(score = S_i, B = B, Sbar = Sbar)
}


score_glm_likelihood <- function(fit){
  # Returns the score contributions for a propensity score model fitted
  #   using 'glm'.
  # Input:
  #   fit: fitted 'glm' model
  # Output:
  #   A list containing (i) score -- a matrix with score contributions for
  #     each data point, and (ii) B, the estimated Hessian matrix.

  # extract components from the model
  X <- model.matrix(fit)
  e.hat <- as.numeric(fitted(fit))
  A     <- fit$y
  n     <- nrow(X)

  # per-observation score, rows are s_i^T
  S_i <- X * as.numeric(e.hat - A)

  # Hessian/Jacobian B
  W <- as.numeric(e.hat * (1 - e.hat))
  B <- crossprod(X, X * W) / n

  list(score = S_i, B = B)
}

#' Influence functions of IPTW and augmented IPTW estimators
#'
#' ATE estimates and their estimated influence function values for the
#' Horvitz-Thompson (`HTInfluence`), Hajek (`HajekInfluence`) and augmented
#' IPTW (`aIPTWInfluence`) estimators. The influence values account for
#' estimation of the propensity score and are the input to [HAC.se()].
#'
#' @param out outcome vector.
#' @param trt binary treatment vector.
#' @param fit,ps.fit fitted propensity score model: a [spBalance()] fit
#'   (class `"bal"`), or a `glm` or `gam` fit.
#' @param safe.est logical; trim extreme influence values (HT only).
#' @param ate.only logical; if `TRUE`, only the ATE is computed (`phi` is `NULL`).
#' @param trt.name name of the treatment variable in `mu.fit`.
#' @param mu.fit fitted outcome regression model (e.g. from [mgcv::gam()]),
#'   including the treatment as a covariate.
#' @return A list with `tau`, the ATE estimate, and `phi`, the influence
#'   function values (one per observation).
#' @seealso [spbalATE()], [HAC.se()]
#' @export
HTInfluence <- function(out, trt, fit, safe.est = FALSE, ate.only = FALSE){
  # Calculates the influence function values for the Horvitz-Thompson estimator.
  # Input:
  #   out: outcome vector
  #   trt: treatment vector
  #   fit: fitted propensity score model (either 'glm', 'gam', or 'bal')
  # Output:
  #   ATE estimate and associated influence function field.

  type <- class(fit)[1]

  if (type == "bal"){

    # fit structures
    n.samps <- length(trt)
    ps = as.numeric(fit$pi.hat)
    X = fit$X
    P = fit$P
    lambda = fit$lambda
    theta = fit$par


    if (!ate.only){
      # balancing score
      score.mat <- score_bal_loss(
        trt = trt, ps = ps, X = X, P = P,
        lambda = lambda, theta = theta
      )

      # Hessian
      # note: need to divide by n.samps, as the loss uses
      #       l* = lambda / n.samps, _not_ lambda

      # W <- as.numeric(trt * (1 - ps) / ps + (1 - trt) * ps / (1 - ps))
      # H.mat <- crossprod(X, X * W) / n.samps + lambda * P / n.samps

      H.mat <- crossprod(X) / n.samps + lambda * P / n.samps

    }

  } else if (type == "gam"){

    X <- predict(fit, type = "lpmatrix")
    ps <- fit$fitted.values

    if (!ate.only){
      score.obj <- score_gam_likelihood(fit, X)
      score.mat <- t(score.obj$score)
      H.mat <- score.obj$B
    }
  } else if (type == "glm"){

    ps <- fit$fitted.values
    score.obj <- score_glm_likelihood(fit)

    if (!ate.only){
      score.mat <- t(score.obj$score)
      H.mat <- score.obj$B
      X <- model.matrix(fit)
    }
  }

  # tau estimate
  tau <- mean((trt / ps - (1 - trt) / (1 - ps)) * out)

  if (!ate.only){
    # HT influence function:
    psi.k <- (trt / ps - (1 - trt) / (1 - ps)) * out - tau

    # correction for \hat{\pi}
    M.hat <- colMeans(-(out * (trt * (1 - ps) / ps + (1 - trt) * ps / (1 - ps))) * X)
    alpha.k <- crossprod(M.hat, solve(H.mat, score.mat))[1,]

    # combined influence function:
    phi.k <- psi.k - alpha.k
  } else {
    phi.k <- NULL
  }

  # return tau estimate and influence values:
  list(tau = tau, phi = phi.k)
}



#' @rdname HTInfluence
#' @export
HajekInfluence <- function(out, trt, fit, ate.only = FALSE){
  # Calculates the influence function values for the Hajek estimator.
  # Input:
  #   out: outcome vector
  #   trt: treatment vector
  #   fit: fitted propensity score model (either 'glm', 'gam', or 'bal')
  # Output:
  #   ATE estimate and associated influence function field.

  type <- class(fit)[1]

    if (type == "bal"){

    # fit structures
    n.samps <- length(trt)
    ps = as.numeric(fit$pi.hat)
    X = fit$X
    P = fit$P
    lambda = fit$lambda
    theta = fit$par


    if (!ate.only){
      # balancing score
      score.mat <- score_bal_loss(
        trt = trt, ps = ps, X = X, P = P,
        lambda = lambda, theta = theta
      )

      # Hessian
      # note: need to divide by n.samps, as the loss uses
      #       l* = lambda / n.samps, _not_ lambda

      # W <- as.numeric(trt * (1 - ps) / ps + (1 - trt) * ps / (1 - ps))
      # H.mat <- crossprod(X, X * W) / n.samps + lambda * P / n.samps

      H.mat <- crossprod(X) / n.samps + lambda * P / n.samps

    }

  } else if (type == "gam"){

    X <- predict(fit, type = "lpmatrix")
    ps <- fit$fitted.values

    if (!ate.only){
      score.obj <- score_gam_likelihood(fit, X)
      score.mat <- t(score.obj$score)
      H.mat <- score.obj$B
    }
  } else if (type == "glm"){

    ps <- fit$fitted.values
    score.obj <- score_glm_likelihood(fit)

    if (!ate.only){
      score.mat <- t(score.obj$score)
      H.mat <- score.obj$B
      X <- model.matrix(fit)
    }
  }

  # Hajek influence function:
  w1.bar <- mean(trt / ps); w0.bar <- mean((1 - trt) / (1 - ps))
  mu1 <- mean(trt / ps * out / w1.bar)
  mu0 <- mean((1 - trt) / (1 - ps) * out / w0.bar)
  tau <- mu1 - mu0

  if (!ate.only){
    # influence function with known ps
    psi.k <- (trt / ps / w1.bar - (1 - trt) / (1 - ps) / w0.bar) * out - tau

    # correction for \hat{\pi}
    mu1.linearization <- -trt * (1 - ps) / ps / w1.bar * (out - mu1) * X
    mu0.linearization <- -(1 - trt) * ps / (1 - ps) / w0.bar * (out - mu0) * X
    M.hat <- colMeans(mu1.linearization + mu0.linearization)
    alpha.k <- crossprod(M.hat, solve(H.mat, score.mat))[1,]

    # combined influence function:
    phi.k <- psi.k - alpha.k
  } else {
    phi.k <- NULL
  }

  # return tau estimate and influence values:
  list(tau = tau, phi = phi.k)
}


#' @rdname HTInfluence
#' @export
aIPTWInfluence <- function(out, trt, trt.name, ps.fit, mu.fit, ate.only = FALSE){
  # Calculates the influence function values for the aIPTW estimator.
  # Input:
  #   out: outcome vector
  #   trt: treatment vector
  #   trt.name: name of treatment in mu.fit model
  #   ps.fit: fitted propensity score model (either 'glm', 'gam', or 'bal')
  #   mu.fit: fitted outcome regression model ('gam' model)
  # Output:
  #   ATE estimate and associated influence function field.

  # estimated propensity scores
  type <- class(ps.fit)[1]
  if (type == "bal"){
    ps <- as.numeric(ps.fit$pi.hat)
  } else {
    ps <- fitted(ps.fit)
  }

  # estimate y(1)
  data.obj <- mu.fit$model
  trt.data <- data.obj
  trt.data[ ,trt.name] <- 1
  m1 <- predict(mu.fit, newdata = trt.data)

  # estimate y(0)
  ctrl.data <- data.obj
  ctrl.data[ , trt.name] <- 0
  m0 <- predict(mu.fit, newdata = ctrl.data)

  # (i) tau estimate
  tau <- mean((trt * out - (trt - ps) * m1) / ps) -
    mean(((1 - trt) * out + (trt - ps) * m0) / (1 - ps))

  if (!ate.only){
    # augmented IPTW plug-in influence function
    phi.k <- m1 - m0 + (trt * (out - m1) / ps) -
      ((1 - trt) * (out - m0) / (1 - ps)) - tau
  } else {
    phi.k <- NULL
  }

  # return tau estimate and influence values:
  list(tau = tau, phi = as.numeric(phi.k))
}


fit_var_mod <- function(v_emp, v_mod){
  # Capture warning from fit.variogram and indicate if convergence occured.
  # Input
  #   v_emp: empirical variogram
  #   v_mod: variogram model used by gstat::fit.variogram
  # Output
  #   A list with (i) fitted variogram model, (ii) an indicator of convergence,
  #     and (iii) any warning messages.

   warn_msg <- NULL
  v_fit <- withCallingHandlers(
    fit.variogram(v_emp, v_mod),
    warning = function(w) {
      warn_msg <<- w$message
      invokeRestart("muffleWarning")
    }
  )

  # Check if the specific timeout text appeared
  converged <- TRUE
  if (!is.null(warn_msg) && grepl("iterations", warn_msg)) {
    converged <- FALSE
  }

  return(list(model = v_fit, converged = converged, warning = warn_msg))
}


#' Spatial HAC standard error from influence function values
#'
#' Heteroskedasticity- and autocorrelation-consistent (Conley-type) standard
#' error of a mean of spatially correlated influence function values, with the
#' kernel bandwidth chosen from a fitted variogram.
#'
#' @param phi influence function values.
#' @param df data frame with the coordinates of the observations.
#' @param coords names of the coordinate columns in `df`.
#' @param range.mult multiplier used to set the bandwidth from the fitted
#'   variogram range (the default gives an effective correlation of about 0.01
#'   under an exponential model).
#' @param nbins number of bins for the empirical variogram.
#' @param k.type HAC kernel, `"bartlett"` or `"uniform"`.
#' @param safe.est logical; trim extreme influence values?
#' @param corr.tol correlogram threshold used for bandwidth selection.
#' @return The estimated standard error.
#' @export
HAC.se <- function(
    phi, df, coords = c("x", "y"),
    range.mult = -log(0.01),
    nbins = 50, k.type = "bartlett", safe.est = FALSE, corr.tol = 0.01
){
  # Calculates the HAC standard error for a given influence function field.
  # Input:
  #   phi: vector of influence function values
  #   df: data frame that contains coordinates associated with phi
  #   coords: coordinate names, defaults to c("x", "y")
  #   range.mult: multiplicative factor used to determine HAC kernel bandwidth;
  #     defaults to -log(0.01), such that the effective correlation in an
  #     exponential covariance model is ~0.01.
  #   nbins: number of bins used to estimate empirical variogram; defaults to 50
  #   k.type: type of kernel, must be one of "bartlett" (default) or "uniform"
  #   safe.est: whether to trim extreme influence function values; default is FALSE
  #   corr.tol: correlogram tolerance threshold if using Lehner's bandwidth selection
  # Output
  #   Estimated HAC standard error, with bandwidth chosen using empirical
  #     variogram.

  # center influence values
  phi.c <- phi - mean(phi)

  # trim extreme values:
  if (safe.est){
    bounds <- quantile(phi.c, probs = c(0.01, 0.99))
    phi.c[phi.c > bounds[2]] <- bounds[2]
    phi.c[phi.c < bounds[1]] <- bounds[1]
  }

  # data frame for hac calc
  hac.df <- df
  hac.df$phi.c <- phi.c
  coord.loc <- which(colnames(hac.df) %in% coords)
  colnames(hac.df)[coord.loc] <- c("x", "y")

  # distance matrix
  D <- fields::rdist(df[,coords])

  # decide on bins for variogram estimation
  max.d <- max(D)
  maxdist_m <- max.d * 2 / 3
  width_m <- maxdist_m / nbins

  # empirical variogram
  vg <- gstat::variogram(
    phi.c ~ 1,
    data = hac.df,
    locations = ~x+y,
    covariogram = FALSE,
    cutoff = maxdist_m,
    width = width_m,
    cressie = TRUE
  )

  # fit spherical variogram model
  var.mod.sph <- fit_var_mod(vg, vgm("Sph", maxdist_m / 2, vg$gamma[1]))
  sph.fail <- !var.mod.sph$converged

  if (!sph.fail){
    # check if range is negative...
    if (var.mod.sph$model$range[2] < 0) {
      sph.fail = TRUE
    } else {
      dist.band <- ceiling(var.mod.sph$model$range[2])
    }
  }

  if (sph.fail){
    # refit with exponential covariogram model
    var.mod.exp <- fit_var_mod(vg, vgm("Exp", maxdist_m / 2, vg$gamma[1]))
    exp.fail <- !var.mod.exp$converged

    if (!exp.fail){
      # check if range is negative...
      if (var.mod.exp$model$range[2] < 0){
        exp.fail = TRUE
      } else {
        # determine distance band for HAC calculation
        dist.band <- ceiling(var.mod.exp$model$range[2]) * range.mult
        while (dist.band >= max.d){
          warning("Estimated distance band is larger than observed max distance.")
          dist.band <- dist.band * 3/4
        }
      }
    }
  } else {
    exp.fail = FALSE
  }

  if (sph.fail & exp.fail){
    warning("Both spherical and exponential variogram models failed to converge.")

    # fit using covariogram tolerance

    covg <- gstat::variogram(
      phi.c ~ 1,
      data = hac.df,
      locations = ~x+y,
      cutoff = maxdist_m,
      width = width_m,
      cressie = TRUE,
      covariogram = TRUE
    )

    # convert to covariogram
    c0 <- var(hac.df$phi.c, na.rm = TRUE)
    covg$rho <- covg$gamma / c0

    # re-order by distance
    cv <- covg[order(covg$dist), ]

    # determine minimum distance with correlelogram under threshold
    dist.band <- with(cv, dist[which(abs(rho) < corr.tol)[1]])
  }

  # estimate HAC standard error
  hac.df$const <- 1
  V.hat <- conleyreg(
    phi.c ~ const,
    data = hac.df,
    model = "ols",
    dist_mat = D,
    dist_cutoff = dist.band,
    kernel = k.type,
    intercept = FALSE,
    vcov = TRUE,
    verbose = FALSE,
    ncores = 1
  )


  if (V.hat < 0){
    V.hat <- NA
    warning("A negative HAC variance estimate was obtained.
         Consider modifying the distance band or kernel type.")
  }

  # return standard error
  sqrt(V.hat[1,1])
}


#' ATE with HAC and iid standard errors
#'
#' Lower-level interface used by [spbalATE()]: one estimator's ATE, with its
#' spatial HAC standard error and (optionally) the iid standard error.
#'
#' @param outcome outcome vector.
#' @param trt binary treatment vector.
#' @param ps.fit fitted propensity score model ([spBalance()], `glm` or `gam`).
#' @param coord.df data frame of coordinates, one row per observation.
#' @param est estimator: `"HT"`, `"Hajek"` or `"aIPTW"`.
#' @param mu.fit fitted outcome model (required for `"aIPTW"`).
#' @param trt.name name of the treatment variable in `mu.fit`.
#' @param iid logical; also return the iid standard error?
#' @param ... further arguments passed to [HAC.se()].
#' @return A list with `tau`, `se.hac` and (if `iid = TRUE`) `se.iid`.
#' @export
est.HAC <- function(
  outcome, trt, ps.fit, coord.df, est,
  mu.fit = NULL, trt.name = "treat", iid = TRUE, ...
){
  # Estimates ATE and corresponding standard error.
  # Input:
  #   outcome: vector with outcome values
  #   trt: treatment vector
  #   ps.fit: fitted propensity score model (either glm, gam, or bal object)
  #   coord.df: data frame with coordinate values
  #   est: type of estimator, must be one of 'HT', 'Hajek', or 'aIPTW'
  #   mu.fit: fitted outcome regression (default = NULL)
  #   trt.name: name of treatment variable (default = "treat")
  #   iid: whether to also compute the iid s.e. estimate (default = TRUE)
  # Output:
  #   A list with tau estimate, se.hac, and se.iid (if iid = TRUE).

  # 1. calculate influence function

  if (est == "HT"){
    inf.func <- HTInfluence(out = outcome, trt = trt, fit = ps.fit)
  } else if (est == "Hajek"){
    inf.func <- HajekInfluence(out = outcome, trt = trt, fit = ps.fit)
  } else if (est == "aIPTW"){
    inf.func <- aIPTWInfluence(
      out = outcome, trt = trt,
      trt.name = trt.name, ps.fit = ps.fit, mu.fit = mu.fit
    )
  }

  if.df <- data.frame(coord.df)
  coord.names <- colnames(coord.df)
  se.hat <- HAC.se(inf.func$phi, if.df, coord.names, ...)

  # return tau estimate and standard error
  res <- list(tau = inf.func$tau, se.hac = se.hat)

  if (iid){
    res$se.iid <- sqrt(mean(inf.func$phi^2) / length(inf.func$phi))
  }
  return(res)
}


# # ----------------------------------------------#
# # --- DEMO                                   ---#
# # ----------------------------------------------#
#
# # sim. details
# sim.start = 0
# nonlin = FALSE
# data.args = list(grid = 100, cov = "SE", r.type = 3, tau = 0.1)
# car.fit = FALSE
# sparse = 1000
# k = 1
#
# ### 1. simulate data
#
# # set seed
# set.seed(sim.start + k)
#
# # simulate data
# sim.data <- dataGen1(sim.args = data.args, K.mat = "default", nonlin)
#
# # convert to SpatRast object
# sim.rast <- dataGen1.SpatRaster(sim.data)
#
# # create data frame
# sim.df <- data.frame(
#     treat = values(sim.rast$z),
#     outcome = values(sim.rast$y),
#     pi = values(sim.rast$pi),
#     v1 = values(sim.rast$x1),
#     v2 = values(sim.rast$x2),
#     crds(sim.rast$y)
# )
# colnames(sim.df) <- c("treat", "outcome", "pi", "v1", "v2", "x", "y")
#
# if (sparse > 0){
#   keep.k <- sample(nrow(sim.df), sparse)
#   sim.df <- sim.df[keep.k, ]
# }
#
#
# if (sparse > 0){
#   n.knots = 300
# } else {
#   n.knots = 500
# }
#
#
# ### 2. estimate using GLM
#
# # fit glm
# glm.fit <- glm(form = treat ~ v1 + v2 + x + y, family = binomial, data = sim.df)
#
# # HT estimate
# est.HAC(
#   outcome = sim.df$outcome,
#   trt = sim.df$treat,
#   ps.fit = glm.fit,
#   coord.df = sim.df[,c(6,7)],
#   est = "HT"
# )
#
# # Hajek estimate
# est.HAC(
#   outcome = sim.df$outcome,
#   trt = sim.df$treat,
#   ps.fit = glm.fit,
#   coord.df = sim.df[,c(6,7)],
#   est = "Hajek"
# )
#
# ### 3. estimate using GAM
#
# # estimate propensity score via mgcv::gam
# gam.fit <- mgcv::gam(form = treat ~ s(x, y, k = n.knots), family = binomial, data = sim.df)
#
# # HT estimate
# est.HAC(
#   out = sim.df$outcome,
#   trt = sim.df$treat,
#   ps.fit = gam.fit,
#   coord.df = sim.df[,c(6,7)],
#   est = "HT"
# )
#
# # Hajek estimate
# est.HAC(
#   out = sim.df$outcome,
#   trt = sim.df$treat,
#   ps.fit = gam.fit,
#   coord.df = sim.df[,c(6,7)],
#   est = "Hajek"
# )
#
#
# ### 6. estimate using spatial balance
#
# # coefvar method
# spb.cvr <- spBalance(
#   formula = treat ~ s(x, y, k = n.knots),
#   data = sim.df,
#   lambda = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 5, 10),
#   tuning = "coefvar",
#   hide.details = TRUE,
#   opt.params = list(tol = 1e-7, max.iter = 1000, alpha = 0.5, beta = 0.5)
# )
#
# # HT estimate
# est.HAC(
#   out = sim.df$outcome,
#   trt = sim.df$treat,
#   ps.fit = spb.cvr,
#   coord.df = sim.df[,c(6,7)],
#   est = "HT"
# )
#
# # Hajek estimate
# est.HAC(
#   out = sim.df$outcome,
#   trt = sim.df$treat,
#   ps.fit = spb.cvr,
#   coord.df = sim.df[,c(6,7)],
#   est = "Hajek"
# )
#
# ### Augmented IPTW
#
# # outcome regression
# out.fit <- gam(outcome ~ treat + s(x, y, k = n.knots), data = sim.df, family = gaussian())
#
# # aIPTW estimate
# est.HAC(
#   out = sim.df$outcome,
#   trt = sim.df$treat,
#   ps.fit = spb.cvr,
#   mu.fit = out.fit,
#   trt.name = "treat",
#   coord.df = sim.df[,c(6,7)],
#   est = "aIPTW"
# )
#
#



