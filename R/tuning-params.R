### tuning-params.R
### Nathan Wikle

### Functions used to select the tuning parameter, lambda.

cvLambda <- function(kfold, formula, data, opt.params, lambda.vals, type = "cv.score", norm = "L2", ...){
  # Choose lambda via cross validation.
  # Input
  #   kfold: number of cv folds to use
  #   formula: propensity score formula
  #   data: data frame containing all variables used in the analysis
  #   opt.params: arguments passed to 'optim' function call
  #   lambda.vals: a vector of tuning parameter values to consider
  #   type: cv metric to use; options include
  #     cv.score: minimize average balancing score on the holdout sets
  #     cv.grad: minimize the norm of the balancing score gradient on the holdout sets
  #     all: returns estimated lambda using both metrics
  #   norm: type of norm used for cv.grad; options = "L1", "L2", "Linf"
  # Output
  #   A list containing the lambda value resulting in smallest CV error and
  # the associated CV metric (eg, score or gradient norm) when using that lambda.

  n.lambda <- length(lambda.vals)
  cv.results <- matrix(NA_real_, nrow = kfold, ncol = n.lambda)

  folds <- caret::createFolds(data[,1], k = kfold)
  cv.score <- cv.grad <- matrix(NA_real_, nrow = kfold, ncol = n.lambda)

  for (k in 1:kfold){

    # fit model
    fit.k <- spBalFit(
      formula = formula,
      data = data[-folds[[k]],],
      lambda = lambda.vals,
      opt.params = opt.params,
      ...
    )

    if ((type == "cv.score") | (type == "all")){
      cv.score[k,] <- cvScore(fit = fit.k, newdata = data[folds[[k]],])
    }

    if ((type == "cv.grad") | (type == "all")){
      cv.grad[k,] <- cvGrad(fit = fit.k, newdata = data[folds[[k]],], norm)
    }

  }


  if ((type == "cv.score") | (type == "all")){
    cv.score.bar <- colMeans(cv.score)
    lambda.score <- lambda.vals[which(cv.score.bar == min(cv.score.bar))]
  } else {
    lambda.score <- NULL
  }

  if ((type == "cv.grad") | (type == "all")){
    cv.grad.bar <- colMeans(cv.grad)
    lambda.grad <- lambda.vals[which(cv.grad.bar == min(cv.grad.bar))]
  } else {
    lambda.grad <- NULL
  }

  ### compile results
  if (type == "cv.score"){
    keep <- c(1)
  } else if (type == "cv.grad"){
    keep <- c(2)
  } else if (type == "all"){
    keep <- c(1,2)
  }

  res <- list(
    lambda = list(lambda.score, lambda.grad)[keep],
    cv = list(cv.score, cv.grad)[keep]
  )

  return(res)
}


predictSpBal <- function(fit, newdata){
  # Predicts propensity score values for new data.
  # Input
  #   fit: a fitted spbalance object from 'spBalFit'
  #   newdata: data frame used as input for prediction
  # Output
  #   A list with the predicted data's design matrix and the
  # estimated propensity scores.

  # initialize new X matrices
  Xs.pred <- NULL
  Xk.pred <- NULL

  # prediction with GAM components
  if (!is.null(fit$gam)){
    if (is.element("gam", class(fit$gam))){
      Xs.pred <- predict(fit$gam, newdata, type = "lpmatrix")
    } else {
      Xs.pred <- model.matrix(formula(fit$gam$pterms), newdata)
      n.sm <- length(which(unlist(fit$terms) == "smooth"))
      if (n.sm > 0){
        for (k in 1:n.sm){
          Xsm.k <- PredictMat(fit$gam$smooth[[k]], newdata)
          Xs.pred <- cbind(Xs.pred, Xsm.k)
        }
      }
      rownames(Xs.pred) <- NULL
    }
  } else {
    if (fit$form$kcentered){
      # add intercept term
      Xs.pred <- matrix(1, nrow = nrow(newdata), ncol = 1)
      colnames(Xs.pred) <- "intercept"
    }
  }

  # prediction with KRR or RFF components
  if (!is.null(fit$krr)){
    n.krr <- length(which(unlist(fit$terms) == "krr"))
    n.rff <- length(which(unlist(fit$terms) == "rff"))

    for (k in 1:(n.krr + n.rff)){
      # grab krr.k object
      krr.k <- fit$krr[[k]]
      # return prediction basis matrix
      K.pred <- krrPredMat(krr.k, newdata)
      Xk.pred <- cbind(Xk.pred, K.pred)
    }
  }

  # combine into single X matrix
  X.pred <- cbind(Xs.pred, Xk.pred)
  pi.pred <- X.pred %*% fit$par

  # return fitted model
  res <- list(
    X = X.pred,
    pi.hat = pi.pred
  )
  return(res)
}


krrPredMat <- function(fit, newdata){
  # Creates a design matrix used when predicting with a KRR model term.
  # Input
  #   fit: a fitted spbalance object from 'spBalFit'
  #   newdata: data frame used as input for prediction
  # Output
  #   The KRR design matrix for new locations.

  # restrict newdata to appropriate terms
  if (is.null(fit$krr.obj$term)){
    X.unscaled <- newdata[,fit$krr.obj$x.inds]
  } else {
    X.unscaled <- newdata[,fit$krr.obj$term]
  }

  # re-scale new data
  if (fit$krr.obj$scale){
    mu <- attr(fit$X, "scaled:center")
    sd <- attr(fit$X, "scaled:scale")
    X.sc <- t((t(X.unscaled) - mu) / sd)
  } else {
    X.sc <- as.matrix(X.unscaled)
  }

  # check if object is 'krr' or 'rff'
  if (fit$krr.obj$rff){
    # rff -> create features matrix

    # create features
    n.features <- length(fit$b)
    f.pred <- X.sc %*% t(fit$w) + rep(fit$b, each = nrow(X.sc))
    Z.pred <- sqrt(2) * cos(f.pred) / sqrt(n.features)

    if (fit$krr.obj$cent){
      mu.z <- attr(fit$Z, "scaled:center")
      # sd.z <- attr(fit$Z, "scaled:scale")
      sigma.p <- t((t(Z.pred) - mu.z))
    } else {
      sigma.p <- Z.pred
    }

  } else {
    # krr -> create covariance matrix
    d.mat <- as.matrix(dist(rbind(fit$X, X.sc), diag = TRUE, upper = TRUE))
    n.fit <- nrow(fit$X)
    n.pred <- nrow(X.sc)

    if (fit$krr.obj$cent){

      s.full <- kernelEval(
        d.mat,
        tuning = fit$krr.obj$k.params,
        type = fit$krr.obj$kernel,
        centered = FALSE
      )

      s.fit <- s.full[1:n.fit, 1:n.fit]
      s.pred <- s.full[n.fit + 1:n.pred, 1:n.fit]

      # centered prediction matrix

      # one.mf <- matrix(1 / n.fit, nrow = n.fit, ncol = n.fit)
      # one.mp <- matrix(1 / n.fit, nrow = n.pred, ncol = n.fit)
      # sigma <- s.pred - one.mp %*% s.fit - s.pred %*% one.mf + one.mp %*% (s.fit %*% one.mf)
      sigma.p <- s.pred -
        matrix(colMeans(s.fit), nrow = n.pred, ncol = n.fit, byrow = TRUE) -
        matrix(rowMeans(s.pred), nrow = n.pred, ncol = n.fit, byrow = FALSE) +
        matrix(mean(rowMeans(s.fit)), nrow = n.pred, ncol = n.fit)

    } else {
      sigma.p <- kernelEval(
        d.mat[n.fit + 1:n.pred, 1:n.fit],
        tuning = fit$krr.obj$k.params,
        type = fit$krr.obj$kernel,
        centered = FALSE
      )
    }
  }

  return(sigma.p)
}


cvScore <- function(fit, newdata){
  # Evalute the balancing score of a fitted model on a CV holdout set.
  # Input
  #   fit: a fitted spbalance object from 'spBalFit'
  #   newdata: data frame used as input for prediction
  # Output
  #   The balancing score, evaluated on the holdout set.

  # get predictions on newdata
  newfit <- predictSpBal(fit, newdata)

  # treatments from newdata
  y.pred <- newdata[,all.vars(formula(fit$form$smf))[1]]

  # calculate balancing score on holdout set
  if (is.matrix(fit$par)){
    cv.score <- apply(fit$par, 2, FUN = function(theta){
      balanceObjFun(y = y.pred, X = newfit$X, theta = theta, lambda = 0)
    })
  } else {
    cv.score <- balanceObjFun(y = y.pred, X = newfit$X, theta = fit$par, lambda = 0)
  }

  return(cv.score)
}


cvGrad <- function(fit, newdata, norm = "L2"){
  # Evalute the norm of the balancing score's gradient on a CV holdout set.
  # Input
  #   fit: a fitted spbalance object from 'spBalFit'
  #   newdata: data frame used as input for prediction
  #   norm: type of norm to use; one of "L1", "L2", "Linf"
  # Output
  #   The norm of the balancing score's gradient, evaluated on the holdout set.

  # get predictions on newdata
  newfit <- predictSpBal(fit, newdata)
  # treatments from newdata
  y.pred <- newdata[,all.vars(formula(fit$form$smf))[1]]
  # pi prediction
  pi.pred <- newfit$pi.hat

  # gradient
  new.z <- (y.pred - pi.pred) / pi.pred / (1 - pi.pred)
  grad.calc <- -crossprod(newfit$X, new.z) / nrow(newfit$X)

  if (norm == "L1"){
    grad.norm <- colSums(abs(grad.calc))
  } else if (norm == "L2"){
    grad.norm <- sqrt(apply(grad.calc, 2, crossprod))
  } else if (norm == "Linf"){
    grad.norm <- apply(abs(grad.calc), 2, max)
  }

  return(grad.norm)
}

coefVar <- function(z, pi.vals, lambda, ratio = 0.9){
  # Choose lambda using the balancing weights coefficient of variation.
  # Input
  #   z: treatment vector
  #   pi.vals: matrix of propensity score (correpsonding to tuning params)
  #   lambda: vector of tuning parameters
  #   ratio: coefficient of variation ratio --- find max lambda such that
  #     coefvar is >= ratio * max(coefvar); default is 0.9
  # Output
  #   A list with the chosen tuning parameter and the associated
  # coefficient of variation ratio.

  n.lambda <- length(lambda)

  if (n.lambda > 1){
    # normalizing weights
    w1.norm <- colSums(z / pi.vals)
    w2.norm <- colSums((1 - z) / (1 - pi.vals))

    # final weights calc
    weights.m <- matrix(NA_real_, nrow = length(z), ncol = n.lambda)
    for (k in 1:n.lambda){
      weights.m[,k] <- (z / pi.vals[,k] / w1.norm[k]) +
        ((1-z) / (1 - pi.vals[,k]) / w2.norm[k])
    }

    # coefficient of variation
    coefvar <- apply(weights.m, 2, FUN = function(x){sd(x)/mean(x)})

  } else {
    # normalizing weights
    w1.norm <- sum(z / pi.vals)
    w2.norm <- sum((1 - z) / (1 - pi.vals))
    # final weights calc
    weights <- (z / pi.vals / w1.norm) + ((1-z) / (1 - pi.vals) / w2.norm)
    # coefficient of variation
    coefvar <- sd(weights) / mean(weights)
  }

  # choice of lambda
  max.cv <- max(coefvar)
  lambda.opt <- max(lambda[coefvar >= ratio * max.cv])

  return(
    list(
      lambda = lambda.opt,
      cv.ratio = coefvar / max.cv
    )
  )
}


maxBal <- function(fit, z, data, max.diff = 0.1, return.fit = FALSE){
  # Find the tuning largest tuning parameter such that the standardized
  #   difference in means for all model terms is less than some threshold.
  # Input
  #   fit: a fitted spbalance object from 'spBalFit'
  #   z: treatment assignment vector
  #   data: data frame containing all variables used in the analysis
  #   max.diff: st. diff-in-means threshold; default is 0.1
  #   return.fit: Boolean indicating if parameters and propensity scores
  #     should be returned, or just best lambda and max.sdm values.
  # Output
  #   A list containing best lambda, the associated max st. difference-in-means,
  # an (optionally) the parameter vec and prop. scores corresponding to lambda.

  # covariates
  model.pred <- predictSpBal(fit, data)
  X.m <- model.pred$X

  # estimated spatial smooth
  fixed.terms <- which(unlist(fit$terms) == "fixed")
  fixed.i <- length(fixed.terms)
  if (fixed.i > 0){
    n.fixed <- length(fit$dims[[fixed.terms]])
  } else {
    n.fixed <- 0
  }
  sm.terms <- which(unlist(fit$terms) == "smooth")
  n.sm <- length(sm.terms)
  krr.terms <- which(unlist(fit$terms) == "krr")
  n.krr <- length(krr.terms)
  rff.terms <- which(unlist(fit$terms) == "rff")
  n.rff <- length(rff.terms)

  # create matrix with "balancing" functions
  n.lambda <- ncol(fit$par)

  if (fixed.i > 0){
    n.bal <- n.fixed +
      n.sm + n.krr + n.rff
  } else {
    n.bal <- n.sm + n.krr + n.rff
  }


  bal.m <- array(data = NA_real_, dim = c(nrow(X.m), n.bal, n.lambda))

  if (fixed.i > 0){
    bal.m[, 1:n.fixed, ] <- X.m[,fit$dims[[fixed.terms]]]
  }

  if (n.sm > 0){
    for (k in 1:n.sm){
      index.k <- fit$dims[[sm.terms[k]]]
      bal.m[, n.fixed + k, ] <- X.m[, index.k] %*% fit$par[index.k,]
    }
  }

  if (n.krr > 0){
    for (k in 1:n.krr){
      index.k <- fit$dims[[krr.terms[k]]]
      bal.m[, n.fixed + n.sm + k, ] <- X.m[, index.k] %*% fit$par[index.k,]
    }
  }

  if (n.rff > 0){
    for (k in 1:n.rff){
      index.k <- fit$dims[[rff.terms[k]]]
      bal.m[, n.fixed + n.sm + n.krr + k, ] <- X.m[, index.k] %*% fit$par[index.k,]
    }
  }

  trt <- which(z==1)
  ctrl <- which(z==0)
  diff.vals <- rep(0, n.lambda)

  for (k in 1:n.lambda){
    # normalizing weights
    w1.norm <- sum(z / fit$pi.hat[,k])
    w2.norm <- sum((1 - z) / (1 - fit$pi.hat[,k]))
    # final weights calc
    weights <- (z / fit$pi.hat[,k] / w1.norm) + ((1-z) / (1 - fit$pi.hat[,k]) / w2.norm)
    diff.vals[k] <- max(calcBalance(bal.m[,,k], trt1 = trt, trt2 = ctrl, weights, d.type = 2)[,2])
  }

  # optimal choice of lambda
  best.lambda <- which(diff.vals < max.diff)
  if (length(best.lambda) > 0){
    lambda.i <- which(fit$lambda == max(fit$lambda[best.lambda]))
  } else {
    lambda.i <- which(diff.vals == min(diff.vals))
  }

  if (return.fit){
    # return best fitting pi.hat and theta.hat
    list(
      theta.hat = fit$par[,lambda.i],
      pi.hat = fit$pi.hat[,lambda.i],
      lambda = fit$lambda[lambda.i],
      max.smd = diff.vals
    )
  } else {
    # return optimal lambda and max standardized difference
    list(
      lambda = fit$lambda[lambda.i],
      max.smd = diff.vals
    )
  }
}


calcBalance <- function(x.mat, trt1, trt2, weights, d.type = 2){
  # Calculate the absolute standardized difference in means between trt/ctrl.
  # Input:
  #   x.mat: covariate matrix
  #   trt1: vector indicating indices of treated units
  #   trt2: vector indicating indices of untreated units
  #   weights: propensity score weights
  #   d.type: type of standard deviation used (default = 2: pooled standard deviation)
  # Output:
  #   Matrix with estimated absolute mean differences for each covariate.

  if (is.matrix(x.mat)){
    diff <- balanceEst(
      m1 = x.mat[trt1, ],
      wt1 = weights[trt1],
      m2 = x.mat[trt2, ],
      wt2 = weights[trt2],
      d.type = 2
    )
  } else {
    diff <- balanceEst(
      m1 = x.mat[trt1],
      wt1 = weights[trt1],
      m2 = x.mat[trt2],
      wt2 = weights[trt2],
      d.type = 2
    )
  }

  # return mean differences
  return(diff)
}


balanceEst <- function(m1, m2, wt1, wt2, d.type){
  # Estimate covariate balance (absolute difference-in-means).
  # Input:
  #   m1: covariate matrix for treated units
  #   wt1: propensity score weights for treated units
  #   m2: covariate matrix for untreated units
  #   wt2: propensity score weights for untreated units
  #   d.type: standard deviation used
  #     (1 = sd of treated group, 2 = pooled standard deviation)
  # Output:
  #   Matrix with estimated absolute mean difference between groups

  if (d.type == 1){
    # standard deviation of the treated group
    denom <- apply(m1, 2, sd)
  } else {

    if (is.matrix(m1)){
      # pooled sd (square root of the mean of the group variances)
      trt.var <- apply(m1, 2, var)
      notrt.var <- apply(m2, 2, var)
      denom <- sqrt(rowMeans(cbind(trt.var, notrt.var)))
      n.variables <- ncol(m1)
    } else {
      # pooled sd (square root of the mean of the group variances)
      trt.var <- var(m1)
      notrt.var <- var(m2)
      denom <- sqrt(mean(c(trt.var, notrt.var)))
      n.variables <- 1
    }
  }

  diff.mat <- matrix(NA_real_, nrow = n.variables, ncol = 2)
  colnames(diff.mat) <- c("Diff.Un", "Diff.Adj")
  rownames(diff.mat) <- colnames(m1)

  if (n.variables > 1){
    diff.mat[,1] <- abs(colMeans(m1) - colMeans(m2)) / denom
    m1.wt <- apply(m1, 2, function(x) {wt1 * x / sum(wt1)})
    m2.wt <- apply(m2, 2, function(x) {wt2 * x / sum(wt2)})
    diff.mat[,2] <- abs(colSums(m1.wt) - colSums(m2.wt)) / denom
  } else {
    diff.mat[,1] <- abs(mean(m1) - mean(m2)) / denom
    m1.wt <- wt1 * m1 / sum(wt1)
    m2.wt <- wt2 * m2 / sum(wt2)
    diff.mat[,2] <- abs(sum(m1.wt) - sum(m2.wt)) / denom
  }

  return(diff.mat)
}

