### ate-est.R
### Nathan Wikle

### Functions to estimate the average treatment effect (ATE).

ateUnadjust <- function(out, trt){
  # Estimates the ATE using unadjusted difference in means.
  # Input
  #   out: outcome vector
  #   trt: treatment vector
  # Output
  #   ATE estimate

  # average of treated vs untreated
  estimand <- mean(out[which(trt == 1)]) - mean(out[which(trt == 0)])
  return(estimand)
}

ateGLM <- function(form, data, out, est.type){
  # Estimates the ATE using IPTW with a logistic regression propensity
  #   score model.
  # Input
  #   form: formula
  #   data: data frame with all data used in the analysis
  #   out: outcome vector
  #   est.type: type of IPTW estimand; options are
  #     'HT' (Horvitz-Thompson), 'Hajek', 'OW' (overlap weights), and 'all'
  # Output
  #   ATE estimate.

  # estimate propensity score via glm
  glm.fit <- glm(form, family = binomial, data = data)
  pi.hat <- glm.fit$fitted.values

  # estimate ATE
  trt <- glm.fit$y
  ate.hat <- ateEst(out, trt, pi.hat)

  # determine which estimand to return
  estimand <-  ate.hat$ate
  names(estimand) <- c("HT", "Hajek", "OW")

  if (est.type == "all"){
    return(estimand)
  } else if (type == "HT"){
    return(estimand[1])
  } else if (type == "Hajek"){
    return(estimand[2])
  } else if (type == "OW"){
    return(estimand[3])
  }

}

ateOracle <- function(out, trt, pi.0, est.type){
  # Estimates the ATE using IPTW with the "true" propensity score.
  # Input
  #   out: outcome vector
  #   trt: treatment vector
  #   pi.0: true propensity scores
  #   est.type: type of IPTW estimand; options are
  #     'HT' (Horvitz-Thompson), 'Hajek', 'OW' (overlap weights), and 'all'
  # Output
  #   ATE estimate.

  # estimate ate using true propensity score
  ate.hat <- ateEst(out, trt, pi.0)

  # determine which estimand to return
  estimand <-  ate.hat$ate
  names(estimand) <- c("HT", "Hajek", "OW")

  if (est.type == "all"){
    return(estimand)
  } else if (type == "HT"){
    return(estimand[1])
  } else if (type == "Hajek"){
    return(estimand[2])
  } else if (type == "OW"){
    return(estimand[3])
  }
}


ateBART <- function(form, data, out, trt, est.type, ...){
  # Estimates the ATE using IPTW with the propensity score estimated
  #   using BART.
  # Input
  #   form: formula (passed to bart object)
  #   data: data frame containing all data used in analysis
  #   out: outcome vector
  #   trt: treatment vector
  #   est.type: type of IPTW estimand; options are
  #     'HT' (Horvitz-Thompson), 'Hajek', 'OW' (overlap weights), and 'all'
  # Output
  #   ATE estimate.

  # design matrix
  xmat <- model.matrix(form, data)

  # estimate propensity score with bart
  prop.model <- bart(
    x.train = xmat[,-1], y.train = trt, ...
  )
  pi.hat <- colMeans(pnorm(prop.model$yhat.train))

  # estimate ATE
  ate.hat <- ateEst(out, trt, pi.hat)

  # determine which estimand to return
  estimand <-  ate.hat$ate
  names(estimand) <- c("HT", "Hajek", "OW")

  if (est.type == "all"){
    return(estimand)
  } else if (type == "HT"){
    return(estimand[1])
  } else if (type == "Hajek"){
    return(estimand[2])
  } else if (type == "OW"){
    return(estimand[3])
  }
}

ateGAM <- function(form, data, out, est.type){
  # Estimates the ATE via IPTW with a generalized additive propensity
  #   score model. This allows for the inclusion of thin-plate spline regression
  #   terms, e.g., to account a latent spatial confounder.
  # Input
  #   form: formula (passed to gam)
  #   data: data frame containing all data used in analysis
  #   out: outcome vector
  #   est.type: type of IPTW estimand; options are
  #     'HT' (Horvitz-Thompson), 'Hajek', 'OW' (overlap weights), and 'all'
  # Output
  #   ATE estimate.

  # estimate propensity score via mgcv::gam
  gam.fit <- mgcv::gam(form, family = binomial, data = data)
  # gam.check(gam.fit)
  pi.hat <- gam.fit$fitted.values

  # estimate ATE
  trt <- gam.fit$y
  ate.hat <- ateEst(out, trt, pi.hat)

  # determine which estimand to return
  estimand <-  ate.hat$ate
  names(estimand) <- c("HT", "Hajek", "OW")


  if (est.type == "all"){
    ate <- estimand
  } else if (type == "HT"){
    ate <- estimand[1]
  } else if (type == "Hajek"){
    ate <- estimand[2]
  } else if (type == "OW"){
    ate <- estimand[3]
  }

  results <- list(
    ate = estimand,
    gam.fit = gam.fit,
    pi.hat = pi.hat
  )
}


ateGLMM <- function(form, data, out, est.type, ...){
  # Estimates the ATE via IPTW with a generalized linear mixed model for the
  #   propensity score. This allows for the inclusion of a spatial random
  #   effect in the propensity score model, e.g., to account a latent
  #   spatial confounder.
  # Input
  #   form: formula (passed to spglm)
  #   data: data frame containing all data used in analysis
  #   out: outcome vector
  #   est.type: type of IPTW estimand; options are
  #     'HT' (Horvitz-Thompson), 'Hajek', 'OW' (overlap weights), and 'all'
  # Output
  #   ATE estimate.


  # fit model with spatial random effect
  spglmm.fit < spglm(formula = treat ~ v1 + v2, family = binomial, data = data, ...)
  pi.hat <- spglmm.fit$fitted$response

  # estimate ATE
  trt <- spglmm.fit$y
  ate.hat <- ateEst(out, trt, pi.hat)

  # determine which estimand to return
  estimand <-  ate.hat$ate
  names(estimand) <- c("HT", "Hajek", "OW")

  if (est.type == "all"){
    return(estimand)
  } else if (type == "HT"){
    return(estimand[1])
  } else if (type == "Hajek"){
    return(estimand[2])
  } else if (type == "OW"){
    return(estimand[3])
  }
}

ateBAL <- function(form, data, out, est.type, opt.params, lambda.vals, tuning = "max.bal", folds = NULL){
  # Estimates the ATE via IPTW with a covariate balancing propensity score.
  #   This follows the approach of Zhao (2019); spatial effects can
  #   be represented with a TPRS or KRR term.
  # Input
  #   form: formula (passed to spglm)
  #   data: data frame containing all data used in analysis
  #   out: outcome vector
  #   est.type: type of IPTW estimand; options are
  #     'HT' (Horvitz-Thompson), 'Hajek', 'OW' (overlap weights), and 'all'
  #   opt.params: optimization arguments passed to 'optim'
  #   lambda.vals: penalty tuning parameter values to consider
  #   tuning: method used to select tuning parameter; options include
  #       'cv.score', 'cv.grad', 'coefvar', 'max.bal', 'min', and 'all'
  #   folds: number of CV folds; default is NULL
  # Output
  #   A list with the ATE estimate as well as the selected tuning parameter.

  # estimate balancing propensity score using full data
  full.fit <- spatialBalanceOptim(
    formula = form,
    data = data,
    lambda = lambda.vals,
    opt.params = opt.params
  )
  trt <- full.fit$gam.obj$y

  if (tuning == "all"){
    ate.vals <- list()
    lambda.final <- rep(NA_real_, 6)
  }

  # select lambda via CV
  if ((tuning == "cv.score") | (tuning == "cv.grad") | (tuning == "all")){
    cv.res <- cvLambda(
      kfold = folds,
      formula = form,
      data = data,
      opt.params = opt.params,
      lambda.vals = lambda.vals,
      type = tuning
    )

    if (tuning == "cv.score"){
      l.cv.score <- which(lambda.vals == cv.res$lambda[[1]])
      ate.vals <- ateEst(out, trt, pr.trt = full.fit$pi.hat[,l.cv.score])
      lambda.final <- cv.res$lambda[[1]]
    } else if (tuning == "cv.grad"){
      l.cv.grad <- which(lambda.vals == cv.res$lambda[[2]])
      ate.vals <- ateEst(out, trt, pr.trt = full.fit$pi.hat[,l.cv.grad])
      lambda.final <- cv.res$lambda[[2]]
    } else if (tuning == "all"){
      l.cv.score <- which(lambda.vals == cv.res$lambda[[1]])
      ate.vals[[1]] <- ateEst(out, trt, pr.trt = full.fit$pi.hat[,l.cv.score])
      lambda.final[1] <- cv.res$lambda[[1]]
      l.cv.grad <- which(lambda.vals == cv.res$lambda[[2]])
      ate.vals[[2]] <- ateEst(out, trt, pr.trt = full.fit$pi.hat[,l.cv.grad])
      lambda.final[2] <- cv.res$lambda[[2]]
    }
  }

  # select lambda via coefficent of variation
  if ((tuning == "coefvar") | (tuning == "all")){

    coefv.9 <- coefVar(
      z = trt,
      pi.vals = full.fit$pi.hat,
      lambda = lambda.vals,
      ratio = 0.9
    )

    coefv.5 <- coefVar(
      z = trt,
      pi.vals = full.fit$pi.hat,
      lambda = lambda.vals,
      ratio = 0.5
    )

    if (tuning == "coefvar"){
      l.cv9 <- which(lambda.vals == coefv.9$lambda)
      ate.vals <- ateEst(out, trt, pr.trt = full.fit$pi.hat[,l.cv9])
      lambda.final <- coefv.9$lambda
    } else if (tuning == "all"){
      l.cv9 <- which(lambda.vals == coefv.9$lambda)
      ate.vals[[3]] <- ateEst(out, trt, pr.trt = full.fit$pi.hat[,l.cv9])
      lambda.final[3] <- coefv.9$lambda
      l.cv5 <- which(lambda.vals == coefv.5$lambda)
      ate.vals[[4]] <- ateEst(out, trt, pr.trt = full.fit$pi.hat[,l.cv5])
      lambda.final[4] <- coefv.5$lambda
    }
  }

  # select lambda via setting balance to a specific threshold
  if ((tuning == "max.bal") | (tuning == "all")){

    mbal <- maxBal(fit = full.fit, z = trt, data = data, max.diff = 0.1, return.fit = FALSE)

    if (tuning == "max.bal"){
      l.bal <- which(lambda.vals == mbal$lambda)
      ate.vals <- ateEst(out, trt, pr.trt = full.fit$pi.hat[,l.bal])
      lambda.final <- mbal$lambda
    } else if (tuning == "all"){
      l.bal <- which(lambda.vals == mbal$lambda)
      ate.vals[[5]] <- ateEst(out, trt, pr.trt = full.fit$pi.hat[,l.bal])
      lambda.final[5] <- mbal$lambda
    }
  }

  # use minimum lambda value
  if (tuning == "min"){
    l.min <- which(lambda.vals == min(lambda.vals))
    ate.vals <- ateEst(out, trt, pr.trt = full.fit$pi.hat[,l.min])
    lambda.final <- lambda.vals[l.min]
  } else if (tuning == "all"){
    l.min <- which(lambda.vals == min(lambda.vals))
    ate.vals[[6]] <- ateEst(out, trt, pr.trt = full.fit$pi.hat[,l.min])
    lambda.final[6] <- lambda.vals[l.min]
  }

  # compile results
  if (tuning == "all"){
    # return all estimands
    estimand <- matrix(data = NA_real_, nrow = 3, ncol = 6)
    rownames(estimand) <- c("HT", "Hajek", "OW")
    colnames(estimand) <- c("cv.score", "cv.grad", "coefv9", "coefv5", "maxbal", "min")
    for (k in 1:6){
      estimand[,k] <- ate.vals[[k]][[2]]
    }
    # determine which estimand to return
    if (est.type == "all"){
      final.est <- estimand
    } else if (type == "HT"){
      final.est <- estimand[1,]
    } else if (type == "Hajek"){
      final.est <- estimand[2,]
    } else if (type == "OW"){
      final.est <- estimand[3,]
    }
  } else {
    # return a specific estimand
    estimand <- ate.vals[[2]]
    names(estimand) <- c("HT", "Hajek", "OW")
    # determine which estimand to return
    if (est.type == "all"){
      final.est <- estimand
    } else if (type == "HT"){
      final.est <- estimand[1]
    } else if (type == "Hajek"){
      final.est <- estimand[2]
    } else if (type == "OW"){
      final.est <- estimand[3]
    }
  }

  # return final estimate and choice of lambda
  res <- list(
    estimate = final.est,
    lambda = lambda.final
  )
  return(res)
}

ateEst <- function(out, trt, pr.trt){
  # Estimates the ATE using IPTW for a given propenisty score model.
  # Input
  #   out: outcome vector
  #   trt: treatment vector
  #   pr.trt: propensity score vector
  # Output
  #   Returns the HT (Horvitz-Thompson), Hajek-type, and overlap weight
  # estimates of the ATE.

  # HT-type estimate
  theta.ipw1 <- c(
    sum(((1 - trt) * out) / (1 - pr.trt)) / length(trt),
    sum((trt * out) / pr.trt) / length(trt)
  )

  # Hajek-type estimate

  w1 <- sum(trt / pr.trt)
  w2 <- sum((1 - trt) / (1 - pr.trt))

  theta.ipw2 <- c(
    sum(((1 - trt) * out) / (1 - pr.trt)) / w2,
    sum((trt * out) / pr.trt) / w1
  )

  # overlap weights estimate
  ow1 <- (1 - pr.trt) * trt
  ow2 <- pr.trt * (1 - trt)
  theta.ow <- c(
    sum(((1 - trt) * out) * ow2) / sum(ow2),
    sum((trt * out) * ow1) / sum(ow1)
  )

  theta.est <- rbind(
    theta.ipw1, theta.ipw2, theta.ow
  )
  rownames(theta.est) <- c("HT", "Hajek", "OW")
  colnames(theta.est) <- c("mu0", "mu1")
  ate.est <- apply(theta.est, 1, base::diff)

  # return estimates
  list(mu = theta.est, ate = ate.est)
}


