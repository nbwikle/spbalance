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


ateBAL <- function(form, data, out, est.type, opt.params, lambda.vals, tuning = "max.bal", coefvar.r = 0.9, folds = 10){
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
  #   coefvar.r: ratio used in "coefvar" tuning parameter selection; scalar or vector valued
  #   folds: number of CV folds; default is 10
  # Output
  #   A list with the ATE estimate as well as the selected tuning parameter.


  # estimate balancing propensity scores
  fit <- spBalance(
    formula = form,
    data = data,
    lambda = lambda.vals,
    opt.params = opt.params,
    tuning = tuning,
    hide.detail = TRUE,
    folds = folds,
    coefvar.r = coefvar.r
  )

  if (tuning == "all"){
    ate.vals <- matrix(data = NA_real_, nrow = length(fit$lambda), ncol = 3)
    colnames(ate.vals) <- c("HT", "Hajek", "OW")
    rownames(ate.vals) <- names(fit$lambda)
    for (k in 1:length(fit$lambda)){
      ate.vals[k,] <- ateEst(out = out, trt = fit$gam$y, pr.trt = fit$pi.hat[,k])[[2]]
    }
  } else {
    ate.vals <- ateEst(out = out, trt = fit$gam$y, pr.trt = fit$pi.hat)[[2]]
  }


  # compile results
  if (tuning == "all"){

    # determine which estimand to return
    if (est.type == "all"){
      final.est <- ate.vals
    } else if (type == "HT"){
      final.est <- ate.vals[,1]
    } else if (type == "Hajek"){
      final.est <- ate.vals[,2]
    } else if (type == "OW"){
      final.est <- ate.vals[,3]
    }
  } else {
    # determine which estimand to return
    if (est.type == "all"){
      final.est <- ate.vals
    } else if (type == "HT"){
      final.est <- ate.vals[1]
    } else if (type == "Hajek"){
      final.est <- ate.vals[2]
    } else if (type == "OW"){
      final.est <- ate.vals[3]
    }
  }

  # return final estimate and choice of lambda
  res <- list(
    ate = final.est,
    lambda = fit$lambda[[1]]
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


