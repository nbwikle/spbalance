### spbalance.R
### Nathan Wikle

### Functions used to create balancing weights for settings with unobserved
###   spatial confounding.

parseFormula <- function(form){
  # Helper function to parse a formula object and determine what
  #   model components to "balance".
  # Input
  #   form: a formulat object
  # Output
  #   A list splitting the formula into a parametric, spline, and krr component.

  ### 1. Parse formula

  # specials attribute indicates which terms are smooth / krr
  tf <- terms.formula(form, specials=c("s", "krr"))
  # labels of the model terms
  terms <- attr(tf, "term.labels")

  nt <- length(terms) # how many terms?
  pf <- "~"          # start of parametric model formula
  if (attr(tf,"response") > 0){  # start the replacement formula
    rf <- paste(as.character(attr(tf, "variables")[2]), "~", sep = "")
    krrf <- paste(as.character(attr(tf, "variables")[2]), "~", sep = "")
  } else {
    rf <- "~"
    krrf <- "~"
  }

  ### tells you which formula terms match "s" or "krr"
  sp <- attr(tf, "specials")$s  # array of indices of smooth terms
  krrp <- attr(tf, "specials")$krr  # array of indices of krr terms

  ### deals with offset (not applicable to balancing weights)
  off <- attr(tf, "offset") # location of offset in formula
  if (!is.null(off)){
    # have to remove the offset from this index list
    sp[sp > off] <- sp[sp > off] - 1
    krrp[krrp > off] <- krrp[krrp > off] - 1
  }

  ns <- length(sp) # number of smooths
  nkrr <- length(krrp)  # number of krr terms
  ks <- 1; kp <- 1; kkern <- 1 # counters for terms in the 2 formulae
  kcentered <- TRUE

  if (nt){

    for (i in 1:nt){
      if (ks <= ns && sp[ks] == i + 1){
        # it's a smooth
        st <- eval(parse(text = terms[i]))
        if (ks > 1 || kp > 1) {
          # add to smooth formula
          rf <- paste(rf, "+", terms[i], sep = "")
        } else {
          rf <- paste(rf, terms[i], sep = "")
        }
        ks <- ks + 1
      } else if (kkern <= nkrr && krrp[kkern] == i + 1){
        # it's a krr term
        krrt <- eval(parse(text = terms[i]))

        if (!krrt$cent){
          kcentered <- FALSE
        }

        if (kkern > 1) {
          # add to krr formula
          krrf <- paste(krrf, "+", terms[i], sep = "")
        } else {
          krrf <- paste(krrf, terms[i], sep = "")
        }
        kkern <- kkern + 1
      } else {
        # parametric term
        if (kp > 1) {
          # add to parametric formula
          pf <- paste(pf, "+", terms[i], sep = "")
        } else {
          pf <- paste(pf, terms[i], sep = "")
        }
        if (ks > 1 || kp > 1) {
          # add to smooth formula
          rf <- paste(rf, "+", terms[i], sep = "")
        } else {
          rf <- paste(rf, terms[i], sep="")
        }
        kp <- kp + 1
      }
    }
  }

  if (nkrr == nt){
    pf <- paste(pf, "1", sep = "")
    rf <- paste(rf, "1", sep = "")
  } else if (ns == nt){
    pf <- paste(pf, "1", sep = "")
    krrf <- paste(krrf, "1", sep = "")
  } else if (ns + nkrr == nt){
    pf <- paste(pf, "1", sep = "")
  }

  ### 2. Store parametric, krr, and smooth terms in a list

  # return parsed formula
  ret <- list(
    pftext = pf,
    pf = as.formula(pf), # parametric forumla
    smtext = rf,
    smf = as.formula(rf), # parametric + spline formula
    krrtext = krrf,
    krrf = as.formula(krrf), # krr formula only
    kcentered = TRUE
  )
  class(ret) <- "split.formula"
  ret
}


krr <- function(..., kern = "SE", kp = 1, var.inds = NULL, centered = TRUE, scale = TRUE, rff = FALSE, nf = 0){
  # Formats a kernel ridge regression (KRR) formula term into a 'krr.spec' object.
  # Input
  #   kern: type of RKHS kernel
  #   kp: kernel parameter (range parameter)
  #   var.inds: index of which variables should be included in KRR term
  #   centered: Boolean indicating if KRR should be centered at zero
  #   scale: should KRR term be scaled to have max = 1?
  #   rff: Boolean indicating if kernel should be approximated with random fourier features
  #   nf: number of RFF features
  # Output
  #   An object of class 'krr.spec'

  if (is.null(var.inds)){
    # variable names included in formula
    vars <- as.list(substitute(list(...)))[-1]
    d <- length(vars)
    # collect term names
    term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
    if (term[1] == "."){
      if (length(d) > 1){
        stop("krr(.) cannot be combined with other terms")
      }
    }
    if (d > 1){
      for (i in 2:d) {
        term[i] <- deparse(vars[[i]], backtick = TRUE, width.cutoff = 500)
        if (term[i] == "."){
          stop("krr(.) cannot be combined with other terms")
        }
      }
    }
    for (i in 1:d){
      term[i] <- attr(terms(reformulate(term[i])), "term.labels")
    }

    # check that there are no repeated terms
    if (length(unique(term)) != d){
      stop("Repeated variables as arguments of krr are not permitted")
    }

    # recreate krr call
    full.call <- paste("krr(", term[1], sep = "")
    if (d > 1){
      for (i in 2:d){
        full.call <- paste(full.call, ",", term[i], sep = "")
      }
    }
    label <- paste(full.call, ")", sep = "")
  } else {
    term <- NULL
    label <- NULL
  }

  # check that kernel function is correctly specified
  if ((kern != "SE") && (kern != "Exp")){
    stop("Kernel must be one of 'SE' or 'Exp'")
  }

  # check that kernel parameter(s) are allowed
  if (kp <= 0){
    stop("Kernel parameter must be > 0")
  }

  # return as a list
  ret <- list(
    term = term,
    kernel = kern,
    k.params = kp,
    x.inds = var.inds,
    cent = centered,
    scale = scale,
    rff = rff,
    n.features = nf,
    label = label
  )
  class(ret) <- paste("krr.spec", sep = "")
  ret
}


spBalance <- function(
  formula, data, lambda, init.params = NULL, fit.gam = FALSE, pen.int = FALSE,
  tuning = "none", folds = 10, grad.norm = "L2", coefvar.r = 0.9, bal.diff = 0.1,
  hide.details = TRUE,
  opt.params = list(tol = 1e-7, max.iter = 100, alpha = 0.5, beta = 0.5),
  ...
){
  # Estimates balancing weights for a given formula, including parametric,
  #   spline, and/or KRR terms.
  # Input
  #   formula: a formula term, eg, 'y ~ x1 + s(x2)', where 's(...)' adopts the
  #     spline formula's used in the mgcv package. Kernel ridge regression
  #     terms are also permissible, using
  #     'krr(
  #       x1, ..., xp, # variables to include in kernel evaluation
  #       kern = "SE", # type of kernel; must be one of 'SE' or 'Exp'
  #       kp = p # number of variables included before "kern" input,
  #       centered = FALSE # boolean indicating if the KRR term should be centered at zero
  #     )'
  #   data: data frame containing all variables used in the analysis
  #   lambda: a vector of tuning possible tuning parameter values
  #   init.params: initial parameter estimates; default is 'NULL'
  #   fit.gam: Boolean indicating if GAM model should be fitted in mgcv; default is 'FALSE'
  #   pen.int: Boolean indicating if the intercept should be penalized; default is 'FALSE'
  #   tuning: indicates how the tuning parameter should be chosen; options include
  #     `cv.score`, `cv.grad`, `coefvar`, `max.bal`, `min`, `all`
  #   folds: number of CV folds (if using a CV method for tuning parameter selection)
  #   grad.norm: type of Lp norm used when evaluating CV with gradient of score;
  #       options are 'L1', 'L2', or 'Linf'
  #   coefvar.r: ratio used in "coefvar" tuning parameter selection; scalar or vector valued
  #   bal.diff: maximum standardized mean difference allowed when choosing tuning param wiht "max.bal"
  #   hide.details: Boolean indicating if tuning parameter selection details should be saved
  #   opt.params: arguments passed to "optim" function call
  # Output
  #   A list containing estimated parameter values, estimated propensity scores,
  # chosen tuning parameter (lambda), and tuning parameter selection details.

  # fit model with specific lambda
  full.fit <- spBalFit(
    formula = formula,
    data = data,
    lambda = lambda,
    opt.params = opt.params,
    init.params = init.params,
    fit.gam = fit.gam,
    pen.int = pen.int
  )

  # begin compiling results
  res <- full.fit
  res$tuning <- tuning
  res$tuning.details <- NULL

  # treatment
  treat <- data[,all.vars(formula)[1]]

  # select lambda via CV
  if ((tuning == "cv.score") | (tuning == "cv.grad") | (tuning == "all")){

    cv.res <- cvLambda(
      kfold = folds,
      formula = formula,
      data = data,
      opt.params = opt.params,
      lambda.vals = lambda,
      type = tuning,
      norm = grad.norm,
      fit.gam = fit.gam,
      pen.int = pen.int
    )

    if ((tuning == "cv.score") | (tuning == "cv.grad")){
      l.choice <- which(lambda == cv.res$lambda[1])
      res$par <- full.fit$par[,l.choice]
      res$pi.hat <- full.fit$pi.hat[,l.choice]
      res$lambda <- cv.res$lambda[1]
      if (!hide.details){
        res$tuning.details <- colMeans(cv.res$cv[[1]])
      }
    } else if (tuning == "all"){
      l.cv.all <- c(
        which(lambda == cv.res$lambda[[1]]),
        which(lambda == cv.res$lambda[2])
      )
      l.cv.details <- lapply(cv.res$cv, colMeans)
      names(l.cv.details) <- c("cv.score", "cv.grad")
    }
  }

  # select lambda via coefficent of variation
  if ((tuning == "coefvar") | (tuning == "all")){

    coefv.res <- list()

    for (j in 1:length(coefvar.r)){
      coefv.res[[j]] <- coefVar(
        z = treat,
        pi.vals = full.fit$pi.hat,
        lambda = lambda,
        ratio = coefvar.r[j]
      )
    }

    l.coef.all <- rep(0, length(coefvar.r))
    l.coef.details <- coefv.res[[1]][[2]]
    for (k in 1:length(coefvar.r)){
      l.coef.all[k] <- which(lambda == coefv.res[[k]][[1]])
    }
    names(l.coef.all) <- paste("cvr", coefvar.r, sep = "")

    if (tuning == "coefvar"){
      res$par <- full.fit$par[,l.coef.all]
      res$pi.hat <- full.fit$pi.hat[,l.coef.all]
      res$lambda <- full.fit$lambda[l.coef.all]
      if (!hide.details){
        res$tuning.details <- l.coef.details
      }
    }
  }

  # select lambda via setting balance to a specific threshold
  if ((tuning == "max.bal") | (tuning == "all")){

    mbal <- maxBal(
      fit = full.fit, z = treat, data = data,
      max.diff = bal.diff, return.fit = FALSE
    )

    l.bal <- which(lambda == mbal$lambda)
    l.bal.details <- mbal$max.smd

    if (tuning == "max.bal"){
      res$par <- full.fit$par[,l.bal]
      res$pi.hat <- full.fit$pi.hat[,l.bal]
      res$lambda <- mbal$lambda
      if (!hide.details){
        res$tuning.details <- l.bal.details
      }
    }
  }

  # use minimum lambda value
  if (tuning == "min"){
    l.min <- which(lambda == min(lambda))
    res$par <-full.fit$par[,l.min]
    res$pi.hat <- full.fit$pi.hat[,l.min]
    res$lambda <- lambda[l.min]
    if (!hide.details){
      res$tuning.details <- lambda
    }
  } else if (tuning == "all"){
    l.min <- which(lambda == min(lambda))
    l.min.details <- lambda
  }

  # compile results
  if (tuning == "all"){

    l.vals <- c(l.cv.all, l.coef.all, l.bal, l.min)
    names(l.vals) <- c("cv.score", "cv.grad", paste("cvr-", coefvar.r, sep = ""), "max.bal", "min")

    # results
    res$par <- full.fit$par[,l.vals]
    res$pi.hat <- full.fit$pi.hat[,l.vals]
    res$lambda <- lambda[l.vals]
    colnames(res$par) <- names(l.vals)
    colnames(res$pi.hat) <- names(l.vals)
    names(res$lambda) <- names(l.vals)

    # return details (for debugging)
    if (!hide.details){
      res$tuning.details <- list(
        cv = l.cv.details,
        coefvar = l.coef.details,
        max.bal = l.bal.details,
        min = l.min.details
      )
    }

  }
  # return fitted model
  return(res)
}


spBalFit <- function(
  formula, data = list(), lambda, init.params = NULL, fit.gam = FALSE, pen.int = FALSE,
  opt.params = list(tol = 1e-7, max.iter = 100, alpha = 0.5, beta = 0.5),...
){
  # Estimate covariate balancing propensity scores for a given formula and
  #   vector of tuning parameter values (lambda).
  # Input
  #   formula: a formula object similar to that used by 'lm' or 'mgcv';
  #     acceptible terms include variable names, splines 's(.)', and
  #     kernel ridge regression 'krr(.)'
  #   data: data frame containing all variables used in the analysis
  #   lambda: a vector of possible tuning parameter values
  #   init.params: initial parameter estimates; default is 'NULL'
  #   fit.gam: Boolean indicating if GAM model should be fitted in mgcv; default is 'FALSE'
  #   pen.int: Boolean indicating if the intercept should be penalized; default is 'FALSE'
  #   opt.params: optimization parameters passed to 'optim'
  # Output
  #   Returns a list with parameter estimates, fitted propensity scores,
  # tuning parameters (lambda), optim convergence results, number of iterations
  # of optimization procedure, a parsed formula object, 'gam' terms (if used),
  # 'krr' terms (if used), dimensionality of each model term, and the type
  # of terms specified in the formula (intercept, fixed, smooth, krr).

  ### 1. parse formula for "krr" and "smooth" (i.e., "s") terms:
  parseform <- parseFormula(formula)

  # treatment (aka, the response)
  treat <- data[, all.vars(formula)[1]]
  # spline formula
  spline.f <- terms.formula(parseform$smf, specials=c("s"))
  # krr formula
  krr.f <- terms.formula(parseform$krrf, specials = c("k"))
  # labels of the model terms
  spline.terms <- attr(spline.f, "term.labels")
  krr.terms <- attr(krr.f, "term.labels")

  ### 2. Construct basis and penalty matrices for each formula term
  basis <- list()
  penalty <- list()
  term.type <- list()
  s.dims <- list()

  # (i) spline terms
  if (length(spline.terms) > 0){
    # extract spline info using mgcv::gam()
    gam.obj <- gam(spline.f, family = binomial, data = data, fit = fit.gam, ...)

    # basis functions
    if (fit.gam){
      spline.basis <- predict(gam.obj, type = "lpmatrix")
      # basis_func <- predict(gam.obj, newdata = data, type = "lpmatrix")
    } else {
      spline.basis <- gam.obj$X # predict(gam.fit, type = "lpmatrix")
      colnames(spline.basis) <- gam.obj$term.names
    }

    ### collect basis functions into a list

    # number of fixed effects (including intercept)
    n.fixed <- length(all.vars(gam.obj$pterms))
    n.comb <- ncol(spline.basis)

    # intercept
    basis[[1]] <- matrix(spline.basis[,1], ncol = 1)
    if (pen.int){
      penalty[[1]] <- matrix(1, nrow = 1, ncol = 1)
    } else {
      penalty[[1]] <- matrix(0, nrow = 1, ncol = 1)
    }
    s.dims[[1]] <- 1
    term.type[[1]] <- "intercept"

    # fixed effects
    if (n.fixed > 1){
      basis[[2]] <- matrix(spline.basis[,2:(n.fixed)], ncol = n.fixed - 1)
      penalty[[2]] <- diag(1, n.fixed - 1)
      s.dims[[2]] <- 2:n.fixed
      term.type[[2]] <- "fixed"
    }

    # spline terms
    if (n.comb > n.fixed){
      # number of parametric terms
      n.p <- length(basis)
      # number of smooth terms
      n.smooths <- length(gam.obj$smooth)
      for (k in 1:n.smooths){
        # grab smooth term
        sm.k <- gam.obj$smooth[[k]]
        # scaling factor
        alpha.k <- sm.k$S.scale
        # penalty matrix
        S.tilde.k <- sm.k$S[[1]]
        S.k <- S.tilde.k * alpha.k
        penalty[[k + n.p]] <- S.k
        s.dims[[k + n.p]] <- max(unlist(s.dims), na.rm = TRUE) + 1:nrow(S.k)
        basis[[k + n.p]] <- matrix(spline.basis[,s.dims[[k + n.p]]], ncol = ncol(S.k))
        term.type[[k + n.p]] <- "smooth"
      }
    }
  } else {
    gam.obj <- NULL
    if (parseform$kcentered){
      # add intercept term
      basis[[1]] <- matrix(0, nrow = length(treat), ncol = 1)
      if (pen.int){
        penalty[[1]] <- matrix(1, nrow = 1, ncol = 1)
      } else {
        penalty[[1]] <- matrix(0, nrow = 1, ncol = 1)
      }
      s.dims[[1]] <- 1
      term.type[[1]] <- "intercept"
    }
  }

  # (ii) krr terms

  # number of krr terms
  n.krr <- length(krr.terms)
  n.s <- length(basis)
  krr.obj <- list()

  if (n.krr > 0){

    for (k in 1:n.krr){
      # grab krr term
      krr.t <- eval(parse(text = krr.terms[k]))
      # create krr structures
      krr.k <- krrCreation(krr.obj = krr.t, data = data)
      krr.obj[[k]] <- krr.k

      # index
      s.dims[[k + n.s]] <- max(unlist(s.dims), na.rm = TRUE) + 1:ncol(krr.k$basis)
      # basis functions
      basis[[k + n.s]] <- krr.k$basis
      # penalty matrix
      penalty[[k + n.s]] <- krr.k$penalty
      # type of term
      if (krr.t$rff){
        term.type[[k + n.s]] <- "rff"
      } else {
        term.type[[k + n.s]] <- "krr"
      }
    }
  } else{
    krr.obj <- NULL
  }

  ### 3. Create final basis and penalty matrices

  # na.col <- which(is.na(unlist(s.dims)))
  basis.full <- do.call(cbind, basis)
  n.theta <- ncol(basis.full)

  penalty.full <- matrix(0, nrow = n.theta, ncol = n.theta)
  for (k in 1:length(penalty)){
    dims.k <- s.dims[[k]]
    penalty.full[dims.k, dims.k] <- penalty[[k]]
  }

  ### 4. Estimate theta.hat for given basis and penalty matrices
  all.types <- unique(unlist(term.type))
  krr.only <- FALSE
  if (length(all.types) == 1){
    if (all.types[[1]] == "krr"){
      krr.only <- TRUE
    }
  } else if (length(all.types) == 2){
    if ((all.types[[1]] == "intercept") & (all.types[[2]] == "krr")){
      krr.only <- TRUE
    }
  }

  # repeat estimation for all given values of lambda
  n.l <- length(lambda)
  theta.m <- matrix(data = NA_real_, nrow = n.theta, ncol = n.l)
  pi.hat <- matrix(data = NA_real_, nrow = nrow(basis.full), ncol = n.l)
  colnames(theta.m) <- colnames(pi.hat) <- paste("l=", lambda, sep = "")
  convergence <- rep(0, n.l)
  counts <- list()

  theta.k <- init.params
  for (k in 1:n.l){
    fit.k <- spbalWeights(
      theta = theta.k,
      z = treat, X = basis.full, P = penalty.full, lambda = lambda[k],
      opt.params = opt.params, krr.only = krr.only,
      centered = parseform$kcentered
    )
    theta.k <- fit.k$theta
    theta.m[,k] <- theta.k
    pi.hat[,k] <- fit.k$pi.hat
    convergence[k] <- fit.k$convergence
    counts[[k]] <- fit.k$counts
  }

  # return results for different lambda values
  results <- list(
    par = theta.m,
    pi.hat = pi.hat,
    lambda = lambda,
    convergence = convergence,
    counts = counts,
    form = parseform,
    gam = gam.obj,
    krr = krr.obj,
    dims = s.dims,
    terms = term.type
  )

  return(results)
}


spbalWeights <- function(theta = NULL, z, X, P, lambda, opt.params, krr.only = FALSE, centered = TRUE){
  # Estimates balancing weights for a given set of basis functions; uses
  #   R's 'optim' function to perform optimization.
  # Input
  #   theta: initial parameter values; default is NULL
  #   z: treatment vector (assumed to be binary)
  #   X: design matrix
  #   P: penalty matrix
  #   lambda: tuning parameter
  #   opt.params: arguments for 'optim' function call
  #   krr.only: Boolean indicating if there is only a KRR term in the model;
  #     if so, will use more efficient 'balanceKRRalg(.)' function
  #   centered: Boolean indicating if KRR term has been centered at zero
  # Output
  #   A list containing parameter estimates, propensity score estimates,
  # optim convergence results, and number of optimization iterations until
  # until convergence ('counts').

  # initialize theta
  if (is.null(theta)){
    theta.0 <- rep(0, ncol(X))
  } else {
    theta.0 <- theta
  }


  if (krr.only){
    # use KRR algorithm if KRR is only term in model
    fit <- balanceKRRalg(
      theta = theta.0, z = z, X = X, P = P,
      lambda = lambda,
      centered = centered,
      params = opt.params
    )
  } else {
    # use optim for all other models
    fit <- optim(
      par = theta.0, fn = lossFunction, gr = gradFunction, method = "BFGS",
      z = z, X = X, P = P, lambda = lambda,
      control = list(maxit = opt.params$max.iter)
    )
  }

  # return fitted model
  res <- list(
    theta = fit$par,
    pi.hat = expit(X %*% fit$par),
    convergence = fit$convergence,
    counts = fit$counts
  )
  res
}


balanceKRRalg <- function(
  theta, z, X, P, lambda, centered = TRUE,
  params = list(tol = 1e-7, max.iter = 20, alpha = 0.5, beta = 0.5, pen.fixed = FALSE)
){
  # Implements iterative Newton's method for fitting balancing weights to
  #   a single KRR term.
  # Input
  #   theta: initial parameter values; default is NULL
  #   z: treatment vector (assumed to be binary)
  #   X: design matrix
  #   P: penalty matrix
  #   lambda: tuning parameter
  #   centered: Boolean indicatin if KRR term is centered at zero
  #   params: arguments for optimization procedure
  #   pen.fixed: Boolean indicating if fixed effects should be penalized
  # Output
  #   A list containing parameter estimates, convergence results, and number
  # of Fisher scoring iterations ('counts').

  # initialize variables
  n.r <- nrow(X)
  change <- Inf
  k <- 1

  # initialize theta
  theta.k <- theta

  # rescale lambda
  lambda.t <- lambda / n.r

  if (centered){
    Sigma <- P[-1,-1]
  }

  # objective function
  L.k <- balanceObjFun(z, X, theta.k, lambda.t, P)

  while(change > params$tol){
    if (k == params$max.iter){
      change <- 0
    } else {
      if (centered){
        pi.val <- expit(X %*% theta.k)
        tp <- theta.k
        tp[1] <- mean(((z - pi.val) / pi.val / (1 - pi.val))) + theta.k[1]
        M.val <- Sigma + diag(n.r * lambda.t, nrow(Sigma))
        z.new <- as.vector(Sigma %*% theta.k[-1] + ((z - pi.val) / pi.val / (1 - pi.val)))
        tp[-1] <- solve(M.val, z.new)
      } else {
        # IRLS:
        pi.val <- expit(X %*% theta.k)
        z.new <- as.vector(X %*% theta.k + ((z - pi.val) / pi.val / (1 - pi.val)))
        M.val <- X + diag(n.r * lambda.t, nrow(X))
        tp <- solve(M.val, z.new)
      }

      # Determine if objective was met:
      L.new <- balanceObjFun(z, X, tp, lambda.t, P)
      L.diff <- L.new - L.k

      D.mat <- crossprod(X, ((z - pi.val) / pi.val / (1 - pi.val))) / n.r - lambda.t * (P %*% theta.k)
      penalty <- crossprod(D.mat, (tp - theta.k))

      if (L.diff > -params$alpha * penalty){
        linesearch <- TRUE
        t.val <- params$beta
      } else {
        theta.new <- tp
        linesearch <- FALSE
      }

      # Backtracking line search
      while(linesearch){
        theta.new <- (1 - t.val) * theta.k + t.val * tp
        L.new <- balanceObjFun(z, X, theta.new, lambda.t, P)
        L.diff <- L.new - L.k

        if (L.diff > -params$alpha * t.val * penalty){
          t.val <- t.val * params$beta
        } else {
          linesearch <- FALSE
        }
      }

      # Assess convergence:
      change <- sqrt(sum((theta.new - theta.k)^2))
      theta.k <- theta.new
      L.k <- L.new
      k <- k + 1
    }
  }

  if (k < params$max.iter){
    list(
      par = c(theta.k),
      counts = paste("Fisher scoring iterations: ", k, sep = ""),
      convergence = 0
    )
  } else {
    list(
      par = c(theta.k),
      counts = paste("Max iter (", params$max.iter, ") reached.", sep = ""),
      convergence = 1
    )
  }
}


expit <- function(x){
  # Inverse link-function.
  # Input
  #   x: probability on logit-scale
  # Output
  #   Probability
  exp(-log1p(exp(-x)))
}


balanceObjFun <- function(y, X, theta, lambda = 0, P = NULL){
  # The balancing weights objective function from Zhao (2019).
  # Input
  #   y: treatment vector (assumed to be binary)
  #   X: design matrix (balance obtained with respect to column vectors)
  #   theta: parameter vector
  #   lambda: tuning parameter; default is 0 (no penalty)
  #   P: penalty matrix; default is NULL (no penalty)
  # Output
  #   The estimated balancing score (i.e., objective function).

  n <- length(y)
  X.theta <- X %*% theta
  score <- sum(y * (X.theta - exp(-X.theta)) - (1-y)*(X.theta + exp(X.theta)) - 1) / n
  # objective function
  if (lambda > 0){
    -score + lambda * crossprod(theta, P %*% theta) / 2
  } else{
    -score
  }
}

lossFunction <- function(theta, z, X, P, lambda = 0){
  # Evaluates loss function; returns Inf if any propensity scores are estimated
  #   to be exactly 0 or 1.
  # Input
  #   theta: parameter vector
  #   z: treatment vector (assumed to be binary)
  #   X: design matrix (balance obtained with respect to column vectors)
  #   P: penalty matrix
  #   lambda: tuning parameter; default is 0 (no penalty)
  # Output
  #   Loss function evaluated for specific theta vector.

  # initialize variables
  n.r <- nrow(X)
  # rescale lambda
  lambda.t <- lambda / n.r
  # make sure theta is within appropriate bounds
  if (max(expit(X %*% theta) == 1)){
    L.k <- Inf
  } else if (max(expit(X %*% theta) == 0)){
    L.k <- Inf
  } else {
    L.k <- balanceObjFun(z, X, theta, lambda.t, P)
  }

  return(L.k)
}

gradFunction <- function(theta, z, X, P, lambda = 0){
  # The gradient of the loss function (used by 'optim').
  # Input
  #   theta: parameter vector
  #   z: treatment vector (assumed to be binary)
  #   X: design matrix (balance obtained with respect to column vectors)
  #   P: penalty matrix
  #   lambda: tuning parameter; default is 0 (no penalty)
  # Output
  #   Loss function gradient for a given theta.

  # initialize variables
  n.r <- nrow(X)
  # rescale lambda
  lambda.t <- lambda / n.r

  pi.hat <- expit(X %*% theta)
  new.z <- (z - pi.hat) / pi.hat / (1 - pi.hat)
  grad.calc <- -crossprod(X, new.z) / n.r + lambda.t * (P %*% theta)

  return(grad.calc)
}


krrCreation <- function(krr.obj, data){
  # Takes a krr object and returns covariance matrix or RFF approximation.
  # Input
  #   krr.obj: a kernel ridge regression (krr) object; created when formula
  #     was originally parsed.
  #   data: data frame containing all variables used in the analysis
  # Output
  #   A list with the KRR basis vectors, penalty matrix, standardized
  # covariate matrix (X), Z (sampled random fourier features),
  # w (frequency samples for RFF implementation), b (shift parameters in RFF
  # implementation).

  # 1. create X matrix
  if (is.null(krr.obj$term)){
    # use x.index formulation
    X <- data[,krr.obj$x.inds]
  } else {
    # use terms to create X matrix
    X <- data[,krr.obj$term]
  }

  # 2. standardize X
  if (krr.obj$scale){
    x.m <- scale(X)
  } else {
    x.m <- as.matrix(X)
  }

  # 3. distance matrix
  d.mat <- as.matrix(dist(x.m, diag = TRUE, upper = TRUE))
  Z.m <- NULL; omega <- NULL; b <- NULL

  # 4. create covariance matrix
  if (krr.obj$rff){
    # uses random fourier features
    r.features <- featureGen(
      n.features = krr.obj$n.features, X = x.m,
      k.params = krr.obj$k.params, type = krr.obj$kernel
    )
    Z.m <- r.features$Z
    omega <- r.features$w
    b <- r.features$b

    # basis functions
    if (krr.obj$cent){
      # center kernel
      basis <- base::scale(Z.m, scale = FALSE)
    } else {
      basis <- Z.m
    }
    # penalty matrix
    penalty <- diag(ncol(basis))
  } else {
    # full covariance matrix
    basis <- kernelEval(
      d.mat, k.params = krr.obj$k.params,
      type = krr.obj$kernel, centered = krr.obj$cent
    )

    # penaltymatrix
    penalty <- basis
  }

  # return kernel evaluation
  ret <- list(
    basis = basis,
    penalty = penalty,
    krr.obj = krr.obj,
    X = x.m,
    Z = basis,
    w = omega,
    b = b
  )
  ret
}

kernelEval <- function(d, k.params, type = "SE", centered = TRUE){
  # Evaluates kernel matrix for a given distance matrix.
  # Input
  #   d: distance matrix
  #   k.params: kernel parameters (range (aka length scale) parameter for
  #     "SE" and "Exp"; range and smoothness for "Matern")
  #   type: one of "SE", "Exp", or "Matern"
  #   centered: whether to center the KRR term at 0.
  # Output
  #   the kernel matrix


  if (type == "SE"){
    K.m <- exp(-d^2 / 2 / k.params)
  } else if (type == "Exp"){
    K.m <- exp(-abs(d) / 2 / k.params)
  } else if (type == "Matern"){
    # rho <- tuning[1]
    # nu <- tuning[2]
    #
    # term1 <- 2^(1-nu) / gamma(nu)
    # term2 <- (2 * sqrt(nu) * abs(d) / rho)^(nu)
    # term3 <- besselK(2 * sqrt(nu) * abs(d) / rho, nu = nu)
    # K.m <- term1 * term2 * term3
    # diag(K.m) <- 1
  }
  if (centered){
    mat.dim <- nrow(K.m)
    K.m <- K.m -
      matrix(colMeans(K.m), nrow = mat.dim, ncol = mat.dim, byrow = TRUE) -
      matrix(colMeans(K.m), nrow = mat.dim, ncol = mat.dim) +
      matrix(mean(colMeans(K.m)), nrow = mat.dim, ncol = mat.dim)
    # One.m <- matrix(1, nrow = mat.dim, ncol = mat.dim) / One.m
    # K.m <- K.m - One.m %*% K.m - K.m %*% One.m + One.m %*% (K.m %*% One.m)
  }
  return(K.m)
}


featureGen <- function(n.features, X, k.params, type = "SE"){
  # Generates random Fourier features (rff) for a given covariance specification.
  # Input
  #   n.features: number of rff's to generate
  #   X: covariate matrix used to evaluate distance
  #   k.params: kernel parameters (range (aka length scale) parameter for
  #     "SE" and "Exp"
  #   type: one of "SE" or "Exp"
  # Output
  #   A list with a matrix of random features (Z), sampled frequencies (w),
  # and sampled shift parameters (b).

  # sample omega from Fourier density
  if (is.matrix(X)){
    x.dim <- ncol(X)
    X.m <- X
  } else {
    X.m <- matrix(X, nrow = length(X), ncol = 1)
    x.dim = 1
  }

  if (type == "SE"){
    fourier.var <- 1 / k.params[1]
    omega <- matrix(rnorm(n = n.features * x.dim, mean = 0, sd = sqrt(fourier.var)),
         nrow = n.features, ncol = x.dim, byrow = TRUE)
  } else if (type == "Exp"){
    # sample from multivariate normal
    rho <- k.params[1] / 2
    mvn.sd <- sqrt(2) / pi / rho / 2
    mvn.samps <- matrix(rnorm(n = n.features * x.dim, mean = 0, sd = mvn.sd),
                        nrow = n.features, ncol = x.dim, byrow = TRUE)
    # sample from chi-square
    chi.samps <- rchisq(n = n.features, df = 1)
    # create frequency samples from multivariate t distribution with 1 degree of freedom
    omega <- mvn.samps / sqrt(chi.samps)
  }

  # sample shift parameter, b
  b <- runif(n = n.features, min = 0, max = 2 * pi)

  # create features
  f.mu <- X.m %*% t(omega) + rep(b, each = nrow(X.m))
  Z.x <- sqrt(2) * cos(f.mu) / sqrt(n.features)

  # return features (and sampled omega and b values)
  list(
    Z = Z.x,
    w = omega,
    b = b
  )
}












