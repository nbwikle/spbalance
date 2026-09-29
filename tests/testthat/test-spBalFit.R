factor_data <- function(n.side = 15, seed = 2) {
  set.seed(seed)
  d <- expand.grid(x = seq_len(n.side), y = seq_len(n.side))
  d$x1 <- rnorm(nrow(d))
  d$f  <- factor(sample(c("a", "b", "c", "d"), nrow(d), replace = TRUE))
  d$trt <- rbinom(nrow(d), 1, plogis(0.5 * d$x1 + (d$f == "b") - sin(d$x / 4)))
  d
}

test_that("factor covariates are split correctly into fixed and smooth blocks", {
  d <- factor_data()
  fit <- spBalance(trt ~ x1 + f + s(x, y, k = 20), data = d, lambda = 1)
  G <- mgcv::gam(trt ~ x1 + f + s(x, y, k = 20), family = binomial, data = d, fit = FALSE)

  # intercept, x1, and three dummy columns for f are fixed effects
  expect_equal(fit$dims[[2]], 2:G$nsdf)
  # the smooth block is exactly mgcv's smooth coefficients
  expect_equal(fit$dims[[3]], G$smooth[[1]]$first.para:G$smooth[[1]]$last.para)
  expect_equal(length(fit$par), ncol(G$X))
})

test_that("a factor covariate gives the same fit as its numeric dummy columns", {
  d <- factor_data()
  dm <- model.matrix(~ f, d)[, -1]
  colnames(dm) <- c("fb", "fc", "fd")
  d <- cbind(d, dm)

  fit.f <- spBalance(trt ~ x1 + f + s(x, y, k = 20), data = d, lambda = 1)
  fit.d <- spBalance(trt ~ x1 + fb + fc + fd + s(x, y, k = 20), data = d, lambda = 1)
  expect_equal(as.numeric(fit.f$pi.hat), as.numeric(fit.d$pi.hat), tolerance = 1e-6)
})

pen_data <- function(n.side = 15, seed = 3) {
  set.seed(seed)
  d <- expand.grid(x = seq_len(n.side), y = seq_len(n.side))
  d$x1 <- rnorm(nrow(d))
  d$x2 <- runif(nrow(d))
  d$trt <- rbinom(nrow(d), 1, plogis(0.5 * d$x1 - d$x2 + sin(d$x / 4)))
  d
}

pen_form <- trt ~ x1 + x2 + s(x, y, k = 20)

test_that("default penalties are the identity ridge and mgcv's spline scaling", {
  d <- pen_data()
  fit <- spBalance(pen_form, data = d, lambda = 1)
  G <- mgcv::gam(pen_form, family = binomial, data = d, fit = FALSE)
  sm <- G$smooth[[1]]
  expect_equal(fit$P[fit$dims[[2]], fit$dims[[2]]], diag(2))
  expect_equal(fit$P[fit$dims[[3]], fit$dims[[3]]], sm$S[[1]] * sm$S.scale)
  expect_equal(fit$spec$pen.fixed, TRUE)
  expect_equal(fit$spec$pen.scale, "mgcv")
})

test_that("pen.fixed = FALSE leaves measured confounders unpenalised and exactly balanced", {
  d <- pen_data()
  opt <- list(tol = 1e-10, max.iter = 1000, alpha = 0.5, beta = 0.5)
  fit <- spBalance(pen_form, data = d, lambda = 10, pen.fixed = FALSE, opt.params = opt)
  expect_equal(fit$P[fit$dims[[2]], fit$dims[[2]]], matrix(0, 2, 2))
  p <- as.numeric(fit$pi.hat)
  for (v in c("x1", "x2")) {
    m1 <- sum(d$trt * d[[v]] / p) / sum(d$trt / p)
    m0 <- sum((1 - d$trt) * d[[v]] / (1 - p)) / sum((1 - d$trt) / (1 - p))
    # exact up to the optimizer's convergence tolerance
    expect_equal(m1, m0, tolerance = 1e-4)
  }
})

test_that("pen.scale = 'equal' is invariant to covariate units", {
  d <- pen_data()
  d2 <- d
  d2$x1 <- 1000 * d$x1
  lam <- c(0.1, 1, 10)
  f1 <- spBalance(pen_form, data = d, lambda = lam, pen.scale = "equal")
  f2 <- spBalance(pen_form, data = d2, lambda = lam, pen.scale = "equal")
  expect_equal(f1$pi.hat, f2$pi.hat, tolerance = 1e-8)
  # the default scaling is not
  g1 <- spBalance(pen_form, data = d, lambda = lam)
  g2 <- spBalance(pen_form, data = d2, lambda = lam)
  expect_gt(max(abs(g1$pi.hat - g2$pi.hat)), 1e-2)
})

test_that("pen.scale = 'equal' is approximately invariant to coordinate units", {
  # not exact: mgcv's truncated thin plate basis changes slightly with the
  #   scale of the coordinates
  d <- pen_data()
  d2 <- d
  d2$x <- 50 * d$x
  d2$y <- 50 * d$y
  lam <- c(0.1, 1, 10)
  f1 <- spBalance(pen_form, data = d, lambda = lam, pen.scale = "equal")
  f2 <- spBalance(pen_form, data = d2, lambda = lam, pen.scale = "equal")
  expect_lt(max(abs(f1$pi.hat - f2$pi.hat)), 0.01)
  g1 <- spBalance(pen_form, data = d, lambda = lam)
  g2 <- spBalance(pen_form, data = d2, lambda = lam)
  expect_gt(max(abs(g1$pi.hat - g2$pi.hat)), 0.05)
})

test_that("pen.scale = 'equal' gives the spline unit average prior variance", {
  d <- pen_data()
  fit <- spBalance(pen_form, data = d, lambda = 1, pen.scale = "equal")
  B <- fit$X[, fit$dims[[3]]]
  P <- fit$P[fit$dims[[3]], fit$dims[[3]]]
  # prior covariance of g = B theta, with theta ~ N(0, P^{-1})
  expect_equal(mean(diag(B %*% solve(P, t(B)))), 1 + 2, tolerance = 1e-6)
  # parametric block is the covariate variances
  expect_equal(diag(fit$P[fit$dims[[2]], fit$dims[[2]]]), c(var(d$x1), var(d$x2)))
})

test_that("penalty options are kept when a fit is refitted", {
  d <- pen_data()
  fit <- spBalance(pen_form, data = d, lambda = 1, pen.fixed = FALSE, pen.scale = "equal")
  expect_equal(spbalance:::refitBalance(fit, d)$par, fit$par)
  expect_equal(spbalance:::refitBalance(fit, d)$P, fit$P)
})
