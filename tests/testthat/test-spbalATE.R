sim_data <- function(n.side = 15, seed = 1) {
  set.seed(seed)
  d <- expand.grid(x = seq_len(n.side), y = seq_len(n.side))
  d$u   <- sin(d$x / 4) + cos(d$y / 3)
  d$x1  <- rnorm(nrow(d))
  d$trt <- rbinom(nrow(d), 1, plogis(0.5 * d$x1 + d$u))
  d$out <- 1 + 0.5 * d$trt + d$x1 + 2 * d$u + rnorm(nrow(d))
  d
}

fit_data <- function(d) spBalance(trt ~ x1 + s(x, y, k = 20), data = d, lambda = 1)

test_that("spBalance stores what the influence functions and refits need", {
  d <- sim_data()
  fit <- fit_data(d)
  expect_s3_class(fit, "bal")
  expect_equal(dim(fit$X), c(nrow(d), length(fit$par)))
  expect_equal(dim(fit$P), c(length(fit$par), length(fit$par)))
  expect_equal(deparse(fit$spec$formula), "trt ~ x1 + s(x, y, k = 20)")
  # refitting with the stored settings reproduces the fit
  expect_equal(spbalance:::refitBalance(fit, d)$par, fit$par)
})

test_that("spbalATE agrees with est.HAC for every estimator", {
  d <- sim_data()
  fit <- fit_data(d)
  mu.f <- out ~ trt + x1 + s(x, y, k = 20)
  a <- suppressWarnings(spbalATE(fit, d, "out", estimator = c("HT", "Hajek", "aIPTW"),
                                 se = c("iid", "hac"), mu.formula = mu.f))
  expect_s3_class(a, "spbalATE")
  expect_equal(a$estimates$estimator, c("HT", "Hajek", "aIPTW"))
  mu.fit <- mgcv::gam(mu.f, data = d)
  for (i in 1:3) {
    e <- a$estimates$estimator[i]
    h <- suppressWarnings(est.HAC(d$out, d$trt, fit, d[, c("x", "y")], est = e,
                                  mu.fit = mu.fit, trt.name = "trt"))
    expect_equal(a$estimates$estimate[i], h$tau)
    expect_equal(a$estimates$se.hac[i], h$se.hac)
    expect_equal(a$estimates$se.iid[i], h$se.iid)
  }
  expect_true(all(a$estimates$se.iid > 0))
})

test_that("spbalATE returns only the requested standard errors", {
  d <- sim_data()
  fit <- fit_data(d)
  a <- spbalATE(fit, d, "out", estimator = "HT", se = "iid")
  expect_named(a$estimates, c("estimator", "estimate", "se.iid"))
  expect_output(print(a), "Average treatment effect of `trt` on `out`")
})

test_that("spbalATE checks its inputs", {
  d <- sim_data()
  fit <- fit_data(d)
  expect_error(spbalATE(fit, d, "out", estimator = "aIPTW"), "mu.formula")
  expect_error(spbalATE(fit, d, "out", estimator = "aIPTW", mu.formula = out ~ x1),
               "must include the treatment")
  expect_error(spbalATE(fit, d[-1, ], "out"), "number of rows differs")
  expect_error(spbalATE(fit, d, "nope"), "not found")
  expect_error(spbalATE(list(), d, "out"), "spBalance\\(\\) fit")
  fit.all <- spBalance(trt ~ x1 + s(x, y, k = 20), data = d, lambda = c(0.5, 1, 2),
                       tuning = "all", coefvar.r = 0.9, folds = 3)
  expect_error(spbalATE(fit.all, d, "out"), "several tuning choices")
})

test_that("block bootstrap standard errors come from spbootstrap", {
  skip_if_not_installed("spbootstrap")
  d <- sim_data(12)
  fit <- fit_data(d)
  set.seed(3)
  a <- spbalATE(fit, d, "out", estimator = c("HT", "Hajek"), se = c("iid", "boot"),
                boot = list(n.boot = 8, block.l = 3))
  expect_true(all(is.finite(a$estimates$se.boot) & a$estimates$se.boot > 0))
  expect_equal(unname(a$block.l), c(3, 3))
  expect_equal(nrow(a$boot$t), 8)
  expect_output(print(a), "Block bootstrap: 8 replicates")
})

test_that("ateEst reproduces the known-propensity IPTW estimates", {
  set.seed(4)
  z <- rbinom(500, 1, 0.4); y <- 2 * z + rnorm(500)
  e <- ateEst(y, z, rep(0.4, 500))
  expect_equal(e$ate[["HT"]], mean(z * y / 0.4) - mean((1 - z) * y / 0.6))
  expect_equal(e$ate[["Hajek"]], mean(y[z == 1]) - mean(y[z == 0]))
})

test_that("cross-validated tuning works with a krr() term", {
  set.seed(5)
  d <- data.frame(x1 = rnorm(120), x2 = rnorm(120))
  d$trt <- rbinom(120, 1, plogis(d$x1 - d$x2))
  fit <- spBalance(trt ~ krr(x1, x2, kern = "SE", kp = 2), data = d,
                   lambda = c(0.5, 1, 2), tuning = "cv.score", folds = 3)
  expect_s3_class(fit, "bal")
  expect_true(fit$lambda %in% c(0.5, 1, 2))
})

test_that("formula terms can use variables from the calling environment", {
  d <- sim_data(12)
  f <- function(n.knots) spBalance(trt ~ x1 + s(x, y, k = n.knots), data = d, lambda = 1)
  fit <- f(15)
  expect_s3_class(fit, "bal")
  expect_equal(ncol(fit$X), length(fit$par))
  # the bootstrap refit (a different call frame) still finds n.knots
  expect_equal(spbalance:::refitBalance(fit, d)$par, fit$par)
  # s() works without mgcv attached
  expect_false("package:mgcv" %in% search())
})
