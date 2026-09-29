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
