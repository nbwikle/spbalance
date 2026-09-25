### spbalATE.R
### Average treatment effect estimates and standard errors from a spBalance fit.

#' Average treatment effect with standard errors
#'
#' Estimates the average treatment effect (ATE) by inverse probability of
#' treatment weighting with propensity scores from [spBalance()], together
#' with iid, spatial HAC and/or spatial block bootstrap standard errors.
#'
#' The iid and HAC standard errors are based on the estimated influence
#' functions ([HTInfluence()]), which account for estimation of the propensity
#' score. The block bootstrap refits the propensity score model (and outcome
#' model, for `"aIPTW"`) on each replicate, with resampled observations moved
#' to their new locations, using \pkg{spbootstrap}.
#'
#' @param fit a [spBalance()] fit (with a single tuning parameter, i.e. not
#'   `tuning = "all"`).
#' @param data the data frame used to fit `fit`, including the outcome and
#'   coordinate columns.
#' @param outcome name of the outcome column in `data`.
#' @param estimator one or more of `"HT"` (Horvitz-Thompson), `"Hajek"` and
#'   `"aIPTW"` (augmented IPTW; requires `mu.formula`).
#' @param se one or more of `"iid"`, `"hac"` and `"boot"`.
#' @param coords names of the x and y coordinate columns in `data`.
#' @param mu.formula outcome model formula for `"aIPTW"`, including the
#'   treatment, e.g. `y ~ trt + s(x, y, k = 100)`; fitted with [mgcv::gam()].
#' @param mu.family family for the outcome model.
#' @param hac.args list of further arguments to [HAC.se()], e.g.
#'   `list(k.type = "uniform")`.
#' @param boot list of block bootstrap settings (used when `"boot"` is in
#'   `se`):
#'   * `n.boot`: number of replicates (default 200);
#'   * `block.l`: block side length (grid cells for `type = "grid"`, coordinate
#'     units otherwise); if `NULL` (default), chosen for each estimator by
#'     `spbootstrap::blockLength()` with `n.pilot` replicates per pilot block
#'     length (default 100);
#'   * `type`: `"grid"` (default; data on a complete regular grid), `"point"`
#'     (points in a rectangle `box.x` x `box.y`) or `"polygon"` (points in the
#'     polygon `boundary`);
#'   * `shift`, `apply.fun`, `box.x`, `box.y`, `boundary`, `c1`, `c2`: passed
#'     to `spbootstrap::spboot()` / `spbootstrap::blockLength()`.
#' @return An object of class `"spbalATE"`: a list with `estimates`, a data
#'   frame with one row per estimator (`estimate` and the requested standard
#'   errors), and, when the bootstrap is used, `boot` (the
#'   `spbootstrap::spboot()` result) and `block.l`.
#' @examples
#' \dontrun{
#' fit <- spBalance(trt ~ x1 + s(x, y, k = 100), data = df,
#'                  lambda = c(0.1, 1, 10), tuning = "coefvar")
#' spbalATE(fit, df, outcome = "y", estimator = c("HT", "Hajek"),
#'          se = c("iid", "hac"))
#' }
#' @export
spbalATE <- function(
  fit, data, outcome, estimator = c("HT", "Hajek"), se = c("iid", "hac"),
  coords = c("x", "y"), mu.formula = NULL, mu.family = stats::gaussian(),
  hac.args = list(), boot = list()
) {

  ### 1. Check inputs

  if (!inherits(fit, "bal")) stop("`fit` must be a spBalance() fit.", call. = FALSE)
  if (is.null(fit$spec)) {
    stop("`fit` was created by an older version of spBalance(); please refit it.", call. = FALSE)
  }
  if (is.matrix(fit$pi.hat) && ncol(fit$pi.hat) > 1) {
    stop("`fit` holds several tuning choices (tuning = \"all\"); refit with a single tuning method.",
         call. = FALSE)
  }
  estimator <- match.arg(estimator, c("HT", "Hajek", "aIPTW"), several.ok = TRUE)
  se <- match.arg(se, c("iid", "hac", "boot"), several.ok = TRUE)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  missing.cols <- setdiff(c(outcome, coords), names(data))
  if (length(missing.cols) > 0) {
    stop(sprintf("Column(s) not found in `data`: %s", paste(missing.cols, collapse = ", ")),
         call. = FALSE)
  }
  if (length(fit$pi.hat) != nrow(data)) {
    stop("`data` must be the data used to fit `fit` (the number of rows differs).", call. = FALSE)
  }
  trt.name <- all.vars(fit$spec$formula)[1]
  if ("aIPTW" %in% estimator) {
    if (is.null(mu.formula)) {
      stop("The \"aIPTW\" estimator needs an outcome model: supply `mu.formula`.", call. = FALSE)
    }
    if (!trt.name %in% all.vars(mu.formula)) {
      stop(sprintf("`mu.formula` must include the treatment, `%s`.", trt.name), call. = FALSE)
    }
  }

  ### 2. Estimates with iid and HAC standard errors

  y <- data[[outcome]]
  z <- data[[trt.name]]
  mu.fit <- if ("aIPTW" %in% estimator) fitOutcome(mu.formula, mu.family, data) else NULL

  rows <- lapply(estimator, function(e) {
    inf <- influence(e, y, z, fit, mu.fit, trt.name, ate.only = FALSE)
    out <- data.frame(estimator = e, estimate = inf$tau)
    if ("iid" %in% se) out$se.iid <- sqrt(mean(inf$phi^2) / length(inf$phi))
    if ("hac" %in% se) {
      out$se.hac <- do.call(HAC.se, c(
        list(phi = inf$phi, df = data.frame(data[, coords]), coords = coords), hac.args
      ))
    }
    out
  })
  estimates <- do.call(rbind, rows)

  result <- list(estimates = estimates, n = nrow(data), trt.name = trt.name,
                 outcome = outcome)

  ### 3. Block bootstrap standard errors

  if ("boot" %in% se) {
    b <- bootATE(fit, data, outcome, trt.name, estimator, coords, mu.formula, mu.family, boot)
    result$estimates$se.boot <- unname(b$fit$se[estimator])
    result$boot    <- b$fit
    result$block.l <- b$block.l
  }

  rownames(result$estimates) <- NULL
  structure(result, class = "spbalATE")
}


#' @export
print.spbalATE <- function(x, digits = 4, ...) {
  cat(sprintf("Average treatment effect of `%s` on `%s` (n = %d)\n\n",
              x$trt.name, x$outcome, x$n))
  print(x$estimates, digits = digits, row.names = FALSE, ...)
  if (!is.null(x$block.l)) {
    cat(sprintf("\nBlock bootstrap: %d replicates, block length %s\n", x$boot$n.boot,
                paste(sprintf("%s = %g", names(x$block.l), x$block.l), collapse = ", ")))
  }
  invisible(x)
}


influence <- function(estimator, y, z, fit, mu.fit, trt.name, ate.only) {
  # Influence function (and ATE) of one estimator.
  switch(estimator,
    HT    = HTInfluence(out = y, trt = z, fit = fit, ate.only = ate.only),
    Hajek = HajekInfluence(out = y, trt = z, fit = fit, ate.only = ate.only),
    aIPTW = aIPTWInfluence(out = y, trt = z, trt.name = trt.name, ps.fit = fit,
                           mu.fit = mu.fit, ate.only = ate.only)
  )
}


fitOutcome <- function(mu.formula, mu.family, data) {
  mgcv::gam(mu.formula, family = mu.family, data = data)
}


refitBalance <- function(fit, data) {
  # Refit a spBalance model, with the same settings, to new data.
  spec <- fit$spec
  do.call(spBalance, c(list(data = data), spec[setdiff(names(spec), "dots")], spec$dots))
}


bootATE <- function(fit, data, outcome, trt.name, estimator, coords,
                    mu.formula, mu.family, boot) {
  # Spatial block bootstrap of the ATE estimators via spbootstrap.
  if (!requireNamespace("spbootstrap", quietly = TRUE)) {
    stop("Bootstrap standard errors need the spbootstrap package; please install it.",
         call. = FALSE)
  }
  opts <- utils::modifyList(
    list(n.boot = 200, block.l = NULL, n.pilot = 100, type = "grid", shift = TRUE,
         apply.fun = lapply, box.x = NULL, box.y = NULL, boundary = NULL,
         c1 = 1, c2 = 0.5),
    boot
  )

  # re-estimate everything on a replicate
  statistic <- function(d) {
    fit.d <- refitBalance(fit, d)
    mu.d  <- if ("aIPTW" %in% estimator) fitOutcome(mu.formula, mu.family, d) else NULL
    vapply(estimator, function(e) {
      influence(e, d[[outcome]], d[[trt.name]], fit.d, mu.d, trt.name, ate.only = TRUE)$tau
    }, numeric(1))
  }

  common <- list(data = data, statistic = statistic, type = opts$type, coord.names = coords,
                 box.x = opts$box.x, box.y = opts$box.y, boundary = opts$boundary,
                 shift = opts$shift, apply.fun = opts$apply.fun)

  block.l <- opts$block.l
  if (is.null(block.l)) {
    bl <- do.call(spbootstrap::blockLength,
                  c(common, list(n.boot = opts$n.pilot, c1 = opts$c1, c2 = opts$c2)))
    block.l <- bl$block.l
  }
  fit.b <- do.call(spbootstrap::spboot, c(common, list(n.boot = opts$n.boot, block.l = block.l)))
  list(fit = fit.b, block.l = fit.b$block.l)
}
