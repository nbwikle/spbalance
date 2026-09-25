#' spbalance: covariate balancing propensity scores with spatial confounding
#'
#' Estimates propensity scores with a loss that prioritises covariate balance
#' ([spBalance()]), allowing spline and kernel terms in space to adjust for
#' unmeasured spatial confounding, and estimates the average treatment effect
#' with iid, spatial HAC or spatial block bootstrap standard errors
#' ([spbalATE()]).
#'
#' @keywords internal
#' @importFrom stats as.formula binomial coef dist fitted formula gaussian
#'   model.matrix nobs optim predict quantile rchisq reformulate rnorm runif sd
#'   terms terms.formula var
#' @importFrom mgcv gam PredictMat s
#' @importFrom gstat variogram fit.variogram vgm
#' @importFrom conleyreg conleyreg
#' @importFrom units set_units
"_PACKAGE"

# conleyreg::conleyreg() calls `units<-` on the distance matrix without loading
# the units package; importing from units ensures its S3 methods are registered
# whenever spbalance is loaded.
