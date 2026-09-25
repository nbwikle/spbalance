### ate-est.R
### IPTW estimates of the ATE from estimated propensity scores.

#' IPTW estimates of the average treatment effect
#'
#' Horvitz-Thompson, Hajek and overlap-weight estimates of the treated and
#' control means and the ATE, from given propensity scores.
#'
#' @param out outcome vector.
#' @param trt binary treatment vector.
#' @param pr.trt estimated propensity scores.
#' @return A list with `mu`, a matrix of estimated control (`mu0`) and treated
#'   (`mu1`) means for each estimator, and `ate`, the corresponding ATE
#'   estimates (`HT`, `Hajek`, `OW`).
#' @seealso [spbalATE()] for estimates with standard errors.
#' @export
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
