# spbalance

*Author: Nathan Wikle*

`spbalance` estimates propensity scores with a loss that prioritises covariate
balance (the tailored loss of Zhao, 2019), using spline and kernel ridge
regression terms in space to adjust for unmeasured spatial confounding. It then
estimates the average treatment effect (ATE) by inverse probability of
treatment weighting, with **iid**, **spatial HAC** or **spatial block bootstrap**
standard errors.

The package accompanies the paper "Balancing weights in settings with
unobserved spatial confounding" (working title).

## Installation

``` r
# from the repository root
devtools::install()
```

Block bootstrap standard errors additionally need the
`spbootstrap` package.

## Usage

``` r
library(spbalance)

# 1. balancing propensity scores; s(x, y) balances a spatial smooth, which
#    adjusts for a spatially smooth unmeasured confounder
fit <- spBalance(
  formula = treat ~ x1 + x2 + s(x, y, k = 300),
  data = df,
  lambda = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 5, 10),
  tuning = "coefvar",
  opt.params = list(tol = 1e-7, max.iter = 1000, alpha = 0.5, beta = 0.5)
)

# 2. ATE with standard errors
spbalATE(
  fit, df, outcome = "outcome",
  estimator = c("HT", "Hajek", "aIPTW"),
  se = c("iid", "hac", "boot"),
  mu.formula = outcome ~ treat + x1 + x2 + s(x, y, k = 300),   # for aIPTW
  boot = list(n.boot = 200, type = "grid")                   # block length chosen automatically
)
```

### Estimators

| `estimator` | Description |
|---|---|
| `"HT"` | Horvitz-Thompson IPTW |
| `"Hajek"` | Hajek (normalised) IPTW |
| `"aIPTW"` | augmented IPTW, with an outcome model given by `mu.formula` (fitted with `mgcv::gam`) |

### Standard errors

| `se` | Method |
|---|---|
| `"iid"` | influence-function standard error assuming independent observations |
| `"hac"` | spatial HAC (Conley) standard error of the influence functions, with bandwidth chosen from a fitted variogram (`HAC.se()`; options via `hac.args`) |
| `"boot"` | spatial block bootstrap: the propensity score (and outcome) model is refitted on each replicate, with resampled blocks moved to their new locations (via `spbootstrap`). Use `boot$type = "grid"` for complete gridded data, or `"point"` / `"polygon"` for point data in a rectangle or polygon study area. |

Lower-level functions are also exported: `ateEst()` (point estimates from
propensity scores), `HTInfluence()`, `HajekInfluence()` and `aIPTWInfluence()`
(influence functions), `HAC.se()`, `est.HAC()` and `calcBalance()`.

## Tuning parameter selection

Ideally, $\lambda$ should be selected in a data-driven fashion. The following selection approaches have been implemented.

-   `tuning = cv.score`: returns the $\lambda$ which minimizes the cross-validated balancing score,

$$ \lambda^{*} = \text{arg min}_{\lambda} \frac{1}{J} \sum_{j = 1}^{J} -S_n\big( \mathbf{Z}^{(j)} ; \hat{\boldsymbol{\alpha}}^{-(j)}_{\lambda} \big) $$

across $J$ CV folds. The number of CV folds, $J$, is specified using `folds`.

-   `tuning = cv.grad`: returns the $\lambda$ which minimizes the $L_p$ norm of the *gradient* of the balancing score across $J$ CV folds,

$$ \lambda^{*} = \text{arg min}_{\lambda} \frac{1}{J} \sum_{j = 1}^{J} \Vert \nabla S_n\big( \mathbf{Z}^{(j)} ; \hat{\boldsymbol{\alpha}}^{-(j)}_{\lambda} \big) \Vert_{p}. $$

Again, specify the number of CV folds using `folds` and the type of $L_p$ norm using `grad.norm` (options are $p = 1$, $2$, or $\infty$).

-   `tuning = coefvar`: returns the $\lambda$ with respect to the coefficient of variation of the balancing weights. In particular, define the coefficient of variation of the weights,

$$ \hat{w}_{\lambda} \equiv w(A_i, \hat{e}_{i, \lambda}) = \frac{A_i}{\hat{e}_i} + \frac{1 - A_i}{1 - \hat{e}_i},$$

as

$$ \text{CV}(\lambda) = \frac{sd(\hat{w}_{\lambda}) }{\text{mean}(\hat{w})}. $$

Choose the largest $\lambda$ such that coefficient of variation of its associated weights is greater than or equal to some specified proportion of the maximum coefficient of variation across all $\lambda$ values. In other words, choose

$$ \lambda^{*} = \text{max} \{ \lambda : CV(\lambda) \geq \rho CV_{max} \}, $$

where $CV_{max} = \text{max}_{\lambda} CV(\lambda)$ and $\rho \in (0,1)$ controls the desired coefficient of variation ratio. The choice of $\rho$ is can be specified with `coefvar.r`; the default is $\rho = 0.9$.

-   `tuning = max.bal`: returns the largest tuning parameter value such that the standardized difference in means for all model terms is less than some threshold. This threshold is specified using `bal.diff`.

-   `tuning = min`: returns the propensity score for the smallest $\lambda$ that was specified in `lambda`.

-   `tuning = all`: returns propensity score estimates using each of the previously mentioned methods. This is useful when comparing the performance of the selection method on simulated data or assessing the sensitivity of the ATE estimate under different selection methods.


## Repository layout

| Path | Contents |
|---|---|
| `R/`, `man/`, `tests/`, `DESCRIPTION`, `NAMESPACE` | the `spbalance` R package |
| `demo/` | `demo("spbalATE-sims", package = "spbalance")`: simulated examples with raster data, points in a box, and points in a polygon |
