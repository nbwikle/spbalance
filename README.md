# Balancing weights with unobserved spatial confounding

*Author: Nathan Wikle*

## Summary

This repository contains the functions and source code necessary to reproduce the main results in the paper, "Balancing weights in settings with unobserved spatial confounding" (WORKING TITLE). A short introduction to the `spBalance` function and a corresponding simulation study are described below.

## Getting Started

### Packages

The code has been tested with R version 4.4.1, "Race for Your Life." The following R packages must be installed before the code will run successfully (package version used during testing are shown in parantheses):

-   `geoR` (1.9.4)
-   `here` (1.0.1)
-   `MASS` (7.3.61)
-   `Matrix` (1.7.0)
-   `matrixStats` (1.5.0)
-   `mgcv` (1.9.1)
-   `raster` (3.6.30)
-   `RColorBrewer` (1.1.3)
-   `reticulate` (1.37.0)
-   `sf` (1.0.16)
-   `spmodel` (0.8.0)
-   `terra` (1.8.54)

### `Reticulate` Package

The Tec et al. (2023) simulation study uses functions from the `scipy.signal` module to perform image convolutions. Performing the simulation study in R requires the installation of the `reticulate` package; see [here](https://rstudio.github.io/reticulate/) for installation instructions.

**Important:** Several R scripts (`spbal-tprs-demo.R`, `spbal-krr-demo.R`) require the existence of a local Python virtual environment named `r-reticulate`. To use these functions, either create this virtual environment or modify instances of the following code to match your preferred virtual environment.

``` r
# indicate that we want to use a specific Python virtualenv
use_virtualenv("r-reticulate")
```

### Covariance matrices

The Tec et al. simulation study repeatedly simulates from a Gaussian process (GP), which becomes computationally intensive as the number of data grow. To limit computational overhead, I've created an R script which saves a spectral decomposition of three common GP covariance matrices (exponential, squared exponential (aka, Gaussian), and Matern), with spatial dependence decaying at three different rates, evaluated on regular grids of size $10\times10$, $25\times25$, $50\times50$, and $100\times100$. The matrices are saved in `./data/cov-mats/` and can be quickly loaded into R when performing the Tec et al. simulation.

To create the matrices, run the following R code **one time**.

``` r
# create and save SVD of GP covariance matrices
library("here")
source(here::here("R", "dgp", "cov-mat-creation.R"))
```

**Important:** creating and saving all the matrices requires \~6.5 GB of space and takes \~11 hours on my laptop.

## Balancing Weights with Spatial Confounding

### Background

The primary focus of this work is to develop inverse probability of treatment weighting (IPTW) estimators that can be used in settings with unmeasured spatial confounding. Let $(Y_s, A_s, X_s)$ denote data observed at location $s \in \mathcal{D}$, where $Y_s$ is the outcome variable, $A_s$ is a binary treatment variable, and $X_s$ is a vector of observed confounders. Similarly, let $Y_s(1)$ and $Y_s(0)$ denote the *potential outcomes* that would have been observed at location $s$ had that unit been assigned treatment or control, respectively. We will make standard identifying assumptions, SUTVA and positivity, however, we assume that exchangeability holds conditional on the observed confounders, $X_s$, as well as an unobserved confounder, $U_s$. That is,

$$
\{Y(1), Y(0)\} \perp A | X, U.
$$

Because $U_s$ is unobserved, estimating a causal effect in such settings is challenging. However, if we are willing to assume that $U_s$ has some additional structure, in particular, that it can be represented as a smooth function of space, we may have some hope in estimating an average causal effect.

In particular, we will focus on estimating the average treatment effect (ATE), $\tau = E\big[Y(1) - Y(0)\big]$, using an IPTW estimator. Let $e_s(x) = Pr(A_s = 1 | X_s = x)$ denote the propensity score at location $s \in \mathcal{D}$. The standard IPTW estimator is given as

$$
\hat{\tau}_{iptw} = \frac{1}{n} \sum_{s = 1}^{n} \bigg( \frac{A_s}{\hat{e}_s}Y_s - \frac{1 - A_s}{1 - \hat{e}_s} Y_s \bigg),
$$

where $\hat{e}_s$ is the estimated propensity score at location $s \in \mathcal{D}$. If we are willing to assume that $U_s$ is a smooth function of space, we might consider estimating $\hat{e}_s$ with a model that explicitly accounts for spatial dependence. For example, we could model $\hat{e}_s$ using a generalized linear mixed model with a spatial random effect, or as a generalized additive model that includes a spatial smooth. Interestingly, such an approach does not always result in well-behaved estimation (see, e.g., Tec et al. (2023)), even when a reasonable spatial propensity score model is used.

Alternatively, we might consider estimating inverse propensity score weights that *explicitly* balances both the observed covariates, $X_s$, as well as the unobserved confounder, $U_s$. To do so, we assume that $U_s$ can be well-approximated by a spatial function, $f(s)$, where $f$ is assumed to belong to some class of spatially smooth smooth functions, $\mathcal{F}$. A variety of function spaces might be considered, including Gaussian processes and reproducing kernel Hilbert spaces (RKHS). Importantly, in many of these cases, $f(s)$ can be represented as a linear combination of *spatial basis functions*, $\{\phi_j(\cdot) : j = 1, 2, \dots \}$, where $f \in \mathcal{F}$ implies that $f(s) = \sum_{j \geq 1} a_j \phi_j(s)$. Each basis function $\phi_j(\cdot)$ generates a particular spatially-varying feature describing some component of the underlying variation in $f(s)$, and $f(s)$ is simply a linear combination of these features. Consequently, if we can find IPT weights that empirically balance the observed confounders, $X_s$, *and* the basis functions, $\phi_j(\cdot)$, then we will have empirically balanced both $X_s$ and $U_s$.

We can accomplish this task by fitting $\hat{e}_s$ using the tailored loss function of Zhao (2019). To do so, let $\tilde{\mathbf{x}}_s = (1, X_s, \phi_1(s), \dots, \phi_J(s))$ denote the collection of observed confounders and spatial basis functions evaluated at location $s \in \mathcal{D}$. The propensity score model, $e_s(\boldsymbol{\alpha})$, is then assumed to have the following functional form, $\text{logit}(e_s) = \tilde{\mathbf{x}}_s' \boldsymbol{\alpha}$, and the model coefficients are estimated using the following objective function:

$$ \mathcal{M}(\boldsymbol{\alpha}) = S_n(\boldsymbol{\alpha}) - \lambda J_{\boldsymbol{\alpha}}(e),$$

where

$$ S_n(\boldsymbol{\alpha}) =  \frac{1}{n} \sum_{i =1 }^{n}\bigg(\log \frac{e_s(\boldsymbol{\alpha})}{1 - e_s(\boldsymbol{\alpha})} - \frac{1}{e_s(\boldsymbol{\alpha})} \bigg) A_s + \bigg( \log \frac{1 - e_s(\boldsymbol{\alpha})}{e_s(\boldsymbol{\alpha})} - \frac{1}{1 - e_s(\boldsymbol{\alpha})} \bigg) (1 - A_s) $$

is the covariate balancing tailored loss proposed by Zhao (2019), $J_{\boldsymbol{\alpha}}(e)$ is a penalty on the function class of $e(\boldsymbol{\alpha})$, and $\lambda \geq 0$ is a tuning parameter controlling the degree of penalization. Notably, it can be shown that when $\lambda = 0$, the weighted difference in means of $\tilde{\mathbf{x}}_s$ between treatment and control locations is exactly zero, implying that the *any* function in the span of $\tilde{\mathbf{X}}$ has been empirically balanced, including $f(s)$. To prevent overfitting, $\lambda$ can be chosen to control the degree of smoothness of $\text{logit}(e_s)$.

### Using `spBalance`

Spatial balancing weights can be estimated using the `spBalance` function (the source code for this function is largely contained in `./R/spbalance.R`). The function contains several important arguments, including

1.  `formula`

-   A formula term relating treatment variable, $a$, to observed covariates, $x$, and some spatial effect, $f(s)$. At the moment, spatial effects can be specified using terms as a spatial smooth term, `s(...)`, from the [`mgcv` package](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/s.html), or as a kernel ridge regression term, `krr(...)`.

2.  `lambda`

-   A vector of possible tuning parameter values that should be considered during the model fit.

3.  `tuning`

-   The method used to select the tuning parameter, $\lambda$, from the list of plausible values specified with `lambda`. Options include `cv.score`, `cv.grad`, `coefvar`, `max.bal`, `min`, and `all`. A description of each method is included below.

Two code demos illustrate the use of `spBalance` on simulated data sets. The first demo, `R/demo/spbal-tprs-demo.R`, simulates data from the Tec et al. (2023) setting and compare average treatment effect estimates using weights where the propensity score is modeled as (i) a generalized additive model with a thin-plate regression spline (TPRS) to account for the unmeasured spatial confounder, and (ii) a covariate balancing propensity score estimate in which balance is obtained with respect to the TPRS basis functions. The key function call is shown below. Notice that the TPRS is specified with `treat ~ s(x, y, k = 500)`.

``` r
# select tuning parameter using coefficient of variation method
spb.cvr <- spBalance(
  formula = treat ~ s(x, y, k = 500),
  data = sim.df,
  lambda = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 5, 10),
  tuning = "coefvar",
  hide.details = TRUE,
  opt.params = list(tol = 1e-7, max.iter = 1000, alpha = 0.5, beta = 0.5)
)

# chosen lambda:
spb.cvr$lambda

# ATE estimate (true value = 0.1)
ateEst(out = sim.df$outcome, trt = sim.df$treat, pr.trt = spb.cvr$pi.hat)$ate
}
```

The second demo, `R/demo/spbal-krr-demo.R`, simulates data similar to the procedure in Zhao (2019). Here, balance is obtained with respect to an RKHS kernel basis function, $k(\mathbf{x}_i, \mathbf{x}_j)$. The demo code illustrates how this can be accomplished using `spBalance`. In particular, `trt ~ krr(x, kern = "SE", kp = 5, centered = FALSE)` specifies that we want to model $f(x)$ using a kernel ridge regression with a squared exponential (i.e., gaussian) kernel function, with lengh-scale (i.e., range) parameter set to $5$; the function is not centered at 0, so we do not need an intercept.

``` r
# Fit model
fit.krr <- spBalance(
  formula = trt ~ krr(x1, x2, x3, x4, x5, kern = "SE", kp = 5, centered = FALSE),
  data = sim.k,
  lambda = c(0.1, 0.25, 0.5, 0.75, 1, 2, 5, 10),
  tuning = "coefvar", coefvar.r = 0.5,
  hide.details = TRUE,
  opt.params = list(tol = 1e-7, max.iter = 1000, alpha = 0.5, beta = 0.5)
)
fit.krr$lambda
ateEst(out = zhao.sim$out, trt = zhao.sim$trt, pr.trt = fit.krr$pi.hat)$ate
```

The KRR optimization is reasonably fast when the data size $n$ is small, however, it becomes computationaly prohibitive as $n$ grows large. To contend with this issue, `spBalance` includes a low-rank approximation of the KRR kernel using random Fourier features (RFF, see Rahimi and Recht (2007)). The `krr` term is the same as before, but with it now includes two additional inputs: `rff = TRUE`, indicating that the RFF approximation should be used, and `nf`, which controls the number of random features used in the approximation.

``` r
fit.rff <- spBalance(
  formula = trt ~ krr(x1, x2, x3, x4, x5, kern = "SE", kp = 5, 
                      centered = FALSE, rff = TRUE, nf = 250),
  data = sim.k,
  lambda = c(0.1, 0.25, 0.5, 0.75, 1, 2, 5, 10),
  tuning = "coefvar", coefvar.r = 0.5,
  hide.details = TRUE,
  opt.params = list(tol = 1e-7, max.iter = 1000, alpha = 0.5, beta = 0.5)
)
fit.rff$lambda
ateEst(out = zhao.sim$out, trt = zhao.sim$trt, pr.trt = fit.rff$pi.hat)$ate
```

### Tuning parameter selection

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

The implementation of these methods can be found in `R/tuning-params.R`.

## Simulation Studies

### Data Generating Process 1 (Tec et al. (2023))

We considered the performance of the spatial balancing weight estimators using data simulated according to the data generating process described in Tec et al. (2023). A full description of the simulation steps can be found in Tec et al. (2023). For convenience, we provide a brief overview of the simulation procedure.

**Overview**

-   Simulate $L_s \sim GP(0, k^{(1)}_{\theta_1})$, where $k^{(1)}_{\theta_1}$ is a squared exponential covariance function with range parameter $\theta_1 = 5.1$:

$$ k^{(1)}_{\theta_1}(d) = \exp\bigg( \frac{-d^2}{2 \theta_1^2} \bigg).$$

-   Let $(X_{1}, X_{2})_s = \nabla L_s$. These are meant to represent observed "wind vector fields".

-   Define $\boldsymbol{\mu} = K_1 * \mathbf{X}_1 + K_2 * \mathbf{X}_2$, where $K_j$, $j = 1,2$ are specially designed convolution kernels. In other words, $\mu_s$ is a function of both the locally observed $X_s$ as well as a smoothed version of the covariates in the region surrounding $s$.

-   Simulate treatment assignment from $A_s \sim \text{Bernoulli}(p_s)$, where $p_s = \text{expit}(\mu_s)$.

-   Simulate outcome observations as

$$ y_s = \tau A_s - \sqrt{0.5} \mu_s + \eta_s + \epsilon_s, $$

where

$$ \eta_s \sim GP(0, \frac{1}{4} k^{(2)}_{\theta_2}), \text{ } \epsilon_s \sim N(0, \sigma^2 = \frac{1}{4}), \text{ and } \tau = 0.1.$$

In particular, the error consists of a spatial random effect, $\eta_s$, and a measurement error term, $\epsilon_s$. The covariance function of the spatial random effect is user-specified, and can be one of squared exponential, exponential, or Matern:

$$ k^{\text{SE}}_{\theta_2}(d) = \exp\bigg( \frac{-d^2}{2 \theta_2^2} \bigg) $$

$$ k^{\text{Exp}}_{\theta_2}(d) = \exp\bigg( \frac{-d}{\theta_2} \bigg) $$

$$ k^{\text{Matern}}_{\theta_2}(d) = \frac{2^{1 - \nu}}{\Gamma(\nu)} \bigg( \sqrt{8\nu} \frac{d}{\rho} \bigg)^{\nu} K_{\nu}\bigg( \sqrt{8\nu} \frac{d}{\rho} \bigg), $$

where $K_{\nu}$ is the modified Bessel function of the second kind (order $\nu$), and $\theta_2 = (\nu, \rho)$ consists of the smoothness and range parameters, respectively.

-   Finally, a variation of the above procedure considers a nonlinear transformation of $\mathbf{X}_1$, $\mathbf{X_2}$:

$$\boldsymbol{\mu} = K_1 * \text{sign}(\mathbf{X}_1) + K_2 * \text{sign}(\mathbf{X}_2) $$

**Implementation**

A simulation study was performed to compare the performance of several estimators of the ATE on data generated from DGP1, including a naive difference-in-means estimator, IPTW estimators using propensity scores fitted with via maximum likelihood estimation or with the balancing score objective function, an outcome regression estimator of the ATE, and augmented IPTW (aIPTW) estimators that efficiently combine the outcome regression and weighting estimators.

The simulation study is implemented in `R/dgp/dgp1-sim-study.R`. The R script can be called via terminal as follows (note: this assumes that your current working directory is set to the repo root directory):

``` console
Rscript ./R/dgp/dgp1-sim-study.R arg1 arg2 arg3 arg4 arg5 arg6 arg7
```

There are several command line arguments that must be specified:

-   `arg1`: type of outcome covariance function; must be one of `SE`, `Matern`, or `Expo`

-   `arg2`: covariance function parameterization; must be one of `1`, `2`, or `3`. These correspond to covariance functions where $k(d) \approx 0.1$ at distances $d = 3$, $5$, and $10$.

-   `arg3`: the $n \times n$ grid size; must be one of `10`, `25`, `50`, or `100`

-   `arg4`: where to seed the random number generator

-   `arg5`: the number of simulations to perform

-   `arg6`: whether the simulated $\boldsymbol{\mu}$ is nonlinear (`T`) or not (`F`)

-   `arg7`: an option argument specifying the number of random samples to keep from the $n \times n$ grid. For example, if this is 1000, then only 1000 random samples are kept from the $n \times n$ grid. Tec et al. refer to this as the "sparse" setting. **Do not include** this argument if you want to keep data for the entire $n \times n$ grid..

As an example,

``` console
Rscript ./R/dgp/dgp1-sim-study.R SE 3 100 0 25 F
```

would perform a simulation study where data are simulated from DGP1 on a $100 \times 100$ grid (resolution $1 \times 1$), with the outcome dependence simulated as a GP with a squared exponential covariance that decays to a correlation of $0.1$ when the distance between two units is 10. The random number generator would be set to `0 + sim.i`, where `sim.i` is the $i$th simulation. The number of simulations performed is 25; for each simulation, the ATE is estimated using a variety of possible estimators. Finally, we simulate from the setting with a linear $\boldsymbol{\mu}$.

The **results** from the simulation study are saved as an `.RDS` file in `/output/dgp1/`. The file size will depend on the number of simulations that are performed, but typically the size is small ($<0.5$ MB). Implementing the simulation study on the $100 \times 100$ grid can take a considerable amount of time. It is recommended you parallelize this task on an HPC, if available.
