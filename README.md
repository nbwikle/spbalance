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

**Important:** Several R scripts (`spbal-tprs-demo.R`, `spbal-krr-demo.R`) there exists a local Python virtual environment named `r-reticulate`. To use these functions, either create this virtual environment or modify instances of the following code to match your preferred virtual environment.

```{r}
# indicate that we want to use a specific Python virtualenv
use_virtualenv("r-reticulate")
```

### Covariance matrices

The Tec et al. simulation study repeatedly simulates from a Gaussian process (GP), which becomes computationally intensive as the number of data grow. To limit computational overhead, I've created an R script which saves a spectral decomposition of three common GP covariance matrices (exponential, squared exponential (aka, Gaussian), and Matern), with spatial dependence decaying at three different rates, evaluated on regular grids of size $10\times10$, $25\times25$, $50\times50$, and $100\times100$. The matrices are saved in `./data/cov-mats/` and can be quickly loaded into R when performing the Tec et al. simulation.

To create the matrices, run the following R code **one time**.

```{r}
# create and save SVD of GP covariance matrices
library("here")
source(here::here("R", "dgp", "cov-mat-creation.R"))
```

**Important:** creating and saving all the matrices requires \~6.5 GB of space and takes \~11 hours on my laptop.

## Balancing Weights with Spatial Confounding

### Background

The primary focus of this work is to develop inverse probability of treatment weighting (IPTW) estimators that can be used in settings with unmeasured spatial confounding. Let $(Y_s, A_s, X_s)$ denote data observed at location $s \in \mathcal{D}$, where $Y_s$ is the outcome variable, $A_s$ is a binary treatment variable, and $X_s$ is a vector of observed confounders. Similarly, let $Y_s(1)$ and $Y_s(0)$ denote the *potential outcomes* that would have been observed at location $s$ had that unit been assigned treatment or control, respectively. We will make standard identifying assumptions, SUTVA and positivity, however, we assume that exchangeability holds conditional on the observed confounders, $X_s$, as well as an unobserved confounder, $U_s$. That is,$$ \{Y(1), Y(0)\} \perp\!\!\!\!\perp A \; | \; X, U $$Because $U_s$ is unobserved, estimating a causal effect in such settings is challenging. However, if we are willing to assume that $U_s$ has some additional structure, in particular, that it can be represented as a smooth function of space, we may have some hope in estimating an average causal effect.

In particular, we will focus on estimating the average treatment effect (ATE), $\tau = E\big[Y(1) - Y(0)\big]$, using an IPTW estimator. Let $e_s(x) = Pr(A_s = 1 | X_s = x)$ denote the propensity score at location $s \in \mathcal{D}$. The standard IPTW estimator is given as

$$
\hat{\tau}_{iptw} = \frac{1}{n} \sum_{s = 1}^{n} \bigg( \frac{A_s}{\hat{e}_s}Y_s - \frac{1 - A_s}{1 - \hat{e}_s} Y_s \bigg),
$$

where $\hat{e}_s$ is the estimated propensity score at location $s \in \mathcal{D}$. If we are willing to assume that $U_s$ is a smooth function of space, we might consider estimating $\hat{e}_s$ with a model that explicitly accounts for spatial dependence. For example, we could model $\hat{e}_s$ using a generalized linear mixed model with a spatial random effect, or as a generalized additive model that includes a spatial smooth. Interestingly, such an approach does not always result in well-behaved estimation (see, e.g., Tec et al. (2023)), even when a reasonable spatial propensity score model is used.

Alternatively, we might consider estimating inverse propensity score weights that *explicitly* balances both the observed covariates, $X_s$, as well as the unobserved confounder, $U_s$. To do so, we assume that $U_s$ can be well-approximated by a spatial function, $f(s)$, where $f$ is assumed to belong to some class of spatially smooth smooth functions, $\mathcal{F}$. A variety of function spaces might be considered, including Gaussian processes and reproducing kernel Hilbert spaces (RKHS). Importantly, in many of these cases, $f(s)$ can be represented as a linear combination of *spatial basis functions*, $\{\phi_j(\cdot) : j = 1, 2, \dots \}$, where $f \in \mathcal{F}$ implies that $f(s) = \sum_{j \geq 1} a_j \phi_j(s)$. Each basis function $\phi_j(\cdot)$ generates a particular spatially-varying feature describing some component of the underlying variation in $f(s)$, and $f(s)$ is simply a linear combination of these features. Consequently, if we can find IPT weights that empirically balance the observed confounders, $X_s$, *and* the basis functions, $\phi_j(\cdot)$, then we will have empirically balanced both $X_s$ and $U_s$.

We can accomplish this task by fitting $\hat{e}_s$ using the tailored loss function of Zhao (2019). To do so, let $\tilde{\mathbf{x}}_s = (1, X_s, \phi_1(s), \dots, \phi_J(s))$ denote the collection of observed confounders and spatial basis functions evaluated at location $s \in \mathcal{D}$. The propensity score model, $e_s(\boldsymbol{\alpha})$, is then assumed to have the following functional form, $\text{logit}(e_s) = \tilde{\mathbf{x}}_s' \boldsymbol{\alpha}$, and the model coefficients are estimated using the following objective function: $$ \mathcal{M}(\boldsymbol{\alpha}) = S_n(\boldsymbol{\alpha}) - \lambda J_{\boldsymbol{\alpha}}(e),$$ where $$ S_n(\boldsymbol{\alpha}) =  \bigg(\log \frac{e_s(\boldsymbol{\alpha})}{1 - e_s(\boldsymbol{\alpha})} - \frac{1}{e_s(\boldsymbol{\alpha})} \bigg) A_s + \bigg( \log \frac{1 - e_s(\boldsymbol{\alpha})}{e_s(\boldsymbol{\alpha})} - \frac{1}{1 - e_s(\boldsymbol{\alpha})} \bigg) (1 - A_s) $$is the covariate balancing tailored loss proposed by Zhao (2019), $J_{\boldsymbol{\alpha}}(e)$ is a penalty on the function class of $e(\boldsymbol{\alpha})$, and $\lambda \geq 0$ is a tuning parameter controlling the degree of penalization. Notably, it can be shown that when $\lambda = 0$, the weighted mean difference of $\tilde{\mathbf{x}}_s$ between treatment and control locations is exactly zero, implying that the any function in the span of $\tilde{\mathbf{x}}_s$ has been empirically balanced---including $f(s) = \sum_{j \geq 1} a_j \phi_j(s)$. To prevent overfitting, $\lambda$ and $J_{\boldsymbol{\alpha}}(e)$ can be chosen to control the degree of smoothness of $\text{logit}(e_s)$.

### Using `spBalance`

Spatial IPTW balancing weights can be estimated using the `spBalance` function (the source code for this function is largely contained in `spbalance.R`). The function contains several important arguments, including

1. `formula`

  - A formula term relating treatment variable, $a$, to observed covariates, $x$, and some spatial effect, $f(s)$. At the moment, spatial effects can be specified using terms as a spatial smooth term, `s(...)`, from the [`mgcv` package](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/s.html), or as a kernel ridge regression term, `krr(...)`.
  
2. `lambda`

  - A vector of possible tuning parameter values that should be considered during the model fit.
  
3. `tuning`

  - The 

```{r}
spBalance <- function(
  formula, data, lambda, init.params = NULL, fit.gam = FALSE, pen.int = FALSE,
  tuning = "none", folds = 10, coefvar.r = 0.9, bal.diff = 0.1, hide.details = TRUE,
  opt.params = list(tol = 1e-7, max.iter = 100, alpha = 0.5, beta = 0.5
){
  ...
  
  'krr(
  #       x1, ..., xp, # variables to include in kernel evaluation
  #       kern = "SE", # type of kernel; must be one of 'SE' or 'Exp'
  #       kp = p # number of variables included before "kern" input,
  #       centered = FALSE # boolean indicating if the KRR term should be centered at zero
  #     )'
}
```

### Tuning parameter selection

## Simulation Studies

All data have been downloaded from publicly available databases. To aid in reproducability efforts, the data have been archived for download on Zenodo (doi = ); they require X GB of space within the local directory.

**Intervention Data**
