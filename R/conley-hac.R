### conley-hac.R
### Nathan Wikle

### Compute the Conley HAC spatial variance quantity
###   V = sum_i sum_j K(d_ij / h) * e_hat_i * e_hat_j
### for Bartlett, Uniform, or Wendland kernel K.
###
### Two implementations are provided:
###   conleyHAC_grid -- for data on a regular Euclidean grid (uses 2D FFT,
###                     O(n log n) in the number of observations n)
###   conleyHAC_sf   -- for arbitrary point locations supplied as an sf object
###                     (uses STRtree spatial indexing, O(n log n + k) where k
###                     is the number of pairs within the bandwidth h)
###
### The wrapper conleyHAC dispatches to the appropriate method.

#--------------------------------------#
#--- 1. Kernel evaluation helper    ---#
#--------------------------------------#

.kernel_eval <- function(u, kernel) {
  # Evaluate spatial kernel K at u = d/h (elementwise).
  # Input
  #   u: non-negative numeric vector (distance / bandwidth)
  #   kernel: "bartlett", "uniform", or "wendland"
  # Output
  #   Numeric vector of kernel weights; zero wherever u > 1.
  #
  # Kernels:
  #   Bartlett  K(u) = max(0, 1 - u)
  #   Uniform   K(u) = 1(u <= 1)
  #   Wendland  K(u) = max(0, 1-u)^4 * (4u + 1)
  #     -- the Wendland C^2 function, positive definite in R^d for d <= 3.
  #        Unlike Bartlett/Uniform, it is guaranteed to produce a positive
  #        semi-definite weight matrix W in 2D, so V = e' W e >= 0 always.

  if (kernel == "bartlett") {
    pmax(0, 1 - u)
  } else if (kernel == "wendland") {
    w <- 1 - u
    w[w < 0] <- 0
    w^4 * (4 * u + 1)
  } else {
    as.numeric(u <= 1)
  }
}

#------------------------------------------------------#
#--- 2. Grid-based Conley HAC via 2D FFT            ---#
#------------------------------------------------------#

conleyHAC_grid <- function(e_hat, coords, h,
                           kernel = c("bartlett", "uniform", "wendland")) {
  # Compute V = sum_i sum_j K(d_ij/h) e_hat_i e_hat_j for data on a
  # regular Euclidean grid using 2D FFT-based autocorrelation.
  #
  # The double sum decomposes as
  #   V = sum_{delta_x, delta_y} K(dist(delta)/h) * C(delta_x, delta_y)
  # where C(delta_x, delta_y) = sum_{k,l} E[k,l] * E[k+delta_x, l+delta_y]
  # is the 2D sample autocovariance at lag (delta_x, delta_y).  All lags
  # C(*) are computed simultaneously in O(n log n) via the FFT power-spectrum
  # identity: C = IFFT(|FFT(E)|^2).
  #
  # Input
  #   e_hat  : numeric vector of residuals, length n = m_x * m_y
  #   coords : n x 2 numeric matrix (or data frame) of grid coordinates;
  #            columns are (x, y) and values must lie on a complete regular
  #            grid (uniform spacing in each dimension)
  #   h      : bandwidth (same units as coords)
  #   kernel : "bartlett" (default), "uniform", or "wendland"
  # Output
  #   Scalar sum_i sum_j K(d_ij/h) e_hat_i e_hat_j

  kernel <- match.arg(kernel)
  coords <- as.matrix(coords)
  n      <- length(e_hat)

  ### 1. detect grid structure
  x_vals <- sort(unique(coords[, 1]))
  y_vals <- sort(unique(coords[, 2]))
  m_x    <- length(x_vals)
  m_y    <- length(y_vals)

  if (m_x * m_y != n) {
    stop(paste0(
      "Grid dimensions (", m_x, " x ", m_y, " = ", m_x * m_y,
      ") do not match length(e_hat) = ", n, ". ",
      "Ensure coords defines a complete regular grid."
    ))
  }

  dx <- if (m_x > 1L) x_vals[2L] - x_vals[1L] else 1
  dy <- if (m_y > 1L) y_vals[2L] - y_vals[1L] else 1

  ### 2. fill grid matrix (robust to arbitrary input ordering)
  row_idx <- match(coords[, 1], x_vals)
  col_idx <- match(coords[, 2], y_vals)
  E <- matrix(0, m_x, m_y)
  E[cbind(row_idx, col_idx)] <- e_hat

  ### 3. zero-pad to avoid circular wrap-around; use next efficient FFT size
  n_pad_x <- nextn(2L * m_x - 1L, factors = c(2, 3, 5))
  n_pad_y <- nextn(2L * m_y - 1L, factors = c(2, 3, 5))
  E_pad <- matrix(0, n_pad_x, n_pad_y)
  E_pad[seq_len(m_x), seq_len(m_y)] <- E

  ### 4. 2D autocorrelation via FFT power spectrum
  # C[delta_x, delta_y] = sum_{k,l} E[k,l] * E[k + delta_x, l + delta_y]
  # Indexing convention after IFFT:
  #   lag delta >= 0  -->  row/col index  delta + 1
  #   lag delta <  0  -->  row/col index  n_pad + delta + 1
  fft_E  <- fft(E_pad)
  acorr  <- Re(fft(Mod(fft_E)^2, inverse = TRUE)) / (n_pad_x * n_pad_y)

  ### 5. enumerate lags within the bandwidth
  max_lag_x <- min(m_x - 1L, ceiling(h / dx))
  max_lag_y <- min(m_y - 1L, ceiling(h / dy))
  lag_grid  <- expand.grid(lx = seq.int(-max_lag_x, max_lag_x),
                           ly = seq.int(-max_lag_y, max_lag_y))

  lag_dist <- sqrt((lag_grid$lx * dx)^2 + (lag_grid$ly * dy)^2)
  k_vals   <- .kernel_eval(lag_dist / h, kernel)

  ### 6. discard zero-weight lags, then look up and accumulate
  active   <- k_vals > 0
  lag_grid <- lag_grid[active, , drop = FALSE]
  k_vals   <- k_vals[active]

  ix <- ifelse(lag_grid$lx >= 0L,
               lag_grid$lx + 1L,
               n_pad_x + lag_grid$lx + 1L)
  iy <- ifelse(lag_grid$ly >= 0L,
               lag_grid$ly + 1L,
               n_pad_y + lag_grid$ly + 1L)

  sum(k_vals * acorr[cbind(ix, iy)])
}

#------------------------------------------------------#
#--- 3. sf-based Conley HAC via spatial indexing    ---#
#------------------------------------------------------#

conleyHAC_sf <- function(e_hat, sf_obj, h,
                         kernel = c("bartlett", "uniform", "wendland")) {
  # Compute V = sum_i sum_j K(d_ij/h) e_hat_i e_hat_j for point data
  # stored as an sf object.
  #
  # Algorithm:
  #   1. Find all pairs (i, j) with d_ij <= h using sf's STRtree spatial
  #      index (sf::st_is_within_distance), O(n log n + k).
  #   2. Compute Euclidean distances from raw coordinates for those pairs.
  #   3. Apply kernel weights and accumulate the sum (fully vectorised).
  #
  # For large n with small h (sparse pairs), this is far more efficient than
  # the naive O(n^2) approach.  For moderate n, a dense O(n^2) alternative
  # using e_hat %*% K %*% e_hat (where K is the full kernel matrix) may be
  # faster due to BLAS.
  #
  # Input
  #   e_hat  : numeric vector of residuals, length n
  #   sf_obj : sf object with POINT geometry, length n; should be in a
  #            projected (planar) CRS so that Euclidean distances are valid
  #   h      : bandwidth in the CRS units of sf_obj
  #   kernel : "bartlett" (default), "uniform", or "wendland"
  # Output
  #   Scalar sum_i sum_j K(d_ij/h) e_hat_i e_hat_j

  kernel <- match.arg(kernel)
  n      <- nrow(sf_obj)

  if (length(e_hat) != n) {
    stop("length(e_hat) must equal the number of features in sf_obj.")
  }

  # warn if geographic (lon/lat) CRS: distances below are Euclidean on raw
  # degree coordinates, not geodesic
  if (!is.na(sf::st_crs(sf_obj)) && isTRUE(sf::st_is_longlat(sf_obj))) {
    warning(
      "sf_obj has a geographic (lon/lat) CRS. Distances are computed as ",
      "Euclidean on raw degree coordinates, not geodesic. ",
      "Transform to a projected CRS (e.g. sf::st_transform) for ",
      "physically meaningful distances."
    )
  }

  ### 1. spatial-index range query: all (i, j) with d_ij <= h
  nb    <- sf::st_is_within_distance(sf_obj, sf_obj, dist = h, sparse = TRUE)
  i_idx <- rep(seq_len(n), lengths(nb))
  j_idx <- unlist(nb, use.names = FALSE)

  if (length(i_idx) == 0L) return(0)

  ### 2. Euclidean distances from raw coordinates (avoids sf overhead per pair)
  coords  <- sf::st_coordinates(sf_obj)[, 1:2, drop = FALSE]
  dx_ij   <- coords[i_idx, 1L] - coords[j_idx, 1L]
  dy_ij   <- coords[i_idx, 2L] - coords[j_idx, 2L]
  d_ij    <- sqrt(dx_ij^2 + dy_ij^2)

  ### 3. kernel weights and weighted sum
  k_vals  <- .kernel_eval(d_ij / h, kernel)
  sum(k_vals * e_hat[i_idx] * e_hat[j_idx])
}

#------------------------------------------------------#
#--- 4. Dispatcher wrapper                          ---#
#------------------------------------------------------#

conleyHAC <- function(e_hat, coords, h,
                      kernel = c("bartlett", "uniform", "wendland"),
                      method = NULL) {
  # Compute the Conley HAC spatial variance quantity
  #   V = sum_i sum_j K(d_ij / h) * e_hat_i * e_hat_j
  # dispatching to the grid (FFT) or sf (spatial-index) implementation.
  #
  # Input
  #   e_hat  : numeric vector of residuals, length n
  #   coords : n x 2 coordinate matrix for method = "grid" (evenly-spaced
  #            regular grid), or an sf object for method = "sf"
  #   h      : bandwidth
  #   kernel : "bartlett" (default), "uniform", or "wendland"
  #   method : "grid" or "sf"; if NULL (default), detected automatically from
  #            the class of coords
  # Output
  #   Scalar sum_i sum_j K(d_ij/h) e_hat_i e_hat_j

  kernel <- match.arg(kernel)

  if (is.null(method)) {
    method <- if (inherits(coords, "sf")) "sf" else "grid"
  }

  if (method == "sf") {
    conleyHAC_sf(e_hat, sf_obj = coords, h = h, kernel = kernel)
  } else {
    conleyHAC_grid(e_hat, coords = coords, h = h, kernel = kernel)
  }
}
