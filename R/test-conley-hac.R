### test-conley-hac.R
### Nathan Wikle

### Correctness tests for conleyHAC_grid, conleyHAC_sf, and conleyHAC.
### Run with: testConleyHAC()

testConleyHAC <- function(tol = 1e-8, verbose = TRUE) {
  # Test suite for the Conley HAC variance estimator.
  # Compares FFT and sf outputs against a brute-force O(n^2) reference,
  # verifies known analytical results, and checks structural properties.
  # Input
  #   tol    : absolute tolerance for floating-point comparisons (default 1e-8)
  #   verbose: if TRUE, prints each test result (default TRUE)
  force(tol)   # evaluate the promise now so .near() can close over it
  # Output
  #   Invisibly returns a list:
  #     results: named logical vector (TRUE = pass, FALSE = fail)
  #     n_pass : integer, number of tests passed
  #     n_total: integer, total tests run

  results <- logical(0)

  .record <- function(passed, label, details = NULL) {
    results[[label]] <<- isTRUE(passed)
    if (verbose) {
      tag <- if (isTRUE(passed)) "PASS" else "FAIL"
      cat(sprintf("  [%s] %s\n", tag, label))
      if (!isTRUE(passed) && !is.null(details))
        cat(sprintf("       %s\n", details))
    }
  }

  .near <- function(a, b) abs(a - b) < tol

  # Brute-force reference: O(n^2) direct double sum
  .brute <- function(e, coords, h, kernel) {
    D <- as.matrix(dist(coords))
    u <- D / h
    # Note: pmax() strips the dim attribute from a matrix, so we restore it
    # afterwards. Using k * (k > 0) in place of pmax(0, k) is an alternative
    # that preserves dims, but the explicit restoration is clearer.
    K <- switch(kernel,
      bartlett = pmax(0, 1 - u),
      uniform  = (u <= 1) * 1.0,
      wendland = { w <- pmax(0, 1 - u); w^4 * (4 * u + 1) }
    )
    dim(K) <- dim(u)
    drop(e %*% K %*% e)
  }

  set.seed(42)
  kernels <- c("bartlett", "uniform", "wendland")

  #-----------------------------------------------------------------#
  #--- 1. Grid FFT vs brute force (5x5 grid, h = 2.5)           ---#
  #-----------------------------------------------------------------#

  if (verbose) cat("--- Grid FFT vs brute force (5x5, h = 2.5) ---\n")

  m5      <- 5L
  coords5 <- as.matrix(expand.grid(x = seq_len(m5), y = seq_len(m5)))
  e5      <- rnorm(m5 * m5)
  h5      <- 2.5

  for (kern in kernels) {
    bf  <- .brute(e5, coords5, h5, kern)
    fft <- conleyHAC_grid(e5, coords5, h5, kern)
    .record(.near(fft, bf),
            paste0("grid FFT vs brute force [", kern, "]"),
            sprintf("FFT = %.10g, BF = %.10g, diff = %.2e", fft, bf, abs(fft - bf)))
  }

  #-----------------------------------------------------------------#
  #--- 2. sf vs brute force (same data as above)                 ---#
  #-----------------------------------------------------------------#

  if (requireNamespace("sf", quietly = TRUE)) {
    if (verbose) cat("--- sf vs brute force (5x5, h = 2.5) ---\n")
    pts5 <- sf::st_as_sf(as.data.frame(coords5),
                         coords = c("x", "y"), crs = NA_crs_)
    for (kern in kernels) {
      bf   <- .brute(e5, coords5, h5, kern)
      sf_v <- suppressWarnings(conleyHAC_sf(e5, pts5, h5, kern))
      .record(.near(sf_v, bf),
              paste0("sf vs brute force [", kern, "]"),
              sprintf("sf = %.10g, BF = %.10g, diff = %.2e", sf_v, bf, abs(sf_v - bf)))
    }
  } else {
    if (verbose) cat("  [SKIP] sf tests (package not available)\n")
  }

  #-----------------------------------------------------------------#
  #--- 3. Analytical results on a 5x5 unit grid                 ---#
  #                                                                 #
  # Residuals: e_i = 1 for all i.  Bandwidth: h = 1.5.             #
  # With constant residuals V = sum_{i,j} K(d_{ij}/h), so the      #
  # exact answer depends only on how many ordered pairs fall in     #
  # each distance class.                                            #
  #                                                                 #
  # Pair counts on a 5x5 unit grid for d <= 1.5:                   #
  #   d = 0:       25  (self-pairs, i = j)                         #
  #   d = 1:       80  (horizontal/vertical adjacent,              #
  #                     4 * 4 + 4 * (4+1) + ... = 2*(4*5+5*4) )   #
  #   d = sqrt(2): 64  (diagonal adjacent, 4*(m-1)^2 = 4*16)       #
  #                     sqrt(2) ~ 1.414 < 1.5, so included         #
  #   d = 2:       60  (2-step h/v), 2 > 1.5, not included         #
  #                                                                 #
  # Expected: V = 25*K(0) + 80*K(1/1.5) + 64*K(sqrt(2)/1.5)       #
  #   Uniform:  25 + 80 + 64         = 169 (exact integer)         #
  #   Bartlett: 25 + 80*(1/3) + 64*(1 - sqrt(2)/1.5)              #
  #   Wendland: 25 + 80*K(2/3) + 64*K(sqrt(2)*2/3)                #
  #-----------------------------------------------------------------#

  if (verbose) cat("--- Analytical results on 5x5 grid (e_i=1, h=1.5) ---\n")

  h_a    <- 1.5
  e_ones <- rep(1, m5 * m5)

  # Pair counts and distances derived from grid geometry
  pair_dists  <- c(0, 1, sqrt(2))
  pair_counts <- c(25L, 80L, 64L)

  for (kern in kernels) {
    exp_v <- sum(pair_counts * .kernel_eval(pair_dists / h_a, kern))
    V     <- conleyHAC_grid(e_ones, coords5, h_a, kern)
    .record(.near(V, exp_v),
            paste0("5x5 grid analytical [", kern, "]"),
            sprintf("V = %.10g, expected = %.10g, diff = %.2e",
                    V, exp_v, abs(V - exp_v)))
  }

  #-----------------------------------------------------------------#
  #--- 4. Small h -> sum(e^2)                                    ---#
  #                                                                 #
  # When h < grid spacing, no off-diagonal pairs fall within the   #
  # bandwidth, so only the i = j terms survive: V = sum(e^2).      #
  #-----------------------------------------------------------------#

  if (verbose) cat("--- Small h -> sum(e^2) (h = 0.1, grid spacing = 1) ---\n")

  h_small    <- 0.1
  expected_s <- sum(e5^2)
  for (kern in kernels) {
    V <- conleyHAC_grid(e5, coords5, h_small, kern)
    .record(.near(V, expected_s),
            paste0("small h -> sum(e^2) [", kern, "]"),
            sprintf("V = %.10g, sum(e^2) = %.10g", V, expected_s))
  }

  #-----------------------------------------------------------------#
  #--- 5. Large h with Uniform kernel -> (sum e)^2               ---#
  #                                                                 #
  # When h exceeds the maximum pairwise distance, the Uniform       #
  # kernel assigns weight 1 to every pair, so:                      #
  #   V = sum_ij e_i * e_j = (sum e)^2.                            #
  #-----------------------------------------------------------------#

  if (verbose) cat("--- Large h, Uniform -> (sum e)^2 ---\n")

  h_large   <- 1e4
  exp_large <- sum(e5)^2
  V_large   <- conleyHAC_grid(e5, coords5, h_large, "uniform")
  .record(.near(V_large, exp_large),
          "large h, uniform -> (sum e)^2",
          sprintf("V = %.10g, (sum e)^2 = %.10g", V_large, exp_large))

  #-----------------------------------------------------------------#
  #--- 6. Input ordering invariance                              ---#
  #                                                                 #
  # Permuting rows of coords and the corresponding entries of       #
  # e_hat must not change the result.                               #
  #-----------------------------------------------------------------#

  if (verbose) cat("--- Input ordering invariance ---\n")

  perm <- sample(nrow(coords5))
  for (kern in kernels) {
    V1 <- conleyHAC_grid(e5,        coords5,                  h5, kern)
    V2 <- conleyHAC_grid(e5[perm],  coords5[perm, , drop = FALSE], h5, kern)
    .record(.near(V1, V2),
            paste0("ordering invariance [", kern, "]"),
            sprintf("original = %.10g, permuted = %.10g", V1, V2))
  }

  #-----------------------------------------------------------------#
  #--- 7. Non-square grid (3 x 7) and anisotropic spacing       ---#
  #       (dx = 1, dy = 2)                                         ---#
  #-----------------------------------------------------------------#

  if (verbose) cat("--- Non-square grid (3x7, dx=1 dy=2) ---\n")

  coords_ns <- as.matrix(expand.grid(x = 1:3, y = seq(0, 12, by = 2)))
  e_ns      <- rnorm(nrow(coords_ns))
  h_ns      <- 3.0
  for (kern in kernels) {
    bf  <- .brute(e_ns, coords_ns, h_ns, kern)
    fft <- conleyHAC_grid(e_ns, coords_ns, h_ns, kern)
    .record(.near(fft, bf),
            paste0("non-square anisotropic grid [", kern, "]"),
            sprintf("FFT = %.10g, BF = %.10g, diff = %.2e", fft, bf, abs(fft - bf)))
  }

  #-----------------------------------------------------------------#
  #--- 8. Wendland non-negativity across random inputs           ---#
  #                                                                 #
  # W(h) is PSD for Wendland, so e' W(h) e >= 0 must hold for      #
  # every e and every h.                                            #
  #-----------------------------------------------------------------#

  if (verbose) cat("--- Wendland V >= 0 (50 random inputs) ---\n")

  n_rnd <- 50L
  all_nonneg <- all(replicate(n_rnd, {
    e_r <- rnorm(m5 * m5)
    h_r <- runif(1, 0.1, 15)
    conleyHAC_grid(e_r, coords5, h_r, "wendland") >= -tol
  }))
  .record(all_nonneg, paste0("Wendland V >= 0 (", n_rnd, " random inputs)"))

  #-----------------------------------------------------------------#
  #--- 9. Wrapper auto-dispatch                                  ---#
  #-----------------------------------------------------------------#

  if (verbose) cat("--- Wrapper auto-dispatch ---\n")

  for (kern in kernels) {
    V_direct <- conleyHAC_grid(e5, coords5, h5, kern)
    V_wrap   <- conleyHAC(e5, coords5, h5, kern)
    .record(.near(V_direct, V_wrap),
            paste0("wrapper -> grid [", kern, "]"))
  }

  if (requireNamespace("sf", quietly = TRUE)) {
    pts5 <- sf::st_as_sf(as.data.frame(coords5),
                         coords = c("x", "y"), crs = NA_crs_)
    for (kern in kernels) {
      V_direct <- suppressWarnings(conleyHAC_sf(e5, pts5, h5, kern))
      V_wrap   <- suppressWarnings(conleyHAC(e5, pts5, h5, kern))
      .record(.near(V_direct, V_wrap),
              paste0("wrapper -> sf [", kern, "]"))
    }
  }

  #-----------------------------------------------------------------#
  #--- Summary                                                   ---#
  #-----------------------------------------------------------------#

  n_pass  <- sum(results)
  n_total <- length(results)
  if (verbose)
    cat(sprintf("\n%d / %d tests passed\n", n_pass, n_total))

  invisible(list(results = results, n_pass = n_pass, n_total = n_total))
}
