library(spatialHVI)

set.seed(17)
y <- matrix(rnorm(6 * 40), 6, 40)
x <- matrix(rnorm(40 * 70), 40, 70)
K <- tcrossprod(x) / 70
A <- matrix(c(0, 1, .3, 1, 0, .6, .3, .6, 0), 3)

fit <- fit_spatial_hvi(y, K, A, tol = 1e-7)
stopifnot(fit$converged)
points <- heritability_from_vi(fit)
vi <- vi_heritability_intervals(fit, draws = 1000, seed = 42)
lrvb <- lrvb_heritability_intervals(fit, verify_derivatives = TRUE)

# Explicit S17 coefficients illustrate the manuscript interval calculation.
# Coverage for a new design requires its own independent validation.
calibrated <- calibrate_heritability_intervals(lrvb, manuscript_interval_calibration())
stopifnot(identical(vi$h2, lrvb$intervals$h2), identical(vi$h2, calibrated$h2))
print(points[c("h2", "bilateral_mean", "H2")])
print(lrvb$diagnostics)
print(calibrated)
