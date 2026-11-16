make_manuscript_test_fit <- function() {
  set.seed(17)
  y <- matrix(rnorm(6 * 40), 6, 40)
  x <- matrix(rnorm(40 * 70), 40, 70)
  A <- matrix(c(0, 1, .3, 1, 0, .6, .3, .6, 0), 3)
  fit_spatial_hvi(y, tcrossprod(x) / 70, A, tol = 1e-7)
}

test_that("PCH uses the generalized eigenvalue and is invariant to pair order", {
  toy <- list(rho_grid = .8, rho_w = 1, E_sigma2 = c(.03, .2, .7, .06, .4, .9),
    nu = 10, S = diag(2), E_Sigma = matrix(c(1, .35, .35, .8), 2),
    mu = matrix(0, 6, 3), Sig_diag = 0, Sig_logdet = 0, Sig_trace_invK = 0,
    A = matrix(c(0, 1, .3, 1, 0, .6, .3, .6, 0), 3))
  original <- heritability_from_vi(toy)
  for (p in list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))) {
    idx <- c(p, p + 3)
    permuted <- toy; permuted$A <- toy$A[p, p]; permuted$E_sigma2 <- toy$E_sigma2[idx]
    h <- heritability_from_vi(permuted)
    G <- h$pieces$Sigma_g; E <- h$pieces$Sigma_e + diag(1e-10, 6)
    lambda <- max(Re(eigen(solve(E, G), only.values = TRUE)$values))
    expect_equal(h$H2, lambda / (1 + lambda), tolerance = 1e-12)
    expect_equal(h$H2, original$H2, tolerance = 1e-12)
    expect_equal(h$h2, original$h2[idx], tolerance = 1e-12)
  }
})

test_that("manuscript point summaries use mean inverse precision", {
  fit <- make_manuscript_test_fit()
  expect_true(fit$converged)
  expect_equal(fit$prior$a0, rep(2.1, 6))
  expect_equal(fit$prior$d0, rep(.33, 6))
  expect_equal(fit$prior$S0, diag(2, 2))
  summary <- heritability_from_vi(fit)
  qi <- Reduce(`+`, lapply(seq_along(fit$rho_grid), function(i)
    fit$rho_w[i] * solve(fit$D_A - fit$rho_grid[i] * fit$A)))
  residual <- rep(diag(qi), 2) * rep(diag(fit$E_Sigma), each = 3)
  expect_equal(summary$h2, fit$E_sigma2 / (fit$E_sigma2 + residual), tolerance = 1e-7)
  lambda <- max(Re(eigen(solve(summary$pieces$Sigma_e, summary$pieces$Sigma_g),
                         only.values = TRUE)$values))
  expect_equal(summary$H2, lambda / (1 + lambda), tolerance = 1e-12)
})

test_that("LRVB derivatives and bilateral covariance agree and centers stay fixed", {
  fit <- make_manuscript_test_fit()
  before <- serialize(fit, NULL)
  h <- heritability_from_vi(fit)
  lr <- lrvb_heritability_intervals(fit, verify_derivatives = TRUE)
  expect_identical(serialize(fit, NULL), before)
  expect_identical(lr$intervals$h2, c(h$h2, h$bilateral_mean))
  expect_identical(lr$original_rho_grid, fit$rho_grid)
  expect_lt(lr$diagnostics[['finite_difference_error']], 5e-4)
  expect_lt(lr$diagnostics[['tilt_refit_error']], 5e-4)
  W <- cbind(diag(3), diag(3)) / 2
  expect_equal(lr$bilateral_h_covariance, W %*% lr$h_covariance %*% t(W), tolerance = 1e-12)
  independent_sides <- (diag(lr$h_covariance)[1:3] + diag(lr$h_covariance)[4:6]) / 4
  expect_gt(max(abs(diag(lr$bilateral_h_covariance) - independent_sides)), 1e-8)
  cal <- calibrate_heritability_intervals(lr, manuscript_interval_calibration())
  expect_identical(cal$h2, lr$intervals$h2)
  expect_equal(cal$lower, plogis(qlogis(cal$h2) - cal$b * cal$sd_logit))
  expect_equal(cal$upper, plogis(qlogis(cal$h2) - cal$a * cal$sd_logit))
  expect_true(all(cal$lower <= cal$upper))
})

test_that("joint variational intervals preserve RNG and use each joint pair draw", {
  fit <- make_manuscript_test_fit()
  set.seed(112); before <- .Random.seed
  out <- vi_heritability_intervals(fit, draws = 100, seed = 42)
  expect_identical(.Random.seed, before)
  expect_identical(out, vi_heritability_intervals(fit, draws = 100, seed = 42))
  set.seed(42)
  index <- sample.int(length(fit$rho_grid), 100, TRUE, prob = fit$rho_w)
  reference <- matrix(0, 100, 9)
  for (b in 1:100) {
    g <- 1 / rgamma(6, shape = fit$a, rate = fit$d)
    Sigma <- chol2inv(chol(rWishart(1, fit$nu, chol2inv(chol(fit$S)))[, , 1]))
    qi <- chol2inv(chol(fit$D_A - fit$rho_grid[index[b]] * fit$A))
    e <- rep(diag(qi), 2) * rep(diag(Sigma), each = 3)
    h <- g / (g + e)
    reference[b, ] <- c(h, (h[1:3] + h[4:6]) / 2)
  }
  limits <- t(apply(reference, 2, quantile, c(.025, .975), type = 7, names = FALSE))
  expect_equal(out$lower, limits[, 1], tolerance = 1e-14)
  expect_equal(out$upper, limits[, 2], tolerance = 1e-14)
  out8 <- vi_heritability_intervals(fit, draws = 100, seed = 42, quantile_type = 8)
  expect_equal(out8$lower, apply(reference, 2, quantile, .025, type = 8, names = FALSE), tolerance = 1e-14)
})

test_that("calibration uses fixed full precision coefficients and development errors", {
  c <- manuscript_interval_calibration()
  expect_identical(c$coefficients$a, c(-2.8135915281262953, -3.7005660188107603))
  expect_identical(c$coefficients$b, c(9.18743593573288, 13.273540095862245))
  d <- data.frame(dataset_id = rep(1:3, each = 2),
    target = rep(c("trait", "bilateral_mean"), 3),
    h2 = .2, truth = c(.05, .1, .2, .25, .4, .5), sd_logit = .5)
  learned <- learn_interval_calibration(d)
  for (i in 1:2) {
    s <- d[d$target == learned$coefficients$target[i], ]
    q <- quantile((qlogis(s$h2) - qlogis(s$truth)) / s$sd_logit, c(.025, .975), type = 7, names = FALSE)
    expect_equal(unlist(learned$coefficients[i, c('a', 'b')], use.names = FALSE), q)
  }
  expect_error(learn_interval_calibration(rbind(d, d[1, ])), "equal target counts")
  d$truth[1] <- 0
  expect_error(learn_interval_calibration(d), "interior")
})

test_that("fixed rho ablation and invalid inputs are handled explicitly", {
  set.seed(19)
  y <- matrix(rnorm(6 * 20), 6, 20)
  K <- diag(20)
  A <- matrix(c(0, 1, .3, 1, 0, .6, .3, .6, 0), 3)
  f <- fit_spatial_hvi(y, K, A, rho_grid = 0, max_iter = 2)
  expect_identical(f$rho_grid, 0)
  expect_equal(f$D_A, diag(rowSums(A)))
  expect_false(isTRUE(all.equal(f$D_A, diag(3))))
  expect_error(lrvb_heritability_intervals(f), "not converged")
  colnames(y) <- paste0('s', 1:20)
  dimnames(K) <- list(rev(colnames(y)), rev(colnames(y)))
  expect_error(fit_spatial_hvi(y, K, A), "subject order")
  A[1, 1] <- 1
  expect_error(fit_spatial_hvi(unname(y), unname(K), A), "zero diagonal")
})
