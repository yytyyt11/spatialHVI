test_that("run_vi_hetero_from_mats returns a valid fit", {
  set.seed(1)
  x <- matrix(rnorm(48), nrow = 8)
  y <- matrix(rnorm(32), nrow = 4)

  fit <- run_vi_hetero_from_mats(
    y = y,
    k_mat = build_kinship_matrix(x),
    rho_grid = seq(0, 0.2, by = 0.1),
    max_iter = 2,
    tol = 0,
    verbose = FALSE
  )

  expect_true(is.list(fit))
  expect_length(fit$E_sigma2, 4)
  expect_true(all(is.finite(fit$E_sigma2)))
  expect_true(fit$E_rho >= 0)

  h <- heritability_from_vi(fit)
  expect_length(h$h2, 4)
  expect_true(h$H2 >= 0 && h$H2 <= 1)
})

test_that("run_vi_homo_from_mats validates dimensions", {
  y <- matrix(rnorm(15), nrow = 3)
  expect_error(
    run_vi_homo_from_mats(y = y, k_mat = diag(5), max_iter = 1, verbose = FALSE),
    "even number of rows"
  )
})
