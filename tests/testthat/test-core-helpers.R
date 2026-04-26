test_that("build_kinship_matrix returns a symmetric matrix", {
  set.seed(1)
  snp <- matrix(sample(0:2, 24, replace = TRUE), nrow = 6)
  k_mat <- build_kinship_matrix(snp)

  expect_equal(dim(k_mat), c(6, 6))
  expect_equal(k_mat, t(k_mat))
})

test_that("build_a_from_y_halves validates trait count", {
  y <- matrix(rnorm(15), nrow = 3)
  expect_error(build_a_from_y_halves(y), "even number of rows")
})

test_that("build_a_from_y_halves can reorder alternating trait pairs", {
  set.seed(2)
  y <- matrix(rnorm(32), nrow = 4)
  rownames(y) <- c("Left_A", "Right_A", "Left_B", "Right_B")

  a_mat <- build_a_from_y_halves(y, k = 1, trait_order = "auto")

  expect_equal(dim(a_mat), c(2, 2))
  expect_equal(a_mat, t(a_mat))
})

test_that("simulate_halves_batch writes expected simulation files", {
  root_dir <- file.path(tempdir(), paste0("spatialHVI-sim-", Sys.getpid(), "-", as.integer(stats::runif(1, 1, 1e6))))
  dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)
  group_dir <- file.path(root_dir, "group_a")
  dir.create(group_dir, showWarnings = FALSE)

  snp <- matrix(sample(0:2, 48, replace = TRUE), nrow = 8)
  utils::write.csv(snp, file.path(group_dir, "SNP1.csv"), row.names = FALSE)

  cfg <- list(
    group_a = list(
      c_traits = 4,
      rho = 0.3,
      knn_k = 1,
      knn_sigma = 1,
      seed_a = 42,
      sigma2_pairs = c(1, 2)
    )
  )

  summary_df <- simulate_halves_batch(
    group_cfgs = cfg,
    sigma_block = diag(c(1, 1)),
    base_dir = root_dir,
    summary_path = file.path(root_dir, "generation_summary.csv")
  )

  expect_equal(nrow(summary_df), 1)
  expect_true(file.exists(file.path(group_dir, "SNP1", "simulated_phenotypes.csv")))
  expect_true(file.exists(file.path(group_dir, "SNP1", "A_used.csv")))
})
