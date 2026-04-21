test_that("run_simulation_batch_halves fits a generated dataset", {
  root_dir <- file.path(tempdir(), paste0("mypackage-fit-", Sys.getpid(), "-", as.integer(stats::runif(1, 1, 1e6))))
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

  simulate_halves_batch(
    group_cfgs = cfg,
    sigma_block = diag(c(1, 1)),
    base_dir = root_dir,
    summary_path = NULL
  )

  summary_df <- run_simulation_batch_halves(
    group_dir = group_dir,
    rho_grid = seq(0, 0.2, by = 0.1),
    max_iter = 1,
    tol = 0,
    verbose = FALSE
  )

  expect_equal(nrow(summary_df), 1)
  expect_true(file.exists(file.path(group_dir, "SNP1", "vi_halves_fast_full_result.rds")))
})

test_that("run_realdata_halves_from_files loads the genotype object and fits", {
  root_dir <- file.path(tempdir(), paste0("mypackage-real-", Sys.getpid(), "-", as.integer(stats::runif(1, 1, 1e6))))
  dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)

  brain_path <- file.path(root_dir, "brain.csv")
  brain_mat <- matrix(rnorm(32), nrow = 4)
  rownames(brain_mat) <- c("Left_A", "Left_B", "Right_A", "Right_B")
  utils::write.csv(as.data.frame(brain_mat), brain_path, row.names = TRUE)

  FreeSurfer_list_Data <- data.frame(
    PTID = paste0("ID", 1:8),
    rs1 = sample(0:2, 8, replace = TRUE),
    rs2 = sample(0:2, 8, replace = TRUE),
    rs3 = sample(0:2, 8, replace = TRUE),
    rs4 = sample(0:2, 8, replace = TRUE),
    rs5 = sample(0:2, 8, replace = TRUE),
    rs6 = sample(0:2, 8, replace = TRUE),
    Left_A = as.numeric(brain_mat["Left_A", ]),
    Right_A = as.numeric(brain_mat["Right_A", ]),
    stringsAsFactors = FALSE
  )
  rdata_path <- file.path(root_dir, "FreeSurfer_Data.RData")
  save(FreeSurfer_list_Data, file = rdata_path)

  fit <- run_realdata_halves_from_files(
    brain_csv = brain_path,
    rdata_path = rdata_path,
    use_k_from = "X",
    x_object = "FreeSurfer_list_Data",
    x_id_col = "PTID",
    brain_row_names = TRUE,
    rho_grid = seq(0, 0.2, by = 0.1),
    max_iter = 1,
    tol = 0,
    verbose = FALSE,
    prefix = NULL,
    make_plot = FALSE,
    show_summary = FALSE
  )

  expect_true(is.list(fit))
  expect_true(all(c("y", "k_mat", "a_mat", "res", "h") %in% names(fit)))
  expect_equal(fit$alignment_method, "shared_traits")
  expect_equal(fit$sample_ids, paste0("ID", 1:8))
})
