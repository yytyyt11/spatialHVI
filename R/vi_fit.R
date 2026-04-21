#' Fit the Heteroscedastic Halves VI Model
#'
#' Fits the fast exact variational inference model for halves-ordered traits.
#'
#' @param y A numeric matrix with traits in rows and subjects in columns, ordered
#'   as `(L1, ..., LJ, R1, ..., RJ)`.
#' @param k_mat A subject-by-subject kinship matrix.
#' @param a_mat Optional `J` by `J` adjacency matrix. When `NULL`, it is estimated
#'   from `y`.
#' @param rho_grid Candidate rho values.
#' @param max_iter Maximum number of VI iterations.
#' @param tol Absolute ELBO tolerance used for convergence.
#' @param verbose Logical; if `TRUE`, prints ELBO progress.
#'
#' @return A list containing posterior summaries, ELBO history, and fitted
#'   hyperparameters.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(48), nrow = 8)
#' k_mat <- build_kinship_matrix(x)
#' y <- matrix(rnorm(32), nrow = 4)
#' fit <- run_vi_hetero_from_mats(
#'   y = y,
#'   k_mat = k_mat,
#'   rho_grid = seq(0, 0.2, by = 0.1),
#'   max_iter = 2,
#'   tol = 0,
#'   verbose = FALSE
#' )
#' fit$E_rho
run_vi_hetero_from_mats <- function(
  y,
  k_mat,
  a_mat = NULL,
  rho_grid = seq(0, 0.99, by = 0.01),
  max_iter = 5000,
  tol = 1e-4,
  verbose = TRUE
) {
  run_vi_core(
    y = y,
    k_mat = k_mat,
    a_mat = a_mat,
    rho_grid = rho_grid,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose,
    model = "hetero"
  )
}

#' Fit the Homoscedastic Halves VI Model
#'
#' Fits the fast homoscedastic variational inference model for halves-ordered
#' traits.
#'
#' @inheritParams run_vi_hetero_from_mats
#'
#' @return A list containing posterior summaries, ELBO history, and fitted
#'   hyperparameters.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(48), nrow = 8)
#' k_mat <- build_kinship_matrix(x)
#' y <- matrix(rnorm(32), nrow = 4)
#' fit <- run_vi_homo_from_mats(
#'   y = y,
#'   k_mat = k_mat,
#'   rho_grid = seq(0, 0.2, by = 0.1),
#'   max_iter = 2,
#'   tol = 0,
#'   verbose = FALSE
#' )
#' fit$model
run_vi_homo_from_mats <- function(
  y,
  k_mat,
  a_mat = NULL,
  rho_grid = seq(0, 0.99, by = 0.01),
  max_iter = 5000,
  tol = 1e-4,
  verbose = TRUE
) {
  run_vi_core(
    y = y,
    k_mat = k_mat,
    a_mat = a_mat,
    rho_grid = rho_grid,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose,
    model = "homo"
  )
}

#' Fit One Simulated Halves Dataset
#'
#' Runs the heteroscedastic VI model for a single simulated dataset and returns
#' both the fitted variational object and heritability summaries.
#'
#' @param y A halves-ordered phenotype matrix.
#' @param k_mat A subject-by-subject kinship matrix.
#' @param a_true The adjacency matrix used to simulate the data.
#' @param rho_grid Candidate rho values.
#' @param max_iter Maximum number of VI iterations.
#' @param tol Absolute ELBO tolerance used for convergence.
#' @param verbose Logical; if `TRUE`, prints ELBO progress.
#' @param prefix Optional file prefix passed to [save_vi_result()]. Use `NULL` to
#'   skip writing files.
#'
#' @return A list with elements `res` and `h`.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(48), nrow = 8)
#' y <- matrix(rnorm(32), nrow = 4)
#' a_true <- matrix(c(0, 1, 1, 0), nrow = 2)
#' out <- fit_simulation_case_halves(
#'   y = y,
#'   k_mat = build_kinship_matrix(x),
#'   a_true = a_true,
#'   rho_grid = seq(0, 0.2, by = 0.1),
#'   max_iter = 2,
#'   tol = 0,
#'   verbose = FALSE
#' )
#' out$h$H2
fit_simulation_case_halves <- function(
  y,
  k_mat,
  a_true,
  rho_grid = seq(0, 0.99, by = 0.01),
  max_iter = 5000,
  tol = 1e-4,
  verbose = TRUE,
  prefix = NULL
) {
  res <- run_vi_hetero_from_mats(
    y = y,
    k_mat = k_mat,
    a_mat = a_true,
    rho_grid = rho_grid,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose
  )
  h <- heritability_from_vi(res, a_mat = a_true)

  if (!is.null(prefix)) {
    save_vi_result(res, prefix = prefix)
  }

  list(res = res, h = h)
}

build_batch_summary_row <- function(snp_name, res, a_path, result_path) {
  a_true <- as.matrix(utils::read.csv(a_path, header = FALSE))
  storage.mode(a_true) <- "double"
  h <- heritability_from_vi(res, a_mat = a_true)

  data.frame(
    snp = snp_name,
    E_rho = res$E_rho,
    H2 = h$H2,
    mean_h2 = mean(h$h2),
    sd_h2 = stats::sd(h$h2),
    result_path = result_path,
    stringsAsFactors = FALSE
  )
}

empty_batch_summary <- function() {
  data.frame(
    snp = character(0),
    E_rho = numeric(0),
    H2 = numeric(0),
    mean_h2 = numeric(0),
    sd_h2 = numeric(0),
    result_path = character(0),
    stringsAsFactors = FALSE
  )
}

#' Fit a Batch of Simulated Halves Datasets from Disk
#'
#' Iterates over `SNP*.csv` files in one simulation group directory, reads the
#' generated phenotype and adjacency files, runs the heteroscedastic VI model,
#' and writes a summary CSV.
#'
#' @param group_dir Path to a simulation group directory.
#' @param resume_from_snp Integer SNP id used as a resume point.
#' @param skip_completed Logical; if `TRUE`, reuse existing saved results when
#'   available.
#' @param rho_grid Candidate rho values.
#' @param max_iter Maximum number of VI iterations.
#' @param tol Absolute ELBO tolerance used for convergence.
#' @param verbose Logical; if `TRUE`, prints ELBO progress.
#' @param result_stub File-name stem used for saved result files inside each SNP
#'   output folder.
#' @param summary_path Optional CSV path for the summary table. Use `NULL` to skip
#'   writing the summary to disk.
#' @param snp_pattern Regular expression used to detect SNP input files.
#'
#' @return A data frame summarizing the fitted SNP datasets.
#' @export
#'
#' @examples
#' \dontrun{
#' run_simulation_batch_halves("output_SNP1000495")
#' }
run_simulation_batch_halves <- function(
  group_dir,
  resume_from_snp = 1L,
  skip_completed = FALSE,
  rho_grid = seq(0, 0.99, by = 0.01),
  max_iter = 5000,
  tol = 1e-4,
  verbose = TRUE,
  result_stub = "vi_halves_fast",
  summary_path = file.path(group_dir, paste0(result_stub, "_summary.csv")),
  snp_pattern = "^SNP[0-9]+\\.csv$"
) {
  if (!dir.exists(group_dir)) {
    stop(sprintf("Directory not found: %s", group_dir), call. = FALSE)
  }

  resume_from_snp <- validate_positive_count(resume_from_snp, arg = "resume_from_snp")
  snp_files <- order_snp_files(group_dir, pattern = snp_pattern)
  if (length(snp_files) == 0L) {
    warning(sprintf("No SNP files matched in `%s`.", group_dir), call. = FALSE)
    summary_df <- empty_batch_summary()
    if (!is.null(summary_path)) {
      ensure_parent_dir(summary_path)
      utils::write.csv(summary_df, summary_path, row.names = FALSE)
    }
    return(summary_df)
  }

  summary_rows <- list()
  row_count <- 0L

  for (snp_path in snp_files) {
    snp_id <- extract_snp_id(snp_path)
    if (snp_id < resume_from_snp) {
      next
    }

    snp_name <- tools::file_path_sans_ext(basename(snp_path))
    data_dir <- file.path(group_dir, snp_name)
    y_path <- file.path(data_dir, "simulated_phenotypes.csv")
    a_path <- file.path(data_dir, "A_used.csv")
    result_path <- file.path(data_dir, paste0(result_stub, "_full_result.rds"))

    if (!file.exists(y_path) || !file.exists(a_path)) {
      warning(sprintf("Skipping `%s` because required inputs are missing.", snp_name), call. = FALSE)
      next
    }

    if (isTRUE(skip_completed) && file.exists(result_path)) {
      res <- readRDS(result_path)
    } else {
      y <- read_numeric_csv_matrix(y_path, row_names = FALSE)
      genotype <- read_snp_matrix(snp_path)
      k_mat <- build_kinship_matrix(genotype)
      a_true <- as.matrix(utils::read.csv(a_path, header = FALSE))
      storage.mode(a_true) <- "double"

      fit <- fit_simulation_case_halves(
        y = y,
        k_mat = k_mat,
        a_true = a_true,
        rho_grid = rho_grid,
        max_iter = max_iter,
        tol = tol,
        verbose = verbose,
        prefix = file.path(data_dir, result_stub)
      )
      res <- fit$res
    }

    row_count <- row_count + 1L
    summary_rows[[row_count]] <- build_batch_summary_row(
      snp_name = snp_name,
      res = res,
      a_path = a_path,
      result_path = result_path
    )
  }

  summary_df <- if (length(summary_rows) > 0L) {
    do.call(rbind, summary_rows)
  } else {
    empty_batch_summary()
  }

  if (nrow(summary_df) > 0L) {
    summary_df <- summary_df[order(as.integer(sub("^SNP", "", summary_df$snp))), , drop = FALSE]
  }

  if (!is.null(summary_path)) {
    ensure_parent_dir(summary_path)
    utils::write.csv(summary_df, summary_path, row.names = FALSE)
  }

  summary_df
}
