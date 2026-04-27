#' Fit the Heteroscedastic Halves VI Model
#'
#' Fits the fast exact variational inference model for halves-ordered traits.
#'
#' @param y Numeric phenotype matrix \eqn{Y} with traits in rows and subjects in
#'   columns. In the manuscript notation this is the \eqn{C \times n} phenotype
#'   matrix, where \eqn{C = 2J} is even and \eqn{n} is the number of subjects.
#'   The rows must already be in halves order `(L1, ..., LJ, R1, ..., RJ)`;
#'   this fitting function does not reorder `y`. The column order must be the
#'   same subject order used in `k_mat`. If traits are supplied in columns,
#'   rows are in alternating pair order, or subjects are not aligned with
#'   `k_mat`, the fitted genetic and residual covariance quantities will refer
#'   to the wrong model axes or the function will stop on a dimension check.
#' @param k_mat Numeric kinship matrix \eqn{K} with dimension \eqn{n \times n}.
#'   Rows and columns represent subjects in exactly the same order as the
#'   columns of `y`. The matrix is symmetrized and scaled internally to have
#'   mean diagonal one. Use [build_kinship_matrix()] to construct it from SNP
#'   dosages, or provide a precomputed genomic relationship matrix. A dimension
#'   mismatch stops the fit; a silently permuted subject order produces invalid
#'   heritability estimates.
#' @param a_mat Optional numeric pair-level spatial adjacency matrix \eqn{A} of
#'   dimension \eqn{J \times J}, where \eqn{J = C/2}. Entry `(i, j)` encodes the
#'   neighborhood weight between paired ROIs `i` and `j`; the diagonal is treated
#'   as zero and the matrix is spectrally normalized internally. When `NULL`, an
#'   adjacency matrix is estimated from `y` using [build_a_from_y_halves()] with
#'   its defaults. Provide `a_mat` when anatomical or simulation adjacency is
#'   known and should define the residual precision \eqn{D_A - \rho A}. A wrong
#'   pair order or wrong dimension stops the fit or assigns spatial dependence
#'   to the wrong ROI pairs.
#' @param rho_grid Numeric vector of candidate spatial dependence values for
#'   \eqn{\rho}. The variational posterior for \eqn{\rho} is represented on
#'   this discrete grid with a uniform prior over values that make
#'   \eqn{D_A - \rho A} positive definite. The default `seq(0, 0.99, by = 0.01)`
#'   searches nonnegative spatial dependence below one. Use a denser or narrower
#'   grid when prior scientific knowledge or computation time warrants it. Empty,
#'   non-finite, or entirely invalid grids stop the fit; overly coarse grids can
#'   blur the posterior mean `E_rho`.
#' @param max_iter Positive integer maximum number of coordinate-ascent VI
#'   iterations. The default `5000` favors convergence for routine analyses.
#'   Increase it if the ELBO has not stabilized; decrease it for smoke tests or
#'   examples. Very small values may return an intentionally under-converged
#'   posterior.
#' @param tol Non-negative numeric scalar ELBO tolerance. Iteration stops when
#'   the absolute change in ELBO is smaller than `tol`; the default `1e-4` is a
#'   practical convergence threshold. Use `0` to force the algorithm to run
#'   until `max_iter` in tests or timing experiments. Larger values stop earlier
#'   but may reduce posterior accuracy.
#' @param verbose Logical scalar. When `TRUE` (default), the function prints the
#'   ELBO at each iteration; set to `FALSE` for batch runs, tests, or scripted
#'   pipelines where per-iteration messages would be noisy.
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
#' @param y Numeric phenotype matrix \eqn{Y} with traits in rows and subjects in
#'   columns. In manuscript notation this is the \eqn{C \times n} phenotype
#'   matrix, where \eqn{C = 2J} is even and \eqn{n} is the number of subjects.
#'   Rows must already be in halves order `(L1, ..., LJ, R1, ..., RJ)`;
#'   this function does not reorder traits. In the homoscedastic model, all
#'   trait rows share one genetic variance component, so a wrong trait layout
#'   affects the common genetic variance estimate as well as the residual
#'   spatial covariance. The column order must match the subject order in
#'   `k_mat`.
#' @param k_mat Numeric kinship matrix \eqn{K} with dimension \eqn{n \times n}.
#'   Rows and columns represent subjects in exactly the same order as the
#'   columns of `y`. The matrix is symmetrized and scaled internally to have
#'   mean diagonal one. Use [build_kinship_matrix()] for SNP-derived kinship or
#'   provide a precomputed genomic relationship matrix. A dimension mismatch
#'   stops the fit; a permuted subject order invalidates the genetic random
#'   effect covariance.
#' @param a_mat Optional numeric pair-level spatial adjacency matrix \eqn{A}
#'   with dimension \eqn{J \times J}, where \eqn{J = C/2}. It defines
#'   neighborhood weights among ROI pairs in the residual precision
#'   \eqn{D_A - \rho A}; the diagonal is treated as zero and the matrix is
#'   spectrally normalized internally. When `NULL`, `A` is estimated from `y`
#'   using [build_a_from_y_halves()]. Provide `a_mat` when an anatomical,
#'   simulation, or otherwise fixed pair-level adjacency should be used.
#' @param rho_grid Numeric vector of candidate values for the spatial dependence
#'   parameter \eqn{\rho}. The variational posterior for \eqn{\rho} is stored as
#'   weights over this discrete grid after excluding values that make
#'   \eqn{D_A - \rho A} non-positive-definite. The default
#'   `seq(0, 0.99, by = 0.01)` searches nonnegative spatial dependence. Use a
#'   narrower or denser grid for sensitivity analysis or faster focused fits.
#' @param max_iter Positive integer maximum number of coordinate-ascent VI
#'   iterations. The default `5000` is intended for routine matrix fits.
#'   Increase it if ELBO diagnostics still change at the limit; reduce it for
#'   examples or smoke tests.
#' @param tol Non-negative numeric ELBO convergence tolerance. Iteration stops
#'   when the absolute ELBO change is below `tol`; the default is `1e-4`. Use
#'   `0` to force running to `max_iter`. Larger values stop earlier but may
#'   return a less stable homoscedastic posterior.
#' @param verbose Logical scalar. If `TRUE` (default), print per-iteration ELBO
#'   progress. Set `FALSE` for scripts, tests, or large batches.
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
#' @param y Numeric simulated phenotype matrix \eqn{Y} with dimension
#'   \eqn{C \times n}: traits in rows, subjects in columns, and rows in halves
#'   order `(L1, ..., LJ, R1, ..., RJ)`. This should usually be the
#'   `simulated_phenotypes.csv` output generated by [simulate_halves_batch()].
#'   Its subject order must match `k_mat`. If a transposed matrix or pairs-order
#'   matrix is supplied, the VI model and heritability summary will be fit to
#'   the wrong trait structure or fail dimension checks.
#' @param k_mat Numeric kinship matrix \eqn{K} with dimension \eqn{n \times n}
#'   in the same subject order as the columns of `y`. In simulation workflows it
#'   is usually computed from the SNP file by [build_kinship_matrix()]. The fit
#'   stops if `nrow(k_mat) != ncol(y)`.
#' @param a_true Numeric adjacency matrix \eqn{A} with dimension
#'   \eqn{J \times J}, where \eqn{J = C/2}. This is the pair-level adjacency used
#'   to generate the simulated residual covariance, typically read from
#'   `A_used.csv`. It is passed as the fixed `a_mat` for fitting and for
#'   [heritability_from_vi()], so its pair order must match the halves ordering
#'   of `y`.
#' @param rho_grid Numeric vector of candidate \eqn{\rho} values used for the
#'   discrete variational posterior. The default `seq(0, 0.99, by = 0.01)` is a
#'   broad nonnegative search grid. Narrow it around the simulation truth for
#'   faster experiments, or widen/refine it when assessing sensitivity.
#' @param max_iter Positive integer maximum number of VI iterations; default
#'   `5000`. Increase when the ELBO continues changing at the limit; reduce for
#'   quick simulation smoke tests.
#' @param tol Non-negative numeric ELBO convergence tolerance; default `1e-4`.
#'   Use `0` to disable early stopping before `max_iter`. A large tolerance can
#'   stop before the simulated posterior has stabilized.
#' @param verbose Logical scalar controlling ELBO messages. Default `TRUE`;
#'   set `FALSE` when fitting many simulated SNP datasets.
#' @param prefix Optional character file prefix passed to [save_vi_result()].
#'   When `NULL` (default), no result files are written. When supplied, files
#'   such as `"<prefix>_full_result.rds"` are created and existing files with
#'   those names are overwritten by the underlying write calls.
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
#' @param group_dir Character scalar path to one simulation group directory,
#'   such as a directory containing `SNP1.csv`, `SNP2.csv`, and per-SNP
#'   subdirectories produced by [simulate_halves_batch()]. For each matched SNP
#'   file, the function expects `group_dir/SNP*/simulated_phenotypes.csv` with
#'   \eqn{Y} in \eqn{C \times n} layout and `group_dir/SNP*/A_used.csv` with
#'   the corresponding \eqn{J \times J} adjacency matrix. The genotype SNP file
#'   itself is read to rebuild \eqn{K}. Missing required files cause that SNP to
#'   be skipped.
#' @param resume_from_snp Positive integer SNP id used as a resume point. The
#'   default `1L` fits all matched SNP files. Set, for example,
#'   `resume_from_snp = 100` to skip `SNP1.csv` through `SNP99.csv` after a
#'   partial batch run.
#' @param skip_completed Logical scalar. If `TRUE`, an existing
#'   `"<result_stub>_full_result.rds"` file inside a SNP output directory is read
#'   instead of refitting that SNP. This is useful for resuming long batches, but
#'   risky if fitting settings, package version, `rho_grid`, or inputs changed
#'   since the saved result was created. The default `FALSE` refits from inputs.
#' @param rho_grid Numeric candidate grid for the spatial dependence parameter
#'   \eqn{\rho}; passed to [fit_simulation_case_halves()]. The default
#'   `seq(0, 0.99, by = 0.01)` is broad. Use the same grid across SNPs when
#'   comparing posterior `E_rho` values.
#' @param max_iter Positive integer maximum VI iterations per SNP dataset.
#'   Default `5000`. Reduce for debugging; increase for difficult batches whose
#'   ELBO traces do not stabilize.
#' @param tol Non-negative numeric ELBO tolerance per SNP fit. Default `1e-4`;
#'   set to `0` for fixed-iteration simulation experiments.
#' @param verbose Logical scalar controlling per-iteration ELBO messages from
#'   each SNP fit. Default `TRUE`; set `FALSE` for production batch logs.
#' @param result_stub Character scalar file-name stem for result files written
#'   inside each per-SNP output folder. The default `"vi_halves_fast"` creates
#'   files such as `vi_halves_fast_full_result.rds` and
#'   `vi_halves_fast_rho_posterior.csv`. Changing this lets multiple analyses
#'   coexist in the same SNP directories.
#' @param summary_path Optional character path for the batch summary CSV. The
#'   default writes `"<result_stub>_summary.csv"` inside `group_dir`; `NULL`
#'   returns the data frame without writing a summary. Existing files are
#'   overwritten.
#' @param snp_pattern Character regular expression used to select SNP input
#'   files in `group_dir`. The default `"^SNP[0-9]+\\.csv$"` expects names such
#'   as `SNP1.csv`. Change it only when your simulation files use a different
#'   naming convention; the numeric SNP id must still be recoverable by the
#'   package's `SNP<integer>.csv` parser.
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
