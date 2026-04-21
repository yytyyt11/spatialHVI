# Internal simulation helpers.

resolve_group_config <- function(cfg) {
  c_traits <- cfg$c_traits %||% cfg$C
  seed_a <- cfg$seed_a %||% cfg$seedA

  missing_fields <- c(
    if (is.null(c_traits)) "c_traits/C",
    if (is.null(cfg$rho)) "rho",
    if (is.null(cfg$knn_k)) "knn_k",
    if (is.null(cfg$knn_sigma)) "knn_sigma",
    if (is.null(seed_a)) "seed_a/seedA",
    if (is.null(cfg$sigma2_pairs)) "sigma2_pairs"
  )

  if (length(missing_fields) > 0L) {
    stop(
      sprintf("Each group config must define: %s.", paste(missing_fields, collapse = ", ")),
      call. = FALSE
    )
  }

  list(
    c_traits = validate_positive_count(c_traits, arg = "c_traits"),
    rho = as.numeric(cfg$rho),
    knn_k = validate_positive_count(cfg$knn_k, arg = "knn_k"),
    knn_sigma = as.numeric(cfg$knn_sigma),
    seed_a = as.integer(seed_a),
    sigma2_pairs = as.numeric(cfg$sigma2_pairs)
  )
}

generate_one_halves_dataset <- function(
  snp_path,
  out_dir,
  c_traits,
  a_mat,
  d_a,
  rho,
  sigma_block,
  sigma2_gc,
  seed_base = 2025,
  run_id = 1
) {
  c_traits <- validate_positive_count(c_traits, arg = "c_traits")
  if (c_traits %% 2L != 0L) {
    stop("`c_traits` must be even.", call. = FALSE)
  }

  a_mat <- validate_square_matrix(a_mat, arg = "a_mat")
  d_a <- validate_square_matrix(d_a, arg = "d_a")
  sigma_block <- validate_square_matrix(sigma_block, arg = "sigma_block")
  sigma2_gc <- as.numeric(sigma2_gc)

  j_pairs <- c_traits / 2L
  if (nrow(a_mat) != j_pairs || nrow(d_a) != j_pairs) {
    stop("Dimensions of `a_mat` and `d_a` must match `c_traits / 2`.", call. = FALSE)
  }
  if (nrow(sigma_block) != 2L) {
    stop("`sigma_block` must be a 2 x 2 matrix.", call. = FALSE)
  }
  if (length(sigma2_gc) != c_traits) {
    stop("`sigma2_gc` must have length equal to `c_traits`.", call. = FALSE)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  genotype <- read_snp_matrix(snp_path)
  subject_count <- nrow(genotype)
  k_mat <- build_kinship_matrix(genotype)
  k_mat <- 0.5 * (k_mat + t(k_mat)) + diag(1e-8, subject_count)

  if (!is.null(seed_base)) {
    set.seed(seed_base + as.integer(run_id) * 11L + c_traits)
  }

  b_mat <- matrix(0, nrow = c_traits, ncol = subject_count)
  for (idx in seq_len(c_traits)) {
    b_mat[idx, ] <- as.numeric(
      mvtnorm::rmvnorm(1, sigma = sigma2_gc[idx] * k_mat)
    )
  }

  spatial_cov <- solve(d_a - rho * a_mat)
  full_cov <- kronecker(sigma_block, spatial_cov)
  full_cov <- 0.5 * (full_cov + t(full_cov)) + diag(1e-8, c_traits)

  epsilon <- matrix(0, nrow = c_traits, ncol = subject_count)
  for (idx in seq_len(subject_count)) {
    epsilon[, idx] <- as.numeric(mvtnorm::rmvnorm(1, sigma = full_cov))
  }

  y_mat <- b_mat + epsilon

  utils::write.csv(
    y_mat,
    file.path(out_dir, "simulated_phenotypes.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    b_mat,
    file.path(out_dir, "random_effects.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    epsilon,
    file.path(out_dir, "residuals.csv"),
    row.names = FALSE
  )
  utils::write.table(
    a_mat,
    file.path(out_dir, "A_used.csv"),
    row.names = FALSE,
    col.names = FALSE,
    sep = ","
  )

  invisible(
    list(
      n_subjects = subject_count,
      c_traits = c_traits,
      output_dir = out_dir
    )
  )
}

empty_generation_summary <- function() {
  data.frame(
    group = character(0),
    snp_file = character(0),
    out_dir = character(0),
    c_traits = integer(0),
    j_pairs = integer(0),
    rho_used = numeric(0),
    time_sec = numeric(0),
    stringsAsFactors = FALSE
  )
}

#' Generate Halves-Ordered Simulation Datasets in Batch
#'
#' Reads SNP files from one or more group directories, constructs a shared
#' adjacency matrix for each group, and writes simulated phenotypes and latent
#' components for every SNP file that matches `snp_pattern`.
#'
#' @param group_cfgs A named list of group configurations. Each element must
#'   define `c_traits` (or `C`), `rho`, `knn_k`, `knn_sigma`, `seed_a` (or
#'   `seedA`), and `sigma2_pairs`.
#' @param sigma_block A 2 by 2 residual covariance matrix shared across groups.
#' @param base_dir Directory containing the group folders.
#' @param summary_path Optional CSV path for the batch summary. Use `NULL` to
#'   skip writing the summary to disk.
#' @param snp_pattern Regular expression used to find SNP CSV files inside each
#'   group directory.
#' @param seed_offset Integer seed offset used when generating replicate-level
#'   randomness.
#'
#' @return A data frame summarizing all generated datasets.
#' @export
#'
#' @examples
#' \dontrun{
#' cfg <- list(
#'   output_SNP_demo = list(
#'     c_traits = 4,
#'     rho = 0.5,
#'     knn_k = 1,
#'     knn_sigma = 1,
#'     seed_a = 42,
#'     sigma2_pairs = c(1, 2)
#'   )
#' )
#' simulate_halves_batch(cfg, sigma_block = diag(c(1, 1)))
#' }
simulate_halves_batch <- function(
  group_cfgs,
  sigma_block,
  base_dir = ".",
  summary_path = file.path(base_dir, "batch_generation_summary.csv"),
  snp_pattern = "^SNP[0-9]+\\.csv$",
  seed_offset = 777
) {
  if (!is.list(group_cfgs) || length(group_cfgs) == 0L || is.null(names(group_cfgs))) {
    stop("`group_cfgs` must be a named, non-empty list.", call. = FALSE)
  }
  sigma_block <- validate_square_matrix(sigma_block, arg = "sigma_block")
  if (nrow(sigma_block) != 2L) {
    stop("`sigma_block` must be a 2 x 2 matrix.", call. = FALSE)
  }

  summary_rows <- list()
  row_count <- 0L

  for (group_name in names(group_cfgs)) {
    cfg <- resolve_group_config(group_cfgs[[group_name]])
    if (cfg$c_traits %% 2L != 0L) {
      stop(sprintf("Group `%s` has an odd `c_traits` value.", group_name), call. = FALSE)
    }

    group_path <- file.path(base_dir, group_name)
    if (!dir.exists(group_path)) {
      warning(sprintf("Group directory not found and was skipped: %s", group_path), call. = FALSE)
      next
    }

    snp_files <- order_snp_files(group_path, pattern = snp_pattern)
    if (length(snp_files) == 0L) {
      warning(sprintf("No SNP files matched in `%s`.", group_path), call. = FALSE)
      next
    }

    j_pairs <- cfg$c_traits / 2L
    a_raw <- make_a_knn(
      j_pairs = j_pairs,
      k = cfg$knn_k,
      sigma = cfg$knn_sigma,
      seed = cfg$seed_a
    )
    a_mat <- normalize_a_spectral(a_raw)
    d_a <- make_degree_matrix(a_mat)
    rho_used <- ensure_spd_rho(d_a = d_a, a_mat = a_mat, rho_init = cfg$rho)
    sigma2_gc <- expand_sigma2_pairs_halves(cfg$sigma2_pairs, cfg$c_traits)

    group_index <- match(group_name, names(group_cfgs))
    for (idx in seq_along(snp_files)) {
      snp_path <- snp_files[idx]
      snp_name <- tools::file_path_sans_ext(basename(snp_path))
      out_dir <- file.path(group_path, snp_name)

      start_time <- proc.time()[3]
      generate_one_halves_dataset(
        snp_path = snp_path,
        out_dir = out_dir,
        c_traits = cfg$c_traits,
        a_mat = a_mat,
        d_a = d_a,
        rho = rho_used,
        sigma_block = sigma_block,
        sigma2_gc = sigma2_gc,
        seed_base = seed_offset + group_index * 1000L,
        run_id = idx
      )
      elapsed <- proc.time()[3] - start_time

      row_count <- row_count + 1L
      summary_rows[[row_count]] <- data.frame(
        group = group_name,
        snp_file = basename(snp_path),
        out_dir = out_dir,
        c_traits = cfg$c_traits,
        j_pairs = j_pairs,
        rho_used = rho_used,
        time_sec = round(elapsed, 4),
        stringsAsFactors = FALSE
      )
    }
  }

  summary_df <- if (length(summary_rows) > 0L) {
    do.call(rbind, summary_rows)
  } else {
    empty_generation_summary()
  }

  if (!is.null(summary_path)) {
    ensure_parent_dir(summary_path)
    utils::write.csv(summary_df, summary_path, row.names = FALSE)
  }

  summary_df
}
