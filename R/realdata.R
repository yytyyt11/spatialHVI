#' Run the Real-Data Halves Workflow from a Brain Measure CSV
#'
#' Reads a brain-measure CSV file, aligns it with genotype or kinship inputs,
#' fits the heteroscedastic VI model, and optionally saves posterior summaries.
#'
#' @param brain_csv Path to a brain-measure CSV file.
#' @param use_k_from Either `"X"` to construct the kinship matrix from a genotype
#'   object or `"K"` to use a precomputed kinship matrix.
#' @param X Optional genotype matrix or data frame used when `use_k_from = "X"`.
#'   Subjects must be in rows. When a data frame is supplied, non-SNP columns are
#'   allowed as long as the SNP columns can be identified.
#' @param K Optional kinship matrix or data frame used when `use_k_from = "K"`.
#' @param a_mat Optional adjacency matrix. When `NULL`, it is estimated from the
#'   phenotype matrix.
#' @param brain_row_names Logical; if `TRUE`, treat the first CSV column as row
#'   names when reading `brain_csv`.
#' @param phenotype_layout Layout of `brain_csv`. Use `"traits_in_rows"` for the
#'   common layout with traits in rows and samples in columns, `"samples_in_rows"`
#'   for the transpose, or `"auto"` to infer the layout from trait names.
#' @param phenotype_id_col Optional sample ID column in `brain_csv` when
#'   `phenotype_layout = "samples_in_rows"`.
#' @param phenotype_sample_ids Optional sample IDs for `brain_csv` when traits are
#'   stored in rows and the sample IDs are not present in the file.
#' @param x_ids Optional sample IDs for `X` when they are not stored in the row
#'   names or an ID column.
#' @param k_ids Optional sample IDs for `K` when they are not stored in row and
#'   column names.
#' @param x_id_col Optional ID column name inside `X`.
#' @param x_snp_cols Optional SNP columns to use from `X`, supplied as column names
#'   or indices after removing `x_id_col`.
#' @param x_snp_pattern Regular expression used to detect SNP columns when
#'   `x_snp_cols` is `NULL`.
#' @param x_annotation Optional SNP annotation data frame used to identify SNP
#'   columns in `X`.
#' @param annotation_rsid_col Column of `x_annotation` containing SNP identifiers.
#' @param trait_order Trait ordering in the phenotype matrix. Use `"halves"` for
#'   `(L1, ..., LJ, R1, ..., RJ)`, `"pairs"` for `(L1, R1, L2, R2, ...)`, or
#'   `"auto"` to infer the order from trait names.
#' @param allow_row_order_alignment Logical; if `TRUE`, fall back to row-order
#'   alignment when no usable sample IDs are available.
#' @param na_action How to handle missing phenotype values. `"error"` stops with a
#'   clear message; `"zero"` replaces them with zero before fitting.
#' @param zscore_rows Logical; if `TRUE`, z-score each phenotype row after sample
#'   alignment.
#' @param k_knn Number of neighbors used when constructing `a_mat` from the
#'   phenotype matrix.
#' @param rho_grid Candidate rho values.
#' @param max_iter Maximum number of VI iterations.
#' @param tol Absolute ELBO tolerance used for convergence.
#' @param verbose Logical; if `TRUE`, prints ELBO progress.
#' @param prefix Optional output prefix passed to [save_vi_result()]. Use `NULL`
#'   to skip writing files.
#' @param make_plot Logical; if `TRUE`, plot the ELBO trace.
#' @param show_summary Logical; if `TRUE`, print a concise text summary.
#'
#' @return A list containing the aligned phenotype matrix, kinship matrix,
#'   adjacency matrix, fitted VI result, heritability summary, sample IDs, trait
#'   names, and alignment metadata.
#' @export
#'
#' @examples
#' \dontrun{
#' run_realdata_brain_hetero_halves(
#'   brain_csv = "path/to/BrainMeasure_data.csv",
#'   use_k_from = "K",
#'   K = diag(10),
#'   phenotype_layout = "traits_in_rows"
#' )
#' }
run_realdata_brain_hetero_halves <- function(
  brain_csv,
  use_k_from = c("X", "K"),
  X = NULL,
  K = NULL,
  a_mat = NULL,
  brain_row_names = TRUE,
  phenotype_layout = c("auto", "traits_in_rows", "samples_in_rows"),
  phenotype_id_col = NULL,
  phenotype_sample_ids = NULL,
  x_ids = NULL,
  k_ids = NULL,
  x_id_col = NULL,
  x_snp_cols = NULL,
  x_snp_pattern = "^rs",
  x_annotation = NULL,
  annotation_rsid_col = "RSID",
  trait_order = c("auto", "halves", "pairs"),
  allow_row_order_alignment = TRUE,
  na_action = c("error", "zero"),
  zscore_rows = TRUE,
  k_knn = 2,
  rho_grid = seq(0, 0.99, by = 0.01),
  max_iter = 10000,
  tol = 1e-4,
  verbose = TRUE,
  prefix = NULL,
  make_plot = FALSE,
  show_summary = FALSE
) {
  use_k_from <- match.arg(use_k_from)
  phenotype_layout <- match.arg(phenotype_layout)
  trait_order <- match.arg(trait_order)
  na_action <- match.arg(na_action)

  phenotype <- read_brain_measure_input(
    brain_csv = brain_csv,
    brain_row_names = brain_row_names,
    phenotype_layout = phenotype_layout,
    phenotype_id_col = phenotype_id_col,
    phenotype_sample_ids = phenotype_sample_ids,
    trait_order = trait_order,
    na_action = na_action
  )

  y_raw <- phenotype$y
  alignment_method <- NULL
  sample_ids <- phenotype$sample_ids
  x_snp_names <- NULL
  shared_trait_validation <- NULL

  if (use_k_from == "X") {
    if (is.null(X)) {
      stop("`X` must be supplied when `use_k_from = \"X\"`.", call. = FALSE)
    }

    x_prepared <- prepare_realdata_x(
      x = X,
      x_ids = x_ids,
      x_id_col = x_id_col,
      x_snp_cols = x_snp_cols,
      x_snp_pattern = x_snp_pattern,
      annotation = x_annotation,
      annotation_rsid_col = annotation_rsid_col
    )

    shared_trait_validation <- validate_shared_trait_alignment(
      y = y_raw,
      trait_names = phenotype$trait_names,
      shared_numeric = x_prepared$shared_numeric
    )

    aligned <- align_realdata_samples(
      y = y_raw,
      phenotype_sample_ids = phenotype$sample_ids,
      x_mat = x_prepared$X,
      x_ids = x_prepared$sample_ids,
      allow_row_order_alignment = allow_row_order_alignment,
      shared_validation = shared_trait_validation
    )

    y_aligned <- aligned$y
    x_aligned <- aligned$x_mat
    sample_ids <- aligned$sample_ids %||% x_prepared$sample_ids
    alignment_method <- aligned$method
    x_snp_names <- x_prepared$snp_names
    k_use <- build_kinship_matrix(x_aligned)
  } else {
    if (is.null(K)) {
      stop("`K` must be supplied when `use_k_from = \"K\"`.", call. = FALSE)
    }

    k_prepared <- prepare_realdata_k(K = K, k_ids = k_ids)
    aligned <- align_realdata_samples(
      y = y_raw,
      phenotype_sample_ids = phenotype$sample_ids,
      k_mat = k_prepared$K,
      k_ids = k_prepared$sample_ids,
      allow_row_order_alignment = allow_row_order_alignment
    )

    y_aligned <- aligned$y
    k_use <- aligned$k_mat
    sample_ids <- aligned$sample_ids %||% k_prepared$sample_ids
    alignment_method <- aligned$method
  }

  y_fit <- if (isTRUE(zscore_rows)) safe_scale_rows(y_aligned) else validate_even_trait_matrix(y_aligned, arg = "y_aligned")

  a_use <- if (is.null(a_mat)) {
    build_a_from_y_halves(y_fit, k = k_knn, trait_order = "halves")
  } else {
    normalize_a_spectral(a_mat)
  }

  res <- run_vi_hetero_from_mats(
    y = y_fit,
    k_mat = k_use,
    a_mat = a_use,
    rho_grid = rho_grid,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose
  )
  h <- heritability_from_vi(res, a_mat = a_use)

  if (isTRUE(make_plot)) {
    graphics::plot(
      seq_along(res$elbo),
      res$elbo,
      type = "l",
      main = "ELBO (real data, heteroscedastic, halves)",
      xlab = "Iteration",
      ylab = "ELBO"
    )
  }

  if (isTRUE(show_summary)) {
    message(format_heritability_summary(res, h))
  }

  if (!is.null(prefix)) {
    save_vi_result(res, prefix = prefix)
  }

  list(
    y = y_fit,
    k_mat = k_use,
    a_mat = a_use,
    res = res,
    h = h,
    sample_ids = sample_ids,
    trait_names = rownames(y_fit),
    alignment_method = alignment_method,
    phenotype_layout = phenotype$phenotype_layout,
    trait_order_detected = phenotype$trait_order_detected,
    trait_order_applied = phenotype$trait_order_applied,
    x_snp_names = x_snp_names,
    shared_trait_validation = shared_trait_validation
  )
}

#' Run the Real-Data Halves Workflow from CSV and RData Files
#'
#' Loads a brain-measure CSV and an `.RData` file, selects a genotype object or
#' kinship object, optionally subsets it, and then calls
#' [run_realdata_brain_hetero_halves()].
#'
#' @param brain_csv Path to a brain-measure CSV file.
#' @param rdata_path Path to an `.RData` file containing genotype, kinship, or
#'   annotation objects.
#' @param use_k_from Whether to use a genotype object (`"X"`), a kinship object
#'   (`"K"`), or `"auto"` detection.
#' @param x_object Optional object name in `rdata_path` to use as the genotype
#'   input.
#' @param k_object Optional object name in `rdata_path` to use as the kinship
#'   input.
#' @param genotype_object Backward-compatible alias for `x_object`.
#' @param annotation_object Optional object name in `rdata_path` containing SNP
#'   annotations such as an `RSID` column.
#' @param genotype_rows Optional row indices used to subset the selected genotype
#'   or kinship object.
#' @param genotype_cols Optional SNP column indices or names used to subset the
#'   selected genotype object after ID columns are removed.
#' @param brain_row_names Logical; if `TRUE`, treat the first CSV column as row
#'   names.
#' @param phenotype_layout Layout of `brain_csv`. See
#'   [run_realdata_brain_hetero_halves()].
#' @param phenotype_id_col Optional sample ID column in `brain_csv` when samples
#'   are stored in rows.
#' @param phenotype_sample_ids Optional sample IDs for `brain_csv` when traits are
#'   stored in rows.
#' @param x_id_col Optional ID column name inside `x_object`.
#' @param x_snp_cols Optional SNP columns to use from `x_object`.
#' @param x_snp_pattern Regular expression used to detect SNP columns when
#'   `x_snp_cols` is `NULL`.
#' @param annotation_rsid_col Column of `annotation_object` containing SNP
#'   identifiers.
#' @param trait_order Trait ordering in the phenotype matrix.
#' @param allow_row_order_alignment Logical; if `TRUE`, fall back to row-order
#'   alignment when no usable sample IDs are available.
#' @param na_action How to handle missing phenotype values.
#' @param zscore_rows Logical; if `TRUE`, z-score each phenotype row after
#'   alignment.
#' @param a_mat Optional adjacency matrix. When `NULL`, it is estimated from the
#'   phenotype matrix.
#' @param k_knn Number of neighbors used when constructing `a_mat`.
#' @param rho_grid Candidate rho values.
#' @param max_iter Maximum number of VI iterations.
#' @param tol Absolute ELBO tolerance used for convergence.
#' @param verbose Logical; if `TRUE`, prints ELBO progress.
#' @param prefix Optional output prefix passed to [save_vi_result()]. Use `NULL`
#'   to skip writing files.
#' @param make_plot Logical; if `TRUE`, plot the ELBO trace.
#' @param show_summary Logical; if `TRUE`, print a concise text summary.
#'
#' @return A list containing the aligned phenotype matrix, kinship matrix,
#'   adjacency matrix, fitted VI result, heritability summary, sample IDs, trait
#'   names, and alignment metadata.
#' @export
#'
#' @examples
#' \dontrun{
#' run_realdata_halves_from_files(
#'   brain_csv = "path/to/BrainMeasure_data.csv",
#'   rdata_path = "path/to/FreeSurfer_Data.RData",
#'   use_k_from = "X",
#'   x_object = "FreeSurfer_list_Data"
#' )
#' }
run_realdata_halves_from_files <- function(
  brain_csv,
  rdata_path,
  use_k_from = c("auto", "X", "K"),
  x_object = NULL,
  k_object = NULL,
  genotype_object = NULL,
  annotation_object = NULL,
  genotype_rows = NULL,
  genotype_cols = NULL,
  brain_row_names = TRUE,
  phenotype_layout = c("auto", "traits_in_rows", "samples_in_rows"),
  phenotype_id_col = NULL,
  phenotype_sample_ids = NULL,
  x_id_col = NULL,
  x_snp_cols = NULL,
  x_snp_pattern = "^rs",
  annotation_rsid_col = "RSID",
  trait_order = c("auto", "halves", "pairs"),
  allow_row_order_alignment = TRUE,
  na_action = c("error", "zero"),
  zscore_rows = TRUE,
  a_mat = NULL,
  k_knn = 2,
  rho_grid = seq(0, 0.99, by = 0.01),
  max_iter = 10000,
  tol = 1e-4,
  verbose = TRUE,
  prefix = NULL,
  make_plot = FALSE,
  show_summary = FALSE
) {
  use_k_from <- match.arg(use_k_from)
  phenotype_layout <- match.arg(phenotype_layout)
  trait_order <- match.arg(trait_order)
  na_action <- match.arg(na_action)

  if (is.null(x_object) && !is.null(genotype_object)) {
    x_object <- genotype_object
  }
  if (is.null(x_snp_cols) && !is.null(genotype_cols)) {
    x_snp_cols <- genotype_cols
  }

  load_info <- load_rdata_env(rdata_path)
  env <- load_info$env
  object_names <- load_info$objects

  if (is.null(annotation_object)) {
    annotation_object <- detect_rdata_object_name(env, object_names, kind = "annotation")
  }
  if (is.null(x_object) && use_k_from %in% c("auto", "X")) {
    x_object <- detect_rdata_object_name(env, object_names, kind = "x")
  }
  if (is.null(k_object) && use_k_from %in% c("auto", "K")) {
    k_object <- detect_rdata_object_name(env, object_names, kind = "k")
  }

  if (use_k_from == "auto") {
    use_k_from <- if (!is.null(k_object)) "K" else "X"
  }

  x_input <- NULL
  k_input <- NULL
  x_annotation <- NULL

  if (!is.null(annotation_object)) {
    if (!annotation_object %in% object_names) {
      stop(sprintf("Object `%s` was not found in `%s`.", annotation_object, rdata_path), call. = FALSE)
    }
    x_annotation <- get(annotation_object, envir = env, inherits = FALSE)
  }

  if (use_k_from == "X") {
    if (is.null(x_object) || !x_object %in% object_names) {
      stop("No genotype object could be identified in the supplied `.RData` file.", call. = FALSE)
    }
    x_input <- get(x_object, envir = env, inherits = FALSE)
    if (!is.null(genotype_rows)) {
      x_input <- x_input[genotype_rows, , drop = FALSE]
    }
  } else {
    if (is.null(k_object) || !k_object %in% object_names) {
      stop("No kinship object could be identified in the supplied `.RData` file.", call. = FALSE)
    }
    k_input <- get(k_object, envir = env, inherits = FALSE)
    if (!is.null(genotype_rows)) {
      k_input <- k_input[genotype_rows, genotype_rows, drop = FALSE]
    }
  }

  run_realdata_brain_hetero_halves(
    brain_csv = brain_csv,
    use_k_from = use_k_from,
    X = x_input,
    K = k_input,
    a_mat = a_mat,
    brain_row_names = brain_row_names,
    phenotype_layout = phenotype_layout,
    phenotype_id_col = phenotype_id_col,
    phenotype_sample_ids = phenotype_sample_ids,
    x_ids = NULL,
    k_ids = NULL,
    x_id_col = x_id_col,
    x_snp_cols = x_snp_cols,
    x_snp_pattern = x_snp_pattern,
    x_annotation = x_annotation,
    annotation_rsid_col = annotation_rsid_col,
    trait_order = trait_order,
    allow_row_order_alignment = allow_row_order_alignment,
    na_action = na_action,
    zscore_rows = zscore_rows,
    k_knn = k_knn,
    rho_grid = rho_grid,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose,
    prefix = prefix,
    make_plot = make_plot,
    show_summary = show_summary
  )
}
