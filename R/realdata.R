#' Run the Real-Data Halves Workflow from a Brain Measure CSV
#'
#' Reads a brain-measure CSV file, aligns it with genotype or kinship inputs,
#' fits the heteroscedastic VI model, and optionally saves posterior summaries.
#'
#' @param brain_csv Character scalar path to the brain-measure phenotype CSV.
#'   The file must contain the observed phenotype matrix \eqn{Y}. With
#'   `phenotype_layout = "traits_in_rows"`, rows are traits and columns are
#'   subjects; with `"samples_in_rows"`, rows are subjects and columns are
#'   traits before the function transposes to the required \eqn{C \times n}
#'   layout. Trait names should allow left/right pairing when `trait_order =
#'   "auto"`. Non-numeric phenotype columns, wrong orientation, or hidden ID
#'   columns treated as traits will cause errors or invalid model axes.
#' @param use_k_from Character scalar choosing the source of the kinship matrix
#'   \eqn{K}. `"X"` constructs \eqn{K} from genotype dosages `X` using
#'   [build_kinship_matrix()]. `"K"` uses the precomputed kinship object `K`.
#'   Choose `"X"` when SNP dosages are available and choose `"K"` when a trusted
#'   genomic relationship matrix has already been computed. The selected source
#'   must have the same subjects as `brain_csv`.
#' @param X Optional numeric matrix or data frame used when `use_k_from = "X"`.
#'   Subjects must be in rows and SNP markers in columns. A data frame may also
#'   contain one sample ID column and non-SNP columns, but SNP columns must be
#'   identifiable by `x_snp_cols`, `x_annotation`, `x_snp_pattern`, or numeric
#'   fallback. Do not let columns such as `PTID` enter the SNP set; that would
#'   either fail numeric validation or distort \eqn{K}.
#' @param K Optional numeric matrix or data frame used when `use_k_from = "K"`.
#'   It must represent a square subject-by-subject kinship matrix \eqn{K} with
#'   rows and columns in the same sample order. Sample IDs may be supplied by
#'   row names, matching row/column names, an ID column detected from a data
#'   frame, or `k_ids`. A non-square matrix or mismatched dimension stops the
#'   workflow; a wrong order without IDs can silently invalidate estimates if
#'   row-order alignment is allowed.
#' @param a_mat Optional numeric pair-level adjacency matrix \eqn{A} with
#'   dimension \eqn{J \times J}, where \eqn{J = C/2}. It encodes spatial
#'   neighborhood weights among left/right ROI pairs for the residual precision
#'   \eqn{D_A - \rho A}. When `NULL` (default), it is estimated from the aligned
#'   phenotype matrix using [build_a_from_y_halves()] and `k_knn`. Provide an
#'   anatomical or externally estimated adjacency when available.
#' @param brain_row_names Logical scalar. If `TRUE` (default), the first column
#'   of `brain_csv` is read as row names. This is appropriate for CSV files where
#'   the first column stores trait names under `traits_in_rows` or sample IDs
#'   under `samples_in_rows`. Set `FALSE` when the CSV has no row-name column;
#'   otherwise the first data column may be removed accidentally.
#' @param phenotype_layout Character scalar describing the orientation of
#'   `brain_csv`. `"auto"` (default) tries to infer orientation from left/right
#'   trait names in rows or columns and otherwise falls back to row names when
#'   present. `"traits_in_rows"` means the file already stores \eqn{Y} as traits
#'   by subjects. `"samples_in_rows"` means the file stores subjects by traits
#'   and will be transposed. Set this manually when auto-detection could confuse
#'   subject IDs with trait names.
#' @param phenotype_id_col Optional character scalar naming the sample ID column
#'   in `brain_csv` when `phenotype_layout = "samples_in_rows"`. That column is
#'   removed before numeric phenotype conversion and used to align samples with
#'   `X` or `K`. If omitted, row names are used when informative. Supplying a
#'   wrong column name stops; omitting IDs may force risky row-order alignment.
#' @param phenotype_sample_ids Optional character vector of sample IDs for
#'   `brain_csv` when phenotypes are stored as traits in rows and the file
#'   columns do not contain usable sample IDs. Its length must equal `ncol(Y)`.
#'   Use it to align phenotype columns with genotype rows or kinship rows when
#'   CSV column names are generic. Incorrect IDs stop alignment or subset to the
#'   wrong common samples.
#' @param x_ids Optional character vector of sample IDs for `X` when they are
#'   not available from row names or `x_id_col`. Length must equal `nrow(X)`.
#'   Use this to make genotype-to-phenotype alignment explicit and avoid
#'   row-order fallback.
#' @param k_ids Optional character vector of sample IDs for `K` when they are
#'   not available from row/column names or an ID column. Length should equal
#'   `nrow(K)`. Use this for precomputed kinship matrices stored without names.
#' @param x_id_col Optional character scalar naming the sample ID column inside
#'   `X`, for example `"PTID"`. The column is removed before SNP selection and
#'   used only for sample alignment. If `NULL`, the function tries to detect a
#'   unique non-numeric ID column. Passing the wrong column or leaving an ID
#'   column inside the SNP set can cause SNP selection errors or invalid
#'   kinship estimates.
#' @param x_snp_cols Optional SNP column selection for `X`, supplied as column
#'   names or column indices after removing `x_id_col`. Use this when SNP
#'   columns do not follow `x_snp_pattern` or when the data frame also contains
#'   numeric phenotype/covariate columns. A wrong selection changes \eqn{K};
#'   unknown names or non-numeric selected columns stop.
#' @param x_snp_pattern Character regular expression used to detect SNP columns
#'   when `x_snp_cols` is `NULL` and annotation matching is unavailable. The
#'   default `"^rs"` targets common rsID-style marker names. Change it for other
#'   marker naming conventions. If it misses all SNPs, the function falls back
#'   to all numeric columns, which can be unsafe when numeric covariates or brain
#'   traits are present.
#' @param x_annotation Optional data frame containing SNP annotation used to
#'   identify SNP columns in `X`. When `annotation_rsid_col` is present, columns
#'   of `X` whose names match annotation rsIDs are selected before applying
#'   `x_snp_pattern`. Use this to avoid treating non-SNP numeric columns as
#'   markers.
#' @param annotation_rsid_col Character scalar naming the SNP identifier column
#'   in `x_annotation`; default `"RSID"`. The values should match column names in
#'   `X`. A wrong name disables annotation-based SNP selection and may trigger
#'   pattern or numeric fallback.
#' @param trait_order Character scalar describing the row order of paired traits
#'   after `brain_csv` is oriented as \eqn{Y}. `"auto"` (default) detects
#'   left/right prefixes and reorders alternating pairs if needed. `"halves"`
#'   means rows are `(L1, ..., LJ, R1, ..., RJ)`. `"pairs"` means rows are
#'   `(L1, R1, L2, R2, ...)` and will be reordered to halves order before
#'   fitting. Supplying the wrong value pairs left/right regions incorrectly and
#'   changes both \eqn{A} and heritability summaries.
#' @param allow_row_order_alignment Logical scalar controlling a safety fallback
#'   when no usable sample IDs are available. If `TRUE` (default), the function
#'   may align phenotype columns to rows of `X` or rows/columns of `K` by their
#'   existing order, with a warning. Set `FALSE` for safer real analyses unless
#'   you have verified row order externally. Allowing row-order alignment with
#'   misordered data silently invalidates \eqn{K}-to-\eqn{Y} correspondence.
#' @param na_action Character scalar controlling missing phenotype values.
#'   `"error"` (default) stops if `brain_csv` contains `NA`, which is safest for
#'   model fitting. `"zero"` replaces missing values with zero before optional
#'   row-wise z-scoring; use it only when zero-imputation is scientifically
#'   justified, because it can alter trait covariance and adjacency estimates.
#' @param zscore_rows Logical scalar. If `TRUE` (default), each phenotype row is
#'   centered and scaled across aligned subjects before fitting, matching the
#'   standardized phenotype convention in the manuscript. Set `FALSE` only when
#'   the input \eqn{Y} has already been standardized or when raw-scale modeling
#'   is intentionally required.
#' @param k_knn Positive integer number of neighbors used when `a_mat = NULL`
#'   and \eqn{A} is estimated from the phenotype matrix. The default `2`
#'   produces a sparse pair-level graph. Increase for denser phenotype-derived
#'   spatial structure; decrease for very local structure.
#' @param rho_grid Numeric vector of candidate \eqn{\rho} values for the
#'   discrete variational posterior. The default `seq(0, 0.99, by = 0.01)`
#'   searches nonnegative spatial dependence. Use a narrower or denser grid for
#'   sensitivity analysis or when anatomical prior knowledge suggests a range.
#' @param max_iter Positive integer maximum number of VI iterations for the
#'   real-data fit. The default `10000` allows more iterations than examples or
#'   simulations. Reduce for smoke tests; increase if convergence diagnostics
#'   show the ELBO still changing at the limit.
#' @param tol Non-negative numeric ELBO convergence tolerance. The default
#'   `1e-4` stops when the ELBO stabilizes; use `0` to force running to
#'   `max_iter`. Larger values trade accuracy for speed.
#' @param verbose Logical scalar. If `TRUE` (default), print per-iteration ELBO
#'   progress from the VI optimizer. Set `FALSE` for batch or reproducible logs.
#' @param prefix Optional character output prefix passed to [save_vi_result()].
#'   `NULL` (default) skips writing files. If supplied, posterior summaries and
#'   the full result are written using this prefix, and existing files with the
#'   same names are overwritten.
#' @param make_plot Logical scalar. If `TRUE`, plot the ELBO trace after fitting
#'   using base graphics. Default `FALSE` avoids interactive graphics in batch
#'   workflows.
#' @param show_summary Logical scalar. If `TRUE`, print a concise text summary
#'   containing posterior `E[rho]`, broad-sense \eqn{H^2}, and narrow-sense
#'   \eqn{h_c^2} summaries. Default `FALSE` keeps the function quiet except for
#'   requested verbose optimizer output.
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
#' @inheritParams run_realdata_brain_hetero_halves
#' @param rdata_path Character scalar path to an `.RData` file containing one or
#'   more objects used by the real-data workflow. Depending on `use_k_from`, the
#'   file should contain either a genotype matrix/data frame \eqn{X}, a
#'   precomputed kinship matrix \eqn{K}, and optionally a SNP annotation data
#'   frame. Objects are loaded into an isolated environment; a missing file or
#'   missing requested object stops the workflow.
#' @param use_k_from Character scalar choosing how to obtain \eqn{K} from the
#'   `.RData` file. `"X"` uses `x_object` and constructs \eqn{K} from genotype
#'   dosages. `"K"` uses `k_object` as a precomputed kinship matrix. `"auto"`
#'   (default) first tries to detect a kinship-like object and otherwise falls
#'   back to a genotype-like object. Manual `"X"` or `"K"` is safer when the file
#'   contains multiple matrix-like objects.
#' @param x_object Optional character scalar object name in `rdata_path` to use
#'   as genotype input \eqn{X}. The object must be a numeric matrix or data frame
#'   with subjects in rows and SNP columns identifiable after removing any ID
#'   column. If `NULL`, the function attempts to detect a genotype-like object.
#'   Specify this manually when multiple genotype/covariate tables are present.
#' @param k_object Optional character scalar object name in `rdata_path` to use
#'   as a precomputed kinship matrix \eqn{K}. The object must be square numeric
#'   matrix-like data with subject order alignable to `brain_csv`. If `NULL`, the
#'   function attempts to detect a symmetric square numeric object.
#' @param genotype_object Optional backward-compatible alias for `x_object`. Use
#'   `x_object` in new code. If both are supplied, `x_object` takes precedence.
#'   A wrong object name stops with an object-not-found or no-genotype error.
#' @param annotation_object Optional character scalar object name in `rdata_path`
#'   containing SNP annotations, commonly a data frame with an `RSID` column.
#'   When supplied or auto-detected, it helps select SNP columns in `x_object`
#'   before pattern or numeric fallback. If the annotation object is named
#'   incorrectly, SNP selection may fall back to less safe rules.
#' @param genotype_rows Optional row indices used to subset the selected
#'   genotype object \eqn{X} or kinship object \eqn{K} before sample alignment.
#'   For `use_k_from = "X"`, rows subset subjects in the genotype table. For
#'   `"K"`, the same indices subset both rows and columns of the kinship matrix.
#'   Use only when the subset is known to match the phenotype subjects; an
#'   incorrect subset breaks \eqn{Y}/\eqn{K} alignment.
#' @param genotype_cols Optional backward-compatible SNP column selection passed
#'   to `x_snp_cols` when `x_snp_cols` is `NULL`. Values may be names or indices
#'   after removing `x_id_col`. Use `x_snp_cols` in new code. A wrong selection
#'   changes the constructed kinship matrix or stops if selected columns are not
#'   numeric.
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
