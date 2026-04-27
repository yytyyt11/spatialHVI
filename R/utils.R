# Internal helpers shared across the package.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

validate_numeric_matrix <- function(x, arg = "x") {
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be a numeric matrix.", arg), call. = FALSE)
  }
  if (nrow(x) == 0L || ncol(x) == 0L) {
    stop(sprintf("`%s` must have at least one row and one column.", arg), call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

validate_square_matrix <- function(x, arg = "x") {
  x <- validate_numeric_matrix(x, arg = arg)
  if (nrow(x) != ncol(x)) {
    stop(sprintf("`%s` must be a square matrix.", arg), call. = FALSE)
  }
  x
}

validate_even_trait_matrix <- function(x, arg = "y") {
  x <- validate_numeric_matrix(x, arg = arg)
  if (nrow(x) %% 2L != 0L) {
    stop(
      sprintf("`%s` must have an even number of rows in halves order.", arg),
      call. = FALSE
    )
  }
  x
}

validate_positive_count <- function(x, arg = "x") {
  if (length(x) != 1L || !is.finite(x) || x < 1 || x != as.integer(x)) {
    stop(sprintf("`%s` must be a positive integer.", arg), call. = FALSE)
  }
  as.integer(x)
}

make_degree_matrix <- function(a_mat, fallback = 1) {
  a_mat <- validate_square_matrix(a_mat, arg = "a_mat")
  if (nrow(a_mat) == 1L) {
    return(matrix(fallback, nrow = 1L, ncol = 1L))
  }
  diag(rowSums(a_mat), nrow = nrow(a_mat))
}

ensure_parent_dir <- function(path) {
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

safe_scale_cols <- function(x) {
  x <- validate_numeric_matrix(x, arg = "x")
  xs <- scale(x)
  xs[!is.finite(xs)] <- 0
  as.matrix(xs)
}

safe_scale_rows <- function(y) {
  y <- validate_numeric_matrix(y, arg = "y")
  ys <- t(apply(y, 1, scale))
  ys[!is.finite(ys)] <- 0
  as.matrix(ys)
}

#' Build a Kinship-Like Matrix from SNP Data
#'
#' Centers and scales SNP columns and returns \eqn{XX^T / L}, where \eqn{L} is
#' the number of SNP columns.
#'
#' @param snp_matrix Numeric SNP dosage matrix \eqn{X} with subjects in rows and
#'   SNPs or genetic markers in columns. In the model workflow this matrix is
#'   used only to construct the kinship matrix \eqn{K = XX^T / L}, where
#'   \eqn{L} is the number of SNP columns after centering and scaling each
#'   column. The row order becomes the subject order of \eqn{K} and must match
#'   the column order of the phenotype matrix \eqn{Y} used later in
#'   [run_vi_hetero_from_mats()] or [run_vi_homo_from_mats()]. Do not include
#'   sample ID columns such as `PTID` or non-SNP annotations; because the input
#'   is coerced to a numeric matrix, such columns can cause conversion errors or
#'   contaminate the genetic similarity calculation.
#' @param normalize_diag Logical scalar. If `FALSE` (default), return the raw
#'   centered-and-scaled relationship matrix \eqn{XX^T / L}. If `TRUE`, divide
#'   the matrix by its mean diagonal so the average self-relatedness is one.
#'   The VI fitting functions normalize `k_mat` internally, so the default is
#'   sufficient for package workflows; set `TRUE` when exporting or inspecting a
#'   standalone kinship matrix on the common mean-diagonal-one scale.
#'
#' @return A symmetric numeric matrix.
#' @export
#'
#' @examples
#' set.seed(1)
#' snp <- matrix(sample(0:2, 24, replace = TRUE), nrow = 6)
#' build_kinship_matrix(snp)
build_kinship_matrix <- function(snp_matrix, normalize_diag = FALSE) {
  snp_matrix <- validate_numeric_matrix(snp_matrix, arg = "snp_matrix")
  scaled <- safe_scale_cols(snp_matrix)
  snp_count <- ncol(scaled)
  if (snp_count < 1L) {
    stop("`snp_matrix` must contain at least one SNP column.", call. = FALSE)
  }

  kinship <- (scaled %*% t(scaled)) / snp_count
  kinship <- 0.5 * (kinship + t(kinship))

  if (normalize_diag) {
    diag_mean <- mean(diag(kinship))
    if (!is.finite(diag_mean) || diag_mean <= 0) {
      stop("Mean diagonal of the kinship matrix must be positive.", call. = FALSE)
    }
    kinship <- kinship / diag_mean
  }

  kinship
}

normalize_kinship_matrix <- function(k_mat) {
  k_mat <- validate_square_matrix(k_mat, arg = "k_mat")
  diag_mean <- mean(diag(k_mat))
  if (!is.finite(diag_mean) || diag_mean <= 0) {
    stop("Mean diagonal of `k_mat` must be positive.", call. = FALSE)
  }
  0.5 * (k_mat + t(k_mat)) / diag_mean
}

read_numeric_csv_matrix <- function(path, row_names = TRUE) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path), call. = FALSE)
  }

  df <- if (row_names) {
    utils::read.csv(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      row.names = 1
    )
  } else {
    utils::read.csv(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  keep <- vapply(df, is.numeric, logical(1))
  if (!all(keep)) {
    df <- df[, keep, drop = FALSE]
  }

  mat <- as.matrix(df)
  storage.mode(mat) <- "double"
  mat[!is.finite(mat)] <- 0
  mat
}

read_snp_matrix <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path), call. = FALSE)
  }

  snp_df <- tryCatch(
    utils::read.table(path, header = FALSE, sep = "", check.names = FALSE),
    error = function(e) NULL
  )

  if (is.null(snp_df) || ncol(snp_df) == 1L) {
    snp_df <- utils::read.csv(path, header = FALSE, check.names = FALSE)
  }

  snp_mat <- as.matrix(snp_df)
  suppressWarnings(storage.mode(snp_mat) <- "double")

  if (anyNA(snp_mat)) {
    numeric_df <- as.data.frame(
      lapply(snp_df, function(col) suppressWarnings(as.numeric(col))),
      stringsAsFactors = FALSE
    )
    keep_cols <- vapply(numeric_df, function(col) any(!is.na(col)), logical(1))
    numeric_df <- numeric_df[, keep_cols, drop = FALSE]
    snp_mat <- as.matrix(numeric_df)
    suppressWarnings(storage.mode(snp_mat) <- "double")

    if (nrow(snp_mat) > 1L && all(is.na(snp_mat[1, ]))) {
      snp_mat <- snp_mat[-1, , drop = FALSE]
    }
  }

  if (!is.numeric(snp_mat) || anyNA(snp_mat)) {
    stop("SNP input must be numeric.", call. = FALSE)
  }

  snp_mat
}

extract_snp_id <- function(path) {
  as.integer(sub("^SNP([0-9]+)\\.csv$", "\\1", basename(path)))
}

order_snp_files <- function(group_dir, pattern = "^SNP[0-9]+\\.csv$") {
  files <- list.files(group_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) {
    return(character(0))
  }

  ids <- vapply(files, extract_snp_id, integer(1))
  files[order(ids, basename(files))]
}

load_rdata_object <- function(path, object_name) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path), call. = FALSE)
  }

  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  if (!object_name %in% loaded) {
    stop(
      sprintf("Object `%s` was not found in `%s`.", object_name, path),
      call. = FALSE
    )
  }

  get(object_name, envir = env, inherits = FALSE)
}
