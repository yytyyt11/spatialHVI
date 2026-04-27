# Internal adjacency and precision helpers.

normalize_a_spectral <- function(a_mat, eps = 1e-12) {
  a_mat <- validate_square_matrix(a_mat, arg = "a_mat")
  a_mat <- 0.5 * (a_mat + t(a_mat))
  diag(a_mat) <- 0

  deg <- rowSums(a_mat)
  d_half_inv <- diag(1 / sqrt(pmax(deg, eps)), nrow = nrow(a_mat))
  scaled <- d_half_inv %*% a_mat %*% d_half_inv
  lam_max <- max(eigen(scaled, symmetric = TRUE, only.values = TRUE)$values)

  if (!is.finite(lam_max) || lam_max <= 0) {
    lam_max <- 1
  }

  a_mat / lam_max
}

ensure_spd_rho <- function(d_a, a_mat, rho_init = 0.8) {
  d_a <- validate_square_matrix(d_a, arg = "d_a")
  a_mat <- validate_square_matrix(a_mat, arg = "a_mat")

  spd_ok <- function(rho) {
    ev <- eigen(d_a - rho * a_mat, symmetric = TRUE, only.values = TRUE)$values
    is.finite(min(ev)) && min(ev) > 0
  }

  rho <- min(as.numeric(rho_init), 0.95)
  if (spd_ok(rho)) {
    return(rho)
  }

  for (factor in c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)) {
    candidate <- rho * factor
    if (spd_ok(candidate)) {
      return(candidate)
    }
  }

  stop("Could not find a positive-definite rho for the supplied adjacency matrix.", call. = FALSE)
}

expand_sigma2_pairs_halves <- function(sigma2_pairs, c_traits) {
  c_traits <- validate_positive_count(c_traits, arg = "c_traits")
  if (c_traits %% 2L != 0L) {
    stop("`c_traits` must be even.", call. = FALSE)
  }

  j_pairs <- c_traits / 2L
  if (length(sigma2_pairs) != j_pairs) {
    warning(
      sprintf(
        "`sigma2_pairs` has length %d but expected %d; values were recycled.",
        length(sigma2_pairs),
        j_pairs
      ),
      call. = FALSE
    )
    sigma2_pairs <- rep(sigma2_pairs, length.out = j_pairs)
  }

  as.numeric(c(sigma2_pairs, sigma2_pairs))
}

load_a_for_vi <- function(a_path, c_traits, renormalize = TRUE) {
  c_traits <- validate_positive_count(c_traits, arg = "c_traits")
  j_pairs <- c_traits / 2L

  a_mat <- as.matrix(utils::read.csv(a_path, header = FALSE))
  storage.mode(a_mat) <- "double"
  if (nrow(a_mat) != j_pairs || ncol(a_mat) != j_pairs) {
    stop("Adjacency matrix dimensions do not match the number of trait pairs.", call. = FALSE)
  }

  a_mat <- 0.5 * (a_mat + t(a_mat))
  diag(a_mat) <- 0

  if (renormalize) {
    a_mat <- normalize_a_spectral(a_mat)
  }

  a_mat
}

make_uniform_rho_prior_valid <- function(a_mat, d_a, rho_grid) {
  rho_grid <- as.numeric(rho_grid)
  valid <- vapply(
    rho_grid,
    function(rho) {
      candidate <- d_a - rho * a_mat
      is.null(tryCatch(chol(0.5 * (candidate + t(candidate))), error = function(e) NULL)) == FALSE
    },
    logical(1)
  )

  if (!any(valid)) {
    stop("No valid rho values were found for the supplied adjacency matrix.", call. = FALSE)
  }

  prior <- rep(0, length(rho_grid))
  prior[valid] <- 1 / sum(valid)
  prior
}

#' Create a k-Nearest-Neighbor Adjacency Matrix
#'
#' Generates a symmetric adjacency matrix from random two-dimensional coordinates.
#'
#' @param j_pairs Positive integer number of paired regions of interest
#'   \eqn{J = C/2}. The returned adjacency matrix has one row and column per
#'   left/right ROI pair, not one row per individual left or right trait. This
#'   must match the pair dimension used by halves-ordered phenotype matrices
#'   with \eqn{C = 2J} rows.
#' @param k Positive integer number of nearest neighbors retained for each ROI
#'   pair when constructing the pair-level adjacency \eqn{A}. The default `2`
#'   creates a sparse local graph. Increase `k` for denser spatial dependence or
#'   decrease it for very local dependence; values larger than `j_pairs - 1` are
#'   effectively capped by the available neighbors.
#' @param sigma Positive numeric bandwidth for converting squared distances
#'   between simulated two-dimensional coordinates into edge weights through a
#'   Gaussian kernel. The default `1` gives moderate decay. Smaller values make
#'   non-nearest edges weaker; larger values make retained edges more similar.
#'   Non-positive or non-finite values can produce invalid adjacency weights.
#' @param seed Optional integer random seed used only for generating the
#'   temporary two-dimensional coordinates. The default `NULL` leaves the
#'   current RNG state unchanged. Set a seed when a simulation group must use a
#'   reproducible adjacency matrix.
#'
#' @return A symmetric `j_pairs` by `j_pairs` adjacency matrix.
#' @export
#'
#' @examples
#' make_a_knn(j_pairs = 3, k = 1, seed = 1)
make_a_knn <- function(j_pairs, k = 2, sigma = 1, seed = NULL) {
  j_pairs <- validate_positive_count(j_pairs, arg = "j_pairs")
  k <- validate_positive_count(max(1L, as.integer(k)), arg = "k")

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (j_pairs == 1L) {
    return(matrix(0, nrow = 1L, ncol = 1L))
  }

  coords <- matrix(stats::rnorm(j_pairs * 2L), ncol = 2L)
  distance_sq <- as.matrix(stats::dist(coords))^2
  weights <- exp(-distance_sq / (2 * sigma^2))
  diag(weights) <- 0
  weights <- 0.5 * (weights + t(weights))

  k_eff <- min(k, max(1L, j_pairs - 1L))
  a_mat <- matrix(0, nrow = j_pairs, ncol = j_pairs)
  for (idx in seq_len(j_pairs)) {
    nn_idx <- order(weights[idx, ], decreasing = TRUE)[seq_len(k_eff)]
    a_mat[idx, nn_idx] <- weights[idx, nn_idx]
  }

  a_mat <- pmax(a_mat, t(a_mat))
  diag(a_mat) <- 0

  zero_degree <- which(rowSums(a_mat) == 0)
  if (length(zero_degree) > 0L) {
    for (idx in zero_degree) {
      best <- which.max(weights[idx, ])
      if (best != idx && weights[idx, best] > 0) {
        a_mat[idx, best] <- weights[idx, best]
        a_mat[best, idx] <- weights[idx, best]
      }
    }
  }

  a_mat
}

#' Build an Adjacency Matrix from Halves-Ordered Traits
#'
#' Constructs a k-nearest-neighbor adjacency matrix from paired left and right
#' traits. By default the function auto-detects whether traits are already in
#' halves order `(L1, ..., LJ, R1, ..., RJ)` or in alternating pair order
#' `(L1, R1, L2, R2, ...)`.
#'
#' @param y Numeric phenotype matrix \eqn{Y} with paired traits in rows and
#'   subjects in columns. The matrix must have an even number of rows
#'   \eqn{C = 2J}. This function uses the similarity of left-hemisphere rows and
#'   right-hemisphere rows across the same subjects to construct the pair-level
#'   adjacency \eqn{A}. Subject columns must therefore be in a consistent order
#'   within all rows, although no kinship matrix is used here. If traits are in
#'   columns instead of rows, the inferred adjacency will be meaningless or the
#'   even-row check will fail.
#' @param k Positive integer number of neighbors retained for each ROI pair in
#'   the estimated adjacency graph. The default `2` produces a sparse
#'   correlation-based graph. Increase it when the expected residual spatial
#'   dependence is broader; decrease it when only the strongest local
#'   relationships should be retained.
#' @param trait_order Character scalar describing how paired trait rows are
#'   arranged before building \eqn{A}. `"auto"` (default) inspects row names with
#'   left/right prefixes and treats detected alternating names as `"pairs"`;
#'   otherwise it assumes `"halves"`. `"halves"` means rows are already
#'   `(L1, ..., LJ, R1, ..., RJ)`. `"pairs"` means rows are alternating
#'   `(L1, R1, L2, R2, ...)` and are internally reordered to halves order for
#'   adjacency construction. Supplying the wrong order connects the wrong ROI
#'   pairs and can bias the spatial residual precision \eqn{D_A - \rho A}.
#'
#' @return A normalized `J` by `J` adjacency matrix.
#' @export
#'
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(24), nrow = 4)
#' rownames(y) <- c("Left_A", "Right_A", "Left_B", "Right_B")
#' build_a_from_y_halves(y, k = 1, trait_order = "auto")
build_a_from_y_halves <- function(y, k = 2, trait_order = c("auto", "halves", "pairs")) {
  y <- validate_even_trait_matrix(y, arg = "y")
  k <- validate_positive_count(max(1L, as.integer(k)), arg = "k")
  trait_order <- match.arg(trait_order)

  y <- reorder_traits_to_halves(y, trait_order = trait_order)$y

  c_traits <- nrow(y)
  j_pairs <- c_traits / 2L
  if (j_pairs == 1L) {
    return(matrix(0, nrow = 1L, ncol = 1L))
  }

  y_left <- y[seq_len(j_pairs), , drop = FALSE]
  y_right <- y[(j_pairs + 1L):c_traits, , drop = FALSE]

  w_left <- suppressWarnings(stats::cor(t(y_left)))
  w_right <- suppressWarnings(stats::cor(t(y_right)))
  w_left[!is.finite(w_left)] <- 0
  w_right[!is.finite(w_right)] <- 0

  weights <- abs(0.5 * (w_left + w_right))
  diag(weights) <- 0

  a_mat <- matrix(0, nrow = j_pairs, ncol = j_pairs)
  for (idx in seq_len(j_pairs)) {
    candidates <- setdiff(seq_len(j_pairs), idx)
    ordered <- candidates[order(weights[idx, candidates], decreasing = TRUE)]
    nn_idx <- ordered[seq_len(min(k, length(ordered)))]
    a_mat[idx, nn_idx] <- weights[idx, nn_idx]
  }

  a_mat <- pmax(a_mat, t(a_mat))
  diag(a_mat) <- 0

  normalize_a_spectral(a_mat)
}
