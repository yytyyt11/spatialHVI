#' Model Priors Used by the Manuscript Interface
#'
#' Creates inverse-gamma genetic-variance and inverse-Wishart bilateral
#' residual-covariance priors. These are model priors, separate from initial
#' variational distributions.
#' @param n_traits Positive even number of halves-ordered traits.
#' @param ig_shape Inverse-gamma shape, greater than one.
#' @param ig_mean Prior mean of each genetic variance.
#' @param iw_df Inverse-Wishart degrees of freedom, greater than three.
#' @param iw_mean Positive definite 2 by 2 prior mean of the bilateral covariance.
#' @return A list with `a0`, `d0`, `v0`, `S0` and prior summaries.
#' @export
spatial_hvi_prior <- function(n_traits, ig_shape = 2.1, ig_mean = 0.30,
                              iw_df = 8, iw_mean = diag(c(0.40, 0.40))) {
  .manuscript_engine$make_vi_prior(n_traits, iw_df, iw_mean, ig_shape, ig_mean)
}

#' Fit the Manuscript Spatial MFVI Model
#'
#' Uses the numerical updates and initialization from the manuscript's frozen
#' experiment implementation. Point estimation is kept separate from LRVB and
#' interval calibration. The older `run_vi_*` functions remain available with
#' their original initialization and summary conventions.
#'
#' @param y Numeric C by n phenotype matrix, in row order
#'   `(L1, ..., LJ, R1, ..., RJ)`, where C = 2J. Covariate adjustment and any
#'   phenotype standardization must already have been performed.
#' @param k_mat Symmetric positive semidefinite n by n GRM. Its row and column
#'   order must match the columns of `y`. Named inputs are checked for alignment.
#' @param a_mat Symmetric nonnegative J by J adjacency with zero diagonal and
#'   positive degrees. A supplied graph is used as given. If NULL, the original
#'   pair-mean correlation k-nearest-neighbor graph (k = 2) is estimated from `y`.
#' @param rho_grid Original discrete spatial-coupling support. Use `0` for the
#'   fixed-rho ablation, whose pair precision is the degree matrix, not identity.
#' @param prior A list produced by [spatial_hvi_prior()], or NULL for its defaults.
#' @param normalize_k Divide the GRM by its mean diagonal before regularization.
#' @param k_ridge Diagonal GRM regularization, applied once after normalization.
#' @param max_iter Maximum number of coordinate-ascent iterations.
#' @param tol Positive absolute ELBO-change stopping tolerance.
#' @param init_h2 Initial genetic-variance allocation.
#' @param init_ig_shape Initial inverse-gamma shape.
#' @param init_iw_df Initial inverse-Wishart degrees of freedom.
#' @param verbose Print iteration progress.
#' @return A `spatialhvi_fit` list containing the original MFVI estimates,
#'   convergence information, prior, graph, GRM eigendecomposition and input
#'   phenotypes. The latter are retained for LRVB; fitted objects from restricted
#'   data must be handled with the same access restrictions as their inputs.
#' @export
fit_spatial_hvi <- function(y, k_mat, a_mat = NULL,
                            rho_grid = seq(0, 0.99, by = 0.01), prior = NULL,
                            normalize_k = TRUE, k_ridge = 1e-6,
                            max_iter = 5000L, tol = 1e-4,
                            init_h2 = 0.15, init_ig_shape = 3, init_iw_df = 8,
                            verbose = FALSE) {
  e <- .manuscript_engine
  e$.assert_matrix(y, "y")
  if (nrow(y) < 4L || nrow(y) %% 2L || ncol(y) < 2L) {
    stop("y needs at least two bilateral pairs and two subjects.", call. = FALSE)
  }
  n <- ncol(y); j <- nrow(y) / 2L
  e$.assert_matrix(k_mat, "k_mat", n, n, symmetric = TRUE)
  e$.assert_psd(k_mat, "k_mat")
  if (!is.null(rownames(k_mat)) && !is.null(colnames(k_mat)) &&
      !identical(rownames(k_mat), colnames(k_mat))) {
    stop("GRM row and column subject names differ.", call. = FALSE)
  }
  for (ids in list(rownames(k_mat), colnames(k_mat))) {
    if (!is.null(ids) && !is.null(colnames(y)) && !identical(ids, colnames(y))) {
      stop("Phenotype and GRM subject order differs.", call. = FALSE)
    }
  }
  if (length(normalize_k) != 1L || is.na(normalize_k) || !is.logical(normalize_k))
    stop("normalize_k must be TRUE or FALSE.", call. = FALSE)
  if (length(k_ridge) != 1L || !is.finite(k_ridge) || k_ridge < 0)
    stop("k_ridge must be nonnegative.", call. = FALSE)
  k_scale <- if (normalize_k) mean(diag(k_mat)) else 1
  if (!is.finite(k_scale) || k_scale <= 0) stop("GRM mean diagonal must be positive.")
  kp <- e$make_K_eig_pack(k_mat / k_scale, eps = k_ridge)
  if (any(kp$lam <= 0)) stop("Regularized GRM must be positive definite.")
  graph_origin <- if (is.null(a_mat)) "estimated_from_y" else "supplied"
  if (is.null(a_mat)) a_mat <- e$build_A_from_Y_halves(y, k = 2)
  e$.assert_matrix(a_mat, "a_mat", j, j, symmetric = TRUE)
  if (any(a_mat < 0) || any(diag(a_mat) != 0) || any(rowSums(a_mat) <= 0))
    stop("a_mat needs nonnegative weights, zero diagonal and positive degrees.")
  if (!is.numeric(rho_grid) || !length(rho_grid) || any(!is.finite(rho_grid)) ||
      any(rho_grid < 0 | rho_grid >= 1) || anyDuplicated(rho_grid))
    stop("rho_grid must contain distinct finite values in [0,1).")
  if (!is.null(prior)) {
    for (key in c("a0", "d0"))
      if (!length(prior[[key]]) %in% c(1L, nrow(y)))
        stop("Prior a0 and d0 must be scalar or have one entry per trait.")
  }
  fit <- e$run_vi_custom(y, kp, a_mat, diag(rowSums(a_mat)), rho_grid = rho_grid,
    prior = prior, max_iter = max_iter, tol = tol, verbose = verbose,
    init_h2 = init_h2, init_ig_shape = init_ig_shape, init_iw_df = init_iw_df)
  fit$input_y <- y
  fit$input_settings <- list(k_normalization = k_scale, k_ridge = k_ridge,
    graph_origin = graph_origin, max_iter = max_iter, tol = tol)
  fit$trait_names <- if (is.null(rownames(y))) paste0("trait_", seq_len(nrow(y))) else rownames(y)
  class(fit) <- c("spatialhvi_fit", class(fit))
  fit
}

manuscript_heritability_summary <- function(fit) {
  qi <- .manuscript_engine$weighted_Rinv_from_rho(fit$A, fit$D_A,
                                                fit$rho_grid, fit$rho_w)
  j <- nrow(fit$A)
  g <- as.numeric(fit$E_sigma2)
  residual <- rep(diag(qi), 2) * rep(diag(fit$E_Sigma), each = j)
  h <- g / pmax(g + residual, 1e-12)
  G <- diag(g)
  E <- kronecker(fit$E_Sigma, qi)
  Uinv <- backsolve(chol((E + t(E)) / 2), diag(length(g)))
  M <- t(Uinv) %*% G %*% Uinv
  lambda <- max(eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  list(h2 = h, H2 = lambda / (1 + lambda),
    bilateral_mean = (h[seq_len(j)] + h[j + seq_len(j)]) / 2,
    pieces = list(pairing = "halves", E_sigma2 = g, var_e = residual,
      R_bar = fit$D_A - fit$E_rho * fit$A, R_inv = qi, Sigma_g = G, Sigma_e = E,
      residual_summary = "mean_inverse"))
}
