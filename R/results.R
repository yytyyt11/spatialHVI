validate_vi_result <- function(res) {
  needed <- c(
    "rho_grid",
    "rho_w",
    "E_sigma2",
    "nu",
    "S",
    "E_Sigma",
    "mu",
    "Sig_diag",
    "Sig_logdet",
    "Sig_trace_invK"
  )
  missing <- setdiff(needed, names(res))
  if (length(missing) > 0L) {
    stop(
      sprintf("The result object is missing required fields: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(res)
}

format_heritability_summary <- function(res, h) {
  sprintf(
    paste(
      "[Summary]",
      "  E[rho] = %.4f",
      "  Broad H^2 = %.4f",
      "  mean(narrow h^2) = %.4f (sd = %.4f)",
      sep = "\n"
    ),
    res$E_rho,
    h$H2,
    mean(h$h2),
    stats::sd(h$h2)
  )
}

#' Save Posterior Summaries from a VI Fit
#'
#' Writes posterior summaries derived from a fitted VI object to disk using a
#' shared file-name prefix.
#'
#' @param res Fitted VI result list returned by [run_vi_hetero_from_mats()] or
#'   [run_vi_homo_from_mats()]. The object must contain posterior summaries for
#'   the model quantities in the manuscript, including the grid and weights for
#'   \eqn{\rho}, genetic variance summaries `E_sigma2`, residual covariance
#'   summaries `E_Sigma`, posterior mean random effects `mu`, and covariance
#'   diagnostics. Passing a partial or modified list stops with a missing-field
#'   error.
#' @param prefix Character scalar output file prefix, optionally including a
#'   directory. The default `"vi_result"` writes files such as
#'   `vi_result_rho_posterior.csv`, `vi_result_sigma2_posterior.csv`, and
#'   `vi_result_full_result.rds` in the current working directory. Parent
#'   directories are created automatically. This is a prefix rather than an
#'   output directory; existing files with the same names are overwritten by the
#'   underlying CSV/RDS writers.
#'
#' @return An invisible named character vector of written file paths.
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' x <- matrix(rnorm(48), nrow = 8)
#' y <- matrix(rnorm(32), nrow = 4)
#' fit <- run_vi_hetero_from_mats(
#'   y = y,
#'   k_mat = build_kinship_matrix(x),
#'   rho_grid = seq(0, 0.2, by = 0.1),
#'   max_iter = 2,
#'   tol = 0,
#'   verbose = FALSE
#' )
#' save_vi_result(fit, prefix = file.path(tempdir(), "example_fit"))
#' }
save_vi_result <- function(res, prefix = "vi_result") {
  validate_vi_result(res)
  ensure_parent_dir(prefix)

  rho_path <- sprintf("%s_rho_posterior.csv", prefix)
  sigma2_path <- sprintf("%s_sigma2_posterior.csv", prefix)
  sigma_path <- sprintf("%s_Sigma_posterior.rds", prefix)
  mu_path <- sprintf("%s_mu_mean.csv", prefix)
  cov_path <- sprintf("%s_mu_cov_summaries.rds", prefix)
  full_path <- sprintf("%s_full_result.rds", prefix)

  rho_df <- data.frame(rho = res$rho_grid, weight = res$rho_w)
  utils::write.csv(rho_df, rho_path, row.names = FALSE)

  c_traits <- length(res$E_sigma2)
  a_vec <- as.numeric(res$a)
  d_vec <- as.numeric(res$d)
  sigma2_mean <- if (length(a_vec) == c_traits) {
    d_vec / pmax(a_vec - 1, 1e-8)
  } else {
    rep(d_vec / pmax(a_vec - 1, 1e-8), c_traits)
  }

  sigma_df <- data.frame(
    trait = seq_len(c_traits),
    a = rep(a_vec, length.out = c_traits),
    d = rep(d_vec, length.out = c_traits),
    mean = sigma2_mean,
    E_sigma2 = res$E_sigma2
  )
  utils::write.csv(sigma_df, sigma2_path, row.names = FALSE)

  saveRDS(
    list(nu = res$nu, S = res$S, E_Sigma = res$E_Sigma),
    sigma_path
  )

  utils::write.csv(as.matrix(res$mu), mu_path, row.names = FALSE)

  saveRDS(
    list(
      Sig_diag = res$Sig_diag,
      Sig_logdet = res$Sig_logdet,
      Sig_trace_invK = res$Sig_trace_invK
    ),
    cov_path
  )
  saveRDS(res, full_path)

  invisible(
    c(
      rho = rho_path,
      sigma2 = sigma2_path,
      sigma = sigma_path,
      mu = mu_path,
      covariance = cov_path,
      full = full_path
    )
  )
}

reconstruct_full_sig_list <- function(res) {
  if (is.null(res$Kpack) || is.null(res$final_alpha) || is.null(res$final_beta)) {
    stop("The result object does not contain enough information to reconstruct `Sig`.", call. = FALSE)
  }

  u_mat <- res$Kpack$U
  inv_lambda <- res$Kpack$inv_lambda

  lapply(seq_along(res$final_alpha), function(idx) {
    diag_terms <- 1 / (res$final_alpha[idx] * inv_lambda + res$final_beta[idx])
    u_mat %*% (diag(diag_terms, nrow = length(diag_terms))) %*% t(u_mat)
  })
}

#' Estimate Heritability from a VI Fit
#'
#' Computes broad-sense and trait-specific narrow-sense heritability from a fitted
#' halves-ordered VI model.
#'
#' @param res Fitted VI result list returned by [run_vi_hetero_from_mats()] or
#'   [run_vi_homo_from_mats()]. It supplies posterior estimates of
#'   \eqn{\sigma^2_{gc}}, \eqn{\Sigma}, \eqn{\rho}, and the fitted adjacency
#'   needed to compute trait-specific narrow-sense heritability \eqn{h_c^2} and
#'   the global principal-component heritability summary \eqn{H^2}. The result
#'   must correspond to an even number of halves-ordered traits.
#' @param a_mat Optional numeric pair-level adjacency matrix \eqn{A} with
#'   dimension \eqn{J \times J}, where \eqn{J = length(res$E_sigma2) / 2}. When
#'   `NULL` (default), the fitted adjacency `res$A` is used. Supply this
#'   explicitly when heritability should be computed with the known simulation or
#'   anatomical adjacency rather than the stored value. The matrix is
#'   normalized internally; a wrong pair order changes the residual variance
#'   decomposition, and a wrong dimension stops.
#' @param pairing Character scalar specifying how the \eqn{C} traits map to
#'   left/right ROI pairs. Currently only `"halves"` is supported, meaning the
#'   fitted result uses rows `(L1, ..., LJ, R1, ..., RJ)`. Other layouts should
#'   be reordered before fitting; unsupported values are rejected by
#'   [base::match.arg()].
#'
#' @return A list with elements `h2`, `H2`, and `pieces`.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(48), nrow = 8)
#' y <- matrix(rnorm(32), nrow = 4)
#' fit <- run_vi_hetero_from_mats(
#'   y = y,
#'   k_mat = build_kinship_matrix(x),
#'   rho_grid = seq(0, 0.2, by = 0.1),
#'   max_iter = 2,
#'   tol = 0,
#'   verbose = FALSE
#' )
#' heritability_from_vi(fit)$H2
heritability_from_vi <- function(res, a_mat = NULL, pairing = c("halves")) {
  pairing <- match.arg(pairing)
  validate_vi_result(res)

  c_traits <- length(res$E_sigma2)
  if (c_traits %% 2L != 0L) {
    stop("The result object must correspond to an even number of traits.", call. = FALSE)
  }
  j_pairs <- c_traits / 2L

  a_use <- if (is.null(a_mat)) res$A else normalize_a_spectral(a_mat)
  a_use <- validate_square_matrix(a_use, arg = "a_mat")
  if (nrow(a_use) != j_pairs) {
    stop("`a_mat` must have dimension equal to `length(res$E_sigma2) / 2`.", call. = FALSE)
  }

  diag(a_use) <- 0
  d_a <- make_degree_matrix(a_use)
  r_bar <- weighted_precision_matrix(res$rho_grid, res$rho_w, d_a, a_use)
  r_inv <- inv_spd(r_bar)

  e_sigma2 <- as.numeric(res$E_sigma2)
  e_sigma <- res$E_Sigma
  var_e <- numeric(c_traits)
  diag_r_inv <- diag(r_inv)
  var_e[seq_len(j_pairs)] <- diag_r_inv * e_sigma[1, 1]
  var_e[(j_pairs + 1L):c_traits] <- diag_r_inv * e_sigma[2, 2]

  h2_narrow <- e_sigma2 / pmax(e_sigma2 + var_e, 1e-12)

  sigma_g <- diag(e_sigma2, nrow = c_traits, ncol = c_traits)
  sigma_e <- kronecker(e_sigma, r_inv)
  chol_sigma_e <- chol(0.5 * (sigma_e + t(sigma_e)) + diag(1e-10, c_traits))
  sigma_e_inv_half <- backsolve(chol_sigma_e, diag(c_traits))
  m_mat <- sigma_e_inv_half %*% sigma_g %*% t(sigma_e_inv_half)
  lambda_max <- max(eigen(0.5 * (m_mat + t(m_mat)), symmetric = TRUE, only.values = TRUE)$values)
  lambda_max <- max(lambda_max, 0)

  list(
    h2 = h2_narrow,
    H2 = lambda_max / (1 + lambda_max),
    pieces = list(
      pairing = pairing,
      E_sigma2 = e_sigma2,
      var_e = var_e,
      R_bar = r_bar,
      R_inv = r_inv,
      Sigma_g = sigma_g,
      Sigma_e = sigma_e
    )
  )
}
