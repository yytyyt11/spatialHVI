## ============================  vi_halves_unified_fast.R  ============================ ##
## Fast exact version for HALVES ordering.
##
## Compared with the dense prototype:
##   - SAME model / SAME updates / SAME targets
##   - Faster q(b_c): use ONE eigendecomposition of K, then exact closed forms
##   - No explicit dense inverse of alpha*K^{-1} + beta*I per trait
##   - By default, full q(b_c) covariance matrices are NOT stored; only exact summaries
##
## Main user-visible difference:
##   - res$Sig is not a list of full n x n covariance matrices by default
##   - instead, store res$Sig_diag, res$Sig_logdet, res$Sig_trace_invK
##   - posterior means / rho / H2 / narrow h2 should match the dense version
##     up to numerical tolerance

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(stats)
})

## --------------------------- Utilities --------------------------- ##
safe_scale_cols <- function(X) {
  Xs <- scale(X)
  Xs[!is.finite(Xs)] <- 0
  as.matrix(Xs)
}

safe_scale_rows <- function(Y) {
  Ys <- t(apply(Y, 1, scale))
  Ys[!is.finite(Ys)] <- 0
  as.matrix(Ys)
}

read_numeric_csv_matrix <- function(path, row_names = TRUE) {
  if (row_names) {
    df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)
  } else {
    df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  is_num <- vapply(df, is.numeric, logical(1))
  if (!all(is_num)) {
    df <- df[, is_num, drop = FALSE]
  }
  M <- as.matrix(df)
  storage.mode(M) <- "double"
  M[!is.finite(M)] <- 0
  M
}

inv_spd <- function(M, jitter = 1e-8, max_tries = 6) {
  p <- nrow(M)
  for (tt in 0:max_tries) {
    Mt <- 0.5 * (M + t(M)) + diag(jitter * (10^tt), p)
    co <- tryCatch(chol(Mt), error = function(e) NULL)
    if (!is.null(co)) return(chol2inv(co))
  }
  stop("inv_spd: matrix not SPD even after jitter.")
}

logdet_spd <- function(M, jitter = 1e-8, max_tries = 6) {
  p <- nrow(M)
  for (tt in 0:max_tries) {
    Mt <- 0.5 * (M + t(M)) + diag(jitter * (10^tt), p)
    co <- tryCatch(chol(Mt), error = function(e) NULL)
    if (!is.null(co)) return(2 * sum(log(diag(co))))
  }
  stop("logdet_spd: matrix not SPD even after jitter.")
}

is_spd <- function(M) {
  Mt <- 0.5 * (M + t(M))
  !is.null(tryCatch(chol(Mt), error = function(e) NULL))
}

make_K_eig_pack <- function(K, eps = 1e-6) {
  Kj <- 0.5 * (K + t(K)) + diag(eps, nrow(K))
  eg <- eigen(Kj, symmetric = TRUE)
  lam <- as.numeric(eg$values)
  U   <- eg$vectors
  inv_lam <- 1 / lam
  U2 <- U * U
  list(
    K = Kj,
    U = U,
    U2 = U2,
    lam = lam,
    inv_lam = inv_lam,
    logdetK = sum(log(lam))
  )
}

normalize_A_spectral <- function(A) {
  A <- as.matrix(A)
  J <- nrow(A)
  stopifnot(nrow(A) == ncol(A))
  A <- 0.5 * (A + t(A))
  diag(A) <- 0

  eps <- 1e-12
  deg <- rowSums(A)
  Dhi <- diag(1 / sqrt(pmax(deg, eps)), J)
  B   <- Dhi %*% A %*% Dhi
  lam_max <- max(eigen(B, symmetric = TRUE, only.values = TRUE)$values)

  if (!is.finite(lam_max) || lam_max <= 0) lam_max <- 1
  A / lam_max
}

load_A_for_vi <- function(A_path, C, renormalize = TRUE) {
  J <- C / 2
  A <- as.matrix(read.csv(A_path, header = FALSE))
  stopifnot(nrow(A) == J, ncol(A) == J)
  A <- 0.5 * (A + t(A))
  diag(A) <- 0
  if (renormalize) A <- normalize_A_spectral(A)
  A
}

make_uniform_rho_prior_valid <- function(A, D_A, rho_grid) {
  valid <- vapply(rho_grid, function(rho) {
    is_spd(D_A - rho * A)
  }, logical(1))

  if (!any(valid)) {
    stop("No valid SPD rho grid points found. Check A scaling or rho_grid.")
  }

  prior <- rep(0, length(rho_grid))
  prior[valid] <- 1 / sum(valid)
  prior
}

log_multigamma <- function(a, p) {
  (p * (p - 1) / 4) * log(pi) + sum(lgamma(a + (1 - (1:p)) / 2))
}

E_logdet_Sigma_IW <- function(nu, S, p = 2) {
  logdet_spd(S) - sum(digamma((nu + 1 - (1:p)) / 2)) - p * log(2)
}

E_Sinv_IW <- function(nu, S) {
  nu * inv_spd(S)
}

entropy_inv_gamma <- function(a, b) {
  a + log(b) + lgamma(a) - (1 + a) * digamma(a)
}
E_log_sigma2_IG <- function(a, b) {
  log(b) - digamma(a)
}
E_inv_sigma2_IG <- function(a, b) {
  a / b
}

## --------------------------- Build A from Y (halves) --------------------------- ##
build_A_from_Y_halves_knn <- function(Y, k = 2) {
  C <- nrow(Y)
  J <- C / 2
  stopifnot(C %% 2 == 0)

  YL <- Y[1:J, , drop = FALSE]
  YR <- Y[(J + 1):(2 * J), , drop = FALSE]

  WL <- suppressWarnings(cor(t(YL)))
  WR <- suppressWarnings(cor(t(YR)))
  WL[!is.finite(WL)] <- 0
  WR[!is.finite(WR)] <- 0

  W <- abs(0.5 * (WL + WR))
  diag(W) <- 0

  A <- matrix(0, J, J)
  for (i in 1:J) {
    cand <- setdiff(seq_len(J), i)
    ord  <- cand[order(W[i, cand], decreasing = TRUE)]
    nn   <- ord[1:min(k, J - 1)]
    A[i, nn] <- W[i, nn]
  }
  A <- pmax(A, t(A))
  diag(A) <- 0

  normalize_A_spectral(A)
}

build_A_from_Y_halves <- function(Y, k = 2) {
  build_A_from_Y_halves_knn(Y, k = k)
}

## --------- halves reshaping helpers --------- ##
vec_to_2xJ_halves <- function(v) {
  C <- length(v)
  J <- C / 2
  rbind(v[1:J], v[(J + 1):(2 * J)])
}

mat_to_mu_list <- function(m_mat) {
  lapply(seq_len(nrow(m_mat)), function(i) matrix(m_mat[i, ], ncol = 1))
}

mu_to_mat <- function(mu_obj) {
  if (is.matrix(mu_obj)) return(mu_obj)
  do.call(rbind, lapply(mu_obj, function(x) as.numeric(x[, 1])))
}

## ========= Helper: save VI results ========= ##
save_vi_result <- function(res, prefix = "vi_result") {
  rho_df <- data.frame(
    rho    = res$rho_grid,
    weight = res$rho_w
  )
  write.csv(rho_df, sprintf("%s_rho_posterior.csv", prefix), row.names = FALSE)

  C <- length(res$E_sigma2)
  a_vec <- as.numeric(res$a)
  d_vec <- as.numeric(res$d)

  if (length(a_vec) == C) {
    sigma_df <- data.frame(
      trait    = 1:C,
      a        = a_vec,
      d        = d_vec,
      mean     = d_vec / pmax(a_vec - 1, 1e-8),
      E_sigma2 = res$E_sigma2
    )
  } else {
    sigma_df <- data.frame(
      trait    = 1:C,
      a        = rep(a_vec, C),
      d        = rep(d_vec, C),
      E_sigma2 = res$E_sigma2
    )
  }
  write.csv(sigma_df, sprintf("%s_sigma2_posterior.csv", prefix), row.names = FALSE)

  Sigma_list <- list(
    nu      = res$nu,
    S       = res$S,
    E_Sigma = res$E_Sigma
  )
  saveRDS(Sigma_list, sprintf("%s_Sigma_posterior.rds", prefix))

  mu_mat <- mu_to_mat(res$mu)
  write.csv(mu_mat, sprintf("%s_mu_mean.csv", prefix), row.names = FALSE)

  ## Fast version stores posterior covariance summaries, not full matrices
  Sig_summary <- list(
    Sig_diag = res$Sig_diag,
    Sig_logdet = res$Sig_logdet,
    Sig_trace_invK = res$Sig_trace_invK
  )
  saveRDS(Sig_summary, sprintf("%s_mu_cov_summaries.rds", prefix))

  saveRDS(res, sprintf("%s_full_result.rds", prefix))
}

## Optional: reconstruct full final q(b_c) covariance matrices (expensive, on demand)
reconstruct_full_Sig_list <- function(res) {
  if (is.null(res$Kpack) || is.null(res$final_alpha) || is.null(res$final_beta)) {
    stop("Result object does not contain enough information to reconstruct full Sig.")
  }
  U <- res$Kpack$U
  inv_lam <- res$Kpack$inv_lam
  lapply(seq_along(res$final_alpha), function(c) {
    d <- 1 / (res$final_alpha[c] * inv_lam + res$final_beta[c])
    U %*% (diag(d, length(d))) %*% t(U)
  })
}

## --------- S sufficient stats, HALVES pairing --------- ##
S_sufficient_halves_fast <- function(y, m_mat, V_mat, R) {
  C <- nrow(y)
  M <- ncol(y)
  J <- C / 2

  acc <- matrix(0, 2, 2)
  for (l in 1:M) {
    y2 <- vec_to_2xJ_halves(y[, l])
    m2 <- vec_to_2xJ_halves(m_mat[, l])
    Rm <- y2 - m2

    acc <- acc + Rm %*% R %*% t(Rm)

    for (j in 1:J) {
      acc <- acc + R[j, j] * diag(c(V_mat[j, l], V_mat[j + J, l]), 2, 2)
    }
  }
  acc
}

## --------------------------- Fast q(b) update --------------------------- ##
## Exact because:
##   prec_c = alpha_c * K^{-1} + beta_c * I
## and K has been eigendecomposed once.
update_q_b_fast <- function(y, Kpack, R, H, m_mat, alpha_vec) {
  C <- nrow(y)
  n <- ncol(y)

  G  <- kronecker(H, R)
  Gy <- G %*% y

  U <- Kpack$U
  U2 <- Kpack$U2
  inv_lam <- Kpack$inv_lam

  V_mat <- matrix(0, C, n)
  logdet_Sig <- numeric(C)
  trace_invK_Sig <- numeric(C)
  quad_mean <- numeric(C)
  quad_total <- numeric(C)
  beta_vec <- diag(G)

  for (c in 1:C) {
    g_row <- G[c, ]
    beta  <- beta_vec[c]

    rhs <- as.numeric(Gy[c, ] - (g_row %*% m_mat - beta * m_mat[c, ]))
    u_rhs <- as.numeric(crossprod(U, rhs))

    denom <- alpha_vec[c] * inv_lam + beta
    d <- 1 / denom

    u_m <- d * u_rhs
    m_c <- as.numeric(U %*% u_m)

    m_mat[c, ] <- m_c
    V_mat[c, ] <- as.numeric(U2 %*% d)

    logdet_Sig[c] <- sum(log(d))
    trace_invK_Sig[c] <- sum(inv_lam * d)
    quad_mean[c] <- sum(inv_lam * (u_m^2))
    quad_total[c] <- quad_mean[c] + trace_invK_Sig[c]
  }

  list(
    m_mat = m_mat,
    V_mat = V_mat,
    logdet_Sig = logdet_Sig,
    trace_invK_Sig = trace_invK_Sig,
    quad_mean = quad_mean,
    quad_total = quad_total,
    beta_vec = beta_vec
  )
}

## --------------------------- Sigma & rho updates --------------------------- ##
update_q_Sigma_halves_fast <- function(y, m_mat, V_mat, R, v0, S0) {
  C <- nrow(y)
  J <- C / 2
  M <- ncol(y)

  nu_t <- v0 + M * J
  S_t  <- S0 + S_sufficient_halves_fast(y, m_mat, V_mat, R)
  S_t  <- 0.5 * (S_t + t(S_t)) + diag(1e-8, 2)

  list(nu_t = nu_t, S_t = S_t)
}

update_q_rho_fast <- function(y, m_mat, V_mat, A, D_A, rho_values, prior_w, H) {
  M <- ncol(y)
  logw <- rep(-Inf, length(rho_values))

  for (tt in seq_along(rho_values)) {
    if (!is.finite(prior_w[tt]) || prior_w[tt] <= 0) next

    rho <- rho_values[tt]
    D_t <- D_A - rho * A
    D_t <- 0.5 * (D_t + t(D_t))

    co <- tryCatch(chol(D_t), error = function(e) NULL)
    if (is.null(co)) next

    ld <- 2 * sum(log(diag(co)))
    val <- M * ld

    S_suff_t <- S_sufficient_halves_fast(y, m_mat, V_mat, D_t)
    val <- val - 0.5 * sum(diag(H %*% S_suff_t))
    val <- val + log(prior_w[tt])

    logw[tt] <- val
  }

  finite_idx <- which(is.finite(logw))
  if (length(finite_idx) == 0) {
    w <- prior_w
  } else {
    mlp <- max(logw[finite_idx])
    w <- exp(logw - mlp)
    w[!is.finite(w)] <- 0
    s <- sum(w)
    if (!is.finite(s) || s <= 0) {
      w <- prior_w
    } else {
      w <- w / s
    }
  }

  list(w = w, logw = logw)
}

## --------------------------- Fast ELBO --------------------------- ##
elbo_core_fast <- function(y, A, D_A,
                           m_mat, V_mat,
                           logdet_Sig, quad_total,
                           v0, S0, nu_t, S_t,
                           rho_values, w_rho, logdetK,
                           mode = c("hetero", "homo"),
                           a_t, d_t, a0, d0,
                           prior_w = NULL,
                           debug = FALSE) {
  mode <- match.arg(mode)

  C <- nrow(y)
  M <- ncol(y)
  J <- C / 2
  p <- 2

  R_bar <- Reduce(`+`, lapply(seq_along(rho_values),
                              function(tt) w_rho[tt] * (D_A - rho_values[tt] * A)))

  H <- E_Sinv_IW(nu_t, S_t)
  E_log_Sigma <- E_logdet_Sigma_IW(nu_t, S_t, p)

  term_like_const <- M * (- C / 2 * log(2 * pi) - (J / 2) * E_log_Sigma)

  E_log_R <- 0
  for (tt in which(w_rho > 0)) {
    D_t <- D_A - rho_values[tt] * A
    if (is_spd(D_t)) {
      E_log_R <- E_log_R + w_rho[tt] * logdet_spd(D_t)
    }
  }
  term_like_det  <- (p / 2) * M * E_log_R

  S_suff_bar <- S_sufficient_halves_fast(y, m_mat, V_mat, R_bar)
  term_like_quad <- -0.5 * sum(diag(H %*% S_suff_bar))
  L_like <- term_like_const + term_like_det + term_like_quad

  if (mode == "hetero") {
    L_prior_b <- sum(
      - M / 2 * log(2 * pi)
      - 0.5 * (logdetK + M * E_log_sigma2_IG(a_t, d_t))
      - 0.5 * E_inv_sigma2_IG(a_t, d_t) * quad_total
    )

    L_prior_sig2 <- sum(
      a0 * log(d0) - lgamma(a0)
      - (a0 + 1) * E_log_sigma2_IG(a_t, d_t)
      - d0 * E_inv_sigma2_IG(a_t, d_t)
    )
  } else {
    Elog_s2 <- E_log_sigma2_IG(a_t, d_t)
    Einv_s2 <- E_inv_sigma2_IG(a_t, d_t)
    quad_sum <- sum(quad_total)

    L_prior_b <- C * (- M / 2 * log(2 * pi) - 0.5 * (logdetK + M * Elog_s2)) -
      0.5 * Einv_s2 * quad_sum

    L_prior_sig2 <- a0 * log(d0) - lgamma(a0) - (a0 + 1) * Elog_s2 - d0 * Einv_s2
  }

  const_IW <- (v0 / 2) * logdet_spd(S0) - (v0 * p / 2) * log(2) - log_multigamma(v0 / 2, p)
  L_prior_S <- const_IW - ((v0 + p + 1) / 2) * E_log_Sigma - 0.5 * sum(diag(S0 %*% H))

  if (is.null(prior_w)) {
    L_prior_rho <- -log(length(rho_values))
  } else {
    L_prior_rho <- sum(ifelse(w_rho > 0 & prior_w > 0, w_rho * log(prior_w), 0))
  }

  H_b <- sum(0.5 * logdet_Sig + 0.5 * M * (1 + log(2 * pi)))
  H_sig2 <- if (mode == "hetero") sum(entropy_inv_gamma(a_t, d_t)) else entropy_inv_gamma(a_t, d_t)

  const_IW_q <- (nu_t / 2) * logdet_spd(S_t) - (nu_t * p / 2) * log(2) - log_multigamma(nu_t / 2, p)
  E_log_q_S <- const_IW_q -
    ((nu_t + p + 1) / 2) * E_logdet_Sigma_IW(nu_t, S_t, p) -
    0.5 * sum(diag(S_t %*% E_Sinv_IW(nu_t, S_t)))
  H_S <- -E_log_q_S

  H_rho <- -sum(ifelse(w_rho > 0, w_rho * log(w_rho), 0))

  L_total <- L_like + L_prior_b + L_prior_sig2 + L_prior_S + L_prior_rho + H_b + H_sig2 + H_S + H_rho

  if (debug) {
    cat("chk E_log_Sigma =", E_log_Sigma, "\n")
    cat("chk term_like_det =", term_like_det, "\n")
    cat("chk term_like_quad =", term_like_quad, "\n")
    flush.console()
  }

  L_total
}

## --------------------------- Fast runners --------------------------- ##
run_vi_hetero <- function(SNP_path,
                          Y_path,
                          A_mat = NULL,
                          rho_grid = seq(0.00, 0.99, by = 0.01),
                          max_iter = 200,
                          tol = 1e-3,
                          verbose = TRUE) {
  seq_i <- read.table(SNP_path, header = FALSE, check.names = FALSE)
  X <- safe_scale_cols(seq_i)
  L <- ncol(seq_i)

  K <- (X %*% t(X)) / L
  diag_mean <- mean(diag(K))
  if (!is.finite(diag_mean) || diag_mean <= 0) stop("Invalid K: mean(diag(K)) must be positive.")
  K <- K / diag_mean
  Kpack <- make_K_eig_pack(K, eps = 1e-6)

  Y <- read_numeric_csv_matrix(Y_path, row_names = FALSE)
  y <- Y
  C <- nrow(y)
  M <- ncol(y)
  stopifnot(C %% 2 == 0)

  A <- if (is.null(A_mat)) build_A_from_Y_halves(y, k = 2) else normalize_A_spectral(A_mat)
  diag(A) <- 0
  D_A <- diag(rowSums(A))

  a0 <- rep(2.1, C)
  m0 <- rep(0.30, C)
  d0 <- m0 * (a0 - 1)

  v0 <- 8
  Sigma0_mean <- diag(c(0.40, 0.40))
  S0 <- Sigma0_mean * (v0 - 3)

  rho_prior <- make_uniform_rho_prior_valid(A, D_A, rho_grid)

  m_mat <- matrix(0, C, M)
  V_mat <- matrix(1, C, M)
  logdet_Sig <- rep(0, C)
  trace_invK_Sig <- rep(0, C)
  quad_total <- rep(0, C)

  a_t <- a0 + M / 2
  d_t <- d0 + 1
  nu_t <- v0 + M * (C / 2)
  S_t  <- S0 + diag(1, 2, 2)
  rho_w <- rho_prior
  elbo_hist <- c()
  beta_vec_final <- rep(NA_real_, C)

  for (iter in 1:max_iter) {
    R_bar <- Reduce(`+`, lapply(seq_along(rho_grid),
                                function(tt) rho_w[tt] * (D_A - rho_grid[tt] * A)))
    H <- E_Sinv_IW(nu_t, S_t)

    alpha_vec <- a_t / d_t
    qb <- update_q_b_fast(y, Kpack, R_bar, H, m_mat, alpha_vec)
    m_mat <- qb$m_mat
    V_mat <- qb$V_mat
    logdet_Sig <- qb$logdet_Sig
    trace_invK_Sig <- qb$trace_invK_Sig
    quad_total <- qb$quad_total
    beta_vec_final <- qb$beta_vec

    a_t <- a0 + M / 2
    d_t <- d0 + 0.5 * quad_total

    qS <- update_q_Sigma_halves_fast(y, m_mat, V_mat, R_bar, v0, S0)
    nu_t <- qS$nu_t
    S_t <- qS$S_t

    qr <- update_q_rho_fast(y, m_mat, V_mat, A, D_A, rho_grid, rho_prior, E_Sinv_IW(nu_t, S_t))
    rho_w <- qr$w

    Lcur <- elbo_core_fast(
      y, A, D_A,
      m_mat, V_mat, logdet_Sig, quad_total,
      v0, S0, nu_t, S_t, rho_grid, rho_w, Kpack$logdetK,
      mode = "hetero",
      a_t = a_t, d_t = d_t, a0 = a0, d0 = d0,
      prior_w = rho_prior,
      debug = (iter <= 2)
    )

    elbo_hist <- c(elbo_hist, Lcur)
    if (verbose) cat(sprintf("Iter %4d | ELBO = %.6f\n", iter, Lcur))
    if (iter >= 2) {
      delta <- elbo_hist[iter] - elbo_hist[iter - 1]
      if (is.finite(delta) && abs(delta) < tol) break
    }
  }

  E_Sigma  <- S_t / (nu_t - 3)
  E_rho    <- sum(rho_grid * rho_w)
  E_sigma2 <- as.numeric(d_t) / pmax(as.numeric(a_t) - 1, 1e-8)

  list(
    elbo = elbo_hist,
    mu = m_mat,
    Sig = NULL,
    Sig_diag = V_mat,
    Sig_logdet = logdet_Sig,
    Sig_trace_invK = trace_invK_Sig,
    a = a_t, d = d_t, a0 = a0, d0 = d0,
    nu = nu_t, S = S_t,
    E_Sigma = E_Sigma,
    E_sigma2 = E_sigma2,
    rho_grid = rho_grid, rho_w = rho_w, E_rho = E_rho,
    A = A, D_A = D_A, model = "hetero",
    Kpack = Kpack,
    final_alpha = a_t / d_t,
    final_beta = beta_vec_final
  )
}

run_vi_homo <- function(SNP_path,
                        Y_path,
                        A_mat = NULL,
                        rho_grid = seq(0.00, 0.99, by = 0.01),
                        max_iter = 200,
                        tol = 1e-3,
                        verbose = TRUE) {
  seq_i <- read.table(SNP_path, header = FALSE, check.names = FALSE)
  X <- safe_scale_cols(seq_i)
  L <- ncol(seq_i)

  K <- (X %*% t(X)) / L
  diag_mean <- mean(diag(K))
  if (!is.finite(diag_mean) || diag_mean <= 0) stop("Invalid K: mean(diag(K)) must be positive.")
  K <- K / diag_mean
  Kpack <- make_K_eig_pack(K, eps = 1e-6)

  Y <- read_numeric_csv_matrix(Y_path, row_names = FALSE)
  y <- Y
  C <- nrow(y)
  M <- ncol(y)
  stopifnot(C %% 2 == 0)

  A <- if (is.null(A_mat)) build_A_from_Y_halves(y, k = 2) else normalize_A_spectral(A_mat)
  diag(A) <- 0
  D_A <- diag(rowSums(A))

  a0 <- 2.1
  m0 <- 0.30
  d0 <- m0 * (a0 - 1)

  v0 <- 8
  Sigma0_mean <- diag(c(0.40, 0.40))
  S0 <- Sigma0_mean * (v0 - 3)

  rho_prior <- make_uniform_rho_prior_valid(A, D_A, rho_grid)

  m_mat <- matrix(0, C, M)
  V_mat <- matrix(1, C, M)
  logdet_Sig <- rep(0, C)
  trace_invK_Sig <- rep(0, C)
  quad_total <- rep(0, C)

  a_t <- a0 + (C * M) / 2
  d_t <- d0 + 1
  nu_t <- v0 + M * (C / 2)
  S_t  <- S0 + diag(1, 2, 2)
  rho_w <- rho_prior
  elbo_hist <- c()
  beta_vec_final <- rep(NA_real_, C)

  for (iter in 1:max_iter) {
    R_bar <- Reduce(`+`, lapply(seq_along(rho_grid),
                                function(tt) rho_w[tt] * (D_A - rho_grid[tt] * A)))
    H <- E_Sinv_IW(nu_t, S_t)

    alpha_scalar <- as.numeric(a_t / d_t)
    qb <- update_q_b_fast(y, Kpack, R_bar, H, m_mat, rep(alpha_scalar, C))
    m_mat <- qb$m_mat
    V_mat <- qb$V_mat
    logdet_Sig <- qb$logdet_Sig
    trace_invK_Sig <- qb$trace_invK_Sig
    quad_total <- qb$quad_total
    beta_vec_final <- qb$beta_vec

    a_t <- a0 + (C * M) / 2
    d_t <- d0 + 0.5 * sum(quad_total)

    qS <- update_q_Sigma_halves_fast(y, m_mat, V_mat, R_bar, v0, S0)
    nu_t <- qS$nu_t
    S_t <- qS$S_t

    qr <- update_q_rho_fast(y, m_mat, V_mat, A, D_A, rho_grid, rho_prior, E_Sinv_IW(nu_t, S_t))
    rho_w <- qr$w

    Lcur <- elbo_core_fast(
      y, A, D_A,
      m_mat, V_mat, logdet_Sig, quad_total,
      v0, S0, nu_t, S_t, rho_grid, rho_w, Kpack$logdetK,
      mode = "homo",
      a_t = a_t, d_t = d_t, a0 = a0, d0 = d0,
      prior_w = rho_prior,
      debug = (iter <= 2)
    )

    elbo_hist <- c(elbo_hist, Lcur)
    if (verbose) cat(sprintf("Iter %4d | ELBO = %.6f\n", iter, Lcur))
    if (iter >= 2) {
      delta <- elbo_hist[iter] - elbo_hist[iter - 1]
      if (is.finite(delta) && abs(delta) < tol) break
    }
  }

  E_Sigma  <- S_t / (nu_t - 3)
  E_rho    <- sum(rho_grid * rho_w)
  E_sigma2 <- rep(as.numeric(d_t) / pmax(as.numeric(a_t) - 1, 1e-8), C)

  list(
    elbo = elbo_hist,
    mu = m_mat,
    Sig = NULL,
    Sig_diag = V_mat,
    Sig_logdet = logdet_Sig,
    Sig_trace_invK = trace_invK_Sig,
    a = a_t, d = d_t, a0 = a0, d0 = d0,
    nu = nu_t, S = S_t,
    E_Sigma = E_Sigma,
    E_sigma2 = E_sigma2,
    rho_grid = rho_grid, rho_w = rho_w, E_rho = E_rho,
    A = A, D_A = D_A, model = "homo",
    Kpack = Kpack,
    final_alpha = rep(a_t / d_t, C),
    final_beta = beta_vec_final
  )
}

run_vi_hetero_from_mats <- function(Y, K, A_mat = NULL,
                                    rho_grid = seq(0.00, 0.99, by = 0.01),
                                    max_iter = 5000,
                                    tol = 1e-4,
                                    verbose = TRUE) {
  K <- as.matrix(K)
  diag_mean <- mean(diag(K))
  if (!is.finite(diag_mean) || diag_mean <= 0) stop("Invalid K: mean(diag(K)) must be positive.")
  K <- K / diag_mean
  Kpack <- make_K_eig_pack(K, eps = 1e-6)

  y <- as.matrix(Y)
  C <- nrow(y)
  M <- ncol(y)
  stopifnot(C %% 2 == 0)

  A <- if (is.null(A_mat)) build_A_from_Y_halves(y, k = 2) else normalize_A_spectral(A_mat)
  diag(A) <- 0
  D_A <- diag(rowSums(A))

  a0 <- rep(2.1, C)
  m0 <- rep(0.30, C)
  d0 <- m0 * (a0 - 1)

  v0 <- 8
  Sigma0_mean <- diag(c(0.40, 0.40))
  S0 <- Sigma0_mean * (v0 - 3)

  rho_prior <- make_uniform_rho_prior_valid(A, D_A, rho_grid)

  m_mat <- matrix(0, C, M)
  V_mat <- matrix(1, C, M)
  logdet_Sig <- rep(0, C)
  trace_invK_Sig <- rep(0, C)
  quad_total <- rep(0, C)

  a_t <- a0 + M / 2
  d_t <- d0 + 1
  nu_t <- v0 + M * (C / 2)
  S_t  <- S0 + diag(1, 2, 2)
  rho_w <- rho_prior
  elbo_hist <- c()
  beta_vec_final <- rep(NA_real_, C)

  for (iter in 1:max_iter) {
    R_bar <- Reduce(`+`, lapply(seq_along(rho_grid),
                                function(tt) rho_w[tt] * (D_A - rho_grid[tt] * A)))
    H <- E_Sinv_IW(nu_t, S_t)

    alpha_vec <- a_t / d_t
    qb <- update_q_b_fast(y, Kpack, R_bar, H, m_mat, alpha_vec)
    m_mat <- qb$m_mat
    V_mat <- qb$V_mat
    logdet_Sig <- qb$logdet_Sig
    trace_invK_Sig <- qb$trace_invK_Sig
    quad_total <- qb$quad_total
    beta_vec_final <- qb$beta_vec

    a_t <- a0 + M / 2
    d_t <- d0 + 0.5 * quad_total

    qS <- update_q_Sigma_halves_fast(y, m_mat, V_mat, R_bar, v0, S0)
    nu_t <- qS$nu_t
    S_t <- qS$S_t

    qr <- update_q_rho_fast(y, m_mat, V_mat, A, D_A, rho_grid, rho_prior, E_Sinv_IW(nu_t, S_t))
    rho_w <- qr$w

    Lcur <- elbo_core_fast(
      y, A, D_A,
      m_mat, V_mat, logdet_Sig, quad_total,
      v0, S0, nu_t, S_t, rho_grid, rho_w, Kpack$logdetK,
      mode = "hetero",
      a_t = a_t, d_t = d_t, a0 = a0, d0 = d0,
      prior_w = rho_prior,
      debug = (iter <= 2)
    )

    elbo_hist <- c(elbo_hist, Lcur)
    if (verbose) cat(sprintf("Iter %4d | ELBO = %.6f\n", iter, Lcur))
    if (iter >= 2) {
      delta <- elbo_hist[iter] - elbo_hist[iter - 1]
      if (is.finite(delta) && abs(delta) < tol) break
    }
  }

  E_Sigma  <- S_t / (nu_t - 3)
  E_rho    <- sum(rho_grid * rho_w)
  E_sigma2 <- as.numeric(d_t) / pmax(as.numeric(a_t) - 1, 1e-8)

  list(
    elbo = elbo_hist,
    mu = m_mat,
    Sig = NULL,
    Sig_diag = V_mat,
    Sig_logdet = logdet_Sig,
    Sig_trace_invK = trace_invK_Sig,
    a = a_t, d = d_t, a0 = a0, d0 = d0,
    nu = nu_t, S = S_t,
    E_Sigma = E_Sigma,
    E_sigma2 = E_sigma2,
    rho_grid = rho_grid, rho_w = rho_w, E_rho = E_rho,
    A = A, D_A = D_A, model = "hetero",
    Kpack = Kpack,
    final_alpha = a_t / d_t,
    final_beta = beta_vec_final
  )
}

## --------------------------- One-call wrappers --------------------------- ##
run_sim_from_mats_hetero <- function(Y, K, A_true,
                                     rho_grid = seq(0.00, 0.99, by = 0.01),
                                     max_iter = 5000,
                                     tol = 1e-4,
                                     verbose = TRUE,
                                     prefix = NULL) {
  res <- run_vi_hetero_from_mats(
    Y = Y,
    K = K,
    A_mat = A_true,
    rho_grid = rho_grid,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose
  )
  h <- heritability_from_vi(res, A_mat = A_true)

  out <- list(res = res, h = h)
  if (!is.null(prefix)) save_vi_result(res, prefix = prefix)
  out
}

run_realdata_brain_hetero_halves <- function(brain_csv,
                                             use_k_from = c("X", "K"),
                                             X = NULL,
                                             K = NULL,
                                             A_mat = NULL,
                                             row_names = TRUE,
                                             zscore_rows = TRUE,
                                             k_knn = 2,
                                             rho_grid = seq(0.00, 0.99, by = 0.01),
                                             max_iter = 10000,
                                             tol = 1e-4,
                                             verbose = TRUE,
                                             prefix = "brain_real_hetero_halves",
                                             make_plot = TRUE) {
  use_k_from <- match.arg(use_k_from)

  BrainMeasures <- read_numeric_csv_matrix(brain_csv, row_names = row_names)
  Y <- if (zscore_rows) safe_scale_rows(BrainMeasures) else BrainMeasures

  C <- nrow(Y)
  M <- ncol(Y)
  stopifnot(C %% 2 == 0)

  if (use_k_from == "X") {
    if (is.null(X)) stop("When use_k_from = 'X', you must provide X.")
    X <- as.matrix(X)
    if (nrow(X) != M) stop("nrow(X) must equal ncol(Y): sample sizes do not match.")
    std_X <- safe_scale_cols(X)
    L <- ncol(std_X)
    K_use <- (std_X %*% t(std_X)) / L
  } else {
    if (is.null(K)) stop("When use_k_from = 'K', you must provide K.")
    K_use <- as.matrix(K)
    if (nrow(K_use) != M || ncol(K_use) != M) stop("K must be M x M, where M = ncol(Y).")
  }

  A_use <- if (is.null(A_mat)) build_A_from_Y_halves(Y, k = k_knn) else normalize_A_spectral(A_mat)

  res <- run_vi_hetero_from_mats(
    Y = Y,
    K = K_use,
    A_mat = A_use,
    rho_grid = rho_grid,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose
  )

  h <- heritability_from_vi(res, A_mat = A_use)

  if (make_plot) {
    plot(seq_along(res$elbo), res$elbo, type = "l",
         main = "ELBO (real data, hetero, halves)",
         xlab = "Iteration", ylab = "ELBO")
  }

  cat(sprintf(
    "\n[Summary]\n  E[rho] = %.4f\n  Broad H^2 = %.4f\n  mean(narrow h^2) = %.4f (sd=%.4f)\n",
    res$E_rho, h$H2, mean(h$h2), sd(h$h2)
  ))
  cat("  First 10 narrow h^2:\n")
  print(round(head(h$h2, 10), 4))

  if (!is.null(prefix)) save_vi_result(res, prefix = prefix)

  list(Y = Y, K = K_use, A = A_use, res = res, h = h)
}

## --------------------------- Heritability --------------------------- ##
heritability_from_vi <- function(res, A_mat = NULL, pairing = c("halves")) {
  pairing <- match.arg(pairing)

  C <- length(res$E_sigma2)
  J <- C / 2
  stopifnot(C %% 2 == 0)

  E_sigma2 <- as.numeric(res$E_sigma2)
  E_Sigma  <- res$E_Sigma

  A <- if (is.null(A_mat)) res$A else normalize_A_spectral(A_mat)
  diag(A) <- 0
  D_A <- diag(rowSums(A))

  R_bar <- Reduce(`+`, lapply(seq_along(res$rho_grid),
                              function(tt) res$rho_w[tt] * (D_A - res$rho_grid[tt] * A)))
  R_inv <- inv_spd(R_bar)

  var_e <- numeric(C)
  rjj <- diag(R_inv)
  var_e[1:J]       <- rjj * E_Sigma[1, 1]
  var_e[(J + 1):C] <- rjj * E_Sigma[2, 2]

  h2_narrow <- E_sigma2 / pmax(E_sigma2 + var_e, 1e-12)

  Sigma_g <- diag(E_sigma2, nrow = C, ncol = C)
  Sigma_e <- kronecker(E_Sigma, R_inv)

  L <- chol(0.5 * (Sigma_e + t(Sigma_e)) + diag(1e-10, C))
  Se_inv_half <- backsolve(L, diag(C))
  Mmat <- Se_inv_half %*% Sigma_g %*% t(Se_inv_half)

  lam_max <- max(eigen(0.5 * (Mmat + t(Mmat)), symmetric = TRUE, only.values = TRUE)$values)
  lam_max <- max(lam_max, 0)

  H2_broad <- lam_max / (1 + lam_max)

  list(
    h2 = h2_narrow,
    H2 = H2_broad,
    pieces = list(
      E_sigma2 = E_sigma2,
      var_e = var_e,
      R_bar = R_bar,
      R_inv = R_inv,
      Sigma_g = Sigma_g,
      Sigma_e = Sigma_e
    )
  )
}
