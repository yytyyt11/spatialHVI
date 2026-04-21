# Internal mathematical helpers for the variational inference code.

inv_spd <- function(mat, jitter = 1e-8, max_tries = 6L) {
  mat <- validate_square_matrix(mat, arg = "mat")
  for (attempt in 0:max_tries) {
    candidate <- 0.5 * (mat + t(mat)) + diag(jitter * (10^attempt), nrow(mat))
    chol_factor <- tryCatch(chol(candidate), error = function(e) NULL)
    if (!is.null(chol_factor)) {
      return(chol2inv(chol_factor))
    }
  }
  stop("Failed to invert a symmetric positive-definite matrix.", call. = FALSE)
}

logdet_spd <- function(mat, jitter = 1e-8, max_tries = 6L) {
  mat <- validate_square_matrix(mat, arg = "mat")
  for (attempt in 0:max_tries) {
    candidate <- 0.5 * (mat + t(mat)) + diag(jitter * (10^attempt), nrow(mat))
    chol_factor <- tryCatch(chol(candidate), error = function(e) NULL)
    if (!is.null(chol_factor)) {
      return(2 * sum(log(diag(chol_factor))))
    }
  }
  stop("Failed to compute a log-determinant for a positive-definite matrix.", call. = FALSE)
}

is_spd <- function(mat) {
  mat <- validate_square_matrix(mat, arg = "mat")
  candidate <- 0.5 * (mat + t(mat))
  !is.null(tryCatch(chol(candidate), error = function(e) NULL))
}

make_k_eig_pack <- function(k_mat, eps = 1e-6) {
  k_mat <- validate_square_matrix(k_mat, arg = "k_mat")
  k_jittered <- 0.5 * (k_mat + t(k_mat)) + diag(eps, nrow(k_mat))
  eig <- eigen(k_jittered, symmetric = TRUE)
  lambda <- as.numeric(eig$values)
  u_mat <- eig$vectors

  list(
    K = k_jittered,
    U = u_mat,
    U2 = u_mat * u_mat,
    lambda = lambda,
    inv_lambda = 1 / lambda,
    logdetK = sum(log(lambda))
  )
}

log_multigamma <- function(a, p) {
  (p * (p - 1) / 4) * log(pi) + sum(lgamma(a + (1 - (1:p)) / 2))
}

e_logdet_sigma_iw <- function(nu, s_mat, p = 2L) {
  logdet_spd(s_mat) - sum(digamma((nu + 1 - (1:p)) / 2)) - p * log(2)
}

e_sinv_iw <- function(nu, s_mat) {
  nu * inv_spd(s_mat)
}

entropy_inv_gamma <- function(a, b) {
  a + log(b) + lgamma(a) - (1 + a) * digamma(a)
}

e_log_sigma2_ig <- function(a, b) {
  log(b) - digamma(a)
}

e_inv_sigma2_ig <- function(a, b) {
  a / b
}

vec_to_2xj_halves <- function(x) {
  c_traits <- length(x)
  j_pairs <- c_traits / 2L
  rbind(x[seq_len(j_pairs)], x[(j_pairs + 1L):c_traits])
}

s_sufficient_halves_fast <- function(y, m_mat, v_mat, r_mat) {
  c_traits <- nrow(y)
  subject_count <- ncol(y)
  j_pairs <- c_traits / 2L
  acc <- matrix(0, nrow = 2L, ncol = 2L)

  for (subject_idx in seq_len(subject_count)) {
    y_pair <- vec_to_2xj_halves(y[, subject_idx])
    m_pair <- vec_to_2xj_halves(m_mat[, subject_idx])
    residual_pair <- y_pair - m_pair

    acc <- acc + residual_pair %*% r_mat %*% t(residual_pair)
    for (pair_idx in seq_len(j_pairs)) {
      acc <- acc + r_mat[pair_idx, pair_idx] * diag(
        c(v_mat[pair_idx, subject_idx], v_mat[pair_idx + j_pairs, subject_idx]),
        nrow = 2L
      )
    }
  }

  acc
}

update_q_b_fast <- function(y, k_pack, r_mat, h_mat, m_mat, alpha_vec) {
  c_traits <- nrow(y)
  subject_count <- ncol(y)

  g_mat <- kronecker(h_mat, r_mat)
  gy_mat <- g_mat %*% y

  u_mat <- k_pack$U
  u_sq <- k_pack$U2
  inv_lambda <- k_pack$inv_lambda

  v_mat <- matrix(0, nrow = c_traits, ncol = subject_count)
  logdet_sig <- numeric(c_traits)
  trace_inv_k_sig <- numeric(c_traits)
  quad_total <- numeric(c_traits)
  beta_vec <- diag(g_mat)

  for (trait_idx in seq_len(c_traits)) {
    g_row <- g_mat[trait_idx, ]
    beta_val <- beta_vec[trait_idx]

    rhs <- as.numeric(gy_mat[trait_idx, ] - (g_row %*% m_mat - beta_val * m_mat[trait_idx, ]))
    rhs_u <- as.numeric(crossprod(u_mat, rhs))

    denom <- alpha_vec[trait_idx] * inv_lambda + beta_val
    diag_terms <- 1 / denom

    mean_u <- diag_terms * rhs_u
    mean_trait <- as.numeric(u_mat %*% mean_u)

    m_mat[trait_idx, ] <- mean_trait
    v_mat[trait_idx, ] <- as.numeric(u_sq %*% diag_terms)

    logdet_sig[trait_idx] <- sum(log(diag_terms))
    trace_inv_k_sig[trait_idx] <- sum(inv_lambda * diag_terms)
    quad_total[trait_idx] <- sum(inv_lambda * (mean_u^2)) + trace_inv_k_sig[trait_idx]
  }

  list(
    m_mat = m_mat,
    v_mat = v_mat,
    logdet_sig = logdet_sig,
    trace_inv_k_sig = trace_inv_k_sig,
    quad_total = quad_total,
    beta_vec = beta_vec
  )
}

update_q_sigma_halves_fast <- function(y, m_mat, v_mat, r_mat, v0, s0) {
  subject_count <- ncol(y)
  j_pairs <- nrow(y) / 2L

  nu_t <- v0 + subject_count * j_pairs
  s_t <- s0 + s_sufficient_halves_fast(y, m_mat, v_mat, r_mat)
  s_t <- 0.5 * (s_t + t(s_t)) + diag(1e-8, nrow = 2L)

  list(nu_t = nu_t, s_t = s_t)
}

update_q_rho_fast <- function(y, m_mat, v_mat, a_mat, d_a, rho_grid, prior_w, h_mat) {
  subject_count <- ncol(y)
  log_w <- rep(-Inf, length(rho_grid))

  for (idx in seq_along(rho_grid)) {
    if (!is.finite(prior_w[idx]) || prior_w[idx] <= 0) {
      next
    }

    rho <- rho_grid[idx]
    precision <- d_a - rho * a_mat
    precision <- 0.5 * (precision + t(precision))
    chol_factor <- tryCatch(chol(precision), error = function(e) NULL)
    if (is.null(chol_factor)) {
      next
    }

    logdet_precision <- 2 * sum(log(diag(chol_factor)))
    s_suff <- s_sufficient_halves_fast(y, m_mat, v_mat, precision)
    log_w[idx] <- subject_count * logdet_precision -
      0.5 * sum(diag(h_mat %*% s_suff)) +
      log(prior_w[idx])
  }

  finite_idx <- which(is.finite(log_w))
  if (length(finite_idx) == 0L) {
    return(list(w = prior_w, log_w = log_w))
  }

  max_log_w <- max(log_w[finite_idx])
  w <- exp(log_w - max_log_w)
  w[!is.finite(w)] <- 0
  if (sum(w) <= 0) {
    w <- prior_w
  } else {
    w <- w / sum(w)
  }

  list(w = w, log_w = log_w)
}

weighted_precision_matrix <- function(rho_grid, rho_w, d_a, a_mat) {
  precision <- matrix(0, nrow = nrow(a_mat), ncol = ncol(a_mat))
  for (idx in seq_along(rho_grid)) {
    precision <- precision + rho_w[idx] * (d_a - rho_grid[idx] * a_mat)
  }
  precision
}

elbo_core_fast <- function(
  y,
  a_mat,
  d_a,
  m_mat,
  v_mat,
  logdet_sig,
  quad_total,
  v0,
  s0,
  nu_t,
  s_t,
  rho_grid,
  rho_w,
  logdet_k,
  mode = c("hetero", "homo"),
  a_t,
  d_t,
  a0,
  d0,
  prior_w = NULL
) {
  mode <- match.arg(mode)

  c_traits <- nrow(y)
  subject_count <- ncol(y)
  j_pairs <- c_traits / 2L
  p_dim <- 2L

  r_bar <- weighted_precision_matrix(rho_grid, rho_w, d_a, a_mat)
  h_mat <- e_sinv_iw(nu_t, s_t)
  e_log_sigma <- e_logdet_sigma_iw(nu_t, s_t, p = p_dim)

  term_like_const <- subject_count * (-c_traits / 2 * log(2 * pi) - (j_pairs / 2) * e_log_sigma)

  e_log_r <- 0
  for (idx in which(rho_w > 0)) {
    precision <- d_a - rho_grid[idx] * a_mat
    if (is_spd(precision)) {
      e_log_r <- e_log_r + rho_w[idx] * logdet_spd(precision)
    }
  }

  term_like_det <- (p_dim / 2) * subject_count * e_log_r
  s_suff_bar <- s_sufficient_halves_fast(y, m_mat, v_mat, r_bar)
  term_like_quad <- -0.5 * sum(diag(h_mat %*% s_suff_bar))
  l_like <- term_like_const + term_like_det + term_like_quad

  if (mode == "hetero") {
    l_prior_b <- sum(
      -subject_count / 2 * log(2 * pi) -
        0.5 * (logdet_k + subject_count * e_log_sigma2_ig(a_t, d_t)) -
        0.5 * e_inv_sigma2_ig(a_t, d_t) * quad_total
    )

    l_prior_sigma2 <- sum(
      a0 * log(d0) - lgamma(a0) -
        (a0 + 1) * e_log_sigma2_ig(a_t, d_t) -
        d0 * e_inv_sigma2_ig(a_t, d_t)
    )
  } else {
    e_log_s2 <- e_log_sigma2_ig(a_t, d_t)
    e_inv_s2 <- e_inv_sigma2_ig(a_t, d_t)
    quad_sum <- sum(quad_total)

    l_prior_b <- c_traits * (
      -subject_count / 2 * log(2 * pi) - 0.5 * (logdet_k + subject_count * e_log_s2)
    ) - 0.5 * e_inv_s2 * quad_sum

    l_prior_sigma2 <- a0 * log(d0) - lgamma(a0) - (a0 + 1) * e_log_s2 - d0 * e_inv_s2
  }

  const_iw <- (v0 / 2) * logdet_spd(s0) -
    (v0 * p_dim / 2) * log(2) -
    log_multigamma(v0 / 2, p_dim)
  l_prior_sigma <- const_iw -
    ((v0 + p_dim + 1) / 2) * e_log_sigma -
    0.5 * sum(diag(s0 %*% h_mat))

  if (is.null(prior_w)) {
    l_prior_rho <- -log(length(rho_grid))
  } else {
    l_prior_rho <- sum(ifelse(rho_w > 0 & prior_w > 0, rho_w * log(prior_w), 0))
  }

  h_b <- sum(0.5 * logdet_sig + 0.5 * subject_count * (1 + log(2 * pi)))
  h_sigma2 <- if (mode == "hetero") {
    sum(entropy_inv_gamma(a_t, d_t))
  } else {
    entropy_inv_gamma(a_t, d_t)
  }

  const_iw_q <- (nu_t / 2) * logdet_spd(s_t) -
    (nu_t * p_dim / 2) * log(2) -
    log_multigamma(nu_t / 2, p_dim)
  e_log_q_sigma <- const_iw_q -
    ((nu_t + p_dim + 1) / 2) * e_logdet_sigma_iw(nu_t, s_t, p = p_dim) -
    0.5 * sum(diag(s_t %*% e_sinv_iw(nu_t, s_t)))
  h_sigma <- -e_log_q_sigma

  h_rho <- -sum(ifelse(rho_w > 0, rho_w * log(rho_w), 0))

  l_like + l_prior_b + l_prior_sigma2 + l_prior_sigma + l_prior_rho +
    h_b + h_sigma2 + h_sigma + h_rho
}

run_vi_core <- function(
  y,
  k_mat,
  a_mat = NULL,
  rho_grid = seq(0, 0.99, by = 0.01),
  max_iter = 5000,
  tol = 1e-4,
  verbose = TRUE,
  model = c("hetero", "homo")
) {
  model <- match.arg(model)
  y <- validate_even_trait_matrix(y, arg = "y")
  k_mat <- normalize_kinship_matrix(k_mat)

  if (nrow(k_mat) != ncol(y) || ncol(k_mat) != ncol(y)) {
    stop("`k_mat` must be square with dimension equal to the number of subjects in `y`.", call. = FALSE)
  }

  rho_grid <- sort(unique(as.numeric(rho_grid)))
  if (length(rho_grid) == 0L || any(!is.finite(rho_grid))) {
    stop("`rho_grid` must contain at least one finite numeric value.", call. = FALSE)
  }

  max_iter <- validate_positive_count(max_iter, arg = "max_iter")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("`tol` must be a non-negative number.", call. = FALSE)
  }
  verbose <- isTRUE(verbose)

  k_pack <- make_k_eig_pack(k_mat, eps = 1e-6)

  c_traits <- nrow(y)
  j_pairs <- c_traits / 2L
  a_use <- if (is.null(a_mat)) {
    build_a_from_y_halves(y, k = 2)
  } else {
    normalize_a_spectral(a_mat)
  }

  if (nrow(a_use) != j_pairs || ncol(a_use) != j_pairs) {
    stop("`a_mat` must be a square matrix with dimension equal to `nrow(y) / 2`.", call. = FALSE)
  }

  diag(a_use) <- 0
  d_a <- make_degree_matrix(a_use)
  rho_prior <- make_uniform_rho_prior_valid(a_use, d_a, rho_grid)

  subject_count <- ncol(y)
  m_mat <- matrix(0, nrow = c_traits, ncol = subject_count)
  v_mat <- matrix(1, nrow = c_traits, ncol = subject_count)
  logdet_sig <- rep(0, c_traits)
  trace_inv_k_sig <- rep(0, c_traits)
  quad_total <- rep(0, c_traits)
  beta_vec_final <- rep(NA_real_, c_traits)

  if (model == "hetero") {
    a0 <- rep(2.1, c_traits)
    d0 <- rep(0.30 * (2.1 - 1), c_traits)
    a_t <- a0 + subject_count / 2
    d_t <- d0 + 1
  } else {
    a0 <- 2.1
    d0 <- 0.30 * (a0 - 1)
    a_t <- a0 + (c_traits * subject_count) / 2
    d_t <- d0 + 1
  }

  v0 <- 8
  s0 <- diag(c(0.40, 0.40)) * (v0 - 3)
  nu_t <- v0 + subject_count * j_pairs
  s_t <- s0 + diag(1, nrow = 2L)
  rho_w <- rho_prior
  elbo_hist <- numeric(0)

  for (iter in seq_len(max_iter)) {
    r_bar <- weighted_precision_matrix(rho_grid, rho_w, d_a, a_use)
    h_mat <- e_sinv_iw(nu_t, s_t)

    alpha_vec <- if (model == "hetero") {
      a_t / d_t
    } else {
      rep(as.numeric(a_t / d_t), c_traits)
    }

    q_b <- update_q_b_fast(y, k_pack, r_bar, h_mat, m_mat, alpha_vec)
    m_mat <- q_b$m_mat
    v_mat <- q_b$v_mat
    logdet_sig <- q_b$logdet_sig
    trace_inv_k_sig <- q_b$trace_inv_k_sig
    quad_total <- q_b$quad_total
    beta_vec_final <- q_b$beta_vec

    if (model == "hetero") {
      a_t <- a0 + subject_count / 2
      d_t <- d0 + 0.5 * quad_total
    } else {
      a_t <- a0 + (c_traits * subject_count) / 2
      d_t <- d0 + 0.5 * sum(quad_total)
    }

    q_sigma <- update_q_sigma_halves_fast(y, m_mat, v_mat, r_bar, v0, s0)
    nu_t <- q_sigma$nu_t
    s_t <- q_sigma$s_t

    q_rho <- update_q_rho_fast(y, m_mat, v_mat, a_use, d_a, rho_grid, rho_prior, e_sinv_iw(nu_t, s_t))
    rho_w <- q_rho$w

    elbo_now <- elbo_core_fast(
      y = y,
      a_mat = a_use,
      d_a = d_a,
      m_mat = m_mat,
      v_mat = v_mat,
      logdet_sig = logdet_sig,
      quad_total = quad_total,
      v0 = v0,
      s0 = s0,
      nu_t = nu_t,
      s_t = s_t,
      rho_grid = rho_grid,
      rho_w = rho_w,
      logdet_k = k_pack$logdetK,
      mode = model,
      a_t = a_t,
      d_t = d_t,
      a0 = a0,
      d0 = d0,
      prior_w = rho_prior
    )

    elbo_hist <- c(elbo_hist, elbo_now)
    if (verbose) {
      message(sprintf("Iter %4d | ELBO = %.6f", iter, elbo_now))
    }

    if (iter >= 2L) {
      delta <- elbo_hist[iter] - elbo_hist[iter - 1L]
      if (is.finite(delta) && abs(delta) < tol) {
        break
      }
    }
  }

  e_sigma2 <- if (model == "hetero") {
    as.numeric(d_t) / pmax(as.numeric(a_t) - 1, 1e-8)
  } else {
    rep(as.numeric(d_t) / pmax(as.numeric(a_t) - 1, 1e-8), c_traits)
  }

  list(
    elbo = elbo_hist,
    mu = m_mat,
    Sig = NULL,
    Sig_diag = v_mat,
    Sig_logdet = logdet_sig,
    Sig_trace_invK = trace_inv_k_sig,
    a = a_t,
    d = d_t,
    a0 = a0,
    d0 = d0,
    nu = nu_t,
    S = s_t,
    E_Sigma = s_t / (nu_t - 3),
    E_sigma2 = e_sigma2,
    rho_grid = rho_grid,
    rho_w = rho_w,
    E_rho = sum(rho_grid * rho_w),
    A = a_use,
    D_A = d_a,
    model = model,
    Kpack = k_pack,
    final_alpha = alpha_vec,
    final_beta = beta_vec_final
  )
}
