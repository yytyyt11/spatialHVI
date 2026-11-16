# Internal numerical functions extracted from the manuscript reference implementation.
# See inst/provenance/manuscript_code.tsv. Only return-path metadata is adapted.
.manuscript_engine <- local({
  cor <- stats::cor
  plogis <- stats::plogis
  qlogis <- stats::qlogis
  qnorm <- stats::qnorm
  median <- stats::median
  var <- stats::var
  clock <- function() proc.time()[["elapsed"]]
  cfg <- list(revalidation = list(maximum_updates = 40L, minimum_updates = 3L, stagnation_window = 6L, per_solve_timeout_seconds = 15))
  .reviewer_v3_sha256 <- "a6fc865730665dd106d8dcc6acbd219ca617821f2e52676534d988b1e9bc7dbb"
  .reviewer_library_sha256 <- "a66437b72f523e5a081accf33f1726357bccc55df4e5792839239b9e12984f91"

  .assert_matrix <- function(M, name, expected_nrow = NULL, expected_ncol = NULL, symmetric = FALSE, tol = 1e-10) {
      if (!is.matrix(M) || !is.numeric(M)) {
          stop(name, " must be a numeric matrix.", call. = FALSE)
      }
      if (!is.null(expected_nrow) && nrow(M) != expected_nrow) {
          stop(name, " has ", nrow(M), " rows; expected ", expected_nrow, ".", call. = FALSE)
      }
      if (!is.null(expected_ncol) && ncol(M) != expected_ncol) {
          stop(name, " has ", ncol(M), " columns; expected ", expected_ncol, ".", call. = FALSE)
      }
      if (any(!is.finite(M)))
          stop(name, " contains non-finite values.", call. = FALSE)
      if (symmetric && nrow(M) != ncol(M)) {
          stop(name, " must be square.", call. = FALSE)
      }
      if (symmetric && max(abs(M - t(M))) > tol) {
          stop(name, " is not symmetric within tolerance ", tol, ".", call. = FALSE)
      }
      invisible(TRUE)
  }

  .assert_psd <- function(M, name = deparse(substitute(M)), tol = 1e-09) {
      .assert_matrix(M, name, symmetric = TRUE, tol = tol)
      ev <- eigen(.symmetrize(M), symmetric = TRUE, only.values = TRUE)$values
      scale <- max(1, max(abs(ev)))
      if (min(ev) < -tol * scale) {
          stop(name, " is not positive semidefinite; min eigenvalue = ", format(min(ev), digits = 8), ".",
              call. = FALSE)
      }
      invisible(ev)
  }

  .assert_spd <- function(M, name = deparse(substitute(M)), tol = 1e-10) {
      .assert_matrix(M, name, symmetric = TRUE, tol = tol)
      co <- tryCatch(chol(.symmetrize(M)), error = function(e) NULL)
      if (is.null(co))
          stop(name, " is not strictly positive definite.", call. = FALSE)
      invisible(co)
  }

  .elbo_core_custom <- function(y, A, D_A, m_mat, V_mat, logdet_Sig, quad_total, prior, nu_t, S_t, rho_values,
      w_rho, rho_pack, logdetK, a_t, d_t) {
      C <- nrow(y)
      M <- ncol(y)
      J <- C/2L
      p <- 2L
      R_bar <- D_A - sum(rho_values * w_rho) * A
      H <- E_Sinv_IW(nu_t, S_t)
      E_log_Sigma <- E_logdet_Sigma_IW(nu_t, S_t, p)
      term_like_const <- M * (-C/2 * log(2 * pi) - (J/2) * E_log_Sigma)
      E_log_R <- sum(w_rho[rho_pack$valid] * rho_pack$logdet[rho_pack$valid])
      term_like_det <- (p/2) * M * E_log_R
      S_suff_bar <- S_sufficient_halves_vectorized(y, m_mat, V_mat, R_bar)
      term_like_quad <- -0.5 * sum(diag(H %*% S_suff_bar))
      L_like <- term_like_const + term_like_det + term_like_quad
      L_prior_b <- sum(-M/2 * log(2 * pi) - 0.5 * (logdetK + M * E_log_sigma2_IG(a_t, d_t)) - 0.5 * E_inv_sigma2_IG(a_t,
          d_t) * quad_total)
      L_prior_sig2 <- sum(prior$a0 * log(prior$d0) - lgamma(prior$a0) - (prior$a0 + 1) * E_log_sigma2_IG(a_t,
          d_t) - prior$d0 * E_inv_sigma2_IG(a_t, d_t))
      const_IW <- (prior$v0/2) * logdet_spd(prior$S0) - (prior$v0 * p/2) * log(2) - log_multigamma(prior$v0/2,
          p)
      L_prior_S <- const_IW - ((prior$v0 + p + 1)/2) * E_log_Sigma - 0.5 * sum(diag(prior$S0 %*% H))
      rho_idx <- which(w_rho > 0)
      L_prior_rho <- sum(w_rho[rho_idx] * log(rho_pack$rho_prior[rho_idx]))
      H_b <- sum(0.5 * logdet_Sig + 0.5 * M * (1 + log(2 * pi)))
      H_sig2 <- sum(entropy_inv_gamma(a_t, d_t))
      const_IW_q <- (nu_t/2) * logdet_spd(S_t) - (nu_t * p/2) * log(2) - log_multigamma(nu_t/2, p)
      E_log_q_S <- const_IW_q - ((nu_t + p + 1)/2) * E_logdet_Sigma_IW(nu_t, S_t, p) - 0.5 * sum(diag(S_t %*%
          E_Sinv_IW(nu_t, S_t)))
      H_S <- -E_log_q_S
      H_rho <- -sum(w_rho[rho_idx] * log(w_rho[rho_idx]))
      L_like + L_prior_b + L_prior_sig2 + L_prior_S + L_prior_rho + H_b + H_sig2 + H_S + H_rho
  }

  .normalize_probability <- function(w, name = "probability vector") {
      w <- as.numeric(w)
      if (!length(w) || any(!is.finite(w)) || any(w < 0)) {
          stop(name, " must contain finite nonnegative values.", call. = FALSE)
      }
      sw <- sum(w)
      if (!is.finite(sw) || sw <= 0)
          stop(name, " has zero total mass.", call. = FALSE)
      w/sw
  }

  .symmetrize <- function(M) 0.5 * (M + t(M))

  .validate_vi_prior <- function(prior, C) {
      if (is.null(prior))
          prior <- make_vi_prior(C)
      required <- c("a0", "d0", "v0", "S0")
      if (!is.list(prior) || !all(required %in% names(prior))) {
          stop("prior must be a list containing a0, d0, v0, and S0.")
      }
      prior$a0 <- rep(as.numeric(prior$a0), length.out = C)
      prior$d0 <- rep(as.numeric(prior$d0), length.out = C)
      if (any(!is.finite(prior$a0)) || any(prior$a0 <= 0) || any(!is.finite(prior$d0)) || any(prior$d0 <=
          0)) {
          stop("prior a0 and d0 must be finite and positive.")
      }
      prior$v0 <- as.numeric(prior$v0)
      if (length(prior$v0) != 1L || !is.finite(prior$v0) || prior$v0 <= 3) {
          stop("prior v0 must be one finite value > 3.")
      }
      .assert_matrix(prior$S0, "prior$S0", expected_nrow = 2L, expected_ncol = 2L, symmetric = TRUE)
      .assert_spd(prior$S0, "prior$S0")
      prior
  }

  E_Sinv_IW <- function(nu, S) {
      nu * inv_spd(S)
  }

  E_inv_sigma2_IG <- function(a, b) {
      a/b
  }

  E_log_sigma2_IG <- function(a, b) {
      log(b) - digamma(a)
  }

  E_logdet_Sigma_IW <- function(nu, S, p = 2) {
      logdet_spd(S) - sum(digamma((nu + 1 - (1:p))/2)) - p * log(2)
  }

  S_sufficient_halves_fast <- function(y, m_mat, V_mat, R) {
      C <- nrow(y)
      M <- ncol(y)
      J <- C/2
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

  S_sufficient_halves_vectorized <- function(y, m_mat, V_mat, R) {
      .assert_matrix(y, "y")
      C <- nrow(y)
      M <- ncol(y)
      if (C%%2L != 0L)
          stop("nrow(y) must be even.")
      .assert_matrix(m_mat, "m_mat", expected_nrow = C, expected_ncol = M)
      .assert_matrix(V_mat, "V_mat", expected_nrow = C, expected_ncol = M)
      if (any(V_mat < -1e-12))
          stop("V_mat contains negative variances.")
      J <- C/2L
      .assert_matrix(R, "R", expected_nrow = J, expected_ncol = J, symmetric = TRUE)
      residual <- y - m_mat
      eL <- residual[seq_len(J), , drop = FALSE]
      eR <- residual[J + seq_len(J), , drop = FALSE]
      ReL <- R %*% eL
      ReR <- R %*% eR
      diag_R <- diag(R)
      a11 <- sum(eL * ReL) + sum(diag_R * rowSums(V_mat[seq_len(J), , drop = FALSE]))
      a22 <- sum(eR * ReR) + sum(diag_R * rowSums(V_mat[J + seq_len(J), , drop = FALSE]))
      a12 <- sum(eL * ReR)
      matrix(c(a11, a12, a12, a22), 2L, 2L)
  }

  build_A_from_Y_halves <- function(Y, k = 2) {
      build_A_from_Y_halves_knn(Y, k = k)
  }

  build_A_from_Y_halves_knn <- function(Y, k = 2) {
      C <- nrow(Y)
      J <- C/2
      stopifnot(C%%2 == 0)
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
          ord <- cand[order(W[i, cand], decreasing = TRUE)]
          nn <- ord[1:min(k, J - 1)]
          A[i, nn] <- W[i, nn]
      }
      A <- pmax(A, t(A))
      diag(A) <- 0
      normalize_A_spectral(A)
  }

  entropy_inv_gamma <- function(a, b) {
      a + log(b) + lgamma(a) - (1 + a) * digamma(a)
  }

  increment_map <- function(delta, ctx, x, base, tilt, full = FALSE) {
      C <- ctx$C
      ii <- 1:C
      hh <- C + 1:3
      ff <- C + 3L + seq_len(ctx$k)
      dtau <- delta[ii]
      dH <- lr_mat(delta[hh])
      dr <- delta[ctx$p]
      dQ <- -dr * ctx$A
      dO <- kronecker(dH, base$Q) + kronecker(base$H, dQ) + kronecker(dH, dQ)
      rhs <- dO %*% base$residual - (dtau * base$m)/rep(ctx$lambda, each = C)
      dm <- base$solveP(rhs)
      ok <- FALSE
      for (i in 1:40) {
          next_m <- base$solveP(rhs - dO %*% dm - (dtau * dm)/rep(ctx$lambda, each = C))
          error <- max(abs(next_m - dm))
          scale <- max(abs(next_m), 1e-300)
          dm <- next_m
          if (error <= 32 * .Machine$double.eps * scale) {
              ok <- TRUE
              break
          }
      }
      if (!ok)
          stop("Increment latent-mean resolvent did not converge")
      delta_p <- outer(dtau, 1/ctx$lambda) + diag(dO)
      dv <- -base$v^2 * delta_p/(1 + base$v * delta_p)
      dC <- -tcrossprod(dm, base$residual) - tcrossprod(base$residual, dm) + tcrossprod(dm) + diag(rowSums(dv))
      dd <- as.numeric((base$m * dm + 0.5 * (dm * dm + dv)) %*% (1/ctx$lambda))
      st <- matrix(0, 2, 2)
      for (k in 1:3) st <- st + tilt[hh[k]] * lr_tilt_scale(k)
      dS <- lr_sym(lr_contract(base$Q, dC) + lr_contract(dQ, base$Cmat) + lr_contract(dQ, dC) - st)
      Si <- base$Hnext/ctx$nu
      dHnext <- -ctx$nu * Si %*% dS %*% solve(diag(2) + Si %*% dS, Si)
      dtau_next <- (-tilt[ii] * base$d - base$a * dd)/(base$d * (base$d + dd))
      dSA <- lr_contract(ctx$A, dC)
      dk <- 0.5 * (sum(dH * base$SA) + sum(base$H * dSA) + sum(dH * dSA))
      dl <- ctx$grid * dk + as.numeric(ctx$features %*% tilt[ff])
      anchor <- which.max(base$w)
      dl <- dl - dl[anchor]
      em <- expm1(dl)
      cc <- sum(base$w * em)
      dw <- base$w * (em - cc)/(1 + cc)
      features <- sweep(ctx$features, 2, ctx$features[anchor, ], "-")
      dfeatures <- as.numeric(crossprod(dw, features))
      out <- c(dtau_next, lr_vech(dHnext), dfeatures[ctx$k])
      if (!full)
          return(out)
      list(state_change = out, output_features_change = dfeatures, latent_delta = dm, latent_iterations = i,
          weight_change = dw)
  }

  increment_output_difference <- function(plus, minus, ctx, x, tp, tm) {
      C <- ctx$C
      ap <- ctx$a - tp[1:C]
      am <- ctx$a - tm[1:C]
      dlog <- log1p((ap - am)/am) - log1p((plus$delta[1:C] - minus$delta[1:C])/(x[1:C] + minus$delta[1:C])) -
          (digamma(ap) - digamma(am))
      p <- increment_map(plus$delta, ctx, x, plus$base, tp, TRUE)
      m <- increment_map(minus$delta, ctx, x, minus$base, tm, TRUE)
      c(dlog, plus$delta[C + 1:3] - minus$delta[C + 1:3], p$output_features_change - m$output_features_change)
  }

  increment_refine <- function(initial, ctx, tilt, numeric_jacobian = FALSE) {
      x <- initial
      b <- lr_map(x, ctx, full = TRUE)
      p <- ctx$p
      sc <- pmax(1, abs(x))
      start <- clock()
      if (numeric_jacobian) {
          J <- matrix(0, p, p)
          eps <- 2e-05
          for (k in 1:p) {
              xp <- xm <- x
              xp[k] <- xp[k] + eps * sc[k]
              xm[k] <- xm[k] - eps * sc[k]
              J[, k] <- (lr_map(xp, ctx) - lr_map(xm, ctx))/(2 * eps * sc[k])
          }
      }
      else J <- lr_jacobian(x, ctx, base = b)$jacobian
      M <- (diag(p) - J) * rep(sc, each = p)/sc
      invM <- solve(M)
      base_residual <- b$state - x
      delta <- numeric(p)
      history <- list()
      stagnant <- 0L
      status <- "maximum_updates"
      for (it in 1:cfg$revalidation$maximum_updates) {
          residual <- base_residual + increment_map(delta, ctx, x, b, tilt) - delta
          step <- sc * as.numeric(invM %*% (residual/sc))
          candidate <- delta + step
          residual2 <- base_residual + increment_map(candidate, ctx, x, b, tilt) - candidate
          err <- max(abs(residual2)/sc)
          relative_delta <- max(abs(residual2)/sc)/max(max(abs(candidate)/sc), 1e-300)
          history[[it]] <- data.frame(update = it, residual_before = max(abs(residual)/sc), residual_after = err,
              step_scaled = max(abs(step)/sc), rate = 1)
          stagnant <- if (identical(candidate, delta))
              stagnant + 1L
          else 0L
          delta <- candidate
          if (it >= cfg$revalidation$minimum_updates && relative_delta < 1e-11) {
              status <- "increment_precision_reached"
              break
          }
          if (stagnant >= cfg$revalidation$stagnation_window) {
              status <- "increment_stagnation"
              break
          }
          if (clock() - start > cfg$revalidation$per_solve_timeout_seconds) {
              status <- "solver_time_limit"
              break
          }
      }
      final <- x + delta
      list(x = final, delta = delta, base = b, updates = it, iterations = it, initial_residual = lr_relative_residual(x,
          lr_map(x, ctx, tilt)), residual = lr_relative_residual(final, lr_map(final, ctx, tilt)), increment_residual = err,
          relative_increment_residual = relative_delta, displacement_norm = sqrt(sum(delta^2)), unchanged = all(delta ==
              0), status = status, history = do.call(rbind, history), seconds = clock() - start)
  }

  inv_spd <- function(M, jitter = 1e-08, max_tries = 6) {
      p <- nrow(M)
      for (tt in 0:max_tries) {
          Mt <- 0.5 * (M + t(M)) + diag(jitter * (10^tt), p)
          co <- tryCatch(chol(Mt), error = function(e) NULL)
          if (!is.null(co))
              return(chol2inv(co))
      }
      stop("inv_spd: matrix not SPD even after jitter.")
  }

  log_multigamma <- function(a, p) {
      (p * (p - 1)/4) * log(pi) + sum(lgamma(a + (1 - (1:p))/2))
  }

  logdet_spd <- function(M, jitter = 1e-08, max_tries = 6) {
      p <- nrow(M)
      for (tt in 0:max_tries) {
          Mt <- 0.5 * (M + t(M)) + diag(jitter * (10^tt), p)
          co <- tryCatch(chol(Mt), error = function(e) NULL)
          if (!is.null(co))
              return(2 * sum(log(diag(co))))
      }
      stop("logdet_spd: matrix not SPD even after jitter.")
  }

  lr_basis <- function(k) {
      if (k == 1L)
          return(matrix(c(1, 0, 0, 0), 2))
      if (k == 2L)
          return(matrix(c(0, 1, 1, 0), 2))
      matrix(c(0, 0, 0, 1), 2)
  }

  lr_context <- function(Y, Kpack, A, D, fit) {
      C <- nrow(Y)
      n <- ncol(Y)
      J <- C/2L
      stopifnot(C == 2L * J, length(Kpack$lam) == n, all(Kpack$lam > 0), max(abs(D - diag(rowSums(A)))) <
          1e-10)
      prior <- fit$prior
      if (is.null(prior))
          stop("Original prior must be supplied explicitly")
      prior$a0 <- rep(prior$a0, length.out = C)
      prior$d0 <- rep(prior$d0, length.out = C)
      rg <- fit$rho_grid
      pw <- fit$rho_prior
      if (is.null(pw))
          pw <- rep(1/length(rg), length(rg))
      qi <- vector("list", length(rg))
      ld <- rep(NA_real_, length(rg))
      ok <- rep(FALSE, length(rg))
      for (k in seq_along(rg)) {
          ch <- tryCatch(chol(D - rg[k] * A), error = function(e) NULL)
          if (!is.null(ch) && pw[k] > 0) {
              ok[k] <- TRUE
              qi[[k]] <- diag(chol2inv(ch))
              ld[k] <- 2 * sum(log(diag(ch)))
          }
      }
      stopifnot(!any(fit$rho_w[!ok] > 0))
      qdiag <- do.call(rbind, qi[ok])
      grid <- rg[ok]
      pw <- pw[ok]/sum(pw[ok])
      features <- cbind(log(qdiag), rho = grid)
      list(C = C, n = n, J = J, Z = Y %*% Kpack$U, lambda = Kpack$lam, A = A, D = D, prior = prior, a = prior$a0 +
          n/2, nu = prior$v0 + n * J, grid = grid, full_grid = rg, valid = ok, prior_w = pw, logdet = ld[ok],
          qdiag = qdiag, features = features, p = C + 4L, k = ncol(features), target_dim = C + 3L + ncol(features),
          iw_numerical_ridge = 1e-08)
  }

  lr_contract <- function(Q, Cmat) {
      j <- nrow(Q)
      out <- matrix(0, 2, 2)
      for (s in 1:2) for (t in 1:2) {
          ii <- (s - 1L) * j + seq_len(j)
          jj <- (t - 1L) * j + seq_len(j)
          out[s, t] <- sum(Q * Cmat[ii, jj, drop = FALSE])
      }
      lr_sym(out)
  }

  lr_covariance <- function(x, ctx) {
      C <- ctx$C
      p <- ctx$p
      kf <- ctx$k
      q <- ctx$target_dim
      ii <- seq_len(C)
      hh <- C + 1:3
      ff <- C + 3L + seq_len(kf)
      jt <- lr_jacobian(x, ctx)
      b <- jt$base
      Jmap <- jt$jacobian
      vf <- b$feature_covariance
      cr <- vf[, ctx$k]
      B <- matrix(0, p, q)
      B[cbind(ii, ii)] <- -1/b$d
      for (j in 1:3) B[hh, hh[j]] <- lr_vech(b$Hnext %*% lr_tilt_scale(j) %*% b$Hnext/ctx$nu)
      B[p, ff] <- cr
      G <- matrix(0, q, p)
      G[cbind(ii, ii)] <- -1/x[ii]
      G[cbind(hh, hh)] <- 1
      G[ff, ] <- outer(cr, jt$dkappa)
      direct <- matrix(0, q, q)
      direct[cbind(ii, ii)] <- trigamma(ctx$a) - 1/ctx$a
      direct[ff, ff] <- vf
      sc <- pmax(1, abs(x))
      lhs <- (diag(p) - Jmap) * rep(sc, each = p)/sc
      dx <- sc * solve(lhs, B/sc)
      raw <- G %*% dx + direct
      cov <- lr_sym(raw)
      diag_scale <- sqrt(pmax(abs(diag(cov)), 1e-24))
      asym_scaled <- max(abs((raw - t(raw))/outer(diag_scale, diag_scale)))
      eig <- eigen(cov, symmetric = TRUE, only.values = TRUE)$values
      residual <- max(abs((diag(p) - Jmap) %*% dx - B))/max(1, max(abs(B)))
      list(cov = cov, raw_cov = raw, dxdt = dx, base = b, jacobian = Jmap, output_derivative = G, direct_response = direct,
          perturbation_derivative = B, dkappa = jt$dkappa, scaled_asymmetry = asym_scaled, min_eigenvalue = min(eig),
          max_eigenvalue = max(eig), scaled_system_rcond = rcond(lhs), linear_solve_residual = residual,
          fixed_point_residual = lr_relative_residual(x, b$state))
  }

  lr_feasible <- function(x, ctx) {
      if (length(x) != ctx$p || any(!is.finite(x)) || any(x[seq_len(ctx$C)] <= 0))
          return(FALSE)
      H <- lr_mat(x[ctx$C + 1:3])
      r <- x[ctx$p]
      if (r < min(ctx$grid) - 1e-12 || r > max(ctx$grid) + 1e-12)
          return(FALSE)
      !is.null(tryCatch(chol(H), error = function(e) NULL)) && !is.null(tryCatch(chol(ctx$D - r * ctx$A),
          error = function(e) NULL))
  }

  lr_feature_moments <- function(ctx, w) {
      anchor <- ctx$features[which.max(w), ]
      shifted <- sweep(ctx$features, 2, anchor, "-")
      ms <- as.numeric(crossprod(w, shifted))
      centered <- sweep(shifted, 2, ms, "-")
      list(mean = anchor + ms, centered = centered, cov = crossprod(centered, w * centered))
  }

  lr_intervals <- function(lr, x, ctx, archived_h) {
      C <- ctx$C
      J <- ctx$J
      h <- as.numeric(archived_h)
      if (any(!is.finite(h) | h <= 0 | h >= 1))
          stop("Original point estimate outside open unit interval; do not clip")
      if (lr$min_eigenvalue < -1e-08 * max(1, lr$max_eigenvalue))
          stop("LRVB covariance not PSD; no eigenvalue repair")
      H <- lr_mat(x[C + 1:3])
      Sstar <- solve(H)
      L <- matrix(0, C, ctx$target_dim)
      L[cbind(seq_len(C), seq_len(C))] <- 1
      for (c in seq_len(C)) {
          side <- if (c <= J)
              1L
          else 2L
          j <- (c - 1L)%%J + 1L
          L[c, C + 1:3] <- c(Sstar[side, 1]^2, 2 * Sstar[side, 1] * Sstar[side, 2], Sstar[side, 2]^2)/Sstar[side,
              side]
          L[c, C + 3L + j] <- -1
      }
      vz <- lr_sym(L %*% lr$cov %*% t(L))
      if (any(diag(vz) < 0))
          stop("Negative logit variance; do not clip")
      dh <- h * (1 - h)
      vh <- vz * outer(dh, dh)
      W <- cbind(diag(J), diag(J))/2
      pair_h <- as.numeric(W %*% h)
      pair_v <- lr_sym(W %*% vh %*% t(W))
      pvz <- pair_v/outer(pair_h * (1 - pair_h), pair_h * (1 - pair_h))
      if (any(diag(pvz) < 0))
          stop("Negative bilateral-mean variance; do not clip")
      make <- function(center, v, target) {
          sd <- sqrt(diag(v))
          z <- qlogis(center)
          data.frame(target = target, index = seq_along(center), h2 = center, sd_logit = sd, lower = plogis(z -
              qnorm(0.975) * sd), upper = plogis(z + qnorm(0.975) * sd))
      }
      list(intervals = rbind(make(h, vz, "trait"), make(pair_h, pvz, "bilateral_mean")), logit_covariance = vz,
          h_covariance = vh, bilateral_h_covariance = pair_v, delta_jacobian = L, point_estimates_unchanged = identical(h,
              as.numeric(archived_h)))
  }

  lr_jacobian <- function(x, ctx, tilt = rep(0, ctx$target_dim), base = NULL) {
      if (is.null(base))
          base <- lr_map(x, ctx, tilt, TRUE)
      C <- ctx$C
      p <- ctx$p
      n <- ctx$n
      ii <- seq_len(C)
      jac <- matrix(0, p, p)
      dkappa <- numeric(p)
      rho_var <- base$feature_covariance[ctx$k, ctx$k]
      for (k in seq_len(p)) {
          dtau <- numeric(C)
          dH <- matrix(0, 2, 2)
          dr <- 0
          if (k <= C)
              dtau[k] <- 1
          else if (k <= C + 3L)
              dH <- lr_basis(k - C)
          else dr <- 1
          dQ <- -dr * ctx$A
          dOmega <- kronecker(dH, base$Q) + kronecker(base$H, dQ)
          dm <- base$solveP(dOmega %*% base$residual - (dtau * base$m)/rep(ctx$lambda, each = C))
          dv <- -base$v^2 * (outer(dtau, 1/ctx$lambda) + diag(dOmega))
          dd <- as.numeric((base$m * dm + 0.5 * dv) %*% (1/ctx$lambda))
          dCmat <- -tcrossprod(dm, base$residual) - tcrossprod(base$residual, dm) + diag(rowSums(dv))
          dS <- lr_contract(dQ, base$Cmat) + lr_contract(base$Q, dCmat)
          jac[ii, k] <- -base$a/base$d^2 * dd
          jac[C + 1:3, k] <- lr_vech(-base$Hnext %*% dS %*% base$Hnext/ctx$nu)
          dkappa[k] <- 0.5 * (sum(dH * base$SA) + sum(base$H * lr_contract(ctx$A, dCmat)))
          jac[p, k] <- rho_var * dkappa[k]
      }
      list(jacobian = jac, dkappa = dkappa, base = base)
  }

  lr_map <- function(x, ctx, tilt = rep(0, ctx$target_dim), full = FALSE) {
      C <- ctx$C
      n <- ctx$n
      J <- ctx$J
      ii <- seq_len(C)
      hh <- C + 1:3
      ff <- C + 3L + seq_len(ctx$k)
      tau <- x[ii]
      H <- lr_mat(x[hh])
      r <- x[ctx$p]
      if (any(tau <= 0))
          stop("Nonpositive genetic precision")
      chol(H)
      Q <- ctx$D - r * ctx$A
      chol(Q)
      Omega <- kronecker(H, Q)
      root <- sqrt(tau)
      eg <- eigen(Omega/outer(root, root), symmetric = TRUE)
      if (any(eg$values <= 0))
          stop("Nonpositive transformed residual precision")
      den <- outer(eg$values, 1/ctx$lambda, "+")
      solveP <- function(rhs) (eg$vectors %*% (crossprod(eg$vectors, rhs/root)/den))/root
      m <- solveP(Omega %*% ctx$Z)
      v <- 1/(outer(tau, 1/ctx$lambda) + diag(Omega))
      residual <- ctx$Z - m
      Cmat <- tcrossprod(residual) + diag(rowSums(v))
      d <- ctx$prior$d0 + 0.5 * as.numeric((m * m + v) %*% (1/ctx$lambda))
      a <- ctx$a - tilt[ii]
      if (any(a <= 1))
          stop("Tilt outside finite-mean inverse-gamma family")
      scale_tilt <- matrix(0, 2, 2)
      for (k in 1:3) scale_tilt <- scale_tilt + tilt[hh[k]] * lr_tilt_scale(k)
      S <- lr_sym(ctx$prior$S0 - scale_tilt + lr_contract(Q, Cmat)) + diag(ctx$iw_numerical_ridge, 2)
      chol(S)
      Hnext <- ctx$nu * solve(S)
      SA <- lr_contract(ctx$A, Cmat)
      kappa <- 0.5 * sum(H * SA)
      logw <- log(ctx$prior_w) + ctx$n * ctx$logdet + ctx$grid * kappa + as.numeric(ctx$features %*% tilt[ff])
      w <- lr_softmax(logw)
      fm <- lr_feature_moments(ctx, w)
      rnext <- fm$mean[ctx$k]
      state <- c(a/d, lr_vech(Hnext), rnext)
      if (!full)
          return(state)
      list(state = state, a = a, d = d, H = H, Hnext = Hnext, Q = Q, Omega = Omega, m = m, v = v, residual = residual,
          Cmat = Cmat, S = S, SA = SA, kappa = kappa, w = w, logw = logw, solveP = solveP, mean_features = fm$mean,
          feature_covariance = fm$cov)
  }

  lr_mat <- function(x) matrix(c(x[1], x[2], x[2], x[3]), 2, 2)

  lr_output_expectation <- function(x, ctx, tilt = rep(0, ctx$target_dim)) {
      b <- lr_map(x, ctx, tilt, TRUE)
      a <- ctx$a - tilt[seq_len(ctx$C)]
      c(log(a/x[seq_len(ctx$C)]) - digamma(a), x[ctx$C + 1:3], b$mean_features)
  }

  lr_plugin_at_state <- function(x, ctx) {
      b <- lr_map(x, ctx, full = TRUE)
      genetic <- (ctx$a/x[seq_len(ctx$C)])/(ctx$a - 1)
      Sigma <- ctx$nu/(ctx$nu - 3) * solve(lr_mat(x[ctx$C + 1:3]))
      qd <- as.numeric(crossprod(b$w, ctx$qdiag))
      residual <- rep(diag(Sigma), each = ctx$J) * rep(qd, 2)
      genetic/(genetic + residual)
  }

  lr_refine <- function(initial, ctx, tilt = rep(0, ctx$target_dim), tol = 1e-10, max_iter = 20000L) {
      x <- initial
      trace <- list()
      for (it in seq_len(max_iter)) {
          b <- lr_map(x, ctx, tilt, TRUE)
          err <- lr_relative_residual(x, b$state)
          if (err < tol)
              return(list(x = x, iterations = it, residual = err, trace = trace))
          mode <- "damped_fixed_point"
          next_x <- 0.5 * x + 0.5 * b$state
          if (it == 1L || it%%5L == 0L) {
              jj <- lr_jacobian(x, ctx, tilt, b)$jacobian
              sc <- pmax(1, abs(x))
              lhs <- (diag(ctx$p) - jj) * rep(sc, each = ctx$p)/sc
              step <- tryCatch(sc * as.numeric(solve(lhs, (b$state - x)/sc)), error = function(e) NULL)
              if (!is.null(step) && all(is.finite(step)))
                  for (rate in 2^-(0:16)) {
                    cand <- x + rate * step
                    if (!lr_feasible(cand, ctx))
                      next
                    val <- tryCatch(lr_map(cand, ctx, tilt), error = function(e) NULL)
                    if (!is.null(val) && lr_relative_residual(cand, val) < err) {
                      next_x <- cand
                      mode <- paste0("Newton_", rate)
                      break
                    }
                  }
          }
          trace[[length(trace) + 1L]] <- c(iteration = it, error = err, mode = mode)
          x <- next_x
      }
      stop("Same-family numerical stationary-point check did not converge")
  }

  lr_relative_residual <- function(x, y) max(abs(y - x)/pmax(1, abs(x)))

  lr_softmax <- function(logw) {
      if (any(!is.finite(logw)))
          stop("Non-finite log weights on valid original rho support")
      x <- exp(logw - max(logw))
      x/sum(x)
  }

  lr_state <- function(fit) {
      c(as.numeric(fit$a/fit$d), lr_vech(fit$nu * solve(fit$S)), sum(fit$rho_grid * fit$rho_w))
  }

  lr_sym <- function(x) (x + t(x))/2

  lr_tilt_scale <- function(k) {
      if (k == 2L)
          lr_basis(k)
      else 2 * lr_basis(k)
  }

  lr_vech <- function(x) c(x[1, 1], x[1, 2], x[2, 2])

  make_K_eig_pack <- function(K, eps = 1e-06) {
      Kj <- 0.5 * (K + t(K)) + diag(eps, nrow(K))
      eg <- eigen(Kj, symmetric = TRUE)
      lam <- as.numeric(eg$values)
      U <- eg$vectors
      inv_lam <- 1/lam
      U2 <- U * U
      list(K = Kj, U = U, U2 = U2, lam = lam, inv_lam = inv_lam, logdetK = sum(log(lam)))
  }

  make_rho_precision_pack <- function(A, D_A, rho_grid, rho_prior = NULL, tol = 1e-10) {
      .assert_matrix(A, "A", symmetric = TRUE, tol = tol)
      J <- nrow(A)
      .assert_matrix(D_A, "D_A", expected_nrow = J, expected_ncol = J, symmetric = TRUE, tol = tol)
      rho_grid <- as.numeric(rho_grid)
      if (!length(rho_grid) || any(!is.finite(rho_grid)) || anyDuplicated(rho_grid)) {
          stop("rho_grid must contain distinct finite values.")
      }
      valid <- logical(length(rho_grid))
      logdet <- rep(-Inf, length(rho_grid))
      min_eigenvalue <- rep(NA_real_, length(rho_grid))
      for (ii in seq_along(rho_grid)) {
          Q <- .symmetrize(D_A - rho_grid[ii] * A)
          ev <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
          min_eigenvalue[ii] <- min(ev)
          co <- tryCatch(chol(Q), error = function(e) NULL)
          if (!is.null(co)) {
              valid[ii] <- TRUE
              logdet[ii] <- 2 * sum(log(diag(co)))
          }
      }
      if (!any(valid))
          stop("No SPD-valid rho grid points.")
      if (is.null(rho_prior)) {
          rho_prior <- as.numeric(valid)
      }
      else {
          rho_prior <- as.numeric(rho_prior)
          if (length(rho_prior) != length(rho_grid) || any(!is.finite(rho_prior)) || any(rho_prior < 0)) {
              stop("rho_prior must be a finite nonnegative vector matching rho_grid.")
          }
          if (any(rho_prior[!valid] > 0)) {
              stop("rho_prior assigns positive mass to a non-SPD grid point.")
          }
      }
      rho_prior[!valid] <- 0
      rho_prior <- .normalize_probability(rho_prior, "rho_prior")
      list(A = A, D_A = D_A, rho_grid = rho_grid, rho_prior = rho_prior, valid = valid, logdet = logdet,
          min_eigenvalue = min_eigenvalue)
  }

  make_vi_init_halves <- function(y, A, D_A, rho_grid, rho_prior, h2_guess = 0.15, ig_shape = 3, iw_df = 8,
      ridge = 1e-06) {
      C <- nrow(y)
      M <- ncol(y)
      J <- C/2
      stopifnot(C%%2 == 0)
      h2_guess <- min(max(h2_guess, 0.01), 0.8)
      ig_shape <- max(ig_shape, 2.05)
      iw_df <- max(iw_df, 3.2)
      a0 <- rep(2.1, C)
      prior_sigma2_mean <- rep(0.3, C)
      d0 <- prior_sigma2_mean * (a0 - 1)
      v0 <- 8
      Sigma0_mean <- diag(c(0.4, 0.4))
      S0 <- Sigma0_mean * (v0 - 3)
      trait_var <- apply(y, 1, var)
      fallback_var <- median(trait_var[is.finite(trait_var) & trait_var > ridge], na.rm = TRUE)
      if (!is.finite(fallback_var) || fallback_var <= ridge)
          fallback_var <- 1
      trait_var[!is.finite(trait_var) | trait_var <= ridge] <- fallback_var
      pair_var_base <- pmin(trait_var[1:J], trait_var[(J + 1):C])
      pair_var_base[!is.finite(pair_var_base) | pair_var_base <= ridge] <- fallback_var
      sigma2_init <- pmax(h2_guess * c(pair_var_base, pair_var_base), ridge)
      a_init <- rep(ig_shape, C)
      d_init <- sigma2_init * (a_init - 1)
      rho_w0 <- rho_prior
      R0 <- weighted_R_from_rho(A, D_A, rho_grid, rho_w0)
      Sigma_init <- pair_block_cov_halves(y, R0)
      Sigma_init <- 0.5 * (Sigma_init + t(Sigma_init))
      diag(Sigma_init) <- pmax(diag(Sigma_init), ridge)
      Sigma_init <- Sigma_init + diag(ridge, 2)
      nu_init <- iw_df
      S_init <- Sigma_init * (nu_init - 3)
      H_init <- E_Sinv_IW(nu_init, S_init)
      zero_mat <- matrix(0, C, M)
      rho_w0 <- update_q_rho_fast(y = y, m_mat = zero_mat, V_mat = zero_mat, A = A, D_A = D_A, rho_values = rho_grid,
          prior_w = rho_prior, H = H_init)$w
      R0 <- weighted_R_from_rho(A, D_A, rho_grid, rho_w0)
      Sigma_init <- pair_block_cov_halves(y, R0)
      Sigma_init <- 0.5 * (Sigma_init + t(Sigma_init))
      diag(Sigma_init) <- pmax(diag(Sigma_init), ridge)
      Sigma_init <- Sigma_init + diag(ridge, 2)
      S_init <- Sigma_init * (nu_init - 3)
      list(a0 = a0, d0 = d0, v0 = v0, S0 = S0, a_init = a_init, d_init = d_init, nu_init = nu_init, S_init = S_init,
          rho_w0 = rho_w0, sigma2_init = sigma2_init, Sigma_init = Sigma_init, trait_var = trait_var, pair_var_base = pair_var_base,
          h2_guess = h2_guess)
  }

  make_vi_prior <- function(C, iw_df = 8, iw_mean = diag(c(0.4, 0.4)), ig_shape = 2.1, ig_mean = 0.3) {
      C <- as.integer(C)
      if (length(C) != 1L || is.na(C) || C < 2L || C%%2L != 0L) {
          stop("C must be one positive even integer.")
      }
      if (!is.finite(iw_df) || iw_df <= 3) {
          stop("iw_df must exceed p + 1 = 3 so the IW mean exists.")
      }
      .assert_matrix(iw_mean, "iw_mean", expected_nrow = 2L, expected_ncol = 2L, symmetric = TRUE)
      .assert_spd(iw_mean, "iw_mean")
      if (!is.finite(ig_shape) || ig_shape <= 1 || !is.finite(ig_mean) || ig_mean <= 0) {
          stop("ig_shape must exceed 1 and ig_mean must be positive.")
      }
      list(a0 = rep(ig_shape, C), d0 = rep(ig_mean * (ig_shape - 1), C), v0 = as.numeric(iw_df), S0 = iw_mean *
          (iw_df - 3), iw_mean = iw_mean, iw_df = as.numeric(iw_df), ig_shape = ig_shape, ig_mean = ig_mean)
  }

  normalize_A_spectral <- function(A) {
      A <- as.matrix(A)
      J <- nrow(A)
      stopifnot(nrow(A) == ncol(A))
      A <- 0.5 * (A + t(A))
      diag(A) <- 0
      eps <- 1e-12
      deg <- rowSums(A)
      Dhi <- diag(1/sqrt(pmax(deg, eps)), J)
      B <- Dhi %*% A %*% Dhi
      lam_max <- max(eigen(B, symmetric = TRUE, only.values = TRUE)$values)
      if (!is.finite(lam_max) || lam_max <= 0)
          lam_max <- 1
      A/lam_max
  }

  pair_block_cov_halves <- function(y, R) {
      C <- nrow(y)
      M <- ncol(y)
      J <- C/2
      acc <- matrix(0, 2, 2)
      for (l in 1:M) {
          y2 <- vec_to_2xJ_halves(y[, l])
          acc <- acc + y2 %*% R %*% t(y2)
      }
      out <- acc/max(M * J, 1)
      0.5 * (out + t(out))
  }

  run_vi_custom <- function(Y, Kpack, A, D_A, rho_grid = seq(0, 0.99, by = 0.01), rho_prior = NULL, prior = NULL,
      max_iter = 1000L, tol = 1e-04, verbose = FALSE, init_h2 = 0.15, init_ig_shape = 3, init_iw_df = 8,
      fit_id = "spatial_vi_custom") {
      y <- as.matrix(Y)
      storage.mode(y) <- "double"
      .assert_matrix(y, "Y")
      C <- nrow(y)
      M <- ncol(y)
      if (C%%2L != 0L)
          stop("Y must have an even number of trait rows.")
      J <- C/2L
      .assert_matrix(A, "A", expected_nrow = J, expected_ncol = J, symmetric = TRUE)
      .assert_matrix(D_A, "D_A", expected_nrow = J, expected_ncol = J, symmetric = TRUE)
      required_k <- c("K", "U", "U2", "lam", "inv_lam", "logdetK")
      if (!is.list(Kpack) || !all(required_k %in% names(Kpack))) {
          stop("Kpack must be a v3 make_K_eig_pack() result.")
      }
      .assert_matrix(Kpack$K, "Kpack$K", expected_nrow = M, expected_ncol = M, symmetric = TRUE)
      .assert_matrix(Kpack$U, "Kpack$U", expected_nrow = M, expected_ncol = M)
      .assert_matrix(Kpack$U2, "Kpack$U2", expected_nrow = M, expected_ncol = M)
      if (length(Kpack$lam) != M || length(Kpack$inv_lam) != M || any(!is.finite(Kpack$lam)) || any(Kpack$lam <=
          0)) {
          stop("Kpack eigenvalue fields are invalid.")
      }
      max_iter <- as.integer(max_iter)
      if (length(max_iter) != 1L || is.na(max_iter) || max_iter < 1L) {
          stop("max_iter must be a positive integer.")
      }
      if (!is.finite(tol) || tol <= 0)
          stop("tol must be finite and positive.")
      prior <- .validate_vi_prior(prior, C)
      rho_pack <- make_rho_precision_pack(A, D_A, rho_grid, rho_prior)
      rho_grid <- rho_pack$rho_grid
      rho_prior <- rho_pack$rho_prior
      init <- make_vi_init_halves(y = y, A = A, D_A = D_A, rho_grid = rho_grid, rho_prior = rho_prior,
          h2_guess = init_h2, ig_shape = init_ig_shape, iw_df = init_iw_df)
      m_mat <- matrix(0, C, M)
      V_mat <- matrix(0, C, M)
      logdet_Sig <- rep(0, C)
      trace_invK_Sig <- rep(0, C)
      quad_total <- rep(0, C)
      a_t <- init$a_init
      d_t <- init$d_init
      nu_t <- init$nu_init
      S_t <- init$S_init
      rho_w <- init$rho_w0
      elbo_hist <- numeric(max_iter)
      beta_vec_final <- rep(NA_real_, C)
      converged <- FALSE
      start <- proc.time()
      for (iter in seq_len(max_iter)) {
          R_bar <- D_A - sum(rho_grid * rho_w) * A
          .assert_spd(R_bar, "posterior-mean pair precision")
          H <- E_Sinv_IW(nu_t, S_t)
          qb <- update_q_b_fast(y, Kpack, R_bar, H, m_mat, a_t/d_t)
          m_mat <- qb$m_mat
          V_mat <- qb$V_mat
          logdet_Sig <- qb$logdet_Sig
          trace_invK_Sig <- qb$trace_invK_Sig
          quad_total <- qb$quad_total
          beta_vec_final <- qb$beta_vec
          a_t <- prior$a0 + M/2
          d_t <- prior$d0 + 0.5 * quad_total
          if (any(!is.finite(d_t)) || any(d_t <= 0))
              stop("Non-finite q(sigma_g^2) update.")
          nu_t <- prior$v0 + M * J
          S_t <- .symmetrize(prior$S0 + S_sufficient_halves_vectorized(y, m_mat, V_mat, R_bar)) + diag(1e-08,
              2L)
          .assert_spd(S_t, "q(Sigma) scale")
          qrho <- update_q_rho_vectorized(y, m_mat, V_mat, A, D_A, rho_grid, rho_prior, E_Sinv_IW(nu_t,
              S_t), rho_pack = rho_pack)
          rho_w <- qrho$w
          Lcur <- .elbo_core_custom(y, A, D_A, m_mat, V_mat, logdet_Sig, quad_total, prior, nu_t, S_t,
              rho_grid, rho_w, rho_pack, Kpack$logdetK, a_t, d_t)
          if (!is.finite(Lcur))
              stop("ELBO became non-finite at iteration ", iter, ".")
          elbo_hist[iter] <- Lcur
          if (verbose)
              cat(sprintf("Iter %4d | ELBO = %.6f\n", iter, Lcur))
          if (iter >= 2L && abs(elbo_hist[iter] - elbo_hist[iter - 1L]) < tol) {
              converged <- TRUE
              break
          }
      }
      elapsed <- proc.time() - start
      elbo_hist <- elbo_hist[seq_len(iter)]
      E_Sigma <- S_t/(nu_t - 3)
      E_sigma2 <- as.numeric(d_t)/(as.numeric(a_t) - 1)
      out <- list(elbo = elbo_hist, mu = m_mat, Sig = NULL, Sig_diag = V_mat, Sig_logdet = logdet_Sig,
          Sig_trace_invK = trace_invK_Sig, a = a_t, d = d_t, a0 = prior$a0, d0 = prior$d0, nu = nu_t, S = S_t,
          E_Sigma = E_Sigma, E_sigma2 = E_sigma2, rho_grid = rho_grid, rho_w = rho_w, E_rho = sum(rho_grid *
              rho_w), A = A, D_A = D_A, model = "hetero", Kpack = Kpack, final_alpha = a_t/d_t, final_beta = beta_vec_final,
          init = init, prior = prior, rho_prior = rho_prior, rho_valid = rho_pack$valid, rho_min_eigenvalue = rho_pack$min_eigenvalue,
          fit_id = fit_id, converged = converged, iterations = iter, elbo_tail_delta = if (iter >= 2L) elbo_hist[iter] -
              elbo_hist[iter - 1L] else NA_real_, timing = c(user_seconds = unname(elapsed[["user.self"]]),
              system_seconds = unname(elapsed[["sys.self"]]), cpu_seconds = unname(elapsed[["user.self"]] +
                  elapsed[["sys.self"]]), wall_seconds = unname(elapsed[["elapsed"]])), dimensions = c(C = C,
              J = J, n = M), trait_order = "block", source_v3 = "spatialHVI::manuscript_engine", source_v3_sha256 = .reviewer_v3_sha256,
          source_library_sha256 = .reviewer_library_sha256)
      class(out) <- c("reviewer_vi_fit", "list")
      out
  }

  sample_iw2 <- function(nu, S) {
      wishart_draw <- stats::rWishart(1L, df = nu, Sigma = chol2inv(chol(S)))[, , 1L]
      chol2inv(chol(wishart_draw))
  }

  update_q_b_fast <- function(y, Kpack, R, H, m_mat, alpha_vec) {
      C <- nrow(y)
      n <- ncol(y)
      G <- kronecker(H, R)
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
          beta <- beta_vec[c]
          rhs <- as.numeric(Gy[c, ] - (g_row %*% m_mat - beta * m_mat[c, ]))
          u_rhs <- as.numeric(crossprod(U, rhs))
          denom <- alpha_vec[c] * inv_lam + beta
          d <- 1/denom
          u_m <- d * u_rhs
          m_c <- as.numeric(U %*% u_m)
          m_mat[c, ] <- m_c
          V_mat[c, ] <- as.numeric(U2 %*% d)
          logdet_Sig[c] <- sum(log(d))
          trace_invK_Sig[c] <- sum(inv_lam * d)
          quad_mean[c] <- sum(inv_lam * (u_m^2))
          quad_total[c] <- quad_mean[c] + trace_invK_Sig[c]
      }
      list(m_mat = m_mat, V_mat = V_mat, logdet_Sig = logdet_Sig, trace_invK_Sig = trace_invK_Sig, quad_mean = quad_mean,
          quad_total = quad_total, beta_vec = beta_vec)
  }

  update_q_rho_fast <- function(y, m_mat, V_mat, A, D_A, rho_values, prior_w, H) {
      M <- ncol(y)
      logw <- rep(-Inf, length(rho_values))
      for (tt in seq_along(rho_values)) {
          if (!is.finite(prior_w[tt]) || prior_w[tt] <= 0)
              next
          rho <- rho_values[tt]
          D_t <- D_A - rho * A
          D_t <- 0.5 * (D_t + t(D_t))
          co <- tryCatch(chol(D_t), error = function(e) NULL)
          if (is.null(co))
              next
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
      }
      else {
          mlp <- max(logw[finite_idx])
          w <- exp(logw - mlp)
          w[!is.finite(w)] <- 0
          s <- sum(w)
          if (!is.finite(s) || s <= 0) {
              w <- prior_w
          }
          else {
              w <- w/s
          }
      }
      list(w = w, logw = logw)
  }

  update_q_rho_vectorized <- function(y, m_mat, V_mat, A, D_A, rho_values, prior_w, H, rho_pack = NULL) {
      .assert_matrix(y, "y")
      C <- nrow(y)
      M <- ncol(y)
      if (C%%2L != 0L)
          stop("nrow(y) must be even.")
      J <- C/2L
      .assert_matrix(H, "H", expected_nrow = 2L, expected_ncol = 2L, symmetric = TRUE)
      .assert_spd(H, "H")
      if (is.null(rho_pack)) {
          rho_pack <- make_rho_precision_pack(A, D_A, rho_values, prior_w)
      }
      else {
          if (!isTRUE(all.equal(as.numeric(rho_values), rho_pack$rho_grid, tolerance = 0, check.attributes = FALSE))) {
              stop("rho_pack grid does not match rho_values.")
          }
          if (!isTRUE(all.equal(as.numeric(prior_w), rho_pack$rho_prior, tolerance = 1e-14, check.attributes = FALSE))) {
              stop("rho_pack prior does not match prior_w.")
          }
      }
      S_D <- S_sufficient_halves_vectorized(y, m_mat, V_mat, D_A)
      S_A <- S_sufficient_halves_vectorized(y, m_mat, V_mat, A)
      trace_HD <- sum(diag(H %*% S_D))
      trace_HA <- sum(diag(H %*% S_A))
      logw <- rep(-Inf, length(rho_values))
      ok <- which(rho_pack$valid & prior_w > 0)
      logw[ok] <- M * rho_pack$logdet[ok] - 0.5 * (trace_HD - rho_values[ok] * trace_HA) + log(prior_w[ok])
      finite <- which(is.finite(logw))
      if (!length(finite)) {
          w <- .normalize_probability(prior_w, "prior_w fallback")
      }
      else {
          w <- exp(logw - max(logw[finite]))
          w[!is.finite(w) | !rho_pack$valid] <- 0
          w <- .normalize_probability(w, "q(rho)")
      }
      if (any(w[!rho_pack$valid] != 0))
          stop("Invalid rho point received posterior mass.")
      list(w = w, logw = logw, S_D = S_D, S_A = S_A, trace_HD = trace_HD, trace_HA = trace_HA)
  }

  vec_to_2xJ_halves <- function(v) {
      C <- length(v)
      J <- C/2
      rbind(v[1:J], v[(J + 1):(2 * J)])
  }

  weighted_R_from_rho <- function(A, D_A, rho_grid, rho_w) {
      Reduce(`+`, lapply(seq_along(rho_grid), function(tt) {
          rho_w[tt] * (D_A - rho_grid[tt] * A)
      }))
  }

  weighted_Rinv_from_rho <- function(A, D_A, rho_grid, rho_w) {
      J <- nrow(A)
      acc <- matrix(0, J, J)
      for (tt in seq_along(rho_grid)) {
          w <- rho_w[tt]
          if (!is.finite(w) || w <= 0)
              next
          R_t <- 0.5 * (D_A - rho_grid[tt] * A + t(D_A - rho_grid[tt] * A))
          acc <- acc + w * inv_spd(R_t)
      }
      0.5 * (acc + t(acc))
  }
  environment()
})
