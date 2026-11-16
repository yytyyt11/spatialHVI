require_manuscript_fit <- function(fit) {
  if (!inherits(fit, "spatialhvi_fit"))
    stop("Use a fit returned by fit_spatial_hvi().", call. = FALSE)
  if (!isTRUE(fit$converged))
    stop("The MFVI fit has not converged; intervals are unavailable.", call. = FALSE)
  invisible(fit)
}

with_interval_seed <- function(seed, expression) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  })
  set.seed(seed)
  force(expression)
}

#' Original MFVI Heritability Intervals
#'
#' Draws jointly from the original variational factors, then forms each left/right
#' mean within the same draw. These are uncorrected variational intervals.
#' @param fit A converged object from [fit_spatial_hvi()].
#' @param draws Number of joint variational draws (1000 in the manuscript).
#' @param seed Random seed. The caller's RNG state is restored afterwards.
#' @param quantile_type R quantile convention: 7 for the simulation analyses;
#'   8 for the archived ADNI analysis.
#' @param level Central interval probability, default 0.95.
#' @return A data frame for individual traits and bilateral means. `h2` is the
#'   unchanged MFVI plug-in estimate, not the mean of posterior draws.
#' @export
vi_heritability_intervals <- function(fit, draws = 1000L, seed = 1L,
                                       quantile_type = 7L, level = 0.95) {
  require_manuscript_fit(fit)
  if (length(draws) != 1L || !is.finite(draws) || draws < 2 || draws != as.integer(draws))
    stop("draws must be an integer of at least two.")
  if (length(quantile_type) != 1L || !quantile_type %in% 1:9)
    stop("quantile_type must be an integer from 1 to 9.")
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be between zero and one.")
  C <- length(fit$E_sigma2); J <- C / 2L
  hdraw <- with_interval_seed(seed, {
    index <- sample.int(length(fit$rho_grid), draws, TRUE, prob = fit$rho_w)
    h <- matrix(NA_real_, draws, C)
    for (b in seq_len(draws)) {
      g <- 1 / stats::rgamma(C, shape = fit$a, rate = fit$d)
      Sigma <- .manuscript_engine$sample_iw2(fit$nu, fit$S)
      Q <- fit$D_A - fit$rho_grid[index[b]] * fit$A
      qi <- chol2inv(chol((Q + t(Q)) / 2))
      residual <- rep(diag(qi), 2) * rep(diag(Sigma), each = J)
      h[b, ] <- g / pmax(g + residual, 1e-12)
    }
    h
  })
  pair <- (hdraw[, seq_len(J), drop = FALSE] + hdraw[, J + seq_len(J), drop = FALSE]) / 2
  bounds <- t(apply(cbind(hdraw, pair), 2, stats::quantile,
    probs = c((1 - level) / 2, (1 + level) / 2), type = quantile_type, names = FALSE))
  point <- manuscript_heritability_summary(fit)
  data.frame(target = rep(c("trait", "bilateral_mean"), c(C, J)),
    index = c(seq_len(C), seq_len(J)), h2 = c(point$h2, point$bilateral_mean),
    lower = bounds[, 1], upper = bounds[, 2], method = "MFVI", level = level)
}

#' LRVB Covariance and Heritability Intervals
#'
#' Evaluates the original discrete-grid covariance response, including latent
#' genetic effects and covariance parameters. A stationary solution of the same
#' MFVI equations is used only to calculate covariance. Reported interval centers
#' remain the original plug-in estimates. Bilateral means use the full covariance.
#' Stops when numerical validity or stationary-point sensitivity checks fail;
#' no covariance eigenvalue repair or empirical inflation is performed here.
#' @param fit A converged object from [fit_spatial_hvi()].
#' @param verify_derivatives Also compare the analytic Jacobian to central finite
#'   differences and response columns to perturbed fixed-point refits.
#' @return A list with `intervals`, `diagnostics`, the trait and bilateral
#'   covariances, and the sufficient-statistic response covariance. The intervals
#'   are nominal 95 percent logit-delta intervals; a posterior covariance check
#'   does not establish repeated-sampling coverage.
#' @export
lrvb_heritability_intervals <- function(fit, verify_derivatives = FALSE) {
  require_manuscript_fit(fit)
  if (!is.logical(verify_derivatives) || length(verify_derivatives) != 1L || is.na(verify_derivatives))
    stop("verify_derivatives must be TRUE or FALSE.")
  e <- .manuscript_engine
  ctx <- e$lr_context(fit$input_y, fit$Kpack, fit$A, fit$D_A, fit)
  x0 <- e$lr_state(fit)
  rr <- e$lr_refine(x0, ctx, tol = 1e-10, max_iter = 20000L)
  lr <- e$lr_covariance(rr$x, ctx)
  h <- manuscript_heritability_summary(fit)$h2
  ci <- e$lr_intervals(lr, rr$x, ctx, h)
  lr0 <- e$lr_covariance(x0, ctx)
  ci0 <- e$lr_intervals(lr0, x0, ctx, h)
  shift <- e$lr_plugin_at_state(rr$x, ctx) - h
  diagnostics <- c(fixed_point_residual = lr$fixed_point_residual,
    scaled_asymmetry = lr$scaled_asymmetry, linear_solve_residual = lr$linear_solve_residual,
    minimum_covariance_eigenvalue = lr$min_eigenvalue,
    system_reciprocal_condition = lr$scaled_system_rcond,
    center_drift_in_posterior_sd = max(abs(shift) / sqrt(diag(ci$h_covariance))),
    relative_logit_variance_change = max(abs(diag(ci$logit_covariance) / diag(ci0$logit_covariance) - 1)),
    interval_endpoint_change = max(abs(ci$intervals$lower - ci0$intervals$lower),
                                   abs(ci$intervals$upper - ci0$intervals$upper)))
  thresholds <- c(fixed_point_residual = 1e-10, scaled_asymmetry = 1e-5,
    linear_solve_residual = 1e-8, center_drift_in_posterior_sd = 0.1,
    relative_logit_variance_change = 0.05, interval_endpoint_change = 0.002)
  bad <- names(thresholds)[!is.finite(diagnostics[names(thresholds)]) |
                           diagnostics[names(thresholds)] >= thresholds]
  if (length(bad)) stop("LRVB numerical checks failed: ", paste(bad, collapse = ", "), call. = FALSE)
  if (isTRUE(verify_derivatives)) {
    x <- rr$x; scale <- pmax(1, abs(x))
    qb <- e$update_q_b_fast(fit$input_y, fit$Kpack, lr$base$Q, lr$base$H,
                            lr$base$m %*% t(fit$Kpack$U), x[seq_len(ctx$C)])
    gaussian_mean_error <- max(abs(qb$m_mat - lr$base$m %*% t(fit$Kpack$U)))
    gaussian_variance_error <- max(abs(qb$V_mat - lr$base$v %*% t(fit$Kpack$U2)))
    if (max(gaussian_mean_error, gaussian_variance_error) >= 1e-8)
      stop("LRVB latent-effect elimination validation failed.")
    categorical_error <- numeric()
    for (j in c(1L, ctx$k)) {
      eps <- 1e-4
      plus <- e$lr_softmax(lr$base$logw + eps * ctx$features[, j])
      minus <- e$lr_softmax(lr$base$logw - eps * ctx$features[, j])
      derivative <- as.numeric(crossprod((plus - minus) / (2 * eps), ctx$features))
      categorical_error <- c(categorical_error, max(abs(derivative - lr$base$feature_covariance[, j])))
    }
    if (max(categorical_error) >= 1e-7) stop("Discrete rho response validation failed.")
    analytic <- lr$jacobian * rep(scale, each = ctx$p) / scale
    fd_error <- numeric()
    for (eps in c(2e-5, 1e-5)) {
      numeric_j <- matrix(0, ctx$p, ctx$p)
      for (j in seq_len(ctx$p)) {
        xp <- xm <- x; xp[j] <- xp[j] + eps * scale[j]; xm[j] <- xm[j] - eps * scale[j]
        numeric_j[, j] <- (e$lr_map(xp, ctx) - e$lr_map(xm, ctx)) / (2 * eps * scale[j])
      }
      numeric_j <- numeric_j * rep(scale, each = ctx$p) / scale
      fd_error <- c(fd_error, max(abs(numeric_j - analytic)) / max(1, max(abs(analytic))))
    }
    tilt_error <- numeric()
    columns <- unique(c(1L, ctx$J, ctx$C + 1L, ctx$C + 2L, ctx$C + 4L, ctx$target_dim))
    for (j in columns) for (eps in c(1e-3, 5e-4)) {
      tp <- tn <- numeric(ctx$target_dim); tp[j] <- eps; tn[j] <- -eps
      plus <- e$increment_refine(x, ctx, tp)
      minus <- e$increment_refine(x, ctx, tn)
      if (plus$status != "increment_precision_reached" || minus$status != "increment_precision_reached")
        stop("Perturbed fixed-point solve did not resolve the response.")
      numeric_response <- e$increment_output_difference(plus, minus, ctx, x, tp, tn) / (2 * eps)
      expected <- lr$raw_cov[, j]
      tilt_error <- c(tilt_error, sqrt(sum((numeric_response - expected)^2)) /
                       max(1e-7, sqrt(sum(expected^2))))
    }
    diagnostics <- c(diagnostics, finite_difference_error = max(fd_error), tilt_refit_error = max(tilt_error),
      gaussian_mean_elimination_error = gaussian_mean_error,
      gaussian_variance_elimination_error = gaussian_variance_error,
      categorical_derivative_error = max(categorical_error))
    if (any(!is.finite(c(fd_error, tilt_error))) || max(fd_error, tilt_error) >= 5e-4)
      stop("LRVB derivative validation failed.", call. = FALSE)
  }
  ci$intervals$method <- "LRVB_logit_delta"
  ci$intervals$level <- 0.95
  stopifnot(identical(ci$intervals$h2, c(h, (h[seq_len(ctx$J)] + h[ctx$J + seq_len(ctx$J)]) / 2)))
  structure(list(intervals = ci$intervals, diagnostics = diagnostics,
    h_covariance = ci$h_covariance, logit_covariance = ci$logit_covariance,
    bilateral_h_covariance = ci$bilateral_h_covariance,
    response_covariance = lr$cov, rho_variance = lr$cov[ctx$target_dim, ctx$target_dim],
    original_rho_grid = fit$rho_grid, refinement_iterations = rr$iterations,
    point_estimates_unchanged = TRUE), class = "spatialhvi_lrvb")
}

#' Fixed Interval Calibration from Manuscript Table S17
#'
#' These coefficients were learned from 80 development data sets before the
#' 160-data-set confirmation. Their performance is design-dependent. They are
#' provided for reproducing the manuscript; they are not universal 95 percent
#' coverage constants for a new study, prior, graph, or sample size.
#' @return A calibration list with separate trait and bilateral-mean coefficients.
#' @export
manuscript_interval_calibration <- function() {
  structure(list(coefficients = data.frame(target = c("trait", "bilateral_mean"),
    a = c(-2.8135915281262953, -3.7005660188107603),
    b = c(9.18743593573288, 13.273540095862245)), level = 0.95,
    source = "S17C; 80 independent development data sets; R quantile type 7",
    scope = "Manuscript S10 design; ADNI application is an uncertainty-sensitivity analysis"),
    class = "spatialhvi_calibration")
}

#' Learn Interval Calibration on Separate Development Simulations
#' @param development Data frame with `dataset_id`, `target` (trait or
#'   bilateral_mean), `h2`, `truth`, and `sd_logit`. Each data set must contribute
#'   the same number of targets within each target class, so type-7 empirical
#'   quantiles give equal data-set weights. Use development simulations only;
#'   do not learn these coefficients from validation data or observed ADNI traits.
#' @return A calibration list to freeze before independent validation.
#' @export
learn_interval_calibration <- function(development) {
  cols <- c("dataset_id", "target", "h2", "truth", "sd_logit")
  if (!is.data.frame(development) || !all(cols %in% names(development)) || !nrow(development))
    stop("development is missing required calibration columns.")
  d <- development
  if (anyNA(d[cols]) || any(!d$target %in% c("trait", "bilateral_mean")) ||
      any(!is.finite(d$h2) | !is.finite(d$truth) | !is.finite(d$sd_logit)) ||
      any(d$h2 <= 0 | d$h2 >= 1 | d$truth <= 0 | d$truth >= 1 | d$sd_logit <= 0))
    stop("Calibration needs finite interior estimates/truths and positive logit SDs.")
  out <- lapply(c("trait", "bilateral_mean"), function(target) {
    z <- d[d$target == target, ]
    counts <- table(z$dataset_id)
    if (length(counts) < 2L || length(unique(as.integer(counts))) != 1L)
      stop("Use at least two development data sets with equal target counts per class.")
    errors <- (stats::qlogis(z$h2) - stats::qlogis(z$truth)) / z$sd_logit
    q <- stats::quantile(errors, c(0.025, 0.975), type = 7, names = FALSE)
    data.frame(target = target, a = q[1], b = q[2])
  })
  structure(list(coefficients = do.call(rbind, out), level = 0.95,
    source = "User-supplied development simulations; R quantile type 7",
    development_dataset_ids = unique(d$dataset_id)), class = "spatialhvi_calibration")
}

#' Apply Frozen Simulation Calibration to LRVB Intervals
#' @param lrvb Result of [lrvb_heritability_intervals()].
#' @param calibration Explicit calibration from [learn_interval_calibration()] or
#'   [manuscript_interval_calibration()]. There is no automatic default.
#' @return A data frame with unchanged `h2` point estimates and asymmetric
#'   simulation-calibrated limits `plogis(logit(h2) - b * sd_logit)` and
#'   `plogis(logit(h2) - a * sd_logit)`.
#' @export
calibrate_heritability_intervals <- function(lrvb, calibration) {
  if (!inherits(lrvb, "spatialhvi_lrvb") || !inherits(calibration, "spatialhvi_calibration"))
    stop("Supply an LRVB result and an explicit frozen calibration.")
  x <- lrvb$intervals; coefficients <- calibration$coefficients
  if (anyDuplicated(coefficients$target)) stop("Calibration target classes must be unique.")
  match_target <- match(x$target, coefficients$target)
  a <- coefficients$a[match_target]; b <- coefficients$b[match_target]
  if (anyNA(match_target) || any(!is.finite(a) | !is.finite(b) | a > b) ||
      any(!is.finite(x$sd_logit) | x$sd_logit <= 0))
    stop("Missing or invalid calibration coefficients or LRVB standard deviations.")
  x$lower <- stats::plogis(stats::qlogis(x$h2) - b * x$sd_logit)
  x$upper <- stats::plogis(stats::qlogis(x$h2) - a * x$sd_logit)
  x$method <- "simulation_calibrated_LRVB"
  x$a <- a; x$b <- b
  x
}
