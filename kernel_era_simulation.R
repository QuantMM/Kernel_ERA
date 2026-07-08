# =============================================================================
# Simulation framework for Kernel ERA
# Requires: source("kernel_era_core.R")
# =============================================================================

# -----------------------------------------------------------------------------
# Data-generating mechanisms
# -----------------------------------------------------------------------------

.sim_standardize_train_apply <- function(z, train_index, tol = 1e-12) {
  z <- as.matrix(z)
  center <- colMeans(z[train_index, , drop = FALSE])
  zc_train <- sweep(z[train_index, , drop = FALSE], 2L, center, "-")
  scale <- sqrt(colMeans(zc_train * zc_train))
  if (any(scale <= tol)) stop("Degenerate true component in simulation.", call. = FALSE)
  sweep(sweep(z, 2L, center, "-"), 2L, scale, "/")
}

.sim_generate_predictors <- function(n,
                                     p_per_set,
                                     rho_within = 0.4,
                                     rho_between = 0.1) {
  p_per_set <- as.integer(p_per_set)
  K <- length(p_per_set)
  if (any(p_per_set < 1L)) stop("p_per_set must be positive.", call. = FALSE)
  if (!is.finite(rho_between) || !is.finite(rho_within) ||
      rho_between < 0 || rho_within < rho_between || rho_within >= 1) {
    stop("Require 0 <= rho_between <= rho_within < 1.", call. = FALSE)
  }

  global <- stats::rnorm(n)
  blocks <- vector("list", K)
  for (k in seq_len(K)) {
    block_latent <- stats::rnorm(n)
    noise <- matrix(stats::rnorm(n * p_per_set[k]), n, p_per_set[k])
    blocks[[k]] <- sqrt(rho_between) * global +
      sqrt(rho_within - rho_between) * block_latent +
      sqrt(1 - rho_within) * noise
    colnames(blocks[[k]]) <- paste0("X", k, "_", seq_len(p_per_set[k]))
  }
  blocks
}

.sim_get_col <- function(X, j) {
  X[, min(j, ncol(X))]
}

.sim_linear_component <- function(X) {
  q <- min(3L, ncol(X))
  w <- c(0.80, -0.60, 0.40)[seq_len(q)]
  as.numeric(X[, seq_len(q), drop = FALSE] %*% w)
}

.sim_curvature_interaction_component <- function(X) {
  x1 <- .sim_get_col(X, 1L)
  x2 <- .sim_get_col(X, 2L)
  x3 <- .sim_get_col(X, 3L)
  0.80 * (x1^2 - 1) +
    0.55 * (x2^2 - 1) +
    0.85 * x1 * x2 +
    0.35 * sin(1.3 * x3)
}

.sim_threshold_saturation_component <- function(X) {
  x1 <- .sim_get_col(X, 1L)
  x2 <- .sim_get_col(X, 2L)
  x3 <- .sim_get_col(X, 3L)
  0.90 * tanh(1.2 * x1) +
    0.75 * (as.numeric(x2 > 0) - 0.5) +
    0.55 * x1 * x3
}

.sim_radial_component <- function(X) {
  x1 <- .sim_get_col(X, 1L)
  x2 <- .sim_get_col(X, 2L)
  x3 <- .sim_get_col(X, 3L)
  radius2 <- x1^2 + x2^2
  1.10 * exp(-0.65 * radius2) + 0.45 * sin(x3)
}

.sim_smooth_component <- function(X) {
  x1 <- .sim_get_col(X, 1L)
  x2 <- .sim_get_col(X, 2L)
  x3 <- .sim_get_col(X, 3L)
  0.85 * sin(1.25 * x1) + 0.65 * tanh(x2) + 0.45 * cos(0.8 * x3)
}

.sim_true_components <- function(blocks,
                                 scenario = c("linear", "mixed", "nonlinear", "weak_nonlinear"),
                                 nonlinearity_strength = 0.35) {
  scenario <- match.arg(scenario)
  K <- length(blocks)
  out <- matrix(NA_real_, nrow(blocks[[1L]]), K)

  for (k in seq_len(K)) {
    linear <- .sim_linear_component(blocks[[k]])
    nonlinear <- switch(as.character(((k - 1L) %% 3L) + 1L),
      "1" = .sim_smooth_component(blocks[[k]]),
      "2" = .sim_curvature_interaction_component(blocks[[k]]),
      "3" = .sim_threshold_saturation_component(blocks[[k]])
    )

    out[, k] <- switch(scenario,
      linear = linear,
      mixed = {
        if (k == 1L) linear else if (k %% 3L == 2L) {
          .sim_curvature_interaction_component(blocks[[k]])
        } else {
          .sim_threshold_saturation_component(blocks[[k]])
        }
      },
      nonlinear = {
        if (k %% 3L == 1L) .sim_radial_component(blocks[[k]]) else nonlinear
      },
      weak_nonlinear = (1 - nonlinearity_strength) * linear +
        nonlinearity_strength * nonlinear
    )
  }
  colnames(out) <- paste0("true_component_", seq_len(K))
  out
}

simulate_kernel_era_data <- function(n_train = 300L,
                                     n_test = 1000L,
                                     p_per_set = c(5L, 5L, 5L),
                                     scenario = c("linear", "mixed", "nonlinear", "weak_nonlinear"),
                                     rho_within = 0.4,
                                     rho_between = 0.1,
                                     target_r2 = 0.50,
                                     beta = NULL,
                                     nonlinearity_strength = 0.35,
                                     seed = 1) {
  scenario <- match.arg(scenario)
  n_train <- as.integer(n_train)
  n_test <- as.integer(n_test)
  p_per_set <- as.integer(p_per_set)
  if (n_train < 20L || n_test < 1L) stop("Use n_train >= 20 and n_test >= 1.", call. = FALSE)
  if (target_r2 <= 0 || target_r2 >= 1) stop("target_r2 must lie in (0, 1).", call. = FALSE)

  set.seed(seed)
  n_total <- n_train + n_test
  blocks <- .sim_generate_predictors(n_total, p_per_set,
                                     rho_within = rho_within,
                                     rho_between = rho_between)
  X <- do.call(cbind, blocks)
  raw_components <- .sim_true_components(blocks, scenario,
                                         nonlinearity_strength = nonlinearity_strength)
  train_index <- seq_len(n_train)
  test_index <- n_train + seq_len(n_test)
  true_components <- .sim_standardize_train_apply(raw_components, train_index)

  K <- length(p_per_set)
  if (is.null(beta)) {
    beta <- seq(0.80, 0.50, length.out = K)
  }
  beta <- as.numeric(beta)
  if (length(beta) != K) stop("beta must have one value per predictor set.", call. = FALSE)

  signal <- as.numeric(true_components %*% beta)
  signal_variance <- mean((signal[train_index] - mean(signal[train_index]))^2)
  error_variance <- signal_variance * (1 - target_r2) / target_r2
  error_sd <- sqrt(error_variance)
  y <- signal + stats::rnorm(n_total, sd = error_sd)

  list(
    X_train = X[train_index, , drop = FALSE],
    y_train = y[train_index],
    X_test = X[test_index, , drop = FALSE],
    y_test = y[test_index],
    true_components_train = true_components[train_index, , drop = FALSE],
    true_components_test = true_components[test_index, , drop = FALSE],
    signal_train = signal[train_index],
    signal_test = signal[test_index],
    beta = beta,
    nvar = p_per_set,
    scenario = scenario,
    target_r2 = target_r2,
    error_sd = error_sd,
    rho_within = rho_within,
    rho_between = rho_between,
    seed = seed
  )
}

# -----------------------------------------------------------------------------
# Metrics and extraction
# -----------------------------------------------------------------------------

.sim_prediction_metrics <- function(y_train, y_test, pred_test) {
  y_train <- as.numeric(y_train)
  y_test <- as.numeric(y_test)
  pred_test <- as.numeric(pred_test)
  sse <- sum((y_test - pred_test)^2)
  denominator <- sum((y_test - mean(y_train))^2)
  data.frame(
    test_mse = mean((y_test - pred_test)^2),
    test_rmse = sqrt(mean((y_test - pred_test)^2)),
    test_mae = mean(abs(y_test - pred_test)),
    test_r2 = 1 - sse / denominator,
    prediction_correlation = .kera_safe_cor(y_test, pred_test),
    stringsAsFactors = FALSE
  )
}

.sim_component_metrics <- function(F_train_est,
                                   F_test_est,
                                   F_train_true,
                                   F_test_true,
                                   method) {
  K <- ncol(F_train_true)
  do.call(rbind, lapply(seq_len(K), function(k) {
    sign_k <- 1
    rr <- .kera_safe_cor(F_train_est[, k], F_train_true[, k])
    if (is.finite(rr) && rr < 0) sign_k <- -1
    est_train <- sign_k * F_train_est[, k]
    est_test <- sign_k * F_test_est[, k]

    # A least-squares affine calibration is used only for the RMSE metric;
    # correlation remains invariant to scale and location.
    calibration <- stats::lm.fit(cbind(1, est_train), F_train_true[, k])$coefficients
    if (any(!is.finite(calibration))) calibration <- c(0, 0)
    calibrated_test <- calibration[1L] + calibration[2L] * est_test

    data.frame(
      method = method,
      predictor_set = k,
      train_component_correlation = .kera_safe_cor(est_train, F_train_true[, k]),
      test_component_correlation = .kera_safe_cor(est_test, F_test_true[, k]),
      test_component_rmse_calibrated = sqrt(mean((calibrated_test - F_test_true[, k])^2)),
      orientation_sign = sign_k,
      stringsAsFactors = FALSE
    )
  }))
}

.sim_condition_columns <- function(condition) {
  out <- as.data.frame(condition, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

# -----------------------------------------------------------------------------
# One Monte Carlo replication
# -----------------------------------------------------------------------------

run_one_kernel_era_replication <- function(condition,
                                           replication,
                                           lambda_grid = 10^seq(-4, 2, by = 1),
                                           sigma_grid = c(0.75, 1.5, 3, 6, 12, 24, 48),
                                           cv_folds = 5L,
                                           cv_selection = c("one_se", "minimum"),
                                           cv_require_all_folds_converged = TRUE,
                                           gaussian_denominator = c("two_sigma_squared", "sigma_squared"),
                                           penalty = c("rkhs", "dual_l2"),
                                           final_n_starts = 3L,
                                           max_iter = 500L,
                                           tol = 1e-6,
                                           objective_tol = 1e-9,
                                           interpolation_fit_threshold = 0.98,
                                           bootstrap_reps = 0L,
                                           include_matlab_gcv = FALSE,
                                           verbose = FALSE) {
  cv_selection <- match.arg(cv_selection)
  gaussian_denominator <- match.arg(gaussian_denominator)
  penalty <- match.arg(penalty)
  cv_require_all_folds_converged <- isTRUE(cv_require_all_folds_converged)
  interpolation_fit_threshold <- as.numeric(interpolation_fit_threshold)
  if (length(interpolation_fit_threshold) != 1L ||
      !is.finite(interpolation_fit_threshold) ||
      interpolation_fit_threshold <= 0 ||
      interpolation_fit_threshold >= 1) {
    stop("interpolation_fit_threshold must lie strictly between 0 and 1.",
         call. = FALSE)
  }
  elapsed_start <- proc.time()[["elapsed"]]
  condition <- as.list(condition)
  condition_id <- as.integer(condition$condition_id %||% 1L)
  seed <- as.integer(condition$seed_base %||% 20260707L) +
    100000L * condition_id + as.integer(replication)

  p_value <- condition$p_per_set %||% 5L
  K <- as.integer(condition$n_sets %||% 3L)
  p_per_set <- if (length(p_value) == 1L) rep(as.integer(p_value), K) else as.integer(p_value)

  dat <- simulate_kernel_era_data(
    n_train = as.integer(condition$n_train %||% 300L),
    n_test = as.integer(condition$n_test %||% 1000L),
    p_per_set = p_per_set,
    scenario = condition$scenario %||% "mixed",
    rho_within = as.numeric(condition$rho_within %||% 0.4),
    rho_between = as.numeric(condition$rho_between %||% 0.1),
    target_r2 = as.numeric(condition$target_r2 %||% 0.5),
    nonlinearity_strength = as.numeric(condition$nonlinearity_strength %||% 0.35),
    seed = seed
  )

  linear_fit <- linear_era_fit(
    X = dat$X_train,
    Y = dat$y_train,
    nvar = dat$nvar,
    standardize = TRUE,
    init = "ones",
    seed = seed + 10L,
    max_iter = max_iter,
    tol = tol,
    verbose = FALSE
  )
  linear_fit <- align_linear_era_signs(linear_fit, dat$true_components_train)
  linear_pred <- predict(linear_fit, dat$X_test, type = "response")
  linear_F_test <- predict(linear_fit, dat$X_test, type = "components")

  cv <- kernel_era_cv(
    X = dat$X_train,
    Y = dat$y_train,
    nvar = dat$nvar,
    lambda_grid = lambda_grid,
    sigma_grid = sigma_grid,
    v = cv_folds,
    seed = seed + 20L,
    kernel = "gaussian",
    gaussian_denominator = gaussian_denominator,
    penalty = penalty,
    solver = "spectral",
    linear_init = "ones",
    max_iter = max_iter,
    tol = tol,
    objective_tol = objective_tol,
    selection = cv_selection,
    require_all_folds_converged = cv_require_all_folds_converged,
    verbose = FALSE
  )

  selected_cv_row <- cv$selected_row
  unrestricted_cv_row <- cv$unrestricted_minimum_row

  kernel_fit <- kernel_era_fit(
    X = dat$X_train,
    Y = dat$y_train,
    nvar = dat$nvar,
    lambda = cv$best_lambda,
    kernel = "gaussian",
    sigma = cv$best_sigma,
    gaussian_denominator = gaussian_denominator,
    penalty = penalty,
    solver = "spectral",
    standardize = TRUE,
    linear_init = "ones",
    n_starts = final_n_starts,
    seed = seed + 30L,
    max_iter = max_iter,
    tol = tol,
    objective_tol = objective_tol,
    verbose = FALSE
  )
  kernel_fit <- align_kernel_era_signs(kernel_fit, dat$true_components_train)
  kernel_pred <- predict(kernel_fit, dat$X_test, type = "response")
  kernel_F_test <- predict(kernel_fit, dat$X_test, type = "components")

  linear_test_metrics <-
    .sim_prediction_metrics(dat$y_train, dat$y_test, linear_pred)

  linear_perf <- cbind(
    data.frame(method = "Linear ERA", stringsAsFactors = FALSE),
    linear_test_metrics,
    data.frame(
      train_fit = linear_fit$fit,
      train_test_fit_gap = linear_fit$fit - linear_test_metrics$test_r2,
      converged = linear_fit$converged,
      iterations = linear_fit$iterations,
      lambda = NA_real_,
      sigma = NA_real_,
      cv_mse = NA_real_,
      cv_se = NA_real_,
      cv_convergence_rate = NA_real_,
      cv_n_successful_folds = NA_integer_,
      cv_eligible = NA,
      cv_selection_rule = NA_character_,
      final_n_starts = 1L,
      selected_start = 1L,
      interpolation_flag = is.finite(linear_fit$fit) &&
        linear_fit$fit > interpolation_fit_threshold,
      stringsAsFactors = FALSE
    )
  )

  kernel_test_metrics <-
    .sim_prediction_metrics(dat$y_train, dat$y_test, kernel_pred)

  kernel_perf <- cbind(
    data.frame(method = "Gaussian Kernel ERA", stringsAsFactors = FALSE),
    kernel_test_metrics,
    data.frame(
      train_fit = kernel_fit$fit,
      train_test_fit_gap = kernel_fit$fit - kernel_test_metrics$test_r2,
      converged = kernel_fit$converged,
      iterations = kernel_fit$iterations,
      lambda = cv$best_lambda,
      sigma = cv$best_sigma,
      cv_mse = selected_cv_row$cv_mse,
      cv_se = selected_cv_row$cv_se,
      cv_convergence_rate = selected_cv_row$convergence_rate,
      cv_n_successful_folds = selected_cv_row$n_successful_folds,
      cv_eligible = selected_cv_row$eligible,
      cv_selection_rule = cv_selection,
      final_n_starts = kernel_fit$n_starts,
      selected_start = kernel_fit$selected_start,
      interpolation_flag = is.finite(kernel_fit$fit) &&
        kernel_fit$fit > interpolation_fit_threshold,
      stringsAsFactors = FALSE
    )
  )

  oracle_perf <- cbind(
    data.frame(method = "Oracle signal", stringsAsFactors = FALSE),
    .sim_prediction_metrics(dat$y_train, dat$y_test, dat$signal_test),
    data.frame(
      train_fit = NA_real_,
      train_test_fit_gap = NA_real_,
      converged = TRUE,
      iterations = 0L,
      lambda = NA_real_,
      sigma = NA_real_,
      cv_mse = NA_real_,
      cv_se = NA_real_,
      cv_convergence_rate = NA_real_,
      cv_n_successful_folds = NA_integer_,
      cv_eligible = NA,
      cv_selection_rule = NA_character_,
      final_n_starts = NA_integer_,
      selected_start = NA_integer_,
      interpolation_flag = NA,
      stringsAsFactors = FALSE
    )
  )

  performance <- rbind(linear_perf, kernel_perf, oracle_perf)
  performance$condition_id <- condition_id
  performance$replication <- replication
  performance$scenario <- dat$scenario
  performance$n_train <- nrow(dat$X_train)
  performance$n_test <- nrow(dat$X_test)
  performance$p_per_set <- paste(dat$nvar, collapse = ":")
  performance$rho_within <- dat$rho_within
  performance$rho_between <- dat$rho_between
  performance$target_r2 <- dat$target_r2
  performance$penalty <- ifelse(performance$method == "Gaussian Kernel ERA", penalty, NA_character_)
  performance$gaussian_denominator <- ifelse(
    performance$method == "Gaussian Kernel ERA", gaussian_denominator, NA_character_)

  components <- rbind(
    .sim_component_metrics(linear_fit$F, linear_F_test,
                           dat$true_components_train, dat$true_components_test,
                           "Linear ERA"),
    .sim_component_metrics(kernel_fit$F, kernel_F_test,
                           dat$true_components_train, dat$true_components_test,
                           "Gaussian Kernel ERA")
  )
  components$condition_id <- condition_id
  components$replication <- replication
  components$scenario <- dat$scenario

  true_B_std <- dat$beta / kernel_fit$y_scaler$scale[1L]
  coefficients <- rbind(
    data.frame(method = "Linear ERA", predictor_set = seq_along(dat$beta),
               estimate = as.numeric(linear_fit$B[, 1L]),
               truth_standardized_scale = dat$beta / linear_fit$y_scaler$scale[1L]),
    data.frame(method = "Gaussian Kernel ERA", predictor_set = seq_along(dat$beta),
               estimate = as.numeric(kernel_fit$B[, 1L]),
               truth_standardized_scale = true_B_std)
  )
  coefficients$bias <- coefficients$estimate - coefficients$truth_standardized_scale
  coefficients$condition_id <- condition_id
  coefficients$replication <- replication
  coefficients$scenario <- dat$scenario

  tuning <- cv$table
  tuning$condition_id <- condition_id
  tuning$replication <- replication
  tuning$scenario <- dat$scenario
  tuning$selection_rule <- cv_selection
  tuning$require_all_folds_converged <-
    cv_require_all_folds_converged

  inference <- NULL
  if (bootstrap_reps > 0L) {
    boot <- kernel_era_bootstrap(
      model = kernel_fit,
      X = dat$X_train,
      Y = dat$y_train,
      n_boot = bootstrap_reps,
      seed = seed + 40L,
      retune = FALSE,
      max_iter = max_iter,
      tol = tol,
      objective_tol = objective_tol,
      verbose = FALSE
    )
    inference <- data.frame(
      predictor_set = seq_along(dat$beta),
      estimate = as.numeric(kernel_fit$B[, 1L]),
      truth_standardized_scale = true_B_std,
      bootstrap_se = as.numeric(boot$se[, 1L]),
      ci_lower = as.numeric(boot$lower[, 1L]),
      ci_upper = as.numeric(boot$upper[, 1L]),
      covered = true_B_std >= as.numeric(boot$lower[, 1L]) &
        true_B_std <= as.numeric(boot$upper[, 1L]),
      bootstrap_convergence_rate = boot$convergence_rate,
      condition_id = condition_id,
      replication = replication,
      scenario = dat$scenario,
      stringsAsFactors = FALSE
    )
  }

  matlab_gcv <- NULL
  if (isTRUE(include_matlab_gcv)) {
    matlab_gcv <- kernel_era_matlab_gcv(
      X = dat$X_train,
      Y = dat$y_train,
      nvar = dat$nvar,
      lambda_grid = lambda_grid,
      sigma_grid = sigma_grid,
      gaussian_denominator = gaussian_denominator,
      penalty = penalty,
      seed = seed + 50L,
      max_iter = max_iter,
      tol = tol,
      objective_tol = objective_tol,
      verbose = FALSE
    )
  }

  elapsed_seconds <- proc.time()[["elapsed"]] - elapsed_start
  performance$elapsed_seconds <- elapsed_seconds

  if (verbose) {
    cat("Completed condition", condition_id, "replication", replication,
        "in", round(elapsed_seconds, 2), "seconds\n")
  }

  list(
    performance = performance,
    components = components,
    coefficients = coefficients,
    tuning = tuning,
    inference = inference,
    matlab_gcv = matlab_gcv,
    selected = data.frame(
      condition_id = condition_id,
      scenario = dat$scenario,
      replication = replication,
      lambda = cv$best_lambda,
      sigma = cv$best_sigma,
      selected_cv_mse = selected_cv_row$cv_mse,
      selected_cv_se = selected_cv_row$cv_se,
      selected_cv_convergence_rate = selected_cv_row$convergence_rate,
      selected_cv_n_successful_folds =
        selected_cv_row$n_successful_folds,
      selected_cv_eligible = selected_cv_row$eligible,
      cv_selection_rule = cv_selection,
      cv_require_all_folds_converged =
        cv_require_all_folds_converged,
      n_eligible_candidates = cv$n_eligible_candidates,
      n_total_candidates = cv$n_total_candidates,
      unrestricted_min_lambda = unrestricted_cv_row$lambda,
      unrestricted_min_sigma = unrestricted_cv_row$sigma,
      unrestricted_min_cv_mse = unrestricted_cv_row$cv_mse,
      unrestricted_min_cv_se = unrestricted_cv_row$cv_se,
      unrestricted_min_convergence_rate =
        unrestricted_cv_row$convergence_rate,
      unrestricted_min_eligible = unrestricted_cv_row$eligible,
      unrestricted_min_was_selected =
        isTRUE(unrestricted_cv_row$selected),
      lambda_at_lower_boundary =
        cv$best_lambda == min(lambda_grid),
      lambda_at_upper_boundary =
        cv$best_lambda == max(lambda_grid),
      sigma_at_lower_boundary =
        cv$best_sigma == min(sigma_grid),
      sigma_at_upper_boundary =
        cv$best_sigma == max(sigma_grid),
      final_converged = kernel_fit$converged,
      final_iterations = kernel_fit$iterations,
      final_train_fit = kernel_fit$fit,
      final_test_r2 = kernel_test_metrics$test_r2,
      final_train_test_fit_gap =
        kernel_fit$fit - kernel_test_metrics$test_r2,
      final_interpolation_flag =
        is.finite(kernel_fit$fit) &&
        kernel_fit$fit > interpolation_fit_threshold,
      final_n_starts = kernel_fit$n_starts,
      final_selected_start = kernel_fit$selected_start,
      any_start_converged = kernel_fit$any_start_converged,
      elapsed_seconds = elapsed_seconds,
      stringsAsFactors = FALSE
    )
  )
}

# -----------------------------------------------------------------------------
# Recommended design grids
# -----------------------------------------------------------------------------

make_kernel_era_quick_grid <- function() {
  data.frame(
    condition_id = 1:3,
    scenario = c("linear", "mixed", "nonlinear"),
    n_train = 150L,
    n_test = 400L,
    n_sets = 3L,
    p_per_set = 5L,
    rho_within = 0.4,
    rho_between = 0.1,
    target_r2 = 0.5,
    nonlinearity_strength = 0.35,
    seed_base = 20260707L,
    stringsAsFactors = FALSE
  )
}

make_kernel_era_main_grid <- function() {
  grid <- expand.grid(
    scenario = c("linear", "mixed", "nonlinear"),
    n_train = c(150L, 300L),
    p_per_set = c(5L, 20L),
    target_r2 = c(0.30, 0.60),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$condition_id <- seq_len(nrow(grid))
  grid$n_test <- 1000L
  grid$n_sets <- 3L
  grid$rho_within <- 0.4
  grid$rho_between <- 0.1
  grid$nonlinearity_strength <- 0.35
  grid$seed_base <- 20260707L
  grid[, c("condition_id", setdiff(names(grid), "condition_id"))]
}

make_kernel_era_correlation_sensitivity_grid <- function() {
  grid <- expand.grid(
    scenario = c("linear", "mixed", "nonlinear"),
    rho_within = c(0.10, 0.70),
    rho_between = c(0.00, 0.20),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid <- grid[grid$rho_between <= grid$rho_within, , drop = FALSE]
  grid$condition_id <- seq_len(nrow(grid))
  grid$n_train <- 300L
  grid$n_test <- 1000L
  grid$n_sets <- 3L
  grid$p_per_set <- 10L
  grid$target_r2 <- 0.50
  grid$nonlinearity_strength <- 0.35
  grid$seed_base <- 20260707L + 5000000L
  grid[, c("condition_id", setdiff(names(grid), "condition_id"))]
}

# -----------------------------------------------------------------------------
# Resumable simulation driver
# -----------------------------------------------------------------------------

.combine_result_element <- function(results, element) {
  pieces <- lapply(results, `[[`, element)
  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) == 0L) return(NULL)
  do.call(rbind, pieces)
}

run_kernel_era_simulation <- function(design,
                                      n_rep = 100L,
                                      results_dir = "kernel_era_simulation_results",
                                      lambda_grid = 10^seq(-4, 2, by = 1),
                                      sigma_grid = c(0.75, 1.5, 3, 6, 12, 24, 48),
                                      cv_folds = 5L,
                                      cv_selection = c("one_se", "minimum"),
                                      cv_require_all_folds_converged = TRUE,
                                      gaussian_denominator = c("two_sigma_squared", "sigma_squared"),
                                      penalty = c("rkhs", "dual_l2"),
                                      final_n_starts = 3L,
                                      max_iter = 500L,
                                      tol = 1e-6,
                                      objective_tol = 1e-9,
                                      interpolation_fit_threshold = 0.98,
                                      bootstrap_reps = 0L,
                                      include_matlab_gcv = FALSE,
                                      n_cores = 1L,
                                      overwrite = FALSE,
                                      verbose = TRUE) {
  cv_selection <- match.arg(cv_selection)
  gaussian_denominator <- match.arg(gaussian_denominator)
  penalty <- match.arg(penalty)
  cv_require_all_folds_converged <-
    isTRUE(cv_require_all_folds_converged)
  design <- as.data.frame(design, stringsAsFactors = FALSE)
  if (!"condition_id" %in% names(design)) design$condition_id <- seq_len(nrow(design))
  n_rep <- as.integer(n_rep)
  n_cores <- as.integer(n_cores)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  rep_dir <- file.path(results_dir, "replications")
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

  tasks <- do.call(rbind, lapply(seq_len(nrow(design)), function(i) {
    data.frame(design_row = i,
               condition_id = design$condition_id[i],
               replication = seq_len(n_rep))
  }))

  run_task <- function(task_row) {
    i <- task_row$design_row
    r <- task_row$replication
    cid <- task_row$condition_id
    file <- file.path(rep_dir, sprintf("condition_%03d_rep_%05d.rds", cid, r))
    if (file.exists(file) && !overwrite) return(readRDS(file))

    result <- try(run_one_kernel_era_replication(
      condition = design[i, , drop = FALSE],
      replication = r,
      lambda_grid = lambda_grid,
      sigma_grid = sigma_grid,
      cv_folds = cv_folds,
      cv_selection = cv_selection,
      cv_require_all_folds_converged =
        cv_require_all_folds_converged,
      gaussian_denominator = gaussian_denominator,
      penalty = penalty,
      final_n_starts = final_n_starts,
      max_iter = max_iter,
      tol = tol,
      objective_tol = objective_tol,
      interpolation_fit_threshold =
        interpolation_fit_threshold,
      bootstrap_reps = bootstrap_reps,
      include_matlab_gcv = include_matlab_gcv,
      verbose = FALSE
    ), silent = TRUE)

    if (inherits(result, "try-error")) {
      result <- list(
        error = data.frame(condition_id = cid, replication = r,
                           message = as.character(result), stringsAsFactors = FALSE)
      )
    }
    saveRDS(result, file)
    result
  }

  task_list <- split(tasks, seq_len(nrow(tasks)))
  if (n_cores > 1L && requireNamespace("future.apply", quietly = TRUE) &&
      requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = n_cores)
    results <- future.apply::future_lapply(
      task_list,
      function(z) run_task(z[1L, , drop = FALSE]),
      future.seed = TRUE,
      future.packages = character(0)
    )
  } else {
    if (n_cores > 1L && verbose) {
      message("future and future.apply are unavailable; running sequentially.")
    }
    results <- lapply(seq_along(task_list), function(j) {
      if (verbose && (j %% max(1L, floor(length(task_list) / 100L)) == 0L)) {
        message("Completed/loaded ", j, " of ", length(task_list), " tasks")
      }
      run_task(task_list[[j]][1L, , drop = FALSE])
    })
  }

  elements <- c("performance", "components", "coefficients", "tuning",
                "inference", "selected", "error")
  combined <- setNames(lapply(elements, function(x) .combine_result_element(results, x)),
                       elements)

  utils::write.csv(design, file.path(results_dir, "simulation_design.csv"), row.names = FALSE)
  for (nm in names(combined)) {
    if (!is.null(combined[[nm]])) {
      utils::write.csv(combined[[nm]], file.path(results_dir, paste0(nm, ".csv")),
                       row.names = FALSE)
    }
  }
  saveRDS(combined, file.path(results_dir, "combined_results.rds"))

  invisible(combined)
}

# -----------------------------------------------------------------------------
# Summary helpers
# -----------------------------------------------------------------------------

summarize_kernel_era_simulation <- function(results) {
  if (is.character(results) && length(results) == 1L) {
    results <- readRDS(results)
  }

  perf <- results$performance
  comp <- results$components
  selected <- results$selected

  mean_or_na <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    mean(x)
  }

  sd_or_na <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) return(NA_real_)
    stats::sd(x)
  }

  perf_summary <- NULL
  if (!is.null(perf)) {
    group <- interaction(
      perf$condition_id,
      perf$scenario,
      perf$method,
      drop = TRUE
    )

    perf_summary <- do.call(
      rbind,
      lapply(split(perf, group), function(d) {
        data.frame(
          condition_id = d$condition_id[1L],
          scenario = d$scenario[1L],
          method = d$method[1L],
          n_rep = nrow(d),
          mean_test_mse = mean_or_na(d$test_mse),
          sd_test_mse = sd_or_na(d$test_mse),
          mean_test_r2 = mean_or_na(d$test_r2),
          sd_test_r2 = sd_or_na(d$test_r2),
          mean_prediction_correlation =
            mean_or_na(d$prediction_correlation),
          convergence_rate = mean(d$converged, na.rm = TRUE),
          mean_train_fit = mean_or_na(d$train_fit),
          mean_train_test_fit_gap =
            mean_or_na(d$train_test_fit_gap),
          interpolation_rate = if (
            all(is.na(d$interpolation_flag))
          ) {
            NA_real_
          } else {
            mean(d$interpolation_flag, na.rm = TRUE)
          },
          stringsAsFactors = FALSE
        )
      })
    )
  }

  component_summary <- NULL
  if (!is.null(comp)) {
    group <- interaction(
      comp$condition_id,
      comp$scenario,
      comp$method,
      comp$predictor_set,
      drop = TRUE
    )

    component_summary <- do.call(
      rbind,
      lapply(split(comp, group), function(d) {
        data.frame(
          condition_id = d$condition_id[1L],
          scenario = d$scenario[1L],
          method = d$method[1L],
          predictor_set = d$predictor_set[1L],
          mean_test_component_correlation =
            mean_or_na(d$test_component_correlation),
          sd_test_component_correlation =
            sd_or_na(d$test_component_correlation),
          mean_test_component_rmse =
            mean_or_na(d$test_component_rmse_calibrated),
          stringsAsFactors = FALSE
        )
      })
    )
  }

  tuning_summary <- NULL
  if (!is.null(selected)) {
    group <- interaction(
      selected$condition_id,
      selected$scenario,
      drop = TRUE
    )

    tuning_summary <- do.call(
      rbind,
      lapply(split(selected, group), function(d) {
        data.frame(
          condition_id = d$condition_id[1L],
          scenario = d$scenario[1L],
          n_rep = nrow(d),
          selected_eligibility_rate =
            mean(d$selected_cv_eligible, na.rm = TRUE),
          selected_cv_full_convergence_rate =
            mean(d$selected_cv_convergence_rate == 1,
                 na.rm = TRUE),
          unrestricted_minimum_exclusion_rate =
            mean(!d$unrestricted_min_eligible, na.rm = TRUE),
          final_convergence_rate =
            mean(d$final_converged, na.rm = TRUE),
          final_interpolation_rate =
            mean(d$final_interpolation_flag, na.rm = TRUE),
          lambda_lower_boundary_rate =
            mean(d$lambda_at_lower_boundary),
          lambda_upper_boundary_rate =
            mean(d$lambda_at_upper_boundary),
          sigma_lower_boundary_rate =
            mean(d$sigma_at_lower_boundary),
          sigma_upper_boundary_rate =
            mean(d$sigma_at_upper_boundary),
          median_lambda = stats::median(d$lambda),
          median_sigma = stats::median(d$sigma),
          median_selected_cv_mse =
            stats::median(d$selected_cv_mse),
          stringsAsFactors = FALSE
        )
      })
    )
  }

  list(
    performance = perf_summary,
    components = component_summary,
    tuning = tuning_summary
  )
}