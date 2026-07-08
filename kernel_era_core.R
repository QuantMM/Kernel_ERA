# =============================================================================
# Kernel Extended Redundancy Analysis (Kernel ERA)
# Core estimation, tuning, prediction, and bootstrap utilities
#
# Default specification follows the current Kernel ERA manuscript:
#   * centered predictor-set Gram matrices
#   * Gaussian kernel exp(-||x-z||^2 / (2 sigma^2))
#   * shared sigma across predictor sets
#   * RKHS-norm penalty lambda * sum_k alpha_k' K_k alpha_k
#   * component normalization f_k'f_k / N = 1
#   * Gauss-Seidel alternating least squares updates
#
# The code uses only base R. Optional parallel execution is handled in the
# separate simulation driver.
# =============================================================================

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# -----------------------------------------------------------------------------
# Basic validation and numerical helpers
# -----------------------------------------------------------------------------

.kera_as_matrix <- function(x, name = deparse(substitute(x))) {
  if (is.data.frame(x)) x <- data.matrix(x)
  if (is.vector(x)) x <- matrix(x, ncol = 1L)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop(name, " must be a finite numeric vector or matrix.", call. = FALSE)
  }
  x
}

.kera_validate_nvar <- function(nvar, p) {
  nvar <- as.integer(nvar)
  if (length(nvar) < 1L || anyNA(nvar) || any(nvar <= 0L)) {
    stop("nvar must contain positive integers.", call. = FALSE)
  }
  if (sum(nvar) != p) {
    stop("sum(nvar) must equal ncol(X).", call. = FALSE)
  }
  nvar
}

.kera_block_indices <- function(nvar) {
  ends <- cumsum(nvar)
  starts <- c(1L, head(ends, -1L) + 1L)
  Map(seq.int, starts, ends)
}

.kera_split_blocks <- function(X, nvar) {
  idx <- .kera_block_indices(nvar)
  lapply(idx, function(ii) X[, ii, drop = FALSE])
}

.kera_safe_cor <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y) || length(x) < 2L) return(NA_real_)
  sx <- sqrt(mean((x - mean(x))^2))
  sy <- sqrt(mean((y - mean(y))^2))
  if (!is.finite(sx) || !is.finite(sy) || sx <= 0 || sy <= 0) return(NA_real_)
  mean((x - mean(x)) * (y - mean(y))) / (sx * sy)
}

.kera_pinv_solve <- function(A, B, tol = sqrt(.Machine$double.eps)) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  sv <- svd(A)
  if (length(sv$d) == 0L || max(sv$d) <= 0) {
    return(matrix(0, nrow = ncol(A), ncol = ncol(B)))
  }
  keep <- sv$d > tol * max(sv$d)
  if (!any(keep)) {
    return(matrix(0, nrow = ncol(A), ncol = ncol(B)))
  }
  tmp <- crossprod(sv$u[, keep, drop = FALSE], B)
  tmp <- sweep(tmp, 1L, sv$d[keep], "/")
  sv$v[, keep, drop = FALSE] %*% tmp
}

.kera_solve_psd <- function(A, B, tol = sqrt(.Machine$double.eps)) {
  A <- (as.matrix(A) + t(as.matrix(A))) / 2
  B <- as.matrix(B)

  ch <- try(chol(A), silent = TRUE)
  if (!inherits(ch, "try-error")) {
    z <- forwardsolve(t(ch), B)
    return(backsolve(ch, z))
  }
  .kera_pinv_solve(A, B, tol = tol)
}

.kera_lm_solve <- function(F, Y, tol = sqrt(.Machine$double.eps)) {
  .kera_pinv_solve(crossprod(F), crossprod(F, Y), tol = tol)
}

.kera_frobenius_sq <- function(A) sum(A * A)

.kera_normalize_component <- function(coef, score, n, tol = 1e-12) {
  score_norm <- sqrt(sum(score * score))
  if (!is.finite(score_norm) || score_norm <= tol) return(NULL)
  multiplier <- sqrt(n) / score_norm
  list(coef = as.numeric(coef) * multiplier,
       score = as.numeric(score) * multiplier)
}

# -----------------------------------------------------------------------------
# Divisor-N standardization
# Equivalent to MATLAB: sqrt(N) * zscore(x) / sqrt(N - 1)
# -----------------------------------------------------------------------------

kera_fit_scaler <- function(x, tol = 1e-12) {
  x <- .kera_as_matrix(x)
  center <- colMeans(x)
  xc <- sweep(x, 2L, center, "-")
  scale <- sqrt(colMeans(xc * xc))
  bad <- !is.finite(scale) | scale <= tol
  if (any(bad)) {
    stop("Zero or near-zero variance column(s): ",
         paste(which(bad), collapse = ", "), call. = FALSE)
  }
  structure(list(center = center, scale = scale), class = "kera_scaler")
}

kera_apply_scaler <- function(x, scaler) {
  x <- .kera_as_matrix(x)
  if (ncol(x) != length(scaler$center)) {
    stop("The number of columns does not match the fitted scaler.", call. = FALSE)
  }
  sweep(sweep(x, 2L, scaler$center, "-"), 2L, scaler$scale, "/")
}

kera_inverse_scaler <- function(z, scaler) {
  z <- .kera_as_matrix(z)
  if (ncol(z) != length(scaler$center)) {
    stop("The number of columns does not match the fitted scaler.", call. = FALSE)
  }
  sweep(sweep(z, 2L, scaler$scale, "*"), 2L, scaler$center, "+")
}

# -----------------------------------------------------------------------------
# Kernels and feature-space centering
# -----------------------------------------------------------------------------

kera_sqdist <- function(A, B = A) {
  A <- .kera_as_matrix(A)
  B <- .kera_as_matrix(B)
  if (ncol(A) != ncol(B)) stop("A and B must have the same columns.", call. = FALSE)
  d2 <- outer(rowSums(A * A), rowSums(B * B), "+") - 2 * tcrossprod(A, B)
  pmax(d2, 0)
}

kera_raw_kernel <- function(A,
                            B = A,
                            kernel = c("gaussian", "polynomial", "linear"),
                            sigma = 1,
                            degree = 2,
                            offset = 1,
                            gaussian_denominator = c("two_sigma_squared", "sigma_squared")) {
  kernel <- match.arg(kernel)
  gaussian_denominator <- match.arg(gaussian_denominator)
  A <- .kera_as_matrix(A)
  B <- .kera_as_matrix(B)

  if (kernel == "gaussian") {
    if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
      stop("sigma must be one positive finite number.", call. = FALSE)
    }
    denom <- if (gaussian_denominator == "two_sigma_squared") {
      2 * sigma^2
    } else {
      sigma^2
    }
    return(exp(-kera_sqdist(A, B) / denom))
  }

  if (kernel == "polynomial") {
    if (length(degree) != 1L || degree < 1 || degree != as.integer(degree)) {
      stop("degree must be a positive integer.", call. = FALSE)
    }
    return((tcrossprod(A, B) + offset)^as.integer(degree))
  }

  tcrossprod(A, B)
}

kera_center_train_kernel <- function(K_raw) {
  K_raw <- .kera_as_matrix(K_raw)
  if (nrow(K_raw) != ncol(K_raw)) stop("Training kernel must be square.", call. = FALSE)
  row_mean <- rowMeans(K_raw)
  col_mean <- colMeans(K_raw)
  grand_mean <- mean(K_raw)
  Kc <- sweep(K_raw, 1L, row_mean, "-")
  Kc <- sweep(Kc, 2L, col_mean, "-")
  Kc <- Kc + grand_mean
  Kc <- (Kc + t(Kc)) / 2
  list(K = Kc,
       train_row_mean = row_mean,
       train_col_mean = col_mean,
       train_grand_mean = grand_mean)
}

kera_center_cross_kernel <- function(K_cross, centering) {
  K_cross <- .kera_as_matrix(K_cross)
  if (ncol(K_cross) != length(centering$train_col_mean)) {
    stop("Cross-kernel columns must equal the training sample size.", call. = FALSE)
  }
  test_row_mean <- rowMeans(K_cross)
  Kc <- sweep(K_cross, 1L, test_row_mean, "-")
  Kc <- sweep(Kc, 2L, centering$train_col_mean, "-")
  Kc + centering$train_grand_mean
}

.kera_kernel_spectral <- function(K, tol = 1e-10) {
  ee <- eigen((K + t(K)) / 2, symmetric = TRUE)
  cutoff <- max(1, max(abs(ee$values))) * tol
  keep <- ee$values > cutoff
  list(values = ee$values[keep],
       vectors = ee$vectors[, keep, drop = FALSE],
       rank = sum(keep),
       tol = tol)
}

.kera_kernel_pc_fallback <- function(spectral, n) {
  if (spectral$rank < 1L) return(NULL)
  d1 <- spectral$values[1L]
  u1 <- spectral$vectors[, 1L]
  score <- sqrt(n) * u1 / sqrt(sum(u1^2))
  alpha <- score / d1
  list(coef = as.numeric(alpha), score = as.numeric(score))
}

.kera_update_alpha <- function(K,
                               spectral,
                               rhs,
                               q,
                               lambda,
                               penalty = c("rkhs", "dual_l2"),
                               solver = c("spectral", "direct"),
                               tol = 1e-10) {
  penalty <- match.arg(penalty)
  solver <- match.arg(solver)
  rhs <- as.numeric(rhs)
  n <- length(rhs)

  if (solver == "spectral" && spectral$rank > 0L) {
    U <- spectral$vectors
    d <- spectral$values
    urhs <- as.numeric(crossprod(U, rhs))

    if (penalty == "rkhs") {
      denom <- q * d + lambda
      good <- abs(denom) > tol
      coef_eig <- numeric(length(d))
      coef_eig[good] <- urhs[good] / denom[good]
    } else {
      # (q K^2 + lambda I) alpha = K rhs
      denom <- q * d^2 + lambda
      good <- abs(denom) > tol
      coef_eig <- numeric(length(d))
      coef_eig[good] <- d[good] * urhs[good] / denom[good]
    }
    return(as.numeric(U %*% coef_eig))
  }

  I_n <- diag(n)
  if (penalty == "rkhs") {
    .kera_solve_psd(q * K + lambda * I_n, rhs, tol = tol)[, 1L]
  } else {
    .kera_solve_psd(q * (K %*% K) + lambda * I_n,
                    K %*% rhs, tol = tol)[, 1L]
  }
}

# -----------------------------------------------------------------------------
# Linear ERA warm start and comparator
# -----------------------------------------------------------------------------

.kera_initialize_linear_weights <- function(X_blocks,
                                            init = c("ones", "random", "pca"),
                                            seed = NULL,
                                            tol = 1e-12) {
  init <- match.arg(init)
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X_blocks[[1L]])

  lapply(X_blocks, function(Xk) {
    p <- ncol(Xk)
    w <- switch(init,
      ones = rep(1, p),
      random = stats::rnorm(p),
      pca = {
        vv <- svd(Xk, nu = 0L, nv = 1L)$v[, 1L]
        as.numeric(vv)
      }
    )
    score <- as.numeric(Xk %*% w)
    normed <- .kera_normalize_component(w, score, n, tol = tol)
    if (is.null(normed)) {
      vv <- svd(Xk, nu = 0L, nv = 1L)$v[, 1L]
      score <- as.numeric(Xk %*% vv)
      normed <- .kera_normalize_component(vv, score, n, tol = tol)
    }
    if (is.null(normed)) stop("Unable to initialize a linear ERA component.", call. = FALSE)
    normed
  })
}

kera_linear_era_core <- function(X_blocks,
                                 Y,
                                 init = c("ones", "random", "pca"),
                                 seed = 1,
                                 max_iter = 200L,
                                 tol = 1e-8,
                                 numerical_tol = 1e-10,
                                 verbose = FALSE) {
  init <- match.arg(init)
  Y <- .kera_as_matrix(Y)
  n <- nrow(Y)
  K <- length(X_blocks)
  if (K < 1L) stop("At least one predictor block is required.", call. = FALSE)
  if (any(vapply(X_blocks, nrow, integer(1L)) != n)) {
    stop("All predictor blocks and Y must have the same rows.", call. = FALSE)
  }

  initialized <- .kera_initialize_linear_weights(X_blocks, init, seed)
  W <- lapply(initialized, `[[`, "coef")
  F <- do.call(cbind, lapply(initialized, `[[`, "score"))
  B <- .kera_lm_solve(F, Y, tol = numerical_tol)

  objective_history <- numeric(max_iter + 1L)
  objective_history[1L] <- .kera_frobenius_sq(Y - F %*% B)
  converged <- FALSE
  iterations <- 0L

  for (it in seq_len(max_iter)) {
    B_old <- B

    for (k in seq_len(K)) {
      fitted_current <- F %*% B
      Rk <- Y - fitted_current + tcrossprod(F[, k], B[k, ])
      bk <- as.numeric(B[k, ])
      qk <- sum(bk * bk)
      if (!is.finite(qk) || qk <= numerical_tol) next

      rhs <- crossprod(X_blocks[[k]], Rk %*% bk)
      lhs <- qk * crossprod(X_blocks[[k]])
      w_tilde <- .kera_solve_psd(lhs, rhs, tol = numerical_tol)[, 1L]
      score_tilde <- as.numeric(X_blocks[[k]] %*% w_tilde)
      normed <- .kera_normalize_component(w_tilde, score_tilde, n,
                                          tol = numerical_tol)
      if (!is.null(normed)) {
        W[[k]] <- normed$coef
        F[, k] <- normed$score
      }
    }

    B <- .kera_lm_solve(F, Y, tol = numerical_tol)
    objective_history[it + 1L] <- .kera_frobenius_sq(Y - F %*% B)
    delta <- sum(abs(B - B_old))
    iterations <- it

    if (verbose) {
      cat(sprintf("Linear ERA iter %d: RSS=%.8g, delta_B=%.4g\n",
                  it, objective_history[it + 1L], delta))
    }
    if (is.finite(delta) && delta <= tol) {
      converged <- TRUE
      break
    }
  }

  objective_history <- objective_history[seq_len(iterations + 1L)]
  list(W = W,
       F = F,
       B = B,
       fitted = F %*% B,
       residuals = Y - F %*% B,
       rss = .kera_frobenius_sq(Y - F %*% B),
       objective_history = objective_history,
       converged = converged,
       iterations = iterations)
}

linear_era_fit <- function(X,
                           Y,
                           nvar,
                           standardize = TRUE,
                           init = c("ones", "random", "pca"),
                           seed = 1,
                           max_iter = 200L,
                           tol = 1e-8,
                           numerical_tol = 1e-10,
                           verbose = FALSE) {
  init <- match.arg(init)
  X <- .kera_as_matrix(X)
  Y <- .kera_as_matrix(Y)
  if (nrow(X) != nrow(Y)) stop("X and Y must have the same rows.", call. = FALSE)
  nvar <- .kera_validate_nvar(nvar, ncol(X))
  block_idx <- .kera_block_indices(nvar)

  if (standardize) {
    x_scalers <- lapply(block_idx, function(ii) kera_fit_scaler(X[, ii, drop = FALSE]))
    X_blocks <- Map(function(ii, sc) kera_apply_scaler(X[, ii, drop = FALSE], sc),
                    block_idx, x_scalers)
    y_scaler <- kera_fit_scaler(Y)
    Y_std <- kera_apply_scaler(Y, y_scaler)
  } else {
    x_scalers <- lapply(block_idx, function(ii) {
      structure(list(center = rep(0, length(ii)), scale = rep(1, length(ii))),
                class = "kera_scaler")
    })
    X_blocks <- .kera_split_blocks(X, nvar)
    y_scaler <- structure(list(center = rep(0, ncol(Y)), scale = rep(1, ncol(Y))),
                          class = "kera_scaler")
    Y_std <- Y
  }

  core <- kera_linear_era_core(X_blocks, Y_std, init = init, seed = seed,
                               max_iter = max_iter, tol = tol,
                               numerical_tol = numerical_tol,
                               verbose = verbose)
  ans <- c(core, list(
    call = match.call(),
    n = nrow(X),
    q = ncol(Y),
    nvar = nvar,
    block_indices = block_idx,
    x_scalers = x_scalers,
    y_scaler = y_scaler,
    X_train_std = X_blocks,
    Y_train_std = Y_std,
    standardize = standardize,
    component_correlation = crossprod(core$F) / nrow(X),
    fit = 1 - core$rss / .kera_frobenius_sq(Y_std)
  ))
  class(ans) <- "linear_era"
  ans
}

predict.linear_era <- function(object,
                               newX = NULL,
                               type = c("response", "standardized", "components", "all"),
                               ...) {
  type <- match.arg(type)
  if (is.null(newX)) {
    F_new <- object$F
  } else {
    newX <- .kera_as_matrix(newX)
    if (ncol(newX) != sum(object$nvar)) stop("newX has incorrect columns.", call. = FALSE)
    X_blocks <- Map(function(ii, sc) {
      kera_apply_scaler(newX[, ii, drop = FALSE], sc)
    }, object$block_indices, object$x_scalers)
    F_new <- do.call(cbind, Map(function(Xk, wk) as.numeric(Xk %*% wk),
                                X_blocks, object$W))
  }
  Y_std <- F_new %*% object$B
  Y_response <- kera_inverse_scaler(Y_std, object$y_scaler)

  if (type == "components") return(F_new)
  if (type == "standardized") return(Y_std)
  if (type == "response") return(if (ncol(Y_response) == 1L) as.numeric(Y_response) else Y_response)
  list(response = if (ncol(Y_response) == 1L) as.numeric(Y_response) else Y_response,
       standardized = Y_std,
       components = F_new)
}

print.linear_era <- function(x, ...) {
  cat("Linear Extended Redundancy Analysis\n")
  cat("  N:", x$n, "  predictor sets:", length(x$nvar), "  outcomes:", x$q, "\n")
  cat("  converged:", x$converged, "  iterations:", x$iterations, "\n")
  cat("  FIT:", format(x$fit, digits = 5), "\n")
  invisible(x)
}

# -----------------------------------------------------------------------------
# Kernel ERA data preparation
# -----------------------------------------------------------------------------

kera_prepare_data <- function(X,
                              Y,
                              nvar,
                              kernel = c("gaussian", "polynomial", "linear"),
                              sigma = 1,
                              degree = 2,
                              offset = 1,
                              gaussian_denominator = c("two_sigma_squared", "sigma_squared"),
                              standardize = TRUE,
                              spectral_tol = 1e-10) {
  kernel <- match.arg(kernel)
  gaussian_denominator <- match.arg(gaussian_denominator)
  X <- .kera_as_matrix(X)
  Y <- .kera_as_matrix(Y)
  if (nrow(X) != nrow(Y)) stop("X and Y must have the same rows.", call. = FALSE)
  nvar <- .kera_validate_nvar(nvar, ncol(X))
  Ksets <- length(nvar)
  block_idx <- .kera_block_indices(nvar)

  if (length(sigma) == 1L) sigma <- rep(sigma, Ksets)
  if (length(sigma) != Ksets) {
    stop("sigma must be scalar (shared) or have one value per predictor set.", call. = FALSE)
  }

  if (standardize) {
    x_scalers <- lapply(block_idx, function(ii) kera_fit_scaler(X[, ii, drop = FALSE]))
    X_blocks <- Map(function(ii, sc) kera_apply_scaler(X[, ii, drop = FALSE], sc),
                    block_idx, x_scalers)
    y_scaler <- kera_fit_scaler(Y)
    Y_std <- kera_apply_scaler(Y, y_scaler)
  } else {
    x_scalers <- lapply(block_idx, function(ii) {
      structure(list(center = rep(0, length(ii)), scale = rep(1, length(ii))),
                class = "kera_scaler")
    })
    X_blocks <- .kera_split_blocks(X, nvar)
    y_scaler <- structure(list(center = rep(0, ncol(Y)), scale = rep(1, ncol(Y))),
                          class = "kera_scaler")
    Y_std <- Y
  }

  kernel_objects <- vector("list", Ksets)
  for (k in seq_len(Ksets)) {
    raw <- kera_raw_kernel(X_blocks[[k]], X_blocks[[k]],
                           kernel = kernel, sigma = sigma[k],
                           degree = degree, offset = offset,
                           gaussian_denominator = gaussian_denominator)
    centered <- kera_center_train_kernel(raw)
    kernel_objects[[k]] <- c(centered,
      list(spectral = .kera_kernel_spectral(centered$K, tol = spectral_tol)))
  }

  list(
    X_raw = X,
    Y_raw = Y,
    X_blocks = X_blocks,
    Y_std = Y_std,
    n = nrow(X),
    q = ncol(Y),
    nvar = nvar,
    block_indices = block_idx,
    x_scalers = x_scalers,
    y_scaler = y_scaler,
    standardize = standardize,
    kernel = kernel,
    sigma = sigma,
    degree = degree,
    offset = offset,
    gaussian_denominator = gaussian_denominator,
    kernel_objects = kernel_objects,
    spectral_tol = spectral_tol
  )
}

# -----------------------------------------------------------------------------
# Kernel ERA estimation
# -----------------------------------------------------------------------------

.kera_penalty_value <- function(alpha, kernel_objects,
                                penalty = c("rkhs", "dual_l2")) {
  penalty <- match.arg(penalty)
  vals <- vapply(seq_along(alpha), function(k) {
    ak <- as.numeric(alpha[[k]])
    if (penalty == "rkhs") {
      as.numeric(crossprod(ak, kernel_objects[[k]]$K %*% ak))
    } else {
      sum(ak * ak)
    }
  }, numeric(1L))
  sum(vals)
}

.kera_fit_one_start <- function(prepared,
                                lambda,
                                penalty = c("rkhs", "dual_l2"),
                                solver = c("spectral", "direct"),
                                linear_init = c("ones", "random", "pca"),
                                seed = 1,
                                max_iter = 200L,
                                tol = 1e-7,
                                objective_tol = 1e-9,
                                numerical_tol = 1e-10,
                                warm_start = NULL,
                                verbose = FALSE) {
  penalty <- match.arg(penalty)
  solver <- match.arg(solver)
  linear_init <- match.arg(linear_init)
  n <- prepared$n
  Ksets <- length(prepared$nvar)
  Y <- prepared$Y_std

  warm <- warm_start
  if (is.null(warm)) {
    warm <- kera_linear_era_core(prepared$X_blocks, Y,
                                 init = linear_init,
                                 seed = seed,
                                 max_iter = max_iter,
                                 tol = tol,
                                 numerical_tol = numerical_tol,
                                 verbose = FALSE)
  }
  F <- warm$F
  B <- warm$B
  alpha <- vector("list", Ksets)

  # Kernel-compatible starting coefficients. These are used only as a robust
  # fallback when a block coefficient is numerically zero at an update.
  for (k in seq_len(Ksets)) {
    fallback <- .kera_kernel_pc_fallback(prepared$kernel_objects[[k]]$spectral, n)
    if (is.null(fallback)) {
      alpha[[k]] <- rep(0, n)
    } else {
      alpha[[k]] <- fallback$coef
    }
  }

  objective_history <- rep(NA_real_, max_iter)
  delta_B_history <- rep(NA_real_, max_iter)
  converged <- FALSE
  iterations <- 0L

  for (it in seq_len(max_iter)) {
    B_old <- B

    # Step 1: Gauss-Seidel block updates for alpha_k, holding B fixed.
    for (k in seq_len(Ksets)) {
      Kk <- prepared$kernel_objects[[k]]$K
      spectral <- prepared$kernel_objects[[k]]$spectral
      fitted_current <- F %*% B
      Rk <- Y - fitted_current + tcrossprod(F[, k], B[k, ])
      bk <- as.numeric(B[k, ])
      qk <- sum(bk * bk)
      rhs <- as.numeric(Rk %*% bk)

      alpha_tilde <- .kera_update_alpha(
        K = Kk,
        spectral = spectral,
        rhs = rhs,
        q = qk,
        lambda = lambda,
        penalty = penalty,
        solver = solver,
        tol = numerical_tol
      )
      score_tilde <- as.numeric(Kk %*% alpha_tilde)
      normed <- .kera_normalize_component(alpha_tilde, score_tilde, n,
                                          tol = numerical_tol)

      if (is.null(normed)) {
        fallback <- .kera_kernel_pc_fallback(spectral, n)
        if (!is.null(fallback)) normed <- fallback
      }
      if (!is.null(normed)) {
        alpha[[k]] <- normed$coef
        F[, k] <- normed$score
      }
    }

    # Step 2: least-squares update for B.
    B <- .kera_lm_solve(F, Y, tol = numerical_tol)
    residual <- Y - F %*% B
    penalty_value <- .kera_penalty_value(alpha, prepared$kernel_objects, penalty)
    objective <- .kera_frobenius_sq(residual) + lambda * penalty_value
    delta_B <- sum(abs(B - B_old))

    objective_history[it] <- objective
    delta_B_history[it] <- delta_B
    iterations <- it

    if (verbose) {
      cat(sprintf("Kernel ERA iter %d: objective=%.8g, delta_B=%.4g\n",
                  it, objective, delta_B))
    }

    relative_obj_change <- if (it == 1L || !is.finite(objective_history[it - 1L])) {
      Inf
    } else {
      abs(objective_history[it - 1L] - objective) /
        max(1, abs(objective_history[it - 1L]))
    }

    if (is.finite(delta_B) &&
        (delta_B <= tol ||
         (delta_B <= 10 * tol && relative_obj_change <= objective_tol))) {
      converged <- TRUE
      break
    }
  }

  objective_history <- objective_history[seq_len(iterations)]
  delta_B_history <- delta_B_history[seq_len(iterations)]
  residual <- Y - F %*% B
  rss <- .kera_frobenius_sq(residual)
  penalty_value <- .kera_penalty_value(alpha, prepared$kernel_objects, penalty)
  objective <- rss + lambda * penalty_value

  list(alpha = alpha,
       F = F,
       B = B,
       fitted = F %*% B,
       residuals = residual,
       rss = rss,
       penalty_value = penalty_value,
       objective = objective,
       objective_history = objective_history,
       delta_B_history = delta_B_history,
       converged = converged,
       iterations = iterations,
       warm_start = warm)
}

kera_fit_prepared <- function(prepared,
                              lambda,
                              penalty = c("rkhs", "dual_l2"),
                              solver = c("spectral", "direct"),
                              linear_init = c("ones", "random", "pca"),
                              n_starts = 1L,
                              seed = 1,
                              max_iter = 200L,
                              tol = 1e-7,
                              objective_tol = 1e-9,
                              numerical_tol = 1e-10,
                              warm_start = NULL,
                              verbose = FALSE) {
  penalty <- match.arg(penalty)
  solver <- match.arg(solver)
  linear_init <- match.arg(linear_init)
  if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    stop("lambda must be one nonnegative finite number.", call. = FALSE)
  }
  n_starts <- as.integer(n_starts)
  if (n_starts < 1L) stop("n_starts must be positive.", call. = FALSE)

  # The first start uses the requested linear ERA initialization. Additional
  # starts use random linear ERA initializations so that n_starts > 1 actually
  # explores distinct solutions even when linear_init = "ones".
  start_initializations <- rep("random", n_starts)
  start_initializations[1L] <- linear_init

  fits <- lapply(seq_len(n_starts), function(s) {
    .kera_fit_one_start(
      prepared = prepared,
      lambda = lambda,
      penalty = penalty,
      solver = solver,
      linear_init = start_initializations[s],
      seed = seed + s - 1L,
      max_iter = max_iter,
      tol = tol,
      objective_tol = objective_tol,
      numerical_tol = numerical_tol,
      warm_start = if (s == 1L) warm_start else NULL,
      verbose = verbose
    )
  })

  objectives <- vapply(fits, `[[`, numeric(1L), "objective")
  start_converged <- vapply(fits, `[[`, logical(1L), "converged")
  start_iterations <- vapply(fits, `[[`, integer(1L), "iterations")

  # Prefer converged solutions. Only when every start fails to converge do we
  # retain the finite solution with the smallest objective for diagnostics.
  converged_candidates <- which(start_converged & is.finite(objectives))
  if (length(converged_candidates) > 0L) {
    best <- converged_candidates[
      which.min(objectives[converged_candidates])
    ]
  } else {
    finite_candidates <- which(is.finite(objectives))
    if (length(finite_candidates) == 0L) {
      stop("All Kernel ERA starts produced non-finite objectives.", call. = FALSE)
    }
    best <- finite_candidates[
      which.min(objectives[finite_candidates])
    ]
  }
  core <- fits[[best]]

  ans <- c(core, list(
    call = match.call(),
    n = prepared$n,
    q = prepared$q,
    nvar = prepared$nvar,
    block_indices = prepared$block_indices,
    x_scalers = prepared$x_scalers,
    y_scaler = prepared$y_scaler,
    X_train_std = prepared$X_blocks,
    Y_train_std = prepared$Y_std,
    standardize = prepared$standardize,
    kernel = prepared$kernel,
    sigma = prepared$sigma,
    degree = prepared$degree,
    offset = prepared$offset,
    gaussian_denominator = prepared$gaussian_denominator,
    kernel_objects = prepared$kernel_objects,
    lambda = lambda,
    penalty = penalty,
    solver = solver,
    n_starts = n_starts,
    selected_start = best,
    start_initializations = start_initializations,
    all_start_objectives = objectives,
    all_start_converged = start_converged,
    all_start_iterations = start_iterations,
    any_start_converged = any(start_converged),
    component_correlation = crossprod(core$F) / prepared$n,
    fit = 1 - core$rss / .kera_frobenius_sq(prepared$Y_std)
  ))
  class(ans) <- "kernel_era"
  ans
}

kernel_era_fit <- function(X,
                           Y,
                           nvar,
                           lambda,
                           kernel = c("gaussian", "polynomial", "linear"),
                           sigma = 1,
                           degree = 2,
                           offset = 1,
                           gaussian_denominator = c("two_sigma_squared", "sigma_squared"),
                           penalty = c("rkhs", "dual_l2"),
                           solver = c("spectral", "direct"),
                           standardize = TRUE,
                           linear_init = c("ones", "random", "pca"),
                           n_starts = 1L,
                           seed = 1,
                           max_iter = 200L,
                           tol = 1e-7,
                           objective_tol = 1e-9,
                           numerical_tol = 1e-10,
                           spectral_tol = 1e-10,
                           verbose = FALSE) {
  kernel <- match.arg(kernel)
  gaussian_denominator <- match.arg(gaussian_denominator)
  penalty <- match.arg(penalty)
  solver <- match.arg(solver)
  linear_init <- match.arg(linear_init)

  prepared <- kera_prepare_data(
    X = X, Y = Y, nvar = nvar,
    kernel = kernel, sigma = sigma,
    degree = degree, offset = offset,
    gaussian_denominator = gaussian_denominator,
    standardize = standardize,
    spectral_tol = spectral_tol
  )
  fit <- kera_fit_prepared(
    prepared = prepared,
    lambda = lambda,
    penalty = penalty,
    solver = solver,
    linear_init = linear_init,
    n_starts = n_starts,
    seed = seed,
    max_iter = max_iter,
    tol = tol,
    objective_tol = objective_tol,
    numerical_tol = numerical_tol,
    verbose = verbose
  )
  fit$call <- match.call()
  fit
}

predict.kernel_era <- function(object,
                               newX = NULL,
                               type = c("response", "standardized", "components", "all"),
                               ...) {
  type <- match.arg(type)
  Ksets <- length(object$nvar)

  if (is.null(newX)) {
    F_new <- object$F
  } else {
    newX <- .kera_as_matrix(newX)
    if (ncol(newX) != sum(object$nvar)) stop("newX has incorrect columns.", call. = FALSE)
    X_new_blocks <- Map(function(ii, sc) {
      kera_apply_scaler(newX[, ii, drop = FALSE], sc)
    }, object$block_indices, object$x_scalers)

    F_new <- matrix(NA_real_, nrow(newX), Ksets)
    for (k in seq_len(Ksets)) {
      K_raw_cross <- kera_raw_kernel(
        A = X_new_blocks[[k]],
        B = object$X_train_std[[k]],
        kernel = object$kernel,
        sigma = object$sigma[k],
        degree = object$degree,
        offset = object$offset,
        gaussian_denominator = object$gaussian_denominator
      )
      K_cross <- kera_center_cross_kernel(K_raw_cross, object$kernel_objects[[k]])
      F_new[, k] <- as.numeric(K_cross %*% object$alpha[[k]])
    }
  }

  Y_std <- F_new %*% object$B
  Y_response <- kera_inverse_scaler(Y_std, object$y_scaler)

  if (type == "components") return(F_new)
  if (type == "standardized") return(Y_std)
  if (type == "response") return(if (ncol(Y_response) == 1L) as.numeric(Y_response) else Y_response)
  list(response = if (ncol(Y_response) == 1L) as.numeric(Y_response) else Y_response,
       standardized = Y_std,
       components = F_new)
}

print.kernel_era <- function(x, ...) {
  shared_sigma <- length(unique(x$sigma)) == 1L
  cat("Kernel Extended Redundancy Analysis\n")
  cat("  N:", x$n, "  predictor sets:", length(x$nvar), "  outcomes:", x$q, "\n")
  cat("  kernel:", x$kernel, "  lambda:", format(x$lambda, digits = 5), "\n")
  if (x$kernel == "gaussian") {
    cat("  sigma:", paste(format(x$sigma, digits = 5), collapse = ", "),
        if (shared_sigma) "(shared)" else "(set-specific)", "\n")
    cat("  Gaussian convention:", x$gaussian_denominator, "\n")
  }
  cat("  penalty:", x$penalty, "  solver:", x$solver, "\n")
  cat("  converged:", x$converged, "  iterations:", x$iterations, "\n")
  cat("  FIT:", format(x$fit, digits = 5),
      "  objective:", format(x$objective, digits = 7), "\n")
  invisible(x)
}

# -----------------------------------------------------------------------------
# Sign orientation helpers
# -----------------------------------------------------------------------------

align_kernel_era_signs <- function(model, reference_components, X_reference = NULL) {
  if (!inherits(model, "kernel_era")) stop("model must be a kernel_era object.", call. = FALSE)
  reference_components <- .kera_as_matrix(reference_components)
  estimated <- if (is.null(X_reference)) model$F else predict(model, X_reference, "components")
  if (!all(dim(reference_components) == dim(estimated))) {
    stop("reference_components must match the estimated component matrix.", call. = FALSE)
  }

  signs <- rep(1, ncol(estimated))
  for (k in seq_len(ncol(estimated))) {
    rr <- .kera_safe_cor(estimated[, k], reference_components[, k])
    if (is.finite(rr) && rr < 0) signs[k] <- -1
  }
  for (k in seq_along(signs)) {
    if (signs[k] < 0) {
      model$alpha[[k]] <- -model$alpha[[k]]
      model$F[, k] <- -model$F[, k]
      model$B[k, ] <- -model$B[k, ]
    }
  }
  model$fitted <- model$F %*% model$B
  model$residuals <- model$Y_train_std - model$fitted
  model$component_correlation <- crossprod(model$F) / model$n
  model$orientation_signs <- signs
  model
}

align_linear_era_signs <- function(model, reference_components, X_reference = NULL) {
  if (!inherits(model, "linear_era")) stop("model must be a linear_era object.", call. = FALSE)
  reference_components <- .kera_as_matrix(reference_components)
  estimated <- if (is.null(X_reference)) model$F else predict(model, X_reference, "components")
  if (!all(dim(reference_components) == dim(estimated))) {
    stop("reference_components must match the estimated component matrix.", call. = FALSE)
  }

  signs <- rep(1, ncol(estimated))
  for (k in seq_len(ncol(estimated))) {
    rr <- .kera_safe_cor(estimated[, k], reference_components[, k])
    if (is.finite(rr) && rr < 0) signs[k] <- -1
  }
  for (k in seq_along(signs)) {
    if (signs[k] < 0) {
      model$W[[k]] <- -model$W[[k]]
      model$F[, k] <- -model$F[, k]
      model$B[k, ] <- -model$B[k, ]
    }
  }
  model$fitted <- model$F %*% model$B
  model$residuals <- model$Y_train_std - model$fitted
  model$component_correlation <- crossprod(model$F) / model$n
  model$orientation_signs <- signs
  model
}

# -----------------------------------------------------------------------------
# K-fold cross-validation tuning
# -----------------------------------------------------------------------------

kera_make_folds <- function(n, v = 5L, seed = 1) {
  v <- as.integer(v)
  if (v < 2L || v > n) stop("v must be between 2 and n.", call. = FALSE)
  set.seed(seed)
  sample(rep(seq_len(v), length.out = n))
}

kernel_era_cv <- function(X,
                          Y,
                          nvar,
                          lambda_grid,
                          sigma_grid = 1,
                          v = 5L,
                          fold_id = NULL,
                          seed = 1,
                          kernel = c("gaussian", "polynomial", "linear"),
                          degree = 2,
                          offset = 1,
                          gaussian_denominator = c("two_sigma_squared", "sigma_squared"),
                          penalty = c("rkhs", "dual_l2"),
                          solver = c("spectral", "direct"),
                          linear_init = c("ones", "random", "pca"),
                          max_iter = 200L,
                          tol = 1e-7,
                          objective_tol = 1e-9,
                          numerical_tol = 1e-10,
                          spectral_tol = 1e-10,
                          selection = c("one_se", "minimum"),
                          require_all_folds_converged = TRUE,
                          verbose = FALSE) {
  kernel <- match.arg(kernel)
  gaussian_denominator <- match.arg(gaussian_denominator)
  penalty <- match.arg(penalty)
  solver <- match.arg(solver)
  linear_init <- match.arg(linear_init)
  selection <- match.arg(selection)
  require_all_folds_converged <- isTRUE(require_all_folds_converged)

  X <- .kera_as_matrix(X)
  Y <- .kera_as_matrix(Y)
  nvar <- .kera_validate_nvar(nvar, ncol(X))
  if (nrow(X) != nrow(Y)) stop("X and Y must have the same rows.", call. = FALSE)

  lambda_grid <- sort(unique(as.numeric(lambda_grid)))
  sigma_grid <- sort(unique(as.numeric(sigma_grid)))
  if (any(!is.finite(lambda_grid)) || any(lambda_grid < 0)) {
    stop("lambda_grid must contain nonnegative finite values.", call. = FALSE)
  }
  if (kernel == "gaussian" &&
      (any(!is.finite(sigma_grid)) || any(sigma_grid <= 0))) {
    stop("sigma_grid must contain positive finite values.", call. = FALSE)
  }
  if (kernel != "gaussian") sigma_grid <- sigma_grid[1L]

  n <- nrow(X)
  if (is.null(fold_id)) fold_id <- kera_make_folds(n, v = v, seed = seed)
  if (length(fold_id) != n || anyNA(fold_id)) {
    stop("Invalid fold_id.", call. = FALSE)
  }
  fold_levels <- sort(unique(fold_id))
  v <- length(fold_levels)

  grid <- expand.grid(
    lambda = lambda_grid,
    sigma = sigma_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  fold_mse <- matrix(NA_real_, nrow(grid), v)
  convergence <- matrix(FALSE, nrow(grid), v)

  for (ff in seq_along(fold_levels)) {
    val <- which(fold_id == fold_levels[ff])
    trn <- setdiff(seq_len(n), val)
    if (verbose) cat("CV fold", ff, "of", v, "
")
    warm_fold <- NULL

    for (ss in seq_along(sigma_grid)) {
      sigma_value <- sigma_grid[ss]
      prepared <- kera_prepare_data(
        X = X[trn, , drop = FALSE],
        Y = Y[trn, , drop = FALSE],
        nvar = nvar,
        kernel = kernel,
        sigma = sigma_value,
        degree = degree,
        offset = offset,
        gaussian_denominator = gaussian_denominator,
        standardize = TRUE,
        spectral_tol = spectral_tol
      )

      if (is.null(warm_fold)) {
        warm_fold <- kera_linear_era_core(
          prepared$X_blocks,
          prepared$Y_std,
          init = linear_init,
          seed = seed + 5000L * ff,
          max_iter = max_iter,
          tol = tol,
          numerical_tol = numerical_tol,
          verbose = FALSE
        )
      }

      rows <- which(grid$sigma == sigma_value)
      for (rr in rows) {
        fit <- try(
          kera_fit_prepared(
            prepared = prepared,
            lambda = grid$lambda[rr],
            penalty = penalty,
            solver = solver,
            linear_init = linear_init,
            n_starts = 1L,
            seed = seed + 10000L * ff + rr,
            max_iter = max_iter,
            tol = tol,
            objective_tol = objective_tol,
            numerical_tol = numerical_tol,
            warm_start = warm_fold,
            verbose = FALSE
          ),
          silent = TRUE
        )

        if (!inherits(fit, "try-error")) {
          pred <- try(
            predict(fit, X[val, , drop = FALSE], type = "response"),
            silent = TRUE
          )
          if (!inherits(pred, "try-error")) {
            pred <- .kera_as_matrix(pred)
            truth <- Y[val, , drop = FALSE]
            fold_mse[rr, ff] <- mean((truth - pred)^2)
            convergence[rr, ff] <- isTRUE(fit$converged)
          }
        }
      }
    }
  }

  grid$cv_mse <- apply(fold_mse, 1L, function(z) {
    z <- z[is.finite(z)]
    if (length(z) == 0L) return(Inf)
    mean(z)
  })
  grid$cv_se <- apply(fold_mse, 1L, function(z) {
    z <- z[is.finite(z)]
    if (length(z) <= 1L) return(NA_real_)
    stats::sd(z) / sqrt(length(z))
  })
  grid$n_successful_folds <- rowSums(is.finite(fold_mse))
  grid$n_converged_folds <- rowSums(convergence)
  grid$convergence_rate <- grid$n_converged_folds / v
  grid$complete_folds <- grid$n_successful_folds == v

  # Candidates with incomplete validation results are never comparable.
  # By default, a candidate is eligible only when all folds also converged.
  grid$eligible <- grid$complete_folds
  if (require_all_folds_converged) {
    grid$eligible <- grid$eligible & grid$n_converged_folds == v
  }

  comparable <- which(grid$complete_folds & is.finite(grid$cv_mse))
  if (length(comparable) == 0L) {
    stop(
      "No tuning combination produced finite validation errors in every fold.",
      call. = FALSE
    )
  }

  unrestricted_minimum_index <- comparable[
    which.min(grid$cv_mse[comparable])
  ]

  eligible_indices <- which(grid$eligible & is.finite(grid$cv_mse))
  if (length(eligible_indices) == 0L) {
    best_rate <- max(grid$convergence_rate[grid$complete_folds], na.rm = TRUE)
    stop(
      paste0(
        "No tuning combination satisfied the CV eligibility rule. ",
        "The highest all-candidate convergence rate was ",
        format(best_rate, digits = 4),
        ". Increase max_iter, relax tol only if scientifically justified, ",
        "or inspect the tuning grid."
      ),
      call. = FALSE
    )
  }

  eligible_minimum_index <- eligible_indices[
    which.min(grid$cv_mse[eligible_indices])
  ]
  selected_index <- eligible_minimum_index
  one_se_threshold <- NA_real_

  if (selection == "one_se" &&
      is.finite(grid$cv_se[eligible_minimum_index])) {
    one_se_threshold <-
      grid$cv_mse[eligible_minimum_index] +
      grid$cv_se[eligible_minimum_index]

    within_one_se <- which(
      grid$eligible &
        is.finite(grid$cv_mse) &
        grid$cv_mse <= one_se_threshold
    )

    # Among statistically comparable eligible candidates, prefer stronger
    # ridge regularization and then a smoother/larger Gaussian bandwidth.
    ordering <- order(
      -grid$lambda[within_one_se],
      -grid$sigma[within_one_se]
    )
    selected_index <- within_one_se[ordering[1L]]
  }

  grid$unrestricted_minimum <- FALSE
  grid$unrestricted_minimum[unrestricted_minimum_index] <- TRUE
  grid$eligible_minimum <- FALSE
  grid$eligible_minimum[eligible_minimum_index] <- TRUE
  grid$selected <- FALSE
  grid$selected[selected_index] <- TRUE
  grid$within_one_se <- if (is.finite(one_se_threshold)) {
    grid$eligible & is.finite(grid$cv_mse) &
      grid$cv_mse <= one_se_threshold
  } else {
    grid$eligible_minimum
  }

  ordered_table <- grid[
    order(
      !grid$selected,
      !grid$eligible,
      grid$cv_mse,
      -grid$lambda,
      -grid$sigma
    ),
    ,
    drop = FALSE
  ]
  rownames(ordered_table) <- NULL

  selected_row <- grid[selected_index, , drop = FALSE]
  unrestricted_row <- grid[unrestricted_minimum_index, , drop = FALSE]
  eligible_minimum_row <- grid[eligible_minimum_index, , drop = FALSE]
  rownames(selected_row) <- NULL
  rownames(unrestricted_row) <- NULL
  rownames(eligible_minimum_row) <- NULL

  ans <- list(
    best_lambda = grid$lambda[selected_index],
    best_sigma = grid$sigma[selected_index],
    selected_index = selected_index,
    minimum_index = eligible_minimum_index,
    unrestricted_minimum_index = unrestricted_minimum_index,
    selection = selection,
    require_all_folds_converged = require_all_folds_converged,
    one_se_threshold = one_se_threshold,
    selected_row = selected_row,
    minimum_row = eligible_minimum_row,
    unrestricted_minimum_row = unrestricted_row,
    n_eligible_candidates = length(eligible_indices),
    n_total_candidates = nrow(grid),
    table = ordered_table,
    fold_mse = fold_mse,
    fold_convergence = convergence,
    fold_id = fold_id,
    kernel = kernel,
    penalty = penalty,
    gaussian_denominator = gaussian_denominator
  )
  class(ans) <- "kernel_era_cv"
  ans
}

print.kernel_era_cv <- function(x, ...) {
  cat("Kernel ERA cross-validation
")
  cat("  selected lambda:", format(x$best_lambda, digits = 5), "
")
  if (x$kernel == "gaussian") {
    cat("  selected sigma:", format(x$best_sigma, digits = 5), "
")
  }
  cat("  selection rule:", x$selection, "
")
  cat(
    "  eligible candidates:",
    x$n_eligible_candidates,
    "of",
    x$n_total_candidates,
    "
"
  )
  cat(
    "  selected CV MSE:",
    format(x$selected_row$cv_mse, digits = 6),
    "  convergence rate:",
    format(x$selected_row$convergence_rate, digits = 4),
    "
"
  )

  if (!isTRUE(x$unrestricted_minimum_row$eligible)) {
    cat(
      "  unrestricted minimum was excluded:",
      "lambda =", format(x$unrestricted_minimum_row$lambda, digits = 5),
      ", sigma =", format(x$unrestricted_minimum_row$sigma, digits = 5),
      ", convergence rate =",
      format(x$unrestricted_minimum_row$convergence_rate, digits = 4),
      "
"
    )
  }

  print(utils::head(x$table, 10L), row.names = FALSE)
  invisible(x)
}

# -----------------------------------------------------------------------------
# MATLAB-style projection GCV, included for replication/diagnosis
# -----------------------------------------------------------------------------

kernel_era_matlab_gcv <- function(X,
                                   Y,
                                   nvar,
                                   lambda_grid,
                                   sigma_grid,
                                   kernel = c("gaussian", "polynomial", "linear"),
                                   degree = 2,
                                   offset = 1,
                                   gaussian_denominator = c("two_sigma_squared", "sigma_squared"),
                                   penalty = c("rkhs", "dual_l2"),
                                   solver = c("spectral", "direct"),
                                   seed = 1,
                                   max_iter = 200L,
                                   tol = 1e-7,
                                   objective_tol = 1e-9,
                                   numerical_tol = 1e-10,
                                   spectral_tol = 1e-10,
                                   verbose = FALSE) {
  kernel <- match.arg(kernel)
  gaussian_denominator <- match.arg(gaussian_denominator)
  penalty <- match.arg(penalty)
  solver <- match.arg(solver)
  X <- .kera_as_matrix(X)
  Y <- .kera_as_matrix(Y)
  if (ncol(Y) != 1L) {
    warning("The MATLAB routine was written for one outcome; multivariate GCV is averaged here.")
  }

  grid <- expand.grid(lambda = sort(unique(lambda_grid)),
                      sigma = sort(unique(sigma_grid)),
                      KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  grid$gcv <- NA_real_
  grid$df_projection <- NA_real_
  grid$converged <- FALSE
  warm_full <- NULL

  for (ss in unique(grid$sigma)) {
    prepared <- kera_prepare_data(X, Y, nvar,
                                  kernel = kernel, sigma = ss,
                                  degree = degree, offset = offset,
                                  gaussian_denominator = gaussian_denominator,
                                  standardize = TRUE,
                                  spectral_tol = spectral_tol)
    if (is.null(warm_full)) {
      warm_full <- kera_linear_era_core(
        prepared$X_blocks, prepared$Y_std, init = "ones", seed = seed,
        max_iter = max_iter, tol = tol, numerical_tol = numerical_tol,
        verbose = FALSE
      )
    }
    rows <- which(grid$sigma == ss)
    for (rr in rows) {
      fit <- kera_fit_prepared(prepared, lambda = grid$lambda[rr],
                               penalty = penalty, solver = solver,
                               linear_init = "ones", n_starts = 1L,
                               seed = seed + rr,
                               max_iter = max_iter, tol = tol,
                               objective_tol = objective_tol,
                               numerical_tol = numerical_tol,
                               warm_start = warm_full,
                               verbose = FALSE)
      S <- fit$F %*% .kera_pinv_solve(crossprod(fit$F), t(fit$F))
      df <- sum(diag(S))
      denom <- 1 - df / fit$n
      grid$gcv[rr] <- mean((fit$Y_train_std - fit$fitted)^2) / denom^2
      grid$df_projection[rr] <- df
      grid$converged[rr] <- fit$converged
      if (verbose) cat("GCV row", rr, "of", nrow(grid), "\n")
    }
  }

  best <- which.min(grid$gcv)
  ans <- list(best_lambda = grid$lambda[best],
              best_sigma = grid$sigma[best],
              selected_index = best,
              table = grid[order(grid$gcv), , drop = FALSE],
              note = paste(
                "This reproduces the MATLAB projection-GCV formula.",
                "Because F is estimated from Y and tr{F(F'F)^+F'} is usually",
                "approximately the number of components, use K-fold CV as the",
                "primary tuning method for simulation and out-of-sample claims."
              ))
  class(ans) <- "kernel_era_matlab_gcv"
  ans
}

print.kernel_era_matlab_gcv <- function(x, ...) {
  cat("MATLAB-style Kernel ERA projection GCV\n")
  cat("  selected lambda:", format(x$best_lambda, digits = 5), "\n")
  cat("  selected sigma:", format(x$best_sigma, digits = 5), "\n")
  cat("  ", x$note, "\n", sep = "")
  print(utils::head(x$table, 10L), row.names = FALSE)
  invisible(x)
}

# -----------------------------------------------------------------------------
# Bootstrap inference with sign alignment
# -----------------------------------------------------------------------------

kernel_era_bootstrap <- function(model,
                                  X,
                                  Y,
                                  n_boot = 500L,
                                  seed = 1,
                                  conf_level = 0.95,
                                  retune = FALSE,
                                  lambda_grid = NULL,
                                  sigma_grid = NULL,
                                  cv_folds = 5L,
                                  max_iter = 200L,
                                  tol = 1e-7,
                                  objective_tol = 1e-9,
                                  numerical_tol = 1e-10,
                                  verbose = FALSE) {
  if (!inherits(model, "kernel_era")) stop("model must be a kernel_era object.", call. = FALSE)
  X <- .kera_as_matrix(X)
  Y <- .kera_as_matrix(Y)
  if (nrow(X) != model$n || nrow(Y) != model$n) {
    stop("X and Y must be the original training sample used for model.", call. = FALSE)
  }
  n_boot <- as.integer(n_boot)
  if (n_boot < 2L) stop("n_boot must be at least 2.", call. = FALSE)
  if (!is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must lie in (0, 1).", call. = FALSE)
  }
  if (retune && (is.null(lambda_grid) || is.null(sigma_grid))) {
    stop("lambda_grid and sigma_grid are required when retune=TRUE.", call. = FALSE)
  }

  set.seed(seed)
  indices <- replicate(n_boot, sample.int(model$n, model$n, replace = TRUE),
                       simplify = FALSE)
  B_boot <- array(NA_real_, dim = c(n_boot, nrow(model$B), ncol(model$B)))
  convergence <- rep(FALSE, n_boot)
  tuning <- data.frame(lambda = rep(NA_real_, n_boot), sigma = rep(NA_real_, n_boot))
  reference_F <- model$F

  for (b in seq_len(n_boot)) {
    idx <- indices[[b]]
    lambda_b <- model$lambda
    sigma_b <- model$sigma[1L]

    if (retune) {
      cv <- kernel_era_cv(
        X = X[idx, , drop = FALSE],
        Y = Y[idx, , drop = FALSE],
        nvar = model$nvar,
        lambda_grid = lambda_grid,
        sigma_grid = sigma_grid,
        v = cv_folds,
        seed = seed + 100000L + b,
        kernel = model$kernel,
        degree = model$degree,
        offset = model$offset,
        gaussian_denominator = model$gaussian_denominator,
        penalty = model$penalty,
        solver = model$solver,
        max_iter = max_iter,
        tol = tol,
        objective_tol = objective_tol,
        numerical_tol = numerical_tol,
        selection = "minimum"
      )
      lambda_b <- cv$best_lambda
      sigma_b <- cv$best_sigma
    }

    fit_b <- try(kernel_era_fit(
      X = X[idx, , drop = FALSE],
      Y = Y[idx, , drop = FALSE],
      nvar = model$nvar,
      lambda = lambda_b,
      kernel = model$kernel,
      sigma = sigma_b,
      degree = model$degree,
      offset = model$offset,
      gaussian_denominator = model$gaussian_denominator,
      penalty = model$penalty,
      solver = model$solver,
      standardize = TRUE,
      linear_init = "ones",
      n_starts = 1L,
      seed = seed + b,
      max_iter = max_iter,
      tol = tol,
      objective_tol = objective_tol,
      numerical_tol = numerical_tol,
      verbose = FALSE
    ), silent = TRUE)

    if (!inherits(fit_b, "try-error")) {
      # Evaluate bootstrap components on the original observations, then orient
      # each predictor-set component to the original fitted solution.
      F_on_original <- predict(fit_b, X, type = "components")
      signs <- rep(1, ncol(F_on_original))
      for (k in seq_len(ncol(F_on_original))) {
        rr <- .kera_safe_cor(F_on_original[, k], reference_F[, k])
        if (is.finite(rr) && rr < 0) signs[k] <- -1
      }
      B_oriented <- fit_b$B * signs
      B_boot[b, , ] <- B_oriented
      convergence[b] <- fit_b$converged
      tuning$lambda[b] <- lambda_b
      tuning$sigma[b] <- sigma_b
    }
    if (verbose && (b %% max(1L, floor(n_boot / 20L)) == 0L)) {
      cat("Bootstrap", b, "of", n_boot, "\n")
    }
  }

  alpha_tail <- (1 - conf_level) / 2
  safe_sd <- function(z) {
    z <- z[is.finite(z)]
    if (length(z) < 2L) return(NA_real_)
    stats::sd(z)
  }
  safe_quantile <- function(z, prob) {
    z <- z[is.finite(z)]
    if (length(z) < 1L) return(NA_real_)
    as.numeric(stats::quantile(z, probs = prob, names = FALSE, type = 7))
  }
  B_se <- apply(B_boot, c(2L, 3L), safe_sd)
  B_lower <- apply(B_boot, c(2L, 3L), safe_quantile, prob = alpha_tail)
  B_upper <- apply(B_boot, c(2L, 3L), safe_quantile, prob = 1 - alpha_tail)
  target_dim <- c(nrow(model$B), ncol(model$B))
  dim(B_se) <- target_dim
  dim(B_lower) <- target_dim
  dim(B_upper) <- target_dim

  list(
    estimate = model$B,
    bootstrap_estimates = B_boot,
    se = B_se,
    lower = B_lower,
    upper = B_upper,
    conf_level = conf_level,
    convergence_rate = mean(convergence),
    convergence = convergence,
    tuning = tuning,
    retune = retune
  )
}
