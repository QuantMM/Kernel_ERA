# =============================================================================
# Kernel ERA 30-replication stability pilot
#
# Primary goals:
#   1. Evaluate the revised tuning procedure over 30 replications per scenario.
#   2. Select the minimum-CV solution only among tuning combinations that
#      converge in every cross-validation fold.
#   3. Reassess convergence, interpolation, tuning-grid boundaries,
#      out-of-sample prediction, and component recovery.
#   4. Examine CV profiles for boundary-selected tuning parameters.
#   5. Summarize replication-level paired comparisons between Kernel ERA and
#      Linear ERA.
#
# Place this file in the same folder as:
#   * kernel_era_core.R
#   * kernel_era_simulation.R
# =============================================================================

getwd()
setwd("C:/Users/kims15/Downloads")

.get_script_dir <- function() {
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) {
    return(dirname(normalizePath(ofile)))
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]))))
  }

  getwd()
}

SCRIPT_DIR <- .get_script_dir()

source(file.path(SCRIPT_DIR, "kernel_era_core.R"))
source(file.path(SCRIPT_DIR, "kernel_era_simulation.R"))

# -----------------------------------------------------------------------------
# Analysis settings
# -----------------------------------------------------------------------------

N_REP <- 30L
N_CORES <- 1L
CV_FOLDS <- 5L
CV_SELECTION <- "minimum"

LAMBDA_GRID <- 10^seq(-4, 2, by = 1)
SIGMA_GRID <- c(0.75, 1.5, 3, 6, 12, 24, 48)

INTERPOLATION_THRESHOLD <- 0.98

PILOT_DESIGN <- make_kernel_era_quick_grid()

RESULTS_DIR <- file.path(
  SCRIPT_DIR,
  "results_pilot_stability_30rep_converged_minimum"
)

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Run the simulation
# -----------------------------------------------------------------------------

pilot_results <- run_kernel_era_simulation(
  design = PILOT_DESIGN,
  n_rep = N_REP,
  results_dir = RESULTS_DIR,

  # The lower boundaries are retained from the verified pilot. The original
  # small-sigma/small-lambda pathology is controlled by requiring convergence
  # in every CV fold rather than by removing valid grid values globally.
  lambda_grid = LAMBDA_GRID,
  sigma_grid = SIGMA_GRID,

  cv_folds = CV_FOLDS,
  cv_selection = CV_SELECTION,
  cv_require_all_folds_converged = TRUE,

  gaussian_denominator = "two_sigma_squared",
  penalty = "rkhs",

  # Start 1 uses the linear ERA warm start. Additional starts use distinct
  # random initializations, and the best converged solution is retained.
  final_n_starts = 3L,

  max_iter = 500L,
  tol = 1e-6,
  objective_tol = 1e-9,
  interpolation_fit_threshold = INTERPOLATION_THRESHOLD,

  bootstrap_reps = 0L,
  include_matlab_gcv = FALSE,

  n_cores = N_CORES,
  overwrite = FALSE,
  verbose = TRUE
)

pilot_summary <- summarize_kernel_era_simulation(pilot_results)

# -----------------------------------------------------------------------------
# Small helper functions used in the diagnostics below
# -----------------------------------------------------------------------------

.mcse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

.mean_ci <- function(x, level = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)

  if (n == 0L) {
    return(c(mean = NA_real_, lower = NA_real_, upper = NA_real_, mcse = NA_real_))
  }

  m <- mean(x)
  se <- .mcse(x)

  if (n <= 1L || !is.finite(se)) {
    return(c(mean = m, lower = NA_real_, upper = NA_real_, mcse = se))
  }

  critical <- stats::qt(1 - (1 - level) / 2, df = n - 1L)
  c(
    mean = m,
    lower = m - critical * se,
    upper = m + critical * se,
    mcse = se
  )
}

.add_within_group_delta <- function(data, value, group_names) {
  split_id <- interaction(data[, group_names, drop = FALSE], drop = TRUE)
  group_min <- ave(data[[value]], split_id, FUN = function(z) {
    if (all(!is.finite(z))) return(rep(NA_real_, length(z)))
    rep(min(z[is.finite(z)]), length(z))
  })
  data$delta_from_profile_min <- data[[value]] - group_min
  data
}

# =============================================================================
# 1. Performance summary
# =============================================================================

cat("\n============================================================\n")
cat("1. PERFORMANCE SUMMARY\n")
cat("============================================================\n")
print(pilot_summary$performance)

# =============================================================================
# 2. Component recovery summary
# =============================================================================

cat("\n============================================================\n")
cat("2. COMPONENT RECOVERY SUMMARY\n")
cat("============================================================\n")
print(pilot_summary$components)

# =============================================================================
# 3. Tuning and convergence summary
# =============================================================================

cat("\n============================================================\n")
cat("3. TUNING AND CONVERGENCE SUMMARY\n")
cat("============================================================\n")
print(pilot_summary$tuning)

# =============================================================================
# 4. Selected tuning combinations
# =============================================================================

cat("\n============================================================\n")
cat("4. SELECTED TUNING COMBINATIONS\n")
cat("============================================================\n")

selected_columns_requested <- c(
  "condition_id",
  "scenario",
  "replication",
  "lambda",
  "sigma",
  "selected_cv_mse",
  "selected_cv_se",
  "selected_cv_convergence_rate",
  "selected_cv_n_successful_folds",
  "selected_cv_eligible",
  "n_eligible_candidates",
  "n_total_candidates",
  "unrestricted_min_lambda",
  "unrestricted_min_sigma",
  "unrestricted_min_cv_mse",
  "unrestricted_min_convergence_rate",
  "unrestricted_min_eligible",
  "final_converged",
  "final_iterations",
  "final_train_fit",
  "final_test_r2",
  "final_train_test_fit_gap",
  "final_interpolation_flag",
  "final_selected_start"
)

selected_columns <- intersect(
  selected_columns_requested,
  names(pilot_results$selected)
)

print(
  pilot_results$selected[, selected_columns, drop = FALSE],
  row.names = FALSE
)

# =============================================================================
# 5. Eligibility check
# =============================================================================

cat("\n============================================================\n")
cat("5. ELIGIBILITY CHECK\n")
cat("============================================================\n")

selected_ineligible <- subset(
  pilot_results$selected,
  !selected_cv_eligible |
    selected_cv_convergence_rate < 1 |
    selected_cv_n_successful_folds < CV_FOLDS
)

if (nrow(selected_ineligible) == 0L) {
  cat(
    "All selected tuning combinations were eligible and converged in all CV folds.\n"
  )
} else {
  cat("The following selected combinations failed the eligibility check:\n")
  print(selected_ineligible, row.names = FALSE)
}

# =============================================================================
# 6. Unrestricted CV minima excluded for nonconvergence
# =============================================================================

cat("\n============================================================\n")
cat("6. UNRESTRICTED CV MINIMA EXCLUDED FOR NONCONVERGENCE\n")
cat("============================================================\n")

excluded_unrestricted_minima <- subset(
  pilot_results$selected,
  !unrestricted_min_eligible,
  select = c(
    condition_id,
    scenario,
    replication,
    unrestricted_min_lambda,
    unrestricted_min_sigma,
    unrestricted_min_cv_mse,
    unrestricted_min_convergence_rate,
    lambda,
    sigma,
    selected_cv_mse
  )
)

if (nrow(excluded_unrestricted_minima) == 0L) {
  cat("No unrestricted CV minimum was excluded.\n")
} else {
  print(excluded_unrestricted_minima, row.names = FALSE)
}

# =============================================================================
# 7. Final-fit nonconvergence
# =============================================================================

cat("\n============================================================\n")
cat("7. FINAL-FIT NONCONVERGENCE\n")
cat("============================================================\n")

nonconverged_final <- subset(
  pilot_results$performance,
  method == "Gaussian Kernel ERA" & !converged,
  select = c(
    condition_id,
    scenario,
    replication,
    iterations,
    lambda,
    sigma,
    train_fit,
    test_r2,
    train_test_fit_gap,
    final_n_starts,
    selected_start
  )
)

if (nrow(nonconverged_final) == 0L) {
  cat("All final Kernel ERA fits converged.\n")
} else {
  print(nonconverged_final, row.names = FALSE)
}

# =============================================================================
# 8. Interpolation check
# =============================================================================

cat("\n============================================================\n")
cat(
  "8. INTERPOLATION CHECK: TRAINING FIT > ",
  INTERPOLATION_THRESHOLD,
  "\n",
  sep = ""
)
cat("============================================================\n")

interpolating_final <- subset(
  pilot_results$performance,
  method == "Gaussian Kernel ERA" & interpolation_flag,
  select = c(
    condition_id,
    scenario,
    replication,
    lambda,
    sigma,
    train_fit,
    test_r2,
    train_test_fit_gap,
    converged
  )
)

if (nrow(interpolating_final) == 0L) {
  cat("No final Kernel ERA fit exceeded the interpolation threshold.\n")
} else {
  print(interpolating_final, row.names = FALSE)
}

# =============================================================================
# 9. Tuning-grid boundary cases
# =============================================================================

cat("\n============================================================\n")
cat("9. TUNING-GRID BOUNDARY CASES\n")
cat("============================================================\n")

boundary_cases <- subset(
  pilot_results$selected,
  lambda_at_lower_boundary |
    lambda_at_upper_boundary |
    sigma_at_lower_boundary |
    sigma_at_upper_boundary,
  select = c(
    condition_id,
    scenario,
    replication,
    lambda,
    sigma,
    lambda_at_lower_boundary,
    lambda_at_upper_boundary,
    sigma_at_lower_boundary,
    sigma_at_upper_boundary,
    final_converged,
    final_interpolation_flag
  )
)

if (nrow(boundary_cases) == 0L) {
  cat("No selected tuning combination was on a grid boundary.\n")
} else {
  print(boundary_cases, row.names = FALSE)
}

# =============================================================================
# 10. Paired Kernel ERA versus Linear ERA comparison
# =============================================================================

cat("\n============================================================\n")
cat("10. PAIRED KERNEL ERA VERSUS LINEAR ERA COMPARISON\n")
cat("============================================================\n")

perf <- pilot_results$performance

kernel_perf <- subset(
  perf,
  method == "Gaussian Kernel ERA",
  select = c(
    condition_id,
    scenario,
    replication,
    test_mse,
    test_r2,
    prediction_correlation
  )
)

linear_perf <- subset(
  perf,
  method == "Linear ERA",
  select = c(
    condition_id,
    scenario,
    replication,
    test_mse,
    test_r2,
    prediction_correlation
  )
)

paired <- merge(
  kernel_perf,
  linear_perf,
  by = c("condition_id", "scenario", "replication"),
  suffixes = c("_kernel", "_linear")
)

paired$delta_r2 <-
  paired$test_r2_kernel - paired$test_r2_linear

paired$delta_mse <-
  paired$test_mse_kernel - paired$test_mse_linear

paired$delta_prediction_correlation <-
  paired$prediction_correlation_kernel -
  paired$prediction_correlation_linear

paired$kernel_r2_win <-
  paired$delta_r2 > 0

paired$kernel_mse_win <-
  paired$delta_mse < 0

paired$mse_reduction_proportion <-
  (paired$test_mse_linear - paired$test_mse_kernel) /
  paired$test_mse_linear

paired_summary <- do.call(
  rbind,
  lapply(split(paired, paired$scenario), function(d) {
    delta_r2_ci <- .mean_ci(d$delta_r2)
    delta_mse_ci <- .mean_ci(d$delta_mse)
    delta_cor_ci <- .mean_ci(d$delta_prediction_correlation)

    data.frame(
      scenario = d$scenario[1L],
      n_rep = nrow(d),

      mean_delta_r2 = delta_r2_ci["mean"],
      sd_delta_r2 = stats::sd(d$delta_r2),
      mcse_delta_r2 = delta_r2_ci["mcse"],
      delta_r2_ci_lower = delta_r2_ci["lower"],
      delta_r2_ci_upper = delta_r2_ci["upper"],

      mean_delta_mse = delta_mse_ci["mean"],
      mcse_delta_mse = delta_mse_ci["mcse"],
      delta_mse_ci_lower = delta_mse_ci["lower"],
      delta_mse_ci_upper = delta_mse_ci["upper"],

      mean_delta_prediction_correlation = delta_cor_ci["mean"],
      mcse_delta_prediction_correlation = delta_cor_ci["mcse"],

      kernel_r2_win_rate = mean(d$kernel_r2_win),
      kernel_mse_win_rate = mean(d$kernel_mse_win),
      mean_mse_reduction = mean(d$mse_reduction_proportion),
      stringsAsFactors = FALSE
    )
  })
)

print(paired_summary, row.names = FALSE)

cat("\nReplication-level paired differences:\n")
print(
  paired[
    order(paired$condition_id, paired$replication),
    c(
      "condition_id",
      "scenario",
      "replication",
      "test_r2_kernel",
      "test_r2_linear",
      "delta_r2",
      "test_mse_kernel",
      "test_mse_linear",
      "delta_mse",
      "mse_reduction_proportion"
    )
  ],
  row.names = FALSE
)

# =============================================================================
# 11. Paired component-recovery comparison
# =============================================================================

cat("\n============================================================\n")
cat("11. PAIRED COMPONENT-RECOVERY COMPARISON\n")
cat("============================================================\n")

component_data <- pilot_results$components

kernel_components <- subset(
  component_data,
  method == "Gaussian Kernel ERA",
  select = c(
    condition_id,
    scenario,
    replication,
    predictor_set,
    test_component_correlation,
    test_component_rmse_calibrated
  )
)

linear_components <- subset(
  component_data,
  method == "Linear ERA",
  select = c(
    condition_id,
    scenario,
    replication,
    predictor_set,
    test_component_correlation,
    test_component_rmse_calibrated
  )
)

component_paired <- merge(
  kernel_components,
  linear_components,
  by = c(
    "condition_id",
    "scenario",
    "replication",
    "predictor_set"
  ),
  suffixes = c("_kernel", "_linear")
)

component_paired$delta_component_correlation <-
  component_paired$test_component_correlation_kernel -
  component_paired$test_component_correlation_linear

component_paired$delta_component_rmse <-
  component_paired$test_component_rmse_calibrated_kernel -
  component_paired$test_component_rmse_calibrated_linear

component_paired$kernel_correlation_win <-
  component_paired$delta_component_correlation > 0

component_paired$kernel_rmse_win <-
  component_paired$delta_component_rmse < 0

component_paired_summary <- do.call(
  rbind,
  lapply(
    split(
      component_paired,
      list(
        component_paired$scenario,
        component_paired$predictor_set,
        drop = TRUE
      )
    ),
    function(d) {
      cor_ci <- .mean_ci(d$delta_component_correlation)
      rmse_ci <- .mean_ci(d$delta_component_rmse)

      data.frame(
        scenario = d$scenario[1L],
        predictor_set = d$predictor_set[1L],
        n_rep = nrow(d),

        mean_delta_component_correlation = cor_ci["mean"],
        mcse_delta_component_correlation = cor_ci["mcse"],
        delta_component_correlation_ci_lower = cor_ci["lower"],
        delta_component_correlation_ci_upper = cor_ci["upper"],
        kernel_component_correlation_win_rate =
          mean(d$kernel_correlation_win),

        mean_delta_component_rmse = rmse_ci["mean"],
        mcse_delta_component_rmse = rmse_ci["mcse"],
        delta_component_rmse_ci_lower = rmse_ci["lower"],
        delta_component_rmse_ci_upper = rmse_ci["upper"],
        kernel_component_rmse_win_rate =
          mean(d$kernel_rmse_win),
        stringsAsFactors = FALSE
      )
    }
  )
)

component_paired_summary <- component_paired_summary[
  order(
    component_paired_summary$scenario,
    component_paired_summary$predictor_set
  ),
]

print(component_paired_summary, row.names = FALSE)

# =============================================================================
# 12. Boundary CV profiles
# =============================================================================

cat("\n============================================================\n")
cat("12. BOUNDARY CV PROFILES\n")
cat("============================================================\n")

lambda_profile <- NULL
sigma_profile <- NULL

if (nrow(boundary_cases) == 0L) {
  cat("No boundary-selected cases were available for CV-profile analysis.\n")
} else {
  boundary_keys <- unique(
    boundary_cases[, c("condition_id", "scenario", "replication")]
  )

  tuning_data <- pilot_results$tuning

  if (!"scenario" %in% names(tuning_data)) {
    tuning_data <- merge(
      tuning_data,
      unique(pilot_results$selected[, c(
        "condition_id",
        "scenario",
        "replication"
      )]),
      by = c("condition_id", "replication"),
      all.x = TRUE
    )
  }

  if ("eligible" %in% names(tuning_data)) {
    tuning_data <- tuning_data[tuning_data$eligible, , drop = FALSE]
  } else if ("selected_cv_eligible" %in% names(tuning_data)) {
    tuning_data <- tuning_data[
      tuning_data$selected_cv_eligible,
      ,
      drop = FALSE
    ]
  }

  boundary_tuning <- merge(
    tuning_data,
    boundary_keys,
    by = c("condition_id", "scenario", "replication")
  )

  selected_lookup <- pilot_results$selected[, c(
    "condition_id",
    "scenario",
    "replication",
    "lambda",
    "sigma"
  )]
  names(selected_lookup)[4:5] <- c(
    "selected_lambda",
    "selected_sigma"
  )

  lambda_profile <- stats::aggregate(
    cv_mse ~ condition_id + scenario + replication + lambda,
    data = boundary_tuning,
    FUN = min
  )

  lambda_profile <- merge(
    lambda_profile,
    selected_lookup,
    by = c("condition_id", "scenario", "replication"),
    all.x = TRUE
  )

  lambda_profile$selected_lambda_value <-
    lambda_profile$lambda == lambda_profile$selected_lambda

  lambda_profile <- .add_within_group_delta(
    lambda_profile,
    value = "cv_mse",
    group_names = c("condition_id", "scenario", "replication")
  )

  lambda_profile <- lambda_profile[
    order(
      lambda_profile$condition_id,
      lambda_profile$replication,
      lambda_profile$lambda
    ),
  ]

  sigma_profile <- stats::aggregate(
    cv_mse ~ condition_id + scenario + replication + sigma,
    data = boundary_tuning,
    FUN = min
  )

  sigma_profile <- merge(
    sigma_profile,
    selected_lookup,
    by = c("condition_id", "scenario", "replication"),
    all.x = TRUE
  )

  sigma_profile$selected_sigma_value <-
    sigma_profile$sigma == sigma_profile$selected_sigma

  sigma_profile <- .add_within_group_delta(
    sigma_profile,
    value = "cv_mse",
    group_names = c("condition_id", "scenario", "replication")
  )

  sigma_profile <- sigma_profile[
    order(
      sigma_profile$condition_id,
      sigma_profile$replication,
      sigma_profile$sigma
    ),
  ]

  cat("\nMinimum eligible CV MSE at each lambda:\n")
  print(lambda_profile, row.names = FALSE)

  cat("\nMinimum eligible CV MSE at each sigma:\n")
  print(sigma_profile, row.names = FALSE)
}

# =============================================================================
# 13. Automated stability checks
# =============================================================================

cat("\n============================================================\n")
cat("13. AUTOMATED STABILITY CHECKS\n")
cat("============================================================\n")

stability_checks <- merge(
  pilot_summary$tuning[, c(
    "condition_id",
    "scenario",
    "final_convergence_rate",
    "final_interpolation_rate",
    "lambda_lower_boundary_rate",
    "lambda_upper_boundary_rate",
    "sigma_lower_boundary_rate",
    "sigma_upper_boundary_rate"
  )],
  paired_summary[, c(
    "scenario",
    "mean_delta_r2",
    "mcse_delta_r2",
    "delta_r2_ci_lower",
    "delta_r2_ci_upper",
    "kernel_r2_win_rate",
    "mean_mse_reduction"
  )],
  by = "scenario",
  all = TRUE
)

stability_checks$convergence_at_least_95_percent <-
  stability_checks$final_convergence_rate >= 0.95

stability_checks$interpolation_at_most_5_percent <-
  stability_checks$final_interpolation_rate <= 0.05

stability_checks$boundary_rate_warning <-
  pmax(
    stability_checks$lambda_lower_boundary_rate,
    stability_checks$lambda_upper_boundary_rate,
    stability_checks$sigma_lower_boundary_rate,
    stability_checks$sigma_upper_boundary_rate,
    na.rm = TRUE
  ) >= 0.20

print(stability_checks, row.names = FALSE)

# =============================================================================
# 14. Task errors
# =============================================================================

cat("\n============================================================\n")
cat("14. TASK ERRORS\n")
cat("============================================================\n")

if (is.null(pilot_results$error) || nrow(pilot_results$error) == 0L) {
  cat("No task-level errors were recorded.\n")
} else {
  print(pilot_results$error, row.names = FALSE)
}

# =============================================================================
# 15. Save combined outputs
# =============================================================================

cat("\n============================================================\n")
cat("15. SAVING COMBINED OUTPUTS\n")
cat("============================================================\n")

saveRDS(
  pilot_results,
  file.path(RESULTS_DIR, "pilot_results_30rep_combined.rds")
)

saveRDS(
  pilot_summary,
  file.path(RESULTS_DIR, "pilot_summary_30rep_combined.rds")
)

utils::write.csv(
  pilot_summary$performance,
  file.path(RESULTS_DIR, "performance_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  pilot_summary$components,
  file.path(RESULTS_DIR, "component_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  pilot_summary$tuning,
  file.path(RESULTS_DIR, "tuning_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  pilot_results$selected,
  file.path(RESULTS_DIR, "selected_tuning_by_replication.csv"),
  row.names = FALSE
)

utils::write.csv(
  paired,
  file.path(RESULTS_DIR, "paired_prediction_results_by_replication.csv"),
  row.names = FALSE
)

utils::write.csv(
  paired_summary,
  file.path(RESULTS_DIR, "paired_prediction_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  component_paired,
  file.path(RESULTS_DIR, "paired_component_results_by_replication.csv"),
  row.names = FALSE
)

utils::write.csv(
  component_paired_summary,
  file.path(RESULTS_DIR, "paired_component_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  boundary_cases,
  file.path(RESULTS_DIR, "boundary_cases.csv"),
  row.names = FALSE
)

if (!is.null(lambda_profile)) {
  utils::write.csv(
    lambda_profile,
    file.path(RESULTS_DIR, "boundary_lambda_profiles.csv"),
    row.names = FALSE
  )
}

if (!is.null(sigma_profile)) {
  utils::write.csv(
    sigma_profile,
    file.path(RESULTS_DIR, "boundary_sigma_profiles.csv"),
    row.names = FALSE
  )
}

utils::write.csv(
  stability_checks,
  file.path(RESULTS_DIR, "stability_checks.csv"),
  row.names = FALSE
)

if (nrow(nonconverged_final) > 0L) {
  utils::write.csv(
    nonconverged_final,
    file.path(RESULTS_DIR, "nonconverged_final_fits.csv"),
    row.names = FALSE
  )
}

if (nrow(interpolating_final) > 0L) {
  utils::write.csv(
    interpolating_final,
    file.path(RESULTS_DIR, "interpolating_final_fits.csv"),
    row.names = FALSE
  )
}

if (!is.null(pilot_results$error) && nrow(pilot_results$error) > 0L) {
  utils::write.csv(
    pilot_results$error,
    file.path(RESULTS_DIR, "task_errors.csv"),
    row.names = FALSE
  )
}

cat("All combined outputs were saved to:\n", RESULTS_DIR, "\n")
