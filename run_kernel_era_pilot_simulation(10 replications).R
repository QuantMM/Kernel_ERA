# =============================================================================
# Kernel ERA revised pilot simulation
#
# Purpose:
#   1. Verify that only tuning combinations with complete, converged CV fits
#      can be selected.
#   2. Use the one-standard-error rule to avoid selecting unnecessarily
#      flexible solutions.
#   3. Recheck tuning-grid boundaries, final convergence, interpolation, and
#      out-of-sample performance using the same three pilot scenarios.
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

# Keep 10 replications for this verification pilot so that the results can be
# compared directly with the original 10-replication pilot using the same seeds.
# After the revised procedure is verified, increase N_REP to 30 for a larger
# stability pilot.
N_REP <- 10L
N_CORES <- 1L

PILOT_DESIGN <- make_kernel_era_quick_grid()

RESULTS_DIR <- file.path(
  SCRIPT_DIR,
  "results_pilot_revised_converged_one_se"
)

pilot_results <- run_kernel_era_simulation(
  design = PILOT_DESIGN,
  n_rep = N_REP,
  results_dir = RESULTS_DIR,

  # Do not extend the lower boundaries. The original small-sigma/small-lambda
  # corner produced interpolation and nonconvergence. Extend only the upper
  # boundaries to evaluate whether stronger smoothing or regularization is
  # needed.
  lambda_grid = 10^seq(-4, 2, by = 1),
  sigma_grid = c(0.75, 1.5, 3, 6, 12, 24, 48),

  cv_folds = 5L,
  cv_selection = "one_se",
  cv_require_all_folds_converged = TRUE,

  gaussian_denominator = "two_sigma_squared",
  penalty = "rkhs",

  # The first final-model start uses the linear ERA warm start. The additional
  # starts use random initializations and the best converged solution is kept.
  final_n_starts = 3L,

  max_iter = 500L,
  tol = 1e-6,
  objective_tol = 1e-9,
  interpolation_fit_threshold = 0.98,

  bootstrap_reps = 0L,
  include_matlab_gcv = FALSE,

  n_cores = N_CORES,
  overwrite = FALSE,
  verbose = TRUE
)

pilot_summary <- summarize_kernel_era_simulation(pilot_results)

cat("\n============================================================\n")
cat("1. PERFORMANCE SUMMARY\n")
cat("============================================================\n")
print(pilot_summary$performance)

cat("\n============================================================\n")
cat("2. COMPONENT RECOVERY SUMMARY\n")
cat("============================================================\n")
print(pilot_summary$components)

cat("\n============================================================\n")
cat("3. TUNING AND CONVERGENCE SUMMARY\n")
cat("============================================================\n")
print(pilot_summary$tuning)

cat("\n============================================================\n")
cat("4. SELECTED TUNING COMBINATIONS\n")
cat("============================================================\n")
selected_columns <- c(
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
print(
  pilot_results$selected[, selected_columns, drop = FALSE],
  row.names = FALSE
)

cat("\n============================================================\n")
cat("5. ELIGIBILITY CHECK\n")
cat("============================================================\n")
selected_ineligible <- subset(
  pilot_results$selected,
  !selected_cv_eligible |
    selected_cv_convergence_rate < 1 |
    selected_cv_n_successful_folds < 5
)

if (nrow(selected_ineligible) == 0L) {
  cat("All selected tuning combinations were eligible and converged in all CV folds.\n")
} else {
  cat("The following selected combinations failed the eligibility check:\n")
  print(selected_ineligible)
}

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

cat("\n============================================================\n")
cat("8. INTERPOLATION CHECK: TRAINING FIT > 0.98\n")
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

paired$kernel_r2_win <-
  paired$delta_r2 > 0

paired$mse_reduction_proportion <-
  (paired$test_mse_linear - paired$test_mse_kernel) /
  paired$test_mse_linear

paired_summary <- do.call(
  rbind,
  lapply(split(paired, paired$scenario), function(d) {
    data.frame(
      scenario = d$scenario[1L],
      n_rep = nrow(d),
      mean_delta_r2 = mean(d$delta_r2),
      sd_delta_r2 = stats::sd(d$delta_r2),
      mcse_delta_r2 =
        stats::sd(d$delta_r2) / sqrt(nrow(d)),
      kernel_r2_win_rate =
        mean(d$kernel_r2_win),
      mean_mse_reduction =
        mean(d$mse_reduction_proportion),
      stringsAsFactors = FALSE
    )
  })
)

print(paired_summary, row.names = FALSE)

cat("\n============================================================\n")
cat("11. TASK ERRORS\n")
cat("============================================================\n")
if (is.null(pilot_results$error) || nrow(pilot_results$error) == 0L) {
  cat("No task-level errors were recorded.\n")
} else {
  print(pilot_results$error, row.names = FALSE)
}

cat("\nResults were saved to:\n", RESULTS_DIR, "\n")



# =============================================================================
# Kernel ERA revised pilot simulation
#
# Purpose:
#   1. Verify that only tuning combinations with complete, converged CV fits
#      can be selected.
#   2. Use the one-standard-error rule to avoid selecting unnecessarily
#      flexible solutions.
#   3. Recheck tuning-grid boundaries, final convergence, interpolation, and
#      out-of-sample performance using the same three pilot scenarios.
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

# Keep 10 replications for this verification pilot so that the results can be
# compared directly with the original 10-replication pilot using the same seeds.
# After the revised procedure is verified, increase N_REP to 30 for a larger
# stability pilot.
N_REP <- 10L
N_CORES <- 1L

PILOT_DESIGN <- make_kernel_era_quick_grid()

RESULTS_DIR <- file.path(
  SCRIPT_DIR,
  "results_pilot_revised_converged_minimum"
)

pilot_results <- run_kernel_era_simulation(
  design = PILOT_DESIGN,
  n_rep = N_REP,
  results_dir = RESULTS_DIR,

  # Do not extend the lower boundaries. The original small-sigma/small-lambda
  # corner produced interpolation and nonconvergence. Extend only the upper
  # boundaries to evaluate whether stronger smoothing or regularization is
  # needed.
  lambda_grid = 10^seq(-4, 2, by = 1),
  sigma_grid = c(0.75, 1.5, 3, 6, 12, 24, 48),

  cv_folds = 5L,
  cv_selection = "minimum",
  cv_require_all_folds_converged = TRUE,

  gaussian_denominator = "two_sigma_squared",
  penalty = "rkhs",

  # The first final-model start uses the linear ERA warm start. The additional
  # starts use random initializations and the best converged solution is kept.
  final_n_starts = 3L,

  max_iter = 500L,
  tol = 1e-6,
  objective_tol = 1e-9,
  interpolation_fit_threshold = 0.98,

  bootstrap_reps = 0L,
  include_matlab_gcv = FALSE,

  n_cores = N_CORES,
  overwrite = FALSE,
  verbose = TRUE
)

pilot_summary <- summarize_kernel_era_simulation(pilot_results)

cat("\n============================================================\n")
cat("1. PERFORMANCE SUMMARY\n")
cat("============================================================\n")
print(pilot_summary$performance)

cat("\n============================================================\n")
cat("2. COMPONENT RECOVERY SUMMARY\n")
cat("============================================================\n")
print(pilot_summary$components)

cat("\n============================================================\n")
cat("3. TUNING AND CONVERGENCE SUMMARY\n")
cat("============================================================\n")
print(pilot_summary$tuning)

cat("\n============================================================\n")
cat("4. SELECTED TUNING COMBINATIONS\n")
cat("============================================================\n")
selected_columns <- c(
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
print(
  pilot_results$selected[, selected_columns, drop = FALSE],
  row.names = FALSE
)

cat("\n============================================================\n")
cat("5. ELIGIBILITY CHECK\n")
cat("============================================================\n")
selected_ineligible <- subset(
  pilot_results$selected,
  !selected_cv_eligible |
    selected_cv_convergence_rate < 1 |
    selected_cv_n_successful_folds < 5
)

if (nrow(selected_ineligible) == 0L) {
  cat("All selected tuning combinations were eligible and converged in all CV folds.\n")
} else {
  cat("The following selected combinations failed the eligibility check:\n")
  print(selected_ineligible)
}

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

cat("\n============================================================\n")
cat("8. INTERPOLATION CHECK: TRAINING FIT > 0.98\n")
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

paired$kernel_r2_win <-
  paired$delta_r2 > 0

paired$mse_reduction_proportion <-
  (paired$test_mse_linear - paired$test_mse_kernel) /
  paired$test_mse_linear

paired_summary <- do.call(
  rbind,
  lapply(split(paired, paired$scenario), function(d) {
    data.frame(
      scenario = d$scenario[1L],
      n_rep = nrow(d),
      mean_delta_r2 = mean(d$delta_r2),
      sd_delta_r2 = stats::sd(d$delta_r2),
      mcse_delta_r2 =
        stats::sd(d$delta_r2) / sqrt(nrow(d)),
      kernel_r2_win_rate =
        mean(d$kernel_r2_win),
      mean_mse_reduction =
        mean(d$mse_reduction_proportion),
      stringsAsFactors = FALSE
    )
  })
)

print(paired_summary, row.names = FALSE)

cat("\n============================================================\n")
cat("11. TASK ERRORS\n")
cat("============================================================\n")
if (is.null(pilot_results$error) || nrow(pilot_results$error) == 0L) {
  cat("No task-level errors were recorded.\n")
} else {
  print(pilot_results$error, row.names = FALSE)
}

cat("\nResults were saved to:\n", RESULTS_DIR, "\n")
