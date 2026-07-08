# =============================================================================
# Kernel ERA final main simulation
#
# Primary analysis specification
#   * 5-fold CV
#   * minimum CV error among candidates that converged in every fold
#   * Gaussian kernel exp(-||x-z||^2 / (2 sigma^2))
#   * RKHS-norm ridge penalty
#   * three final-model starts; best converged solution retained
#   * dimension-adjusted sigma grids for P_k = 5 and P_k = 20
#   * resumable condition-replication files
#
# Place this file in the same folder as:
#   * kernel_era_core.R
#   * kernel_era_simulation.R
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Working directory and source files
# -----------------------------------------------------------------------------

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
# 1. User settings
# -----------------------------------------------------------------------------

# Use "grid_check" for a short end-to-end check of all 24 conditions.
# Use "main" for the final 200-replication simulation.
# RUN_MODE <- "main"
RUN_MODE <- "grid_check"

N_REP_GRID_CHECK <- 5L
N_REP_MAIN <- 200L

N_REP <- switch(
  RUN_MODE,
  grid_check = N_REP_GRID_CHECK,
  main = N_REP_MAIN,
  stop("RUN_MODE must be either 'grid_check' or 'main'.", call. = FALSE)
)

# On Windows, parallel execution requires the future and future.apply packages.
# Start conservatively. Increase this only when adequate RAM is available.
N_CORES <- 1L

CV_FOLDS <- 5L
FINAL_N_STARTS <- 3L
MAX_ITER <- 500L
TOL <- 1e-6
OBJECTIVE_TOL <- 1e-9
INTERPOLATION_THRESHOLD <- 0.98

# The pilot showed lower-bound selections at lambda = 1e-4 and upper-bound
# selections at lambda = 1e2. The final grid extends one order in each direction.
# Candidates that fail to converge in any CV fold remain in the diagnostic table
# but are not eligible for selection.
LAMBDA_GRID <- 10^seq(-5, 3, by = 1)

# The P_k = 5 stability pilot frequently selected sigma = 48 under linear truth.
# The base grid therefore extends to 96. For other predictor-set dimensions,
# sigma is multiplied by sqrt(P_k / 5), because Euclidean distances among
# standardized observations grow approximately in proportion to sqrt(P_k).
BASE_P_FOR_SIGMA <- 5L
BASE_SIGMA_GRID <- c(0.75, 1.5, 3, 6, 12, 24, 48, 96)

sigma_grid_for_p <- function(p_per_set) {
  p_per_set <- as.numeric(p_per_set)
  if (length(p_per_set) != 1L || !is.finite(p_per_set) || p_per_set <= 0) {
    stop("p_per_set must be one positive finite number.", call. = FALSE)
  }

  scale_factor <- sqrt(p_per_set / BASE_P_FOR_SIGMA)
  unique(as.numeric(signif(BASE_SIGMA_GRID * scale_factor, digits = 12L)))
}

RESULTS_ROOT <- file.path(
  SCRIPT_DIR,
  if (RUN_MODE == "main") {
    "results_main_final_converged_minimum_dimension_scaled"
  } else {
    "results_main_grid_check_converged_minimum_dimension_scaled"
  }
)

dir.create(RESULTS_ROOT, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Main design
# -----------------------------------------------------------------------------

MAIN_DESIGN <- make_kernel_era_main_grid()
MAIN_DESIGN <- MAIN_DESIGN[order(MAIN_DESIGN$condition_id), , drop = FALSE]
rownames(MAIN_DESIGN) <- NULL

required_design_columns <- c(
  "condition_id",
  "scenario",
  "n_train",
  "n_test",
  "n_sets",
  "p_per_set",
  "rho_within",
  "rho_between",
  "target_r2",
  "seed_base"
)

missing_design_columns <- setdiff(required_design_columns, names(MAIN_DESIGN))
if (length(missing_design_columns)) {
  stop(
    "MAIN_DESIGN is missing: ",
    paste(missing_design_columns, collapse = ", "),
    call. = FALSE
  )
}

if (anyDuplicated(MAIN_DESIGN$condition_id)) {
  stop("condition_id must be unique in MAIN_DESIGN.", call. = FALSE)
}

utils::write.csv(
  MAIN_DESIGN,
  file.path(RESULTS_ROOT, "main_simulation_design.csv"),
  row.names = FALSE
)

# Record the actual sigma grid used for each predictor-set dimension.
p_values <- sort(unique(as.integer(MAIN_DESIGN$p_per_set)))

sigma_grid_table <- do.call(
  rbind,
  lapply(p_values, function(p) {
    values <- sigma_grid_for_p(p)
    data.frame(
      p_per_set = p,
      sigma_index = seq_along(values),
      sigma = values,
      scale_relative_to_p5 = sqrt(p / BASE_P_FOR_SIGMA),
      stringsAsFactors = FALSE
    )
  })
)

utils::write.csv(
  sigma_grid_table,
  file.path(RESULTS_ROOT, "sigma_grids_by_predictor_dimension.csv"),
  row.names = FALSE
)

run_settings <- data.frame(
  setting = c(
    "run_mode",
    "n_rep",
    "n_cores",
    "cv_folds",
    "cv_selection",
    "cv_require_all_folds_converged",
    "final_n_starts",
    "max_iter",
    "tol",
    "objective_tol",
    "interpolation_threshold",
    "gaussian_denominator",
    "penalty",
    "lambda_grid",
    "base_sigma_grid_p5"
  ),
  value = c(
    RUN_MODE,
    as.character(N_REP),
    as.character(N_CORES),
    as.character(CV_FOLDS),
    "minimum",
    "TRUE",
    as.character(FINAL_N_STARTS),
    as.character(MAX_ITER),
    format(TOL, scientific = TRUE),
    format(OBJECTIVE_TOL, scientific = TRUE),
    as.character(INTERPOLATION_THRESHOLD),
    "two_sigma_squared",
    "rkhs",
    paste(LAMBDA_GRID, collapse = ";"),
    paste(BASE_SIGMA_GRID, collapse = ";")
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  run_settings,
  file.path(RESULTS_ROOT, "run_settings.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("KERNEL ERA MAIN SIMULATION\n")
cat("============================================================\n")
cat("Run mode:", RUN_MODE, "\n")
cat("Conditions:", nrow(MAIN_DESIGN), "\n")
cat("Replications per condition:", N_REP, "\n")
cat("Expected condition-replication tasks:", nrow(MAIN_DESIGN) * N_REP, "\n")
cat("Lambda grid:", paste(LAMBDA_GRID, collapse = ", "), "\n")
cat("Results root:", RESULTS_ROOT, "\n")

for (p in p_values) {
  cat(
    "P_k =", p,
    "sigma grid:",
    paste(sigma_grid_for_p(p), collapse = ", "),
    "\n"
  )
}

# -----------------------------------------------------------------------------
# 3. Run the simulation separately for each predictor-set dimension
# -----------------------------------------------------------------------------

# Splitting by p_per_set permits dimension-adjusted sigma grids without changing
# kernel_era_core.R or kernel_era_simulation.R. Each subgroup has its own
# resumable directory, while condition_id remains globally unique.

results_by_p <- vector("list", length(p_values))
names(results_by_p) <- paste0("p", p_values)

for (j in seq_along(p_values)) {
  p <- p_values[j]
  design_p <- MAIN_DESIGN[MAIN_DESIGN$p_per_set == p, , drop = FALSE]
  sigma_grid_p <- sigma_grid_for_p(p)

  results_dir_p <- file.path(
    RESULTS_ROOT,
    sprintf("p_per_set_%02d", p)
  )

  cat("\n============================================================\n")
  cat("RUNNING PREDICTOR-SET DIMENSION P_k =", p, "\n")
  cat("============================================================\n")
  cat("Conditions in this subgroup:", nrow(design_p), "\n")
  cat("Tasks in this subgroup:", nrow(design_p) * N_REP, "\n")
  cat("Sigma grid:", paste(sigma_grid_p, collapse = ", "), "\n")
  cat("Subgroup directory:", results_dir_p, "\n")

  results_by_p[[j]] <- run_kernel_era_simulation(
    design = design_p,
    n_rep = N_REP,
    results_dir = results_dir_p,

    lambda_grid = LAMBDA_GRID,
    sigma_grid = sigma_grid_p,

    cv_folds = CV_FOLDS,
    cv_selection = "minimum",
    cv_require_all_folds_converged = TRUE,

    gaussian_denominator = "two_sigma_squared",
    penalty = "rkhs",

    final_n_starts = FINAL_N_STARTS,

    max_iter = MAX_ITER,
    tol = TOL,
    objective_tol = OBJECTIVE_TOL,
    interpolation_fit_threshold = INTERPOLATION_THRESHOLD,

    bootstrap_reps = 0L,
    include_matlab_gcv = FALSE,

    n_cores = N_CORES,
    overwrite = FALSE,
    verbose = TRUE
  )
}

# -----------------------------------------------------------------------------
# 4. Combine the dimension-specific runs
# -----------------------------------------------------------------------------

.combine_across_runs <- function(run_list, element) {
  pieces <- lapply(run_list, function(x) x[[element]])
  pieces <- Filter(Negate(is.null), pieces)

  if (!length(pieces)) {
    return(NULL)
  }

  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

result_elements <- c(
  "performance",
  "components",
  "coefficients",
  "tuning",
  "inference",
  "selected",
  "error"
)

main_results <- setNames(
  lapply(result_elements, function(element) {
    .combine_across_runs(results_by_p, element)
  }),
  result_elements
)

saveRDS(
  main_results,
  file.path(RESULTS_ROOT, "main_results_combined.rds")
)

for (nm in names(main_results)) {
  if (!is.null(main_results[[nm]])) {
    utils::write.csv(
      main_results[[nm]],
      file.path(RESULTS_ROOT, paste0("combined_", nm, ".csv")),
      row.names = FALSE
    )
  }
}

# -----------------------------------------------------------------------------
# 5. Add design information to result and summary tables
# -----------------------------------------------------------------------------

DESIGN_KEY <- MAIN_DESIGN[, c(
  "condition_id",
  "scenario",
  "n_train",
  "n_test",
  "n_sets",
  "p_per_set",
  "rho_within",
  "rho_between",
  "target_r2",
  "nonlinearity_strength",
  "seed_base"
), drop = FALSE]

.add_design_columns <- function(x) {
  if (is.null(x)) return(NULL)

  add_names <- setdiff(names(DESIGN_KEY), names(x))
  if (!length(add_names)) return(x)

  key <- DESIGN_KEY[, c("condition_id", add_names), drop = FALSE]
  out <- merge(x, key, by = "condition_id", all.x = TRUE, sort = FALSE)

  order_columns <- intersect(
    c("condition_id", "scenario", "n_train", "p_per_set", "target_r2",
      "replication", "method", "predictor_set"),
    names(out)
  )

  if (length(order_columns)) {
    ordering <- do.call(
  order,
  unname(out[order_columns])
)
    out <- out[ordering, , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}

main_summary <- summarize_kernel_era_simulation(main_results)

performance_summary <- .add_design_columns(main_summary$performance)
component_summary <- .add_design_columns(main_summary$components)
tuning_summary <- .add_design_columns(main_summary$tuning)
selected_with_design <- .add_design_columns(main_results$selected)

utils::write.csv(
  performance_summary,
  file.path(RESULTS_ROOT, "summary_performance.csv"),
  row.names = FALSE
)

utils::write.csv(
  component_summary,
  file.path(RESULTS_ROOT, "summary_components.csv"),
  row.names = FALSE
)

utils::write.csv(
  tuning_summary,
  file.path(RESULTS_ROOT, "summary_tuning.csv"),
  row.names = FALSE
)

utils::write.csv(
  selected_with_design,
  file.path(RESULTS_ROOT, "selected_tuning_with_design.csv"),
  row.names = FALSE
)

saveRDS(
  main_summary,
  file.path(RESULTS_ROOT, "main_summary_combined.rds")
)

# -----------------------------------------------------------------------------
# 6. Paired prediction comparison within each design condition
# -----------------------------------------------------------------------------

.mean_ci <- function(x, level = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)

  if (!n) {
    return(c(mean = NA_real_, sd = NA_real_, mcse = NA_real_,
             lower = NA_real_, upper = NA_real_))
  }

  m <- mean(x)
  s <- if (n > 1L) stats::sd(x) else NA_real_
  se <- if (n > 1L) s / sqrt(n) else NA_real_

  if (!is.finite(se)) {
    return(c(mean = m, sd = s, mcse = se,
             lower = NA_real_, upper = NA_real_))
  }

  crit <- stats::qt(1 - (1 - level) / 2, df = n - 1L)

  c(
    mean = m,
    sd = s,
    mcse = se,
    lower = m - crit * se,
    upper = m + crit * se
  )
}

perf <- main_results$performance

kernel_perf <- subset(
  perf,
  method == "Gaussian Kernel ERA",
  select = c(
    condition_id,
    scenario,
    replication,
    test_mse,
    test_r2,
    prediction_correlation,
    train_fit,
    train_test_fit_gap,
    converged,
    interpolation_flag
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
    prediction_correlation,
    train_fit,
    train_test_fit_gap,
    converged
  )
)

paired_prediction <- merge(
  kernel_perf,
  linear_perf,
  by = c("condition_id", "scenario", "replication"),
  suffixes = c("_kernel", "_linear")
)

paired_prediction$delta_r2 <-
  paired_prediction$test_r2_kernel - paired_prediction$test_r2_linear

paired_prediction$delta_mse <-
  paired_prediction$test_mse_kernel - paired_prediction$test_mse_linear

paired_prediction$delta_prediction_correlation <-
  paired_prediction$prediction_correlation_kernel -
  paired_prediction$prediction_correlation_linear

paired_prediction$kernel_r2_win <- paired_prediction$delta_r2 > 0
paired_prediction$kernel_mse_win <- paired_prediction$delta_mse < 0

paired_prediction$mse_reduction_proportion <-
  (paired_prediction$test_mse_linear - paired_prediction$test_mse_kernel) /
  paired_prediction$test_mse_linear

paired_prediction <- .add_design_columns(paired_prediction)

paired_prediction_summary <- do.call(
  rbind,
  lapply(split(paired_prediction, paired_prediction$condition_id), function(d) {
    r2_stats <- .mean_ci(d$delta_r2)
    mse_stats <- .mean_ci(d$delta_mse)
    cor_stats <- .mean_ci(d$delta_prediction_correlation)

    data.frame(
      condition_id = d$condition_id[1L],
      scenario = d$scenario[1L],
      n_rep = nrow(d),

      mean_delta_r2 = r2_stats["mean"],
      sd_delta_r2 = r2_stats["sd"],
      mcse_delta_r2 = r2_stats["mcse"],
      delta_r2_ci_lower = r2_stats["lower"],
      delta_r2_ci_upper = r2_stats["upper"],

      mean_delta_mse = mse_stats["mean"],
      mcse_delta_mse = mse_stats["mcse"],
      delta_mse_ci_lower = mse_stats["lower"],
      delta_mse_ci_upper = mse_stats["upper"],

      mean_delta_prediction_correlation = cor_stats["mean"],
      mcse_delta_prediction_correlation = cor_stats["mcse"],

      kernel_r2_win_rate = mean(d$kernel_r2_win),
      kernel_mse_win_rate = mean(d$kernel_mse_win),
      mean_mse_reduction = mean(d$mse_reduction_proportion),

      mean_kernel_train_test_gap = mean(d$train_test_fit_gap_kernel),
      mean_linear_train_test_gap = mean(d$train_test_fit_gap_linear),
      stringsAsFactors = FALSE
    )
  })
)

paired_prediction_summary <- .add_design_columns(paired_prediction_summary)

utils::write.csv(
  paired_prediction,
  file.path(RESULTS_ROOT, "paired_prediction_by_replication.csv"),
  row.names = FALSE
)

utils::write.csv(
  paired_prediction_summary,
  file.path(RESULTS_ROOT, "paired_prediction_summary_by_condition.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 7. Paired component-recovery comparison within each design condition
# -----------------------------------------------------------------------------

comp <- main_results$components

kernel_comp <- subset(
  comp,
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

linear_comp <- subset(
  comp,
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

paired_components <- merge(
  kernel_comp,
  linear_comp,
  by = c("condition_id", "scenario", "replication", "predictor_set"),
  suffixes = c("_kernel", "_linear")
)

paired_components$delta_component_correlation <-
  paired_components$test_component_correlation_kernel -
  paired_components$test_component_correlation_linear

paired_components$delta_component_rmse <-
  paired_components$test_component_rmse_calibrated_kernel -
  paired_components$test_component_rmse_calibrated_linear

paired_components$kernel_correlation_win <-
  paired_components$delta_component_correlation > 0

paired_components$kernel_rmse_win <-
  paired_components$delta_component_rmse < 0

paired_components <- .add_design_columns(paired_components)

component_split <- interaction(
  paired_components$condition_id,
  paired_components$predictor_set,
  drop = TRUE
)

paired_component_summary <- do.call(
  rbind,
  lapply(split(paired_components, component_split), function(d) {
    cor_stats <- .mean_ci(d$delta_component_correlation)
    rmse_stats <- .mean_ci(d$delta_component_rmse)

    data.frame(
      condition_id = d$condition_id[1L],
      scenario = d$scenario[1L],
      predictor_set = d$predictor_set[1L],
      n_rep = nrow(d),

      mean_delta_component_correlation = cor_stats["mean"],
      mcse_delta_component_correlation = cor_stats["mcse"],
      component_correlation_ci_lower = cor_stats["lower"],
      component_correlation_ci_upper = cor_stats["upper"],
      kernel_component_correlation_win_rate =
        mean(d$kernel_correlation_win),

      mean_delta_component_rmse = rmse_stats["mean"],
      mcse_delta_component_rmse = rmse_stats["mcse"],
      component_rmse_ci_lower = rmse_stats["lower"],
      component_rmse_ci_upper = rmse_stats["upper"],
      kernel_component_rmse_win_rate = mean(d$kernel_rmse_win),
      stringsAsFactors = FALSE
    )
  })
)

paired_component_summary <- .add_design_columns(paired_component_summary)

utils::write.csv(
  paired_components,
  file.path(RESULTS_ROOT, "paired_components_by_replication.csv"),
  row.names = FALSE
)

utils::write.csv(
  paired_component_summary,
  file.path(RESULTS_ROOT, "paired_component_summary_by_condition.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 8. Convergence, interpolation, and tuning-boundary diagnostics
# -----------------------------------------------------------------------------

selected <- selected_with_design

selected$any_boundary <-
  selected$lambda_at_lower_boundary |
  selected$lambda_at_upper_boundary |
  selected$sigma_at_lower_boundary |
  selected$sigma_at_upper_boundary

condition_diagnostics <- do.call(
  rbind,
  lapply(split(selected, selected$condition_id), function(d) {
    data.frame(
      condition_id = d$condition_id[1L],
      scenario = d$scenario[1L],
      n_rep = nrow(d),

      selected_eligibility_rate = mean(d$selected_cv_eligible),
      selected_cv_full_convergence_rate =
        mean(d$selected_cv_convergence_rate == 1),
      unrestricted_minimum_exclusion_rate =
        mean(!d$unrestricted_min_eligible),

      final_convergence_rate = mean(d$final_converged),
      final_interpolation_rate = mean(d$final_interpolation_flag),

      lambda_lower_boundary_rate = mean(d$lambda_at_lower_boundary),
      lambda_upper_boundary_rate = mean(d$lambda_at_upper_boundary),
      sigma_lower_boundary_rate = mean(d$sigma_at_lower_boundary),
      sigma_upper_boundary_rate = mean(d$sigma_at_upper_boundary),
      any_boundary_rate = mean(d$any_boundary),

      median_lambda = stats::median(d$lambda),
      median_sigma = stats::median(d$sigma),
      median_cv_mse = stats::median(d$selected_cv_mse),
      mean_elapsed_seconds = mean(d$elapsed_seconds),
      stringsAsFactors = FALSE
    )
  })
)

condition_diagnostics <- .add_design_columns(condition_diagnostics)

stability_checks <- merge(
  condition_diagnostics,
  paired_prediction_summary[, c(
    "condition_id",
    "mean_delta_r2",
    "mcse_delta_r2",
    "delta_r2_ci_lower",
    "delta_r2_ci_upper",
    "kernel_r2_win_rate",
    "mean_mse_reduction"
  )],
  by = "condition_id",
  all.x = TRUE,
  sort = FALSE
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

stability_checks <- stability_checks[
  order(stability_checks$condition_id),
  ,
  drop = FALSE
]

utils::write.csv(
  condition_diagnostics,
  file.path(RESULTS_ROOT, "condition_diagnostics.csv"),
  row.names = FALSE
)

utils::write.csv(
  stability_checks,
  file.path(RESULTS_ROOT, "stability_checks.csv"),
  row.names = FALSE
)

boundary_cases <- selected[selected$any_boundary, , drop = FALSE]

utils::write.csv(
  boundary_cases,
  file.path(RESULTS_ROOT, "boundary_cases.csv"),
  row.names = FALSE
)

nonconverged_final <- subset(selected, !final_converged)
interpolating_final <- subset(selected, final_interpolation_flag)
excluded_unrestricted_minima <- subset(selected, !unrestricted_min_eligible)

utils::write.csv(
  nonconverged_final,
  file.path(RESULTS_ROOT, "nonconverged_final_fits.csv"),
  row.names = FALSE
)

utils::write.csv(
  interpolating_final,
  file.path(RESULTS_ROOT, "interpolating_final_fits.csv"),
  row.names = FALSE
)

utils::write.csv(
  excluded_unrestricted_minima,
  file.path(RESULTS_ROOT, "excluded_unrestricted_cv_minima.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 9. Console summary
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("1. PERFORMANCE SUMMARY\n")
cat("============================================================\n")
print(performance_summary, row.names = FALSE)

cat("\n============================================================\n")
cat("2. COMPONENT RECOVERY SUMMARY\n")
cat("============================================================\n")
print(component_summary, row.names = FALSE)

cat("\n============================================================\n")
cat("3. TUNING AND CONVERGENCE SUMMARY\n")
cat("============================================================\n")
print(tuning_summary, row.names = FALSE)

cat("\n============================================================\n")
cat("4. PAIRED PREDICTION SUMMARY BY CONDITION\n")
cat("============================================================\n")
print(paired_prediction_summary, row.names = FALSE)

cat("\n============================================================\n")
cat("5. PAIRED COMPONENT SUMMARY BY CONDITION AND SET\n")
cat("============================================================\n")
print(paired_component_summary, row.names = FALSE)

cat("\n============================================================\n")
cat("6. STABILITY CHECKS\n")
cat("============================================================\n")
print(stability_checks, row.names = FALSE)

cat("\n============================================================\n")
cat("7. FINAL-FIT NONCONVERGENCE\n")
cat("============================================================\n")
if (!nrow(nonconverged_final)) {
  cat("All final Kernel ERA fits converged.\n")
} else {
  print(nonconverged_final, row.names = FALSE)
}

cat("\n============================================================\n")
cat("8. INTERPOLATION CHECK\n")
cat("============================================================\n")
if (!nrow(interpolating_final)) {
  cat("No final Kernel ERA fit exceeded the interpolation threshold.\n")
} else {
  print(interpolating_final, row.names = FALSE)
}

cat("\n============================================================\n")
cat("9. TASK ERRORS\n")
cat("============================================================\n")
if (is.null(main_results$error) || !nrow(main_results$error)) {
  cat("No task-level errors were recorded.\n")
} else {
  print(main_results$error, row.names = FALSE)
}

cat("\n============================================================\n")
cat("MAIN SIMULATION COMPLETE\n")
cat("============================================================\n")
cat("Results were saved to:\n", RESULTS_ROOT, "\n", sep = "")
