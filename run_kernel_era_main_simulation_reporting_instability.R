# =============================================================================
# Kernel ERA main simulation with explicit finite-sample stability reporting
#
# Primary analysis specification
#   * 5-fold CV
#   * minimum CV error among candidates that converged in every fold
#   * Gaussian kernel exp(-||x-z||^2 / (2 sigma^2))
#   * RKHS-norm ridge penalty
#   * three final-model starts; best converged solution retained when available
#   * dimension-adjusted sigma grids for P_k = 5 and P_k = 20
#   * resumable condition-replication files
#
# Important reporting principle
#   * Final nonconvergence and near-interpolation are retained as simulation
#     outcomes. They are not automatically replaced by another tuning solution.
#   * Performance is summarized both:
#       (a) for all operational outputs, and
#       (b) conditionally for replications in which the final Kernel ERA fit
#           converged and did not exceed the interpolation threshold.
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

# Use "pilot" for a short end-to-end run of all 24 conditions.
# Use "main" for the final 200-replication simulation.
# RUN_MODE <- "pilot"
RUN_MODE <- "main"

N_REP_PILOT <- 10L
N_REP_MAIN <- 200L

N_REP <- switch(
  RUN_MODE,
  pilot = N_REP_PILOT,
  main = N_REP_MAIN,
  stop("RUN_MODE must be either 'pilot' or 'main'.", call. = FALSE)
)

# On Windows, parallel execution requires the future and future.apply packages.
# Start conservatively. Increase only when adequate RAM is available.
N_CORES <- 1L

CV_FOLDS <- 5L
FINAL_N_STARTS <- 3L
MAX_ITER <- 500L
TOL <- 1e-6
OBJECTIVE_TOL <- 1e-9
INTERPOLATION_THRESHOLD <- 0.98

# The lower bound 1e-5 repeatedly generated a clearly pathological
# small-lambda/small-sigma corner in preliminary checks. The 1e-4 lower bound
# is retained because it was selected in stable, fully converged solutions.
# The upper bound remains 1e3 so that strong regularization can be selected and
# reported rather than ruled out in difficult conditions.
LAMBDA_GRID <- 10^seq(-4, 3, by = 1)

# Keep the small-sigma region in the design. Instability arising there under
# small-N, larger-P, low-signal conditions is a substantive finite-sample result,
# not something to remove after observing the outcome.
#
# Sigma is scaled by sqrt(P_k / 5), because Euclidean distances among
# standardized observations increase approximately in proportion to sqrt(P_k).
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
    "results_main_final_report_instability_dimension_scaled"
  } else {
    "results_main_pilot_report_instability_dimension_scaled"
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
    "stable_replication_definition",
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
    "final_converged AND NOT final_interpolation_flag",
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
# 3. Run separately for each predictor-set dimension
# -----------------------------------------------------------------------------

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
# 4. Combine dimension-specific runs
# -----------------------------------------------------------------------------

.combine_across_runs <- function(run_list, element) {
  pieces <- lapply(run_list, function(x) x[[element]])
  pieces <- Filter(Negate(is.null), pieces)

  if (!length(pieces)) return(NULL)

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
# 5. Shared post-processing helpers
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

  out <- x
  add_names <- setdiff(names(DESIGN_KEY), names(out))

  if (length(add_names)) {
    key <- DESIGN_KEY[, c("condition_id", add_names), drop = FALSE]
    out <- merge(out, key, by = "condition_id", all.x = TRUE, sort = FALSE)
  }

  order_columns <- intersect(
    c(
      "condition_id", "scenario", "n_train", "p_per_set", "target_r2",
      "replication", "method", "predictor_set"
    ),
    names(out)
  )

  if (length(order_columns)) {
    order_args <- unname(as.list(out[order_columns]))
    ordering <- do.call(base::order, order_args)
    out <- out[ordering, , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}

.mean_finite <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

.rate_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

.mean_ci <- function(x, level = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)

  if (!n) {
    return(c(
      mean = NA_real_, sd = NA_real_, mcse = NA_real_,
      lower = NA_real_, upper = NA_real_
    ))
  }

  m <- mean(x)
  s <- if (n > 1L) stats::sd(x) else NA_real_
  se <- if (n > 1L) s / sqrt(n) else NA_real_

  if (!is.finite(se)) {
    return(c(mean = m, sd = s, mcse = se, lower = NA_real_, upper = NA_real_))
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

.filter_by_replication_keys <- function(x, keys) {
  if (is.null(x)) return(NULL)
  if (!nrow(keys)) return(x[0, , drop = FALSE])

  out <- merge(
    x,
    keys,
    by = c("condition_id", "replication"),
    all = FALSE,
    sort = FALSE
  )
  rownames(out) <- NULL
  out
}

# -----------------------------------------------------------------------------
# 6. Define final-fit status and summarize all versus stable replications
# -----------------------------------------------------------------------------

selected_with_design <- .add_design_columns(main_results$selected)

selected_with_design$final_nonconvergence_flag <-
  !selected_with_design$final_converged

selected_with_design$final_stable <-
  selected_with_design$final_converged &
  !selected_with_design$final_interpolation_flag

selected_with_design$final_operational_failure <-
  selected_with_design$final_nonconvergence_flag |
  selected_with_design$final_interpolation_flag

selected_with_design$final_status <- ifelse(
  selected_with_design$final_nonconvergence_flag &
    selected_with_design$final_interpolation_flag,
  "nonconverged_and_interpolating",
  ifelse(
    selected_with_design$final_nonconvergence_flag,
    "nonconverged",
    ifelse(
      selected_with_design$final_interpolation_flag,
      "interpolating",
      "stable"
    )
  )
)

utils::write.csv(
  selected_with_design,
  file.path(RESULTS_ROOT, "selected_tuning_with_design_and_status.csv"),
  row.names = FALSE
)

stable_keys <- selected_with_design[
  selected_with_design$final_stable,
  c("condition_id", "replication"),
  drop = FALSE
]

main_summary_all <- summarize_kernel_era_simulation(main_results)

stable_results <- main_results
stable_results$performance <- .filter_by_replication_keys(
  main_results$performance,
  stable_keys
)
stable_results$components <- .filter_by_replication_keys(
  main_results$components,
  stable_keys
)
stable_results$coefficients <- .filter_by_replication_keys(
  main_results$coefficients,
  stable_keys
)
stable_results$selected <- .filter_by_replication_keys(
  main_results$selected,
  stable_keys
)

main_summary_stable <- summarize_kernel_era_simulation(stable_results)

performance_summary_all <- .add_design_columns(main_summary_all$performance)
performance_summary_all$analysis_sample <- "all_operational_outputs"

performance_summary_stable <- .add_design_columns(main_summary_stable$performance)
performance_summary_stable$analysis_sample <-
  "kernel_converged_and_noninterpolating_replications"

performance_summary <- rbind(
  performance_summary_all,
  performance_summary_stable
)

component_summary_all <- .add_design_columns(main_summary_all$components)
component_summary_all$analysis_sample <- "all_operational_outputs"

component_summary_stable <- .add_design_columns(main_summary_stable$components)
component_summary_stable$analysis_sample <-
  "kernel_converged_and_noninterpolating_replications"

component_summary <- rbind(
  component_summary_all,
  component_summary_stable
)

tuning_summary <- .add_design_columns(main_summary_all$tuning)

utils::write.csv(
  performance_summary_all,
  file.path(RESULTS_ROOT, "summary_performance_all_operational.csv"),
  row.names = FALSE
)
utils::write.csv(
  performance_summary_stable,
  file.path(RESULTS_ROOT, "summary_performance_kernel_stable.csv"),
  row.names = FALSE
)
utils::write.csv(
  performance_summary,
  file.path(RESULTS_ROOT, "summary_performance_combined.csv"),
  row.names = FALSE
)

utils::write.csv(
  component_summary_all,
  file.path(RESULTS_ROOT, "summary_components_all_operational.csv"),
  row.names = FALSE
)
utils::write.csv(
  component_summary_stable,
  file.path(RESULTS_ROOT, "summary_components_kernel_stable.csv"),
  row.names = FALSE
)
utils::write.csv(
  component_summary,
  file.path(RESULTS_ROOT, "summary_components_combined.csv"),
  row.names = FALSE
)

utils::write.csv(
  tuning_summary,
  file.path(RESULTS_ROOT, "summary_tuning.csv"),
  row.names = FALSE
)

saveRDS(
  list(all = main_summary_all, stable = main_summary_stable),
  file.path(RESULTS_ROOT, "main_summaries_all_and_stable.rds")
)

# -----------------------------------------------------------------------------
# 7. Paired prediction comparison: all operational and stable replications
# -----------------------------------------------------------------------------

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

paired_prediction$kernel_final_stable <-
  paired_prediction$converged_kernel &
  !paired_prediction$interpolation_flag

paired_prediction$kernel_operational_failure <-
  !paired_prediction$converged_kernel |
  paired_prediction$interpolation_flag

paired_prediction <- .add_design_columns(paired_prediction)

.summarize_paired_prediction <- function(data, analysis_sample) {
  if (!nrow(data)) return(NULL)

  out <- do.call(
    rbind,
    lapply(split(data, data$condition_id), function(d) {
      r2_stats <- .mean_ci(d$delta_r2)
      mse_stats <- .mean_ci(d$delta_mse)
      cor_stats <- .mean_ci(d$delta_prediction_correlation)

      data.frame(
        condition_id = d$condition_id[1L],
        scenario = d$scenario[1L],
        analysis_sample = analysis_sample,
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

        kernel_r2_win_rate = .rate_or_na(d$kernel_r2_win),
        kernel_mse_win_rate = .rate_or_na(d$kernel_mse_win),
        mean_mse_reduction = .mean_finite(d$mse_reduction_proportion),

        mean_kernel_train_test_gap =
          .mean_finite(d$train_test_fit_gap_kernel),
        mean_linear_train_test_gap =
          .mean_finite(d$train_test_fit_gap_linear),
        stringsAsFactors = FALSE
      )
    })
  )

  .add_design_columns(out)
}

paired_prediction_summary_all <- .summarize_paired_prediction(
  paired_prediction,
  "all_operational_outputs"
)

paired_prediction_summary_stable <- .summarize_paired_prediction(
  paired_prediction[paired_prediction$kernel_final_stable, , drop = FALSE],
  "kernel_converged_and_noninterpolating_replications"
)

paired_prediction_summary <- rbind(
  paired_prediction_summary_all,
  paired_prediction_summary_stable
)

utils::write.csv(
  paired_prediction,
  file.path(RESULTS_ROOT, "paired_prediction_by_replication.csv"),
  row.names = FALSE
)
utils::write.csv(
  paired_prediction_summary_all,
  file.path(RESULTS_ROOT, "paired_prediction_summary_all_operational.csv"),
  row.names = FALSE
)
utils::write.csv(
  paired_prediction_summary_stable,
  file.path(RESULTS_ROOT, "paired_prediction_summary_kernel_stable.csv"),
  row.names = FALSE
)
utils::write.csv(
  paired_prediction_summary,
  file.path(RESULTS_ROOT, "paired_prediction_summary_combined.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 8. Paired component recovery: all operational and stable replications
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

status_key <- selected_with_design[, c(
  "condition_id",
  "replication",
  "final_converged",
  "final_interpolation_flag",
  "final_stable",
  "final_status"
), drop = FALSE]

paired_components <- merge(
  paired_components,
  status_key,
  by = c("condition_id", "replication"),
  all.x = TRUE,
  sort = FALSE
)
paired_components <- .add_design_columns(paired_components)

.summarize_paired_components <- function(data, analysis_sample) {
  if (!nrow(data)) return(NULL)

  component_split <- interaction(
    data$condition_id,
    data$predictor_set,
    drop = TRUE
  )

  out <- do.call(
    rbind,
    lapply(split(data, component_split), function(d) {
      cor_stats <- .mean_ci(d$delta_component_correlation)
      rmse_stats <- .mean_ci(d$delta_component_rmse)

      data.frame(
        condition_id = d$condition_id[1L],
        scenario = d$scenario[1L],
        predictor_set = d$predictor_set[1L],
        analysis_sample = analysis_sample,
        n_rep = nrow(d),

        mean_delta_component_correlation = cor_stats["mean"],
        mcse_delta_component_correlation = cor_stats["mcse"],
        component_correlation_ci_lower = cor_stats["lower"],
        component_correlation_ci_upper = cor_stats["upper"],
        kernel_component_correlation_win_rate =
          .rate_or_na(d$kernel_correlation_win),

        mean_delta_component_rmse = rmse_stats["mean"],
        mcse_delta_component_rmse = rmse_stats["mcse"],
        component_rmse_ci_lower = rmse_stats["lower"],
        component_rmse_ci_upper = rmse_stats["upper"],
        kernel_component_rmse_win_rate =
          .rate_or_na(d$kernel_rmse_win),
        stringsAsFactors = FALSE
      )
    })
  )

  .add_design_columns(out)
}

paired_component_summary_all <- .summarize_paired_components(
  paired_components,
  "all_operational_outputs"
)

paired_component_summary_stable <- .summarize_paired_components(
  paired_components[paired_components$final_stable, , drop = FALSE],
  "kernel_converged_and_noninterpolating_replications"
)

paired_component_summary <- rbind(
  paired_component_summary_all,
  paired_component_summary_stable
)

utils::write.csv(
  paired_components,
  file.path(RESULTS_ROOT, "paired_components_by_replication.csv"),
  row.names = FALSE
)
utils::write.csv(
  paired_component_summary_all,
  file.path(RESULTS_ROOT, "paired_component_summary_all_operational.csv"),
  row.names = FALSE
)
utils::write.csv(
  paired_component_summary_stable,
  file.path(RESULTS_ROOT, "paired_component_summary_kernel_stable.csv"),
  row.names = FALSE
)
utils::write.csv(
  paired_component_summary,
  file.path(RESULTS_ROOT, "paired_component_summary_combined.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 9. Operational stability, convergence, interpolation, and boundaries
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
      n_rep_requested = N_REP,
      n_rep_observed = nrow(d),

      n_final_converged = sum(d$final_converged, na.rm = TRUE),
      n_final_nonconverged = sum(!d$final_converged, na.rm = TRUE),
      n_final_interpolating = sum(d$final_interpolation_flag, na.rm = TRUE),
      n_final_stable = sum(d$final_stable, na.rm = TRUE),
      n_operational_failures = sum(d$final_operational_failure, na.rm = TRUE),

      selected_eligibility_rate = .rate_or_na(d$selected_cv_eligible),
      selected_cv_full_convergence_rate =
        .rate_or_na(d$selected_cv_convergence_rate == 1),
      unrestricted_minimum_exclusion_rate =
        .rate_or_na(!d$unrestricted_min_eligible),

      final_convergence_rate = .rate_or_na(d$final_converged),
      final_nonconvergence_rate = .rate_or_na(!d$final_converged),
      final_interpolation_rate = .rate_or_na(d$final_interpolation_flag),
      final_stable_rate = .rate_or_na(d$final_stable),
      final_operational_failure_rate =
        .rate_or_na(d$final_operational_failure),

      lambda_lower_boundary_rate = .rate_or_na(d$lambda_at_lower_boundary),
      lambda_upper_boundary_rate = .rate_or_na(d$lambda_at_upper_boundary),
      sigma_lower_boundary_rate = .rate_or_na(d$sigma_at_lower_boundary),
      sigma_upper_boundary_rate = .rate_or_na(d$sigma_at_upper_boundary),
      any_boundary_rate = .rate_or_na(d$any_boundary),

      median_lambda = stats::median(d$lambda, na.rm = TRUE),
      median_sigma = stats::median(d$sigma, na.rm = TRUE),
      median_cv_mse = stats::median(d$selected_cv_mse, na.rm = TRUE),
      mean_elapsed_seconds = .mean_finite(d$elapsed_seconds),
      stringsAsFactors = FALSE
    )
  })
)

condition_diagnostics <- .add_design_columns(condition_diagnostics)

# These are reporting flags, not pass/fail rules. Elevated rates in difficult
# conditions are retained as findings rather than used to modify the estimator.
condition_diagnostics$report_nonconvergence <-
  condition_diagnostics$final_nonconvergence_rate > 0
condition_diagnostics$report_interpolation <-
  condition_diagnostics$final_interpolation_rate > 0
condition_diagnostics$report_boundary_concentration <-
  pmax(
    condition_diagnostics$lambda_lower_boundary_rate,
    condition_diagnostics$lambda_upper_boundary_rate,
    condition_diagnostics$sigma_lower_boundary_rate,
    condition_diagnostics$sigma_upper_boundary_rate,
    na.rm = TRUE
  ) >= 0.20

operational_summary <- merge(
  condition_diagnostics,
  paired_prediction_summary_all[, c(
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

stable_prediction_columns <- paired_prediction_summary_stable[, c(
  "condition_id",
  "n_rep",
  "mean_delta_r2",
  "mcse_delta_r2",
  "kernel_r2_win_rate",
  "mean_mse_reduction"
), drop = FALSE]

names(stable_prediction_columns)[-1L] <- paste0(
  names(stable_prediction_columns)[-1L],
  "_stable_only"
)

operational_summary <- merge(
  operational_summary,
  stable_prediction_columns,
  by = "condition_id",
  all.x = TRUE,
  sort = FALSE
)

operational_summary <- operational_summary[
  order(operational_summary$condition_id),
  ,
  drop = FALSE
]

utils::write.csv(
  condition_diagnostics,
  file.path(RESULTS_ROOT, "condition_operational_stability.csv"),
  row.names = FALSE
)
utils::write.csv(
  operational_summary,
  file.path(RESULTS_ROOT, "condition_operational_and_predictive_summary.csv"),
  row.names = FALSE
)

boundary_cases <- selected[selected$any_boundary, , drop = FALSE]
nonconverged_final <- selected[selected$final_nonconvergence_flag, , drop = FALSE]
interpolating_final <- selected[selected$final_interpolation_flag, , drop = FALSE]
operational_failures <- selected[
  selected$final_operational_failure,
  ,
  drop = FALSE
]
excluded_unrestricted_minima <- selected[
  !selected$unrestricted_min_eligible,
  ,
  drop = FALSE
]

utils::write.csv(
  boundary_cases,
  file.path(RESULTS_ROOT, "boundary_cases.csv"),
  row.names = FALSE
)
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
  operational_failures,
  file.path(RESULTS_ROOT, "final_operational_failures.csv"),
  row.names = FALSE
)
utils::write.csv(
  excluded_unrestricted_minima,
  file.path(RESULTS_ROOT, "excluded_unrestricted_cv_minima.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 10. Pipeline-integrity check
# -----------------------------------------------------------------------------

EXPECTED_TASKS <- nrow(MAIN_DESIGN) * N_REP
N_TASK_ERRORS <- if (is.null(main_results$error)) 0L else nrow(main_results$error)
N_SELECTED_ROWS <- if (is.null(main_results$selected)) 0L else nrow(main_results$selected)

required_objects_created <- all(vapply(
  list(
    performance_summary_all,
    performance_summary_stable,
    component_summary_all,
    component_summary_stable,
    tuning_summary,
    paired_prediction_summary_all,
    paired_prediction_summary_stable,
    paired_component_summary_all,
    paired_component_summary_stable,
    condition_diagnostics
  ),
  function(x) !is.null(x),
  logical(1L)
))

pipeline_integrity <- data.frame(
  run_mode = RUN_MODE,
  expected_tasks = EXPECTED_TASKS,
  selected_result_rows = N_SELECTED_ROWS,
  task_error_rows = N_TASK_ERRORS,
  all_expected_selected_rows_present = N_SELECTED_ROWS == EXPECTED_TASKS,
  no_task_level_errors = N_TASK_ERRORS == 0L,
  all_selected_cv_candidates_eligible =
    all(selected$selected_cv_eligible, na.rm = TRUE),
  all_summary_objects_created = required_objects_created,
  pipeline_complete =
    N_SELECTED_ROWS == EXPECTED_TASKS &&
    N_TASK_ERRORS == 0L &&
    required_objects_created,
  stringsAsFactors = FALSE
)

utils::write.csv(
  pipeline_integrity,
  file.path(RESULTS_ROOT, "pipeline_integrity.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 11. Console summary
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("1. PERFORMANCE SUMMARY: ALL OPERATIONAL OUTPUTS\n")
cat("============================================================\n")
print(performance_summary_all, row.names = FALSE)

cat("\n============================================================\n")
cat("2. PERFORMANCE SUMMARY: KERNEL-STABLE REPLICATIONS ONLY\n")
cat("============================================================\n")
print(performance_summary_stable, row.names = FALSE)

cat("\n============================================================\n")
cat("3. COMPONENT SUMMARY: ALL OPERATIONAL OUTPUTS\n")
cat("============================================================\n")
print(component_summary_all, row.names = FALSE)

cat("\n============================================================\n")
cat("4. COMPONENT SUMMARY: KERNEL-STABLE REPLICATIONS ONLY\n")
cat("============================================================\n")
print(component_summary_stable, row.names = FALSE)

cat("\n============================================================\n")
cat("5. TUNING SUMMARY\n")
cat("============================================================\n")
print(tuning_summary, row.names = FALSE)

cat("\n============================================================\n")
cat("6. OPERATIONAL STABILITY BY CONDITION\n")
cat("============================================================\n")
print(condition_diagnostics, row.names = FALSE)

cat("\n============================================================\n")
cat("7. PAIRED PREDICTION: ALL OPERATIONAL OUTPUTS\n")
cat("============================================================\n")
print(paired_prediction_summary_all, row.names = FALSE)

cat("\n============================================================\n")
cat("8. PAIRED PREDICTION: KERNEL-STABLE REPLICATIONS ONLY\n")
cat("============================================================\n")
print(paired_prediction_summary_stable, row.names = FALSE)

cat("\n============================================================\n")
cat("9. FINAL NONCONVERGENCE\n")
cat("============================================================\n")
if (!nrow(nonconverged_final)) {
  cat("No final Kernel ERA nonconvergence was observed.\n")
} else {
  print(nonconverged_final, row.names = FALSE)
}

cat("\n============================================================\n")
cat("10. FINAL INTERPOLATION\n")
cat("============================================================\n")
if (!nrow(interpolating_final)) {
  cat("No final Kernel ERA fit exceeded the interpolation threshold.\n")
} else {
  print(interpolating_final, row.names = FALSE)
}

cat("\n============================================================\n")
cat("11. TASK ERRORS\n")
cat("============================================================\n")
if (is.null(main_results$error) || !nrow(main_results$error)) {
  cat("No task-level errors were recorded.\n")
} else {
  print(main_results$error, row.names = FALSE)
}

cat("\n============================================================\n")
cat("12. PIPELINE INTEGRITY\n")
cat("============================================================\n")
print(pipeline_integrity, row.names = FALSE)

if (isTRUE(pipeline_integrity$pipeline_complete)) {
  cat(
    "\nThe computational pipeline completed successfully.\n",
    "Observed nonconvergence and interpolation are retained as substantive\n",
    "simulation outcomes and do not, by themselves, indicate pipeline failure.\n",
    sep = ""
  )

  if (RUN_MODE == "pilot") {
    cat(
      "After reviewing the condition-specific operational rates, change\n",
      "RUN_MODE <- \"main\" to start the 200-replication simulation.\n",
      sep = ""
    )
  }
} else {
  cat(
    "\nThe computational pipeline did not complete cleanly. Review task errors\n",
    "or missing result rows before starting the main simulation.\n",
    sep = ""
  )
}

cat("\n============================================================\n")
cat("SIMULATION RUN COMPLETE\n")
cat("============================================================\n")
cat("Results were saved to:\n", RESULTS_ROOT, "\n", sep = "")
