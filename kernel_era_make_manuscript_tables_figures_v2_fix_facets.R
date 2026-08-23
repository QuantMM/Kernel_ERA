# =============================================================================
# Kernel ERA manuscript tables and figures
#
# Inputs expected in DATA_DIR:
#   pipeline_integrity.csv
#   condition_operational_stability.csv
#   paired_prediction_summary_combined.csv
#   paired_component_summary_combined.csv
#   summary_tuning.csv
#
# Outputs written to OUT_DIR:
#   manuscript_master_results_by_condition.csv
#   manuscript_component_master_by_condition_set.csv
#   manuscript_aggregate_prediction_by_scenario.csv
#   manuscript_aggregate_component_by_scenario_set.csv
#   manuscript_operational_stability_by_design.csv
#   manuscript_prediction_by_design_cell.csv
#   figures/*.png and figures/*.pdf
# =============================================================================


install.packages(c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "forcats",
  "scales"
))

# -----------------------------------------------------------------------------
# 0. Paths and packages
# -----------------------------------------------------------------------------

DATA_DIR <- "C:/Users/kims15/Desktop/Kernel ERA Claude code/results_main_final_report_instability_dimension_scaled"
# Change this to your results-summary folder if needed.
OUT_DIR  <- "C:/Users/kims15/Desktop/Kernel ERA Claude code/kernel_era_manuscript_outputs_500rep"
FIG_DIR  <- file.path(OUT_DIR, "figures")

required_packages <- c("readr", "dplyr", "tidyr", "ggplot2", "stringr", "forcats", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Please install the following packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(forcats)
library(scales)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Read required files and verify computational integrity
# -----------------------------------------------------------------------------

pipeline_integrity <- read_csv(file.path(DATA_DIR, "pipeline_integrity.csv"), show_col_types = FALSE)
condition_stability <- read_csv(file.path(DATA_DIR, "condition_operational_stability.csv"), show_col_types = FALSE)
paired_prediction <- read_csv(file.path(DATA_DIR, "paired_prediction_summary_combined.csv"), show_col_types = FALSE)
paired_components <- read_csv(file.path(DATA_DIR, "paired_component_summary_combined.csv"), show_col_types = FALSE)
summary_tuning <- read_csv(file.path(DATA_DIR, "summary_tuning.csv"), show_col_types = FALSE)

if (!isTRUE(pipeline_integrity$pipeline_complete[1])) {
  warning("pipeline_complete is not TRUE. Inspect pipeline_integrity.csv before interpreting results.")
}

cat("Pipeline integrity:\n")
print(pipeline_integrity)

# -----------------------------------------------------------------------------
# 2. Helper variables and labels
# -----------------------------------------------------------------------------

scenario_levels <- c("linear", "mixed", "nonlinear")
scenario_labels <- c(
  linear = "Linear truth",
  mixed = "Mixed truth",
  nonlinear = "Nonlinear truth"
)

sample_labels <- c(
  all_operational_outputs = "all",
  kernel_converged_and_noninterpolating_replications = "stable"
)

label_pk <- function(x) paste0("P[k] == ", x)
label_r2 <- function(x) paste0("R^2 == ", x)
label_n <- function(x) paste0("N = ", x)

# -----------------------------------------------------------------------------
# 3. Step 1: condition-level master results table
# -----------------------------------------------------------------------------

prediction_wide <- paired_prediction %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    analysis_sample_short = recode(analysis_sample, !!!sample_labels)
  ) %>%
  select(
    condition_id,
    analysis_sample_short,
    n_rep,
    mean_delta_r2,
    sd_delta_r2,
    mcse_delta_r2,
    delta_r2_ci_lower,
    delta_r2_ci_upper,
    mean_delta_mse,
    mcse_delta_mse,
    delta_mse_ci_lower,
    delta_mse_ci_upper,
    mean_delta_prediction_correlation,
    mcse_delta_prediction_correlation,
    kernel_r2_win_rate,
    kernel_mse_win_rate,
    mean_mse_reduction,
    mean_kernel_train_test_gap,
    mean_linear_train_test_gap
  ) %>%
  pivot_wider(
    id_cols = condition_id,
    names_from = analysis_sample_short,
    values_from = c(
      n_rep,
      mean_delta_r2,
      sd_delta_r2,
      mcse_delta_r2,
      delta_r2_ci_lower,
      delta_r2_ci_upper,
      mean_delta_mse,
      mcse_delta_mse,
      delta_mse_ci_lower,
      delta_mse_ci_upper,
      mean_delta_prediction_correlation,
      mcse_delta_prediction_correlation,
      kernel_r2_win_rate,
      kernel_mse_win_rate,
      mean_mse_reduction,
      mean_kernel_train_test_gap,
      mean_linear_train_test_gap
    ),
    names_glue = "{.value}_{analysis_sample_short}"
  )

master_results_by_condition <- condition_stability %>%
  select(
    condition_id,
    scenario,
    n_train,
    p_per_set,
    target_r2,
    n_test,
    n_sets,
    rho_within,
    rho_between,
    nonlinearity_strength,
    seed_base,
    n_rep_requested,
    n_rep_observed,
    n_final_converged,
    n_final_nonconverged,
    n_final_interpolating,
    n_final_stable,
    n_operational_failures,
    selected_eligibility_rate,
    selected_cv_full_convergence_rate,
    unrestricted_minimum_exclusion_rate,
    final_convergence_rate,
    final_nonconvergence_rate,
    final_interpolation_rate,
    final_stable_rate,
    final_operational_failure_rate,
    lambda_lower_boundary_rate,
    lambda_upper_boundary_rate,
    sigma_lower_boundary_rate,
    sigma_upper_boundary_rate,
    any_boundary_rate,
    median_lambda,
    median_sigma,
    median_cv_mse,
    mean_elapsed_seconds,
    report_nonconvergence,
    report_interpolation,
    report_boundary_concentration
  ) %>%
  left_join(prediction_wide, by = "condition_id") %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    n_train_factor = factor(n_train, levels = c(150, 300)),
    p_per_set_factor = factor(p_per_set, levels = c(5, 20)),
    target_r2_factor = factor(target_r2, levels = c(0.30, 0.60))
  ) %>%
  arrange(target_r2, p_per_set, n_train, scenario)

write_csv(
  master_results_by_condition,
  file.path(OUT_DIR, "manuscript_master_results_by_condition.csv")
)

aggregate_prediction_by_scenario <- master_results_by_condition %>%
  group_by(scenario) %>%
  summarize(
    n_conditions = n(),
    mean_delta_r2_all = mean(mean_delta_r2_all, na.rm = TRUE),
    mean_delta_r2_stable = mean(mean_delta_r2_stable, na.rm = TRUE),
    mean_kernel_r2_win_rate_all = mean(kernel_r2_win_rate_all, na.rm = TRUE),
    mean_mse_reduction_all = mean(mean_mse_reduction_all, na.rm = TRUE),
    mean_final_stable_rate = mean(final_stable_rate, na.rm = TRUE),
    min_final_stable_rate = min(final_stable_rate, na.rm = TRUE),
    mean_final_nonconvergence_rate = mean(final_nonconvergence_rate, na.rm = TRUE),
    max_final_nonconvergence_rate = max(final_nonconvergence_rate, na.rm = TRUE),
    mean_final_interpolation_rate = mean(final_interpolation_rate, na.rm = TRUE),
    mean_unrestricted_minimum_exclusion_rate = mean(unrestricted_minimum_exclusion_rate, na.rm = TRUE),
    mean_any_boundary_rate = mean(any_boundary_rate, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  aggregate_prediction_by_scenario,
  file.path(OUT_DIR, "manuscript_aggregate_prediction_by_scenario.csv")
)

prediction_by_design_cell <- master_results_by_condition %>%
  group_by(scenario, n_train, p_per_set, target_r2) %>%
  summarize(
    n_conditions = n(),
    mean_delta_r2_all = mean(mean_delta_r2_all, na.rm = TRUE),
    mcse_delta_r2_all = mean(mcse_delta_r2_all, na.rm = TRUE),
    kernel_r2_win_rate_all = mean(kernel_r2_win_rate_all, na.rm = TRUE),
    mean_mse_reduction_all = mean(mean_mse_reduction_all, na.rm = TRUE),
    final_stable_rate = mean(final_stable_rate, na.rm = TRUE),
    final_nonconvergence_rate = mean(final_nonconvergence_rate, na.rm = TRUE),
    final_interpolation_rate = mean(final_interpolation_rate, na.rm = TRUE),
    unrestricted_minimum_exclusion_rate = mean(unrestricted_minimum_exclusion_rate, na.rm = TRUE),
    any_boundary_rate = mean(any_boundary_rate, na.rm = TRUE),
    median_lambda = median(median_lambda, na.rm = TRUE),
    median_sigma = median(median_sigma, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(target_r2, p_per_set, n_train, scenario)

write_csv(
  prediction_by_design_cell,
  file.path(OUT_DIR, "manuscript_prediction_by_design_cell.csv")
)

# -----------------------------------------------------------------------------
# 4. Step 2: component-level master table
# -----------------------------------------------------------------------------

component_wide <- paired_components %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    analysis_sample_short = recode(analysis_sample, !!!sample_labels)
  ) %>%
  select(
    condition_id,
    predictor_set,
    analysis_sample_short,
    n_rep,
    mean_delta_component_correlation,
    mcse_delta_component_correlation,
    component_correlation_ci_lower,
    component_correlation_ci_upper,
    kernel_component_correlation_win_rate,
    mean_delta_component_rmse,
    mcse_delta_component_rmse,
    component_rmse_ci_lower,
    component_rmse_ci_upper,
    kernel_component_rmse_win_rate
  ) %>%
  pivot_wider(
    id_cols = c(condition_id, predictor_set),
    names_from = analysis_sample_short,
    values_from = c(
      n_rep,
      mean_delta_component_correlation,
      mcse_delta_component_correlation,
      component_correlation_ci_lower,
      component_correlation_ci_upper,
      kernel_component_correlation_win_rate,
      mean_delta_component_rmse,
      mcse_delta_component_rmse,
      component_rmse_ci_lower,
      component_rmse_ci_upper,
      kernel_component_rmse_win_rate
    ),
    names_glue = "{.value}_{analysis_sample_short}"
  )

component_design <- paired_components %>%
  distinct(
    condition_id,
    predictor_set,
    scenario,
    n_train,
    n_test,
    n_sets,
    p_per_set,
    rho_within,
    rho_between,
    target_r2,
    nonlinearity_strength,
    seed_base
  )

component_master_by_condition_set <- component_design %>%
  left_join(component_wide, by = c("condition_id", "predictor_set")) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    predictor_set = factor(predictor_set, levels = c(1, 2, 3)),
    n_train_factor = factor(n_train, levels = c(150, 300)),
    p_per_set_factor = factor(p_per_set, levels = c(5, 20)),
    target_r2_factor = factor(target_r2, levels = c(0.30, 0.60))
  ) %>%
  arrange(target_r2, p_per_set, n_train, scenario, predictor_set)

write_csv(
  component_master_by_condition_set,
  file.path(OUT_DIR, "manuscript_component_master_by_condition_set.csv")
)

aggregate_component_by_scenario_set <- component_master_by_condition_set %>%
  group_by(scenario, predictor_set) %>%
  summarize(
    n_conditions = n(),
    mean_delta_component_correlation_all = mean(mean_delta_component_correlation_all, na.rm = TRUE),
    mean_kernel_component_correlation_win_rate_all = mean(kernel_component_correlation_win_rate_all, na.rm = TRUE),
    mean_delta_component_rmse_all = mean(mean_delta_component_rmse_all, na.rm = TRUE),
    mean_kernel_component_rmse_win_rate_all = mean(kernel_component_rmse_win_rate_all, na.rm = TRUE),
    mean_delta_component_correlation_stable = mean(mean_delta_component_correlation_stable, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(scenario, predictor_set)

write_csv(
  aggregate_component_by_scenario_set,
  file.path(OUT_DIR, "manuscript_aggregate_component_by_scenario_set.csv")
)

operational_stability_by_design <- master_results_by_condition %>%
  group_by(n_train, p_per_set, target_r2) %>%
  summarize(
    mean_stable_rate = mean(final_stable_rate, na.rm = TRUE),
    min_stable_rate = min(final_stable_rate, na.rm = TRUE),
    mean_nonconvergence_rate = mean(final_nonconvergence_rate, na.rm = TRUE),
    max_nonconvergence_rate = max(final_nonconvergence_rate, na.rm = TRUE),
    mean_interpolation_rate = mean(final_interpolation_rate, na.rm = TRUE),
    max_interpolation_rate = max(final_interpolation_rate, na.rm = TRUE),
    mean_cv_exclusion_rate = mean(unrestricted_minimum_exclusion_rate, na.rm = TRUE),
    max_cv_exclusion_rate = max(unrestricted_minimum_exclusion_rate, na.rm = TRUE),
    mean_any_boundary_rate = mean(any_boundary_rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(target_r2, p_per_set, n_train)

write_csv(
  operational_stability_by_design,
  file.path(OUT_DIR, "manuscript_operational_stability_by_design.csv")
)

# -----------------------------------------------------------------------------
# 5. Step 3: manuscript figure code
# -----------------------------------------------------------------------------

# Common theme ---------------------------------------------------------------
base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey92", color = "grey60"),
    plot.title.position = "plot"
  )

# Figure 1: prediction delta R2 by design condition --------------------------
fig1_data <- master_results_by_condition %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels, labels = scenario_labels),
    n_train_factor = factor(n_train, levels = c(150, 300), labels = c("N = 150", "N = 300")),
    p_panel = factor(p_per_set, levels = c(5, 20), labels = c("P[k] == 5", "P[k] == 20")),
    r2_panel = factor(target_r2, levels = c(0.30, 0.60), labels = c("R^2 == .30", "R^2 == .60"))
  )

figure1_prediction_delta_r2 <- ggplot(
  fig1_data,
  aes(x = scenario, y = mean_delta_r2_all, ymin = delta_r2_ci_lower_all, ymax = delta_r2_ci_upper_all,
      shape = n_train_factor, group = n_train_factor)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = 0.45), linewidth = 0.35) +
  facet_grid(r2_panel ~ p_panel, labeller = label_parsed) +
  labs(
    title = "Prediction advantage of Kernel ERA over Linear ERA",
    subtitle = expression(Delta*R^2 == R^2[Kernel] - R^2[Linear]),
    x = NULL,
    y = expression(Delta*R^2),
    shape = "Training sample size"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# Figure 1b: all-operational versus stable-only sensitivity ------------------
figure1b_all_vs_stable <- ggplot(
  master_results_by_condition,
  aes(x = mean_delta_r2_all, y = mean_delta_r2_stable)
) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4, linetype = "dashed") +
  geom_point(aes(shape = factor(p_per_set), size = final_stable_rate)) +
  facet_wrap(~ scenario, labeller = as_labeller(scenario_labels)) +
  labs(
    title = "All-operational and stable-only summaries give the same conclusions",
    x = expression(Delta*R^2~"(all operational outputs)"),
    y = expression(Delta*R^2~"(Kernel-stable replications only)"),
    shape = expression(P[k]),
    size = "Stable rate"
  ) +
  base_theme

# Figure 2: component recovery aggregate by scenario and predictor set --------
fig2_data <- aggregate_component_by_scenario_set %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels, labels = scenario_labels),
    predictor_set = factor(predictor_set, levels = c(1, 2, 3), labels = c("Set 1", "Set 2", "Set 3"))
  )

figure2_component_recovery_overview <- ggplot(
  fig2_data,
  aes(x = predictor_set, y = mean_delta_component_correlation_all)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
  geom_col(width = 0.65) +
  facet_wrap(~ scenario, nrow = 1) +
  labs(
    title = "Component recovery advantage by true functional form",
    subtitle = expression(Delta*cor == cor(hat(f)[Kernel], f[true]) - cor(hat(f)[Linear], f[true])),
    x = NULL,
    y = expression("Mean "*Delta*" component correlation")
  ) +
  base_theme

# Figure 2b: condition-level component recovery ------------------------------
fig2b_data <- component_master_by_condition_set %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels, labels = scenario_labels),
    predictor_set = factor(predictor_set, levels = c(1, 2, 3), labels = c("Set 1", "Set 2", "Set 3")),
    n_train_factor = factor(n_train, levels = c(150, 300), labels = c("N = 150", "N = 300")),
    p_panel = factor(p_per_set, levels = c(5, 20), labels = c("P[k] == 5", "P[k] == 20")),
    r2_panel = factor(target_r2, levels = c(0.30, 0.60), labels = c("R^2 == .30", "R^2 == .60"))
  )

figure2b_component_recovery_detailed <- ggplot(
  fig2b_data,
  aes(x = predictor_set, y = mean_delta_component_correlation_all,
      ymin = component_correlation_ci_lower_all, ymax = component_correlation_ci_upper_all,
      shape = n_train_factor, group = n_train_factor)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = 0.45), linewidth = 0.30) +
  facet_grid(
    scenario ~ r2_panel + p_panel,
    labeller = labeller(
      scenario = label_value,
      r2_panel = label_parsed,
      p_panel = label_parsed
    )
  ) +
  labs(
    title = "Component recovery by design condition",
    x = NULL,
    y = expression(Delta*" component correlation"),
    shape = "Training sample size"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# Figure 3: operating characteristics and tuning pressure ---------------------
fig3_data <- master_results_by_condition %>%
  select(
    condition_id, scenario, n_train, p_per_set, target_r2,
    final_stable_rate, final_nonconvergence_rate, final_interpolation_rate,
    unrestricted_minimum_exclusion_rate, any_boundary_rate
  ) %>%
  pivot_longer(
    cols = c(
      final_stable_rate,
      final_nonconvergence_rate,
      final_interpolation_rate,
      unrestricted_minimum_exclusion_rate,
      any_boundary_rate
    ),
    names_to = "metric",
    values_to = "rate"
  ) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels, labels = scenario_labels),
    metric = factor(
      metric,
      levels = c(
        "final_stable_rate",
        "final_nonconvergence_rate",
        "final_interpolation_rate",
        "unrestricted_minimum_exclusion_rate",
        "any_boundary_rate"
      ),
      labels = c(
        "Stable final fit",
        "Final nonconvergence",
        "Final interpolation",
        "CV minimum excluded",
        "Any tuning boundary"
      )
    ),
    n_train_factor = factor(n_train, levels = c(150, 300), labels = c("N = 150", "N = 300")),
    p_panel = factor(p_per_set, levels = c(5, 20), labels = c("P[k] == 5", "P[k] == 20")),
    r2_panel = factor(target_r2, levels = c(0.30, 0.60), labels = c("R^2 == .30", "R^2 == .60"))
  )

figure3_operating_characteristics <- ggplot(
  fig3_data,
  aes(x = scenario, y = rate, shape = n_train_factor, group = n_train_factor)
) +
  geom_point(position = position_dodge(width = 0.45), size = 2) +
  facet_grid(
    metric ~ r2_panel + p_panel,
    labeller = labeller(
      metric = label_value,
      r2_panel = label_parsed,
      p_panel = label_parsed
    ),
    scales = "free_y"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Operating characteristics and tuning pressure",
    x = NULL,
    y = "Rate",
    shape = "Training sample size"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# Figure 3b: compact stability/tuning summary by design -----------------------
fig3b_data <- operational_stability_by_design %>%
  select(
    n_train, p_per_set, target_r2,
    mean_stable_rate,
    mean_cv_exclusion_rate,
    mean_any_boundary_rate
  ) %>%
  pivot_longer(
    cols = c(mean_stable_rate, mean_cv_exclusion_rate, mean_any_boundary_rate),
    names_to = "metric",
    values_to = "rate"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("mean_stable_rate", "mean_cv_exclusion_rate", "mean_any_boundary_rate"),
      labels = c("Stable final fit", "CV minimum excluded", "Any tuning boundary")
    ),
    design_cell = paste0("N=", n_train, ", Pk=", p_per_set, ", R2=", target_r2)
  )

figure3b_stability_tuning_compact <- ggplot(
  fig3b_data,
  aes(x = design_cell, y = rate)
) +
  geom_col(width = 0.70) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Compact summary of stability and tuning pressure",
    x = NULL,
    y = "Mean rate across truth conditions"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

# -----------------------------------------------------------------------------
# 6. Save figures
# -----------------------------------------------------------------------------

save_figure <- function(plot, filename, width = 8, height = 5.5, dpi = 300) {
  ggsave(file.path(FIG_DIR, paste0(filename, ".png")), plot = plot, width = width, height = height, dpi = dpi)
  ggsave(file.path(FIG_DIR, paste0(filename, ".pdf")), plot = plot, width = width, height = height)
}

save_figure(figure1_prediction_delta_r2, "figure1_prediction_delta_r2", width = 9, height = 6)
save_figure(figure1b_all_vs_stable, "figure1b_all_vs_stable_delta_r2", width = 8, height = 4.8)
save_figure(figure2_component_recovery_overview, "figure2_component_recovery_overview", width = 8.5, height = 4.5)
save_figure(figure2b_component_recovery_detailed, "figure2b_component_recovery_detailed", width = 12, height = 7.5)
save_figure(figure3_operating_characteristics, "figure3_operating_characteristics", width = 12, height = 10)
save_figure(figure3b_stability_tuning_compact, "figure3b_stability_tuning_compact", width = 9, height = 7)

# -----------------------------------------------------------------------------
# 7. Console summaries
# -----------------------------------------------------------------------------

cat("\nAggregate prediction by scenario:\n")
print(aggregate_prediction_by_scenario)

cat("\nAggregate component recovery by scenario and predictor set:\n")
print(aggregate_component_by_scenario_set)

cat("\nOperational stability by design cell:\n")
print(operational_stability_by_design)

cat("\nOutputs written to: ", normalizePath(OUT_DIR), "\n", sep = "")


# -----------------------------------------------------------------------------
# Selected tuning heatmaps
# -----------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)   # optional, only for combining the 4 panels

RESULTS_ROOT <- "C:/Users/kims15/Desktop/Kernel ERA Claude code/results_main_final_report_instability_dimension_scaled"

stab <- read.csv(file.path(RESULTS_ROOT, "condition_operational_stability.csv"))

# --- Build tidy factors for the axes ---------------------------------------
hm <- stab %>%
  transmute(
    condition_id,
    scenario = factor(scenario, levels = c("linear", "mixed", "nonlinear")),
    # y-axis label: one row per (N, p, R2) design cell
    design_cell = factor(
      sprintf("N=%d, p=%d, R\u00B2=%.1f", n_train, p_per_set, target_r2)
    ),
    median_lambda,
    median_sigma,
    any_boundary_rate,
    unrestricted_minimum_exclusion_rate
  )

# Order the y-axis sensibly (by N, then p, then R2)
cell_order <- hm %>%
  distinct(design_cell, .keep_all = FALSE) %>%
  pull(design_cell)
# (if you want a specific order, set levels explicitly here)

# --- A reusable heatmap function -------------------------------------------
make_heatmap <- function(df, fill_var, title, log_fill = FALSE,
                         fill_lab = NULL, percent = FALSE) {
  p <- ggplot(df, aes(x = scenario, y = design_cell,
                      fill = .data[[fill_var]])) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = if (percent)
                    scales::percent(.data[[fill_var]], accuracy = 1)
                  else
                    formatC(.data[[fill_var]], format = "g", digits = 2)),
              size = 2.8) +
    labs(title = title, x = NULL, y = NULL, fill = fill_lab) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(face = "bold"),
          plot.title  = element_text(face = "bold"))

  if (log_fill) {
    p <- p + scale_fill_viridis_c(trans = "log10", option = "C",
                                  labels = scales::label_number())
  } else {
    p <- p + scale_fill_viridis_c(option = "C",
                                  labels = if (percent) scales::percent else waiver())
  }
  p
}

# --- The four heatmaps ------------------------------------------------------
p_lambda <- make_heatmap(hm, "median_lambda", NULL, log_fill = TRUE) +
  labs(title = expression("Selected " * lambda * " (median)"), fill = expression(lambda))

p_sigma  <- make_heatmap(hm, "median_sigma", NULL, log_fill = TRUE) +
  labs(title = expression("Selected " * sigma * " (median)"), fill = expression(sigma))

p_bound  <- make_heatmap(hm, "any_boundary_rate",
                         "Any-boundary rate", percent = TRUE,
                         fill_lab = "Rate")

p_excl   <- make_heatmap(hm, "unrestricted_minimum_exclusion_rate",
                         "CV-minimum exclusion rate", percent = TRUE,
                         fill_lab = "Rate")

# --- Save individually ------------------------------------------------------
ggsave("heatmap_selected_lambda.pdf", p_lambda, width = 5, height = 5)
ggsave("heatmap_selected_sigma.pdf",  p_sigma,  width = 5, height = 5)
ggsave("heatmap_boundary_rate.pdf",   p_bound,  width = 5, height = 5)
ggsave("heatmap_cv_exclusion_rate.pdf", p_excl, width = 5, height = 5)