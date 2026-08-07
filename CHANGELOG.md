# Changelog

All notable changes to this repository are documented here. Versions
correspond to GitHub releases where applicable.

## [2.0] – 2026-08 – Full simulation (manuscript revision)

### Changed
- Increased Monte Carlo replications from 200 to 500 per condition
  (24 conditions × 500 = 12,000 tasks), in response to a reviewer request.
  This is the only substantive change from the conference / technical-report
  version.
- The full run was performed with `overwrite = FALSE` against the existing
  results directory, so the original 200 replications were reused on the same
  seed scheme and only replications 201–500 were newly computed. The
  resumable, condition-by-replication design made this incremental extension
  possible without recomputing the first 200 replications.
- Added `run_kernel_era_main_simulation_reporting_instability_500rep.R`
  (identical to the 200-replication script except for `N_REP_MAIN <- 500L`).
  The 200-replication script,
  `run_kernel_era_main_simulation_reporting_instability.R`, is retained
  alongside it so the 200 → 500 workflow is transparent and reproducible.

### Not changed
- Data-generating mechanisms, 24-condition design, and seed scheme.
- Estimator (RKHS-penalty Kernel ERA), lambda grid (10^-4 … 10^3),
  dimension-scaled Gaussian sigma grids, five-fold minimum-CV tuning with the
  all-fold-convergence requirement, and interpolation threshold (0.98).
- Reporting structure and all output file names (all-operational vs
  kernel-stable summaries, `pipeline_integrity.csv`,
  `condition_operational_stability.csv`, etc.); these were already in place
  for the technical-report version.

### Results
- Pipeline integrity: all 12,000 tasks completed, 0 task-level errors.
- Operational stability: 17 final nonconvergences and 0 interpolations across
  12,000 final fits; minimum condition-level stable rate 0.982.

## [1.0] – 2026-07 – Technical report / conference version

- Public release underlying the conference presentation and technical report.
- 200 replications per condition, produced with
  `run_kernel_era_main_simulation_reporting_instability.R`.
- Preserved as release `v1.0-techreport-200rep`.
