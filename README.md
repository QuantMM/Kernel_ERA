# Kernel ERA R Implementation and Simulation Framework

This repository implements the current Kernel Extended Redundancy Analysis (Kernel ERA) model and provides a resumable Monte Carlo simulation framework for evaluating prediction, component recovery, tuning behavior, and finite-sample estimation stability.

## Files

### Core implementation

* `kernel_era_core.R`

  * divisor-(N) standardization
  * Linear ERA initialization and comparison model
  * Gaussian, polynomial, and linear kernels
  * double-centering of training Gram matrices
  * training-centered test-to-training cross-kernels
  * Kernel ERA alternating least-squares estimation
  * RKHS-norm penalty used in the current manuscript
  * optional dual-coefficient (L_2) penalty for sensitivity analysis
  * K-fold cross-validation
  * MATLAB-style projection GCV for replication and diagnosis
  * out-of-sample prediction
  * component sign orientation
  * pairs bootstrap with sign alignment

* `kernel_era_simulation.R`

  * linear, mixed, nonlinear, and weakly nonlinear data-generating mechanisms
  * independent training and large test samples
  * Linear ERA versus Gaussian Kernel ERA comparison
  * oracle-signal reference
  * prediction, component-recovery, coefficient, tuning-boundary, convergence, and interpolation metrics
  * optional bootstrap coverage assessment
  * resumable condition-by-replication execution
  * optional Windows-compatible multisession parallelization through `future` and `future.apply`

### Tests and diagnostics

* `run_kernel_era_unit_tests.R`

  * lightweight checks of divisor-(N) scaling, kernel centering, spectral/direct update equivalence, one-observation prediction, stored fitted values, and multivariate outcomes

* `run_kernel_era_quick_check.R`

  * small diagnostic checking standardization, kernel centering, component normalization, prediction, tuning, and one complete replication

* `run_kernel_era_stability_pilot_30rep.R`

  * focused stability pilot for the three core linear, mixed, and nonlinear scenarios
  * evaluates minimum-CV tuning among fully converged candidates
  * checks final convergence, interpolation, tuning boundaries, paired prediction, and component recovery

### Main simulation

* `run_kernel_era_main_simulation_reporting_instability.R`

  * current recommended main simulation script
  * supports both short pilot and full main-simulation modes
  * uses dimension-adjusted Gaussian bandwidth grids
  * retains nonconvergence and interpolation as substantive simulation outcomes
  * produces both all-operational and stable-replication summaries
  * saves resumable condition-by-replication results

Older exploratory simulation scripts may be retained for reproducibility, but the current workflow is based on `run_kernel_era_main_simulation_reporting_instability.R`.

## Current model specification

The current implementation estimates

$$
\min_{\{\alpha_k\},B}
\left\|
Y-\sum_{k=1}^K K_k\alpha_k b_k^\top
\right\|_F^2
+
\lambda\sum_{k=1}^K
\alpha_k^\top K_k\alpha_k
$$

subject to

$$
\frac{1}{N}
\left\|
K_k\alpha_k
\right\|_2^2
=
1,
\qquad
k=1,\ldots,K.
$$

Here, \(K_k\) is the centered Gram matrix for predictor set \(k\), \(\alpha_k\) is the corresponding vector of kernel coefficients, and \(b_k\) is the row of outcome coefficients associated with the \(k\)th predictor-set component.

For fixed \(b_k\), define the partial residual matrix

$$
R_k
=
Y-\sum_{\ell\neq k}
K_\ell\alpha_\ell b_\ell^\top.
$$

The provisional RKHS-penalty update is

$$
\widetilde{\alpha}_k
=
\left[
(b_k^\top b_k)K_k+\lambda I_N
\right]^{-1}
R_kb_k.
$$

This is followed by normalization so that the component \(K_k\alpha_k\) has divisor-\(N\) variance one.

The default Gaussian kernel is

$$
\kappa(x,z)
=
\exp\left\{
-\frac{\|x-z\|^2}{2\sigma^2}
\right\}.
$$

Use

```r
gaussian_denominator = "sigma_squared"
```

only for sensitivity analysis using the alternative parameterization without the factor (2).

## Important implementation choices

### 1. Training-only standardization

Each predictor and outcome variable is standardized using the training-sample mean and divisor-(N) standard deviation.

Validation and test observations reuse the corresponding training-sample quantities. They are never standardized using their own means or standard deviations.

### 2. Correct out-of-sample kernel evaluation

For a new observation, the implementation constructs a test-to-training cross-kernel.

It does not construct a kernel using only the test observations. Cross-kernel centering uses the feature-space centering quantities obtained from the training Gram matrix.

### 3. Shared tuning parameters

The simulation uses one ridge parameter (\lambda) and one Gaussian bandwidth (\sigma) shared across all predictor sets, matching the current manuscript specification.

The fitting functions can accept set-specific bandwidths, but the supplied simulation and tuning procedures deliberately select shared scalar values.

### 4. Cross-validation tuning

The primary tuning procedure is five-fold cross-validation based on validation prediction mean squared error.

The selected tuning combination is the minimum-CV candidate among those that converged in every validation fold:

```r
cv_selection = "minimum"
cv_require_all_folds_converged = TRUE
```

The one-standard-error rule is not used as the primary procedure because ordering Kernel ERA model complexity jointly through (\lambda) and (\sigma) is not straightforward, and preliminary analyses showed that the implemented one-standard-error rule could produce excessive smoothing.

`kernel_era_matlab_gcv()` reproduces the MATLAB projection-GCV calculation. It is retained for replication and methodological comparison, but it is not used as the default simulation tuner.

### 5. Dimension-adjusted Gaussian bandwidth grids

The main simulation contains predictor sets with different dimensions. Because Euclidean distances among standardized observations increase approximately in proportion to (\sqrt{P_k}), the Gaussian bandwidth grid is scaled according to

$$
\sigma(P_k)
===========

\sigma_{\mathrm{base}}
\sqrt{\frac{P_k}{5}}.
$$

The base grid for (P_k=5) is

```r
BASE_SIGMA_GRID <- c(
  0.75, 1.5, 3, 6, 12, 24, 48, 96
)
```

This produces the following grids:

```text
P_k = 5:
0.75, 1.5, 3, 6, 12, 24, 48, 96

P_k = 20:
1.5, 3, 6, 12, 24, 48, 96, 192
```

The current ridge grid is

```r
LAMBDA_GRID <- 10^seq(-4, 3, by = 1)
```

The very small value (\lambda=10^{-5}) was excluded after preliminary checks repeatedly identified it as part of a pathological small-(\lambda), small-(\sigma) interpolation region. Smaller bandwidth values are otherwise retained so that finite-sample instability in difficult conditions can be observed and reported.

### 6. Spectral solver

The default

```r
solver = "spectral"
```

uses the eigendecomposition of each centered Gram matrix.

For the RKHS update, it is algebraically equivalent on the kernel range to solving the corresponding direct linear system, while avoiding repeated (N\times N) factorizations at every alternating least-squares iteration.

Use

```r
solver = "direct"
```

for a literal matrix-solve implementation.

### 7. Multiple final-model starts

The final Kernel ERA model is estimated using three initializations:

```r
final_n_starts = 3L
```

The first start uses the Linear ERA solution as a warm start, and the remaining starts use random initialization. The best converged solution is retained.

The current numerical settings are

```r
max_iter = 500L
tol = 1e-6
objective_tol = 1e-9
```

### 8. Sign orientation

A predictor-set component and its outcome coefficient may both be multiplied by (-1) without changing the fitted contribution.

Simulation recovery and bootstrap summaries therefore align each fitted component to a reference component before reporting coefficients or component-recovery measures.

### 9. Estimation instability as a simulation outcome

Nonconvergence and interpolation in difficult simulation conditions are not automatically removed or replaced by another tuning solution.

A final fit is classified as stable when

```r
final_converged & !final_interpolation_flag
```

where interpolation is defined by a training FIT greater than the prespecified threshold:

```r
interpolation_fit_threshold = 0.98
```

The simulation reports:

* final convergence and nonconvergence rates
* interpolation rates
* stable-estimation rates
* combined operational-failure rates
* exclusion of nominal CV minima because of fold nonconvergence
* tuning-boundary selection rates
* prediction and component recovery for all operational outputs
* prediction and component recovery conditional on stable Kernel ERA estimation

This distinction allows estimation instability under small-sample, relatively high-dimensional, or low-signal conditions to be treated as a substantive finite-sample result rather than hidden through automatic fallback procedures.


## Main output files

### Design and run settings

* `main_simulation_design.csv`
* `run_settings.csv`
* `sigma_grids_by_predictor_dimension.csv`

### Combined raw results

* `combined_performance.csv`
* `combined_components.csv`
* `combined_coefficients.csv`
* `combined_tuning.csv`
* `combined_selected.csv`
* `combined_inference.csv`
* `combined_error.csv`
* `main_results_combined.rds`

### Performance and recovery summaries

* `summary_performance_all_operational.csv`
* `summary_performance_kernel_stable.csv`
* `summary_performance_combined.csv`
* `summary_components_all_operational.csv`
* `summary_components_kernel_stable.csv`
* `summary_components_combined.csv`
* `summary_tuning.csv`
* `main_summaries_all_and_stable.rds`

### Paired Kernel ERA versus Linear ERA comparisons

* `paired_prediction_by_replication.csv`
* `paired_prediction_summary_all_operational.csv`
* `paired_prediction_summary_kernel_stable.csv`
* `paired_prediction_summary_combined.csv`
* `paired_components_by_replication.csv`
* `paired_component_summary_all_operational.csv`
* `paired_component_summary_kernel_stable.csv`
* `paired_component_summary_combined.csv`

### Stability and diagnostic outputs

* `selected_tuning_with_design_and_status.csv`
* `condition_operational_stability.csv`
* `condition_operational_and_predictive_summary.csv`
* `boundary_cases.csv`
* `nonconverged_final_fits.csv`
* `interpolating_final_fits.csv`
* `final_operational_failures.csv`
* `excluded_unrestricted_cv_minima.csv`
* `pipeline_integrity.csv`

`pipeline_integrity.csv` evaluates whether the requested tasks and post-processing steps completed successfully. Observed nonconvergence or interpolation is treated as a simulation result and does not, by itself, constitute a pipeline failure.

## Main design currently encoded

The main simulation contains 24 conditions crossing:

* functional form:

  * linear
  * mixed linear and nonlinear
  * fully nonlinear
* training sample size:

  * (N=150)
  * (N=300)
* predictors per set:

  * (P_k=5)
  * (P_k=20)
* target population signal level:

  * (R^2=.30)
  * (R^2=.60)

All conditions use:

* three prespecified predictor sets
* an independent test sample of 1,000 observations
* within-set predictor correlation of (.40)
* between-set predictor correlation of (.10)
* nonlinearity-strength setting of (.35)

The main script currently uses 200 Monte Carlo replications per condition. The run is resumable at the condition-replication level.
