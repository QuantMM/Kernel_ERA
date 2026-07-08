# Kernel ERA R Implementation and Simulation Framework

This repository implements the current Kernel Extended Redundancy Analysis (Kernel ERA) model and provides a resumable Monte Carlo simulation framework.

## Files

- `kernel_era_core.R`
  - divisor-`N` standardization
  - linear ERA initialization and comparator
  - Gaussian, polynomial, and linear kernels
  - double-centering of training Gram matrices
  - training-centered test-to-training cross-kernels
  - Kernel ERA alternating least-squares estimation
  - RKHS-norm penalty used in the current manuscript
  - optional dual-coefficient `L2` penalty for sensitivity analysis
  - K-fold cross-validation
  - MATLAB-style projection GCV for replication/diagnosis
  - out-of-sample prediction
  - component sign orientation
  - pairs bootstrap with sign alignment

- `kernel_era_simulation.R`
  - linear, mixed, nonlinear, and weakly nonlinear data-generating mechanisms
  - independent training and large test samples
  - Linear ERA versus Gaussian Kernel ERA comparison
  - oracle-signal reference
  - prediction, recovery, tuning-boundary, convergence, and coefficient metrics
  - optional bootstrap coverage assessment
  - resumable condition-by-replication execution
  - optional Windows-compatible multisession parallelization through `future` and `future.apply`

- `run_kernel_era_unit_tests.R`
  - lightweight checks of divisor-`N` scaling, kernel centering, spectral/direct update equivalence, one-observation prediction, stored fitted values, and multivariate outcomes

- `run_kernel_era_quick_check.R`
  - small diagnostic that checks standardization, kernel centering, component normalization, prediction, tuning, and one complete replication

- `run_kernel_era_pilot_simulation.R`
  - ten replications in three core scenarios using a coarse grid and three-fold cross-validation

- `run_kernel_era_main_simulation.R`
  - main-study template using the proposed factorial design, five-fold cross-validation, and resumable output

## Current manuscript specification

The default implementation is

$$
\min_{\{\alpha_k\},B}
\left\|Y-\sum_{k=1}^K K_k\alpha_k b_k^\top\right\|_F^2
+\lambda\sum_{k=1}^K\alpha_k^\top K_k\alpha_k,
$$

subject to

$$
\frac{1}{N}\|K_k\alpha_k\|_2^2=1.
$$

For a fixed outcome-coefficient row `b_k`, the provisional RKHS-penalty update is

$$
\widetilde\alpha_k=
\left[(b_k^\top b_k)K_k+\lambda I_N\right]^{-1}R_kb_k,
$$

followed by normalization to make the component have divisor-`N` variance one.

The default Gaussian kernel is

$$
\kappa(x,z)=\exp\left\{-\frac{\|x-z\|^2}{2\sigma^2}\right\}.
$$

Use `gaussian_denominator = "sigma_squared"` only for a sensitivity analysis using the alternative parameterization without the factor `2`.
## Important implementation choices

### 1. Training-only standardization

Each predictor and outcome is standardized using the training-sample mean and divisor-`N` standard deviation. Validation and test observations reuse those training quantities. They are never standardized using their own means and standard deviations.

### 2. Correct out-of-sample kernel evaluation

For a new observation, the code constructs a test-to-training cross-kernel. It does not construct a kernel only among test observations. Cross-kernel centering uses the feature-space centering quantities estimated from the training Gram matrix.

### 3. Shared tuning

The simulation uses one `lambda` and one Gaussian `sigma` shared across all predictor sets, matching the current manuscript. The fitting function can accept a vector of set-specific bandwidths, but the supplied tuning functions deliberately select a shared scalar.

### 4. K-fold cross-validation as the primary tuning method

`kernel_era_cv()` minimizes validation prediction mean squared error. This is the default for all simulation claims.

`kernel_era_matlab_gcv()` reproduces the MATLAB projection-GCV calculation. It is retained for replication and comparison, not used as the default simulation tuner. In the projection formula, the trace is generally close to the number of fitted components and does not fully represent the adaptive, penalized construction of those components.

### 5. Spectral solver

The default `solver = "spectral"` uses the eigendecomposition of each centered Gram matrix. For the RKHS update, it is algebraically equivalent on the kernel range to solving the MATLAB linear system, while avoiding repeated `N x N` factorizations at every alternating least-squares iteration. Set `solver = "direct"` for a literal matrix-solve implementation.

### 6. Sign orientation

A component and its outcome coefficient may both be multiplied by `-1` without changing fitted values. Simulation recovery and bootstrap summaries therefore align each fitted predictor-set component to a reference component before reporting coefficients or recovery.

## Suggested workflow

Run the unit tests and then the diagnostic:

```r
source("run_kernel_era_unit_tests.R")
source("run_kernel_era_quick_check.R")
```

Then run the pilot:

```r
source("run_kernel_era_pilot_simulation.R")
```

Inspect at least the following pilot outputs:

- convergence rates
- selected `lambda` and `sigma`
- lower- and upper-boundary selection rates
- test prediction performance
- predictor-set component recovery
- elapsed time per replication

Only after the pilot should the main run be launched:

```r
source("run_kernel_era_main_simulation.R")
```

## Resuming a run

Each condition-replication result is saved separately under `results_dir/replications`. With `overwrite = FALSE`, completed replications are loaded rather than recomputed.

## Main output files

- `performance.csv`: test mean squared error, test R-squared, prediction correlation, train FIT, train-test FIT gap, tuning values, and convergence
- `components.csv`: train/test correlation with the true predictor-set component and calibrated test component root mean squared error
- `coefficients.csv`: oriented outcome coefficients and bias on the standardized outcome scale
- `tuning.csv`: complete cross-validation surface for every replication
- `selected.csv`: selected tuning parameters and tuning-grid boundary indicators
- `inference.csv`: bootstrap standard errors, percentile intervals, and coverage when bootstrap is enabled
- `error.csv`: failed replications, if any
- `combined_results.rds`: all combined results

## Main design currently encoded

The main grid crosses:

- truth: linear, mixed linear/nonlinear, fully nonlinear
- training sample size: 150 and 300
- predictors per set: 5 and 20
- target population signal-to-noise level: `R2 = .30` and `.60`

with three prespecified predictor sets, a test sample of 1,000, within-set correlation `.40`, and between-set correlation `.10`.

This design is a defensible starting point rather than a locked manuscript decision. The final number of conditions and Monte Carlo replications should be set after timing and tuning-boundary diagnostics from the pilot.

## Optional penalty sensitivity analysis

The manuscript-matching default is

```r
penalty = "rkhs"
```

An optional dual-coefficient penalty is available as

```r
penalty = "dual_l2"
```

Its update solves

\[
\left[(b_k^\top b_k)K_k^2+\lambda I_N\right]\alpha_k
=K_kR_kb_k.
\]

This option should be treated as a separate method or sensitivity analysis, not mixed with the primary RKHS specification.
