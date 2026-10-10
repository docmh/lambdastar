# lambdastar — Measurement Core, Model Adapter, and Bootstrap

A standalone R package for the measurement and projection core of Lambda Star,
featuring numerical parsimony functional and person/cluster bootstrap with direct
percentile intervals. Version **0.8.2** uses the manuscript estimator by default
for point estimates and all draws. The alternative population estimator is
explicitly selectable. Nuisance effects are jointly accounted for. All intervals
use `ci.method="percentile"`. The [model space learning path](inst/doc/MODEL_SPACE_WORKFLOW.md),
additional user guides, and executable examples are installed with the package.

Methodological background: Mike Hammes (2026),
[*From Null Rejection to Structured Localization: Lambda Star and the Thermodynamics of Frequentist Statistics*](https://doi.org/10.5281/zenodo.22962377)
(Preprint, DOI: 10.5281/zenodo.22962377).

## Installation and Automated Tests

The published version is available on CRAN:

```r
install.packages("lambdastar")
```

The [source code on GitHub](https://github.com/docmh/lambdastar) contains the
published package and its tests. [GitHub Actions](https://github.com/docmh/lambdastar/actions)
checks macOS on Apple Silicon with R release, oldrel, and devel, as well as Linux
with R release and devel, and Windows with R release. The macOS jobs verify
that R runs natively on ARM64. Check logs and platform details are available as
artifacts. The CI check uses `--as-cran --no-manual`; the PDF manual check is
performed separately before CRAN submission.

## Local Installation

Requirements: R >= 4.1 with `stats` and the recommended R package `boot`.
A source archive can be created with `R CMD build .` in the package directory:

```r
install.packages("lambdastar_0.8.2.tar.gz", repos = NULL, type = "source")
library(lambdastar)
packageVersion("lambdastar")
```

Installation uses the first writable standard library or the personal library
offered by R. With `lib = "/path/to/library"`, an alternative destination can
be chosen. The package is licensed under [Apache 2.0](inst/doc/LICENSE-2.0.txt).
Copyright 2026 Mike Hammes. Contact: <mike.hammes@mikehammes.name>.
A copy of the license is installed at `system.file("doc", "LICENSE-2.0.txt",
package="lambdastar")`.

A standalone example and learning guide are installed:

```r
source(system.file("examples", "quickstart.R", package = "lambdastar"))
system.file("doc", "LEARN_LAMBDASTAR.md", package = "lambdastar")
help("lambdastar")
```

## Learning Lambda Star Step by Step

The [executable tutorial](inst/doc/LEARN_LAMBDASTAR.md) explains measurement noise,
λ*, localization in X_H, actual bootstrap percentiles, and F_H using a consistent
example with known structure. Seven figures, prediction questions, and adjustable
parameters make the computation process transparent. An [already executed version](inst/tutorial/WALKTHROUGH.md)
with default values is also included.

```r
source(system.file("examples", "learn_lambdastar.R", package="lambdastar"))
lesson <- lambdastar_tutorial(output_dir="my-lambda-star-lesson")
```

In interactive R, execution pauses for a prediction before each evaluation.
`pause=FALSE, show_plots=FALSE` runs it without dialog and graphics windows.
The output directory contains the illustrated `WALKTHROUGH.md`, CSV tables,
synthetic data, bootstrap draws, and `lesson.rds`. Existing results are not
overwritten. Default: B=5000 and 95%; for short practice runs, e.g., `B=199`
can be explicitly set.

[Advanced examples](inst/doc/ADVANCED_EXAMPLES.md) separately cover population
estimation, location adjustment, school clusters, and the two null cases.
The [model space learning path](inst/doc/MODEL_SPACE_WORKFLOW.md) then leads
from research question to X_H, conditional effects, contrasts, basic recipes, and
F_H comparisons. Its scope table delineates currently supported hypotheses and
reference spaces. A separate example shows available and undefined intervals on two
synthetic E4 datasets:

```r
source(system.file("examples", "intervals_and_targets.R", package="lambdastar"))
intervals <- lambdastar_intervals_example()  # B=5000, 95%, seed=1
intervals$intervals
intervals$fits$weak$bootstrap$draw_reasons
```

## Persons Across Multiple Measurement Time Points

`lambda_panel()` supplements fixed person/timepoint models. Each row contains
a person at a time point and at least two parallel measurement indicators.
Person intercepts are adjusted by default through group centering; additional
target and adjustment terms come from the declared model. The bootstrap draws
complete person trajectories with separate identities for multiply-drawn copies.
[Derivation and learning path](inst/doc/PANEL_MODELS.md).

```r
source(system.file("examples", "panel_models.R", package="lambdastar"))
panel <- lambdastar_panel_example()  # two analyses with B=5000, 95%, seed=1
panel$trend$estimates
confint(panel$interaction)         # Group × Time after person and time adjustment
panel$comparison$delta_F           # Time trend versus curvature with same nuisance space
panel$trend$panel$map               # actual person/observation mapping
```

The measurement bath must be homogeneous and isotropic in the remaining
observation space. Time points are not parallel replicates of an unchanged response.
Unequal numbers of observations are possible; explicit complete-case selection
removes entire persons. T_q remains on the scale of the observed item mean.

## Mean Against Reference and Paired Differences

Since 0.7.0, `lambda_reference()` captures deviation from a specified reference,
including mean direction. `lambda_paired()` forms `after - before` for each
parallel item pair and estimates T from these differences. Derivation, resolution
convention, and example: [Reference means](inst/doc/REFERENCE_MEANS.md).

```r
source(system.file("examples", "reference_means.R", package="lambdastar"))
means <- lambdastar_reference_example()  # B=5000, 95%, seed=1
means$one_sample$estimates              # Mean against 50
means$paired$estimates                  # Mean change against 0
means$paired_reference$estimates        # Change against .5
means$paired$mean_departure             # signed mean deviation
confint(means$paired)
```

With fixed origin and no nuisance, d=N, `~1` describes the free mean direction
and `~0` the fixed null value. A comparison requires the same reference vector
and nuisance space. T_q is explicitly specified on the mean or difference scale.
For pairs, at least two parallel measurement pairs are required; a single before/after
column does not identify a separate measurement bath. Standard centered calls
remain unchanged.

## Brief Example

```r
library(lambdastar)
set.seed(17)
d <- data.frame(x = rnorm(80), z = rnorm(80),
                group = factor(rep(c("a", "b"), 40)),
                condition = factor(rep(c("a", "b", "c", "d"), 20)))
signal <- d$x + 0.5 * d$z
d[paste0("item", 1:4)] <- replicate(4, signal + rnorm(80))
items <- paste0("item", 1:4)
d$y_mean <- rowMeans(d[items])
T_q <- 1/12  # Example assumption on mean scale, no estimated resolution

measurement <- lambda_measure(d, items)
two_group <- lambda_model(d, ~ group, items, T_q)       # Group mean hypothesis
anova <- lambda_model(d, ~ condition, items, T_q)       # Full model
ancova <- lambda_model(d, ~ condition + x, items, T_q)
regression <- lambda_model(d, y_mean ~ x + z, items, T_q)
interaction <- lambda_model(d, ~ group * x, items, T_q)
summary(regression)
regression$parsimony[c("rss", "dimension_cost", "geometry_cost", "free_energy")]

# Same coding via explicit X_H or a checked lm object:
via_matrix <- lambda_model(d, as.matrix(d[c("x", "z")]), items, T_q)
via_lm <- lambda_model(d, lm(y_mean ~ x + z, d), items, T_q)
stopifnot(all.equal(regression$parsimony$free_energy,
                    via_lm$parsimony$free_energy))

# One data object, two models, same single-model computation twice:
comparison <- lambda_compare(d, ~ x, ~ x + z, items, T_q)
print(comparison)
stopifnot(all.equal(comparison$delta_F,
  comparison$model_b$parsimony$free_energy - comparison$model_a$parsimony$free_energy))

# Same space, different coding: projection identical, geometry possibly different.
coding <- lambda_compare(d, ~ condition, ~ condition, items, T_q,
  contrasts_a = list(condition = "contr.treatment"),
  contrasts_b = list(condition = "contr.helmert"))
print(coding)
```

Items or replications must be correctly coded already and measure the same stable
structure. The core assumes parallel measurements with independent, homogeneous
measurement errors. It does not empirically check these measurement assumptions
and does not automatically interpret temporal changes as noise. Selected inputs
must be complete and finite. Exclusion of missing cases must be explicitly
requested.

## Model Adapters and Common Cases

- `lambda_measure(data, indicators)` estimates the measurement state. Only here
  does `indicators = NULL` mean all columns are measurements.
- `lambda_model(data, model, indicators, T_q)` evaluates a complete model.
  `model` is a `lambda_design`, a formula, a suitable `lm`, or a numeric matrix.
  One-sided formulas use the item mean. For two-sided formulas, the response
  must match; deviations trigger an error.
- `lambda_compare(data, model_a, model_b, indicators, T_q)` prepares common
  data and measurement conditions once and calls the same evaluation path twice.
  `delta_F = F_Hb - F_Ha`; a negative delta means smaller F_H for b.
  This is not a significance test and not automatic model selection.

`na.action = "fail"` is the default. `na.action = "complete"` selects complete
persons across items and all variables in both models once. Unneeded data columns
do not affect selection. Excluded positions appear in `$cases$excluded`, used ones
in `$cases$included`; row identities and required variables are also stored.
Two separate single calls are comparable only when missingness is present if
both are evaluated on exactly this common subset.

Case selection occurs **before** transformations like `poly()` or `scale()`.
NA/Inf arising from transformations are rejected, not additionally filtered
per model. Formula variables must come from `data`. An explicit matrix contains
all input rows in data order; existing row names must match. Without row names,
the specified order is used. Missing matrix values enter the common selection.

`lm` requires a stored model frame, an intercept, the same cases in the same
order, a matching mean response, and matching coded predictor values. Stored
contrasts apply even after changes to global R contrast options. Weights, offsets,
GLM, multivariate responses, and other model classes are rejected. Differing
`lm` samples are explicitly refitted; the adapter does not perform a silent refit.

An intercept does not count as a hypothesis direction. `~ 1` is the null model.
Factor levels and actual contrasts are preserved. Explicit contrast matrices may
have fewer than G−1 columns; these exact directions are preserved in the bootstrap.
Rank deficiency, lost levels, or a redundant dummy block after centering like
`~ 0 + factor` cause an error rather than silently changing coding. `$model`
contains formula, columns, term mapping, contrasts, factor levels, rank, and
cases; `$basis$X` contains the centered, possibly nuisance-adjusted coded design.
Only the projection uses the QR basis.

The manuscript path assumes stable structure plus homogeneous isotropic measurement
bath with finite second moments. Neither normality of response nor equal group
variances is required. Group comparisons (including alongside Welch), fixed ANOVA/ANCOVA,
regression, and interactions use the same X_H. Welch test statistics are computed
separately; they replace neither T nor the projection correction. The population
estimator retains its additional assumptions.

`lambda_design()` binds a reusable model to coding and input rows. `lambda_hypothesis()`
selects target and adjustment terms from the full model. Conditional λ* and κ refer
to the structure remaining after adjustment. Guide: [Model spaces](inst/doc/MODEL_SPACES.md).

```r
spec <- lambda_design(d, ~ group * x)
effect <- lambda_hypothesis(d, ~ group * x, target="group:x", adjust=c("group", "x"))
conditional <- lambda_model(d, effect, items, T_q)
```

`lambda_contrast()` encodes named homogeneous constraints `L beta=0`;
`lambda_marginal_contrast()` uses explicitly set predictor values and weights
for simple or averaged effects. `encoding="fixed"` keeps the basis learned after
case selection fixed during bootstrap; `"reevaluate"` is the default. Guide and
executable example: [Contrasts and basic recipes](inst/doc/CONTRASTS_AND_BASES.md).

```r
h <- lambda_contrast(d, ~group*x, c(groupb=1, `groupb:x`=1))
fit <- lambda_model(d, h, items, T_q)  # b-a at x=1, conditional on the null model
source(system.file("examples", "contrasts_and_bases.R", package="lambdastar"))
lesson <- lambdastar_contrasts_example(B=199)  # short practice run; default B=5000
```

Single-sample/paired mean hypotheses use the explicit
[reference API](inst/doc/REFERENCE_MEANS.md). Standard centering still removes
the mean. Fixed person/timepoint models use
[lambda_panel](inst/doc/PANEL_MODELS.md); general covariance and mixed models
remain separate extensions.

`measurement_scale = "mean"` is the default; `"item"` changes only the
noncentrality denominator. The adapter correctly uses the mean response and
`temperature_mean` for F_H in both cases. `T_q` must be explicitly specified
as positive. `capacity`, `numerical_control`, `eta_lookup`, and `response_tol`
are overridable and are stored in `$settings`. Response/lm checking uses pointwise
`abs(a-b) <= sqrt(.Machine$double.eps) * max(1, abs(a), abs(b))` by default.
`print()` and `summary()` show point estimates and F_H decomposition.

## Percentile Intervals and Common Person Draws

`bootstrap=TRUE` uses `ci.method="percentile"`, B=5000, 95%, and seed=1.
The bootstrap draws N complete persons with replacement. Each draw computes T, λ*,
λ_H, λ_perp, and κ_H using the same selected estimator as the original data
(`estimator="paper"` as default). For 95%, the 2.5% and 97.5% percentiles are
used directly, computed with `boot::boot.ci(type="perc")`. Constant distributions
retain identical lower and upper bounds. Level, B, and seed are overridable.

**The selected point estimates are in `$estimates`** and in the corresponding
fields `$lambda_star`, `$lambda_h`, and `$kappa`, also when `bootstrap=FALSE`.
The CI table contains the same points. `$inference` and `$settings$estimator`
name the method. `$paper` preserves manuscript values; `$population` is created
only when explicitly selected.

```r
ci_model <- lambda_model(d, ~ x + z, items, T_q,
  bootstrap = lambda_bootstrap_control(B = 499, conf.level = .95, seed = 17,
                                       keep_draws = TRUE))
ci_model$estimates
confint(ci_model)
confint(ci_model, parm = c("lambda_star", "lambda_h", "kappa"), level = .90)
ci_model$bootstrap$ci

# Common person draws, fixed observed group proportions:
ci_comparison <- lambda_compare(d, ~ group, ~ group + x + z, items, T_q,
  bootstrap = lambda_bootstrap_control(B = 499, strata = "group"))
confint(ci_comparison)
```

`estimator="paper"` uses the bath-corrected projection energies from the
manuscript: on the mean scale λ_H,raw=Q_H/T−p and λ*,raw=Q_total/T−d.
Without additional nuisance variables, d=N−1. Original estimate and each
bootstrap draw use these formulas.

`estimator="population"` explicitly selects the alternative correction for
population with fixed N: λ*ₙ=(N−1)Var(S)/T and λ_H,ₙ from the linear
population projection. It uses `V_S=Var(Y_mean)-T_mean` and
`V_H=Var(Y_mean)-RSS/(N-1-p)` under linear conditional expectation and
homogeneous residual variance. It currently supports intercept adjustment and
person bootstrap. Nuisance or cluster combinations are explicitly rejected with
this option. `help("lambda_inference")` contains the formulas.

### Nuisance Effects and School Clusters as in the Manuscript

```r
# sses: one row per person, five preprocessed WHO-5 items.
# skills: names of 15 predictors; original resolution of item mean.
fit_sses <- lambda_model(sses, as.matrix(sses[skills]), items, T_q = 1/300,
  nuisance = ~ factor(SiteID),
  bootstrap = lambda_bootstrap_control(B = 1000, seed = 2026073116,
    resampling = "cluster", cluster = "school_uid", strata = "SiteID",
    percentile_interpolation = "type7", keep_draws = TRUE))
fit_sses$estimates
confint(fit_sses, parm = "kappa")
```

`nuisance` accepts a one-sided formula or a numeric matrix. An intercept always
belongs to the nuisance space. Measurements and hypothesis matrix are jointly
residualized; d=N−rank(B) determines measurement correction. The nuisance space
is estimated anew in each draw. Coding and rank loss are checked. The model
formula contains the directions to be tested; nuisance variables are separate
in `nuisance`.

Clusters are drawn completely and with replacement within their strata; per
stratum the cluster count remains fixed, the person count can vary. Cluster
labels apply within each stratum. Common case selection also includes nuisance,
cluster, and stratum variables. `percentile_interpolation="type7"` computes
the unchanged R percentiles (`quantile(type=7)`) of the original analysis.
Default `"boot"` still uses `boot::boot.ci(type="perc")`; both are ordinary
percentile intervals.

| Situation | Treatment |
|---|---|
| λ*=0 | valid null value; remains in λ* distribution |
| κ with estimated total signal null | `NA`, reason `undefined_kappa` |
| κ=0 or κ=1 with positive total signal | valid values; remain in distribution |
| Original data with T=0 | singularity; no intervals |
| Bootstrap draw with T=0 | T=0 retained; dependent values undefined |
| Rank loss or missing factor level | affected model values `NA`, own diagnostic reason |

An undefined original estimate receives `undefined_estimate`. Once a draw is
undefined for a measure, its regular interval remains `NA` (`undefined_draws`
for κ without signal, else `failed_draws`). Intervals for other measures are
computed independently. `n_valid`, `n_failed`, `n_undefined`, and failure reasons
are visible; draws are not replaced. At least `min_valid=200` values are needed
for quantiles.

With `conditional_quantiles=TRUE`, quantiles can be additionally requested for
valid draws only. They are stored separately in `$bootstrap$conditional_quantiles`,
with the count of excluded draws, as **conditional descriptive quantiles**.
`confint()` continues to provide only regular confidence intervals.

`strata="group"` maintains observed group sizes during person bootstrap. Only
`estimator="population"` uses the corrected mixture variance here and requires
at least two persons per stratum. A saturated model has no residual degrees of
freedom for population projection; its affected population values are output
as `NA` with justification.

`keep_draws=TRUE` stores **the actual measures** from each draw, including zero
and NA, along with `draw_reasons`. `keep_indices=TRUE` stores original row
positions, for clusters as a variable-length list. `confint(..., level=...)` uses
existing draws and requires `keep_draws=TRUE` for this. All measures and comparison
models share the same person draws. Intervals are pointwise; F_H and ΔF remain
point estimates. Validation reports availability, width, and coverage separately;
with true null signal, coverage for κ is not defined.

## Functions and Scales

- `estimate_lambda_star_items()` provides temperature, raw and non-negative total
  structure, and row identity. Default is `scale = "composite_mean"`.
- `estimate_lambda_partition_items()` supplements projection onto an explicit `x`,
  raw/constrained decomposition, and kappa. The design contains no intercept.
- `hypothesis_basis()` preserves the coded matrix; the normalized QR basis serves
  only projection. `design_volume(logarithm = TRUE)` preserves log(V) even for
  very small volume.
- `eta_numerical()` solves quantized capacity via a logarithmic moment equation.
  `details = TRUE` shows solution and truncation diagnostics.
- `xi_eta_numerical()` integrates directly and adaptively by default. A reusable
  table from `make_eta_numerical_lookup()` accelerates many evaluations; below its
  first positive node, direct integration is used.
- `thermodynamic_model_state()` computes F_H for a single X_H. Default is
  `eta_method = "numerical"`; `"smooth"` is an explicitly marked approximation.
  Both use the same geometric activation.
- `lambda_control()` holds overridable tolerances separately for rank, capacity
  solution, truncated distribution tails, quadrature, and lookup. Lookup errors
  are estimated discretization errors, not formal bounds for the entire computation.
  Gibbs tails are checked via integral limits.

`scale = "item"` changes the denominator of noncentralities to individual
measurement variance. The stored `composite_response` and `stable_profile` are
still means. Those using the mean response in F_H must therefore use
**temperature_mean**. Scaling the response by factor c requires corresponding
scaling of T and T_q by c². The resolution T_q is specified by subject matter.

At T = 0, the core provides `status = "singular_temperature"` and NA for values
not regularly estimable including F_H. For explicit numerical boundary checks,
`boundary = TRUE` is possible: F_H = RSS, status `boundary_only`. This is not
an empirical estimate or comparison state. Kappa remains undetermined at zero
even with positive temperature and estimated total structure. Very small positive
values are not replaced by epsilon.

## Building and Testing

The package imports `stats` and `boot`. The test runner uses base R; documentation
and NAMESPACE are already included.

```sh
R CMD build .
R CMD check --as-cran lambdastar_0.8.2.tar.gz
```

Tests include an existing synthetic simulation and numerically independent
reference values generated with Python Decimal. User examples use synthetic data;
external OECD raw data are not included. The Apache 2.0 license applies to this
package; it does not license the separate manuscript or external OECD data.
