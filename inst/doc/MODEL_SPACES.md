# Reusable designs and conditional hypotheses

lambdastar delivers the centered model-space implementation and its
[learning path and supported-scope table](MODEL_SPACE_WORKFLOW.md). The core input is
the actual hypothesis encoding X_H, together with its nuisance reference and
measurement data. Formula names do not select different lambda-star methods.

## The measurement axiom

The manuscript estimator evaluates lambda_star, lambda_H, kappa_H and numerical
F_H under stable structure plus zero-mean homogeneous, isotropic background
fluctuation with finite second moments. Parallel measurements identify that
background under the independent-replicate assumption. Fitting uses this axiom;
it does not demonstrate it for a particular instrument.

Normality or equal group variances of the observed response are not required.
Different stable within-group dispersions can contribute to lambda_star without
being captured by a group-mean contrast. If stable group means are equal, the
true group-mean lambda_H can be zero while total structure is positive. A
separately calculated Welch comparison can accompany this analysis; its t value
and degrees of freedom are not inputs to the lambda-star calculation.

The explicit `estimator="population"` option keeps its separate linear
conditional-mean/homogeneous-residual assumptions and existing limitations.
This increment does not widen that estimator's supported scope.

## A reusable specification

```r
design <- lambda_design(d, ~group*x + z)
fit <- lambda_model(d, design, items, T_q=1/12)
saveRDS(design, "design.rds")
fit_again <- lambda_model(d, readRDS("design.rds"), items, T_q=1/12)
```

`lambda_design()` accepts an ordinary supported formula/lm or a numeric matrix.
It records schema version, original row identities, the relevant input columns,
coding defaults, and optional nuisance specification. Formula matrices are built
after shared case selection; inspect `fit$encoded_x`, `fit$encoded_nuisance`,
`fit$model` and `fit$basis$X` for the realized encoding and adjusted geometry.
Matrix inputs are already encoded and appear in `design$X_H`.

Changing predictors or row order requires constructing a new design. For a
one-sided formula, new measurement values on unchanged predictors are allowed;
temperature and projections are always reestimated. A two-sided formula's
response must still match the measurement composite. Objects preserve standard
contrast defaults independently of later global options. Custom functions and
custom contrast functions should be evaluated into numeric matrices explicitly.
Supported recipes use deterministic base/stats terms; spline calls can be
namespace-qualified. No caller-local functions or constants are captured.

## Select a target block and its adjustment explicitly

```r
interaction <- lambda_hypothesis(
  d, ~group*x + z, target="group:x", adjust=c("group", "x", "z"))
fit <- lambda_model(d, interaction, items, T_q=1/12,
  bootstrap=lambda_bootstrap_control(keep_draws=TRUE, keep_indices=TRUE))
fit$estimates
fit$model$hypothesis
confint(fit, parm="kappa")
```

Use exact term labels from `attr(terms(~group*x + z), "term.labels")`.
Target and adjustment selections must be disjoint. `adjust=character()` means
no selected adjustment terms; the intercept is always included. Optional
`nuisance` adds a separate common adjustment. Unselected terms are labeled
`omitted`. Their variables remain necessary to evaluate the full recipe, and
participate in common case selection. The helper extracts columns from the
full model; it does not reconstruct a reduced formula with potentially different
factor or interaction coding.

For this conditional analysis, lambda_star is the structure left after nuisance
adjustment and kappa is the fraction of that structure captured by the target.
Effects with different adjustment spaces have different denominators. Their
partial kappas need not add to one. No automatic ANOVA sum-of-squares type or
causal interpretation is assigned. General coefficient restrictions L beta=0
use the separate [contrast constructors](CONTRASTS_AND_BASES.md).

## Compare models under one reference

```r
a <- lambda_design(d, ~x, nuisance=~group)
b <- lambda_design(d, ~x+z, nuisance=~group)
comparison <- lambda_compare(d, a, b, items, T_q=1/12)
comparison$delta_F  # F_b - F_a
```

Both objects must declare the same nuisance recipe (contrast recodings may
also demonstrate an identical coefficient null space); a coincidental match of
fitted nuisance spaces on the original sample is insufficient for replay.
Non-NULL call-level coding/nuisance arguments must agree with the object.
A legacy formula/matrix partner in a mixed comparison uses the object's shared
nuisance space. All models use the same selected cases and bootstrap indices.
Explicitly encoded nuisance matrices must agree in values, row and column coding.

Projection uses QR internally, while F_H retains the actual adjusted coding.
Invertible recoding can preserve lambda_H and change the geometry cost of F_H.
The comparison calls the same single-model evaluator twice and shares the
measurement state and thermal terms. It is not a new significance test.

## Bootstrap and scope of this increment

The defaults remain B=5000, 95% ordinary percentile intervals, and seed=1.
Simple formulas and explicit matrices are resampled by row. Sample-dependent
formula transforms are reevaluated after case selection and in each draw by default.
With `encoding="fixed"`, constructors retain the raw encoding learned after common
selection and resample its rows. The [basis guide](CONTRASTS_AND_BASES.md) explains
the distinction between freezing a learned basis and prespecifying its parameters.
This includes the full encoding behind a conditional block, before adjustment.
Storing a design does not freeze a fitted QR basis or convert case bootstrap to
fixed-predictor inference. The encoding policy is visible in the object/result.

Rank loss, lost levels, T=0 and undefined kappa retain their existing diagnostics;
draws are neither replaced nor used to widen endpoints. F_H remains a point
estimate. An ordinary compiled design still removes the overall mean.
Since 0.7.0, [lambda_reference and lambda_paired](REFERENCE_MEANS.md) explicitly
retain a known reference and support the corresponding mean hypotheses.

Runnable examples: `system.file("examples", "model_spaces.R", package="lambdastar")`.
The source defines `lambdastar_model_spaces_example()`; it does not run analysis
when sourced. The [learning path](MODEL_SPACE_WORKFLOW.md)
connects the model definitions to the examples; the separate reference guide records
the new retained-mean semantics. Version 0.8.0's [panel extension](PANEL_MODELS.md)
adds fixed person-by-occasion models, parallel indicators within observations,
and whole-person resampling with distinct fixed-effect copy identities.
