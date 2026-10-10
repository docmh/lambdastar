# From a scientific question to an encoded model

Version 0.6.1 completed this centered fixed-model learning path. Version 0.7.0
adds the separate [reference-mean and paired extension](REFERENCE_MEANS.md), and
0.8.0 adds [fixed person-by-occasion models](PANEL_MODELS.md). Start with the
[measurement tutorial](LEARN_LAMBDASTAR.md) if lambda_star, T and F_H are new to
you. Here you will choose a hypothesis, inspect its X_H and reference space,
compare explicit codings, and examine actual percentile draws when signal is weak.

All these models use the same calculation under **stable structure plus zero-mean,
homogeneous isotropic measurement fluctuation with finite second moments**.
Parallel measurements identify the bath under the independent-replicate
assumption. Total responses can be non-Gaussian or have unequal group variances.
Fitting a model uses the measurement axiom; it does not establish it from the data.

## 1. Choose the question before the constructor

| Scientific question / model | Encoding or constructor | Supported meaning |
|---|---|---|
| Two independent group means | `~group` with two levels | Centered group-mean space; no conventional t/Welch statistic is calculated |
| One-way or factorial ANOVA | `~group`, `~a*b` | Whole declared model space |
| Regression, ANCOVA, moderation | `~x+z`, `~group+x`, `~group*x+z` | Fixed linear predictor columns and interactions |
| Conditional term or block | `lambda_hypothesis()` | Selected target after explicitly selected adjustment terms |
| Reduced planned contrasts | Numeric X_H or explicit factor contrast matrix | Only the supplied directions, including fewer than G−1 |
| Homogeneous coefficient restrictions | `lambda_contrast(..., L)` | Named L beta=0; restricted null model becomes the reference |
| Simple effects / marginal contrasts | `lambda_marginal_contrast()` | Explicit predictor grid and signed weights; no automatic averaging population |
| Polynomial, spline, Fourier, piecewise basis | Formula recipe or numeric X_H | Fixed linear coefficients; declared knots/periods or explicit basis-learning policy |
| Fixed blocking/site adjustment | `nuisance=~site` | Common centered reference; lambda_star and kappa are conditional on it |
| Clustered sampling | `resampling="cluster", cluster=...` | Resample whole clusters; does not change the measurement covariance model |
| Alternative population estimation | `estimator="population"` | Intercept-only adjustment and case resampling; additional linear-mean/homogeneous-residual assumptions; see `help("lambda_inference")` |
| One-sample and paired mean hypotheses | `lambda_reference()`, `lambda_paired()` | Known reference retained explicitly; T_q supplied on the difference scale; see the E6a guide |
| General affine L beta=c | Explicit reference vector, X_H and nuisance matrix | Algebraic construction in the E6a guide; E3 constructors retain their centered contract |
| Fixed person×occasion models | `lambda_panel()` | Parallel indicators within observations, explicit person adjustment and whole-person bootstrap; see [panel guide](PANEL_MODELS.md) |
| General serial covariance structures | Separate covariance-metric extension | Cluster resampling does not make an anisotropic measurement bath isotropic |
| GLS/WLS, mixed or generalized linear models, penalized/adaptive model search | Separate mathematical extensions | Not supported by the current fixed-projector calculation |

Choose the basis from the substantive hypothesis. X_H retains its column names,
units and coding. A QR basis is an internal projection device; it does not replace
the geometry of the declared model in F_H.

## 2. Inspect a whole model and a conditional interaction

```r
library(lambdastar)
source(system.file("examples", "model_spaces.R", package="lambdastar"))
m <- lambdastar_model_spaces_example()  # one B=5000 bootstrap; B=199 for practice
d <- m$data; items <- m$items
m$full$estimates
m$conditional$estimates
m$conditional$model$hypothesis
m$conditional$encoded_x
m$conditional$encoded_nuisance
m$conditional$basis$X
```

**Predict:** Does interaction kappa divide by the whole model's total structure?

The full design is `~group*x+z`. The conditional design selects `group:x`,
adjusting for `group`, `x` and `z`. Its lambda_star describes the structure
remaining after that adjustment; kappa is the fraction captured from that
remaining structure. Changing the adjustment changes the denominator. Kappas
from differently adjusted effects need not add to one.

`encoded_x` contains raw hypothesis columns, `encoded_nuisance` the reference
coding, and `basis$X` the adjusted hypothesis matrix. Recipes are realized after
common case selection. Saved designs retain their input-row binding:

```r
path <- tempfile(fileext=".rds")
saveRDS(m$design, path)
replayed <- lambda_model(d, readRDS(path), items, T_q=1/12)
stopifnot(isTRUE(all.equal(replayed$estimates, m$full$estimates)))
unlink(path)
```

Changed predictors or row order require a new design. See the
[design guide](MODEL_SPACES.md) for formula/matrix inputs and shared cases.

## 3. Separate dispersion from a mean hypothesis

```r
v <- m$unequal_variance
tapply(v$signal, v$data$group, mean)
tapply(v$signal, v$data$group, var)
v$oracle
v$fit$estimates
v$welch
```

**Predict:** Can there be substantial total structure but no true group-mean capture?

This example centers its known stable profile within each group and gives the
groups different stable dispersions. Its fixed-profile oracle has positive
lambda_star but lambda_H=kappa=0. Its four item errors share one homogeneous
bath. A particular noisy estimate need not be exactly zero. The separately
returned Welch result asks its conventional mean question; none of its standard
errors, t values or degrees of freedom enter lambda-star computations.

Equal *population* means are a different setup. A fresh finite sample can have
different realized latent group means and therefore positive latent projection.
Keep three quantities distinct:

| Quantity | What is held fixed / averaged |
|---|---|
| Fixed-profile captured energy | Project this particular latent S into this particular X_H |
| Population target at fixed N | `(N−1) * Var(population projection of S) / T` |
| Expected paper estimate | Average the actual corrected ratio statistic over the declared sampling process |

In E4's equal-population-means cell the population lambda_H target is zero;
the expected raw paper statistic is 5.42695, and the simulated mean is
5.48542 with MC SE .90388. Incidental finite-sample projection and division by
estimated T explain this distinction. The fixed-profile test agrees with its
zero captured-energy expectation. The package keeps `estimator="paper"`
selected unless you explicitly choose the supported population alternative.

## 4. Specify effects and basis policies

```r
source(system.file("examples", "contrasts_and_bases.R", package="lambdastar"))
c <- lambdastar_contrasts_example(B=199)  # explicit practice budget
c$coefficient_names
c$joint$model$L
c$joint$estimates
c$via_grid$estimates
c$averaged$estimates
c$fixed$model$full_model$grid_recipe$terms
```

**Predict:** Will freezing polynomial coordinates freeze the sampled persons?

The joint example encodes b−a and c−a at x=1. The grid version states those
settings directly; its values agree. The averaged version states a 25%/75%
weighting over two x settings. The [contrast guide](CONTRASTS_AND_BASES.md)
explains the coefficient null space and supplies editable L/grid/weight code.

`encoding="fixed"` freezes the raw encoding learned after common selection for
that call. `"reevaluate"` rebuilds sample-dependent recipes in each draw. Both
resample persons with all predictors and items and rebuild the reference and
projection. Inspect their shared indices and differing coefficient draws:

```r
identical(c$fixed$bootstrap$indices, c$adaptive$bootstrap$indices)
head(c$fixed$bootstrap$draws - c$adaptive$bootstrap$draws)
```

For prespecified coordinates, use literal spline knots/boundaries, a declared
Fourier period or an explicit numeric matrix. A frozen learned basis is different
from prespecifying it independently of the sample.

## 5. Compare the same data and reference, with actual coding

```r
a <- lambda_design(d, ~x, nuisance=~group)
b <- lambda_design(d, ~x+z, nuisance=~group)
comparison <- lambda_compare(d, a, b, items, T_q=1/12)
components <- function(f) unlist(f$parsimony[
  c("rss", "dimension_cost", "geometry_cost", "free_energy")])
rbind(a=components(comparison$model_a), b=components(comparison$model_b))
comparison$delta_F  # F_b - F_a; negative means b has the smaller functional
c$comparison$delta_F  # same contrast space, changed contrast coding
```

F_H is RSS plus dimension and geometry costs. Each model gets its own calculation;
the wrapper shares the data, measurement state and thermal work. Invertible
recoding preserves capture but can change encoded volume and F_H. Comparisons
require compatible reference recipes, the same selected cases, response scale
and T_q. Different conditional null models do not become comparable just because
their fits use the same number of rows. F_H and delta_F remain point estimates.

## 6. Inspect available and unavailable percentile intervals

```r
source(system.file("examples", "intervals_and_targets.R", package="lambdastar"))
i <- lambdastar_intervals_example()  # B=5000, level=.95, seed=1
i$intervals[, c("scenario", "parameter", "population_target", "estimate",
               "lower", "upper", "n_valid", "n_undefined", "status")]
regular <- i$fits$baseline
weak <- i$fits$weak
draws <- weak$bootstrap$draws
table(total_zero=draws[, "lambda_star"]==0,
      kappa_undefined=is.na(draws[, "kappa"]))
weak$bootstrap$diagnostics
```

**Predict:** What happens to a ratio when its estimated denominator is zero?

These are the two prespecified synthetic E4 datasets, shipped with generator
provenance. The lesson uses fresh default-seed draws; E4 used its stored per-cell
seeds. The weak dataset has estimated lambda_star=0 and undefined kappa. Valid
zero lambda-star and lambda_H draws remain in their respective distributions;
undefined kappa draws retain their reasons. Its regular kappa CI is unavailable,
while the other parameter intervals are calculated independently.

The default lesson run gives the following rounded values:

| Dataset | Parameter | Estimate | 95% percentile interval | Undefined draws / 5000 |
|---|---|---:|---|---:|
| Baseline | lambda_star | 627.7241 | [454.5871, 850.1655] | 0 |
| Baseline | lambda_H | 300.8166 | [198.0474, 437.0789] | 0 |
| Baseline | kappa | .479218 | [.363757, .603042] | 0 |
| Weak | lambda_star | 0 | [0, 39.8796] | 0 |
| Weak | lambda_H | 0 | [0, 8.6081] | 0 |
| Weak | kappa | undefined | unavailable | 2967 |

These describe these datasets with seed 1; they are not coverage estimates.

Recreate a draw using the public API and inspect a percentile without dropping NAs:

```r
ix <- regular$bootstrap$indices[1, ]
one <- lambda_model(i$data$baseline[ix, ], ~group, i$items, T_q=1/12)
rbind(stored=regular$bootstrap$draws[1, ], replayed=one$estimates)
v <- regular$bootstrap$draws[, "kappa"]
stopifnot(all(is.finite(v)))
percentile <- boot::boot.ci(structure(list(t0=regular$kappa,
  t=matrix(v, ncol=1), R=length(v), sim="ordinary"), class="boot"),
  type="perc", index=1)$percent[1, 4:5]
rbind(independent=percentile, package=as.numeric(confint(regular, parm="kappa")))
confint(regular, parm="kappa", level=.90)  # reuses the actual stored draws
```

`B`, `conf.level` and `seed` are editable lesson arguments. B=199 uses an explicit
`min_valid=199` practice override. Changing the level or B does not change the
point statistic. No simulation-truth value is supplied to the estimator.

E4's diagnostic M=100/B=199 study observed 91–98% kappa coverage in informative
in-scope cells. The weak cell had 0/100 regular kappa intervals and 85% lambda_H
coverage. These are finite simulation results, with MC uncertainty, not a uniform
95% coverage guarantee or grounds to alter the percentile procedure.
At true lambda_star=0 there is no population kappa coverage target.

## 7. Work at the measured scale

The E4 N=30,000 benchmark's median point-comparison time fell from .856 s in
0.4.1 to .632 s in 0.6.0; process peak RSS rose from 313.5 to 324.1 MiB. These
describe that local workload, not a scaling guarantee. Designs retain inspectable
raw encodings, avoid N-by-N projectors and share measurement/thermal work.
Choose `keep_draws` and `keep_indices` deliberately when retaining large results.

Version 0.6.1 packaged these lessons and compact E4 evidence. Version 0.7.0 adds
E6a's [known-reference calculations](REFERENCE_MEANS.md) while preserving this
centered path. Version 0.8.0 implements E6b's fixed panel observation mapping and
person resampling. NLSS integration remains the later H6 handoff.
