# Fixed models on person-by-occasion observations (E6b)

Version 0.8.0 adds `lambda_panel()`: a constructor for the same X_H/nuisance
representation used by `lambda_model()` and `lambda_compare()`. It adds an
explicit observation map and the resampling rule needed for repeated observations.

Each row is **one person at one occasion**. At least two scored parallel items
or replicates occupy columns within that row. The stable profile may change
between occasions; those changes are part of the structure to be studied.
Occasions are not used as if they were repeated measures of an unchanged profile.
Person IDs are globally unique identities and each person–occasion pair occurs
once. Different persons may have different sets or numbers of observed occasions.

## 1. Retained space and measurement temperature

Let R be the number of observation rows, P the number of persons, and K the
parallel measurements per row. These distinguish person count from observation
count; the result's existing `n`/`N` fields count observation rows.

For default person adjustment, let G contain person-indicator columns,
D=I−P_G, and B_tilde=D B for any additional nuisance columns. Then

```text
A = (I - P_B_tilde) D
d = R - P - rank(B_tilde)
Z_k = A Y_k,  z = mean_k Z_k,  H = A X_H
T_hat = sum_k ||Z_k-z||² / [d (K-1) K]
```

The required measurement assumption is independent parallel errors across k
with `A Cov(error_k) A = sigma² A` on the retained observation space. Consequently,

```text
E[T_hat | S] = sigma²/K = T
E[||z||² - d T_hat | S] = ||A S||²
E[||P_H z||² - p T_hat | S] = ||P_H A S||²
```

The paper-style ratios, constraints and kappa use these corrected energies.
Numerical F_H retains the same RSS, capacity and normalized encoded-volume terms.
T_q is supplied on the composite-mean scale of the observed response. No response
transformation occurs in this constructor, so no automatic resolution change is
introduced. As elsewhere, `measurement_scale="item"` changes the noncentrality
denominator; F_H continues to use the composite-mean response and temperature.

Person demeaning implements D directly using group sums. It constructs neither
an R-by-R projector nor a dummy matrix for default person-intercept adjustment.
Only additional nuisance and hypothesis columns require thin QR decompositions.
Explicit person effects in a target formula do require their encoded columns.

Stable non-Gaussian responses and unequal stable between-group dispersions remain
allowed. General serial measurement covariance is not made isotropic by cluster
bootstrap. Some covariance disappears under the declared nuisance projection:
an error shared by all occasions of a person, separately for each replicate,
lies in G and is removed exactly. This specific case is checked in the validation;
it is not permission to treat arbitrary longitudinal covariance as scalar.

With `person_effects="none"`, ordinary global centering and explicit nuisance
apply instead. Person identity can then be a target term. These are different
retained spaces and have different total-structure/kappa interpretations.

## 2. Work through the installed example

The synthetic lesson has 40 persons, four occasions and four parallel items per
occasion. Two treatment arms differ in their stable time slopes. Person intercepts,
person-specific slopes and curvature also contribute stable structure. Its known
mean-scale bath temperature is `.7^2/4 = .1225`.

```r
library(lambdastar)
source(system.file("examples", "panel_models.R", package="lambdastar"))
lesson <- lambdastar_panel_example()   # two B=5000, 95%, seed-1 analyses
d <- lesson$data
items <- lesson$indicators
head(d)

# Average linear occasion trend after removing person intercepts.
h_time <- lambda_panel(d, ~time, person="id", occasion="time")
lesson$trend$estimates
lesson$trend$panel
confint(lesson$trend)

# Arm difference in slopes, conditional on the common time trend.
h_interaction <- lambda_panel(d, ~arm*time, person="id", occasion="time",
                              target="arm:time", adjust="time")
lesson$interaction$estimates
confint(lesson$interaction)
```

The arm main effect is constant within each person, so person adjustment removes
it. It should not also be included as a target column in that adjusted space.
The interaction uses the full model's declared coding and explicitly selects
the target and adjustment terms; no automatic ANOVA type is inferred.
The two fits above have different nuisance spaces, so their kappas do not form
an additive partition and their F_H values are not a common-space comparison.

```r
# Default model: categorical occasions, allowing any mean occasion pattern.
h_occasions <- lambda_panel(d, person="id", occasion="time")
occasion_fit <- lambda_model(d, h_occasions, items, T_q=1/12)

# Common person adjustment permits a time-versus-curvature comparison.
h_curve <- lambda_panel(d, ~time + I(time^2), person="id", occasion="time")
comparison <- lambda_compare(d, h_time, h_curve, items, T_q=1/12)
comparison$delta_F  # F_curve - F_time

# Person-specific slopes beyond the common trend.
h_slopes <- lambda_panel(d, ~id*time, person="id", occasion="time",
                         target="id:time", adjust="time")
slopes <- lambda_model(d, h_slopes, items, T_q=1/12)

# Alternatively include person intercepts among the target model columns.
h_full <- lambda_panel(d, ~id+time, person="id", occasion="time",
                       person_effects="none")
full <- lambda_model(d, h_full, items, T_q=1/12)
```

`1/12` is this lesson's illustrative scored-response resolution convention,
not a value inferred from continuous observations. Choose T_q for your actual
measurement scale. A numeric `time` in an explicit formula describes a linear
trend; the default constructor instead uses `factor(time)` for a categorical
occasion effect. Arbitrary fixed bases and explicit non-person X_H columns can
be used through the same constructor.

## 3. Whole-person bootstrap and inspectable copies

`bootstrap=TRUE` resolves to B=5000 whole-person cluster draws. A configurable
call can preserve observed treatment-arm person counts:

```r
cfg <- lambda_bootstrap_control(
  B=5000, conf.level=.95, seed=1,
  resampling="cluster", cluster="id", strata="arm",
  keep_draws=TRUE, keep_indices=TRUE)
fit <- lambda_model(d, h_interaction, items, T_q=1/12, bootstrap=cfg)
confint(fit)
fit$bootstrap$draws[1:5, ]
fit$bootstrap$failures

# Reconstruct an actual draw, including separate sampled copies of a person.
j <- 1L
ix <- fit$bootstrap$indices[[j]]
draw <- d[ix, ]
draw$id <- fit$bootstrap$sampling$person_instances[[j]]
h_draw <- lambda_panel(draw, ~arm*time, person="id", occasion="time",
                       target="arm:time", adjust="time")
replay <- lambda_model(draw, h_draw, items, T_q=1/12)
all.equal(unname(replay$estimates), unname(fit$bootstrap$draws[j, ]))
```

Persons and their full selected trajectories are sampled together. If one person
is drawn twice, the copies have separate fixed intercept/slope identities. Simply
reusing the original ID would merge their nuisance columns, change rank and
miscalculate T. `indices` refers to original input row positions;
`person_instances` gives the draw-local identity for every row of each copy.
The same point estimator is recomputed on each draw, including temperature and
retained degrees of freedom. Actual intervals remain ordinary percentiles.

P is held fixed, or fixed within each declared stratum. For unbalanced panels,
R and d can vary with the sampled trajectory lengths. These are empirical
whole-person sampling intervals; this does not implement the separate population
point-estimator correction. It is not occasion-row bootstrap or a conditional
fixed-X procedure. Strata must be constant within each person; explicit controls
requesting case sampling or a different cluster column are rejected.

Person-dependent formulas require `encoding="reevaluate"` because their columns
belong to sampled copies. The original numerical factor coding is retained;
sampled copies occupy draw-local factor slots. Fixed non-person bases retain
their declared raw encoding. Matrix columns are treated as fixed row-bound
directions; encode person intercepts/slopes through the person declaration and
formula rather than a matrix of original-ID dummies.
Person labels serve as fixed-effect identities, not substantive numerical
predictors. Inference about contrasts between particular named individuals would
require a different resampling target from the whole-person procedure here.

## 4. Cases, comparisons and scope

- Missing data fail by default. With `na.action="complete"`, any incomplete
  required row excludes the entire person, using variables from both models in
  a comparison. Unobserved occasions already absent from an unbalanced panel
  are not imputed or silently converted into missing rows.
- At least two persons with at least two observed occasions each, and at least
  two parallel indicators per observation, are required. A panel of one score
  per occasion does not identify its fluctuation through this estimator.
- Comparisons use one common data matrix, observation map, person-adjustment
  policy and nuisance specification. The existing evaluator runs twice.
- Lost rank/levels, zero temperature, zero estimated total signal and undefined
  kappa retain their diagnostics. Draws are retained and never replaced. For a
  saturated retained hypothesis, captured structure equals total structure
  exactly; numerical projection residue is not treated as uncaptured energy.
- At the two-occasion limit, the mean change agrees with E6a's parallel paired
  differences: T_difference=2 T_panel and the lambda values agree. With matching
  resolution scaling, F_difference=2 F_panel. This is a mathematical change of
  response representation, not a license to compare F values on different scales.
- `target`/`adjust` provide full-model term blocks. E3's coefficient/grid contrast
  helpers do not propagate panel maps and reject a panel source; declare the
  wanted within-person contrast as explicit non-person X_H instead.

General mixed models, estimated covariance metrics/GLS and uncertainty of an
estimated external reference remain separate extensions. The
[reference guide](REFERENCE_MEANS.md) covers known
means and paired differences separately.
