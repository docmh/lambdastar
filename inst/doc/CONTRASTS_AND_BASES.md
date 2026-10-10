# Linear contrasts and reproducible bases

Since 0.6.0, lambdastar provides `lambda_contrast()` for named homogeneous coefficient
restrictions, `lambda_marginal_contrast()` for explicitly weighted predictor
settings, and `encoding="fixed"` for the design constructors. They all use the
existing measurement, projection, F_H and percentile-bootstrap calculations.
The [measurement axiom and conditional interpretation](MODEL_SPACES.md) still apply.

## Play through the example

```r
library(lambdastar)
source(system.file("examples", "contrasts_and_bases.R", package="lambdastar"))
lesson <- lambdastar_contrasts_example()  # B=5000, 95%, seed=1 for bootstrap
lesson$coefficient_names
lesson$joint
confint(lesson$joint)
```

Use `B=199` for a short practice run. Sourcing defines the function without
running it; the function restores the caller's random state and writes no files.
The synthetic data have stable group effects, moderation and nonlinear structure,
plus independent parallel measurements of one homogeneous bath. T_q=1/12 is an
explicit example resolution on the composite-mean scale, not an estimated value.

## 1. State the coefficient hypothesis

For a full model `~group*x+z` with reference group a, the b-a effect at x=1 is
`groupb + groupb:x`. The analogous c-a effect can be tested jointly:

```r
d <- lesson$data; items <- lesson$items
L <- rbind(b_minus_a=c(groupb=1,groupc=0,`groupb:x`=1,`groupc:x`=0),
           c_minus_a=c(groupb=0,groupc=1,`groupb:x`=0,`groupc:x`=1))
h <- lambda_contrast(d, ~group*x+z, L)
fit <- lambda_model(d, h, items, T_q=1/12)
fit$model$L
fit$encoded_x
fit$encoded_nuisance
fit$basis$X
```

L's column names refer to the full model's actual non-intercept coefficient
names. Omitted columns receive zeros. A named vector specifies one restriction;
matrix row names and order specify the effect coding. Unknown columns, dependent
restrictions or an unidentified full conditional design give errors. A full-model
`lambda_design` can be supplied instead, retaining its coding and nuisance policy.
Formula/lm and explicit matrix sources are supported. User-supplied nuisance
columns are added before deriving the restricted null model.

## 2. See what is projected

In coefficient space the implementation obtains a right inverse D and a null
basis N, with `L D = I` and `L N = 0`, through scaled QR solves. With full encoded
columns X and existing nuisance B (including the intercept):

```text
restricted null model: span(B, X N)
raw hypothesis coding: X D
evaluated hypothesis:  X D residualized against the restricted null model
```

The solver uses coefficient-sized matrices, without an observation-by-observation
projector or an inverse of normal equations. `fit$model` retains L, D, N and the
full-model metadata.

Lambda_star now describes structure remaining **after the restricted null
model**; kappa divides captured structure by that conditional total. It is not
the fraction of an unadjusted whole-model lambda_star. Different null models
therefore yield different denominators, whose kappas need not sum to one.

Change the ordered L rows with an invertible transformation and inspect
`lesson$comparison`: lambda_H is unchanged, while F_H can change because the
actual contrast geometry changes. `lambda_compare()` permits this recoding only
when the full recipe, original nuisance and coefficient null space agree for
every draw. For grid contrasts, this check requires identical grids and equal
weight row spaces. Different restricted null models are rejected. Two separate
fits do not establish comparison compatibility. Delta_F remains F_b-F_a.

## 3. Specify settings and averaging weights directly

```r
grid <- data.frame(group=c("a","b","c"),x=1,z=0,
                   row.names=c("a1","b1","c1"))
w <- rbind(b_minus_a=c(a1=-1,b1=1,c1=0),
           c_minus_a=c(a1=-1,b1=0,c1=1))
h_grid <- lambda_marginal_contrast(d, ~group*x+z, grid, w)
```

Each grid row contains all required predictors at an explicitly chosen setting.
Each weights row must sum to zero and have named columns matching every grid row.
The helper computes `L = weights %*% X_grid`. It does not infer a reference
population, covariate mean, equal-cell weights or ANOVA sum-of-squares type.
The example's `averaged` result uses an explicit 25%/75% average over x=-1 and x=2.
Changing those weights changes the question.

Polynomial, spline and scale grid values use the coordinates learned on the
selected analysis data, not parameters fitted on the small grid. R's
[prediction recipe mechanism](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/makepredictcall.html)
retains polynomial coefficients and spline knots; lambdastar also records scale
centers and divisors. Unsupported adaptive transformations such as `I(x-mean(x))`
require an explicit coefficient L or a supported prediction recipe. Missing grid
predictors, unknown factor levels and zero/dependent resulting restrictions are errors.

## 4. Choose the basis policy explicitly

| Policy | Original sample | Bootstrap draw |
|---|---|---|
| `encoding="reevaluate"` (default) | Evaluate after common case selection | Reevaluate sample-dependent terms, then rebuild full/null geometry |
| `encoding="fixed"` | Evaluate after common case selection and retain the learned encoding | Resample raw encoded rows; rebuild full/null geometry |
| Explicit numeric X | Use supplied encoding | Resample its rows; rebuild geometry |

This argument is available on `lambda_design()`, `lambda_hypothesis()` and both
contrast constructors. A design supplied to a contrast constructor already owns
its policy. `fixed` applies also to its explicit nuisance recipe. Neither policy
slices a previously residualized X or fitted orthogonal projection basis.
Persons and predictors continue to be resampled together, including whole
clusters when requested. A frozen encoding is not fixed-predictor inference.

`lesson$fixed` and `lesson$adaptive` have identical original estimates and draw
indices, but different draws for the named polynomial coefficient hypothesis.
Inspect `lesson$fixed$model$full_model$grid_recipe$terms` to see the learned
polynomial coefficients. For a whole model, these metadata are in
`fit$model$grid_recipe$terms`. Scale parameters and spline interior/boundary knots
are retained there too. Equivalent full polynomial spans may yield identical
capture under the two policies; coefficient hypotheses and F_H still depend on
the declared coordinates.

A saved design stores a deferred recipe bound to its input rows. A `fixed`
encoding is learned anew on each call's common complete cases and then frozen
within that call's bootstrap. It is not a claim that basis parameters were
prespecified independently of the data. For that, supply explicit parameters:

```r
spline <- lambda_design(d, ~splines::ns(x,knots=c(-.5,.5),Boundary.knots=c(-4,4)))
fourier <- lambda_design(d, ~sin(2*pi*x)+cos(2*pi*x))  # period 1
piecewise <- lambda_design(d, ~x+I(pmax(x-.5,0)))      # knot .5
```

Use literal parameters in portable recipes; custom transforms or caller-local
constants should be evaluated into a numeric matrix explicitly. No knot search
or smoothing selection is added. Saved designs replay under changed global
factor-contrast defaults because the original coding is retained.

## 5. Inspect intervals and scope

`lesson$joint$bootstrap$draws` contains the actual B estimates;
`lesson$joint$bootstrap$indices` records the resampled persons. All CIs remain
ordinary equal-tail percentile intervals. Any invalid draw retains its diagnostic
and suppresses the affected regular CI; it is never replaced. Valid zeros remain
in the distribution. T=0 is singular and kappa is undefined at lambda_star=0.
F_H is still a point estimate, calculated separately for each model.

E3 contrasts use `estimator="paper"`: their null-model adjustment does not extend
the optional population estimator's scope. Nonzero restrictions and intercept,
one-sample or paired-reference hypotheses use the separate
[E6a reference API](REFERENCE_MEANS.md); the centered E3 constructors do not
silently switch their reference space.
The [model-space learning path](MODEL_SPACE_WORKFLOW.md) explains the
population targets and interval interpretation with runnable examples.
