# Advanced Lambda Star examples

Complete [the main lesson](LEARN_LAMBDASTAR.md) first. These examples introduce
different estimation targets and sampling designs separately from that lesson.

```r
library(lambdastar)
source(system.file("examples", "learn_lambdastar.R", package="lambdastar"))
source(system.file("examples", "advanced_lambdastar.R", package="lambdastar"))
advanced <- lambdastar_advanced_examples()  # B=5000, seed=29
# For a short practice run: lambdastar_advanced_examples(B=199)
```

The function prints its results and returns all data and fitted objects. It
does not write files. Save the return value with `saveRDS()` if desired.
Both teaching functions preserve the caller's random-number state.

## A. Choose the population estimator explicitly

The first example generates independent random predictors and an independent
unobserved stable component. Its known random-population targets are based on
Var(S)=1.55 and Var(Pi_H S)=1.30, at N=160 and T=0.36.

```r
rbind(paper=advanced$paper$estimates,
      population=advanced$population$estimates,
      population_target=advanced$population_truth)
advanced$population$inference
confint(advanced$population, parm="kappa")
```

**Predict:** Does giving both estimators a percentile CI make their points equal?

The paper path corrects captured energy for p measurement-noise directions.
The population path instead uses the residual-variance correction for the
linear random-design population projection. Their targets/assumptions differ;
their finite-sample estimates can differ. Each path uses its own selected
statistic consistently on the original data and every draw. This example uses
the same person resamples for both methods so resampling randomness does not
hide the distinction. The sample's realized latent signal and the population
variance target should not be confused.

The population option assumes a linear conditional mean and homogeneous
residual variance. It currently supports intercept adjustment and person
resampling, optionally within fixed strata. It does not support additional
nuisance or cluster combinations. The main lesson and the following examples
use `estimator="paper"`.

## B. Define what remains after nuisance adjustment

The second dataset has 192 students in 24 schools across four sites. School
sizes vary between 6, 8 and 10. Stable site and school effects are shared across
the parallel items; independent item errors are added separately.

```r
d <- advanced$cluster_data
items <- advanced$indicators
adjusted <- lambda_model(d, ~workload+resources, items, T_q=1/12,
                        nuisance=~site)
adjusted$nuisance$rank
adjusted$d
```

**Predict:** Should site-adjusted kappa retain the same total-structure denominator
as the unadjusted analysis?

The intercept and site contrasts span the common nuisance space B. Both the
measurements and hypothesis columns are projected off that space. The nuisance
rank is four, so d=192-4=188. The resulting partition concerns structure remaining
after site adjustment. Its denominator is conditional on that adjustment, so
the unadjusted and adjusted kappas describe different quantities.

For a comparison, provide one data frame and one common nuisance space:

```r
lambda_compare(d, ~workload, ~workload+resources, items, T_q=1/12,
               nuisance=~site)
```

## C. Resample schools within sites

```r
fit <- advanced$clustered
fit$bootstrap$config
confint(fit)
fit$bootstrap$sampling$cluster_counts
range(fit$bootstrap$sampling$draw_sizes)
fit$bootstrap$indices[[1]]
```

**Predict:** Will resampling six whole schools per site always yield 192 students?

School sizes differ, so the number of students can vary. Every selected school
brings all its student rows and their measurements. A school selected twice
brings every student twice. There are six school draws within each site; the
site effects are refitted and residual measurement degrees of freedom are
recomputed in every draw.

This sampling scheme reflects the clustered design and the shared school
effects. It is not selected because it gives wider or narrower intervals.
The output remains an ordinary percentile interval with explicit diagnostics.

## D. Inspect the two distinct zero cases

```r
advanced$zero_signal$estimates
advanced$zero_signal$status
advanced$zero_temperature$estimates
advanced$zero_temperature$status
advanced$zero_temperature$temperature
```

Opposite item profiles have positive disagreement but an identically zero
mean: estimated lambda_star=0 and kappa is undefined. The example exposes an
estimation boundary; these opposite items are not a valid parallel scale.
Identical nonconstant item profiles have estimated T=0: that is the singular
temperature case, with no background fluctuation for regular estimation.
It is not reported as an infinite noncentrality or as perfect localization.

For complete API definitions, use `help("lambda_inference")`.
