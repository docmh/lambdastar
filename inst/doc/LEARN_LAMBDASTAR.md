# Learn Lambda Star by changing a working example

This lesson starts with noisy measurements of a known stable signal. You will
estimate the background fluctuation, find out how much structure a model
captures, inspect its bootstrap distribution, and compare explicitly coded
models using the numerical parsimony functional F_H.

You need only R and an installed `lambdastar`. The plots use base R graphics.
The compact `quickstart.R` remains available as a short API introduction.
You can also [read the completed default run](../tutorial/WALKTHROUGH.md)
with its seven figures before running the lesson yourself.

## Start the lesson

Run this in R or RStudio:

```r
library(lambdastar)
source(system.file("examples", "learn_lambdastar.R", package="lambdastar"))
lesson <- lambdastar_tutorial(
  output_dir="my-lambda-star-lesson",
  pause=TRUE,
  show_plots=TRUE
)
```

Sourcing the file only defines the teaching functions. Calling
`lambdastar_tutorial()` starts the lesson. At each chapter, make a prediction
and press Enter to reveal the calculation and plot. In noninteractive R,
the lesson advances automatically. In RStudio, plots appear in the plot pane.
The script also writes every figure to the output folder.

Open `my-lambda-star-lesson/WALKTHROUGH.md` in your editor's Markdown preview.
It contains the explanations, your actual computed tables, and seven figures.
The same folder holds the synthetic data, CSV tables, all bootstrap draws,
`lesson.rds`, and session information. A nonempty output folder is never
overwritten: choose a new name for each experiment.

The default uses **5,000 draws and 95% percentile intervals**. For a short
practice run, explicitly use `B=199`; that gives less precise empirical tail
quantiles. Computational time depends on the machine. This is one worked
example, not a simulation study of interval coverage.

```r
# Run without pausing or opening a graphics window:
lesson <- lambdastar_tutorial(
  output_dir="my-lambda-star-batch", pause=FALSE, show_plots=FALSE
)

# Reopen saved results without repeating the bootstrap:
lesson <- readRDS("my-lambda-star-lesson/lesson.rds")
d <- lesson$example$data
indicators <- lesson$example$indicators
T_q <- lesson$settings$T_q
```

## 1. See what repeats and what fluctuates

The synthetic outcome is wellbeing in arbitrary continuous score units.
It is not an empirical dataset or a validated questionnaire. Its stable signal is

```text
S = -0.9 * workload + 0.7 * resources + 0.5 * hidden_profile
item_j = S + independent measurement error_j
```

There is also an `unused` predictor which contributes no stable signal.
All four underlying profiles are centered, orthogonal and have sample variance
one. This deliberately controlled design makes the true contributions easy to
inspect; real predictors need not be orthogonal. Once these profiles are
generated, they are treated as the fixed signal/design in this example.

**Predict:** Will the mean of four items track S more closely than one item?

Figure 1 shows the individual items, their mean and the known signal for twelve
people, then compares the observed mean with the known signal for all people.
The known signal is available only in this simulated example. The package
receives the observed items and model variables, not the hidden profile.

## 2. Estimate T and lambda_star

```r
measurement <- lambda_measure(d, indicators)
measurement$estimates
lesson$oracle[c("temperature", "lambda_star")]
```

**Predict:** Should every difference between observed mean scores count as
stable structure?

The example's item-error SD is 1.2. With K=4 independent parallel items, the
generating background variance of their mean is T=1.2²/4=0.36. Lambda Star
estimates that background from item disagreement.

On the default mean-score scale, and with only an intercept removed:

```text
d = N - 1
observed energy = sum of squared centered item means
raw stable energy = observed energy - d*T
raw lambda_star = raw stable energy / T
```

Figure 2 makes that subtraction visible. Estimated stable energy is constrained
to be nonnegative. T is a variance; lambda_star is dimensionless total
noncentrality. **Lambda_star is not a percentage.** With the default settings,
the known stable variance is 1.55 and true lambda_star=159*1.55/0.36=684.5833.
An estimate from one noisy dataset need not equal that value.

## 3. Locate the structure captured by a model

```r
workload_model <- lambda_model(d, ~workload, indicators, T_q)
full_model <- lambda_model(d, ~workload + resources, indicators, T_q)
workload_model$estimates
full_model$estimates
```

**Predict:** Does the second model capture all the stable structure?

The actual encoded predictor columns define X_H. The model projects the
observed profile into that space and corrects the captured and remaining
energies for their measurement-noise budgets. For p hypothesis directions:

```text
raw lambda_H    = captured observed energy / T - p
raw lambda_perp = residual observed energy / T - (d-p)
kappa_H        = estimated captured signal / estimated total signal
```

The returned partition constrains captured signal between zero and total.
At positive estimated total signal, lambda_H + lambda_perp = lambda_star.
If estimated total signal is zero, kappa is undefined.

Workload's known capture fraction is 0.81/1.55, about **52.3%**. Adding
resources raises it to 1.30/1.55, about **83.9%**. The hidden profile accounts
for the remainder. It is stable structure omitted by the hypothesis, not item
measurement noise. Figure 3 shows the estimated captured/remaining partition
for the null, workload, full, and extra-predictor models.
Estimated lambda_H can decrease after adding a direction because the extra
noise correction can exceed the gain in observed captured energy.

## 4. Look inside the percentile bootstrap

```r
ix <- lesson$bootstrap$bootstrap$indices[1, ]
d[ix[1:10], c("id", "workload", "resources", indicators)]
draws <- lesson$bootstrap$bootstrap$draws
lesson$bootstrap$bootstrap$ci
confint(lesson$bootstrap, parm="kappa")

# Inspectability: recreate the first draw through the public API.
first_draw <- lambda_model(d[ix, ], ~workload+resources, indicators, T_q)
rbind(saved=draws[1, ], recalculated=first_draw$estimates)

# Check availability before calculating a percentile of the full distribution.
alpha <- (1-lesson$settings$conf.level)/2
if (all(is.finite(draws[, "kappa"]))) {
  quantile(draws[, "kappa"], c(alpha, 1-alpha), type=7)
}
```

**Predict:** When one person appears twice, which measurements and predictors
must accompany that person?

Every draw resamples N complete rows with replacement and recalculates the
same selected paper estimator. The tutorial explicitly uses type-7
interpolation, so the final `quantile()` call reproduces the endpoints. The
package's default `boot` interpolation is also an ordinary percentile method.
Figure 4 shows the actual endpoints, without widening or replacing them.

Changing the confidence level reuses the retained draws:

```r
confint(lesson$bootstrap, parm="kappa", level=.90)
```

The diagnostics show unavailable statistics explicitly. A zero lambda_star
draw stays in its own distribution; its kappa is undefined. Such draws are not
silently discarded to manufacture a regular interval. Histograms label the
finite-draw count and CI status. A 95% percentile interval is the central 95%
of the bootstrap statistic's distribution; nominal confidence does not mean
that every chosen data-generating design has exactly 95% repeated-sample coverage.

## 5. Understand model comparison through F_H

```r
useful <- lambda_compare(d, ~workload, ~workload+resources, indicators, T_q)
extra <- lambda_compare(d, ~workload+resources,
                       ~workload+resources+unused, indicators, T_q)
useful$delta_F
extra$delta_F
lesson$tables$parsimony
lesson$tables$increments
```

**Predict:** Can a model fit the observed scores better and still have a larger F_H?

```text
F_H = residual energy + dimension cost + geometry cost
    = RSS + p*Xi_eta(T,T_q) - Xi_a(T,T_q)*log(V_H)
delta_F = F_H(model_b) - F_H(model_a)
```

Xi_eta is numerically evaluated by the package. V_H is the normalized volume
of the declared model columns. Figure 5 separates the change in residual energy
from the extra costs of adding the unused predictor. Negative delta_F means
model b has the smaller functional; positive delta_F favors model a by that
criterion. It is not a p-value. A noisy dataset can sometimes favor an
irrelevant extra predictor, so the walkthrough interprets the actual result.

T_q is the chosen resolution variance on the mean-score representation. The
default 1/12 is an explicit illustrative assumption for this continuous-data
lesson; it is not estimated from the data and does not imply an actual Likert
scale. Models within a comparison must share the same data, measurement scale
and resolution. F_H and delta_F remain point estimates.

## 6. Predict, change one argument, rerun

Figure 6 already calculates three controlled series of experiments. Known
values and estimates are distinguished by line style and marker shape.
Use the table for exact values:

```r
lesson$tables$experiments

# Run these separately, after making your prediction.
noisier <- lambdastar_tutorial("noise-2", noise_sd=2)
more_items <- lambdastar_tutorial("eight-items", items=8)
more_people <- lambdastar_tutorial("320-people", n=320)
finer_resolution <- lambdastar_tutorial("finer-resolution", T_q=1/48)
```

- More measurement noise raises T and lowers true lambda_star.
- More independent parallel items lower the mean's T and raise true lambda_star.
- More people raise true lambda_star at fixed stable variance and T:
  **lambda_star=(N-1)*Var(S)/T** in this controlled design.
- The true capture fraction stays 1.30/1.55 in these experiments, although its
  noisy estimate can change.
- Changing only T_q leaves T and the noncentrality partition unchanged while
  changing the resolution-dependent terms in F_H.

These are individual generated datasets, not averaged Monte Carlo curves.
At fixed N, the seed preserves the profiles and underlying item-error draws
when noise or K changes. Changing N generates a new standardized design.
Compare models within a measurement condition; these plots do not justify
comparing F_H across different data matrices or measurement conditions.

## 7. See the role of coding

```r
d$resources_tilted <- d$resources + 2*d$workload
coding <- lambda_compare(d, ~workload+resources,
                        ~workload+resources_tilted, indicators, T_q)
lesson$tables$coding
```

**Predict:** Which quantities change when the hypothesis space stays the same?

The transformation is invertible: both models span the same space. Their
projections, lambda_H and kappa agree. The normalized column volume changes,
so the geometry cost and F_H can change. Figure 7 shows both encoded designs.
The lesson is to state and justify the model coding, not to search for a coding
that produces a desired criterion value. The package uses QR for projection
without substituting that computational basis for the declared geometry.

## Continue with advanced examples

The [model-space learning path](MODEL_SPACE_WORKFLOW.md) adds an explicit
supported-scope table, conditional effects, coefficient/grid contrasts, fixed
and reevaluated bases, and an executable weak-signal interval example. It also
distinguishes this tutorial's fixed-profile truth from a population target at
fixed N and the expectation of the selected paper statistic.

[ADVANCED_EXAMPLES.md](ADVANCED_EXAMPLES.md) and the separately executable
`advanced_lambdastar.R` introduce the population estimator, site adjustment,
whole-school resampling, zero signal, and singular temperature.
