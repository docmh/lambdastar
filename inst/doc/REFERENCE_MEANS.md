# Reference means and paired differences (E6a)

Since 0.8.0, the separate [panel guide](PANEL_MODELS.md) covers fixed
person-by-occasion observations with parallel indicators and person bootstrap.

This extends the manuscript's centered projection identities to a **prespecified
origin and an explicitly retained mean direction**. The existing centered API
is unchanged. The source checked for this derivation is `symmetry_localization.pdf`,
21 September 2026 snapshot: pp. 9 (axiom and centering), 13–14 (projected energy),
17–18 (parallel-measurement moments). This reference-space extension
is derived here; it is not presented as an already implemented manuscript example.

## Reference and retained space

Let r be a known reference vector, B an explicitly specified nuisance matrix,
and A=I-P_B. With no B, A=I. Write d=rank(A)=N-rank(B), U=A(S_absolute-r),
Z_k=A(Y_k-r), and z=mean_k Z_k. The conditional error covariance on this
space is sigma² A, with independent parallel errors across k. Then

```text
T_hat = sum_k ||Z_k-z||² / [d (K-1) K]
E[T_hat | U] = sigma²/K = T
E[||z||²-d T_hat | U] = ||U||²
E[||P_H z||²-p T_hat | U] = ||P_H U||²
```

P_H projects onto the actual A X_H and p is its rank. Division by estimated T
and the existing nonnegative constraints give the same paper-style statistic
as before; the raw energy identities do not imply unbiased constrained ratios.
The parsimony functional remains

```text
F_H = ||z-P_H z||² + p Xi_eta(T,T_q) - Xi_a(T,T_q) log(V_H).
```

V_H is the volume of column-normalized **A X_H**, without additional centering.
All computations use thin decompositions, never an N-by-N projector.
The old centered case is r=0, B=1, d=N-1. A fixed-reference analysis has no
implicit intercept in B; a formula's intercept is retained in X_H.

For a scalar reference mu0 and the mean direction X_H=1, d=N and p=1:

```text
lambda_star,reference = [N (mu-mu0)² + ||S_centered||²] / T
lambda_H,mean         = N (mu-mu0)² / T
kappa_H               = N (mu-mu0)² / [N (mu-mu0)² + ||S_centered||²]
lambda_H,raw          = N (mean(Y_bar)-mu0)² / T_hat - 1
```

These are squared departures. The observed signed mean departure is also
reported, separately from lambda_H. Total structure now includes the departure
from the known reference, so it differs substantively from centered lambda_star.
At estimated total zero, kappa remains undefined; T=0 remains singular.

The reference is prespecified, not estimated from these observations and then
treated as known. The optional population estimator is not extended by E6a.
For population-oriented person bootstrap at fixed N, the reference-energy
comparator is N E[(S_absolute-r)²]/T; a mean-direction population comparator
is N(E[S_absolute-r])²/T. These need not equal the expected finite-sample paper
statistic. No coverage guarantee or alternate endpoint rule is introduced.

## Paired observations and resolution

Rows are persons. `before[j]` and `after[j]` are matching parallel measurement
columns, for j=1,...,K. Form D_ij=Y_after,ij-Y_before,ij **before** estimating
temperature. On a stacked observation vector the transformation is C=[-I,I],
represented by its two coefficients and column pairs, not a dense matrix.
With scalar within-occasion covariances and common within-pair covariance gamma:

```text
Cov(error_difference,k) = (sigma_after² + sigma_before² - 2 gamma) I
T_difference = (sigma_after² + sigma_before² - 2 gamma) / K.
```

Different j pairs must have independent difference errors. Homogeneous
cross-occasion correlation within a matching pair is allowed; temperature is
estimated from the differences and no independent-occasion shortcut is used.
For a general row transform, C Sigma C' must be scalar on the retained space.
Under Sigma=T I, C C'=a I gives transformed T=a T. General transforms or panel
covariances are not automatically justified by this special paired case.

T_q is supplied explicitly on the **mean difference response scale**. Two
different resolution models illustrate why the package does not double it:

- For independent quantization errors with variances q_before and q_after,
  the difference error variance is q_before+q_after; with covariance q_cov,
  it is q_before+q_after-2 q_cov.
- For differences represented on one lattice of spacing h, the declared
  bin-resolution variance is h²/12, even when both original scores used that
  spacing. This is a lattice convention, not the variance of two accumulated
  quantization errors.

Choose the model appropriate to the actual scored difference and its composite
representation. An exact fixed reference shift does not change a resolution
variance. A response rescaling by c requires T and T_q scaled by c².
All resolution variances above refer to the composite response being analyzed.
The illustrative T_q=1/12 below is a declared score-scale resolution, not an
automatically inferred resolution for a multi-item Likert composite.

## Usage and common-reference comparisons

The installed lesson provides synthetic parallel measurements with known
baseline/change profiles and correlated within-pair errors. Source the file,
then run it; B=199 is an explicit short practice budget:

```r
source(system.file("examples","reference_means.R",package="lambdastar"))
lesson <- lambdastar_reference_example()  # B=5000, 95%, seed=1
lesson$one_sample$estimates              # baseline mean against 50
lesson$paired$estimates                  # mean change against zero
lesson$paired_reference$estimates        # mean change against .5
lesson$truth                            # fixed-profile oracles, not fitted values
```

```r
library(lambdastar)
d <- lesson$data
items <- lesson$before
h <- lambda_reference(d, ~1, value=50)
fit <- lambda_model(d, h, items, T_q=1/12, bootstrap=TRUE)
fit$reference
fit$mean_departure
confint(fit)

# Fixed null value versus a free mean, under one common origin.
h0 <- lambda_reference(d, ~0, value=50)
cmp <- lambda_compare(d, h0, h, items, T_q=1/12)
cmp$delta_F

# Positive differences mean after minus before; reference may be nonzero.
paired <- lambda_paired(d, before=paste0("pre",1:4), after=paste0("post",1:4),
  T_q=1/12, value=0, bootstrap=TRUE)
paired$paired
paired$mean_departure
```

`lambda_reference()` returns a row-bound `lambda_design`. Formula intercepts
are retained: `~1` is the mean direction and `~0` is the zero-dimensional null.
Numeric matrices supply exactly their encoded columns, including a constant
column when intended. The value can be scalar or a known vector on the input
rows. Nuisance formulas are literal: `~0+x` adjusts for x alone; `~x` includes
an intercept and therefore removes the mean direction. Confounded target
directions produce the existing rank diagnostics.

Every model in a comparison must use the same reference vector and nuisance
recipe. Mixing a centered design with a fixed-reference design is rejected.
Case selection is common, reference values travel with resampled persons, and
paired columns are never resampled independently. Case/stratified/cluster
percentile settings and diagnostics are retained. F_H remains a point estimate.

## Affine hypotheses and scope

An affine hypothesis L beta=c can be represented explicitly: choose D and N0
with L D=I and L N0=0, then r=X D c, nuisance B=X N0 and target X_H=X D.
Use `lambda_reference(data, X_H, value=r, nuisance=B)` with those independently
declared encodings. This retains the actual target coding and tests departures
from that affine null. Different affine origins are not automatically comparable.
The E3 coefficient/grid constructors retain their centered L beta=0 contract;
they do not silently transform this reference specification.

E6a supports a known mean, paired mean differences and explicit affine designs
under the scalar-bath condition. It still needs at least two parallel items or
replicates per response representation. A single before/after column pair alone
cannot identify a separate stable change profile and its measurement bath.
Fixed within-person/occasion models with parallel indicators are implemented
separately by [lambda_panel](PANEL_MODELS.md). Estimated external-reference
uncertainty and general covariance/mixed models still require their own extensions.
