# lambdastar 0.8.1

- Correct README links to files included in the source package and remove the
  private audit reference.
- Exclude internal audit, validation and delivery reports and result tables
  from the CRAN archive. Retain user guides, runnable examples and package tests.
- Remove documentation references to excluded reports. Statistical calculations
  and percentile semantics are unchanged.

# lambdastar 0.8.0

- Add lambda_panel() for fixed person-by-occasion models with parallel indicators,
  explicit observation maps, conditional term blocks and person-level case selection.
- Implement default person-intercept adjustment by group demeaning; reuse the
  existing hypothesis, temperature and numerical F_H evaluator.
- Bootstrap whole persons with distinct copy identities, rebuilding person-specific
  columns and degrees of freedom while preserving declared non-person encoding.
- Deliver the derivation, an executable panel lesson and independent numerical,
  resampling, fixed-profile and compatibility evidence.
- Enforce the exact saturation identity for captured energy. If boot::boot.ci
  returns no percentile endpoints, report quantile_unavailable without a crash
  or substituted limits. Ordinary percentile semantics remain unchanged.

# lambdastar 0.7.0

- Add lambda_reference() for a known scalar/vector origin and explicit retained
  mean direction, through the existing design, projection and F_H evaluator.
  Use d=N without nuisance, literal nuisance columns and uncentered encoded volume.
- Add lambda_paired() for matching after-minus-before parallel measurement pairs.
  Estimate temperature from the differences and require T_q on their mean scale.
- Carry known reference values with case/stratum/cluster draws; require a common
  origin/reference recipe for comparisons. Retain ordinary percentile semantics,
  zero/singular diagnostics and the existing centered/population behavior.
- Include the derivation, a default-bootstrap lesson, independent fixed-profile
  and paired-covariance checks, plus regression evidence for the centered API.

# lambdastar 0.6.1

- Complete the centered model-space learning path with an explicit supported-scope
  table, conditional-reference interpretation, coding comparisons and finite-N
  population-target examples.
- Include the E4 distribution/resource report, compact evidence and two synthetic
  datasets with a runnable ordinary-percentile lesson for regular and weak signal.
- Make the model-space lessons reproducible under caller RNG and contrast options;
  regenerate help and deliver a checked source archive. Numerical algorithms,
  estimator selection, bootstrap defaults and endpoint rules are unchanged.

# lambdastar 0.6.0

- Add lambda_contrast() for named homogeneous coefficient restrictions L beta=0,
  using stable coefficient-space solves and the explicit restricted null model.
- Add lambda_marginal_contrast() for predictor settings and signed averaging
  weights, with learned polynomial, spline and scale prediction coordinates.
- Add encoding="fixed" to design constructors; retain reevaluation as the default.
  Freeze raw encodings only after common selection, and rebuild full/null geometry
  in every case or cluster draw. Preserve declared L coding for F_H.
- Allow comparisons of invertible contrast recodings with a common null reference.
  Diagnose unidentified full models, redundant restrictions and unsupported grid
  recipes. Keep percentile intervals, estimators and invalid-draw rules unchanged.
- Include a runnable contrasts/bases lesson and independent restricted-fit,
  compatibility, resampling and portable-replay validation.

# lambdastar 0.5.0

- Add lambda_design() for portable, row-bound formula/lm/matrix specifications
  and lambda_hypothesis() for target terms conditional on explicit adjustment.
- Preserve declared coding, factor contrasts, common case selection and existing
  bootstrap transformation policies. Reject stale inputs and incompatible model
  comparison reference recipes. Report conditional targets and omitted terms.
- Reuse the measurement nuisance state and projected residual energy for F_H;
  share thermal calculations within comparisons. Keep the projection algorithms,
  estimators, percentile endpoints and undefined-draw rules unchanged.
- Clarify that homogeneous isotropic measurement fluctuation does not require
  Gaussian responses or equal total group variances. The population estimator
  retains its separate assumptions. General coefficient restrictions and new
  reference/paired-mean spaces remain later increments.

# lambdastar 0.4.1

- Add an executable seven-chapter teaching example with known simulated signal,
  plots, prediction checkpoints, inspectable bootstrap draws and editable noise,
  item-count, sample-size and coding experiments. Preserve the compact quickstart.
- Add separate population, nuisance, school-cluster and zero-case examples.
  Teaching runs preserve RNG state and save results only to a new/empty folder.
- Include a worked default run and its figures in the installed documentation.
  Statistical estimators and numerical algorithms are unchanged from 0.4.0.

# lambdastar 0.4.0

- License the package under Apache 2.0 and use mike.hammes@mikehammes.name
  as the author and maintainer contact.
- Make `estimator="paper"` the default throughout points, returned lambda/kappa
  fields, summaries and bootstrap draws. Keep `estimator="population"` explicit;
  `estimates` and `inference` identify the selected calculation. This changes the
  default bootstrap statistic from 0.3.1; choose population to retain that method.
- Add shared nuisance formula/matrix adjustment to measurement, model and
  comparison APIs, using N-rank(B) and refitting nuisance effects on every draw.
- Add whole-cluster resampling within strata, variable draw sizes and optional
  row-index lists. Add explicit type-7 percentile interpolation for paper replay.
- Retain ordinary percentile intervals, reduced contrast coding, undefined-value
  diagnostics and caller RNG restoration. No CI endpoint alterations.
- Add independent conditional-moment, cluster and estimator checks, plus public
  API reproduction of the OECD manuscript analysis and random-design validation.

# lambdastar 0.3.1

- Preserve explicit factor contrast matrices, including fewer than G-1 columns,
  when rebuilding bootstrap models. Formula transformations, saved lm models
  and comparisons use the supplied coding without adding contrast directions.
- Add regression tests against explicit design matrices and independent OLS for
  reduced contrasts, transformations and interactions.

# lambdastar 0.3.0

- All bootstrap confidence intervals use direct percentiles of the parameter
  draws, computed by boot::boot.ci(type="perc"). Defaults: B=5000, 95%, seed=1.
- Population point estimation is independent of confidence intervals and is
  available also with bootstrap=FALSE. The same estimator runs on every draw.
- Undefined kappa, singular temperature and model failures are counted explicitly.
  Missing draws suppress the affected regular interval; valid zero and one
  values stay in their distributions.
- Optional conditional_quantiles are a separate descriptive output. They never
  substitute for a regular interval.
- ci.method accepts only "percentile", also for temperature intervals.
- Existing finite-person numerical estimates and F_H calculations are preserved.

# lambdastar 0.1.0

- Locally installable measurement core, model adapters, numerical F_H and
  common-data comparison wrapper, with executable examples and validation.

# lambdastar 0.0.3

- Person-bootstrap infrastructure, RNG restoration and diagnostics.

# lambdastar 0.0.2

- Model adapters, common-case comparison and S3 summaries.

# lambdastar 0.0.1

- Measurement/projection core, numerical capacity, temperature integrals,
  encoded geometry and parsimony functional.
