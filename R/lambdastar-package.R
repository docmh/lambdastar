#' Lambda Star: measurement, linear hypotheses and percentile intervals
#'
#' A standalone R package for parallel measurements, nuisance-adjusted linear
#' models, the numerical parsimony functional and ordinary case/cluster percentile
#' intervals. Use [lambda_measure()], [lambda_model()] and [lambda_compare()].
#' The manuscript estimator is the default; an alternative population estimator
#' is explicitly selectable. Points and bootstrap draws use the same selected
#' statistic. See [lambda_inference].
#' [lambda_design()] binds a reusable encoded model to its input rows;
#' [lambda_hypothesis()] selects target terms conditional on explicit adjustment
#' terms. [lambda_contrast()] encodes named homogeneous coefficient restrictions;
#' [lambda_marginal_contrast()] uses explicit predictor settings and weights.
#' Design constructors support fixed or reevaluated bootstrap basis recipes.
#' The manuscript path assumes stable structure plus zero-mean homogeneous,
#' isotropic background fluctuation with finite second moments. Neither Gaussian
#' responses nor equal total group variances are required. The measurement
#' assumption is used, not established by fitting. The optional population
#' estimator has additional assumptions described in [lambda_inference].
#' An installed example is available at
#' `system.file("examples", "quickstart.R", package = "lambdastar")`.
#' A guided lesson with plots, known simulated signal and editable experiments
#' is at `system.file("doc", "LEARN_LAMBDASTAR.md", package = "lambdastar")`.
#' Source `system.file("examples", "learn_lambdastar.R", package = "lambdastar")`
#' and call lambdastar_tutorial() to create a personal worked example.
#' The centered fixed-model learning path and supported-scope table are at
#' `system.file("doc", "MODEL_SPACE_WORKFLOW.md", package = "lambdastar")`.
#' It links reusable designs, conditional effects, coefficient/grid contrasts,
#' basis policies and the distinction between fixed-profile and population targets.
#' Source `system.file("examples", "intervals_and_targets.R", package="lambdastar")`
#' and call lambdastar_intervals_example() for regular and weak-signal percentile
#' examples (B=5000, 95%, seed=1), with actual draws and undefined-value diagnostics.
#' [lambda_reference()] retains a known reference origin and mean direction;
#' [lambda_paired()] evaluates after-minus-before parallel measurement pairs.
#' The derivation and resolution conventions are in installed `doc/REFERENCE_MEANS.md`.
#' [lambda_panel()] encodes fixed person-by-occasion designs with parallel items,
#' explicit person adjustment and whole-person bootstrap with distinct copies.
#' Its derivation and runnable lesson are in installed `doc/PANEL_MODELS.md`.
#' General covariance/mixed models and uncertainty in an estimated reference
#' remain separate extensions. Ordinary centered calls keep their existing meaning.
#' T_q is supplied on the composite-mean response scale. Measurements must be
#' already scored; the package does not perform imputation or validate a scale.
#' @name lambdastar
#' @aliases lambdastar-package
"_PACKAGE"
