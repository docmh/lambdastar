#' Prepare parallel measurements from a common data frame
#'
#' Only complete, already scored parallel measurements are supported. Missing
#' cases are rejected by default. Complete-case selection uses whole persons,
#' never a different number of items per person. Point estimates and bootstrap
#' draws use the same explicitly selected estimator.
#' @param data A data frame or matrix, ordinarily with one row per person.
#'   A lambda_panel model instead uses one row per person-occasion observation.
#'   Row names identify input observations; panel person labels are declared
#'   separately. Otherwise input positions provide the row identities.
#' @param indicators Character names of at least two numeric measurement columns.
#'   For lambda_measure only, NULL selects all columns.
#' @param measurement_scale Noncentrality denominator: "mean" (default) or
#'   "item". The observed response remains the item mean in either case.
#' @param na.action "fail" or explicit whole-person "complete" selection.
#' @param bootstrap FALSE for point estimates (default), TRUE for 5000-draw
#'   percentile intervals, or settings from [lambda_bootstrap_control()].
#'   See [lambda_inference] for estimator and resampling definitions.
#' @param estimator "paper" (default) uses the manuscript's projection bath
#'   correction. "population" explicitly selects the alternative random-design
#'   population correction; it currently supports case resampling and intercept
#'   adjustment only.
#' @param nuisance NULL, a one-sided formula, or a numeric nuisance matrix on
#'   the input rows. An intercept is included in ordinary centered analyses.
#'   Fixed references use only explicit nuisance columns. Measurements and X_H are
#'   adjusted jointly and measurement df are N minus the nuisance rank.
#' @param nuisance_contrasts Optional contrasts for a nuisance formula.
#' @param numerical_control Controls from lambda_control().
#' @param reference For lambda_measure, NULL uses ordinary centering; a known
#'   scalar/vector retains departures from that reference without an implicit
#'   intercept. For models use a [lambda_reference()] design.
#' @return A lambda_measurement containing the core measurement state plus
#'   included/excluded input positions, original row IDs and effective settings.
#'   estimates, inference and settings identify the selected estimator. The
#'   paper estimates remain in paper; population is present only when selected.
#' @export
lambda_measure <- function(data, indicators = NULL,
    measurement_scale = c("mean", "item"), na.action = c("fail", "complete"),
    bootstrap = FALSE, estimator = c("paper", "population"), nuisance = NULL,
    nuisance_contrasts = NULL, numerical_control = lambda_control(), reference = NULL) {
  boot <- .bootstrap_control(bootstrap)
  data <- .lambda_data(data)
  if (is.null(indicators)) indicators <- names(data)
  prepared <- .prepare_lambda(data, indicators, list(),
    match.arg(measurement_scale), match.arg(na.action), c(boot$strata, boot$cluster),
    nuisance, nuisance_contrasts, numerical_control, match.arg(estimator), origin=reference)
  .attach_bootstrap(prepared$measurement, prepared, boot, NULL)
}

#' Evaluate a complete linear hypothesis model
#'
#' Models describe the hypothesis space after adjustment for the common nuisance
#' design, if supplied separately. Intercepts are removed from the encoded
#' design before ordinary centering. Use ~ 1 for the centered null model. Redundant directions
#' (including a full dummy block from ~ 0 + factor) are rejected, not recoded.
#' A [lambda_design()], [lambda_hypothesis()], [lambda_contrast()] or
#' [lambda_marginal_contrast()] object carries its own coding
#' and nuisance recipe. Non-NULL call settings must agree with the object.
#' A [lambda_reference()] design instead retains its known origin and explicit
#' mean direction: ~1 is the mean model, ~0 its fixed null. No nuisance intercept
#' is then added implicitly. Reference values travel with rows during resampling.
#' A [lambda_panel()] design declares person-by-occasion rows with parallel
#' indicators, person-level complete-case selection and whole-person bootstrap.
#' Its default person adjustment uses within-person demeaning; bootstrap copies
#' receive distinct identities. See lambda_panel for the measurement assumptions.
#' Conditional blocks report total structure and kappa after that adjustment.
#' The manuscript path assumes homogeneous isotropic measurement fluctuation,
#' not Gaussian responses or equal total group variances.
#'
#' One-sided formulas use the item mean as response. A two-sided formula or lm
#' must have a numeric response equal to that mean on the selected rows. Formula
#' variables must be columns of data. Transformations are evaluated after common
#' case selection; transformation-created missing or infinite values are errors.
#' Factors retain their levels. The actual model.matrix coding determines volume.
#'
#' Only unweighted, univariate lm objects with an intercept, saved model frame,
#' no offset, and exactly the selected rows in the same order are accepted. The
#' saved response and model matrix are checked against data. The stored contrasts
#' are retained. Objects fitted on different cases must be refitted explicitly.
#' An explicit matrix supplies encoded directions without an intercept, in data
#' row order; when it has row names they must match data before case selection.
#'
#' F_H always uses the composite mean and its temperature, including when
#' measurement_scale = "item" changes the noncentrality denominator. T_q must
#' therefore be supplied on the composite-mean response scale. Zero temperature
#' returns singular_temperature and NA estimates/F_H, never a boundary estimate.
#' @inheritParams lambda_measure
#' @param model A lambda_design, one- or two-sided formula, suitable lm object, or numeric X_H
#'   matrix with the same input rows as data.
#' @param T_q Positive measurement resolution variance on the mean-response scale.
#' @param contrasts Optional model.matrix contrasts.arg list for a formula.
#'   Matrices and lm objects already specify their coding and reject this argument.
#' @param capacity "numerical" (default) or the explicit "smooth" approximation.
#' @param numerical_control Controls from lambda_control().
#' @param eta_lookup Optional numerical integration lookup from
#'   make_eta_numerical_lookup().
#' @param response_tol Nonnegative pointwise tolerance for checking response and
#'   saved lm design: abs(a-b) <= response_tol * max(1, abs(a), abs(b)).
#' @return A lambda_model with raw/constrained noncentralities, basis, case
#'   provenance, model metadata (formula, columns, term assignment, contrasts,
#'   factor levels, rank), and parsimony containing RSS, costs, volume and F_H.
#' @examples
#' set.seed(17)
#' d <- data.frame(x = rnorm(40), group = factor(rep(c("a", "b"), 20)))
#' d$i1 <- d$x + rnorm(40)
#' d$i2 <- d$x + rnorm(40)
#' lambda_model(d, ~ group + x, c("i1", "i2"), T_q = 1/12)
#' @export
lambda_model <- function(data, model, indicators, T_q,
    measurement_scale = c("mean", "item"), na.action = c("fail", "complete"),
    contrasts = NULL, capacity = c("numerical", "smooth"),
    numerical_control = lambda_control(), eta_lookup = NULL,
    response_tol = sqrt(.Machine$double.eps), bootstrap = FALSE,
    estimator = c("paper", "population"), nuisance = NULL, nuisance_contrasts = NULL) {
  boot <- .panel_bootstrap(bootstrap,list(model))
  settings <- .model_settings(T_q, match.arg(capacity), numerical_control,
                              eta_lookup, response_tol)
  prepared <- .prepare_lambda(.lambda_data(data), indicators,
    list(list(model = model, contrasts = contrasts)),
    match.arg(measurement_scale), match.arg(na.action), c(boot$strata, boot$cluster),
    nuisance, nuisance_contrasts, numerical_control, match.arg(estimator))
  result <- .evaluate_lambda_model(prepared, 1L, settings)
  .attach_bootstrap(result, prepared, boot, settings)
}

#' Compare two models on one common measurement matrix
#'
#' This wrapper prepares the measurement state and case selection once and calls
#' the same single-model evaluator twice. No separate F_H formula is used.
#' Complete-case selection spans indicators and variables from BOTH models.
#' delta_F is F_Hb - F_Ha; negative values mean the second model has lower F_H.
#' This is a descriptive difference, not a p-value or automatic model selection.
#' @inheritParams lambda_model
#' @param model_a,model_b Two explicitly supplied models or codings, as for model.
#' @param contrasts_a,contrasts_b Optional contrasts for the respective formulas.
#' @return A lambda_comparison with model_a, model_b, delta_F, shared measurement,
#'   case provenance and effective settings. At T=0 delta_F is NA.
#' @export
lambda_compare <- function(data, model_a, model_b, indicators, T_q,
    measurement_scale = c("mean", "item"), na.action = c("fail", "complete"),
    contrasts_a = NULL, contrasts_b = NULL,
    capacity = c("numerical", "smooth"), numerical_control = lambda_control(),
    eta_lookup = NULL, response_tol = sqrt(.Machine$double.eps), bootstrap = FALSE,
    estimator = c("paper", "population"), nuisance = NULL, nuisance_contrasts = NULL) {
  boot <- .panel_bootstrap(bootstrap,list(model_a,model_b))
  settings <- .model_settings(T_q, match.arg(capacity), numerical_control,
                              eta_lookup, response_tol)
  prepared <- .prepare_lambda(.lambda_data(data), indicators,
    list(list(model = model_a, contrasts = contrasts_a),
         list(model = model_b, contrasts = contrasts_b)),
    match.arg(measurement_scale), match.arg(na.action), c(boot$strata, boot$cluster),
    nuisance, nuisance_contrasts, numerical_control, match.arg(estimator))
  thermal <- .thermal_terms(prepared$measurement$temperature_mean, settings$T_q,
    settings$capacity, settings$eta_lookup, settings$numerical_control)
  a <- .evaluate_lambda_model(prepared, 1L, settings, thermal)
  b <- .evaluate_lambda_model(prepared, 2L, settings, thermal)
  result <- structure(list(status = prepared$measurement$status,
    message = prepared$measurement$message, model_a = a, model_b = b,
    delta_F = b$parsimony$free_energy - a$parsimony$free_energy,
    measurement = prepared$measurement, cases = prepared$measurement$cases,
    settings = a$settings), class = "lambda_comparison")
  .attach_bootstrap(result, prepared, boot, settings)
}

.lambda_data <- function(data) {
  if ((!is.data.frame(data) && !is.matrix(data)) || nrow(data) < 2L)
    stop("data must be a common data frame or matrix with at least two rows.", call. = FALSE)
  if (is.null(colnames(data)) || anyNA(colnames(data)) ||
      any(!nzchar(colnames(data))) || anyDuplicated(colnames(data)))
    stop("data must have unique, nonempty column names.", call. = FALSE)
  if (!is.null(rownames(data)) && (anyNA(rownames(data)) || anyDuplicated(rownames(data))))
    stop("data must have unique, non-missing row names.", call. = FALSE)
  as.data.frame(data, optional = TRUE)
}

.model_settings <- function(T_q, capacity, control, lookup, response_tol) {
  .scalar(T_q, "T_q", positive = TRUE)
  .scalar(response_tol, "response_tol", nonnegative = TRUE)
  if (capacity == "smooth" && !is.null(lookup))
    stop("A lookup cannot be used with smooth capacity.", call. = FALSE)
  list(T_q = T_q, capacity = capacity, numerical_control = .control(control),
       eta_lookup = lookup, response_tol = response_tol)
}

# Specification inspection happens before case selection; no model transformations
# are evaluated until the common set of raw variables is complete.
.model_spec <- function(model, contrasts, data) {
  if (inherits(model, "lambda_design")) {
    .bound_design(model, data)
    original <- if (model$spec$type %in% c("block", "contrast")) model$spec$full else model$spec
    if (!is.null(contrasts) && !identical(contrasts, original$contrasts))
      stop("contrasts conflict with the lambda_design coding.", call. = FALSE)
    spec <- model$spec
    spec$hypothesis <- model$hypothesis
    spec$design_policy <- model$resampling_policy
    return(spec)
  }
  if (is.matrix(model)) {
    if (!is.numeric(model) || is.complex(model) || nrow(model) != nrow(data))
      stop("Explicit model matrix must be numeric with the same input rows as data.", call. = FALSE)
    .rows(data, model)
    if (!is.null(contrasts)) stop("An explicit matrix already specifies its coding.", call. = FALSE)
    return(list(type = "matrix", x = model, variables = character()))
  }
  fit <- NULL
  if (inherits(model, "lm")) {
    if (!identical(class(model), "lm"))
      stop("Only ordinary univariate lm objects are supported; no GLM or other families.", call. = FALSE)
    fit <- model
    if (is.null(fit$model)) stop("lm requires a saved model frame (model = TRUE).", call. = FALSE)
    if (!is.null(fit$weights) || !is.null(fit$call$weights))
      stop("Weighted lm objects are not supported.", call. = FALSE)
    if (!is.null(fit$offset) || !is.null(fit$call$offset))
      stop("Offsets are not supported.", call. = FALSE)
    if (!is.null(contrasts)) stop("An lm object already specifies its contrasts.", call. = FALSE)
    model <- stats::formula(fit)
    contrasts <- fit$contrasts
  }
  if (!inherits(model, "formula"))
    stop("model must be a formula, an ordinary lm object, or a numeric matrix.", call. = FALSE)
  if (any(c("|", "||") %in% all.names(model)))
    stop("Random-effect formulas are not supported.", call. = FALSE)
  terms <- if (is.null(fit)) stats::terms(model, data = data) else fit$terms
  if (length(attr(terms, "offset"))) stop("Offsets are not supported.", call. = FALSE)
  if (!is.null(fit) && attr(terms, "intercept") != 1L)
    stop("An lm object must include an intercept for the centered-space adapter.", call. = FALSE)
  variables <- all.vars(terms)
  # pi is a deterministic basis constant unless explicitly supplied as a column.
  if (!"pi" %in% names(data) && "pi" %in% variables &&
      identical(get0("pi", envir = environment(model), inherits = TRUE), base::pi))
    variables <- setdiff(variables, "pi")
  absent <- setdiff(variables, names(data))
  if (length(absent)) stop("Formula variables must be columns of data: ",
                           paste(absent, collapse = ", "), call. = FALSE)
  if (!is.null(contrasts) && (!is.list(contrasts) || is.null(names(contrasts)) ||
      any(!nzchar(names(contrasts))) || anyDuplicated(names(contrasts))))
    stop("contrasts must be a named list.", call. = FALSE)
  list(type = if (is.null(fit)) "formula" else "lm", formula = model,
       terms = terms, contrasts = contrasts, fit = fit, variables = variables)
}

.prepare_lambda <- function(data, indicators, candidates, scale, na_action,
                            extra_variables = NULL, nuisance = NULL,
                            nuisance_contrasts = NULL, control = lambda_control(),
                            estimator = "paper", origin = NULL) {
  if (!is.character(indicators) || length(indicators) < 2L || anyNA(indicators) ||
      anyDuplicated(indicators) || !all(indicators %in% names(data)))
    stop("indicators must name at least two distinct measurement columns in data.", call. = FALSE)
  if (!all(vapply(data[indicators], function(x) is.numeric(x) &&
                 !is.complex(x) && is.null(dim(x)), logical(1))))
    stop("All indicators must be numeric vectors, already scored.", call. = FALSE)
  specs <- lapply(candidates, function(c) .model_spec(c$model, c$contrasts, data))
  panel <- .panel_candidates(candidates)
  if(!is.null(panel)&&any(c(panel$person,panel$occasion)%in%indicators))
    stop("Person and occasion identity columns cannot also be measurement indicators.",call.=FALSE)
  if(!is.null(panel)&&estimator!="paper")
    stop("Panel estimation currently requires estimator='paper'; no panel population correction is implemented.",call.=FALSE)
  origins <- lapply(candidates,function(c) if(inherits(c$model,"lambda_design")) c$model$reference else NULL)
  has_origin <- !vapply(origins,is.null,logical(1))
  if(any(has_origin)) {
    if(!all(has_origin)||!all(vapply(origins,function(r)identical(r,origins[[1]]),logical(1))))
      stop("Model comparison requires the same fixed reference; centered and reference designs cannot be mixed.",call.=FALSE)
    origin <- origins[[1]]$value
  }
  origin <- .reference_values(origin,nrow(data),rownames(data))
  if(!is.null(origin)&&estimator!="paper")
    stop("Fixed-reference estimation currently uses estimator = 'paper'; the population estimator has a centered target.",call.=FALSE)
  if (inherits(nuisance, "lambda_design"))
    stop("nuisance must be a one-sided formula or numeric matrix, not a lambda_design.", call. = FALSE)
  objects <- vapply(candidates, function(c) inherits(c$model, "lambda_design"), logical(1))
  nuisance_spec <- if (is.null(nuisance)) NULL else if (any(objects))
    .portable_spec(nuisance, nuisance_contrasts, data) else .model_spec(nuisance, nuisance_contrasts, data)
  if (is.null(nuisance) && !is.null(nuisance_contrasts))
    stop("nuisance_contrasts requires a nuisance formula.", call. = FALSE)
  .check_nuisance_spec(nuisance_spec)
  if (any(objects)) {
    declared <- lapply(candidates[objects], function(c) c$model$B)
    reference <- declared[[1L]]
    if (any(!vapply(declared, function(b) .same_reference(b, reference, control), logical(1))))
      stop("Model comparison requires the same nuisance recipe/reference space.", call. = FALSE)
    if (!is.null(nuisance_spec) && !.same_reference(nuisance_spec, reference, control))
      stop("nuisance conflicts with the lambda_design reference space.", call. = FALSE)
    # A legacy partner in a mixed comparison uses the explicitly shared space.
    nuisance_spec <- reference
  }
  if (estimator == "population" && !is.null(nuisance_spec))
    stop("The population estimator currently supports intercept adjustment only; use estimator = 'paper' with nuisance.", call. = FALSE)
  if (!all(extra_variables %in% names(data)))
    stop("strata must name a column in data, as must cluster when supplied.", call. = FALSE)
  variables <- unique(c(indicators, unlist(lapply(specs, `[[`, "variables")),
                         nuisance_spec$variables, extra_variables, panel$person, panel$occasion))
  selected <- data[variables]
  if (any(vapply(selected, function(x) is.list(x) || !is.null(dim(x)), logical(1))))
    stop("Input variables must be vector columns, not list or matrix columns.", call. = FALSE)
  complete <- stats::complete.cases(selected)
  for (s in c(specs, list(nuisance_spec))) complete <- complete & .spec_complete(s, nrow(data))
  if (any(!complete) && na_action == "fail")
    stop("Missing values in required variables; use na.action = 'complete' explicitly.", call. = FALSE)
  if(!is.null(panel)) complete <- complete &
    !as.character(data[[panel$person]])%in%as.character(data[[panel$person]][!complete])
  rows <- which(complete)
  if (length(rows) < 2L) stop("Fewer than two complete cases remain.", call. = FALSE)
  subset <- data[rows, , drop = FALSE]
  if(!is.null(panel)) {
    .panel_map(subset,panel$person,panel$occasion,panel$person_effects)
    subset <- .panel_data(subset,panel)
  }
  items <- as.matrix(subset[indicators])
  encoding_cache <- new.env(parent = emptyenv())
  encoded_nuisance <- if (is.null(nuisance_spec)) NULL else .encode_model(nuisance_spec,
    list(data = subset, rows = rows, encoding_cache = encoding_cache, numerical_control = control,
      retain_intercept = !is.null(origin),
      measurement = list(composite_response = rowMeans(items))), sqrt(.Machine$double.eps))
  state <- estimate_lambda_star_items(items,
    if (scale == "mean") "composite_mean" else "item", rownames(subset),
    .panel_nuisance(encoded_nuisance$x,subset,panel), control,
    reference=if(is.null(origin)) NULL else origin[rows])
  state$input_rows <- rows
  state$cases <- list(input_n = nrow(data), input_row_ids = rownames(data),
    included = rows, excluded = which(!complete), row_ids = rownames(subset),
    variables = variables, na.action = na_action)
  if(!is.null(panel)) state$cases$selection_unit <- "person"
  state$settings <- list(indicators = indicators, measurement_scale = scale, na.action = na_action,
    estimator = estimator, nuisance = encoded_nuisance$metadata, numerical_control = .control(control))
  if(!is.null(origin)) state$settings$reference <- state$reference
  class(state) <- "lambda_measurement"
  list(data = subset, rows = rows, specs = specs, measurement = state, estimator = estimator,
       nuisance_spec = nuisance_spec, nuisance_x = encoded_nuisance$x,
       nuisance_metadata = encoded_nuisance$metadata, encoding_cache = encoding_cache,
       numerical_control = control, retain_intercept = !is.null(origin), panel=panel)
}

.same_values <- function(a, b, tol) {
  is.numeric(a) && is.numeric(b) && length(a) == length(b) &&
    all(is.finite(a)) && all(is.finite(b)) &&
    all(abs(a - b) <= tol * pmax(1, abs(a), abs(b)))
}

.encode_model <- function(spec, prepared, tol) {
  # This cache exists only inside one prepared dataset/draw, never in results.
  cache <- prepared$encoding_cache
  if (!is.null(cache)) for (entry in cache$entries)
    if (identical(entry$spec, spec) && identical(entry$tol, tol)) return(entry$value)
  value <- .encode_model_uncached(spec, prepared, tol)
  if (!is.null(spec$hypothesis)) {
    value$metadata$hypothesis <- spec$hypothesis
    value$metadata$encoding_policy <- spec$design_policy
  }
  if (!is.null(cache)) cache$entries <- c(cache$entries, list(list(spec = spec, tol = tol, value = value)))
  value
}

.encode_model_uncached <- function(spec, prepared, tol) {
  if (spec$type == "contrast") return(.encode_contrast(spec, prepared, tol))
  if (spec$type == "block") return(.encode_block(spec, prepared, tol))
  if (spec$type == "combined") {
    parts <- lapply(spec$parts, .encode_model, prepared = prepared, tol = tol)
    x <- do.call(cbind, lapply(parts, `[[`, "x"))
    return(list(x = x, metadata = list(type = "combined", columns = colnames(x),
      parts = lapply(parts, `[[`, "metadata"))))
  }
  data <- prepared$data
  if (spec$type == "matrix") {
    rows <- if (identical(spec$row_scope, "selected")) prepared$local_rows else prepared$rows
    if (is.null(rows)) rows <- seq_len(nrow(data))
    x <- spec$x[rows, , drop = FALSE]
    if (is.null(colnames(x)) && ncol(x)) colnames(x) <- paste0("X", seq_len(ncol(x)))
    rownames(x) <- rownames(data)
    if (!is.null(spec$metadata)) return(list(x = x, metadata = spec$metadata))
    return(list(x = x, metadata = list(type = "matrix", formula = NULL,
      terms = NULL, columns = colnames(x), term_assign = NULL, contrasts = NULL,
      factor_levels = list(), intercept_removed = FALSE)))
  }
  fit <- spec$fit
  if (!is.null(fit) && !identical(rownames(fit$model), rownames(data)))
    stop("lm cases/order differ from the common selected rows; refit explicitly.", call. = FALSE)
  frame <- stats::model.frame(spec$terms, data = data, na.action = stats::na.pass,
    drop.unused.levels = FALSE, xlev = if (is.null(fit)) spec$xlevels else fit$xlevels)
  if (nrow(frame) != nrow(data) || !identical(rownames(frame), rownames(data)))
    stop("Formula evaluation changed the common rows.", call. = FALSE)
  terms <- attr(frame, "terms")
  if (attr(terms, "response") == 1L) {
    response <- stats::model.response(frame)
    if (!is.null(dim(response)) || !.same_values(response,
        prepared$measurement$composite_response, tol))
      stop("Model response must equal the composite mean on the common rows.", call. = FALSE)
  }
  if (any(!stats::complete.cases(frame)))
    stop("Formula transformation created missing values on the common rows.", call. = FALSE)
  # Prevent model.matrix from merely warning and ignoring a mistyped contrast.
  factors <- names(frame)[vapply(frame, function(x)
    is.factor(x) || is.character(x) || is.logical(x), logical(1))]
  if (!is.null(spec$contrasts) && !all(names(spec$contrasts) %in% factors))
    stop("Contrasts must refer to factor predictors in the model frame.", call. = FALSE)
  coding <- spec$contrasts
  if (!is.null(spec$contrast_defaults)) {
    for (name in factors) {
      f <- frame[[name]]
      if (!is.factor(f)) f <- factor(f)
      choice <- coding[[name]]
      if (is.null(choice)) choice <- attr(f, "contrasts", exact = TRUE)
      if (is.null(choice)) choice <- spec$contrast_defaults[if (is.ordered(f)) 2L else 1L]
      .check_portable_contrast(choice)
      coding[[name]] <- if (is.matrix(choice)) choice else
        get(choice, envir = asNamespace("stats"))(levels(f))
    }
  }
  mm <- stats::model.matrix(terms, data = frame, contrasts.arg = coding)
  if (!is.null(fit)) {
    saved <- stats::model.matrix(fit)
    if (!identical(dimnames(mm), dimnames(saved)) || !.same_values(mm, saved, tol) ||
        !.same_values(stats::model.response(fit$model),
                      prepared$measurement$composite_response, tol))
      stop("lm response or coding does not match the supplied data.", call. = FALSE)
  }
  assignment <- attr(mm, "assign")
  keep <- if(isTRUE(prepared$retain_intercept)) rep(TRUE,length(assignment)) else assignment != 0L
  x <- mm[, keep, drop = FALSE]
  list(x = x, metadata = list(type = spec$type, formula = spec$formula,
    terms = terms, columns = colnames(x), term_assign = assignment[keep],
    predictor_variables = intersect(all.vars(stats::delete.response(terms)), spec$variables),
    grid_recipe = .grid_recipe(terms, frame),
    encoding = if (is.null(spec$encoding)) "reevaluate" else spec$encoding,
    contrasts = attr(mm, "contrasts"),
    factor_levels = lapply(frame[factors], function(x) if (is.factor(x)) levels(x)
      else if (is.logical(x)) c("FALSE", "TRUE") else levels(factor(x))),
    logical_factors = factors[vapply(frame[factors], is.logical, logical(1))],
    intercept_removed = any(!keep)))
}

# Deterministic H3 handoff: after joint row selection, both public entry points
# call this evaluator. It uses exactly the H1 F_H and partition computations.
.evaluate_lambda_model <- function(prepared, candidate, settings, thermal = NULL) {
  encoded <- .encode_model(prepared$specs[[candidate]], prepared, settings$response_tol)
  state <- prepared$measurement
  basis <- .conditional_basis(encoded$x, state$nuisance, settings$numerical_control)
  result <- .partition_measurement(state, basis)
  parsimony <- .parsimony_state(state$stable_profile, basis, state$temperature_mean,
    settings$T_q, "H", settings$numerical_control, settings$capacity,
    settings$eta_lookup, FALSE, state$nuisance$rank, thermal,
    rss = result$residual_observed)
  encoded$metadata$rank <- parsimony$df_h
  encoded$metadata$row_ids <- state$row_ids
  encoded$metadata$input_rows <- state$input_rows
  result$model <- encoded$metadata
  result$encoded_x <- encoded$x
  result$encoded_nuisance <- prepared$nuisance_x
  result$parsimony <- parsimony
  result$settings <- c(state$settings, settings)
  class(result) <- c("lambda_model", "lambda_measurement")
  result
}
