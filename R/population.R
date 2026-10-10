# Population point estimation, independent of confidence interval construction.
# The finite-person core and numerical F_H are retained separately.
.population_variance <- function(y, groups, stratified) {
  if (!stratified) return(stats::var(y))
  sizes <- lengths(groups)
  if (any(sizes < 2L))
    stop("Population inference requires at least two persons per stratum.", call. = FALSE)
  weights <- sizes / length(y)
  means <- vapply(groups, function(g) mean(y[g]), numeric(1))
  variances <- vapply(groups, function(g) stats::var(y[g]), numeric(1))
  # Unbiased variance of the population mixture with fixed observed weights.
  sum(weights * (1 - (1 - weights) / sizes) * variances) +
    sum(weights * (means - sum(weights * means))^2)
}

.population_h <- function(state, basis, variance_y) {
  p <- basis$rank
  if (p == 0L) return(0)
  if (state$d <= p)
    stop("Population inference requires positive residual degrees of freedom.", call. = FALSE)
  residual <- state$stable_profile - .project(state$stable_profile, basis$Q)
  variance_y - .ssq(residual) / (state$d - p)
}

.component_estimates <- function(context) {
  raw <- context$raw
  total <- max(0, raw[["variance_signal"]])
  scale <- context$d / raw[["temperature"]]
  values <- c(temperature = raw[["temperature"]], lambda_star = scale * total)
  for (i in seq_along(context$prefixes)) {
    h <- min(total, max(0, raw[[i + 2L]]))
    v <- c(lambda_h = scale * h, lambda_perp = scale * (total - h),
           kappa = if (is.finite(total) && total > 0) h / total else NA_real_)
    names(v) <- paste0(context$prefixes[i], names(v))
    values <- c(values, v)
  }
  if (context$singular) values[] <- NA_real_
  values
}

.paper_context <- function(result, prepared) {
  state <- if (inherits(result, "lambda_comparison")) result$measurement else result
  models <- .bootstrap_models(result)
  prefixes <- if (inherits(result, "lambda_comparison")) c("model_a.", "model_b.") else
    if (length(models)) "" else character()
  h <- vapply(models, function(m)
    (m$captured_observed - m$df_h * state$temperature_mean) / state$d, numeric(1))
  list(raw = c(temperature = state$temperature,
    variance_signal = state$signal_energy_raw / state$d,
    stats::setNames(h, paste0(prefixes, "variance_h")[seq_along(h)])),
    d = state$d, prefixes = prefixes,
    variance_y = .ssq(state$stable_profile) / state$d,
    singular = state$status != "ok", model_issues = character())
}

.select_estimator <- function(result, prepared, groups, stratified) {
  paper <- .paper_context(result, prepared)
  context <- if (prepared$estimator == "paper") paper else
    .population_context(result, prepared, groups, stratified)
  estimates <- .component_estimates(context)
  raw_paper <- c(temperature = paper$raw[[1]], lambda_star = paper$d*paper$raw[[2]]/paper$raw[[1]])
  for (i in seq_along(paper$prefixes)) {
    h <- paper$d * paper$raw[[i+2L]] / paper$raw[[1]]
    v <- c(lambda_h = h, lambda_perp = unname(raw_paper["lambda_star"])-h)
    names(v) <- paste0(paper$prefixes[i], names(v))
    raw_paper <- c(raw_paper, v)
  }
  if (paper$singular) raw_paper[-1] <- NA_real_
  result$paper <- list(estimates = .component_estimates(paper), raw_components = paper$raw,
                       raw_estimates = raw_paper)
  result$population <- if (prepared$estimator == "population") list(estimates = estimates,
    raw_components = context$raw, target = "population_at_fixed_N",
    variance_y = context$variance_y, model_issues = context$model_issues) else NULL
  result$estimates <- estimates
  result$inference <- list(estimator = prepared$estimator,
    target = if (prepared$estimator == "paper") "paper_projection" else "population_at_fixed_N",
    N = nrow(prepared$data), scale_factor = context$d,
    projection_assumption = if (prepared$estimator == "paper")
      "fixed_signal_projection_with_measurement_bath_correction" else
      "linear_conditional_mean_homoskedastic_residuals")
  result$settings$estimator <- prepared$estimator
  if(!is.null(prepared$measurement$reference)) {
    result$inference$target <- "fixed_reference_projection"
    result$inference$reference <- prepared$measurement$reference
  }
  if (!is.null(result$model$hypothesis))
    result$inference$hypothesis <- result$model$hypothesis
  result <- .panel_result(result,prepared)
  if (inherits(result, "lambda_comparison")) {
    for (name in c("model_a", "model_b"))
      result[[name]] <- .select_estimator(result[[name]], prepared, groups, stratified)$result
  } else {
    for (key in setdiff(names(estimates), "temperature")) result[[key]] <- unname(estimates[key])
    result$rho_star <- result$lambda_star / (result$lambda_star + context$d)
    raw_total <- context$d * context$raw["variance_signal"] / context$raw["temperature"]
    result$lambda_star_raw <- if (context$singular) NA_real_ else unname(raw_total)
    if (inherits(result, "lambda_model")) {
      raw_h <- context$d * context$raw[3] / context$raw["temperature"]
      result$lambda_h_raw <- if (context$singular) NA_real_ else unname(raw_h)
      result$lambda_perp_raw <- if (context$singular) NA_real_ else unname(raw_total - raw_h)
    }
  }
  list(result = result, context = context)
}

.population_context <- function(result, prepared, groups, stratified) {
  state <- prepared$measurement
  models <- .bootstrap_models(result)
  variance_y <- .population_variance(state$composite_response, groups, stratified)
  issues <- character(length(models))
  h <- vapply(seq_along(models), function(i) {
    value <- .bootstrap_attempt(.population_h(state, models[[i]]$basis, variance_y))
    if (inherits(value, "error")) {
      issues[i] <<- conditionMessage(value)
      NA_real_
    } else value
  }, numeric(1))
  prefixes <- if (inherits(result, "lambda_comparison")) c("model_a.", "model_b.") else
    if (length(h)) "" else character()
  raw <- c(temperature = state$temperature,
    variance_signal = variance_y - state$temperature_mean,
    stats::setNames(h, paste0(prefixes, "variance_h")[seq_along(h)]))
  list(raw = raw, d = state$d, prefixes = prefixes, variance_y = variance_y,
       singular = state$status != "ok", model_issues = issues)
}
