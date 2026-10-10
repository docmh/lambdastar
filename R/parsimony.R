#' Evaluate the parsimony functional for one encoded hypothesis
#'
#' The response and supplied temperature must describe the same measurement
#' representation. For a mean of parallel items use its composite_mean
#' temperature, not the single-item temperature. The encoded X determines
#' volume; its orthonormal basis is used only for projection.
#' @param response Finite numeric response vector.
#' @param x Explicit hypothesis matrix on response rows. Omit an intercept for
#'   ordinary centered analyses; retain it explicitly to evaluate a reference mean.
#' @param temperature Nonnegative analysis-scale background variance.
#' @param resolution Positive T_q on the same response scale.
#' @param model Descriptive label for this single model.
#' @param control Numerical controls.
#' @param eta_method Numerical capacity (default) or the labeled smooth approximation.
#' @param eta_lookup Optional lookup for numerical integration.
#' @param boundary Explicitly request the mathematical zero-temperature limit.
#' @param nuisance Optional numeric nuisance design on the response rows.
#'   Response and hypothesis directions are projected off this common space;
#'   an intercept is implicit only for ordinary centered analyses.
#' @param reference NULL uses the centered reference. A known scalar or numeric
#'   vector retains departures from that origin without an implicit intercept.
#' @return List containing status, model, n, df_h, df_perp, rss, log_volume,
#'   volume, temperature, resolution, xi_eta, xi_a, dimension_cost,
#'   geometry_cost, free_energy, numerical diagnostics and encoded basis.
#'   At T=0 regular estimation returns NA for the thermal terms and F_H.
#'   boundary=TRUE instead labels the result boundary_only; it is not inference.
#' @export
thermodynamic_model_state <- function(response, x, temperature, resolution,
    model = "H", control = lambda_control(),
    eta_method = c("numerical", "smooth"), eta_lookup = NULL, boundary = FALSE,
    nuisance = NULL, reference = NULL) {
  response <- .vector(response, "response")
  .scalar(temperature, "temperature", nonnegative = TRUE)
  .scalar(resolution, "resolution", positive = TRUE)
  control <- .control(control)
  eta_method <- match.arg(eta_method)
  if (!is.character(model) || length(model) != 1 || is.na(model))
    stop("model must be one character label.", call. = FALSE)
  if (!is.logical(boundary) || length(boundary) != 1 || is.na(boundary))
    stop("boundary must be TRUE or FALSE.", call. = FALSE)
  if (boundary && temperature > 0)
    stop("boundary=TRUE is reserved for the T=0 mathematical limit.", call. = FALSE)
  origin <- .reference_values(reference,length(response),names(response))
  space <- .nuisance_space(nuisance, length(response), names(response), control,center=is.null(origin))
  basis <- .conditional_basis(x, space, control)
  if (nrow(basis$X) != length(response))
    stop("response and x have incompatible row counts.", call. = FALSE)
  if (!is.null(names(response)) && !is.null(rownames(basis$X)) &&
      !identical(names(response), rownames(basis$X)))
    stop("response and x have different row identities or order.", call. = FALSE)
  z <- .residualize_nuisance(if(is.null(origin)) response else .finite(response-origin), space)
  .parsimony_state(z, basis, temperature, resolution, model, control, eta_method,
    eta_lookup, boundary, space$rank)
}

.thermal_terms <- function(temperature, resolution, eta_method, eta_lookup, control,
                           boundary = FALSE) {
  singular <- temperature == 0
  integral <- list(value = NA_real_, abs_error = NA_real_, method = "not_estimated")
  geometry <- NA_real_
  if (!singular || boundary) {
    if (eta_method == "numerical") {
      integral <- xi_eta_numerical(temperature, resolution, eta_lookup, control, details = TRUE)
    } else {
      if (!is.null(eta_lookup)) stop("A lookup cannot be used with smooth capacity.", call. = FALSE)
      integral <- list(value = xi_eta_app(temperature, resolution),
                       abs_error = NA_real_, method = "smooth_approximation")
    }
    geometry <- xi_a(temperature, resolution)
  }
  list(integral = integral, geometry = geometry)
}

.parsimony_state <- function(z, basis, temperature, resolution, model, control,
    eta_method, eta_lookup, boundary, nuisance_rank, thermal = NULL, rss = NULL) {
  if (is.null(rss)) rss <- .ssq(z - .project(z, basis$Q))
  singular <- temperature == 0
  status <- if (singular) {
    if (boundary) "boundary_only" else "singular_temperature"
  } else "ok"
  if (is.null(thermal)) thermal <- .thermal_terms(temperature, resolution,
    eta_method, eta_lookup, control, boundary)
  integral <- thermal$integral; geometry <- thermal$geometry
  dimension_cost <- basis$rank * integral$value
  geometry_cost <- -geometry * basis$log_volume
  free_energy <- rss + dimension_cost + geometry_cost
  if (!singular || boundary) .finite(c(dimension_cost, geometry_cost, free_energy))
  list(status = status, message = if (singular) .singular_message() else NULL,
       model = model, n = length(z), df_h = basis$rank,
       df_perp = length(z) - nuisance_rank - basis$rank, nuisance_rank = nuisance_rank, rss = rss,
       log_volume = basis$log_volume, volume = basis$volume,
       temperature = temperature, resolution = resolution,
       xi_eta = integral$value, xi_a = geometry, eta_method = eta_method,
       dimension_cost = dimension_cost, geometry_cost = geometry_cost,
       free_energy = free_energy,
       numerical = list(method = integral$method, integral_abs_error = integral$abs_error,
                        control = control), basis = basis)
}
