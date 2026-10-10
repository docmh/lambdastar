.coordinates <- function(temperature, resolution) {
  for (name in c("temperature", "resolution")) {
    v <- get(name)
    if (!is.numeric(v) || is.complex(v) || !length(v) || any(!is.finite(v)))
      stop(name, " must be a nonempty finite real numeric vector.", call. = FALSE)
  }
  n <- max(length(temperature), length(resolution))
  if (!length(temperature) %in% c(1L, n) || !length(resolution) %in% c(1L, n))
    stop("temperature and resolution have incompatible lengths.", call. = FALSE)
  temperature <- rep_len(temperature, n)
  resolution <- rep_len(resolution, n)
  if (any(temperature < 0) || any(resolution <= 0))
    stop("temperature must be nonnegative and resolution positive.", call. = FALSE)
  ratio <- temperature / resolution
  if (any(!is.finite(ratio)) || any(temperature > 0 & ratio == 0))
    stop("Temperature/resolution ratio is outside numerical range.", call. = FALSE)
  list(temperature = temperature, resolution = resolution, ratio = ratio)
}

.log_sum <- function(v) {
  a <- max(v)
  if (a == -Inf) return(-Inf)
  a + log(sum(exp(v - a)))
}

# Upper integral bounds for both omitted positive-bin tails. Bounds start
# after the mode of c * exp(-beta*c), so the summands are decreasing.
.gibbs <- function(log_beta, control) {
  beta <- exp(log_beta)
  limit <- max(8, ceiling(0.5 + sqrt(-log(control$tail_tol) / beta)))
  repeat {
    if (!is.finite(limit) || limit > control$max_support)
      stop("Numerical eta support exceeds max_support; adjust the control or range.",
           call. = FALSE)
    costs <- (seq_len(limit) - 0.5)^2
    lw <- -beta * costs
    log_weights <- .log_sum(lw)
    log_numerator <- log(2) + .log_sum(lw + log(costs))
    log_z <- log1p(exp(log(2) + log_weights))
    u <- limit - 0.5
    log_i0 <- 0.5 * (log(pi) - log_beta) +
      stats::pnorm(sqrt(2 * beta) * u, lower.tail = FALSE, log.p = TRUE)
    log_i2 <- .log_sum(c(log(u) - beta * u^2 - log(2) - log_beta,
                         log_i0 - log(2) - log_beta))
    log_tail_z <- log(2) + log_i0 - log_z
    log_tail_m <- log(2) + log_i2 - log_numerator
    if (max(log_tail_z, log_tail_m) <= log(control$tail_tol)) break
    limit <- limit * 2
  }
  list(log_moment = log_numerator - log_z, log_z = log_z,
       beta = beta, support = limit,
       tail_relative_bound = exp(max(log_tail_z, log_tail_m)))
}

.eta_scalar <- function(ratio, control) {
  if (ratio == 0) return(c(eta = 0, beta = NA_real_, relative_moment_error = 0,
                           support = 0, tail_relative_bound = 0))
  # Keep the target in log space, including below the normal floating range.
  log_target <- log(ratio) - log(12)
  initial <- if (ratio < 1) log(4 * (log(6) - log(ratio))) else log(6) - log(ratio)
  objective <- function(b) .gibbs(b, control)$log_moment - log_target
  lower <- upper <- initial
  for (i in seq_len(128L)) {
    if (objective(lower) >= 0) break
    lower <- lower - log(2)
  }
  for (i in seq_len(128L)) {
    if (objective(upper) <= 0) break
    upper <- upper + log(2)
  }
  if (lower == upper) {
    root <- initial
  } else {
    root <- stats::uniroot(objective, c(lower, upper),
      tol = max(1e-15, control$capacity_tol / (1 + abs(log_target))),
      maxiter = 200L, check.conv = TRUE)$root
  }
  g <- .gibbs(root, control)
  error <- g$log_moment - log_target
  numerical_floor <- 128 * .Machine$double.eps * (1 + abs(log_target))
  if (!is.finite(error) || abs(error) > max(control$capacity_tol, numerical_floor))
    stop("Numerical eta did not satisfy the relative moment tolerance.", call. = FALSE)
  value <- g$log_z + exp(root + log_target)
  if (!is.finite(value) || value <= 0)
    stop("Positive entropy is outside numerical range.", call. = FALSE)
  c(eta = value, beta = g$beta, relative_moment_error = expm1(error),
    support = g$support, tail_relative_bound = g$tail_relative_bound)
}

#' Quantized maximum-entropy capacity
#'
#' Numerically solves the Gibbs distribution over nearest-bin quantizer costs
#' c_j = (abs(j) - 1/2)^2 for nonzero j, with c_0 = 0. The target minimum
#' second moment is T/(12 T_q). This evaluates the defining supremum; it is
#' not the entropy of a fitted Gaussian distribution. A log-moment equation
#' avoids premature convergence at small positive temperatures.
#' @param temperature Nonnegative temperature, scalar or vector.
#' @param resolution Positive T_q, scalar or vector of compatible length.
#' @param control Numerical controls.
#' @param details Return solver diagnostics in a data frame.
#' @return Numeric capacity in nats, or diagnostics with one row per coordinate.
#'   The value at zero is the mathematical boundary, not an inferential estimate.
#' @export
eta_numerical <- function(temperature, resolution, control = lambda_control(),
                          details = FALSE) {
  control <- .control(control)
  z <- .coordinates(temperature, resolution)
  if (!is.logical(details) || length(details) != 1 || is.na(details))
    stop("details must be TRUE or FALSE.", call. = FALSE)
  values <- t(vapply(z$ratio, .eta_scalar, numeric(5), control = control))
  if (!details) return(unname(values[, "eta"]))
  data.frame(ratio = z$ratio, values, row.names = NULL)
}

#' Smooth capacity approximation and geometric activation
#' @param temperature Nonnegative temperature.
#' @param resolution Positive T_q; scalar recycling is supported.
#' @return Numeric values. eta_app is a labeled approximation, xi_eta_app
#'   its primitive, geometry_activation the Hill profile and xi_a its primitive.
#'   Zero-temperature values are mathematical boundary values only.
#' @export
eta_app <- function(temperature, resolution) {
  r <- .coordinates(temperature, resolution)$ratio
  # log1p(a*r) without overflowing a*r at extreme but finite r.
  a <- pi * exp(1) / 6
  value <- numeric(length(r))
  ordinary <- r <= .Machine$double.xmax / a
  value[ordinary] <- 0.5 * log1p(a * r[ordinary])
  value[!ordinary] <- 0.5 * (log(a) + log(r[!ordinary]))
  value
}

#' @rdname eta_app
#' @export
xi_eta_app <- function(temperature, resolution) {
  z <- .coordinates(temperature, resolution)
  r <- (pi * exp(1) / 6) * z$ratio
  small <- r < 1e-3
  factor <- numeric(length(r))
  # [(1+r) log(1+r)-r]/r = r/2-r^2/6+r^3/12-...
  v <- r[small]
  factor[small] <- v * (1/2 + v * (-1/6 + v * (1/12 + v * (-1/20 + v/30))))
  v <- r[!small]
  if (any(!small)) factor[!small] <- (1 + 1/v) *
    (2 * eta_app(z$temperature[!small], z$resolution[!small])) - 1
  .finite(z$temperature * 0.5 * factor)
}

#' @rdname eta_app
#' @export
geometry_activation <- function(temperature, resolution) {
  z <- .coordinates(temperature, resolution)
  1 / (1 + 1 / z$ratio)
}

#' @rdname eta_app
#' @export
xi_a <- function(temperature, resolution) {
  z <- .coordinates(temperature, resolution)
  r <- z$ratio
  small <- r < 1e-3
  factor <- numeric(length(r))
  v <- r[small]
  factor[small] <- v * (1/2 + v * (-1/3 + v * (1/4 + v * (-1/5 + v/6))))
  v <- r[!small]
  factor[!small] <- 1 - log1p(v)/v
  .finite(z$temperature * factor)
}

.xi_direct <- function(temperature, resolution, control) {
  if (temperature == 0) return(c(value = 0, abs_error = 0))
  ratio <- .coordinates(temperature, resolution)$ratio
  top <- .eta_scalar(ratio, control)[["eta"]]
  scale <- temperature * top
  if (!is.finite(scale) || scale == 0)
    stop("Positive temperature primitive is outside numerical range.", call. = FALSE)
  # Scale both axis and ordinate: quadrature is O(1) even for tiny T/T_q.
  integrand <- function(u) vapply(u, function(v) {
    .eta_scalar(ratio * v, control)[["eta"]] / top
  }, numeric(1))
  result <- stats::integrate(integrand, 0, 1,
    rel.tol = control$integration_rel_tol,
    abs.tol = min(1, control$integration_abs_tol / scale),
    subdivisions = control$subdivisions, stop.on.error = TRUE)
  c(value = .finite(scale * result$value), abs_error = scale * result$abs.error)
}

#' Build a reusable capacity integral lookup
#'
#' An adaptively refined mesh approximates eta by linear segments. Their
#' integrals are accumulated and evaluated quadratically inside each segment.
#' Below min_ratio, direct quadrature is used, preventing extrapolation of
#' the first interval into the singular boundary region. Local trapezoidal
#' errors are estimated by refinement; these estimates are not rigorous bounds.
#' @param max_ratio Largest T/T_q covered (strictly positive).
#' @param min_ratio Smallest cached ratio, default max_ratio times 1e-6.
#' @param grid_points Initial mesh size, integer >= 3.
#' @param control Numerical controls.
#' @return An eta_numerical_lookup object with nodes, capacities, cumulative
#'   integrals, estimated errors, and its actual numerical controls.
#' @export
make_eta_numerical_lookup <- function(max_ratio, min_ratio = max_ratio * 1e-6,
                                    grid_points = 65L, control = lambda_control()) {
  control <- .control(control)
  .scalar(max_ratio, "max_ratio", positive = TRUE)
  .scalar(min_ratio, "min_ratio", positive = TRUE)
  .scalar(grid_points, "grid_points", positive = TRUE)
  if (min_ratio >= max_ratio || grid_points < 3 || grid_points != floor(grid_points) ||
      grid_points > control$lookup_max_points)
    stop("Invalid lookup range or grid_points.", call. = FALSE)
  x <- exp(seq(log(min_ratio), log(max_ratio), length.out = grid_points))
  x[c(1, length(x))] <- c(min_ratio, max_ratio)
  y <- eta_numerical(x, 1, control)
  baseline <- .xi_direct(min_ratio, 1, control)
  repeat {
    dx <- diff(x)
    midpoint <- x[-length(x)] + dx / 2
    ym <- eta_numerical(midpoint, 1, control)
    yl <- y[-length(y)]
    yr <- y[-1L]
    coarse <- dx * (yl + yr) / 2
    fine <- dx * (yl + 2 * ym + yr) / 4
    errors <- 4 * abs(fine - coarse) / 3
    refine <- errors > control$lookup_rel_tol * abs(fine) +
      control$lookup_abs_tol * dx / (max_ratio - min_ratio)
    if (!any(refine)) break
    if (length(x) + sum(refine) > control$lookup_max_points)
      stop("Lookup accuracy requires more than lookup_max_points.", call. = FALSE)
    nx <- c(x, midpoint[refine])
    ny <- c(y, ym[refine])
    order <- order(nx)
    x <- nx[order]
    y <- ny[order]
  }
  structure(list(format_version = 1L, min_ratio = min_ratio, max_ratio = max_ratio,
    ratio = x, eta = y, primitive = baseline[["value"]] + c(0, cumsum(coarse)),
    cumulative_error = baseline[["abs_error"]] + c(0, cumsum(errors)),
    interval_error = errors, control = control), class = "eta_numerical_lookup")
}

.check_lookup <- function(lookup, control) {
  if (!inherits(lookup, "eta_numerical_lookup") || !identical(lookup$format_version, 1L))
    stop("Use a lookup created by make_eta_numerical_lookup().", call. = FALSE)
  n <- length(lookup$ratio)
  if (n < 3 || length(lookup$eta) != n || length(lookup$primitive) != n ||
      length(lookup$cumulative_error) != n || length(lookup$interval_error) != n-1 ||
      any(!is.finite(c(lookup$ratio, lookup$eta, lookup$primitive,
                       lookup$cumulative_error, lookup$interval_error))) ||
      any(diff(lookup$ratio) <= 0) || lookup$ratio[1] != lookup$min_ratio ||
      lookup$ratio[n] != lookup$max_ratio)
    stop("Malformed capacity lookup.", call. = FALSE)
  built <- .control(lookup$control)
  for (key in c("capacity_tol", "tail_tol", "lookup_rel_tol", "lookup_abs_tol",
                "integration_rel_tol", "integration_abs_tol"))
    if (built[[key]] > control[[key]])
      stop("Lookup was built with looser ", key, "; rebuild it.", call. = FALSE)
}

#' Integrate numerical capacity over temperature
#' @param temperature Nonnegative temperature, scalar or vector.
#' @param resolution Positive T_q, scalar or vector.
#' @param lookup Optional reusable lookup; otherwise direct adaptive quadrature.
#' @param control Numerical controls. A reused lookup must meet these settings.
#' @param details Return values, estimated absolute errors and evaluation methods.
#' @return Numeric integral values, or a list containing value, abs_error,
#'   method and control. Lookup error estimates describe mesh discretization;
#'   solver/truncation tolerances are recorded separately in control.
#' @export
xi_eta_numerical <- function(temperature, resolution, lookup = NULL,
                             control = lambda_control(), details = FALSE) {
  control <- .control(control)
  z <- .coordinates(temperature, resolution)
  if (!is.logical(details) || length(details) != 1 || is.na(details))
    stop("details must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(lookup)) {
    .check_lookup(lookup, control)
    if (any(z$ratio > lookup$max_ratio))
      stop("Lookup does not cover the requested ratio.", call. = FALSE)
  }
  value <- error <- numeric(length(z$ratio))
  method <- character(length(value))
  for (i in seq_along(value)) {
    r <- z$ratio[i]
    if (is.null(lookup) || r < lookup$min_ratio) {
      direct <- .xi_direct(z$temperature[i], z$resolution[i], control)
      value[i] <- direct[["value"]]
      error[i] <- direct[["abs_error"]]
      method[i] <- if (r == 0) "boundary" else "quadrature"
    } else {
      j <- min(findInterval(r, lookup$ratio), length(lookup$ratio) - 1L)
      dx <- r - lookup$ratio[j]
      slope <- (lookup$eta[j+1] - lookup$eta[j]) /
        (lookup$ratio[j+1] - lookup$ratio[j])
      value[i] <- z$resolution[i] * (lookup$primitive[j] +
                                      lookup$eta[j] * dx + slope * dx^2 / 2)
      error[i] <- z$resolution[i] * (lookup$cumulative_error[j] + lookup$interval_error[j])
      method[i] <- "lookup"
    }
  }
  .finite(value)
  if (details) list(value = value, abs_error = error, method = method, control = control)
  else value
}
