#' Numerical settings for the Lambda Star core
#'
#' @param rank_tol Relative rank tolerance after column normalization.
#' @param capacity_tol Relative target-moment tolerance for the Gibbs solver.
#' @param tail_tol Relative truncation tolerance for partition and moment tails.
#' @param integration_rel_tol Relative quadrature tolerance.
#' @param integration_abs_tol Absolute quadrature tolerance (zero by default).
#' @param subdivisions Maximum adaptive quadrature subdivisions.
#' @param lookup_rel_tol Relative local trapezoidal-error target for lookups.
#' @param lookup_abs_tol Absolute dimensionless lookup-error target.
#' @param lookup_max_points Maximum lookup mesh size.
#' @param max_support Maximum positive quantizer-bin index.
#' @return A validated list of numerical controls. Lookup error estimates are
#'   numerical estimates, not rigorous error certificates.
#' @export
lambda_control <- function(rank_tol = 1e-10, capacity_tol = 1e-10,
                           tail_tol = 1e-14, integration_rel_tol = 1e-8,
                           integration_abs_tol = 0, subdivisions = 200L,
                           lookup_rel_tol = 1e-5, lookup_abs_tol = 0,
                           lookup_max_points = 32769L, max_support = 1000000L) {
  out <- as.list(environment())
  for (key in c("rank_tol", "capacity_tol", "tail_tol", "integration_rel_tol",
                "lookup_rel_tol")) {
    .scalar(out[[key]], key, positive = TRUE)
    if (out[[key]] >= 1 || out[[key]] < 1e-15)
      stop(key, " must be at least 1e-15 and less than one.", call. = FALSE)
  }
  for (key in c("integration_abs_tol", "lookup_abs_tol"))
    .scalar(out[[key]], key, nonnegative = TRUE)
  for (key in c("subdivisions", "lookup_max_points", "max_support")) {
    .scalar(out[[key]], key, positive = TRUE)
    if (out[[key]] != floor(out[[key]]) || out[[key]] > .Machine$integer.max)
      stop(key, " must be an integer in range.", call. = FALSE)
  }
  if (lookup_max_points < 17 || max_support < 8)
    stop("lookup_max_points must be >= 17 and max_support >= 8.", call. = FALSE)
  class(out) <- "lambda_control"
  out
}

.control <- function(control) {
  if (!inherits(control, "lambda_control"))
    stop("control must be created by lambda_control().", call. = FALSE)
  do.call(lambda_control, unclass(control))
}

.scalar <- function(x, name, positive = FALSE, nonnegative = FALSE) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x))
    stop(name, " must be one finite real number.", call. = FALSE)
  if ((positive && x <= 0) || (nonnegative && x < 0))
    stop(name, if (positive) " must be positive." else " must be nonnegative.",
         call. = FALSE)
  invisible(x)
}

.matrix <- function(x, name) {
  x <- as.matrix(x)
  if (!is.numeric(x) || is.complex(x) || nrow(x) < 2L || any(!is.finite(x)))
    stop(name, " must be a finite real numeric matrix with at least two rows.",
         call. = FALSE)
  x
}

.vector <- function(x, name) {
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) ||
      length(x) < 2L || any(!is.finite(x)))
    stop(name, " must be a finite real numeric vector of length >= 2.", call. = FALSE)
  x
}

.finite <- function(x) {
  if (any(!is.finite(x)))
    stop("Numerical range exceeded; rescale the input explicitly.", call. = FALSE)
  x
}

.ssq <- function(x) {
  if (!length(x) || all(x == 0)) return(0)
  a <- max(abs(x))
  value <- .finite(a * sqrt(sum((x / a)^2)))^2
  if (!is.finite(value) || value == 0)
    stop("Squared energy is outside the representable numerical range.", call. = FALSE)
  value
}

.rows <- function(a, b) {
  if (nrow(a) != nrow(b)) stop("Inputs have incompatible row counts.", call. = FALSE)
  if (!is.null(rownames(a)) && !is.null(rownames(b)) &&
      !identical(rownames(a), rownames(b)))
    stop("Inputs have different row identities or row order.", call. = FALSE)
}

.singular_message <- function() {
  "No background fluctuation; regular estimation is not possible."
}
