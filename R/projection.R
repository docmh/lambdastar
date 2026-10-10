#' Center a numeric matrix
#' @param x Numeric matrix, with observations in rows.
#' @param name Name used in diagnostics.
#' @return Matrix centered by column, retaining row and column names.
#' @export
center_columns <- function(x, name = "x") {
  x <- .matrix(x, name)
  .finite(sweep(x, 2L, colMeans(x), "-"))
}

#' Construct a centered hypothesis basis
#'
#' The returned X retains the encoded design. Q is used only for projection.
#' Rank is tested after column normalization, making the decision insensitive
#' to changes in predictor units. Redundant or constant columns are rejected.
#' @param x Finite numeric design without an intercept. Zero columns specify H0.
#' @param control Numerical controls from [lambda_control()].
#' @return List containing X, orthonormal Q, rank, log_volume, volume and control.
#'   log_volume is retained even when exponentiation underflows.
#' @export
hypothesis_basis <- function(x, control = lambda_control()) {
  control <- .control(control)
  x <- center_columns(x)
  .encoded_basis(x, control, nrow(x) - 1L)
}

# Already realized in the retained space: never change the declared origin here.
.encoded_basis <- function(x, control, dimension) {
  p <- ncol(x)
  if (!p) return(list(X = x, Q = x, rank = 0L, log_volume = 0,
                      volume = 1, control = control))
  if (p > dimension)
    stop("x is rank deficient in the retained space.", call. = FALSE)
  norms <- vapply(seq_len(p), function(j) {
    a <- max(abs(x[, j]))
    if (a == 0) return(0)
    .finite(a * sqrt(sum((x[, j] / a)^2)))
  }, numeric(1))
  if (any(norms == 0)) stop("x has a constant or redundant column.", call. = FALSE)
  normalized <- sweep(x, 2L, norms, "/")
  singular <- svd(normalized, nu = 0, nv = 0)$d
  if (min(singular) <= control$rank_tol * max(singular))
    stop("x is rank deficient at the requested tolerance.", call. = FALSE)
  decomposition <- qr(normalized, tol = control$rank_tol, LAPACK = FALSE)
  if (decomposition$rank != p)
    stop("x is rank deficient at the requested tolerance.", call. = FALSE)
  log_volume <- min(0, sum(log(singular)))
  list(X = x, Q = qr.Q(decomposition), rank = p,
       log_volume = log_volume, volume = exp(log_volume), control = control)
}

#' Project a vector using an orthonormal basis
#' @param x Finite numeric vector.
#' @param q Orthonormal matrix with the same number of rows as x.
#' @return Numeric projected vector. No N by N projector is constructed.
#' @export
project_onto_basis <- function(x, q) {
  x <- .vector(x, "x")
  q <- .matrix(q, "q")
  if (nrow(q) != length(x)) stop("q and x have incompatible dimensions.", call. = FALSE)
  if (!ncol(q)) return(rep(0, length(x)))
  if (max(abs(crossprod(q) - diag(ncol(q)))) > 1e-8)
    stop("q must have orthonormal columns.", call. = FALSE)
  .project(x, q)
}

.project <- function(x, q) {
  if (!ncol(q)) return(rep(0, length(x)))
  .finite(drop(q %*% crossprod(q, x)))
}

#' Encoded design volume
#' @param x Centered-space hypothesis design (centering is performed internally).
#' @param logarithm Return log volume instead of its exponential.
#' @param control Numerical controls.
#' @return Scalar volume, or log volume. Null and one-direction models have V=1.
#' @export
design_volume <- function(x, logarithm = FALSE, control = lambda_control()) {
  if (!is.logical(logarithm) || length(logarithm) != 1 || is.na(logarithm))
    stop("logarithm must be TRUE or FALSE.", call. = FALSE)
  basis <- hypothesis_basis(x, control)
  if (logarithm) basis$log_volume else basis$volume
}
