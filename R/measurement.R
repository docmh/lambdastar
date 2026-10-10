#' Estimate intrinsic structure from parallel measurements
#'
#' Measurements must target the same stable profile with independent,
#' homogeneous replicate errors. Homogeneous isotropic fluctuation refers to
#' measurement errors, not the total
#' response variance within groups. Finite second moments are required; a Gaussian
#' response distribution is not required for these moment/projection calculations.
#' Complete, already scored items are required;
#' neither reverse scoring nor missing-value substitution is inferred.
#' @param items Finite numeric N by K matrix, K >= 2.
#' @param scale Noncentrality denominator: composite mean (default) or single item.
#' @param row_ids Optional unique observation identifiers in input order.
#' @param nuisance Optional numeric nuisance design on the same rows. The
#'   intercept is included by default, but not with a fixed reference;
#'   redundant nuisance columns count only once.
#' @param control Numerical controls for nuisance rank determination.
#' @param reference NULL retains ordinary centering. A known scalar or numeric
#'   vector on the input rows retains departures from that reference, with no
#'   implicit nuisance intercept. It must be prespecified, not fitted here.
#' @return Measurement state containing temperature, raw and constrained signal
#'   and noncentrality, retained composite profile and row identity. At zero
#'   temperature, status is singular_temperature and normalized estimates are NA.
#'   composite_response always refers to the mean; item scale only changes
#'   the noncentrality denominator, not this observed response.
#' @export
estimate_lambda_star_items <- function(items,
    scale = c("composite_mean", "item"), row_ids = NULL, nuisance = NULL,
    control = lambda_control(), reference = NULL) {
  scale <- match.arg(scale)
  items <- .matrix(items, "items")
  n <- nrow(items)
  k <- ncol(items)
  if (k < 2L) stop("items must contain at least two replicate columns.", call. = FALSE)
  if (is.null(row_ids)) row_ids <- rownames(items)
  if (is.null(row_ids)) row_ids <- as.character(seq_len(n))
  if (length(row_ids) != n || anyNA(row_ids) || anyDuplicated(row_ids))
    stop("row_ids must contain one unique non-missing ID per row.", call. = FALSE)
  origin <- .reference_values(reference,n,row_ids)
  nuisance_space <- .nuisance_space(nuisance, n, row_ids, control, center=is.null(origin))
  shifted <- if(is.null(origin)) items else .finite(sweep(items,1L,origin,"-"))
  z <- .residualize_nuisance(shifted, nuisance_space)
  z_bar <- rowMeans(z)
  d <- nuisance_space$d
  # Item disagreement is unchanged by a common row profile or reference.
  # Remove the first item before projection and row averaging: identical
  # indicators then give exact zeros even without extended-precision sums.
  # This also keeps a large common profile out of the fluctuation calculation,
  # without rounding genuinely small positive temperatures down to zero.
  item_offsets <- .finite(sweep(items, 1L, items[, 1L], "-"))
  fluctuation_items <- .residualize_nuisance(item_offsets, nuisance_space)
  deviations <- sweep(fluctuation_items, 1L, rowMeans(fluctuation_items), "-")
  sigma2_item <- .ssq(deviations) / (d * (k - 1L))
  temperature_mean <- sigma2_item / k
  if (any(deviations != 0) && temperature_mean == 0)
    stop("Positive temperature underflowed; rescale the input explicitly.", call. = FALSE)
  temperature <- if (scale == "item") sigma2_item else temperature_mean
  signal_energy_raw <- .finite(.ssq(z_bar) - d * temperature_mean)
  signal_energy <- max(0, signal_energy_raw)
  singular <- temperature == 0
  lambda_star_raw <- if (singular) NA_real_ else .finite(signal_energy_raw / temperature)
  lambda_star <- if (singular) NA_real_ else max(0, lambda_star_raw)
  result <- list(status = if (singular) "singular_temperature" else "ok",
       message = if (singular) .singular_message() else NULL,
       n = n, k = k, d = d, scale = scale, nuisance = nuisance_space,
       row_ids = row_ids, input_rows = seq_len(n), item_names = colnames(items),
       sigma2_item = sigma2_item, temperature_mean = temperature_mean,
       temperature = temperature, signal_energy_raw = signal_energy_raw,
       signal_energy = signal_energy, lambda_star_raw = lambda_star_raw,
       lambda_star = lambda_star,
       rho_star = if (singular) NA_real_ else lambda_star / (lambda_star + d),
       stable_profile = z_bar, composite_response = rowMeans(items))
  if(!is.null(origin)) {
    result$reference <- list(type="fixed",value=origin,row_ids=row_ids,implicit_intercept=FALSE)
    result$mean_departure <- mean(z_bar)
  }
  result
}

#' Estimate the noncentrality partition for an explicit design
#' @param items Complete parallel measurement matrix.
#' @param x Explicit encoded hypothesis matrix on the same rows as items.
#' @param scale Noncentrality denominator, composite_mean or item.
#' @param control Numerical controls.
#' @param row_ids Optional observation IDs.
#' @param nuisance Optional numeric nuisance design shared by measurements and
#'   hypotheses. An intercept is implicit only in the ordinary centered case.
#' @inheritParams estimate_lambda_star_items
#' @return Measurement state extended with raw and constrained partitions,
#'   degrees of freedom, observed energies and the encoded basis. The raw
#'   partition preserves additivity; the constrained partition caps captured
#'   structure at the nonnegative total. Kappa is NA when total structure is zero.
#' @export
estimate_lambda_partition_items <- function(items, x,
    scale = c("composite_mean", "item"), control = lambda_control(), row_ids = NULL,
    nuisance = NULL, reference = NULL) {
  items <- .matrix(items, "items")
  x <- .matrix(x, "x")
  .rows(items, x)
  state <- estimate_lambda_star_items(items, scale, row_ids, nuisance, control, reference)
  basis <- .conditional_basis(x, state$nuisance, control)
  .partition_measurement(state, basis)
}

# Shared by the matrix core and the model adapters; measurement is prepared once.
.partition_measurement <- function(state, basis) {
  # A saturated retained space is the identity, not an approximate QR projection.
  fitted <- if(basis$rank==state$d) state$stable_profile else .project(state$stable_profile, basis$Q)
  captured <- .ssq(fitted)
  residual <- .ssq(state$stable_profile - fitted)
  df_h <- basis$rank
  df_perp <- state$d - df_h
  raw_h <- raw_perp <- h <- perp <- kappa <- NA_real_
  if (state$status == "ok") {
    raw_h <- .finite((captured - df_h * state$temperature_mean) / state$temperature)
    raw_perp <- .finite((residual - df_perp * state$temperature_mean) / state$temperature)
    state$lambda_star_raw <- .finite(raw_h + raw_perp)
    state$lambda_star <- max(0, state$lambda_star_raw)
    h <- min(state$lambda_star, max(0, raw_h))
    perp <- state$lambda_star - h
    if (state$lambda_star > 0) kappa <- h / state$lambda_star
    state$rho_star <- state$lambda_star / (state$lambda_star + state$d)
  }
  c(state, list(df_h = df_h, df_perp = df_perp, lambda_h_raw = raw_h,
                lambda_perp_raw = raw_perp, lambda_h = h, lambda_perp = perp,
                kappa = kappa, captured_observed = captured,
                residual_observed = residual, basis = basis))
}
