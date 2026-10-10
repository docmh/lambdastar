#' Settings for case or cluster percentile bootstrap intervals
#'
#' Defaults are 5000 draws, 95 percent two-sided percentile intervals and seed 1.
#' Whole persons and predictors are resampled together. Every draw uses the same
#' selected estimator as the original data. Undefined values are retained as NA;
#' draws are never replaced. Any invalid draw suppresses the affected regular CI.
#' Optional quantiles conditional on valid draws are returned in a separate table.
#' @param B Integer number of draws, at least 2.
#' @param conf.level Confidence level strictly between zero and one.
#' @param ci.sides Currently only "two.sided".
#' @param ci.method Currently only "percentile". This controls interval
#'   construction, not the point estimator.
#' @param resampling "case" for whole persons or "cluster" for whole clusters,
#'   jointly including measurements, predictors and nuisance variables.
#' @param strata NULL or one data column name. Resample within each stratum,
#'   retaining its observed number of persons (case) or clusters (cluster).
#'   At least two persons per stratum are required
#'   for estimation of the fixed-mixture population variance.
#' @param seed Integer from zero to .Machine$integer.max. The caller's random
#'   state and RNG kind are restored, including absence of .Random.seed.
#' @param min_valid Minimum finite draws for quantiles, default 200.
#' @param conditional_quantiles FALSE by default. If TRUE, separately report
#'   quantiles conditional on valid draws when some draws are undefined. These
#'   are not confidence intervals for the original unconditional target.
#' @param keep_draws Retain the actual parameter estimates for all B draws,
#'   including zeros and NA, plus their per-parameter failure reasons.
#' @param keep_indices Retain the B by N matrix of original input row positions.
#'   Cluster resampling instead returns a list of variable-length row vectors.
#' @param cluster Column identifying clusters for resampling="cluster".
#'   Labels are scoped within strata. Each stratum draws its original number
#'   of clusters with replacement, retaining all rows of every selected cluster.
#' @param percentile_interpolation "boot" (default) uses boot::boot.ci's
#'   normal-quantile interpolation; "type7" uses stats::quantile(type=7), as
#'   in the manuscript. Both compute ordinary equal-tail percentile intervals.
#' @return Settings for the bootstrap argument of lambda_measure(),
#'   lambda_model() or lambda_compare().
#' @export
lambda_bootstrap_control <- function(B = 5000L, conf.level = 0.95,
    ci.sides = "two.sided", ci.method = "percentile", resampling = "case",
    strata = NULL, seed = 1L, min_valid = 200L, conditional_quantiles = FALSE,
    keep_draws = FALSE, keep_indices = FALSE, cluster = NULL,
    percentile_interpolation = "boot") {
  out <- as.list(environment())
  for (key in c("B", "seed", "min_valid")) {
    .scalar(out[[key]], key, nonnegative = TRUE)
    if (out[[key]] != floor(out[[key]]) || out[[key]] > .Machine$integer.max ||
        (key != "seed" && out[[key]] < 2))
      stop(key, " must be an integer in range (B and min_valid >= 2).", call. = FALSE)
    out[[key]] <- as.integer(out[[key]])
  }
  .scalar(conf.level, "conf.level", positive = TRUE)
  if (conf.level >= 1) stop("conf.level must be less than one.", call. = FALSE)
  for (key in c("ci.sides", "ci.method")) {
    expected <- c(ci.sides = "two.sided", ci.method = "percentile")
    if (!identical(out[[key]], unname(expected[[key]])))
      stop("Only ", expected[[key]], " is supported for ", key, ".", call. = FALSE)
  }
  out$resampling <- match.arg(resampling, c("case", "cluster"))
  out$percentile_interpolation <- match.arg(percentile_interpolation, c("boot", "type7"))
  if (!is.null(cluster) && (!is.character(cluster) || length(cluster) != 1L ||
      is.na(cluster) || !nzchar(cluster))) stop("cluster must be NULL or one column name.", call. = FALSE)
  if ((out$resampling == "cluster") != !is.null(cluster))
    stop("Supply cluster exactly when resampling = 'cluster'.", call. = FALSE)
  if (!is.null(strata) && (!is.character(strata) || length(strata) != 1L ||
                           is.na(strata) || !nzchar(strata)))
    stop("strata must be NULL or one column name.", call. = FALSE)
  for (key in c("keep_draws", "keep_indices", "conditional_quantiles"))
    if (!is.logical(out[[key]]) || length(out[[key]]) != 1L || is.na(out[[key]]))
      stop(key, " must be TRUE or FALSE.", call. = FALSE)
  class(out) <- "lambda_bootstrap_control"
  out
}

.bootstrap_control <- function(bootstrap) {
  if (identical(bootstrap, FALSE)) return(NULL)
  if (identical(bootstrap, TRUE)) return(lambda_bootstrap_control())
  if (!inherits(bootstrap, "lambda_bootstrap_control"))
    stop("bootstrap must be FALSE, TRUE, or lambda_bootstrap_control().", call. = FALSE)
  do.call(lambda_bootstrap_control, unclass(bootstrap))
}

.bootstrap_models <- function(result) {
  if (inherits(result, "lambda_comparison")) list(result$model_a, result$model_b) else
    if (inherits(result, "lambda_model")) list(result) else list()
}

# Resolve contrast functions once to numeric matrices. Model frames still run
# anew on every draw unless the design explicitly freezes its selected encoding.
.bootstrap_specs <- function(prepared, points) {
  lapply(seq_along(prepared$specs), function(i)
    .bootstrap_spec(prepared$specs[[i]], points[[i]]$model, prepared))
}

.bootstrap_spec <- function(spec, meta, prepared) {
  if (spec$type == "contrast") {
    spec$full <- .bootstrap_spec(spec$full, meta$full_model, prepared)
    if (!is.null(spec$base)) spec$base <- .bootstrap_spec(spec$base, meta$base_model, prepared)
    return(spec)
  }
  if (spec$type == "block") {
    spec$full <- .bootstrap_spec(spec$full, meta$full_model, prepared)
    return(spec)
  }
  if (spec$type == "combined") {
    spec$parts <- Map(function(s,m) .bootstrap_spec(s,m,prepared), spec$parts, meta$parts)
    return(spec)
  }
  if (spec$type == "matrix") return(spec)
  if (identical(spec$encoding, "fixed")) {
    encoded <- .encode_model(spec, prepared, sqrt(.Machine$double.eps))
    return(list(type = "matrix", x = encoded$x, metadata = encoded$metadata,
                row_scope = "selected", variables = character()))
  }
  spec$type <- "formula"
  spec$fit <- NULL
  spec$terms <- stats::terms(spec$formula, data = prepared$data)
  spec$xlevels <- meta$factor_levels[setdiff(names(meta$factor_levels), meta$logical_factors)]
  spec$contrasts <- lapply(names(meta$contrasts), function(name) {
    contrast <- meta$contrasts[[name]]
    # contrasts<- pads reduced matrices to G-1 columns. Keep an explicit
    # coding unchanged so resampling uses exactly the supplied directions.
    if (is.matrix(contrast)) return(contrast)
    f <- factor(meta$factor_levels[[name]], levels = meta$factor_levels[[name]])
    stats::contrasts(f) <- contrast
    stats::contrasts(f)
  })
  names(spec$contrasts) <- names(meta$contrasts)
  if (!length(spec$contrasts)) spec$contrasts <- NULL
  spec
}

.bootstrap_groups <- function(prepared, config) {
  n <- nrow(prepared$data)
  if (is.null(config$strata)) return(list(all = seq_len(n)))
  s <- prepared$data[[config$strata]]
  if ((!is.atomic(s) && !is.factor(s)) || anyNA(s) ||
      (is.numeric(s) && any(!is.finite(s))))
    stop("strata must contain finite, non-missing group labels.", call. = FALSE)
  groups <- split(seq_len(n), factor(as.character(s), levels = unique(as.character(s))))
  if (config$resampling == "case" && all(lengths(groups) == 1L))
    stop("All strata are singletons; no person resampling is possible.", call. = FALSE)
  groups
}

.cluster_groups <- function(prepared, config) {
  if (config$resampling != "cluster") return(NULL)
  data <- prepared$data
  labels <- data[[config$cluster]]
  if (anyNA(labels) || (is.numeric(labels) && any(!is.finite(labels))))
    stop("cluster must contain finite, non-missing labels.", call. = FALSE)
  strata <- if (is.null(config$strata)) rep("all", nrow(data)) else as.character(data[[config$strata]])
  lapply(split(seq_len(nrow(data)), strata), function(g) split(g, as.character(labels[g])))
}

.bootstrap_attempt <- function(code) {
  tryCatch(withCallingHandlers(force(code), warning = function(w)
    stop("Warning during bootstrap evaluation: ", conditionMessage(w), call. = FALSE)),
    error = identity)
}

.bootstrap_error_code <- function(e) {
  message <- conditionMessage(e)
  if (grepl("rank deficient|constant or redundant", message)) "rank_deficient"
  else if (grepl("factor levels|contrasts can be applied|new levels", message)) "factor_levels"
  else if (grepl("coding changed", message)) "coding_changed"
  else "evaluation_error"
}

.percentile_limits <- function(values, level, estimate, interpolation = "boot") {
  # Degenerate empirical distributions have identical lower/upper percentiles.
  if (all(values == values[1]))
    return(list(limits = rep(values[1], 2L), message = "constant_distribution"))
  if (interpolation == "type7") return(list(limits = unname(stats::quantile(values,
    c((1-level)/2, (1+level)/2), type = 7)), message = ""))
  object <- structure(list(t0 = estimate, t = matrix(values, ncol = 1),
    R = length(values), sim = "ordinary"), class = "boot")
  messages <- character()
  interval <- withCallingHandlers(
    boot::boot.ci(object, conf = level, type = "perc", index = 1),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  limits <- unname(interval$percent[1, 4:5])
  if(length(limits)!=2L) {
    limits <- c(NA_real_,NA_real_)
    messages <- c(messages,"boot::boot.ci did not return percentile limits.")
  }
  list(limits = limits,
       message = paste(unique(messages), collapse = "; "))
}

.bootstrap_ci_table <- function(estimates, draws, config, singular = FALSE,
                                reasons = NULL) {
  attempted <- nrow(draws)
  if (is.null(reasons)) reasons <- matrix(NA_character_, nrow(draws), ncol(draws))
  conditional <- list()
  records <- lapply(seq_along(estimates), function(j) {
    values <- draws[, j]
    valid <- is.finite(values)
    n_valid <- sum(valid)
    n_failed <- attempted - n_valid
    n_undefined <- sum(reasons[, j] == "undefined_kappa", na.rm = TRUE)
    limits <- c(NA_real_, NA_real_)
    message <- ""
    status <- if (singular) "singular_temperature" else
      if (!is.finite(estimates[j])) "undefined_estimate" else
      if (n_undefined) "undefined_draws" else
      if (n_failed) "failed_draws" else
      if (n_valid < config$min_valid) "too_few_valid" else "ok"
    if (status == "ok") {
      interval <- .percentile_limits(values, config$conf.level, estimates[j], config$percentile_interpolation)
      limits <- interval$limits
      message <- interval$message
      if (any(!is.finite(limits))) {
        limits[] <- NA_real_
        status <- "quantile_unavailable"
      }
    }
    if (config$conditional_quantiles && !singular && n_failed &&
        n_valid >= config$min_valid) {
      interval <- .percentile_limits(values[valid], config$conf.level, estimates[j], config$percentile_interpolation)
      status_conditional <- "conditional_on_valid"
      if (any(!is.finite(interval$limits))) {
        interval$limits[] <- NA_real_
        status_conditional <- "quantile_unavailable"
      }
      conditional[[length(conditional) + 1L]] <<- data.frame(
        parameter = names(estimates)[j], lower = interval$limits[1], upper = interval$limits[2],
        central_mass = config$conf.level, n_valid = n_valid, n_excluded = n_failed,
        status = status_conditional, quantile_message = interval$message)
    }
    data.frame(parameter = names(estimates)[j], estimate = unname(estimates[j]),
      lower = limits[1], upper = limits[2], conf.level = config$conf.level,
      B = config$B, attempted = attempted, n_valid = n_valid, n_failed = n_failed,
      n_undefined = n_undefined,
      failure_rate = if (attempted) n_failed / attempted else NA_real_,
      n_zero = sum(values[valid] == 0), expected_tail_count = n_valid * (1-config$conf.level)/2,
      status = status, method = "percentile", quantile_message = message, row.names = NULL)
  })
  list(ci = do.call(rbind, records), conditional_quantiles = if (length(conditional))
    do.call(rbind, conditional) else data.frame(parameter = character(), lower = numeric(),
      upper = numeric(), central_mass = numeric(), n_valid = integer(), n_excluded = integer(),
      status = character(), quantile_message = character()))
}

.attach_bootstrap <- function(result, prepared, config, model_settings) {
  started <- proc.time()[["elapsed"]]
  if (prepared$estimator == "population" && !is.null(config) && config$resampling == "cluster")
    stop("The population estimator currently supports case resampling only; use estimator = 'paper' for cluster resampling.", call. = FALSE)
  groups <- if (is.null(config)) list(all = seq_len(nrow(prepared$data))) else .bootstrap_groups(prepared, config)
  selected <- .select_estimator(result, prepared, groups, !is.null(config$strata))
  result <- selected$result; context <- selected$context; estimates <- result$estimates
  if (is.null(config)) return(result)
  clusters <- .cluster_groups(prepared, config)
  panel <- prepared$panel
  if(!is.null(panel)&&!is.null(config$strata)) {
    strata_per_person <- split(as.character(prepared$data[[config$strata]]),prepared$data[[panel$person]])
    if(any(vapply(strata_per_person,function(x)length(unique(x)),integer(1))!=1L))
      stop("Panel bootstrap strata must be constant within each person.",call.=FALSE)
  }
  models <- .bootstrap_models(result)
  state <- prepared$measurement
  attempted <- if (context$singular) 0L else config$B
  draws <- matrix(NA_real_, attempted, length(estimates), dimnames = list(NULL, names(estimates)))
  reasons <- matrix(NA_character_, attempted, length(estimates), dimnames = dimnames(draws))
  indices <- if (!config$keep_indices) NULL else if (is.null(clusters))
    matrix(NA_integer_, attempted, state$n) else vector("list", attempted)
  draw_sizes <- draw_df <- rep(NA_integer_, attempted)
  person_instances <- if(!is.null(panel)&&config$keep_indices) vector("list",attempted) else NULL
  failures <- list()
  record <- function(b, scope, code, message, columns) {
    reasons[b, columns] <<- code
    failures[[length(failures) + 1L]] <<- data.frame(draw = b, scope = scope,
      code = code, message = message, stringsAsFactors = FALSE)
  }
  if (attempted) {
    specs <- .bootstrap_specs(prepared, models)
    cached <- lapply(seq_along(specs), function(i)
      if (.simple_spec(specs[[i]])) models[[i]]$encoded_x else NULL)
    nuisance_spec <- if (is.null(prepared$nuisance_spec)) NULL else
      .bootstrap_spec(prepared$nuisance_spec, prepared$nuisance_metadata, prepared)
    nuisance_cached <- if (!is.null(nuisance_spec) && .simple_spec(nuisance_spec)) prepared$nuisance_x else NULL
    items <- as.matrix(prepared$data[state$settings$indicators])
    control <- state$settings$numerical_control
    saved_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) saved_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      do.call(RNGkind, as.list(saved_kind))
      if (had_seed) assign(".Random.seed", saved_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(config$seed)
    for (b in seq_len(attempted)) {
      if (is.null(clusters)) {
        index <- seq_len(state$n)
        for (g in groups) index[g] <- g[sample.int(length(g), length(g), replace = TRUE)]
      } else {
        picked <- unlist(lapply(clusters,function(g)
          g[sample.int(length(g),length(g),replace=TRUE)]),recursive=FALSE,use.names=FALSE)
        index <- unlist(picked,use.names=FALSE)
      }
      size <- length(index); draw_sizes[b] <- size
      if (config$keep_indices) {
        if (is.null(clusters)) indices[b, ] <- prepared$rows[index]
        else indices[[b]] <- prepared$rows[index]
      }
      data <- prepared$data[index, , drop = FALSE]
      rownames(data) <- as.character(seq_len(size))
      if(!is.null(panel)) {
        instance <- rep(seq_along(picked),lengths(picked))
        labels <- levels(prepared$data[[panel$person]])
        data[[panel$person]] <- factor(labels[instance],levels=labels)
        if(config$keep_indices) person_instances[[b]] <- instance
      }
      sample <- list(data = data, rows = prepared$rows[index], local_rows = index,
        retain_intercept = prepared$retain_intercept,
        numerical_control = control,
        encoding_cache = new.env(parent = emptyenv()),
        measurement = list(composite_response = rowMeans(items[index, , drop = FALSE])))
      s <- .bootstrap_attempt({
        nuisance <- if (is.null(nuisance_spec)) NULL else if (!is.null(nuisance_cached))
          nuisance_cached[index, , drop = FALSE] else {
            enc <- .encode_model(nuisance_spec, sample, sqrt(.Machine$double.eps))
            if (!identical(enc$metadata$columns, prepared$nuisance_metadata$columns))
              stop("Nuisance coding changed during resampling.", call. = FALSE)
            enc$x
          }
        if (!is.null(nuisance)) rownames(nuisance) <- rownames(data)
        value <- estimate_lambda_star_items(items[index, , drop = FALSE], state$scale,
          rownames(data), .panel_nuisance(nuisance,data,panel), control,
          reference=if(is.null(state$reference)) NULL else state$reference$value[index])
        if (value$nuisance$rank != state$nuisance$rank)
          stop("Nuisance space is rank deficient relative to the original sample.", call. = FALSE)
        value
      })
      if (inherits(s, "error")) {
        record(b, "measurement", .bootstrap_error_code(s), conditionMessage(s), seq_along(estimates)); next
      }
      draw_df[b] <- s$d; draws[b, 1] <- s$temperature
      if (s$status != "ok") {
        record(b, "measurement", s$status, s$message, seq.int(2L, length(estimates))); next
      }
      sample$measurement <- s
      v <- if (prepared$estimator == "paper") .ssq(s$stable_profile)/s$d else
        .population_variance(s$composite_response, groups, !is.null(config$strata))
      draw_context <- context; draw_context$d <- s$d
      draw_context$raw <- c(temperature = s$temperature,
        variance_signal = if (prepared$estimator == "paper") s$signal_energy_raw/s$d else v-s$temperature_mean,
        context$raw[-(1:2)] * NA_real_)
      for (i in seq_along(specs)) {
        columns <- 2L + (i-1L)*3L + 1:3
        scope <- if (length(specs) == 1L) "model" else c("model_a", "model_b")[i]
        value <- .bootstrap_attempt({
          if (!is.null(cached[[i]])) x <- cached[[i]][index, , drop = FALSE] else {
            encoded <- .encode_model(specs[[i]], sample, model_settings$response_tol)
            if (!identical(encoded$metadata$columns, models[[i]]$model$columns))
              stop("Model coding changed during resampling.", call. = FALSE)
            x <- encoded$x
          }
          basis <- .conditional_basis(x, s$nuisance, control)
          if (prepared$estimator == "population") .population_h(s, basis, v) else {
            partition <- .partition_measurement(s, basis)
            (partition$captured_observed - basis$rank*s$temperature_mean)/s$d
          }
        })
        if (inherits(value, "error")) record(b, scope, .bootstrap_error_code(value), conditionMessage(value), columns)
        else {
          draw_context$raw[i+2L] <- value
          if (draw_context$raw[2] <= 0)
            record(b, scope, "undefined_kappa", "Estimated total signal is zero; kappa is undefined.", columns[3])
        }
      }
      draws[b, ] <- .component_estimates(draw_context)
    }
  }
  failures <- if (length(failures)) do.call(rbind, failures) else
    data.frame(draw = integer(), scope = character(), code = character(), message = character())
  diagnostics <- if (nrow(failures)) stats::aggregate(list(n = failures$draw),
    failures[c("scope", "code")], length) else data.frame(scope = character(), code = character(), n = integer())
  tables <- .bootstrap_ci_table(estimates, draws, config, context$singular, reasons)
  inference <- c(result$inference, list(method = "percentile", status = "computed",
    design = if (!is.null(clusters)) "clusters_within_strata" else
      if (is.null(config$strata)) "random_persons_and_predictors" else "fixed_observed_stratum_counts",
    percentile_interpolation = config$percentile_interpolation,
    simultaneous = FALSE, F_H_intervals = FALSE, kappa_requires_positive_signal = TRUE))
  result$bootstrap <- list(config = config, ci = tables$ci,
    conditional_quantiles = tables$conditional_quantiles, diagnostics = diagnostics,
    failures = failures, attempted = attempted, skipped = config$B - attempted,
    draws = if (config$keep_draws) draws else NULL,
    draw_reasons = if (config$keep_draws) reasons else NULL, indices = indices,
    sampling = list(input_rows = prepared$rows, row_ids = state$row_ids,
      stratum_sizes = lengths(groups), cluster_counts = if (is.null(clusters)) NULL else lengths(clusters),
      draw_sizes = draw_sizes, draw_df = draw_df,
      rng_kind = c("Mersenne-Twister", "Inversion", "Rejection")),
    inference = inference, elapsed_seconds = proc.time()[["elapsed"]] - started)
  result$settings$bootstrap <- config
  if(!is.null(panel)) {
    result$bootstrap$sampling$person_instances <- person_instances
    result$bootstrap$sampling$unit <- "person"
    result$bootstrap$sampling$person_count <- result$panel$n_persons
  }
  result
}

#' Extract Lambda Star percentile intervals
#'
#' The stored confidence level is used when level is omitted. A different level
#' requires keep_draws=TRUE and reuses the same parameter draws without resampling.
#' The diagnostics attribute includes estimates, valid/undefined draw counts and
#' explicit statuses for missing intervals. Conditional quantiles are available
#' separately in object$bootstrap$conditional_quantiles. See [lambda_inference].
#' @param object A lambda_measurement, lambda_model or lambda_comparison computed
#'   with bootstrap=TRUE or lambda_bootstrap_control().
#' @param parm Optional parameter names or integer positions in the CI table.
#' @param level Confidence level; defaults to the level used during estimation.
#' @param ... Reserved; unsupported arguments produce an error.
#' @return A lambda_ci matrix with lower/upper columns, diagnostics and inference
#'   attributes. Suppressed intervals contain NA and an explicit status.
#' @export
confint.lambda_measurement <- function(object, parm, level = 0.95, ...) {
  if (length(list(...))) stop("Unsupported confint arguments.", call. = FALSE)
  boot <- object$bootstrap
  if (is.null(boot)) stop("No bootstrap computed; refit with bootstrap = TRUE or lambda_bootstrap_control().", call. = FALSE)
  if (missing(level)) level <- boot$config$conf.level
  config <- boot$config
  config$conf.level <- level
  config <- do.call(lambda_bootstrap_control, unclass(config))
  table <- boot$ci
  if (level != boot$config$conf.level) {
    if (is.null(boot$draws)) stop("Changing level requires keep_draws = TRUE.", call. = FALSE)
    estimates <- stats::setNames(table$estimate, table$parameter)
    table <- .bootstrap_ci_table(estimates, boot$draws, config,
      object$status == "singular_temperature", boot$draw_reasons)$ci
  }
  if (!missing(parm)) {
    if (is.character(parm) && !anyNA(parm)) index <- match(parm, table$parameter)
    else if (is.numeric(parm) && !anyNA(parm) && all(is.finite(parm)) &&
             all(parm == floor(parm))) index <- parm
    else stop("parm must contain parameter names or integer positions.", call. = FALSE)
    if (!length(index) || anyNA(index) || any(index < 1 | index > nrow(table)) || anyDuplicated(index))
      stop("Unknown or duplicated parm.", call. = FALSE)
    table <- table[index, , drop = FALSE]
  }
  out <- as.matrix(table[c("lower", "upper")])
  dimnames(out) <- list(table$parameter, paste0(format(100 * c((1-level)/2, (1+level)/2), trim = TRUE), " %"))
  attr(out, "diagnostics") <- table
  attr(out, "inference") <- boot$inference
  class(out) <- c("lambda_ci", "matrix", "array")
  out
}

#' @rdname confint.lambda_measurement
#' @export
confint.lambda_model <- function(object, parm, level = 0.95, ...) {
  confint.lambda_measurement(object, parm, if (missing(level)) object$bootstrap$config$conf.level else level, ...)
}

#' @rdname confint.lambda_measurement
#' @export
confint.lambda_comparison <- function(object, parm, level = 0.95, ...) {
  confint.lambda_measurement(object, parm, if (missing(level)) object$bootstrap$config$conf.level else level, ...)
}

#' @rdname lambda_results
#' @export
print.lambda_ci <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Pointwise percentile bootstrap intervals\n")
  cat("Estimator:", attr(x, "inference")$estimator, "| interpolation:",
      attr(x, "inference")$percentile_interpolation, "\n")
  d <- attr(x, "diagnostics")
  print(d[c("parameter", "estimate", "lower", "upper", "n_valid", "n_failed", "status")],
        digits = digits, row.names = FALSE)
  invisible(x)
}
