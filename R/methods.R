#' Inspect Lambda Star results
#' @param object,x A Lambda Star result or summary.
#' @param digits Number of significant digits for printing.
#' @param ... Reserved for S3 compatibility.
#' @return Summary methods return inspectable summary lists. Print methods return
#'   their input invisibly. Full encoded matrices remain in the original results.
#' @name lambda_results
NULL

#' @rdname lambda_results
#' @export
summary.lambda_measurement <- function(object, ...) {
  estimates <- data.frame(quantity = "lambda_star", raw = object$lambda_star_raw,
                          estimate = object$lambda_star)
  result <- structure(list(status = object$status, message = object$message,
    n = object$n, k = object$k, cases = object$cases, settings = object$settings,
    temperatures = c(item = object$sigma2_item, mean = object$temperature_mean,
                     noncentrality_denominator = object$temperature),
    estimates = estimates, inference = object$inference, nuisance_rank = object$nuisance$rank,
    paper = object$paper, population = object$population,
    bootstrap = .bootstrap_summary(object$bootstrap)),
    class = "lambda_summary")
  if(!is.null(object$reference)) {
    result$reference <- object$reference
    result$mean_departure <- object$mean_departure
  }
  if(!is.null(object$panel)) result$panel <- object$panel
  result
}

#' @rdname lambda_results
#' @export
summary.lambda_model <- function(object, ...) {
  result <- summary.lambda_measurement(object)
  result$estimates <- rbind(result$estimates,
    data.frame(quantity = c("lambda_h", "lambda_perp", "kappa"),
      raw = c(object$lambda_h_raw, object$lambda_perp_raw, NA_real_),
      estimate = c(object$lambda_h, object$lambda_perp, object$kappa)))
  result$model <- object$model
  result$parsimony <- object$parsimony[c("df_h", "df_perp", "rss", "volume",
    "log_volume", "temperature", "resolution", "xi_eta", "xi_a", "eta_method",
    "dimension_cost", "geometry_cost", "free_energy", "numerical")]
  result
}

#' @rdname lambda_results
#' @export
summary.lambda_comparison <- function(object, ...) {
  a <- object$model_a$parsimony
  b <- object$model_b$parsimony
  structure(list(status = object$status, message = object$message,
    cases = object$cases, settings = object$settings, delta_F = object$delta_F,
    inference = object$inference, population = object$population, bootstrap = .bootstrap_summary(object$bootstrap),
    models = data.frame(model = c("a", "b"), df_h = c(a$df_h, b$df_h),
      RSS = c(a$rss, b$rss), dimension_cost = c(a$dimension_cost, b$dimension_cost),
      geometry_cost = c(a$geometry_cost, b$geometry_cost),
      F_H = c(a$free_energy, b$free_energy)),
    model_a = summary(object$model_a), model_b = summary(object$model_b)),
    class = "lambda_comparison_summary")
}

#' @rdname lambda_results
#' @export
print.lambda_measurement <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print(summary(x), digits = digits)
  invisible(x)
}

#' @rdname lambda_results
#' @export
print.lambda_model <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print(summary(x), digits = digits)
  invisible(x)
}

#' @rdname lambda_results
#' @export
print.lambda_comparison <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print(summary(x), digits = digits)
  invisible(x)
}

#' @rdname lambda_results
#' @export
print.lambda_summary <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Lambda Star point estimates | estimator:", x$inference$estimator, "\n")
  cat("Status:", x$status, "| N =", x$n, "| measurements =", x$k,
      "| excluded =", length(x$cases$excluded), "\n")
  if (!is.null(x$message)) cat(x$message, "\n")
  cat("Temperatures:\n")
  print(x$temperatures, digits = digits)
  cat("Noncentrality scale:", x$settings$measurement_scale, "\n")
  cat("Nuisance rank:", x$nuisance_rank, "\n")
  if(!is.null(x$panel)) cat("Panel:",x$panel$n_persons,"persons |",x$panel$n_observations,
    "person-occasion observations | person effects:",x$panel$person_effects,"\n")
  if(!is.null(x$reference)) {
    cat("Fixed reference:",if(length(unique(x$reference$value))==1L)
      x$reference$value[1] else "row-specific", "| no implicit centering\n")
    cat("Signed mean departure after nuisance adjustment:",format(x$mean_departure,digits=digits),"\n")
  }
  print(x$estimates, digits = digits, row.names = FALSE)
  if (!is.null(x$model)) {
    cat("Model:", if (is.null(x$model$formula)) "explicit X_H" else
        paste(deparse(x$model$formula), collapse = " "), "\n")
    if (identical(x$model$hypothesis$type, "conditional_block")) {
      h <- x$model$hypothesis
      cat("Conditional target:", paste(h$target, collapse = ", "),
          "| adjust:", paste(h$adjust, collapse = ", "), "\n")
      cat("Omitted terms:", paste(h$omitted, collapse = ", "), "\n")
      cat("lambda_star and kappa refer to structure after nuisance adjustment.\n")
    }
    if (identical(x$model$hypothesis$type, "linear_contrast")) {
      cat("Linear restrictions L beta = 0 (", x$model$hypothesis$kind, "):\n", sep = "")
      print(x$model$L, digits = digits)
      cat("lambda_star and kappa refer to structure after the restricted null model.\n")
    }
    p <- x$parsimony
    cat(if(is.null(x$reference)) "Centered rank:" else "Retained-space rank:", p$df_h, "| capacity:", p$eta_method, "\n")
    print(unlist(p[c("rss", "volume", "log_volume", "dimension_cost",
                     "geometry_cost", "free_energy")]), digits = digits)
  }
  .print_bootstrap_summary(x$bootstrap, digits)
  invisible(x)
}

#' @rdname lambda_results
#' @export
print.lambda_comparison_summary <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Lambda Star model comparison (point estimates)\n")
  cat("Estimator:", x$inference$estimator, "\n")
  cat("Status:", x$status, "| N =", length(x$cases$included),
      "| excluded =", length(x$cases$excluded), "\n")
  if (!is.null(x$message)) cat(x$message, "\n")
  if(!is.null(x$inference$reference)) {
    r <- x$inference$reference$value
    cat("Common fixed reference:",if(length(unique(r))==1L) r[1] else "row-specific","\n")
  }
  print(x$models, digits = digits, row.names = FALSE)
  cat("delta_F (b - a):", format(x$delta_F, digits = digits), "\n")
  .print_bootstrap_summary(x$bootstrap, digits)
  invisible(x)
}

.bootstrap_summary <- function(bootstrap) {
  if (is.null(bootstrap)) return(NULL)
  bootstrap[c("config", "ci", "conditional_quantiles", "diagnostics", "attempted", "skipped", "inference", "elapsed_seconds")]
}

.print_bootstrap_summary <- function(bootstrap, digits) {
  if (is.null(bootstrap)) return(invisible(NULL))
  cat("Pointwise percentile bootstrap intervals\n")
  cat("Estimator:", bootstrap$inference$estimator, "| resampling:",
      bootstrap$config$resampling, "| interpolation:", bootstrap$config$percentile_interpolation, "\n")
  cat("Level:", bootstrap$config$conf.level, "| draws:", bootstrap$attempted,
      "/", bootstrap$config$B, "| seed:", bootstrap$config$seed, "\n")
  print(bootstrap$ci[c("parameter", "estimate", "lower", "upper", "n_valid", "n_failed", "status")],
        digits = digits, row.names = FALSE)
  if (nrow(bootstrap$conditional_quantiles)) {
    cat("Quantiles conditional on valid draws (not regular confidence intervals):\n")
    print(bootstrap$conditional_quantiles, digits = digits, row.names = FALSE)
  }
  if (nrow(bootstrap$diagnostics)) print(bootstrap$diagnostics, row.names = FALSE)
  invisible(NULL)
}
