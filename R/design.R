#' Bind a reusable encoded model specification to input rows
#'
#' A design keeps the actual coding and its recipe, not just a projection basis.
#' Formula realization is deferred until common case selection in lambda_model()
#' or lambda_compare(). Explicit matrices are already encoded. No temperature,
#' fitted projection or measurement state is cached in this object.
#'
#' Relevant input columns and row identities are bound exactly. New measurements
#' on the same predictors are allowed for one-sided formulas; changed predictors
#' or row order require a new design. Ordinary deterministic base/stats terms
#' and supported namespace-qualified terms are portable; caller-defined functions
#' must be evaluated explicitly into a matrix before construction.
#'
#' Factor contrast defaults are recorded and realized as numeric contrasts.
#' During bootstrap, explicit
#' matrices and simple formulas are sliced by row; sample-dependent formula
#' transforms are reevaluated by default, as with the existing formula interface.
#' With encoding="fixed", the raw encoding learned after common case selection
#' is retained during resampling. Projections are rebuilt in either case.
#' @param data Common input data frame or matrix, with unique row identities.
#' @param model Formula, compatible unweighted lm, or encoded numeric matrix.
#' @param contrasts Explicit formula contrast list, as in lambda_model().
#' @param nuisance Optional common one-sided nuisance formula or numeric matrix.
#' @param nuisance_contrasts Contrast list for the nuisance formula.
#' @param encoding Bootstrap basis policy. `"reevaluate"` (default) evaluates
#'   sample-dependent transforms in each draw. `"fixed"` learns the encoding
#'   after common case selection and resamples its rows, including nuisance
#'   terms. Persons are still resampled; this does not fix the observed X.
#' @return A versioned lambda_design specification. X_H is populated for matrix
#'   inputs; formula matrices are available as encoded_x on the evaluated result.
#'   B stores the nuisance recipe. row_ids and bound_data record the input binding.
#' @seealso [lambda_hypothesis()], [lambda_contrast()], [lambda_marginal_contrast()].
#'   The installed learning path is at
#'   `system.file("doc", "MODEL_SPACE_WORKFLOW.md", package="lambdastar")`.
#' @examples
#' d <- data.frame(g=factor(rep(c("a", "b"), each=10)), x=seq_len(20))
#' design <- lambda_design(d, ~g+x)
#' design
#' @export
lambda_design <- function(data, model, contrasts = NULL, nuisance = NULL,
                          nuisance_contrasts = NULL,
                          encoding = c("reevaluate", "fixed")) {
  encoding <- match.arg(encoding)
  data <- .lambda_data(data)
  if (inherits(model, "lambda_design"))
    stop("model is already a lambda_design; use it directly or construct a new specification.", call. = FALSE)
  if (inherits(nuisance, "lambda_design"))
    stop("nuisance must be a one-sided formula or numeric matrix, not a lambda_design.", call. = FALSE)
  spec <- .portable_spec(model, contrasts, data)
  b <- if (is.null(nuisance)) NULL else .portable_spec(nuisance, nuisance_contrasts, data)
  spec$encoding <- encoding
  if (!is.null(b)) b$encoding <- encoding
  if (is.null(nuisance) && !is.null(nuisance_contrasts))
    stop("nuisance_contrasts requires nuisance.", call. = FALSE)
  .check_nuisance_spec(b)
  .new_design(data, spec, b, list(type = "whole_model"))
}

.new_design <- function(data, spec, b, hypothesis) {
  variables <- unique(c(spec$variables, b$variables))
  structure(list(schema_version = 1L, spec = spec, B = b,
    X_H = if (spec$type == "matrix") spec$x else NULL,
    row_ids = rownames(data), bound_data = data[variables],
    hypothesis = hypothesis,
    resampling_policy = if (.simple_spec(spec)) "subset_rows" else "reevaluate_recipe"),
    class = "lambda_design")
}

.check_nuisance_spec <- function(spec) {
  if (!is.null(spec) && (spec$type == "lm" ||
      (!is.null(spec$terms) && attr(spec$terms, "response") != 0L)))
    stop("nuisance must be a one-sided formula or numeric matrix.", call. = FALSE)
}

.portable_formula <- function(formula) {
  allowed <- c("~", "+", "-", "*", "/", ":", "^", "(", "I", "log", "log1p",
    "log2", "log10", "sqrt", "exp", "abs", "sin", "cos", "tan", "pmax", "pmin",
    "scale", "poly", "factor", "as.factor", "ordered", "as.numeric", "as.integer", "cut",
    "relevel", "interaction", "mean", "sd", "var", "min", "max", "range", "median",
    "quantile", "C", "ns", "bs", "c", "list")
  walk <- function(x) {
    if (!is.call(x)) return(invisible(NULL))
    head <- x[[1L]]
    if (is.call(head) && identical(head[[1L]], as.name("::"))) {
      pkg <- as.character(head[[2L]]); fun <- as.character(head[[3L]])
      if (!pkg %in% c("base", "stats", "splines") || !fun %in% allowed)
        stop("Design recipes require supported deterministic base/stats/splines terms; supply an encoded matrix.", call. = FALSE)
    } else if (!is.symbol(head) || !as.character(head) %in% setdiff(allowed,c("ns","bs")))
      stop("Design recipes require supported deterministic terms; supply an encoded matrix for custom functions.", call. = FALSE)
    for (part in as.list(x)[-1L]) walk(part)
  }
  walk(formula)
  environment(formula) <- asNamespace("stats")
  formula
}

.portable_spec <- function(model, contrasts, data) {
  spec <- .model_spec(model, contrasts, data)
  if (spec$type == "matrix") return(spec)
  spec$formula <- .portable_formula(spec$formula)
  spec$terms <- stats::terms(spec$formula, data = data)
  # Resolve factor coding without evaluating transformations or selecting cases.
  # Transformed factor terms are resolved when their model frame is realized.
  spec$contrast_defaults <- getOption("contrasts")
  for (choice in c(as.list(spec$contrast_defaults), spec$contrasts))
    .check_portable_contrast(choice)
  for (column in data[spec$variables]) if (is.factor(column) && !is.null(attr(column, "contrasts")))
    .check_portable_contrast(attr(column, "contrasts"))
  if (!is.null(spec$fit)) {
    fit <- spec$fit
    attr(fit$model, "terms") <- spec$terms
    spec$fit <- structure(list(model = fit$model, terms = spec$terms,
      contrasts = fit$contrasts, xlevels = fit$xlevels), class = "lm")
  }
  spec
}

.check_portable_contrast <- function(choice) {
  if (is.matrix(choice) && is.numeric(choice) && all(is.finite(choice))) return(invisible(NULL))
  if (is.character(choice) && length(choice) == 1L && !is.na(choice) &&
      choice %in% c("contr.treatment", "contr.sum", "contr.helmert", "contr.poly", "contr.SAS"))
    return(invisible(NULL))
  stop("Portable designs require numeric contrast matrices or standard stats contrast names.", call. = FALSE)
}

.bound_design <- function(object, data) {
  if (!identical(object$schema_version, 1L) || !is.list(object$spec))
    stop("Unsupported lambda_design schema.", call. = FALSE)
  if (!identical(object$row_ids, rownames(data)))
    stop("Design row identities/order changed; construct a new lambda_design.", call. = FALSE)
  if(!is.null(object$reference)) {
    if(!identical(object$reference$type,"fixed")||is.null(object$reference$value)||
       !identical(object$reference$row_ids,rownames(data)))
      stop("Reference specification or row identities changed; construct a new lambda_reference.",call.=FALSE)
    .reference_values(object$reference$value,nrow(data),rownames(data))
  }
  if(!is.null(object$panel)) {
    expected <- .panel_map(data,object$panel$person,object$panel$occasion,object$panel$person_effects)
    if(!identical(object$panel,expected)||!is.null(object$reference))
      stop("Panel observation map changed or conflicts with a fixed reference; construct a new lambda_panel.",call.=FALSE)
  }
  vars <- names(object$bound_data)
  if (!all(vars %in% names(data)) || !identical(object$bound_data, data[vars]))
    stop("Design input columns changed; construct a new lambda_design.", call. = FALSE)
  invisible(object)
}

.spec_signature <- function(spec) {
  if (is.null(spec)) return(NULL)
  if (spec$type == "matrix") return(list(type = "matrix", x = spec$x))
  if (spec$type == "block") return(list(type = "block", full = .spec_signature(spec$full),
    selected = spec$selected))
  if (spec$type == "combined") return(lapply(spec$parts, .spec_signature))
  if (spec$type == "contrast") return(list(type = "contrast", role = spec$role,
    full = .spec_signature(spec$full), base = .spec_signature(spec$base),
    restriction = spec$restriction))
  list(type = "formula", formula = paste(deparse(spec$formula), collapse = ""),
    contrasts = spec$contrasts, defaults = spec$contrast_defaults,
    encoding = if (is.null(spec$encoding)) "reevaluate" else spec$encoding)
}

.simple_spec <- function(spec) {
  if(!is.null(spec$panel_person)) return(FALSE)
  if (spec$type == "matrix") return(TRUE)
  if (spec$type == "block") return(.simple_spec(spec$full))
  if (spec$type == "combined") return(all(vapply(spec$parts, .simple_spec, logical(1))))
  if (spec$type == "contrast") return(FALSE)
  if (identical(spec$encoding, "fixed")) return(TRUE)
  all(all.names(spec$formula, functions = TRUE) %in%
    c(all.vars(spec$formula), "~", "+", "-", "*", "/", ":", "^", "("))
}

.spec_complete <- function(spec, n) {
  if (is.null(spec)) return(rep(TRUE, n))
  if (spec$type == "matrix") return(stats::complete.cases(spec$x))
  if (spec$type == "block") return(.spec_complete(spec$full, n))
  if (spec$type == "combined")
    return(Reduce(`&`, lapply(spec$parts, .spec_complete, n = n)))
  if (spec$type == "contrast")
    return(.spec_complete(spec$full, n) & .spec_complete(spec$base, n))
  rep(TRUE, n)
}

#' @rdname lambda_design
#' @param x A lambda_design object.
#' @param ... Reserved for S3 compatibility.
#' @export
print.lambda_design <- function(x, ...) {
  cat("Lambda Star design | schema", x$schema_version, "|", length(x$row_ids), "input rows\n")
  cat("Hypothesis:", x$hypothesis$type, "| bootstrap encoding:", x$resampling_policy, "\n")
  if(!is.null(x$panel)) cat("Panel:",x$panel$person,"x",x$panel$occasion,
    "| person effects:",x$panel$person_effects,"\n")
  if (x$hypothesis$type == "conditional_block") {
    cat("Target:", paste(x$hypothesis$target, collapse = ", "), "\n")
    cat("Adjust:", paste(x$hypothesis$adjust, collapse = ", "), "\n")
    cat("Omitted:", paste(x$hypothesis$omitted, collapse = ", "), "\n")
  }
  if (x$hypothesis$type == "linear_contrast")
    cat("Restriction:", x$hypothesis$kind, "| source encoding:", x$hypothesis$source_encoding,
        "| reference: restricted null model\n")
  if(!is.null(x$panel)) {
    cat("Nuisance:",if(x$panel$person_effects=="adjust") "person intercepts" else "global intercept",
      if(is.null(x$B)) "" else "+ explicit columns","\n")
  } else if(!is.null(x$reference)) {
    cat("Reference: fixed | no implicit centering | values:",
        if(length(unique(x$reference$value))==1) x$reference$value[1] else "row-specific", "\n")
    cat("Nuisance:",if(is.null(x$B)) "none" else "explicit; no implicit intercept","\n")
  } else cat("Nuisance:", if (is.null(x$B)) "intercept" else "explicit (plus intercept)", "\n")
  invisible(x)
}
