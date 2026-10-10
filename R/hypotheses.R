#' Construct a conditional block from full-model terms
#'
#' Target and adjustment terms are explicit, disjoint selections from the full
#' model's term labels. Columns are taken from that model's actual coding.
#' Unselected terms are reported as omitted, although their variables still
#' participate in full-recipe evaluation and common case selection.
#'
#' The resulting lambda_star and kappa refer to structure remaining after the
#' specified adjustment. Partial kappas with different adjustment spaces are
#' not an additive partition. F_H comparisons require the same nuisance recipe.
#' This is a term-block constructor, not a general coefficient restriction or
#' an automatic Type I/II/III analysis. The intercept is always removed.
#' @inheritParams lambda_design
#' @param target Nonempty character vector of exact full-model term labels.
#' @param adjust Character vector of adjustment term labels. Explicitly use
#'   character() for no adjustment terms. No automatic complement is inferred.
#' @return A lambda_design for the target block with its effective nuisance recipe,
#'   target, adjustment and omitted-term metadata.
#' @examples
#' d <- data.frame(g=factor(rep(c("a","b"), each=10)), x=seq_len(20))
#' lambda_hypothesis(d, ~g*x, target="g:x", adjust=c("g","x"))
#' @export
lambda_hypothesis <- function(data, model, target, adjust, contrasts = NULL,
    nuisance = NULL, nuisance_contrasts = NULL, encoding = c("reevaluate", "fixed")) {
  design <- lambda_design(data, model, contrasts, nuisance, nuisance_contrasts,
                          encoding = match.arg(encoding))
  full <- design$spec
  if (full$type == "matrix")
    stop("Term hypotheses require a formula or lm; use lambda_design for explicit matrix blocks.", call. = FALSE)
  labels <- attr(full$terms, "term.labels")
  for (name in c("target", "adjust")) {
    value <- get(name)
    if (!is.character(value) || anyNA(value) || anyDuplicated(value) ||
        !all(value %in% labels) || (name == "target" && !length(value)))
      stop(name, " must select distinct exact full-model term labels: ",
        paste(labels, collapse = ", "), call. = FALSE)
  }
  if (length(intersect(target, adjust)))
    stop("target and adjust must be disjoint.", call. = FALSE)
  block <- function(selected) list(type = "block", full = full,
    selected = selected, variables = full$variables)
  b <- design$B
  if (length(adjust)) {
    adjustment <- block(adjust)
    b <- if (is.null(b)) adjustment else list(type = "combined", parts = list(b, adjustment),
      variables = unique(c(b$variables, adjustment$variables)))
  }
  .new_design(.lambda_data(data), block(target), b,
    list(type = "conditional_block", target = target, adjust = adjust,
      omitted = setdiff(labels, c(target, adjust)), full_terms = labels,
      total_structure = "after_nuisance_adjustment"))
}

.encode_block <- function(spec, prepared, tol) {
  full <- .encode_model(spec$full, prepared, tol)
  labels <- attr(full$metadata$terms, "term.labels")
  keep <- full$metadata$term_assign %in% match(spec$selected, labels)
  meta <- full$metadata
  meta$full_model <- full$metadata
  meta$type <- "block"
  meta$columns <- meta$columns[keep]
  meta$term_assign <- meta$term_assign[keep]
  meta$selected_terms <- spec$selected
  list(x = full$x[, keep, drop = FALSE], metadata = meta)
}
