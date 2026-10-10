#' Encode a hypothesis relative to a known reference
#'
#' The mean direction is retained: ~1 is a free mean and ~0 is the fixed null
#' value. This extends the centered manuscript calculation to a prespecified
#' affine origin. Formula intercepts are retained, including in an explicitly
#' supplied nuisance formula; use ~0+x for nuisance x without an intercept.
#' Reference values must be known independently of this analysis.
#'
#' All compared models must declare the same reference and nuisance recipe.
#' The ordinary centered API is unchanged. This path uses estimator="paper";
#' it does not extend the population estimator or the centered L beta=0 helpers.
#' See the installed derivation at
#' `system.file("doc","REFERENCE_MEANS.md",package="lambdastar")`.
#' @inheritParams lambda_design
#' @param model Formula or explicit encoded matrix. Formula intercepts are
#'   retained. lm sources are not accepted; declare the reference recipe explicitly.
#' @param value Known scalar reference, or finite numeric vector in input-row
#'   order. If the vector is named, its names must match the row identities.
#' @return A row-bound lambda_design with an explicit fixed reference. After
#'   evaluation, reference, mean_departure and retained df describe the origin
#'   and signed mean departure after any nuisance adjustment. Lambda values
#'   quantify squared departures; F_H uses the uncentered adjusted coding.
#' @examples
#' d <- data.frame(x=seq_len(20))
#' lambda_reference(d, ~1, value=50)
#' lambda_reference(d, ~0, value=50)
#' @export
lambda_reference <- function(data, model=~1, value=0, nuisance=NULL,
    contrasts=NULL, nuisance_contrasts=NULL, encoding=c("reevaluate","fixed")) {
  data <- .lambda_data(data)
  if(is.null(value)) stop("value must specify a known reference; use lambda_design for centering.",call.=FALSE)
  if(inherits(model,"lm")) stop("Reference designs require an explicit formula or matrix, not lm.",call.=FALSE)
  origin <- .reference_values(value,nrow(data),rownames(data))
  design <- lambda_design(data,model,contrasts,nuisance,nuisance_contrasts,match.arg(encoding))
  design$reference <- list(type="fixed",value=origin,row_ids=rownames(data),implicit_intercept=FALSE)
  design
}

.reference_values <- function(value,n,row_ids=NULL) {
  if(is.null(value)) return(NULL)
  if(!is.numeric(value)||is.complex(value)||!is.null(dim(value))||
     !length(value)||any(!is.finite(value))||!length(value)%in%c(1L,n))
    stop("reference must be a finite numeric scalar or one value per input row.",call.=FALSE)
  if(length(value)==n&&!is.null(names(value))&&!is.null(row_ids)&&
     !identical(names(value),as.character(row_ids)))
    stop("reference names must match input row identities and order.",call.=FALSE)
  rep(unname(value),length.out=n)
}

#' Evaluate paired mean differences from parallel measurement pairs
#'
#' Forms after-before within each person and each matching replicate pair,
#' then calls lambda_model with a lambda_reference design. Temperature is
#' estimated from the difference items, retaining homogeneous within-pair
#' error covariance. Different replicate-pair errors must be independent.
#' Persons and both occasions travel together in case or cluster resampling.
#' @inheritParams lambda_reference
#' @param before,after Character vectors of matching scored measurement columns,
#'   in pair order. Each has the same length K >= 2; columns cannot be reused.
#' @param T_q Positive resolution variance explicitly on the composite-mean
#'   difference scale. It is not inferred or doubled from occasion resolutions.
#' @param ... Other lambda_model settings, including bootstrap, na.action,
#'   measurement_scale, capacity and numerical_control.
#' @return A lambda_model plus paired metadata: column pairs, after-minus-before
#'   coefficients and the difference matrix in original input-row order. The
#'   default ~1 estimates captured mean departure from value. T=0 is singular.
#' @export
lambda_paired <- function(data, before, after, T_q, value=0, model=~1,
    nuisance=NULL, contrasts=NULL, nuisance_contrasts=NULL,
    encoding=c("reevaluate","fixed"), ...) {
  data <- .lambda_data(data)
  columns <- c(before,after)
  if(!is.character(before)||!is.character(after)||length(before)<2L||
     length(before)!=length(after)||anyNA(columns)||anyDuplicated(columns)||
     !all(columns%in%names(data)))
    stop("before and after must specify at least two distinct, matching measurement pairs.",call.=FALSE)
  if(!all(vapply(data[columns],function(x)is.numeric(x)&&!is.complex(x)&&is.null(dim(x)),logical(1))))
    stop("Paired measurements must be numeric scored columns.",call.=FALSE)
  differences <- as.matrix(data[after])-as.matrix(data[before])
  names <- make.unique(c(colnames(data),paste0(".lambda_difference_",seq_along(before))))[
    ncol(data)+seq_along(before)]
  colnames(differences) <- names
  data[names] <- as.data.frame(differences)
  design <- lambda_reference(data,model,value,nuisance,contrasts,nuisance_contrasts,match.arg(encoding))
  result <- lambda_model(data,design,names,T_q,...)
  result$paired <- list(before=before,after=after,difference_columns=names,
    coefficients=c(before=-1,after=1),differences=differences,
    temperature="estimated_from_parallel_differences",
    resolution="supplied_on_composite_mean_difference_scale")
  result
}
