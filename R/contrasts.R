#' Encode a homogeneous linear coefficient hypothesis
#'
#' Constructs L beta = 0 in the non-intercept coefficient coding of a full
#' model. Named columns not supplied in L receive zero weights. The full design
#' must be identifiable after common nuisance adjustment, and L must have
#' independent rows. Intercept/nonzero-reference hypotheses require a different
#' reference space and are not part of this constructor.
#'
#' The null-model directions become additional nuisance columns. Consequently,
#' lambda_star and kappa describe structure after adjustment for that null model.
#' The actual ordered contrast coding is retained for F_H, not replaced by the
#' numerical orthogonal basis. Geometry and estimability are rebuilt each draw.
#' @inheritParams lambda_design
#' @param model A whole-model lambda_design, formula, compatible unweighted lm,
#'   or encoded numeric matrix. A supplied design owns its coding, nuisance and
#'   encoding policy; contradictory constructor arguments are rejected.
#' @param L Finite numeric matrix with unique coefficient column names and
#'   optional unique row labels, or a named numeric vector for one restriction.
#' @return A lambda_design for the conditional contrast space. The realized L,
#'   coefficient right inverse, null basis and full-model metadata are available
#'   in the evaluated result's model component.
#' @examples
#' d <- data.frame(x=seq_len(20), z=sin(seq_len(20)))
#' lambda_contrast(d, ~x+z, c(x=1, z=-1))
#' @export
lambda_contrast <- function(data, model, L, contrasts = NULL, nuisance = NULL,
    nuisance_contrasts = NULL, encoding = c("reevaluate", "fixed")) {
  source <- .contrast_source(data, model, contrasts, nuisance, nuisance_contrasts,
                             if (missing(encoding)) NULL else match.arg(encoding))
  L <- .named_restrictions(L, "L")
  if ("(Intercept)" %in% colnames(L))
    stop("L must use non-intercept coefficients; reference-mean hypotheses are not supported.", call. = FALSE)
  .new_contrast(source, .lambda_data(data), list(kind="coefficients", L=L))
}

#' Encode contrasts of explicitly specified predictor settings
#'
#' The grid lists predictor settings; each weights row is a signed contrast
#' over those grid rows and must sum to zero. Supply all averaging weights
#' explicitly (for example -0.5,-0.5,0.5,0.5 for an equally averaged difference).
#' No covariate values, reference population or equal-cell averaging is inferred.
#' The learned polynomial/spline/scale prediction recipe is used to encode the
#' grid in the full model's coordinates, separately in each bootstrap draw.
#' Unsupported sample-dependent transforms require an explicit coefficient L.
#' @inheritParams lambda_contrast
#' @param grid Data frame of at least two complete predictor settings, with
#'   unique row names. It does not contain measurement outcomes.
#' @param weights Named numeric vector or matrix whose columns match the grid
#'   row names exactly. Row order/labels define the contrast coding.
#' @return A lambda_design as for lambda_contrast(), with grid and weights in
#'   its hypothesis metadata. The realized coefficient L is stored on results.
#' @examples
#' d <- data.frame(g=factor(rep(c("a","b"),10)), x=seq_len(20))
#' grid <- data.frame(g=c("a","b"), x=c(10,10), row.names=c("a10","b10"))
#' lambda_marginal_contrast(d, ~g*x, grid, c(a10=-1,b10=1))
#' @export
lambda_marginal_contrast <- function(data, model, grid, weights, contrasts = NULL,
    nuisance = NULL, nuisance_contrasts = NULL, encoding = c("reevaluate", "fixed")) {
  source <- .contrast_source(data, model, contrasts, nuisance, nuisance_contrasts,
                             if (missing(encoding)) NULL else match.arg(encoding))
  if (source$spec$type == "matrix")
    stop("A predictor grid needs a formula/lm recipe; use lambda_contrast for a matrix design.", call. = FALSE)
  grid <- .lambda_data(grid)
  if (any(!stats::complete.cases(grid))) stop("grid must be complete.", call. = FALSE)
  weights <- .named_restrictions(weights, "weights")
  if (!setequal(colnames(weights), rownames(grid)))
    stop("weights columns must match every grid row name exactly.", call. = FALSE)
  weights <- weights[,rownames(grid),drop=FALSE]
  if (any(abs(rowSums(weights)) > 1e-12*rowSums(abs(weights))))
    stop("Each weights row must sum to zero; mean-reference hypotheses are not supported.", call. = FALSE)
  .new_contrast(source, .lambda_data(data), list(kind="grid", grid=grid, weights=weights))
}

.contrast_source <- function(data, model, contrasts, nuisance, nuisance_contrasts, encoding) {
  if(inherits(model,"lambda_design")&&!is.null(model$panel))
    stop("Panel term contrasts must be declared via lambda_panel target/adjust or explicit row coding.",call.=FALSE)
  if(inherits(model,"lambda_design")&&!is.null(model$reference))
    stop("Centered contrast constructors cannot use a fixed reference; supply explicit affine X, reference and nuisance coding via lambda_reference.",call.=FALSE)
  if (!inherits(model,"lambda_design")) return(lambda_design(data,model,contrasts,
    nuisance,nuisance_contrasts,encoding=if(is.null(encoding)) "reevaluate" else encoding))
  .bound_design(model,.lambda_data(data))
  if (!identical(model$hypothesis$type,"whole_model"))
    stop("Use a whole-model lambda_design as the contrast source.",call.=FALSE)
  if (!is.null(contrasts) || !is.null(nuisance) || !is.null(nuisance_contrasts) ||
      (!is.null(encoding) && !identical(encoding,model$spec$encoding)))
    stop("A supplied design already specifies coding, nuisance and encoding policy.",call.=FALSE)
  model
}

.named_restrictions <- function(x, name) {
  if (is.numeric(x) && is.null(dim(x)))
    x <- matrix(x,1L,dimnames=list(NULL,names(x)))
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) || !nrow(x) || !ncol(x) ||
      any(!is.finite(x)) || is.null(colnames(x)) || anyNA(colnames(x)) ||
      any(!nzchar(colnames(x))) || anyDuplicated(colnames(x)))
    stop(name," must be a finite numeric matrix/vector with unique named columns.",call.=FALSE)
  if (is.null(rownames(x))) rownames(x) <- paste0("contrast",seq_len(nrow(x)))
  if (anyNA(rownames(x)) || any(!nzchar(rownames(x))) || anyDuplicated(rownames(x)))
    stop(name," must have unique nonempty row labels.",call.=FALSE)
  if (any(rowSums(abs(x))==0)) stop(name," has a zero restriction row.",call.=FALSE)
  x
}

.new_contrast <- function(source, data, restriction) {
  common <- list(type="contrast",full=source$spec,base=source$B,restriction=restriction,
    variables=unique(c(source$spec$variables,source$B$variables)))
  target <- common;target$role <- "target"
  null <- common;null$role <- "null"
  .new_design(data,target,null,c(list(type="linear_contrast",
    total_structure="after_restricted_null_model",source_encoding=source$resampling_policy),restriction))
}

.align_L <- function(L, columns) {
  if (anyNA(columns) || anyDuplicated(columns)) stop("Full model needs unique coefficient names.",call.=FALSE)
  if (!all(colnames(L) %in% columns))
    stop("Unknown L coefficient names: ",paste(setdiff(colnames(L),columns),collapse=", "),call.=FALSE)
  out <- matrix(0,nrow(L),length(columns),dimnames=list(rownames(L),columns))
  out[,colnames(L)] <- L
  out
}

# Stable coefficient-space solves: L D = I and L N = 0. No normal equations.
.restriction_basis <- function(L, control) {
  r <- nrow(L);p <- ncol(L)
  if (r>p) stop("L restriction rows are rank deficient.",call.=FALSE)
  maxima <- apply(abs(L),1L,max)
  if (any(maxima==0)) stop("L has a zero/nonestimable restriction row.",call.=FALSE)
  normalized <- sweep(L,1L,maxima,"/")
  norms <- sqrt(rowSums(normalized^2))
  normalized <- sweep(normalized,1L,norms,"/")
  singular <- svd(normalized,nu=0,nv=0)$d
  if (min(singular)<=control$rank_tol*max(singular))
    stop("L restriction rows are rank deficient at the requested tolerance.",call.=FALSE)
  q <- qr(t(normalized),tol=control$rank_tol,LAPACK=FALSE)
  if (q$rank!=r) stop("L restriction rows are rank deficient.",call.=FALSE)
  Q <- qr.Q(q,complete=TRUE)
  D <- matrix(0,p,r)
  D[,q$pivot] <- Q[,seq_len(r),drop=FALSE] %*% forwardsolve(t(qr.R(q)),diag(r))
  D <- .finite(sweep(sweep(D,2L,norms,"/"),2L,maxima,"/"))
  N <- Q[,seq_len(p)>r,drop=FALSE]
  dimnames(D) <- list(colnames(L),rownames(L))
  dimnames(N) <- list(colnames(L),if(ncol(N)) paste0("null",seq_len(ncol(N))) else character())
  list(D=D,N=N,rank=r)
}

.contrast_components <- function(spec, prepared, tol) {
  key <- spec;key[c("role","hypothesis","design_policy")] <- NULL
  cache <- prepared$encoding_cache
  if (!is.null(cache)) for(entry in cache$contrasts)
    if (identical(entry$key,key) && identical(entry$tol,tol)) return(entry$value)
  full <- .encode_model(spec$full,prepared,tol)
  base <- if(is.null(spec$base)) NULL else .encode_model(spec$base,prepared,tol)
  control <- prepared$numerical_control
  if(is.null(control)) control <- lambda_control()
  control <- .control(control)
  space <- .nuisance_space(base$x,nrow(full$x),rownames(full$x),control)
  # Validate identifiability of the entire conditional coefficient system.
  .conditional_basis(full$x,space,control)
  L <- if(spec$restriction$kind=="coefficients")
    .align_L(spec$restriction$L,colnames(full$x)) else {
      grid <- .grid_encoding(full$metadata,spec$restriction$grid)
      .finite(spec$restriction$weights %*% grid)
    }
  basis <- .restriction_basis(L,control)
  target <- .finite(full$x %*% basis$D)
  null <- .finite(full$x %*% basis$N)
  if(!is.null(base)) null <- cbind(base$x,null)
  rownames(target) <- rownames(null) <- rownames(full$x)
  meta <- list(type="contrast",full_model=full$metadata,base_model=base$metadata,
    L=L,coefficient_right_inverse=basis$D,coefficient_null_basis=basis$N,
    restriction=spec$restriction)
  value <- list(target=target,null=null,metadata=meta)
  if(!is.null(cache)) cache$contrasts <- c(cache$contrasts,list(list(key=key,tol=tol,value=value)))
  value
}

.encode_contrast <- function(spec, prepared, tol) {
  value <- .contrast_components(spec,prepared,tol)
  x <- value[[spec$role]]
  meta <- value$metadata;meta$columns <- colnames(x);meta$role <- spec$role
  list(x=x,metadata=meta)
}

.same_reference <- function(a,b,control) {
  if (identical(.spec_signature(a),.spec_signature(b))) return(TRUE)
  if (is.null(a)||is.null(b)||a$type!="contrast"||b$type!="contrast"||
      a$role!="null"||b$role!="null") return(FALSE)
  if (!identical(.spec_signature(a$full),.spec_signature(b$full)) ||
      !.same_reference(a$base,b$base,control) ||
      a$restriction$kind!=b$restriction$kind) return(FALSE)
  if(a$restriction$kind=="grid") {
    if(!identical(a$restriction$grid,b$restriction$grid)) return(FALSE)
    L1 <- a$restriction$weights;L2 <- b$restriction$weights
  } else {L1 <- a$restriction$L;L2 <- b$restriction$L}
  columns <- union(colnames(L1),colnames(L2))
  L1 <- .align_L(L1,columns);L2 <- .align_L(L2,columns)
  q1 <- .restriction_basis(L1,control);q2 <- .restriction_basis(L2,control)
  if(q1$rank!=q2$rank) return(FALSE)
  # Equal row spaces imply the same coefficient null model in every draw.
  scaled <- sweep(L2,1L,apply(abs(L2),1L,max),"/")
  max(c(0,abs(scaled %*% q1$N))) <= control$rank_tol
}
