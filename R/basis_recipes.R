# Prediction terms are learned only after common selection. makepredictcall()
# records poly/ns/bs parameters; scale() needs its center/scale added explicitly.
.grid_recipe <- function(terms, frame) {
  out <- stats::delete.response(terms)
  expressions <- as.list(attr(out,"predvars"))
  if (!length(expressions)) expressions <- as.list(attr(out,"variables"))
  frames <- frame
  if(attr(terms,"response")) frames <- frames[-1L]
  function_name <- function(call) {
    head <- call[[1L]]
    if(is.call(head) && identical(head[[1L]],as.name("::"))) as.character(head[[3L]]) else
      if(is.symbol(head)) as.character(head) else "unsupported"
  }
  for(i in seq_along(frames)) {
    call <- expressions[[i+1L]]
    if(is.call(call) && function_name(call)=="scale") {
      # Normalize positional arguments before replacing fitted parameters.
      call <- match.call(base::scale.default,call,expand.dots=TRUE)
      call[[1L]] <- as.name("scale")
      center <- attr(frames[[i]],"scaled:center");scale <- attr(frames[[i]],"scaled:scale")
      call$center <- if(is.null(center)) FALSE else center
      call$scale <- if(is.null(scale)) FALSE else scale
      expressions[[i+1L]] <- call
    }
  }
  unsafe <- function(x) {
    if(!is.call(x)) return(FALSE)
    name <- function_name(x)
    if(name %in% c("mean","sd","var","min","max","range","median","quantile","cut")) return(TRUE)
    if(name=="scale" && (is.null(x$center)||is.null(x$scale))) return(TRUE)
    if(name=="poly" && is.null(x$coefs) && !identical(x$raw,TRUE)) return(TRUE)
    if(name %in% c("ns","bs") && (is.null(x$knots)||is.null(x$Boundary.knots))) return(TRUE)
    any(vapply(as.list(x)[-1L],unsafe,logical(1)))
  }
  if(any(vapply(expressions[-1L],unsafe,logical(1))))
    return(list(terms=NULL,reason="Grid contrasts need pointwise terms or supported poly/spline/scale prediction recipes; use explicit L for other transforms."))
  attr(out,"predvars") <- as.call(expressions)
  list(terms=out,reason=NULL)
}

.grid_encoding <- function(meta, grid) {
  recipe <- meta$grid_recipe
  if(is.null(recipe$terms)) stop(if(is.null(recipe$reason))
    "The source does not provide a predictor-grid recipe." else recipe$reason,call.=FALSE)
  needed <- meta$predictor_variables
  if(!all(needed %in% names(grid)))
    stop("grid must supply all full-model predictors: ",paste(setdiff(needed,names(grid)),collapse=", "),call.=FALSE)
  frame <- stats::model.frame(recipe$terms,data=grid,na.action=stats::na.fail,
    xlev=meta$factor_levels[setdiff(names(meta$factor_levels),meta$logical_factors)],
    drop.unused.levels=FALSE)
  mm <- stats::model.matrix(recipe$terms,frame,contrasts.arg=meta$contrasts)
  x <- mm[,attr(mm,"assign")!=0L,drop=FALSE]
  if(!identical(colnames(x),meta$columns)) stop("Grid coding differs from the full model.",call.=FALSE)
  .finite(x)
}
