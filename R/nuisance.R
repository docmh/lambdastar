# Ordinary nuisance spaces include the intercept. Fixed-reference spaces use
# the explicitly supplied columns only. QR rank determines residual df.
.nuisance_space <- function(nuisance, n, row_ids = NULL, control = lambda_control(), center = TRUE) {
  control <- .control(control)
  if(inherits(nuisance,"lambda_person_nuisance"))
    return(.person_nuisance_space(nuisance,n,row_ids,control))
  if (!center && is.null(nuisance)) {
    return(list(rank=0L,d=n,qr=NULL,columns=character(),control=control,center=FALSE))
  }
  if (is.null(nuisance)) return(list(rank = 1L, d = n - 1L, qr = NULL,
    columns = "(Intercept)", control = control))
  nuisance <- .matrix(nuisance, "nuisance")
  if (nrow(nuisance) != n) stop("nuisance and measurements have incompatible rows.", call. = FALSE)
  if (!is.null(row_ids) && !is.null(rownames(nuisance)) &&
      !identical(as.character(row_ids), rownames(nuisance)))
    stop("nuisance and measurements have different row identities or order.", call. = FALSE)
  if(!center&&!ncol(nuisance)) return(list(rank=0L,d=n,qr=NULL,columns=character(),control=control,center=FALSE))
  b <- if (center) cbind(`(Intercept)` = 1, center_columns(nuisance)) else nuisance
  scales <- apply(abs(b), 2L, max)
  keep <- scales > 0
  if (!any(keep)) return(list(rank=0L,d=n,qr=NULL,columns=colnames(b),control=control,center=FALSE))
  normalized <- sweep(b[, keep, drop = FALSE], 2L, scales[keep], "/")
  normalized <- sweep(normalized, 2L, sqrt(colSums(normalized^2)), "/")
  q <- qr(normalized, tol = control$rank_tol, LAPACK = FALSE)
  if (q$rank >= n) stop("Nuisance space leaves no residual measurement degrees of freedom.", call. = FALSE)
  out <- list(rank = q$rank, d = n - q$rank, qr = q, columns = colnames(b), control = control)
  if (!center) out$center <- FALSE
  out
}

.residualize_nuisance <- function(x, space) {
  if(!is.null(space$person)) {
    x <- .demean_person(x,space$person)
    if(is.null(space$qr)||(!is.null(dim(x))&&!ncol(x))) return(.finite(x))
    return(.finite(qr.resid(space$qr,x)))
  }
  if (is.null(space$qr)) {
    if (identical(space$center,FALSE)) return(.finite(x))
    if (is.null(dim(x))) return(.finite(x - mean(x)))
    return(center_columns(x))
  }
  if (!is.null(dim(x)) && !ncol(x)) return(x)
  .finite(qr.resid(space$qr, x))
}

.conditional_basis <- function(x, space, control = lambda_control()) {
  x <- .matrix(x, "x")
  if (identical(space$center,FALSE)) {
    adjusted <- .residualize_nuisance(x,space)
    for(j in seq_len(ncol(x))) {
      scale <- max(abs(x[,j]))
      if(scale==0 || sqrt(sum((adjusted[,j]/scale)^2)) <=
         control$rank_tol*sqrt(sum((x[,j]/scale)^2)))
        stop("x has a constant or redundant column after nuisance adjustment.",call.=FALSE)
    }
    return(.encoded_basis(adjusted,.control(control),space$d))
  }
  if (is.null(space$qr)&&is.null(space$person)) return(hypothesis_basis(x, control))
  adjusted <- .residualize_nuisance(x, space)
  original <- center_columns(x)
  for (j in seq_len(ncol(x))) {
    scale <- max(abs(original[, j]))
    if (scale == 0 || sqrt(sum((adjusted[, j]/scale)^2)) <=
        control$rank_tol * sqrt(sum((original[, j]/scale)^2)))
      stop("x has a constant or redundant column after nuisance adjustment.", call. = FALSE)
  }
  basis <- hypothesis_basis(adjusted, control)
  if (basis$rank > space$d) stop("x is rank deficient in the nuisance residual space.", call. = FALSE)
  basis
}
