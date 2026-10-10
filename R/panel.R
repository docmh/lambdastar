#' Encode fixed models on person-by-occasion observations
#'
#' Each row is one person-occasion observation and indicators are parallel
#' measurements within that row. Occasions are not treated as parallel repeats
#' of an unchanged response. The measurement bath must be homogeneous and
#' isotropic across the retained observation space, with independent parallel
#' replicate errors. Person bootstrap does not correct serial measurement error.
#'
#' Person intercepts are adjusted by default, using group demeaning without a
#' dense dummy matrix. Additional nuisance columns are adjusted in that space.
#' Set person_effects="none" to encode person effects explicitly in model.
#' The person column is treated as a categorical identity, including numeric IDs.
#' Formula columns involving person are rebuilt for distinct sampled copies;
#' fixed encoding of such formulas is rejected. Other fixed bases keep their
#' declared encoding. Explicit matrices encode row-bound non-person directions.
#'
#' With bootstrap=TRUE in lambda_model/lambda_compare, whole persons are sampled
#' with replacement. An explicit bootstrap control must request cluster resampling
#' with cluster=person. Strata must be constant within persons. Duplicate sampled
#' persons receive distinct draw-level identities. Unbalanced observation counts
#' are allowed, so the number of observation rows and retained df can vary.
#' na.action="complete" excludes an entire person if any required row is incomplete.
#' At least two persons, each with at least two observed occasions, are required.
#'
#' The default paper estimator, numerical F_H and percentile rules are retained.
#' This does not extend the optional population estimator. Comparisons require
#' the same observation map, person adjustment and nuisance specification.
#' @inheritParams lambda_design
#' @param model Formula or encoded matrix. NULL encodes a categorical occasion
#'   effect. For a numeric time trend supply an explicit formula.
#' @param person,occasion Single column names identifying persons and occasions.
#'   Each pair must be unique and nonmissing. Missing observation rows are allowed.
#' @param target Optional exact full-model term labels, as in lambda_hypothesis.
#' @param adjust Explicit adjustment terms when target is supplied.
#' @param person_effects "adjust" (default) removes person intercepts, or "none".
#' @return A row-bound lambda_design with a panel observation map. Evaluation
#'   returns panel metadata; retained bootstrap indices also have person_instances
#'   giving each sampled copy's draw-local ID in matching row order.
#' @examples
#' d <- data.frame(id=rep(1:12,each=3), time=rep(0:2,12))
#' lambda_panel(d, person="id", occasion="time")
#' lambda_panel(d, ~time, person="id", occasion="time")
#' lambda_panel(d, ~id*time, person="id", occasion="time",
#'   target="id:time", adjust="time")
#' @export
lambda_panel <- function(data, model=NULL, person, occasion, target=NULL,
    adjust=character(), person_effects=c("adjust","none"), contrasts=NULL,
    nuisance=NULL, nuisance_contrasts=NULL, encoding=c("reevaluate","fixed")) {
  original <- .lambda_data(data)
  map <- .panel_map(original,person,occasion,match.arg(person_effects))
  data <- .panel_data(original,map)
  if(is.null(model)) model <- stats::reformulate(paste0("factor(`",occasion,"`)"))
  if(inherits(model,c("lambda_design","lm")))
    stop("Panel construction requires an explicit formula or matrix; use target/adjust for term blocks.",call.=FALSE)
  if(is.null(target)&&length(adjust)) stop("adjust requires target.",call.=FALSE)
  design <- if(is.null(target)) lambda_design(data,model,contrasts,nuisance,
    nuisance_contrasts,match.arg(encoding)) else lambda_hypothesis(data,model,target,
      adjust,contrasts,nuisance,nuisance_contrasts,match.arg(encoding))
  design$spec <- .panel_spec(design$spec,person)
  design$B <- .panel_spec(design$B,person)
  variables <- unique(c(names(design$bound_data),person,occasion))
  design$bound_data <- original[variables]
  design$panel <- map
  design$resampling_policy <- "whole_persons_with_distinct_copy_identities"
  design
}

.panel_map <- function(data,person,occasion,person_effects) {
  person_effects <- match.arg(person_effects,c("adjust","none"))
  for(name in c("person","occasion")) {
    column <- get(name)
    if(!is.character(column)||length(column)!=1L||is.na(column)||
       !column%in%names(data)||grepl("`",column,fixed=TRUE))
      stop(name," must name one data column (without backticks).",call.=FALSE)
    x <- data[[column]]
    if(!is.atomic(x)||!is.null(dim(x))||anyNA(x)||
       (is.numeric(x)&&any(!is.finite(x))))
      stop("Person and occasion identities must be finite, nonmissing vector labels.",call.=FALSE)
  }
  if(person==occasion) stop("person and occasion must be different columns.",call.=FALSE)
  id <- as.character(data[[person]]);time <- as.character(data[[occasion]])
  if(anyDuplicated(data.frame(id,time)))
    stop("Each person-occasion pair must be unique; put parallel replicates in indicator columns.",call.=FALSE)
  counts <- table(id)
  if(length(counts)<2L||any(counts<2L))
    stop("A panel requires at least two persons, each with at least two occasions.",call.=FALSE)
  list(schema_version=1L,person=person,occasion=occasion,person_effects=person_effects,
    row_ids=rownames(data),person_ids=id,occasion_ids=time)
}

.panel_data <- function(data,panel) {
  if(is.null(panel)) return(data)
  id <- data[[panel$person]]
  levels <- if(is.factor(id)) levels(droplevels(id)) else sort(unique(as.character(id)))
  data[[panel$person]] <- factor(as.character(id),levels=levels)
  if(!is.null(attr(id,"contrasts",exact=TRUE)))
    attr(data[[panel$person]],"contrasts") <- attr(id,"contrasts",exact=TRUE)
  data
}

.panel_spec <- function(spec,person) {
  if(is.null(spec)) return(NULL)
  if(spec$type=="combined") spec$parts <- lapply(spec$parts,.panel_spec,person=person)
  if(spec$type=="block") spec$full <- .panel_spec(spec$full,person)
  if(person%in%spec$variables) {
    if(identical(spec$encoding,"fixed"))
      stop("Person-dependent recipes require encoding='reevaluate' for sampled copy identities.",call.=FALSE)
    spec$panel_person <- person
  }
  spec
}

.panel_candidates <- function(candidates) {
  maps <- lapply(candidates,function(c) if(inherits(c$model,"lambda_design")) c$model$panel else NULL)
  present <- !vapply(maps,is.null,logical(1))
  if(!any(present)) return(NULL)
  if(!all(present)||!all(vapply(maps,function(m)identical(m,maps[[1L]]),logical(1))))
    stop("Panel comparisons require the same person/occasion map and person adjustment in both designs.",call.=FALSE)
  maps[[1L]]
}

.panel_bootstrap <- function(bootstrap,models) {
  panel <- .panel_candidates(lapply(models,function(model)list(model=model)))
  if(!is.null(panel)&&identical(bootstrap,TRUE))
    bootstrap <- lambda_bootstrap_control(resampling="cluster",cluster=panel$person)
  config <- .bootstrap_control(bootstrap)
  if(!is.null(panel)&&!is.null(config)&&
     (config$resampling!="cluster"||config$cluster!=panel$person))
    stop("Panel bootstrap must resample whole persons: resampling='cluster', cluster=person.",call.=FALSE)
  config
}

.panel_nuisance <- function(x,data,panel) {
  if(is.null(panel)||panel$person_effects!="adjust") return(x)
  structure(list(x=x,person=as.character(data[[panel$person]])),class="lambda_person_nuisance")
}

.demean_person <- function(x,person) {
  vector <- is.null(dim(x))
  if(vector) x <- matrix(x,ncol=1L)
  if(!ncol(x)) return(x)
  means <- rowsum(x,person$group,reorder=FALSE)/person$counts
  result <- x-means[person$group,,drop=FALSE]
  if(vector) drop(result) else result
}

.person_nuisance_space <- function(nuisance,n,row_ids,control) {
  labels <- nuisance$person
  if(length(labels)!=n||anyNA(labels)) stop("Invalid person observation map.",call.=FALSE)
  group <- match(labels,unique(labels));person <- list(group=group,counts=tabulate(group))
  p <- length(person$counts);x <- nuisance$x;q <- NULL;rank <- p
  if(!is.null(x)&&ncol(x)) {
    x <- .matrix(x,"nuisance")
    if(nrow(x)!=n) stop("Nuisance and panel rows differ.",call.=FALSE)
    residual <- .demean_person(x,person)
    scale <- apply(abs(x),2,max)
    keep <- scale>0
    normalized <- sweep(residual[,keep,drop=FALSE],2,scale[keep],"/")
    norms <- sqrt(colSums(normalized^2))
    keep <- norms>control$rank_tol*sqrt(n)
    if(any(keep)) {
      normalized <- sweep(normalized[,keep,drop=FALSE],2,norms[keep],"/")
      q <- qr(normalized,tol=control$rank_tol,LAPACK=FALSE);rank <- rank+q$rank
    }
  }
  if(rank>=n) stop("Nuisance space leaves no residual measurement degrees of freedom.",call.=FALSE)
  list(rank=as.integer(rank),d=as.integer(n-rank),qr=q,person=person,
    columns=c("(Person intercepts)",colnames(x)),control=control)
}

.panel_result <- function(result,prepared) {
  panel <- prepared$panel
  if(is.null(panel)) return(result)
  data <- prepared$data
  result$panel <- c(panel[c("person","occasion","person_effects")],list(
    observation_unit="person_occasion",replicate_columns=prepared$measurement$settings$indicators,
    n_persons=length(unique(data[[panel$person]])),n_observations=nrow(data),
    map=data.frame(input_row=prepared$rows,row_id=rownames(data),
      person=as.character(data[[panel$person]]),occasion=as.character(data[[panel$occasion]])),
    missing_policy="whole_person_complete",bath="isotropic_on_retained_observation_space"))
  result$inference$target <- "panel_projection"
  result$inference$sampling_unit <- "person"
  result$inference$n_persons <- result$panel$n_persons
  result$settings$panel <- result$panel
  result
}
