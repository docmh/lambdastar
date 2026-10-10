# Source this file; then call lambdastar_reference_example(). No work on source.
lambdastar_reference_example <- function(B=5000L,conf.level=.95,seed=1L) {
  had <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  old <- if(had)get(".Random.seed",envir=.GlobalEnv) else NULL;kind <- RNGkind()
  on.exit({do.call(RNGkind,as.list(kind));if(had)assign(".Random.seed",old,envir=.GlobalEnv)
    else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv)},add=TRUE)
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(68061)
  n <- 80L;k <- 4L
  u <- as.numeric(scale(rnorm(n)));v <- as.numeric(scale(rnorm(n)))
  baseline <- 52+1.5*u;change <- 1+.5*v
  d <- data.frame(id=seq_len(n));before <- paste0("pre",1:k);after <- paste0("post",1:k)
  for(j in seq_len(k)) {
    e1 <- rnorm(n);e2 <- rnorm(n)
    d[[before[j]]] <- baseline+e1
    d[[after[j]]] <- baseline+change+1.2*(.6*e1+sqrt(1-.6^2)*e2)
  }
  cfg <- lambdastar::lambda_bootstrap_control(B=B,conf.level=conf.level,seed=seed,
    min_valid=min(200L,B),keep_draws=TRUE,keep_indices=TRUE)
  reference <- lambdastar::lambda_reference(d,~1,value=50)
  one_sample <- lambdastar::lambda_model(d,reference,before,1/12,bootstrap=cfg)
  paired <- lambdastar::lambda_paired(d,before,after,1/12,value=0,bootstrap=cfg)
  paired_reference <- lambdastar::lambda_paired(d,before,after,1/12,value=.5)
  comparison <- lambdastar::lambda_compare(d,lambdastar::lambda_reference(d,~0,value=50),
    reference,before,1/12)
  truth <- function(S,r,T) {
    total <- sum((S-r)^2)/T;captured <- length(S)*(mean(S)-r)^2/T
    c(temperature=T,lambda_star=total,lambda_h=captured,lambda_perp=total-captured,
      kappa=if(total>0)captured/total else NA_real_)
  }
  list(data=d,before=before,after=after,baseline=baseline,change=change,
    design=reference,one_sample=one_sample,paired=paired,paired_reference=paired_reference,
    comparison=comparison,truth=list(one_sample=truth(baseline,50,1/k),
      paired=truth(change,0,(1+1.2^2-2*.6*1.2)/k),
      paired_reference=truth(change,.5,(1+1.2^2-2*.6*1.2)/k)),
    settings=list(B=B,conf.level=conf.level,seed=seed,T_q=1/12,
      pair_error_correlation=.6,difference="after - before"))
}
