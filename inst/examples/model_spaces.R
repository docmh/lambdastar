# Source this file, then call lambdastar_model_spaces_example(B=5000).
lambdastar_model_spaces_example <- function(B=5000L) {
  old <- if (exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) .Random.seed else NULL
  kind <- RNGkind()
  on.exit({do.call(RNGkind,as.list(kind));if(is.null(old)) {
    if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)
  } else assign(".Random.seed",old,envir=.GlobalEnv)},add=TRUE)
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(1202)
  d <- data.frame(group=factor(rep(c("a","b"),c(55,105))),x=rnorm(160),z=rnorm(160))
  s <- 2*d$x + d$z + .8*d$x*(d$group=="b")
  items <- paste0("i",1:4);d[items] <- replicate(4,s+rnorm(160))
  coding <- list(group="contr.treatment")
  design <- lambdastar::lambda_design(d,~group*x+z,contrasts=coding)
  effect <- lambdastar::lambda_hypothesis(d,~group*x+z,"group:x",c("group","x","z"),contrasts=coding)
  cfg <- lambdastar::lambda_bootstrap_control(B=B,min_valid=min(200L,B),
    keep_draws=TRUE,keep_indices=TRUE)
  full <- lambdastar::lambda_model(d,design,items,1/12)
  conditional <- lambdastar::lambda_model(d,effect,items,1/12,bootstrap=cfg)
  # Equal stable group means with unequal stable dispersions, plus one common bath.
  dispersion <- unlist(lapply(c(55,105),function(n) {
    u <- qexp((seq_len(n)-.5)/n);u-mean(u)
  })) * ifelse(d$group=="a",.6,1.8)
  variance_data <- d
  variance_data[items] <- replicate(4,dispersion+rnorm(160))
  variance_fit <- lambdastar::lambda_model(variance_data,~group,items,1/12,contrasts=coding)
  oracle <- c(T=.25,lambda_star=sum(dispersion^2)/.25,lambda_H=0,kappa=0)
  welch <- stats::t.test(rowMeans(variance_data[items])~group,variance_data)
  list(data=d,items=items,signal=s,design=design,hypothesis=effect,full=full,conditional=conditional,
    unequal_variance=list(data=variance_data,signal=dispersion,oracle=oracle,
      fit=variance_fit,welch=welch))
}
