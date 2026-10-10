# Source, then call lambdastar_panel_example(); sourcing does not run analyses.
lambdastar_panel_example <- function(B=5000L,conf.level=.95,seed=1L) {
  had <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  old <- if(had)get(".Random.seed",envir=.GlobalEnv) else NULL;kind <- RNGkind()
  on.exit({do.call(RNGkind,as.list(kind));if(had)assign(".Random.seed",old,envir=.GlobalEnv)
    else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv)},add=TRUE)
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(68261)
  P <- 40L;J <- 4L;K <- 4L
  d <- data.frame(id=rep(seq_len(P),each=J),time=rep(0:3,P))
  d$arm <- factor(ifelse(d$id<=P/2,"control","training"),levels=c("control","training"))
  intercept <- rnorm(P,sd=2);slope <- rnorm(P,sd=.3);curve <- rnorm(P,sd=.2)
  S <- 10+intercept[d$id]+(.5+.65*(d$arm=="training")+slope[d$id])*d$time+
    curve[d$id]*(d$time-1.5)^2
  indicators <- paste0("item",seq_len(K))
  for(k in seq_len(K))d[[indicators[k]]] <- S+rnorm(nrow(d),sd=.7)
  # Preserve treatment-group person counts while resampling whole trajectories.
  cfg <- lambdastar::lambda_bootstrap_control(B=B,conf.level=conf.level,seed=seed,
    resampling="cluster",cluster="id",strata="arm",min_valid=min(B,200L),
    keep_draws=TRUE,keep_indices=TRUE)
  trend_design <- lambdastar::lambda_panel(d,~time,person="id",occasion="time")
  interaction_design <- lambdastar::lambda_panel(d,~arm*time,person="id",occasion="time",
    target="arm:time",adjust="time")
  trend <- lambdastar::lambda_model(d,trend_design,indicators,1/12,bootstrap=cfg)
  interaction <- lambdastar::lambda_model(d,interaction_design,indicators,1/12,bootstrap=cfg)
  curved_design <- lambdastar::lambda_panel(d,~time+I(time^2),person="id",occasion="time")
  comparison <- lambdastar::lambda_compare(d,trend_design,curved_design,indicators,1/12)
  # An explicit model of person-specific slopes, evaluated without bootstrap here.
  slopes_design <- lambdastar::lambda_panel(d,~id*time,person="id",occasion="time",
    target="id:time",adjust="time")
  slopes <- lambdastar::lambda_model(d,slopes_design,indicators,1/12)
  list(data=d,indicators=indicators,stable_profile=S,
    trend_design=trend_design,interaction_design=interaction_design,slopes_design=slopes_design,
    trend=trend,interaction=interaction,slopes=slopes,comparison=comparison,
    settings=list(B=B,conf.level=conf.level,seed=seed,T_q=1/12,
      temperature_truth=.7^2/K,persons=P,occasions=J,parallel_items=K,
      generator_seed=68261,resampling="persons within treatment arms"))
}
