# Source, then call lambdastar_intervals_example(). Sourcing runs no analysis.
# This lesson uses two retained synthetic E4 datasets, with fresh seed=1 draws.
lambdastar_intervals_example <- function(B=5000L, conf.level=.95, seed=1L) {
  cfg <- lambdastar::lambda_bootstrap_control(B=B,conf.level=conf.level,
    seed=seed,min_valid=min(B,200L),keep_draws=TRUE,keep_indices=TRUE)
  directory <- system.file("examples","data",package="lambdastar")
  provenance <- read.csv(file.path(directory,"e4-provenance.csv"))
  items <- paste0("i",1:4)
  fits <- data <- designs <- list()
  for(name in c("baseline","weak")) {
    d <- read.csv(file.path(directory,paste0("e4-",name,"-1001.csv")))
    d$group <- factor(d$group,levels=c("a","b"))
    data[[name]] <- d
    designs[[name]] <- lambdastar::lambda_design(d,~group,contrasts=list(group="contr.treatment"))
    fits[[name]] <- lambdastar::lambda_model(d,designs[[name]],items,T_q=1/12,
      estimator="paper",bootstrap=cfg)
  }
  # Zero total estimates stay in the lambda-star distribution. Kappa is undefined
  # for those same draws; the package records why its interval is unavailable.
  intervals <- do.call(rbind,lapply(names(fits),function(name) {
    f <- fits[[name]]
    truth <- subset(provenance,scenario==name)
    cbind(scenario=name,population_target=truth$truth[match(f$bootstrap$ci$parameter,truth$parameter)],
      f$bootstrap$ci)
  }))
  rownames(intervals) <- NULL
  list(data=data,items=items,designs=designs,fits=fits,intervals=intervals,
    provenance=provenance,settings=list(B=B,conf.level=conf.level,seed=seed,T_q=1/12))
}
