# Load learn_lambdastar.R first. This file defines examples; it runs no analysis
# until lambdastar_advanced_examples() is called. See ADVANCED_EXAMPLES.md.
lambdastar_advanced_examples <- function(B=5000L, seed=29L) {
  if (!exists("lambdastar_tutorial_data",mode="function"))
    stop("Source the installed examples/learn_lambdastar.R first.",call.=FALSE)
  .ls_integer(B,"B",20);.ls_integer(seed,"seed",0)
  config <- lambdastar::lambda_bootstrap_control(B=B,seed=seed+1L,
    min_valid=min(B,200L),keep_draws=TRUE,keep_indices=TRUE,
    percentile_interpolation="type7")
  # Independent random-design population: these profiles are not orthogonalized.
  population_data <- .ls_seed(seed, {
    d <- data.frame(workload=stats::rnorm(160),resources=stats::rnorm(160))
    signal <- -.9*d$workload+.7*d$resources+.5*stats::rnorm(160)
    d[paste0("item",1:4)] <- replicate(4,signal+1.2*stats::rnorm(160))
    d
  })
  indicators <- paste0("item",1:4)
  paper <- lambdastar::lambda_model(population_data,~workload+resources,indicators,1/12,
    estimator="paper",bootstrap=config)
  population <- lambdastar::lambda_model(population_data,~workload+resources,indicators,1/12,
    estimator="population",bootstrap=config)
  cat("\nA. Choose an estimator and its target\n")
  print(rbind(paper=paper$estimates,population=population$estimates))
  print(stats::confint(population))
  population_truth <- c(temperature=1.2^2/4,lambda_star=159*1.55/(1.2^2/4),
    lambda_h=159*1.30/(1.2^2/4),lambda_perp=159*.25/(1.2^2/4),kappa=1.30/1.55)
  cat("Known random-population target at N=160:\n");print(population_truth)
  cat("The population target uses population variances; the realized finite-sample signal is a different object.\n",
      "Both methods use their selected statistic on every draw. They do not become interchangeable because both use percentiles.\n")

  # 24 schools in 4 sites, with unequal school sizes: whole-school draws can
  # contain different numbers of persons. Stable school effects are shared
  # across items and are therefore not treated as item measurement errors.
  sizes <- rep(c(6L,8L,10L),8)
  school <- rep(seq_len(24),sizes)
  site <- factor(rep(rep(1:4,each=6),sizes))
  e <- lambdastar_tutorial_data(n=length(school),seed=seed)
  d <- e$data;d$school <- school;d$site <- site
  school_effect <- .ls_seed(seed+2L,stats::rnorm(24,sd=.6))[school]
  site_effect <- c(-2,-.7,.7,2)[as.integer(site)]
  d[e$indicators] <- d[e$indicators]+site_effect+school_effect
  d$score <- rowMeans(d[e$indicators])
  unadjusted <- lambdastar::lambda_model(d,~workload+resources,e$indicators,1/12)
  adjusted <- lambdastar::lambda_model(d,~workload+resources,e$indicators,1/12,nuisance=~site)
  cat("\nB. Define a common nuisance space\n")
  print(rbind(unadjusted=unadjusted$estimates,within_site=adjusted$estimates))
  cat("N =",nrow(d),"; nuisance rank =",adjusted$nuisance$rank,"; d =",adjusted$d,"\n")
  cat("The adjusted target concerns structure remaining after site effects are removed.\n",
      "Its kappa denominator differs from the unadjusted analysis; these are not two estimates of the same unconditional quantity.\n")

  cluster_config <- lambdastar::lambda_bootstrap_control(B=B,seed=seed+1L,
    min_valid=min(B,200L),keep_draws=TRUE,keep_indices=TRUE,resampling="cluster",
    cluster="school",strata="site",percentile_interpolation="type7")
  clustered <- lambdastar::lambda_model(d,~workload+resources,e$indicators,1/12,
    nuisance=~site,bootstrap=cluster_config)
  cat("\nC. Resample the independent sampling units\n")
  print(stats::confint(clustered))
  cat("Cluster counts by site:\n");print(clustered$bootstrap$sampling$cluster_counts)
  cat("Bootstrap person-count range:",range(clustered$bootstrap$sampling$draw_sizes),"\n")
  cat("Each draw samples six whole schools per site, with replacement. All students and their items travel with the school.\n",
      "Site nuisance effects are refitted and d is recalculated in each draw. Choose this design because of the sampling structure, not to obtain a particular CI width.\n")

  # Deterministic boundary examples distinguish zero signal from zero noise.
  base <- c(-1,0,1,0)
  edge <- data.frame(x=c(0,1,0,-1),i1=base,i2=-base)
  zero_signal <- lambdastar::lambda_model(edge,~x,c("i1","i2"),1/12)
  edge$i2 <- edge$i1
  zero_temperature <- lambdastar::lambda_model(edge,~x,c("i1","i2"),1/12)
  cat("\nD. Distinguish undefined ratios from singular temperature\n")
  print(rbind(zero_signal=zero_signal$estimates,zero_temperature=zero_temperature$estimates))
  cat("Observed T in zero-signal case:",zero_signal$temperature,
      "; observed T in singular case:",zero_temperature$temperature,"\n")
  cat("Opposite item profiles give an observed mean with no stable signal: T is positive, lambda_star=0 and kappa is undefined.\n",
      "Identical measurements give T=0: a singular state, not an infinite regular estimate.\n")
  invisible(list(population_data=population_data,population_truth=population_truth,
    paper=paper,population=population,cluster_data=d,indicators=e$indicators,
    unadjusted=unadjusted,adjusted=adjusted,clustered=clustered,
    zero_signal=zero_signal,zero_temperature=zero_temperature))
}
