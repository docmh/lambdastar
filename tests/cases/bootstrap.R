boot_data <- read.csv(fixture("ideal_costly_model3.csv"))
boot_items <- paste0("y_item_", 1:6)
boot_cfg <- lambda_bootstrap_control(B = 37, conf.level = .8, seed = 27,
                                    min_valid = 20, keep_draws = TRUE, keep_indices = TRUE)
boot_fit <- lambda_compare(boot_data, ~x1, ~x1+x2+x3, boot_items, 1/12, bootstrap = boot_cfg)

test("bootstrap defaults and explicit overrides are resolved and validated", {
  d <- lambda_bootstrap_control()
  equal(d$B, 5000L)
  equal(d$conf.level, .95)
  equal(d$seed, 1L)
  equal(d$ci.method, "percentile")
  equal(d$ci.sides, "two.sided")
  equal(d$resampling, "case")
  equal(boot_fit$bootstrap$config, boot_cfg)
  equal(boot_fit$settings$bootstrap, boot_cfg)
  equal(boot_fit$bootstrap$inference$target, "paper_projection")
  equal(boot_fit$bootstrap$inference$scale_factor, nrow(boot_data)-1L)
  for (bad in list(list(B=1), list(B=2.5), list(seed=NA), list(seed=-1),
    list(conf.level=1), list(conf.level=0), list(ci.sides="one.sided"),
    list(ci.method="bca"), list(resampling="residual"), list(strata=c("x1","x2")),
    list(keep_draws=1), list(conditional_quantiles=1), list(min_valid=1))) {
    fails(do.call(lambda_bootstrap_control, bad), ".")
  }
  fails(lambda_measure(boot_data, boot_items, bootstrap = list(B=20)), "bootstrap must")
})

test("same seed gives identical joint draws and changing seed changes draws", {
  a <- lambda_compare(boot_data, ~x1, ~x1+x2+x3, boot_items, 1/12, bootstrap = boot_cfg)
  equal(a$bootstrap$draws, boot_fit$bootstrap$draws)
  equal(a$bootstrap$indices, boot_fit$bootstrap$indices)
  equal(a$bootstrap$ci, boot_fit$bootstrap$ci)
  cfg <- boot_cfg
  cfg$seed <- 28L
  b <- lambda_compare(boot_data, ~x1, ~x1+x2+x3, boot_items, 1/12, bootstrap = cfg)
  stopifnot(!identical(a$bootstrap$indices, b$bootstrap$indices))
  equal(a$parsimony, NULL) # F_H lives only in the two model results.
  equal(a$model_a$parsimony,
    lambda_model(boot_data, ~x1, boot_items, 1/12)$parsimony)
})

test("bootstrap preserves RNG kind and present or absent caller seed", {
  local({
    old_kind <- RNGkind()
    existed <- exists(".Random.seed", .GlobalEnv, inherits=FALSE)
    if (existed) old_seed <- get(".Random.seed", .GlobalEnv)
    on.exit({ do.call(RNGkind, as.list(old_kind))
      if (existed) assign(".Random.seed", old_seed, .GlobalEnv)
      else if (exists(".Random.seed", .GlobalEnv, inherits=FALSE)) rm(".Random.seed", envir=.GlobalEnv)
    })
    RNGkind("L'Ecuyer-CMRG")
    set.seed(700)
    before <- get(".Random.seed", .GlobalEnv)
    lambda_measure(boot_data, boot_items, bootstrap=lambda_bootstrap_control(B=3,min_valid=2))
    equal(RNGkind()[1], "L'Ecuyer-CMRG")
    equal(get(".Random.seed", .GlobalEnv), before)
    rm(".Random.seed", envir=.GlobalEnv)
    lambda_measure(boot_data, boot_items, bootstrap=lambda_bootstrap_control(B=3,min_valid=2))
    stopifnot(!exists(".Random.seed", .GlobalEnv, inherits=FALSE))
  })
})

test("joint resampling matches independent pairwise moments and lm.fit", {
  draws <- boot_fit$bootstrap$draws
  for (b in c(1, 7, 37)) {
    selected <- boot_fit$bootstrap$indices[b, ]
    d <- boot_data[selected, ]
    y <- rowMeans(d[boot_items])
    k <- length(boot_items)
    pairs <- combn(k,2)
    T_ref <- mean(apply(pairs,2,function(j)
      var(d[[boot_items[j[1]]]]-d[[boot_items[j[2]]]])))/(2*k)
    total <- max(0, (length(y)-1)*(var(y)-T_ref)/T_ref)
    equal(unname(draws[b,"temperature"]), unname(T_ref), tol=1e-9)
    equal(unname(draws[b,"lambda_star"]), total, tol=1e-9)
    for (j in 1:2) {
      cols <- if(j==1) "x1" else c("x1","x2","x3")
      fit <- lm.fit(cbind(1,as.matrix(d[cols])), y)
      captured <- sum((fit$fitted.values-mean(y))^2)
      h <- min(total,max(0,captured/T_ref-(fit$rank-1)))
      prefix <- c("model_a.","model_b.")[[j]]
      equal(unname(draws[b,paste0(prefix,"lambda_h")]),h,tol=1e-9)
      equal(unname(draws[b,paste0(prefix,"kappa")]),h/total,tol=1e-9)
    }
  }
  stopifnot(stats::sd(draws[,"temperature"]) > 0)
})

test("percentile intervals use requested level and joint kappa ratios", {
  ci <- boot_fit$bootstrap$ci
  for (j in seq_len(nrow(ci))) {
    reference <- structure(list(t0=ci$estimate[j],t=matrix(boot_fit$bootstrap$draws[,j]),
      R=boot_cfg$B,sim="ordinary"),class="boot")
    expected <- unname(boot::boot.ci(reference,conf=.8,type="perc")$percent[1,4:5])
    equal(c(ci$lower[j],ci$upper[j]),expected)
  }
  bounds <- confint(boot_fit, parm=c("temperature","model_b.kappa"), level=.9)
  reference <- structure(list(t0=ci$estimate[1],t=matrix(boot_fit$bootstrap$draws[,1]),
    R=boot_cfg$B,sim="ordinary"),class="boot")
  equal(as.numeric(bounds[1,]), unname(boot::boot.ci(reference,conf=.9,type="perc")$percent[1,4:5]))
  equal(attr(bounds,"diagnostics")$conf.level, rep(.9,2))
  equal(confint(boot_fit, parm=1:2), confint(boot_fit, parm=c("temperature","lambda_star")))
  fails(confint(boot_fit,parm="F_H"), "Unknown")
  fails(confint(boot_fit,parm=NA), "parm")
  fails(confint(boot_fit,level=1), "less than")
  fails(confint(boot_fit,method="bca"), "Unsupported")
  fails(confint(lambda_measure(boot_data,boot_items)), "No bootstrap")
})

test("missingness and strata use one original selection and retain row mapping", {
  d <- boot_data
  d$group <- rep(c("a","b","c"),length.out=nrow(d))
  d$group[1] <- NA
  d$y_item_1[3] <- NA
  d$x3[9] <- NA
  cfg <- lambda_bootstrap_control(B=19,min_valid=10,strata="group",keep_indices=TRUE,keep_draws=TRUE)
  fails(lambda_compare(d,~x1,~x3,boot_items,1/12,bootstrap=cfg), "Missing")
  fit <- lambda_compare(d,~x1,~x3,boot_items,1/12,na.action="complete",bootstrap=cfg)
  equal(fit$cases$excluded,c(1L,3L,9L))
  equal(fit$bootstrap$sampling$input_rows,fit$cases$included)
  equal(fit$bootstrap$sampling$row_ids,rownames(d)[fit$cases$included])
  for(b in seq_len(cfg$B)) {
    rows <- fit$bootstrap$indices[b,]
    stopifnot(!any(rows %in% fit$cases$excluded))
    equal(table(d$group[rows]),table(d$group[fit$cases$included]))
    equal(d$group[rows],d$group[fit$cases$included])
  }
  equal(fit$bootstrap$inference$design,"fixed_observed_stratum_counts")
})

test("fixed-mixture variance requires estimable within-stratum variances", {
  d <- boot_data[1:12,]
  d$group <- c(rep("many",11),"single")
  cfg <- lambda_bootstrap_control(B=15,min_valid=10,strata="group",keep_indices=TRUE)
  fails(lambda_measure(d,boot_items,bootstrap=cfg,estimator="population"),"at least two")
  d$group <- seq_len(nrow(d))
  fails(lambda_measure(d,boot_items,bootstrap=cfg),"singletons")
  cfg$strata <- "absent"
  fails(lambda_measure(d,boot_items,bootstrap=cfg),"strata must name")
})

test("comparison shares draws with separate matrix/formula/lm pipelines", {
  cfg <- lambda_bootstrap_control(B=15,min_valid=10,keep_draws=TRUE,keep_indices=TRUE)
  d <- boot_data
  d$group <- factor(rep(letters[1:3],80))
  f <- y_mean~group+x1
  fit <- lm(f,d,contrasts=list(group="contr.sum"))
  a <- lambda_model(d,fit,boot_items,1/12,bootstrap=cfg)
  b <- lambda_model(d,f,boot_items,1/12,contrasts=list(group="contr.sum"),bootstrap=cfg)
  m <- model.matrix(fit)[,-1]
  c <- lambda_model(d,m,boot_items,1/12,bootstrap=cfg)
  cmp <- lambda_compare(d,fit,m,boot_items,1/12,bootstrap=cfg)
  equal(a$bootstrap$draws,b$bootstrap$draws)
  equal(a$bootstrap$draws,c$bootstrap$draws)
  equal(unname(cmp$bootstrap$draws[,3:5]),unname(a$bootstrap$draws[,3:5]))
  equal(unname(cmp$bootstrap$draws[,6:8]),unname(a$bootstrap$draws[,3:5]))
  equal(cmp$bootstrap$indices,a$bootstrap$indices)
  equal(cmp$delta_F,0)
})

test("data-dependent transformations rebuild while factor coding stays fixed", {
  d <- boot_data
  d$group <- rep(letters[1:3],80) # characters are converted to fixed factor levels
  cfg <- lambda_bootstrap_control(B=12,min_valid=10,keep_draws=TRUE,keep_indices=TRUE)
  a <- lambda_model(d,~poly(x1,2)+factor(group),boot_items,1/12,
    contrasts=list("factor(group)"="contr.helmert"),bootstrap=cfg)
  for(b in c(1,12)) {
    rows <- a$bootstrap$indices[b,]
    x <- model.matrix(~poly(x1,2)+factor(group),d[rows,],
      contrasts.arg=list("factor(group)"="contr.helmert"))[,-1]
    ref <- lambda_model(d[rows,],unname(x),boot_items,1/12)
    equal(unname(a$bootstrap$draws[b,"lambda_h"]),unname(ref$estimates["lambda_h"]))
  }
  d$binary <- d$x1>0
  b <- lambda_model(d,~binary,boot_items,1/12,bootstrap=cfg)
  equal(b$bootstrap$ci$n_failed,rep(0L,5))
})

test("lost factor levels fail only the affected model and never get replaced", {
  d <- boot_data[1:24,]
  d$rare <- factor(c(rep("a",23),"b"))
  cfg <- lambda_bootstrap_control(B=80,min_valid=20,keep_draws=TRUE)
  a <- lambda_compare(d,~rare,~x1,boot_items,1/12,bootstrap=cfg)
  stopifnot(any(a$bootstrap$failures$scope=="model_a"),
    any(a$bootstrap$failures$code=="rank_deficient"),a$bootstrap$attempted==80)
  table <- a$bootstrap$ci
  row <- table[table$parameter=="model_a.lambda_h",]
  stopifnot(row$n_failed>0,is.na(row$lower),row$status=="failed_draws",
            row$n_valid+row$n_failed==80)
  stopifnot(all(table$n_failed[table$parameter %in% c("temperature","model_b.lambda_h")]==0))
  cfg$conditional_quantiles <- TRUE
  b <- lambda_compare(d,~rare,~x1,boot_items,1/12,bootstrap=cfg)
  equal(a$bootstrap$draws,b$bootstrap$draws)
  table <- b$bootstrap$ci
  row <- table[table$parameter=="model_a.lambda_h",]
  equal(row$status,"failed_draws")
  stopifnot(is.na(row$lower))
  conditional <- b$bootstrap$conditional_quantiles
  row <- conditional[conditional$parameter=="model_a.lambda_h",]
  equal(row$status,"conditional_on_valid")
  stopifnot(is.finite(row$lower))
})

test("zero temperature skips every CI including temperature without sampling", {
  d <- data.frame(i1=1:8,i2=1:8,x=1:8)
  a <- lambda_model(d,~x,c("i1","i2"),1/12,bootstrap=boot_cfg)
  equal(a$bootstrap$attempted,0L)
  equal(a$bootstrap$skipped,boot_cfg$B)
  stopifnot(all(is.na(a$bootstrap$ci$lower)),all(is.na(a$bootstrap$ci$upper)),
    all(a$bootstrap$ci$status=="singular_temperature"))
  equal(nrow(a$bootstrap$draws),0L)
})

test("zero total structure never yields an invented kappa interval", {
  d <- data.frame(i1=c(-1,1,-1,1),i2=c(1,-1,1,-1),x=c(-1,-1,1,1))
  cfg <- lambda_bootstrap_control(B=80,min_valid=10,conditional_quantiles=TRUE,keep_draws=TRUE)
  a <- lambda_model(d,~x,c("i1","i2"),1/12,bootstrap=cfg)
  k <- a$bootstrap$ci[a$bootstrap$ci$parameter=="kappa",]
  equal(a$lambda_star,0)
  equal(k$status,"undefined_estimate")
  stopifnot(is.na(k$lower),is.na(k$upper),k$n_valid==0,
            any(a$bootstrap$failures$code=="undefined_kappa"),
            any(a$bootstrap$failures$code=="singular_temperature"))
})

test("item scaling transforms both point and interval noncentralities", {
  cfg <- lambda_bootstrap_control(B=12,min_valid=10,keep_draws=TRUE)
  a <- lambda_model(boot_data,~x1,boot_items,1/12,bootstrap=cfg)
  b <- lambda_model(boot_data,~x1,boot_items,1/12,measurement_scale="item",bootstrap=cfg)
  equal(a$bootstrap$draws[,1]*6,b$bootstrap$draws[,1])
  equal(a$bootstrap$draws[,2:4]/6,b$bootstrap$draws[,2:4])
  equal(a$bootstrap$draws[,5],b$bootstrap$draws[,5])
  equal(a$parsimony,b$parsimony)
})

test("retention is optional, thresholds explicit, summaries include diagnostics", {
  cfg <- lambda_bootstrap_control(B=5,min_valid=10)
  a <- lambda_measure(boot_data,boot_items,bootstrap=cfg)
  stopifnot(is.null(a$bootstrap$draws),is.null(a$bootstrap$indices),
    all(a$bootstrap$ci$status=="too_few_valid"),all(is.na(a$bootstrap$ci$lower)))
  fails(confint(a,level=.9),"keep_draws")
  equal(summary(a)$bootstrap$ci,a$bootstrap$ci)
  equal(summary(boot_fit)$bootstrap$ci,boot_fit$bootstrap$ci)
  out <- capture.output(print(confint(boot_fit)))
  stopifnot(any(grepl("percentile",out)))
  out <- capture.output(print(boot_fit))
  stopifnot(any(grepl("percentile",out)))
})
