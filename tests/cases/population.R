pop_data <- read.csv(fixture("ideal_costly_model3.csv"))
pop_items <- paste0("y_item_",1:6)
pop_cfg <- lambda_bootstrap_control(B=99,min_valid=50,seed=37,keep_draws=TRUE,keep_indices=TRUE)
pop_fit <- lambda_model(estimator="population", pop_data,~x1+x2+x3,pop_items,1/12,bootstrap=pop_cfg)

test("population estimator is independent of bootstrap and confidence level", {
  plain <- lambda_model(estimator="population", pop_data,~x1+x2+x3,pop_items,1/12)
  equal(plain$population,pop_fit$population)
  equal(unname(pop_fit$population$estimates),pop_fit$bootstrap$ci$estimate)
  cfg <- pop_cfg; cfg$conf.level <- .8; cfg$B <- 19L; cfg$min_valid <- 10L
  other <- lambda_model(estimator="population", pop_data,~x1+x2+x3,pop_items,1/12,bootstrap=cfg)
  equal(other$population,plain$population)
  equal(pop_cfg$ci.method,"percentile")
})

test("every stored draw matches independent pairwise moments and OLS", {
  for(b in seq_len(pop_cfg$B)) {
    d <- pop_data[pop_fit$bootstrap$indices[b,],]
    y <- rowMeans(d[pop_items]); n <- length(y)
    fit <- lm.fit(cbind(1,as.matrix(d[c("x1","x2","x3")])),y)
    pairs <- combn(6,2)
    T <- mean(apply(pairs,2,function(j) var(d[[pop_items[j[1]]]]-d[[pop_items[j[2]]]])))/12
    total <- max(0,var(y)-T)
    h <- min(total,max(0,var(y)-sum(fit$residuals^2)/(n-fit$rank)))
    expected <- c(T,(n-1)*total/T,(n-1)*h/T,(n-1)*(total-h)/T,h/total)
    equal(unname(pop_fit$bootstrap$draws[b,]),unname(expected),tol=1e-9)
  }
})

test("percentile construction agrees with boot on the identical person draws", {
  # Independent resampling engine; evaluate its rows using the public estimator.
  set.seed(637)
  ref <- boot::boot(pop_data,function(data,indices) {
    f <- lambda_model(estimator="population", data[indices,],~x1+x2+x3,pop_items,1/12,capacity="smooth")
    unname(f$population$estimates)
  },R=399)
  cfg <- lambda_bootstrap_control(B=399,keep_draws=TRUE)
  table <- getFromNamespace(".bootstrap_ci_table","lambdastar")(
    setNames(ref$t0,names(pop_fit$population$estimates)),ref$t,cfg)$ci
  for(j in 1:5) {
    expected <- boot::boot.ci(ref,conf=.95,type="perc",index=j)$percent[1,4:5]
    equal(c(table$lower[j],table$upper[j]),unname(expected))
  }
})

test("weak signal records undefined ratios without substituting an interval", {
  set.seed(2026020001)
  d <- data.frame(x1=rnorm(100),x2=rnorm(100))
  s <- .1*(d$x1+d$x2)+sqrt(.02)*rnorm(100)
  d[paste0("i",1:4)] <- replicate(4,s+rnorm(100))
  cfg <- lambda_bootstrap_control(B=299,min_valid=50,seed=120001,keep_draws=TRUE,
    conditional_quantiles=TRUE)
  fit <- lambda_model(estimator="population", d,~x1+x2,paste0("i",1:4),1/12,bootstrap=cfg)
  ci <- fit$bootstrap$ci; k <- ci[ci$parameter=="kappa",]
  z <- fit$bootstrap$draws[,"lambda_star"]==0
  stopifnot(any(z),any(!z),all(is.na(fit$bootstrap$draws[z,"kappa"])),
    k$n_undefined==sum(z),k$n_failed==sum(z),is.na(k$lower),is.na(k$upper),
    all(is.finite(ci$lower[1:4])))
  stopifnot(k$status %in% c("undefined_estimate","undefined_draws"))
  cond <- fit$bootstrap$conditional_quantiles
  stopifnot(nrow(cond)==1,cond$parameter=="kappa",cond$n_excluded==sum(z))
  fit2 <- fit; fit2$bootstrap$config$conditional_quantiles <- FALSE
  equal(confint(fit2,level=.8),confint(fit,level=.8))
})

test("valid zero/one ratios and degenerate distributions retain their percentiles", {
  table <- getFromNamespace(".bootstrap_ci_table","lambdastar")
  cfg <- lambda_bootstrap_control(B=299,min_valid=100)
  draws <- cbind(lambda_star=rep(10,299),lambda_h=rep(10,299),kappa=rep(1,299))
  x <- table(c(lambda_star=10,lambda_h=10,kappa=1),draws,cfg)$ci
  equal(x$lower,c(10,10,1)); equal(x$upper,c(10,10,1))
  stopifnot(all(x$status=="ok"))
  d <- pop_fit$bootstrap$draws
  d[,"kappa"] <- 0
  x <- table(pop_fit$population$estimates,d,pop_cfg)$ci
  equal(x$lower[5],0); equal(x$upper[5],0)
  # A zero-signal point estimate remains undefined even with finite resamples.
  estimates <- pop_fit$population$estimates; estimates["kappa"] <- NA_real_
  x <- table(estimates,d,pop_cfg)$ci
  stopifnot(x$status[5]=="undefined_estimate",is.na(x$lower[5]))
})

test("fixed-stratum population draws match independent mixture variance", {
  d <- pop_data[1:60,]; d$g <- factor(rep(c("a","b"),c(20,40)))
  cfg <- pop_cfg; cfg$strata <- "g"
  fit <- lambda_model(estimator="population", d,~g,pop_items,1/12,bootstrap=cfg)
  equal(fit$rho_star,fit$lambda_star/(fit$lambda_star+fit$d))
  for(b in c(1,99)) {
    selected <- d[fit$bootstrap$indices[b,],]; y <- rowMeans(selected[pop_items])
    z <- split(y,selected$g); w <- c(1/3,2/3); n <- c(20,40)
    means <- vapply(z,mean,numeric(1)); v <- vapply(z,var,numeric(1))
    vy <- sum(w*v)+sum(w*(means-sum(w*means))^2)-sum(w*(1-w)*v/n)
    f <- lm(y~selected$g)
    T <- fit$bootstrap$draws[b,1]; S <- max(0,vy-T)
    H <- min(S,max(0,vy-sum(residuals(f)^2)/58))
    equal(unname(fit$bootstrap$draws[b,"lambda_h"]),unname(59*H/T))
    equal(selected$g,d$g)
  }
})

test("null models have valid zero localization for positive signal", {
  f <- lambda_model(estimator="population", pop_data,~1,pop_items,1/12,bootstrap=pop_cfg)
  equal(f$population$estimates["lambda_h"],c(lambda_h=0))
  stopifnot(all(f$bootstrap$draws[,"lambda_h"]==0),
    all(f$bootstrap$draws[,"kappa"]==0),all(confint(f,parm="kappa")==0))
})

test("saturated models keep finite-person F_H and identify unestimable population projection", {
  d <- data.frame(i1=c(1,2,3,4),i2=c(2,1,5,3),g=factor(1:4))
  fit <- lambda_model(estimator="population", d,~g,c("i1","i2"),1/12,bootstrap=pop_cfg)
  stopifnot(is.finite(fit$parsimony$free_energy),
    is.na(fit$population$estimates["lambda_h"]),
    grepl("residual degrees",fit$population$model_issues),
    fit$bootstrap$ci$status[3]=="undefined_estimate")
})

test("singular temperature draws preserve temperature zeros and report dependent failures", {
  d <- data.frame(i1=c(0,0,1),i2=c(0,1,0))
  cfg <- lambda_bootstrap_control(B=299,min_valid=50,keep_draws=TRUE)
  f <- lambda_measure(estimator="population", d,names(d),bootstrap=cfg)
  z <- f$bootstrap$draws[,"temperature"]==0
  stopifnot(any(z),all(is.na(f$bootstrap$draws[z,"lambda_star"])),
    f$bootstrap$ci$n_failed[1]==0,f$bootstrap$ci$n_zero[1]==sum(z),
    f$bootstrap$ci$n_failed[2]==sum(z))
})


test("zero total point retains a lambda interval and an undefined kappa point", {
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(160),40,4)); indicators <- names(d)
  d$x <- seq_len(40)
  cfg <- lambda_bootstrap_control(B=299,min_valid=100,keep_draws=TRUE)
  f <- lambda_model(estimator="population", d,~x,indicators,1/12,bootstrap=cfg)
  total <- f$bootstrap$ci[2,]; k <- f$bootstrap$ci[5,]
  stopifnot(total$estimate==0,total$lower==0,total$upper>0,total$status=="ok",
    is.na(k$estimate),is.na(k$lower),k$status=="undefined_estimate",
    all(is.na(f$bootstrap$draws[f$bootstrap$draws[,"lambda_star"]==0,"kappa"])))
})

test("unrepresentable tail probabilities are reported without substituting limits", {
  cfg <- pop_cfg; cfg$conf.level <- 1-.Machine$double.eps/2
  table <- getFromNamespace(".bootstrap_ci_table","lambdastar")(
    pop_fit$population$estimates,pop_fit$bootstrap$draws,cfg)$ci
  stopifnot(all(table$status=="quantile_unavailable"),all(is.na(table$lower)),
    all(is.na(table$upper)),all(nzchar(table$quantile_message)))
})
