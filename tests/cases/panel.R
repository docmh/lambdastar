# E6b: explicit dummy-model references independent of implicit person projection.
panel_data <- function(unbalanced=FALSE) {
  set.seed(6821)
  d <- data.frame(id=rep(1:12,each=4),time=rep(0:3,12))
  d$group <- factor(ifelse(d$id<=6,"a","b"))
  d$z <- rnorm(48);u <- rnorm(12);v <- rnorm(12,sd=.6)
  signal <- 10+u[d$id]+(.8+v[d$id]+.4*(d$group=="b"))*d$time+.3*d$z
  for(k in 1:4)d[[paste0("y",k)]] <- signal+rnorm(48,sd=.5)
  if(unbalanced)d <- d[!(d$id%%2==0&d$time==3),]
  d
}
panel_items <- paste0("y",1:4)
panel_cfg <- function(B=19,...)lambda_bootstrap_control(B=B,min_valid=2,seed=29,
  resampling="cluster",cluster="id",keep_draws=TRUE,keep_indices=TRUE,...)
panel_oracle <- function(d,X,extra=NULL) {
  Y <- as.matrix(d[panel_items]);n <- nrow(Y);k <- ncol(Y)
  B <- cbind(model.matrix(~0+factor(id),d),extra)
  qb <- qr(B);df <- n-qb$rank;Z <- qr.resid(qb,Y);z <- rowMeans(Z)
  T <- sum((Z-z)^2)/(df*(k-1)*k)
  A <- qr.resid(qb,X);p <- ncol(A)
  fitted <- if(p) qr.fitted(qr(A),z) else z*0
  total <- max(0,sum(z^2)/T-df);h <- min(total,max(0,sum(fitted^2)/T-p))
  v <- if(p) sqrt(det(crossprod(sweep(A,2,sqrt(colSums(A^2)),"/")))) else 1
  list(estimates=c(temperature=T,lambda_star=total,lambda_h=h,lambda_perp=total-h,
    kappa=if(total>0)h/total else NA_real_),d=as.integer(df),rank=as.integer(p),
    rss=sum((z-fitted)^2),volume=v,Z=Z)
}

test("panel time trend matches explicit person dummy least squares", {
  d <- panel_data();h <- lambda_panel(d,~time,person="id",occasion="time")
  f <- lambda_model(d,h,panel_items,1/12);o <- panel_oracle(d,matrix(d$time,ncol=1))
  equal(f$estimates,o$estimates);equal(f$d,o$d);equal(f$df_h,o$rank)
  equal(f$parsimony$rss,o$rss);equal(f$panel$n_persons,12L)
  equal(f$nuisance$rank,12L);stopifnot(is.null(f$nuisance$qr))
  equal(f$panel$map$input_row,seq_len(nrow(d)));equal(f$inference$target,"panel_projection")
  shuffled <- d[sample.int(nrow(d)),]
  other <- lambda_model(shuffled,lambda_panel(shuffled,~time,person="id",occasion="time"),panel_items,1/12)
  equal(other$estimates,f$estimates);equal(other$parsimony$free_energy,f$parsimony$free_energy)
})

test("default panel occasion effect is categorical with literal coding", {
  d <- panel_data();f <- lambda_model(d,lambda_panel(d,person="id",occasion="time"),panel_items,1/12)
  X <- model.matrix(~factor(time),d)[,-1,drop=FALSE];o <- panel_oracle(d,X)
  equal(f$estimates,o$estimates);equal(f$df_h,3L);equal(f$parsimony$volume,o$volume)
  equal(f$parsimony$free_energy,o$rss+3*f$parsimony$xi_eta-f$parsimony$xi_a*log(o$volume))
})

test("additional nuisance is projected within persons with correct rank", {
  d <- panel_data();h <- lambda_panel(d,~time,person="id",occasion="time",nuisance=~z+group)
  f <- lambda_model(d,h,panel_items,1/12);o <- panel_oracle(d,matrix(d$time,ncol=1),matrix(d$z,ncol=1))
  equal(f$estimates,o$estimates);equal(f$d,35L);equal(f$nuisance$rank,13L)
  equal(f$parsimony$rss,o$rss)
})

test("panel target blocks express group by occasion effects after person adjustment", {
  d <- panel_data();h <- lambda_panel(d,~group*factor(time),person="id",occasion="time",
    target="group:factor(time)",adjust="factor(time)")
  f <- lambda_model(d,h,panel_items,1/12)
  mm <- model.matrix(~group*factor(time),d);asgn <- attr(mm,"assign")
  o <- panel_oracle(d,mm[,asgn==3,drop=FALSE],mm[,asgn==2,drop=FALSE])
  equal(f$estimates,o$estimates);equal(f$df_h,3L);equal(f$d,33L)
})

test("person slopes rebuild separate copies in every bootstrap draw", {
  d <- panel_data();h <- lambda_panel(d,~id*time,person="id",occasion="time",
    target="id:time",adjust="time")
  f <- lambda_model(d,h,panel_items,1/12,bootstrap=panel_cfg())
  stopifnot(!any(f$bootstrap$failures$code%in%c("rank_deficient","coding_changed","factor_levels")))
  for(j in c(1,8,19)) {
    ix <- f$bootstrap$indices[[j]];s <- d[ix,];s$id <- f$bootstrap$sampling$person_instances[[j]]
    mm <- model.matrix(~factor(id)*time,s);asgn <- attr(mm,"assign")
    o <- panel_oracle(s,mm[,asgn==3,drop=FALSE],mm[,asgn==2,drop=FALSE])
    equal(unname(f$bootstrap$draws[j,]),unname(o$estimates),tol=1e-8)
    equal(f$bootstrap$sampling$draw_df[j],o$d)
  }
})

test("unbalanced panel bootstrap retains each selected person's entire pattern", {
  d <- panel_data(TRUE);f <- lambda_model(d,lambda_panel(d,~time,person="id",occasion="time"),
    panel_items,1/12,bootstrap=panel_cfg())
  stopifnot(length(unique(f$bootstrap$sampling$draw_sizes))>1)
  for(j in seq_len(19)) {
    ix <- f$bootstrap$indices[[j]];instance <- f$bootstrap$sampling$person_instances[[j]]
    s <- d[ix,];s$id <- instance
    for(g in unique(instance)) {
      input <- ix[instance==g];id <- d$id[input[1]]
      equal(input,which(d$id==id))
    }
    o <- panel_oracle(s,matrix(s$time,ncol=1))
    equal(unname(f$bootstrap$draws[j,]),unname(o$estimates),tol=1e-8)
    equal(f$bootstrap$sampling$draw_df[j],as.integer(length(ix)-12L))
  }
})

test("person resampling keeps stratum counts and rejects time-varying strata", {
  d <- panel_data();h <- lambda_panel(d,~time,person="id",occasion="time")
  f <- lambda_model(d,h,panel_items,1/12,bootstrap=panel_cfg(strata="group"))
  for(ix in f$bootstrap$indices)equal(as.integer(table(d$group[ix])),c(24L,24L))
  fails(lambda_model(d,h,panel_items,1/12,bootstrap=panel_cfg(strata="time")),"constant within")
  fails(lambda_model(d,h,panel_items,1/12,bootstrap=lambda_bootstrap_control(B=9)),"whole persons")
  cfg <- getFromNamespace(".panel_bootstrap","lambdastar")(TRUE,list(h))
  equal(cfg$B,5000L);equal(cfg$cluster,"id");equal(cfg$ci.method,"percentile")
})

test("panel complete selection removes whole persons on common comparison variables", {
  d <- panel_data();d$y1[3] <- NA;d$z[10] <- NA
  a <- lambda_panel(d,~time,person="id",occasion="time")
  b <- lambda_panel(d,~time+z,person="id",occasion="time")
  fails(lambda_compare(d,a,b,panel_items,1/12),"Missing values")
  f <- lambda_compare(d,a,b,panel_items,1/12,na.action="complete",bootstrap=panel_cfg(9))
  equal(f$cases$excluded,c(1:4,9:12));equal(f$panel$n_persons,10L)
  equal(f$model_a$cases,f$model_b$cases)
  for(ix in f$bootstrap$indices) stopifnot(!any(ix%in%c(1:4,9:12)))
})

test("panel comparisons reuse the same measurement and standalone numerical F", {
  d <- panel_data();a <- lambda_panel(d,~time,person="id",occasion="time")
  b <- lambda_panel(d,~poly(time,2,raw=TRUE),person="id",occasion="time")
  f <- lambda_compare(d,a,b,panel_items,1/12)
  fa <- lambda_model(d,a,panel_items,1/12);fb <- lambda_model(d,b,panel_items,1/12)
  equal(f$delta_F,fb$parsimony$free_energy-fa$parsimony$free_energy)
  fails(lambda_compare(d,a,~time,panel_items,1/12),"same person/occasion")
  c <- lambda_panel(d,~time,person="id",occasion="time",person_effects="none")
  fails(lambda_compare(d,a,c,panel_items,1/12),"same person/occasion")
})

test("person fixed effects can be targets with distinct bootstrap identities", {
  d <- panel_data();h <- lambda_panel(d,~id+time,person="id",occasion="time",person_effects="none")
  f <- lambda_model(d,h,panel_items,1/12,bootstrap=panel_cfg(9))
  stopifnot(!nrow(f$bootstrap$failures));equal(f$df_h,12L)
  j <- 1;ix <- f$bootstrap$indices[[j]];s <- d[ix,];s$id <- factor(f$bootstrap$sampling$person_instances[[j]])
  g <- lambda_model(s,~id+time,panel_items,1/12)
  equal(unname(f$bootstrap$draws[j,]),unname(g$estimates),tol=1e-8)
})

test("paired two-occasion reduction agrees with E6a lambda values and scaled F", {
  d <- panel_data();d <- d[d$time%in%0:1,]
  f <- lambda_model(d,lambda_panel(d,~time,person="id",occasion="time"),panel_items,1/12)
  wide <- data.frame(d[d$time==0,panel_items],d[d$time==1,panel_items]);names(wide) <- c(paste0("a",1:4),paste0("b",1:4))
  g <- lambda_paired(wide,paste0("a",1:4),paste0("b",1:4),T_q=2/12)
  equal(f$estimates[-1],g$estimates[-1],tol=1e-8)
  equal(f$temperature*2,g$temperature)
  equal(f$parsimony$free_energy*2,g$parsimony$free_energy,tol=1e-8)
})

test("panel scale and person-constant shifts preserve appropriate quantities", {
  d <- panel_data();h <- lambda_panel(d,~time,person="id",occasion="time")
  f <- lambda_model(d,h,panel_items,1/12)
  s <- d;s[panel_items] <- s[panel_items]+s$id*100
  g <- lambda_model(s,lambda_panel(s,~time,person="id",occasion="time"),panel_items,1/12)
  equal(f$estimates,g$estimates,tol=1e-8)
  s[panel_items] <- 3*d[panel_items]
  g <- lambda_model(s,lambda_panel(s,~time,person="id",occasion="time"),panel_items,9/12)
  equal(f$estimates[-1],g$estimates[-1],tol=1e-8);equal(f$temperature*9,g$temperature)
  equal(f$parsimony$free_energy*9,g$parsimony$free_energy,tol=1e-8)
  item <- lambda_model(d,h,panel_items,1/12,measurement_scale="item")
  equal(item$lambda_star*4,f$lambda_star);equal(item$kappa,f$kappa)
})

test("panel validation rejects ambiguous maps and incompatible estimator paths", {
  d <- panel_data();h <- lambda_panel(d,~time,person="id",occasion="time")
  fails(lambda_panel(d[c(1,seq_len(nrow(d))),],~time,person="id",occasion="time"),"pair must be unique")
  s <- d;s$id[1] <- NA;fails(lambda_panel(s,~time,person="id",occasion="time"),"nonmissing")
  fails(lambda_model(d,h,panel_items,1/12,estimator="population"),"no panel population")
  fails(lambda_model(d[nrow(d):1,],h,panel_items,1/12),"row identities")
  fails(lambda_panel(d,~id*time,person="id",occasion="time",encoding="fixed"),"require encoding")
  fails(lambda_contrast(d,h,c(time=1)),"Panel term contrasts")
  fails(lambda_model(d,lambda_panel(d,~group,person="id",occasion="time"),panel_items,1/12),"redundant")
})

test("panel fixed non-person bases and explicit matrices remain usable", {
  d <- panel_data();h <- lambda_panel(d,~poly(time,2),person="id",occasion="time",encoding="fixed")
  f <- lambda_model(d,h,panel_items,1/12,bootstrap=panel_cfg(9))
  X <- unname(poly(d$time,2));g <- lambda_model(d,lambda_panel(d,X,person="id",occasion="time"),panel_items,1/12,bootstrap=panel_cfg(9))
  equal(f$estimates,g$estimates);equal(f$bootstrap$draws,g$bootstrap$draws)
  equal(f$parsimony$free_energy,g$parsimony$free_energy)
})

test("declared person factor contrasts survive panel preparation and resampling", {
  d <- panel_data();d$id <- factor(d$id);C <- contr.helmert(levels(d$id))
  contrasts(d$id) <- C
  h <- lambda_panel(d,~id+time,person="id",occasion="time",person_effects="none")
  f <- lambda_model(d,h,panel_items,1/12,bootstrap=panel_cfg(9))
  equal(f$model$contrasts$id,C)
  plain <- d;attr(plain$id,"contrasts") <- NULL
  g <- lambda_model(plain,lambda_panel(plain,~id+time,person="id",occasion="time",
    person_effects="none",contrasts=list(id=C)),panel_items,1/12,bootstrap=panel_cfg(9))
  equal(f$estimates,g$estimates);equal(f$parsimony$free_energy,g$parsimony$free_energy)
  equal(f$bootstrap$draws,g$bootstrap$draws)
})

test("panel zero temperature and zero signal retain explicit diagnostic meanings", {
  d <- panel_data();for(k in panel_items)d[[k]] <- d$id+d$time
  f <- lambda_model(d,lambda_panel(d,~time,person="id",occasion="time"),panel_items,1/12)
  equal(f$status,"singular_temperature");stopifnot(is.na(f$parsimony$free_energy))
  d <- panel_data();v <- rep(c(-1,1,-1,1),12)
  for(k in 1:4)d[[panel_items[k]]] <- d$id+c(-1,-1,1,1)[k]*v
  f <- lambda_model(d,lambda_panel(d,~time,person="id",occasion="time"),panel_items,1/12)
  equal(f$lambda_star,0);stopifnot(is.na(f$kappa),f$temperature>0)
})

test("saturated person by occasion interactions retain positive bath df", {
  d <- panel_data();h <- lambda_panel(d,~id*factor(time),person="id",occasion="time",
    target="id:factor(time)",adjust="factor(time)")
  f <- lambda_model(d,h,panel_items,1/12,bootstrap=panel_cfg(5))
  equal(f$d,33L);equal(f$df_h,f$d);stopifnot(f$temperature>0)
  equal(f$lambda_perp,0);stopifnot(!any(f$bootstrap$failures$code=="rank_deficient"))
  equal(unname(f$bootstrap$draws[,"lambda_perp"]),rep(0,5))
  equal(unname(f$bootstrap$draws[,"kappa"]),rep(1,5))
})

test("unavailable library percentiles return diagnostics without substituted bounds", {
  limits <- getFromNamespace(".percentile_limits","lambdastar")
  capture.output(x <- limits(c(1,1+.Machine$double.eps),.95,1))
  stopifnot(length(x$limits)==2L,all(is.na(x$limits)),grepl("did not return",x$message))
})
