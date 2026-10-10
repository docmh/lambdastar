ref_items <- c("i1", "i2", "i3")
set.seed(405)
ref_data <- data.frame(site = factor(rep(c("a", "b"), each = 18)),
  school = rep(rep(letters[1:4], c(3,4,5,6)), 2), x = rnorm(36), z = rnorm(36))
ref_signal <- 3*ref_data$x + as.numeric(ref_data$site) + .6*ref_data$z + rnorm(36,sd=.5)
ref_data[ref_items] <- replicate(3, ref_signal + rnorm(36,sd=.6))

# Independent QR/paired-difference reference; no package estimation functions.
ref_projection <- function(data, x, b = matrix(1,nrow(data),1)) {
  y <- as.matrix(data[ref_items]); qb <- qr(b); df <- nrow(y)-qb$rank
  yr <- qr.resid(qb,y); mean_y <- rowMeans(yr)
  pairs <- combn(ncol(y),2)
  T <- mean(apply(pairs,2,function(j) sum((yr[,j[1]]-yr[,j[2]])^2)/df))/(2*ncol(y))
  xr <- qr.resid(qb,x); f <- lm.fit(xr,mean_y)
  total <- max(0,sum(mean_y^2)/T-df)
  raw_h <- sum(f$fitted.values^2)/T-f$rank
  h <- min(total,max(0,raw_h))
  c(temperature=T,lambda_star=total,lambda_h=h,lambda_perp=total-h,kappa=h/total)
}

test("paper is the default and selected points are the CI statistics", {
  cfg <- lambda_bootstrap_control(B=29,min_valid=10,keep_draws=TRUE,keep_indices=TRUE)
  p <- lambda_model(ref_data,~x,ref_items,.1,bootstrap=cfg)
  equal(p$inference$estimator,"paper"); stopifnot(is.null(p$population))
  equal(unname(p$estimates),p$bootstrap$ci$estimate)
  equal(p$lambda_h,unname(p$estimates["lambda_h"]))
  equal(p$estimates,ref_projection(ref_data,as.matrix(ref_data["x"])))
  for(b in seq_len(cfg$B)) {
    d <- ref_data[p$bootstrap$indices[b,],]
    equal(p$bootstrap$draws[b,],ref_projection(d,as.matrix(d["x"])))
  }
  a <- lambda_model(ref_data,~x,ref_items,.1,estimator="population",bootstrap=cfg)
  equal(a$lambda_h,unname(a$population$estimates["lambda_h"]))
  equal(a$estimates,a$population$estimates)
  equal(a$paper$estimates,p$estimates)
  equal(a$parsimony,p$parsimony)
  stopifnot(abs(a$lambda_h-p$lambda_h)>.01)
  equal(attr(confint(a),"inference")$estimator,"population")
})

test("nuisance formula and full matrix reproduce independent conditional moments", {
  B <- model.matrix(~site+z,ref_data); X <- as.matrix(ref_data["x"])
  a <- lambda_model(ref_data,~x,ref_items,.1,nuisance=~site+z)
  b <- lambda_model(ref_data,X,ref_items,.1,nuisance=B)
  c <- estimate_lambda_partition_items(as.matrix(ref_data[ref_items]),X,nuisance=B)
  expected <- ref_projection(ref_data,X,B)
  equal(a$estimates,expected); equal(b$estimates,expected)
  equal(a$lambda_h,c$lambda_h); equal(a$d,33L); equal(a$nuisance$rank,3L)
  equal(a$df_perp,32L); equal(a$parsimony$df_perp,32L)
  fit <- lm.fit(cbind(B,X),rowMeans(ref_data[ref_items]))
  equal(a$parsimony$rss,sum(fit$residuals^2))
  low <- thermodynamic_model_state(rowMeans(ref_data[ref_items]),X,a$temperature,.1,nuisance=B)
  equal(a$parsimony$free_energy,low$free_energy)
  shifted <- ref_data
  shifted[ref_items] <- shifted[ref_items]+7*as.numeric(shifted$site)+3*shifted$z
  equal(lambda_model(shifted,~x,ref_items,.1,nuisance=~site+z)$estimates,a$estimates)
  equal(lambda_measure(ref_data,ref_items,nuisance=~site+z)$lambda_star,a$lambda_star)
  fails(lambda_model(ref_data,~site,ref_items,.1,nuisance=~site),"redundant")
  fails(lambda_model(ref_data,~x,ref_items,.1,nuisance=~site,estimator="population"),"intercept adjustment only")
})

test("nuisance and cluster columns participate in the common case selection", {
  d <- ref_data; d$z[2] <- NA; d$school[4] <- NA; d$x[7] <- NA
  cfg <- lambda_bootstrap_control(B=5,min_valid=2,resampling="cluster",cluster="school",strata="site")
  fails(lambda_compare(d,~x,~1,ref_items,.1,nuisance=~site+z,bootstrap=cfg),"Missing")
  a <- lambda_compare(d,~x,~1,ref_items,.1,nuisance=~site+z,bootstrap=cfg,na.action="complete")
  equal(a$cases$excluded,c(2L,4L,7L))
  equal(a$model_a$d,30L);equal(a$model_b$d,30L)
  equal(a$estimates["model_a.lambda_h"],setNames(a$model_a$lambda_h,"model_a.lambda_h"))
})

test("cluster sampling preserves whole clusters and refits nuisance for every draw", {
  cfg <- lambda_bootstrap_control(B=79,min_valid=20,seed=44,keep_draws=TRUE,
    keep_indices=TRUE,resampling="cluster",cluster="school",strata="site",percentile_interpolation="type7")
  a <- lambda_model(ref_data,~x,ref_items,.1,nuisance=~site+z,bootstrap=cfg)
  stopifnot(is.list(a$bootstrap$indices),length(unique(a$bootstrap$sampling$draw_sizes))>1)
  equal(unname(a$bootstrap$sampling$cluster_counts),c(4L,4L))
  for(j in seq_len(cfg$B)) {
    ix <- a$bootstrap$indices[[j]]; d <- ref_data[ix,]
    for(g in split(seq_len(nrow(ref_data)),interaction(ref_data$site,ref_data$school)))
      stopifnot(length(unique(tabulate(ix,nbins=nrow(ref_data))[g]))==1)
    expected <- ref_projection(d,as.matrix(d["x"]),model.matrix(~site+z,d))
    equal(a$bootstrap$draws[j,],expected,tol=1e-8)
    equal(a$bootstrap$sampling$draw_df[j],nrow(d)-3L)
  }
  for(j in seq_len(ncol(a$bootstrap$draws)))
    equal(as.numeric(confint(a)[j,]),unname(quantile(a$bootstrap$draws[,j],c(.025,.975),type=7)))
  equal(as.numeric(confint(a,parm="kappa",level=.8)),
    unname(quantile(a$bootstrap$draws[,"kappa"],c(.1,.9),type=7)))
  cmp <- lambda_compare(ref_data,~x,~x+z,ref_items,.1,nuisance=~site,bootstrap=cfg)
  one <- lambda_model(ref_data,~x,ref_items,.1,nuisance=~site,bootstrap=cfg)
  equal(cmp$bootstrap$indices,one$bootstrap$indices)
  equal(unname(cmp$bootstrap$draws[,3:5]),unname(one$bootstrap$draws[,3:5]))
  fails(lambda_model(ref_data,~x,ref_items,.1,estimator="population",bootstrap=cfg),"case resampling only")
})

test("nuisance rank loss is diagnosed without replacing bootstrap draws", {
  d <- ref_data; d$rare <- factor(c(rep("common",35),"rare"))
  cfg <- lambda_bootstrap_control(B=79,min_valid=20,keep_draws=TRUE)
  a <- lambda_model(d,~x,ref_items,.1,nuisance=~rare,bootstrap=cfg)
  stopifnot(any(a$bootstrap$failures$scope=="measurement"),
    any(grepl("Nuisance",a$bootstrap$failures$message)),a$bootstrap$attempted==79,
    any(a$bootstrap$ci$n_failed>0),all(is.na(confint(a))))
})

test("nuisance transforms are rebuilt and explicit row ordering is enforced", {
  cfg <- lambda_bootstrap_control(B=9,min_valid=2,keep_draws=TRUE,keep_indices=TRUE)
  a <- lambda_model(ref_data,~x,ref_items,.1,nuisance=~poly(z,2),bootstrap=cfg)
  for(b in seq_len(cfg$B)) {
    d <- ref_data[a$bootstrap$indices[b,],]
    equal(a$bootstrap$draws[b,],ref_projection(d,as.matrix(d["x"]),model.matrix(~poly(z,2),d)))
  }
  B <- model.matrix(~site,ref_data)
  fails(lambda_model(ref_data,~x,ref_items,.1,nuisance=B[36:1,]),"row identities")
  fails(lambda_model(ref_data,~x,ref_items,.1,nuisance=i1~site),"one-sided")
  fails(lambda_model(ref_data,~x,ref_items,.1,nuisance_contrasts=list(site="contr.sum")),"requires")
  fails(lambda_bootstrap_control(resampling="cluster"),"Supply cluster")
  fails(lambda_bootstrap_control(cluster="school"),"Supply cluster")
  fails(lambda_bootstrap_control(percentile_interpolation="invented"),"arg")
})

test("conditional encoded geometry governs F_H after nuisance adjustment", {
  X <- as.matrix(ref_data[c("x","z")]);B <- model.matrix(~site,ref_data)
  residual_x <- qr.resid(qr(B),X)
  volume <- sqrt(det(cor(residual_x)))
  a <- lambda_model(ref_data,X,ref_items,.1,nuisance=B)
  equal(a$parsimony$volume,volume)
  expected <- a$residual_observed+2*xi_eta_numerical(a$temperature,.1)-xi_a(a$temperature,.1)*log(volume)
  equal(a$parsimony$free_energy,expected)
  X[,2] <- X[,2]+2*X[,1]
  b <- lambda_model(ref_data,X,ref_items,.1,nuisance=B)
  equal(a$lambda_h,b$lambda_h)
  stopifnot(abs(a$parsimony$free_energy-b$parsimony$free_energy)>.01)
})

test("paper saturation is supported and nuisance temperature zero stays singular", {
  d <- data.frame(i1=c(-3,-1,1,3),i2=c(-2.9,-1.2,1.2,2.9),g=factor(1:4))
  a <- lambda_model(d,~g,c("i1","i2"),.1)
  equal(a$df_perp,0L);equal(a$kappa,1);equal(a$lambda_h,a$lambda_star)
  d <- ref_data; d[ref_items] <- replicate(3,ref_signal)
  a <- lambda_model(d,~x,ref_items,.1,nuisance=~site,
    bootstrap=lambda_bootstrap_control(B=9,min_valid=2))
  equal(a$status,"singular_temperature");equal(a$bootstrap$attempted,0L)
  stopifnot(all(is.na(confint(a))))
  fails(lambda_measure(ref_data,ref_items,nuisance=diag(36)),"no residual measurement")
})
