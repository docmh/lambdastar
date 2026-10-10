# E6a: independent uncentered and affine least-squares references.
ref_data <- function() {
  set.seed(6801)
  d <- data.frame(x=rnorm(60),z=rnorm(60),group=factor(rep(c("a","b"),each=30)))
  S <- 12+.6*d$x+.4*d$z
  for(k in 1:4) d[[paste0("i",k)]] <- S+rnorm(60)
  d
}
ref_items <- paste0("i",1:4)
ref_cfg <- lambda_bootstrap_control(B=29,min_valid=2,seed=371,keep_draws=TRUE,keep_indices=TRUE)

test("reference mean uses N degrees of freedom and preserves its direction", {
  d <- ref_data();Y <- as.matrix(d[ref_items]);n <- nrow(Y);k <- ncol(Y);r <- 10
  T <- mean(combn(k,2,function(j)sum((Y[,j[1]]-Y[,j[2]])^2)))/(2*n*k)
  f <- lambda_model(d,lambda_reference(d,~1,value=r),ref_items,1/12)
  equal(f$d,n);equal(f$nuisance$rank,0L);equal(f$df_h,1L)
  equal(f$temperature,T);equal(f$lambda_h_raw,n*(mean(Y)-r)^2/T-1)
  equal(f$lambda_star_raw,sum((rowMeans(Y)-r)^2)/T-n)
  equal(f$mean_departure,mean(Y)-r);stopifnot(!f$model$intercept_removed)
  equal(f$parsimony$rss,sum((rowMeans(Y)-mean(Y))^2))
  equal(lambda_measure(d,ref_items,reference=r)$estimates[1:2],f$estimates[1:2])
  far <- lambda_measure(d,ref_items,reference=1e20)
  equal(far$temperature,f$temperature)
})

test("reference null and free mean share measurement and numerical F components", {
  d <- ref_data();Y <- as.matrix(d[ref_items]);r <- 10
  h0 <- lambda_reference(d,~0,value=r);h1 <- lambda_reference(d,~1,value=r)
  cmp <- lambda_compare(d,h0,h1,ref_items,1/12)
  equal(cmp$model_a$df_h,0L);equal(cmp$model_a$lambda_h,0)
  equal(cmp$model_a$parsimony$rss,sum((rowMeans(Y)-r)^2))
  equal(cmp$model_b$parsimony$free_energy,cmp$model_b$parsimony$rss+cmp$model_b$parsimony$xi_eta)
  equal(cmp$delta_F,cmp$model_b$parsimony$free_energy-cmp$model_a$parsimony$free_energy)
  manual <- thermodynamic_model_state(rowMeans(Y),matrix(1,nrow(Y),1),cmp$measurement$temperature,
    1/12,reference=r)
  equal(manual$free_energy,cmp$model_b$parsimony$free_energy)
})

test("reference formula and matrix use uncentered encoded geometry", {
  d <- ref_data();X <- cbind(`(Intercept)`=1,x=d$x,z=d$z)
  f <- lambda_model(d,lambda_reference(d,~1+x+z,value=10),ref_items,1/12)
  g <- lambda_model(d,lambda_reference(d,X,value=10),ref_items,1/12)
  equal(f$estimates,g$estimates);equal(f$parsimony$free_energy,g$parsimony$free_energy)
  U <- sweep(X,2,sqrt(colSums(X^2)),"/")
  equal(f$parsimony$volume,sqrt(det(crossprod(U))))
  X[,3] <- X[,3]+2*X[,2]
  h <- lambda_model(d,lambda_reference(d,X,value=10),ref_items,1/12)
  equal(f$estimates,h$estimates,tol=1e-8)
  stopifnot(abs(f$parsimony$free_energy-h$parsimony$free_energy)>1e-5)
})

test("translation and paired sign preserve magnitudes without hiding direction", {
  d <- ref_data();f <- lambda_model(d,lambda_reference(d,~1,value=10),ref_items,1/12)
  shifted <- d;shifted[ref_items] <- shifted[ref_items]+23
  g <- lambda_model(shifted,lambda_reference(shifted,~1,value=33),ref_items,1/12)
  equal(f$estimates,g$estimates,tol=1e-8);equal(f$parsimony$free_energy,g$parsimony$free_energy,tol=1e-8)
  h <- lambda_model(d,lambda_reference(d,~1,value=12),ref_items,1/12)
  stopifnot(f$lambda_h>h$lambda_h)
})

test("explicit nuisance retains only the declared reference space", {
  d <- ref_data();Y <- as.matrix(d[ref_items]);B <- matrix(d$z,ncol=1)
  f <- lambda_model(d,lambda_reference(d,~1+x,value=10,nuisance=~0+z),ref_items,1/12)
  Z <- qr.resid(qr(B),Y-10);z <- rowMeans(Z);A <- qr.resid(qr(B),cbind(1,d$x))
  T <- sum((Z-z)^2)/((nrow(d)-1)*3*4);fit <- lm.fit(A,z)
  equal(f$d,nrow(d)-1L);equal(f$temperature,T)
  equal(f$lambda_h_raw,sum(fit$fitted.values^2)/T-2)
  equal(f$parsimony$rss,sum(fit$residuals^2))
  fails(lambda_model(d,lambda_reference(d,~1,value=10,nuisance=~z),ref_items,1/12),"redundant")
})

test("explicit affine L beta=c matches independent constrained least squares", {
  d <- ref_data();X <- cbind(1,d$x,d$z);D <- matrix(c(1,2,0)/5,ncol=1)
  N0 <- cbind(c(-2,1,0),c(0,0,1));r <- drop(X%*%D*9);B <- X%*%N0;H <- X%*%D
  f <- lambda_model(d,lambda_reference(d,H,value=r,nuisance=B),ref_items,1/12)
  y <- rowMeans(d[ref_items])-r
  rss0 <- sum(lm.fit(B,y)$residuals^2);rss1 <- sum(lm.fit(X,y)$residuals^2)
  equal(f$lambda_h_raw,(rss0-rss1)/f$temperature-1,tol=1e-8)
  equal(f$parsimony$rss,rss1,tol=1e-8)
})

test("fixed-reference case draws carry row-specific origins and exact percentiles", {
  d <- ref_data();r <- 10+.2*d$x
  f <- lambda_model(d,lambda_reference(d,~1+x,value=r),ref_items,1/12,bootstrap=ref_cfg)
  for(j in c(1,15,29)) {
    ix <- f$bootstrap$indices[j,];s <- d[ix,]
    g <- lambda_model(s,lambda_reference(s,~1+x,value=r[ix]),ref_items,1/12)
    equal(g$estimates,f$bootstrap$draws[j,],tol=1e-8)
  }
  for(j in seq_len(ncol(f$bootstrap$draws))) {
    v <- f$bootstrap$draws[,j]
    q <- if(length(unique(v))==1) rep(v[1],2) else suppressWarnings(boot::boot.ci(
      structure(list(t0=f$estimates[j],t=matrix(v,ncol=1),R=length(v),sim="ordinary"),class="boot"),type="perc",index=1)$percent[1,4:5])
    equal(as.numeric(confint(f)[j,]),unname(q))
  }
})

test("reference cluster draws recompute N and resample the origin vector", {
  d <- ref_data();d$cluster <- rep(seq_len(12),times=rep(c(3,7),6));r <- 10+.2*d$x
  cfg <- lambda_bootstrap_control(B=19,min_valid=2,seed=45,resampling="cluster",cluster="cluster",keep_draws=TRUE,keep_indices=TRUE)
  f <- lambda_model(d,lambda_reference(d,~1,value=r),ref_items,1/12,bootstrap=cfg)
  stopifnot(length(unique(f$bootstrap$sampling$draw_sizes))>1)
  equal(f$bootstrap$sampling$draw_sizes,f$bootstrap$sampling$draw_df)
  for(j in c(1,10,19)) {
    ix <- f$bootstrap$indices[[j]];s <- d[ix,]
    g <- lambda_model(s,lambda_reference(s,~1,value=r[ix]),ref_items,1/12)
    equal(g$estimates,f$bootstrap$draws[j,],tol=1e-8)
  }
})

test("reference comparison rejects changed origins and implicit centered partners", {
  d <- ref_data();h <- lambda_reference(d,~1,value=10)
  fails(lambda_compare(d,h,lambda_reference(d,~1,value=11),ref_items,1/12),"same fixed reference")
  fails(lambda_compare(d,h,~1,ref_items,1/12),"same fixed reference")
  fails(lambda_model(d,h,ref_items,1/12,estimator="population"),"centered target")
  fails(lambda_contrast(d,h,c(x=1)),"Centered contrast constructors")
  fails(lambda_model(d[60:1,],h,ref_items,1/12),"row identities")
  fails(lambda_reference(d,~1,value=NA_real_),"finite numeric")
  fails(lambda_reference(d,~1,value=setNames(rep(0,60),rev(rownames(d)))),"names must match")
  broken <- h;broken$reference$value <- NULL
  fails(lambda_model(d,broken,ref_items,1/12),"Reference specification")
})

test("paired differences are recomputed on the supplied pair coding", {
  d <- ref_data();pre <- ref_items;post <- paste0("post",1:4)
  set.seed(33);for(j in 1:4)d[[post[j]]] <- d[[pre[j]]]+2+.2*d$x+rnorm(60)
  f <- lambda_paired(d,pre,post,1/12,value=.5,bootstrap=ref_cfg)
  delta <- as.matrix(d[post])-as.matrix(d[pre]);temp <- data.frame(delta)
  g <- lambda_model(temp,lambda_reference(temp,~1,value=.5),names(temp),1/12,bootstrap=ref_cfg)
  equal(f$estimates,g$estimates);equal(f$bootstrap$draws,g$bootstrap$draws)
  equal(f$bootstrap$indices,g$bootstrap$indices)
  h <- lambda_paired(d,post,pre,1/12,value=-.5)
  equal(f$estimates,h$estimates);equal(f$mean_departure,-h$mean_departure)
  equal(f$parsimony$free_energy,h$parsimony$free_energy)
  for(j in c(1,29)) {
    ix <- f$bootstrap$indices[j,]
    replay <- lambda_paired(d[ix,],pre,post,1/12,value=.5)
    equal(replay$estimates,f$bootstrap$draws[j,])
  }
})

test("paired selection keeps complete persons and response units", {
  d <- ref_data();post <- paste0("post",1:4)
  set.seed(71);for(j in 1:4)d[[post[j]]] <- d[[ref_items[j]]]+1+rnorm(60)
  d$i1[3] <- NA;d$post4[7] <- NA;d$x[9] <- NA
  fails(lambda_paired(d,ref_items,post,1/12,model=~1+x),"Missing values")
  f <- lambda_paired(d,ref_items,post,1/12,model=~1+x,na.action="complete")
  equal(f$cases$excluded,c(3L,7L,9L));equal(f$d,57L)
  g <- lambda_paired(d,ref_items,post,1/12,model=~1+x,na.action="complete",measurement_scale="item")
  equal(f$lambda_star,g$lambda_star*4);equal(f$kappa,g$kappa)
  equal(f$parsimony$free_energy,g$parsimony$free_energy)
  h <- lambda_paired(d,ref_items,post,1/48,model=~1+x,na.action="complete")
  equal(f$estimates,h$estimates);stopifnot(abs(f$parsimony$free_energy-h$parsimony$free_energy)>1e-5)
  fails(lambda_paired(d,ref_items[1],post[1],1/12),"at least two")
  fails(lambda_paired(d,ref_items,ref_items,1/12),"matching measurement pairs")
})

test("reference zero-signal and zero-temperature cases keep their meaning", {
  d <- data.frame(i1=rep(1,8),i2=rep(-1,8));h <- lambda_reference(d,~1,value=0)
  f <- lambda_model(d,h,c("i1","i2"),1/12,bootstrap=ref_cfg)
  stopifnot(f$temperature>0,f$lambda_star==0,is.na(f$kappa),all(f$bootstrap$draws[,"lambda_star"]==0),
    all(is.na(f$bootstrap$draws[,"kappa"])))
  d$i2 <- d$i1
  g <- lambda_model(d,lambda_reference(d,~1),c("i1","i2"),1/12,bootstrap=ref_cfg)
  stopifnot(g$status=="singular_temperature",is.na(g$parsimony$free_energy),g$bootstrap$attempted==0L)
  d$i2 <- d$i1+c(rep(1e-6,4),rep(-1e-6,4))
  tiny <- lambda_model(d,lambda_reference(d,~1),c("i1","i2"),1/12)
  stopifnot(tiny$temperature>0,tiny$temperature<1e-12,is.finite(tiny$parsimony$free_energy))
})

test("reference and centered saturated ranks have their appropriate dimensions", {
  d <- ref_data()[1:4,];X <- diag(4)
  f <- lambda_model(d,lambda_reference(d,X,value=0),ref_items,1/12)
  equal(f$df_h,4L);equal(f$df_perp,0L);stopifnot(f$residual_observed<1e-20)
  fails(lambda_model(d,X,ref_items,1/12),"rank deficient")
  r <- rep(10,4)
  m <- estimate_lambda_partition_items(as.matrix(d[ref_items]),X,reference=r)
  g <- lambda_model(d,lambda_reference(d,X,value=r),ref_items,1/12)
  equal(m$lambda_h,g$lambda_h);equal(m$temperature,g$temperature)
})

test("reference response rescaling changes temperatures and F but not noncentrality", {
  d <- ref_data();f <- lambda_model(d,lambda_reference(d,~1+x,value=10),ref_items,1/12)
  scaled <- d;scaled[ref_items] <- scaled[ref_items]*2
  g <- lambda_model(scaled,lambda_reference(scaled,~1+x,value=20),ref_items,4/12)
  equal(f$estimates[-1],g$estimates[-1]);equal(g$temperature,4*f$temperature)
  equal(g$parsimony$free_energy,4*f$parsimony$free_energy,tol=1e-8)
})
