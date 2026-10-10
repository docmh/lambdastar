ms <- readRDS(fixture("model-space-baseline.rds"))
ms_data <- ms$data; ms_items <- ms$items; ms_cfg <- ms$config
ms_same <- function(a,b) {
  equal(a$estimates,b$estimates)
  equal(a$parsimony,b$parsimony, tol=1e-9)
  equal(a$bootstrap$draws,b$bootstrap$draws)
  equal(a$bootstrap$indices,b$bootstrap$indices)
  equal(a$bootstrap$ci,b$bootstrap$ci)
}

test("0.4.1 baseline points, geometry, draws and intervals remain compatible", {
  for (name in names(ms$forms)) {
    f <- ms$forms[[name]]
    a <- lambda_model(ms_data,f,ms_items,1/12,bootstrap=ms_cfg)
    ms_same(a,ms$fits[[name]])
    object <- lambda_design(ms_data,f)
    ms_same(lambda_model(ms_data,object,ms_items,1/12,bootstrap=ms_cfg),a)
  }
  for (scale in c("mean","item")) {
    estimator <- if(scale=="mean") "population" else "paper"
    a <- lambda_model(ms_data,lambda_design(ms_data,~x+z),ms_items,1/12,
      estimator=estimator,measurement_scale=scale,bootstrap=ms_cfg)
    ms_same(a,ms$fits[[if(scale=="mean") "population" else "item"]])
  }
  a <- lambda_model(ms_data,lambda_design(ms_data,~g*x,nuisance=~site+z),
    ms_items,1/12,bootstrap=ms_cfg)
  ms_same(a,ms$fits$adjusted)
  cmp <- lambda_compare(ms_data,lambda_design(ms_data,~x),lambda_design(ms_data,~x+z),
    ms_items,1/12,bootstrap=ms_cfg)
  equal(cmp$delta_F,ms$comparison$delta_F)
  equal(cmp$bootstrap$draws,ms$comparison$bootstrap$draws)
})

test("saved designs retain coding independently of global contrasts and caller functions", {
  object <- local({ temporary <- 123; lambda_design(ms_data,~g*poly(x,2)+z) })
  expected <- lambda_model(ms_data,object,ms_items,1/12,bootstrap=ms_cfg)
  path <- tempfile(fileext=".rds");saveRDS(object,path)
  restored <- readRDS(path);unlink(path)
  local({ old <- options(contrasts=c("contr.helmert","contr.poly"));on.exit(options(old))
    ms_same(lambda_model(ms_data,restored,ms_items,1/12,bootstrap=ms_cfg),expected) })
  fails(lambda_design(ms_data,~rnorm(x)),"deterministic")
  f <- function(x) x^2
  fails(lambda_design(ms_data,~f(x)),"encoded matrix")
  fails(lambda_design(ms_data,~g,contrasts=list(g=function(n)contr.sum(n))),"Portable")
  d <- ms_data; contrasts(d$g) <- contr.helmert(3)
  a <- lambda_model(d,~g,ms_items,1/12)
  ms_same(lambda_model(d,lambda_design(d,~g),ms_items,1/12),a)
})

test("lm and matrix design objects preserve their encoded models", {
  d <- ms_data;d$y <- rowMeans(d[ms_items])
  fitted <- lm(y~g*x+z,d,contrasts=list(g="contr.sum"))
  a <- lambda_model(d,fitted,ms_items,1/12,bootstrap=ms_cfg)
  ms_same(lambda_model(d,lambda_design(d,fitted),ms_items,1/12,bootstrap=ms_cfg),a)
  X <- model.matrix(fitted)[,-1,drop=FALSE]
  design <- lambda_design(d,X)
  equal(design$X_H,X)
  ms_same(lambda_model(d,design,ms_items,1/12,bootstrap=ms_cfg),a)
})

test("bound designs reject stale predictors and row order but permit new measurements", {
  object <- lambda_design(ms_data,~g+x)
  fails(lambda_model(ms_data[90:1,],object,ms_items,1/12),"row identities")
  d <- ms_data;d$x[2] <- d$x[2]+.1
  fails(lambda_model(d,object,ms_items,1/12),"input columns changed")
  d <- ms_data;d$i1 <- d$i1+.1*d$x
  ms_same(lambda_model(d,object,ms_items,1/12),lambda_model(d,~g+x,ms_items,1/12))
  object$schema_version <- 99L
  fails(lambda_model(ms_data,object,ms_items,1/12),"schema")
})

# Independent full/restricted fits and item-pair temperature reference.
ms_reference <- function(d,form,target,adjust) {
  mm <- model.matrix(form,d); ass <- attr(mm,"assign");labels <- attr(terms(form),"term.labels")
  X <- mm[,ass %in% match(target,labels),drop=FALSE]
  B <- mm[,ass==0L | ass %in% match(adjust,labels),drop=FALSE]
  Y <- as.matrix(d[ms_items]);qb <- qr(B);df <- nrow(Y)-qb$rank
  pairs <- combn(ncol(Y),2)
  T <- mean(apply(pairs,2,function(j)sum(qr.resid(qb,Y[,j[1]]-Y[,j[2]])^2)/df))/(2*ncol(Y))
  y <- rowMeans(Y)
  restricted <- lm.fit(B,y); full <- lm.fit(cbind(B,X),y)
  total <- max(0,sum(restricted$residuals^2)/T-df)
  captured <- sum(restricted$residuals^2)-sum(full$residuals^2)
  h <- min(total,max(0,captured/T-(full$rank-restricted$rank)))
  x <- qr.resid(qb,X);volume <- if(ncol(x)==1L) 1 else sqrt(det(cor(x)))
  list(values=c(temperature=T,lambda_star=total,lambda_h=h,lambda_perp=total-h,
    kappa=if(total>0)h/total else NA_real_),rss=sum(full$residuals^2),volume=volume,
    d=df,X=X,B=B)
}

test("conditional blocks match restricted fits and actual encoded geometry on unequal cells", {
  form <- ~g*x+z
  object <- lambda_hypothesis(ms_data,form,"g:x",c("g","x","z"))
  a <- lambda_model(ms_data,object,ms_items,1/12,bootstrap=ms_cfg)
  r <- ms_reference(ms_data,form,"g:x",c("g","x","z"))
  equal(a$estimates,r$values);equal(a$d,r$d)
  equal(a$parsimony$rss,r$rss);equal(a$parsimony$volume,r$volume)
  equal(a$parsimony$free_energy,r$rss+2*xi_eta_numerical(r$values[1],1/12)-
    xi_a(r$values[1],1/12)*log(r$volume))
  ms_same(a,lambda_model(ms_data,r$X,ms_items,1/12,nuisance=r$B,bootstrap=ms_cfg))
  equal(a$inference$hypothesis$total_structure,"after_nuisance_adjustment")
  stopifnot(any(grepl("after nuisance",capture.output(print(a)))))
  for (b in seq_len(ms_cfg$B)) {
    ix <- a$bootstrap$indices[b,]
    equal(a$bootstrap$draws[b,],ms_reference(ms_data[ix,],form,"g:x",c("g","x","z"))$values)
  }
})

test("block transformations are rebuilt after shared selection and inside each draw", {
  d <- ms_data;d$i2[4] <- NA; d$z[7] <- NA
  form <- ~g*poly(x,2)+z
  target <- tail(attr(terms(form),"term.labels"),1)
  adjust <- setdiff(attr(terms(form),"term.labels"),target)
  design <- lambda_hypothesis(d,form,target,adjust)
  a <- lambda_model(d,design,ms_items,1/12,na.action="complete",bootstrap=ms_cfg)
  equal(a$cases$excluded,c(4L,7L))
  equal(a$estimates,ms_reference(d[-c(4,7),],form,target,adjust)$values)
  for (b in seq_len(ms_cfg$B))
    equal(a$bootstrap$draws[b,],ms_reference(d[a$bootstrap$indices[b,],],form,target,adjust)$values)
})

test("additional nuisance and conditional cluster bootstrap keep common row provenance", {
  d <- ms_data;d$school <- rep(letters[1:10],9)
  cfg <- lambda_bootstrap_control(B=13,min_valid=2,seed=14,resampling="cluster",
    cluster="school",keep_draws=TRUE,keep_indices=TRUE)
  object <- lambda_hypothesis(d,~g*x+z,"g:x",c("g","x","z"),nuisance=~site)
  a <- lambda_model(d,object,ms_items,1/12,bootstrap=cfg)
  full <- ~site+g*x+z;adjust <- c("site","g","x","z")
  equal(a$estimates,ms_reference(d,full,"g:x",adjust)$values)
  for (b in seq_len(cfg$B)) {
    ix <- a$bootstrap$indices[[b]]
    equal(a$bootstrap$draws[b,],ms_reference(d[ix,],full,"g:x",adjust)$values)
  }
})

test("comparison requires explicit compatible references and preserves paired draws", {
  a <- lambda_design(ms_data,~x,nuisance=~g)
  b <- lambda_design(ms_data,~x+z,nuisance=~g)
  cmp <- lambda_compare(ms_data,a,b,ms_items,1/12,bootstrap=ms_cfg)
  old <- lambda_compare(ms_data,~x,~x+z,ms_items,1/12,nuisance=~g,bootstrap=ms_cfg)
  equal(cmp$delta_F,old$delta_F);equal(cmp$bootstrap$draws,old$bootstrap$draws)
  fails(lambda_compare(ms_data,a,lambda_design(ms_data,~x+z),ms_items,1/12),"same nuisance")
  fails(lambda_model(ms_data,a,ms_items,1/12,nuisance=~z),"conflicts")
  fails(lambda_model(ms_data,~x,ms_items,1/12,nuisance=a),"not a lambda_design")
  fails(lambda_design(ms_data,~x,nuisance=a),"not a lambda_design")
  point <- cmp$model_a;point$bootstrap <- NULL
  ms_same(lambda_model(ms_data,a,ms_items,1/12,nuisance=~g),point)
  fails(lambda_model(ms_data,a,ms_items,1/12,estimator="population"),"intercept adjustment")
  fails(lambda_model(ms_data,lambda_design(ms_data,~g),ms_items,1/12,
    contrasts=list(g="contr.helmert")),"conflict")
})

test("conditional comparisons share nuisance while retaining distinct target columns", {
  form <- ~g+x+z
  a <- lambda_hypothesis(ms_data,form,"x","g")
  b <- lambda_hypothesis(ms_data,form,c("x","z"),"g")
  cmp <- lambda_compare(ms_data,a,b,ms_items,1/12,bootstrap=ms_cfg)
  old <- lambda_compare(ms_data,~x,~x+z,ms_items,1/12,nuisance=~g,bootstrap=ms_cfg)
  equal(cmp$delta_F,old$delta_F);equal(cmp$bootstrap$draws,old$bootstrap$draws)
  mixed <- lambda_compare(ms_data,~x,b,ms_items,1/12,bootstrap=ms_cfg)
  equal(mixed$delta_F,cmp$delta_F)
  equal(mixed$bootstrap$draws,cmp$bootstrap$draws)
})

test("shared point evaluation prepares nuisance and thermal terms only once", {
  local({
    ns <- asNamespace("lambdastar")
    counts <- new.env();counts$nuisance <- counts$thermal <- 0L
    # Tracers are scoped to this test and removed even if a check fails.
    trace(".nuisance_space",tracer=substitute(assign("nuisance", get("nuisance",envir=E)+1L,envir=E),list(E=counts)),
      where=ns,print=FALSE)
    on.exit(untrace(".nuisance_space",where=ns),add=TRUE)
    trace(".thermal_terms",tracer=substitute(assign("thermal", get("thermal",envir=E)+1L,envir=E),list(E=counts)),
      where=ns,print=FALSE)
    on.exit(untrace(".thermal_terms",where=ns),add=TRUE)
    a <- lambda_design(ms_data,~x,nuisance=~g)
    b <- lambda_design(ms_data,~x+z,nuisance=~g)
    invisible(lambda_compare(ms_data,a,b,ms_items,1/12))
    equal(counts$nuisance,1L);equal(counts$thermal,1L)
  })
})

test("term selections, omissions and rank failures have explicit meaning", {
  fails(lambda_hypothesis(ms_data,~g*x,"x","x"),"disjoint")
  fails(lambda_hypothesis(ms_data,~g*x,"missing",character()),"exact")
  fails(lambda_hypothesis(ms_data,~g*x,character(),character()),"exact")
  fails(lambda_hypothesis(ms_data,~g*x,"x"),"adjust")
  object <- lambda_hypothesis(ms_data,~g*x+z,"x",character())
  equal(object$hypothesis$omitted,c("g","z","g:x"))
  ms_same(lambda_model(ms_data,object,ms_items,1/12),lambda_model(ms_data,~x,ms_items,1/12))
  d <- ms_data;d$copy <- d$x
  fails(lambda_model(d,lambda_hypothesis(d,~x+copy,"x","copy"),ms_items,1/12),"redundant")
  d$rare <- factor(c(rep("common",89),"rare"))
  object <- lambda_hypothesis(d,~x+rare,"x","rare")
  a <- lambda_model(d,object,ms_items,1/12,bootstrap=ms_cfg)
  stopifnot(nrow(a$bootstrap$failures)>0,all(is.na(confint(a))),a$bootstrap$attempted==29L)
})

test("designs preserve singular temperatures and raw coding sensitivity", {
  d <- ms_data;d[ms_items] <- replicate(4,2*d$x)
  a <- lambda_model(d,lambda_design(d,~x),ms_items,1/12,bootstrap=ms_cfg)
  equal(a$status,"singular_temperature");equal(a$bootstrap$attempted,0L)
  d <- ms_data;X <- as.matrix(d[c("x","z")]);Y <- cbind(X[,1],X[,2]+3*X[,1])
  a <- lambda_model(d,lambda_design(d,X),ms_items,1/12)
  b <- lambda_model(d,lambda_design(d,Y),ms_items,1/12)
  equal(a$lambda_h,b$lambda_h)
  stopifnot(abs(a$parsimony$free_energy-b$parsimony$free_energy)>.01)
})
