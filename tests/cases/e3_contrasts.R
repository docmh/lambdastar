e3 <- readRDS(fixture("e3-baseline.rds"))
e3_data <- e3$data;e3_data$i2[4] <- mean(e3_data$i2,na.rm=TRUE);e3_data$z[7] <- 0
e3_items <- e3$items;e3_cfg <- e3$config
e3_fit <- function(d,h,boot=FALSE,...) lambda_model(d,h,e3_items,1/12,bootstrap=boot,...)
e3_same <- function(a,b) {
  equal(a$estimates,b$estimates,tol=1e-9)
  pa <- a$parsimony;pb <- b$parsimony
  pa$basis <- pb$basis <- NULL
  equal(pa,pb,tol=1e-9)
  equal(a$bootstrap$draws,b$bootstrap$draws,tol=1e-9)
  equal(a$bootstrap$indices,b$bootstrap$indices)
  equal(a$bootstrap$ci,b$bootstrap$ci,tol=1e-9)
}
# Independent SVD null-model and least-squares fits; no package basis helper.
e3_oracle <- function(d,X,L,B=NULL) {
  Y <- as.matrix(d[e3_items]);p <- ncol(X);r <- nrow(L)
  N <- svd(L,nu=0,nv=p)$v[,seq_len(p)>r,drop=FALSE]
  B <- cbind(1,B,X%*%N);q <- qr(B);df <- nrow(d)-q$rank
  T <- mean(combn(ncol(Y),2,function(j)
    sum(qr.resid(q,Y[,j[1]]-Y[,j[2]])^2)/df))/(2*ncol(Y))
  reduced <- lm.fit(B,rowMeans(Y));full <- lm.fit(cbind(1,X,if(ncol(B)>1) B[,-1,drop=FALSE]),rowMeans(Y))
  rss0 <- sum(reduced$residuals^2);rss1 <- sum(full$residuals^2)
  total <- max(0,rss0/T-df);h <- min(total,max(0,(rss0-rss1)/T-r))
  # For a test oracle only: direct coefficient right inverse via a small solve.
  D <- t(L)%*%solve(L%*%t(L));A <- qr.resid(q,X%*%D)
  volume <- if(r==1) 1 else sqrt(det(cor(A)))
  list(values=c(temperature=T,lambda_star=total,lambda_h=h,lambda_perp=total-h,
    kappa=if(total>0)h/total else NA_real_),rss=rss1,volume=volume,df=df)
}
e3_check_oracle <- function(a,ref) {
  equal(a$estimates,ref$values,tol=1e-8)
  equal(a$parsimony$rss,ref$rss,tol=1e-8)
  equal(a$parsimony$volume,ref$volume,tol=1e-8)
}

test("E3 preserves 0.5.0 whole and conditional adaptive-basis outputs", {
  for(i in seq_along(e3$objects)) e3_same(e3_fit(e3$data,e3$objects[[i]],e3_cfg,
    na.action="complete"),e3$fits[[i]])
})

test("general multirow L agrees with independently fitted restricted models", {
  d <- e3_data;X <- model.matrix(~x+z+g,d)[,-1,drop=FALSE]
  L <- rbind(difference=c(x=1,z=-1,gb=0,gc=0),joint=c(x=.5,z=0,gb=1,gc=-1))
  a <- e3_fit(d,lambda_contrast(d,~x+z+g,L))
  e3_check_oracle(a,e3_oracle(d,X,L))
  equal(unname(a$model$L%*%a$model$coefficient_right_inverse),diag(2),tol=1e-9)
  equal(unname(a$model$L%*%a$model$coefficient_null_basis),matrix(0,2,2),tol=1e-9)
  equal(rownames(a$model$L),rownames(L))
})

test("coefficient, matrix and lm sources retain the same conditional geometry", {
  d <- e3_data;d$y <- rowMeans(d[e3_items]);fitted <- lm(y~x+z+g,d)
  L <- c(x=1,z=-1);a <- e3_fit(d,lambda_contrast(d,~x+z+g,L),e3_cfg)
  for(source in list(fitted,model.matrix(fitted)[,-1,drop=FALSE],lambda_design(d,~x+z+g)))
    e3_same(e3_fit(d,lambda_contrast(d,source,L),e3_cfg),a)
  b <- e3_fit(d,lambda_contrast(d,~x+z+g,L,nuisance=~site))
  X <- model.matrix(~x+z+g,d)[,-1,drop=FALSE];L <- matrix(c(1,-1,0,0),1)
  e3_check_oracle(b,e3_oracle(d,X,L,model.matrix(~site,d)[,-1,drop=FALSE]))
})

test("full coefficient restrictions reduce to a whole-model state", {
  d <- e3_data;L <- diag(2);colnames(L) <- c("x","z")
  a <- e3_fit(d,lambda_contrast(d,~x+z,L),e3_cfg)
  e3_same(a,e3_fit(d,~x+z,e3_cfg))
  equal(ncol(a$model$coefficient_null_basis),0L)
})

test("L units and invertible recoding preserve capture but retain F_H geometry", {
  d <- e3_data;L <- rbind(c(x=1,z=0,gb=0,gc=0),c(x=0,z=1,gb=0,gc=0))
  recoded <- rbind(c(1,2),c(0,1))%*%L
  h1 <- lambda_contrast(d,~x+z+g,L);h2 <- lambda_contrast(d,~x+z+g,recoded)
  cmp <- lambda_compare(d,h1,h2,e3_items,1/12,bootstrap=e3_cfg)
  equal(cmp$model_a$estimates,cmp$model_b$estimates,tol=1e-8)
  stopifnot(abs(cmp$model_a$parsimony$volume-cmp$model_b$parsimony$volume)>.1,
            abs(cmp$delta_F)>.01)
  for(i in seq_len(e3_cfg$B)) {
    sample <- d[cmp$bootstrap$indices[i,],]
    X <- model.matrix(~x+z+g,sample)[,-1,drop=FALSE]
    ref <- e3_oracle(sample,X,L)$values
    equal(unname(cmp$bootstrap$draws[i,1:5]),unname(ref),tol=1e-8)
  }
  scaled <- sweep(L,1,c(1e100,1e-100),"*")
  e3_same(e3_fit(d,lambda_contrast(d,~x+z+g,scaled)),e3_fit(d,h1))
  fails(lambda_compare(d,h1,lambda_contrast(d,~x+z+g,c(x=1)),e3_items,1/12),"same nuisance")
})

test("L input errors and nonidentifiability are explicit", {
  d <- e3_data
  fails(lambda_contrast(d,~x+z,c(1,-1)),"named columns")
  fails(lambda_contrast(d,~x,c(x=0)),"zero restriction")
  fails(lambda_contrast(d,~x,c(`(Intercept)`=1)),"non-intercept")
  fails(e3_fit(d,lambda_contrast(d,~x,c(unknown=1))),"Unknown L")
  fails(e3_fit(d,lambda_contrast(d,~x+z,rbind(c(x=1,z=2),c(x=2,z=4)))),"rank deficient")
  fails(e3_fit(d,lambda_contrast(d,~x+I(2*x),c(x=1))),"rank deficient")
  fails(e3_fit(d,lambda_contrast(d,~x+z,c(x=1),nuisance=~x)),"redundant")
  fails(e3_fit(d,lambda_contrast(d,~x,c(x=1)),estimator="population"),"intercept adjustment")
  fails(lambda_contrast(d,lambda_hypothesis(d,~x+z,"x","z"),c(x=1)),"whole-model")
  fails(lambda_contrast(d,lambda_design(d,~x),c(x=1),encoding="fixed"),"already specifies")
})

test("derived null and target are rebuilt for every complete-case bootstrap draw", {
  d <- e3$data;L <- rbind(c(x=1,z=-1,gb=0,gc=0),c(x=0,z=0,gb=1,gc=-1))
  a <- e3_fit(d,lambda_contrast(d,~x+z+g,L,nuisance=~site),e3_cfg,na.action="complete")
  equal(a$cases$excluded,c(4L,7L))
  for(b in seq_len(e3_cfg$B)) {
    sample <- d[a$bootstrap$indices[b,],]
    X <- model.matrix(~x+z+g,sample)[,-1,drop=FALSE]
    ref <- e3_oracle(sample,X,L,model.matrix(~site,sample)[,-1,drop=FALSE])
    equal(a$bootstrap$draws[b,],ref$values,tol=1e-8)
  }
})

test("contrast cluster draws preserve all rows and variable sample sizes", {
  d <- e3_data;d$cluster <- rep(seq_len(12),length.out=nrow(d))
  cfg <- lambda_bootstrap_control(B=11,min_valid=2,seed=72,resampling="cluster",
    cluster="cluster",keep_draws=TRUE,keep_indices=TRUE)
  a <- e3_fit(d,lambda_contrast(d,~x+z,c(x=1,z=-1)),cfg)
  stopifnot(length(unique(lengths(a$bootstrap$indices)))>1)
  for(b in seq_len(cfg$B)) {
    sample <- d[a$bootstrap$indices[[b]],]
    equal(a$bootstrap$draws[b,],e3_oracle(sample,as.matrix(sample[c("x","z")]),
      matrix(c(1,-1),1))$values,tol=1e-8)
  }
})

test("fixed and reevaluated polynomial/spline blocks reproduce their declared draw estimators", {
  d <- e3$data
  for(form in list(~poly(x,2)+z,~splines::ns(x,df=3)+z)) {
    target <- attr(terms(form),"term.labels")[1]
    fixed <- lambda_hypothesis(d,form,target,"z",encoding="fixed")
    adaptive <- lambda_hypothesis(d,form,target,"z",encoding="reevaluate")
    a <- e3_fit(d,fixed,e3_cfg,na.action="complete")
    b <- e3_fit(d,adaptive,e3_cfg,na.action="complete")
    equal(a$estimates,b$estimates);equal(a$encoded_x,b$encoded_x)
    equal(a$bootstrap$indices,b$bootstrap$indices)
    selected <- a$cases$included
    for(i in seq_len(e3_cfg$B)) {
      ix <- a$bootstrap$indices[i,];j <- match(ix,selected);sample <- d[ix,]
      raw <- a$encoded_x[j,,drop=FALSE];rownames(raw) <- rownames(sample)
      nu <- a$encoded_nuisance[j,,drop=FALSE];rownames(nu) <- rownames(sample)
      reference <- e3_fit(sample,lambda_design(sample,raw,nuisance=nu))
      equal(a$bootstrap$draws[i,],reference$estimates,tol=1e-8)
      reference <- e3_fit(sample,lambda_hypothesis(sample,form,target,"z"))
      equal(b$bootstrap$draws[i,],reference$estimates,tol=1e-8)
    }
  }
})

test("fixed adaptive nuisance and whole-model encoding are learned after selection", {
  d <- e3$data;h <- lambda_design(d,~g*poly(x,2),nuisance=~scale(z),encoding="fixed")
  a <- e3_fit(d,h,e3_cfg,na.action="complete")
  for(i in c(1L,9L,19L)) {
    ix <- a$bootstrap$indices[i,];j <- match(ix,a$cases$included);sample <- d[ix,]
    X <- a$encoded_x[j,,drop=FALSE];B <- a$encoded_nuisance[j,,drop=FALSE]
    rownames(X) <- rownames(B) <- rownames(sample)
    equal(a$bootstrap$draws[i,],e3_fit(sample,lambda_design(sample,X,nuisance=B))$estimates)
  }
})

test("simple and explicitly averaged grid effects agree with named coefficient L", {
  d <- e3_data
  grid <- data.frame(g=c("a","b"),x=c(2,2),z=0,row.names=c("a2","b2"))
  a <- e3_fit(d,lambda_marginal_contrast(d,~g*x+z,grid,c(b2=1,a2=-1)),e3_cfg)
  e3_same(a,e3_fit(d,lambda_contrast(d,~g*x+z,c(gb=1,`gb:x`=2)),e3_cfg))
  grid <- data.frame(g=rep(c("a","b"),each=2),x=rep(c(-1,2),2),z=0,
    row.names=c("a_lo","a_hi","b_lo","b_hi"))
  w <- c(a_lo=-.25,a_hi=-.75,b_lo=.25,b_hi=.75)
  a <- e3_fit(d,lambda_marginal_contrast(d,~g*x+z,grid,w))
  e3_same(a,e3_fit(d,lambda_contrast(d,~g*x+z,c(gb=1,`gb:x`=1.25))))
})

test("grid prediction uses learned poly, spline and scale coordinates", {
  d <- e3_data;d$y <- rowMeans(d[e3_items])
  grid <- data.frame(x=c(-1,.6),z=c(.2,.2),row.names=c("lo","hi"));w <- c(lo=-1,hi=1)
  for(form in list(y~poly(x,3)+z,y~splines::ns(x,df=3)+z,y~splines::bs(x,df=4)+z,y~scale(x)+z)) {
    ordinary <- lm(form,d)
    a <- e3_fit(d,lambda_marginal_contrast(d,form,grid,w))
    equal(drop(a$model$L%*%coef(ordinary)[-1]),
      setNames(sum(w*predict(ordinary,grid)),"contrast1"),tol=1e-8)
    mm <- model.matrix(delete.response(terms(ordinary)),model.frame(delete.response(terms(ordinary)),grid))
    if(!grepl("scale",deparse(form))) equal(unname(a$model$L),unname(matrix(w,1)%*%mm[,-1,drop=FALSE]))
  }
})

test("fixed grid contrasts use original learned basis with draw-specific restricted fits", {
  d <- e3$data;grid <- data.frame(x=c(-1,1),z=0,row.names=c("lo","hi"));w <- c(lo=-1,hi=1)
  form <- ~poly(x,3)+z
  a <- e3_fit(d,lambda_marginal_contrast(d,form,grid,w,encoding="fixed"),e3_cfg,na.action="complete")
  b <- e3_fit(d,lambda_marginal_contrast(d,form,grid,w),e3_cfg,na.action="complete")
  equal(a$estimates,b$estimates)
  X <- model.matrix(form,d[a$cases$included,])[,-1,drop=FALSE]
  for(i in seq_len(e3_cfg$B)) {
    ix <- a$bootstrap$indices[i,];sample <- d[ix,];raw <- X[match(ix,a$cases$included),,drop=FALSE]
    rownames(raw) <- rownames(sample)
    equal(a$bootstrap$draws[i,],e3_fit(sample,lambda_contrast(sample,raw,a$model$L))$estimates,tol=1e-8)
    equal(b$bootstrap$draws[i,],e3_fit(sample,lambda_marginal_contrast(sample,form,grid,w))$estimates,tol=1e-8)
  }
})

test("fixed named polynomial coefficients differ from reevaluated coordinates as declared", {
  d <- e3_data;L <- c(`poly(x, 3)1`=1,`poly(x, 3)2`=-1)
  a <- e3_fit(d,lambda_contrast(d,~poly(x,3)+z,L,encoding="fixed"),e3_cfg)
  b <- e3_fit(d,lambda_contrast(d,~poly(x,3)+z,L),e3_cfg)
  equal(a$estimates,b$estimates)
  stopifnot(max(abs(a$bootstrap$draws-b$bootstrap$draws),na.rm=TRUE)>1)
  for(i in c(1L,9L,19L)) {
    ix <- a$bootstrap$indices[i,];sample <- d[ix,]
    X <- model.matrix(~poly(x,3)+z,d)[ix,-1,drop=FALSE];rownames(X) <- rownames(sample)
    equal(a$bootstrap$draws[i,],e3_fit(sample,lambda_contrast(sample,X,L))$estimates,tol=1e-8)
    equal(b$bootstrap$draws[i,],e3_fit(sample,lambda_contrast(sample,~poly(x,3)+z,L))$estimates,tol=1e-8)
  }
})

test("grid validation prevents implicit weighting and unsupported prediction recipes", {
  d <- e3_data;grid <- data.frame(x=c(-1,1),row.names=c("lo","hi"));w <- c(lo=-1,hi=1)
  fails(lambda_marginal_contrast(d,~x,grid,c(lo=1,hi=1)),"sum to zero")
  fails(lambda_marginal_contrast(d,~x,grid,c(LO=-1,hi=1)),"every grid row")
  fails(e3_fit(d,lambda_marginal_contrast(d,~x+z,grid,w)),"all full-model predictors")
  for(form in list(~I(x-mean(x)),~I(poly(x,2))))
    fails(e3_fit(d,lambda_marginal_contrast(d,form,grid,w)),"prediction recipes")
  grid$g <- c("a","unknown")
  fails(e3_fit(d,lambda_marginal_contrast(d,~g+x,grid,w)),"new level")
  fails(lambda_marginal_contrast(d,as.matrix(d["x"]),grid,w),"formula/lm")
})

test("Fourier and explicitly fixed-knot piecewise/spline recipes equal encoded matrices", {
  d <- e3_data
  for(form in list(~sin(2*pi*x)+cos(2*pi*x),~x+I(pmax(x-.5,0)),
      ~splines::ns(x,knots=c(-.5,.5),Boundary.knots=c(-4,4)))) {
    h <- lambda_design(d,form,encoding="fixed");X <- model.matrix(form,d)[,-1,drop=FALSE]
    e3_same(e3_fit(d,h,e3_cfg),e3_fit(d,X,e3_cfg))
  }
})

test("rank failures remain original draw failures without replacement", {
  d <- e3_data[1:12,];d$x <- c(1,rep(0,11));d$z <- seq_len(12)
  cfg <- lambda_bootstrap_control(B=31,min_valid=2,seed=2,keep_draws=TRUE,keep_indices=TRUE)
  a <- e3_fit(d,lambda_contrast(d,~x+z,c(x=1)),cfg)
  absent <- apply(a$bootstrap$indices,1,function(ix)!1L %in% ix)
  stopifnot(any(absent),all(is.na(a$bootstrap$draws[absent,])),
    nrow(a$bootstrap$draws)==31L,any(grepl("rank|constant",a$bootstrap$failures$message)))
})

test("multirow grid recoding shares one null reference and retains factor coding", {
  d <- e3_data;d$flag <- d$x>0
  grid <- data.frame(g=c("a","b","c"),x=0,z=0,flag=TRUE,row.names=c("a","b","c"))
  w <- rbind(c(a=-1,b=1,c=0),c(a=-1,b=0,c=1))
  a <- lambda_marginal_contrast(d,~g*x+z+flag,grid,w)
  b <- lambda_marginal_contrast(d,~g*x+z+flag,grid,rbind(c(1,2),c(0,1))%*%w)
  cmp <- lambda_compare(d,a,b,e3_items,1/12,bootstrap=e3_cfg)
  equal(cmp$model_a$estimates,cmp$model_b$estimates,tol=1e-8)
  equal(cmp$model_a$model$L[,"flagTRUE"],setNames(c(0,0),c("contrast1","contrast2")))
})

test("contrast singular temperature and undefined kappa retain their meanings", {
  d <- e3_data;d[e3_items] <- rep(list(2*d$x+d$z),4)
  a <- e3_fit(d,lambda_contrast(d,~x+z,c(x=1)),e3_cfg)
  equal(a$status,"singular_temperature");equal(a$temperature,0)
  stopifnot(all(is.na(a$estimates[-1])))
  d[e3_items] <- cbind(d$x,-d$x,d$z,-d$z)
  a <- e3_fit(d,lambda_contrast(d,~x+z,c(x=1)))
  equal(a$lambda_star,0);equal(a$lambda_h,0);stopifnot(is.na(a$kappa))
})
