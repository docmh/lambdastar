# Source, then run lambdastar_contrasts_example(). Sourcing alone runs no analysis.
# B=5000 gives the ordinary default percentile intervals; B=199 is a short lesson.
lambdastar_contrasts_example <- function(B=5000L) {
  old <- if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) .Random.seed else NULL
  kind <- RNGkind()
  on.exit({do.call(RNGkind,as.list(kind));if(is.null(old)) {
    if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)
  } else assign(".Random.seed",old,envir=.GlobalEnv)},add=TRUE)
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(3301)
  d <- data.frame(group=factor(rep(c("a","b","c"),c(35,50,65))),
                  x=rnorm(150),z=rnorm(150))
  signal <- 2*d$x+.7*d$z+1.2*(d$group=="b")+.5*(d$group=="c")+
    .8*d$x*(d$group=="b")+.4*(d$x^2-1)
  items <- paste0("i",1:4);d[items] <- replicate(4,signal+rnorm(150))
  cfg <- lambdastar::lambda_bootstrap_control(B=B,min_valid=min(B,200L),
    keep_draws=TRUE,keep_indices=TRUE)
  evaluate <- function(h,bootstrap=FALSE) lambdastar::lambda_model(d,h,items,1/12,bootstrap=bootstrap)

  # 1. Inspect the full model's actual coefficient names before specifying L.
  full_design <- lambdastar::lambda_design(d,~group*x+z,contrasts=list(group="contr.treatment"))
  full <- evaluate(full_design)
  coefficient_names <- colnames(full$encoded_x)

  # 2. Jointly ask for b-a and c-a at x=1. Unnamed/omitted coefficients get zero.
  L <- rbind(b_minus_a=c(groupb=1,groupc=0,`groupb:x`=1,`groupc:x`=0),
             c_minus_a=c(groupb=0,groupc=1,`groupb:x`=0,`groupc:x`=1))
  hypothesis <- lambdastar::lambda_contrast(d,full_design,L)
  joint <- evaluate(hypothesis,cfg)

  # 3. The same settings can be stated directly; z is held explicitly at zero.
  grid <- data.frame(group=c("a","b","c"),x=1,z=0,row.names=c("a1","b1","c1"))
  weights <- rbind(b_minus_a=c(a1=-1,b1=1,c1=0),c_minus_a=c(a1=-1,b1=0,c1=1))
  via_grid <- evaluate(lambdastar::lambda_marginal_contrast(d,full_design,grid,weights))
  stopifnot(isTRUE(all.equal(joint$estimates,via_grid$estimates,tolerance=1e-9)))

  # 4. An explicit 25%/75% average over x=-1 and x=2, for the b-a difference.
  average_grid <- data.frame(group=rep(c("a","b"),each=2),x=rep(c(-1,2),2),z=0,
    row.names=c("a_lo","a_hi","b_lo","b_hi"))
  averaged <- evaluate(lambdastar::lambda_marginal_contrast(d,full_design,average_grid,
    c(a_lo=-.25,a_hi=-.75,b_lo=.25,b_hi=.75)))

  # 5. A recoding of those same two restrictions has the same null model.
  recoded <- rbind(c(1,2),c(0,1))%*%L
  comparison <- lambdastar::lambda_compare(d,hypothesis,
    lambdastar::lambda_contrast(d,full_design,recoded),items,1/12)
  # Capture agrees; the ordered contrast geometry (hence F_H) can differ.

  # 6. Fix or relearn orthogonal polynomial coordinates under person resampling.
  # These named coefficient restrictions depend on the chosen coordinates.
  polynomial_L <- c(`poly(x, 3)1`=1,`poly(x, 3)2`=-1)
  fixed_design <- lambdastar::lambda_contrast(d,~poly(x,3)+z,polynomial_L,encoding="fixed")
  adaptive_design <- lambdastar::lambda_contrast(d,~poly(x,3)+z,polynomial_L,encoding="reevaluate")
  fixed <- evaluate(fixed_design,cfg);adaptive <- evaluate(adaptive_design,cfg)
  stopifnot(identical(fixed$bootstrap$indices,adaptive$bootstrap$indices))

  # 7. Fully prespecified basis parameters: no sample-learned knot or period.
  recipes <- list(spline=~splines::ns(x,knots=c(-.5,.5),Boundary.knots=c(-4,4)),
    fourier=~sin(2*pi*x)+cos(2*pi*x),piecewise=~x+I(pmax(x-.5,0)))
  bases <- lapply(recipes,function(f)evaluate(lambdastar::lambda_design(d,f,encoding="fixed")))
  list(data=d,items=items,signal=signal,coefficient_names=coefficient_names,
    full=full,hypothesis=hypothesis,joint=joint,grid=grid,weights=weights,via_grid=via_grid,
    averaged=averaged,comparison=comparison,fixed_design=fixed_design,
    adaptive_design=adaptive_design,fixed=fixed,adaptive=adaptive,bases=bases)
}
