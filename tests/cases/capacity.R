test("numerical eta agrees with independent 90-digit Decimal reference", {
  ref <- read.csv(fixture("eta-reference.csv"))
  relative(eta_numerical(ref$ratio,1),ref$eta,5e-10)
})

test("tiny eta follows analytic three-bin limit instead of 0.5", {
  r <- c(1e-300,1e-100,1e-24,1e-16,1e-12,1e-8)
  p <- r/3
  expected <- -(1-p)*log1p(-p)-p*(log(p)-log(2))
  actual <- eta_numerical(r,1,details=TRUE)
  relative(actual$eta,expected,1e-10)
  stopifnot(all(diff(actual$eta)>0),all(abs(actual$relative_moment_error)<1e-9),
            all(actual$tail_relative_bound<=lambda_control()$tail_tol))
  stopifnot(eta_numerical(1e-12,1)<1e-10)
})

test("capacity scaling, vectorization and high-resolution behavior", {
  r <- c(.1,1,10,100,1e6)
  a <- eta_numerical(r,1)
  relative(eta_numerical(r*3,3),a)
  equal(eta_numerical(c(2,20),c(2,2)),eta_numerical(c(1,10),1))
  stopifnot(all(diff(a)>0),all(a>eta_app(r,1)))
  relative(tail(a,1), .5*log(pi*exp(1)/6*1e6), 1e-3)
})

test("thermal boundaries and smooth small-temperature stability", {
  for (f in list(eta_numerical,eta_app,xi_eta_numerical,xi_eta_app,geometry_activation,xi_a))
    equal(f(0,1),0)
  r <- c(1e-16,1e-12,1e-8)
  relative(eta_app(r,1),pi*exp(1)/12*r,1e-8)
  relative(xi_eta_app(r,1),pi*exp(1)/24*r^2,1e-8)
  relative(xi_a(r,1),r^2/2,1e-8)
  relative(geometry_activation(r,1),r,2e-8)
})

test("direct quadrature agrees with known primitive and small-ratio limit", {
  relative(xi_eta_numerical(10,1)/10,1.3253096,2e-7)
  r <- c(1e-16,1e-12,1e-8)
  # Integral of u/3 * [1 + log(6/u)] in the three-bin limit.
  reference <- r^2/6 * (1.5+log(6/r))
  relative(xi_eta_numerical(r,1),reference,1e-8)
  relative(xi_eta_numerical(c(.1,10)*7,7),7*xi_eta_numerical(c(.1,10),1))
})

test("lookup agrees with direct integration across range and off-grid", {
  lookup <- make_eta_numerical_lookup(1e4)
  r <- c(0,1e-12,1e-5,.013,.15,1,3,10,100,3456,1e4)
  a <- xi_eta_numerical(r,1,lookup,details=TRUE)
  b <- xi_eta_numerical(r,1)
  equal(a$value[1],0)
  relative(a$value[-1],b[-1],1.1e-5)
  equal(a$method[2],"quadrature")
  stopifnot(all(a$abs_error>=0),tail(a$method,1)=="lookup")
  fails(xi_eta_numerical(1e5,1,lookup),"does not cover")
  fails(xi_eta_numerical(1,1,lookup,lambda_control(lookup_rel_tol=1e-7)),"looser")
  broken <- lookup; broken$ratio[2] <- broken$ratio[1]
  fails(xi_eta_numerical(1,1,broken),"Malformed")
})

test("numerical controls are validated and limits fail explicitly", {
  fails(lambda_control(rank_tol=0),"positive")
  fails(lambda_control(max_support=8.5),"integer")
  fails(eta_numerical(-1,1),"nonnegative")
  fails(eta_numerical(1,0),"positive")
  fails(eta_numerical(c(1,2),1:3),"lengths")
  fails(eta_numerical(Inf,1),"finite")
  fails(eta_numerical(1e10,1,lambda_control(max_support=8)),"max_support")
  fails(make_eta_numerical_lookup(10,grid_points=3,
    control=lambda_control(lookup_rel_tol=1e-10,lookup_max_points=17)),"lookup_max_points")
})
