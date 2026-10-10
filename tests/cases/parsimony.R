test("single-model F_H agrees with independent lm and determinant", {
  d <- read.csv(fixture("ideal_costly_model3.csv"))
  items <- as.matrix(d[paste0("y_item_",1:6)])
  x <- as.matrix(d[c("x1","x2","x3")])
  measurement <- estimate_lambda_star_items(items)
  state <- thermodynamic_model_state(rowMeans(items),x,measurement$temperature,1/12)
  rss <- sum(residuals(lm(rowMeans(items)~x))^2)
  logv <- .5*as.numeric(determinant(cor(x),logarithm=TRUE)$modulus)
  expected <- rss + 3*xi_eta_numerical(measurement$temperature,1/12) -
    xi_a(measurement$temperature,1/12)*logv
  equal(state$status,"ok"); equal(state$eta_method,"numerical")
  equal(state$rss,rss); equal(state$log_volume,logv)
  equal(state$free_energy,expected)
  equal(state$free_energy,state$rss+state$dimension_cost+state$geometry_cost)
})

test("recoding changes F_H solely through geometry", {
  x <- qr.Q(qr(scale(cbind(sin(1:40),cos(1:40)),scale=FALSE)))
  recoded <- cbind(x[,1],.8*x[,1]+.6*x[,2])
  y <- sin((1:40)*3)+x[,1]
  a <- thermodynamic_model_state(y,x,1,.1)
  b <- thermodynamic_model_state(y,recoded,1,.1)
  equal(a$rss,b$rss); equal(a$dimension_cost,b$dimension_cost)
  equal(b$free_energy-a$free_energy,-xi_a(1,.1)*log(.6))
})

test("singular inference is distinct from the mathematical boundary", {
  x <- matrix(1:10,ncol=1); y <- sin(1:10)
  a <- thermodynamic_model_state(y,x,0,.1)
  equal(a$status,"singular_temperature"); stopifnot(is.na(a$free_energy))
  b <- thermodynamic_model_state(y,x,0,.1,boundary=TRUE)
  equal(b$status,"boundary_only"); equal(b$free_energy,b$rss)
  equal(b$dimension_cost,0); equal(b$geometry_cost,0)
  fails(thermodynamic_model_state(y,x,1,.1,boundary=TRUE),"reserved")
  fails(thermodynamic_model_state(y,x,-1,.1),"nonnegative")
})

test("null model, smooth opt-in, and response-scale covariance", {
  x <- matrix(numeric(),10,0); y <- sin(1:10)
  null <- thermodynamic_model_state(y,x,1,.1)
  equal(null$free_energy,sum((y-mean(y))^2))
  x <- matrix(1:10,ncol=1)
  a <- thermodynamic_model_state(y,x,1,.1)
  b <- thermodynamic_model_state(y*3,x,9,.9)
  relative(b$free_energy,9*a$free_energy)
  smooth <- thermodynamic_model_state(y,x,1,.1,eta_method="smooth")
  stopifnot(smooth$dimension_cost<a$dimension_cost)
  equal(smooth$eta_method,"smooth")
})
