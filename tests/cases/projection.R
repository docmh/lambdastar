test("known orthogonal energy partition is 9 plus 16", {
  x <- matrix(sin(1:20), ncol=1)
  q <- as.numeric(scale(x,scale=FALSE)); q <- q/sqrt(sum(q^2))
  v <- (1:20)^2; v <- v-mean(v); v <- v-sum(v*q)*q; v <- v/sqrt(sum(v^2))
  s <- 3*q+4*v
  basis <- hypothesis_basis(x)
  captured <- project_onto_basis(s,basis$Q)
  equal(sum(s^2),25)
  equal(sum(captured^2),9)
  equal(sum((s-captured)^2),16)
})

test("encoded geometry differs while recoded projection agrees", {
  z <- scale(cbind(sin(1:40),cos(1:40)),scale=FALSE)
  q <- qr.Q(qr(z))
  x <- cbind(q[,1], .8*q[,1]+.6*q[,2])
  a <- hypothesis_basis(x); b <- hypothesis_basis(q)
  equal(a$volume,.6)
  equal(b$volume,1)
  y <- sin((1:40)*2)+q[,1]
  equal(project_onto_basis(y,a$Q),project_onto_basis(y,b$Q))
  c <- hypothesis_basis(sweep(x,2,c(1e-100,-1e100),"*"))
  equal(c$log_volume,a$log_volume)
  equal(project_onto_basis(y,c$Q),project_onto_basis(y,a$Q))
})

test("null, one direction, rank loss, and bad projection", {
  equal(design_volume(matrix(numeric(),10,0)),1)
  equal(design_volume(matrix(1:10,ncol=1)),1)
  equal(design_volume(matrix(numeric(),10,0),logarithm=TRUE),0)
  fails(hypothesis_basis(cbind(1:10,2*(1:10))), "rank deficient")
  fails(hypothesis_basis(cbind(1:10,rep(1,10))), "constant")
  fails(hypothesis_basis(diag(3)), "rank deficient")
  fails(project_onto_basis(1:10,matrix(1,10,1)), "orthonormal")
  fails(project_onto_basis(1:10,matrix(0,9,0)), "dimensions")
})
