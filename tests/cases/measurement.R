test("parallel measurements agree with independent pairwise identities", {
  items <- cbind(c(1, 3, 2, 6, 4), c(2, 2, 4, 5, 3), c(1, 4, 3, 5, 6))
  z <- scale(items, scale = FALSE)
  pairs <- combn(ncol(z), 2)
  signal <- mean(apply(pairs, 2, function(j) sum(z[, j[1]] * z[, j[2]])))
  bath <- mean(apply(pairs, 2, function(j) sum((z[, j[1]] - z[, j[2]])^2))) /
    (2 * (nrow(z) - 1))
  a <- estimate_lambda_star_items(items)
  equal(a$sigma2_item, bath)
  equal(a$signal_energy_raw, signal)
  equal(a$temperature, bath / ncol(items))
  equal(a$lambda_star_raw, signal / a$temperature)
})

test("two-replicate identity and default composite scale", {
  x <- cbind(1:8, c(3, 1, 5, 4, 7, 5, 9, 6))
  z <- scale(x, scale = FALSE)
  a <- estimate_lambda_star_items(x)
  equal(a$signal_energy_raw, sum(z[,1] * z[,2]))
  b <- estimate_lambda_star_items(x, scale = "item")
  equal(a$scale, "composite_mean")
  equal(a$temperature, b$temperature/2)
  equal(a$lambda_star_raw, b$lambda_star_raw * 2)
})

test("partition identity, scaling, and preserved observation identity", {
  data <- read.csv(fixture("ideal_costly_model3.csv"))
  items <- as.matrix(data[paste0("y_item_", 1:6)])
  x <- as.matrix(data[c("x1", "x2")])
  a <- estimate_lambda_partition_items(items, x, row_ids = data$id)
  b <- estimate_lambda_partition_items(items, x, scale = "item")
  equal(a$lambda_star_raw, a$lambda_h_raw + a$lambda_perp_raw)
  equal(a$lambda_star, a$lambda_h + a$lambda_perp)
  equal(a$lambda_star_raw, estimate_lambda_star_items(items)$lambda_star_raw)
  equal(a$kappa, b$kappa)
  equal(a$lambda_h_raw, 6 * b$lambda_h_raw)
  equal(a$row_ids, data$id)
  equal(a$input_rows, seq_len(nrow(data)))
  # Independent least-squares fit: no call to our projection implementation.
  fit <- lm(rowMeans(items) ~ x)
  equal(a$captured_observed, sum((fitted(fit)-mean(rowMeans(items)))^2))
  equal(a$residual_observed, sum(residuals(fit)^2))
})

test("zero bath is singular while zero signal has undefined kappa", {
  items <- cbind(1:10, 1:10, 1:10)
  a <- estimate_lambda_partition_items(items, matrix(1:10, ncol=1))
  equal(a$status, "singular_temperature")
  equal(a$temperature, 0)
  stopifnot(all(is.na(unlist(a[c("lambda_star", "lambda_star_raw", "lambda_h",
                                 "lambda_h_raw", "lambda_perp", "kappa", "rho_star")]))))
  b <- estimate_lambda_partition_items(cbind(1:10, -(1:10)), matrix(1:10, ncol=1))
  equal(b$status, "ok")
  equal(b$lambda_star, 0)
  stopifnot(b$lambda_star_raw < 0, is.na(b$kappa))
})

test("invalid measurements and row mismatches are rejected", {
  fails(estimate_lambda_star_items(matrix(1:5,ncol=1)), "two replicate")
  fails(estimate_lambda_star_items(matrix(c(1,NA,2,3),2)), "finite")
  fails(estimate_lambda_star_items(matrix(c(1,Inf,2,3),2)), "finite")
  fails(estimate_lambda_star_items(cbind(letters[1:5],1:5)), "numeric")
  fails(estimate_lambda_star_items(cbind(1:5,2:6), row_ids=rep("a",5)), "unique")
  a <- cbind(1:5,2:6); rownames(a) <- letters[1:5]
  b <- matrix(1:5,ncol=1); rownames(b) <- rev(letters[1:5])
  fails(estimate_lambda_partition_items(a,b), "row identities")
  fails(estimate_lambda_partition_items(a,b[1:4,,drop=FALSE]), "row counts")
})

test("tiny positive temperature is not clipped to zero", {
  a <- cbind(c(1,2,4,3),c(2,1,3,5))
  normal <- estimate_lambda_star_items(a)
  tiny <- estimate_lambda_star_items(a * 1e-100)
  equal(tiny$status, "ok")
  relative(tiny$temperature, normal$temperature*1e-200)
  relative(tiny$lambda_star_raw, normal$lambda_star_raw)
  fails(estimate_lambda_star_items(a*1e-200), "numerical range")
})
