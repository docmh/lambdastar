test("comparison equals separate evaluations and obeys swap/self identities", {
  a <- lambda_model(d2, ~ x1, items2, 1/12)
  b <- lambda_model(d2, ~ x1 + x2 + x3, items2, 1/12)
  c <- lambda_compare(d2, ~ x1, ~ x1 + x2 + x3, items2, 1/12)
  # Provenance lists the UNION of needed variables for the comparison.
  equal(c$model_a$parsimony, a$parsimony)
  equal(c$model_b$parsimony, b$parsimony)
  equal(c$model_a$lambda_h, a$lambda_h)
  equal(c$model_b$lambda_h, b$lambda_h)
  equal(c$delta_F, b$parsimony$free_energy - a$parsimony$free_energy)
  reverse <- lambda_compare(d2, ~ x1 + x2 + x3, ~ x1, items2, 1/12)
  equal(reverse$delta_F, -c$delta_F)
  equal(lambda_compare(d2, ~x1, ~x1, items2, 1/12)$delta_F, 0)
  s <- summary(c)
  equal(s$models$F_H, c(a$parsimony$free_energy, b$parsimony$free_energy))
  capture.output(returned <- print(c))
  stopifnot(identical(returned, c))
})

test("common missingness selection precedes both data-dependent designs", {
  d <- d2
  d$y_item_1[2] <- NA
  d$x1[4] <- NA
  d$x3[7] <- NA
  d$unused <- NA_real_
  fails(lambda_compare(d, ~ poly(x1, 2), ~ x3, items2, 1/12), "Missing values")
  c <- lambda_compare(d, ~ poly(x1, 2), ~ x3, items2, 1/12, na.action = "complete")
  rows <- setdiff(seq_len(nrow(d)), c(2L, 4L, 7L))
  a <- lambda_model(d[rows, ], ~ poly(x1, 2), items2, 1/12)
  b <- lambda_model(d[rows, ], ~ x3, items2, 1/12)
  equal(c$cases$included, rows)
  equal(c$cases$excluded, c(2L, 4L, 7L))
  equal(c$model_a$row_ids, c$model_b$row_ids)
  equal(c$model_a$temperature, c$model_b$temperature)
  equal(c$model_a$basis$X, a$basis$X)
  equal(c$model_a$parsimony$free_energy, a$parsimony$free_energy)
  equal(c$model_b$parsimony$free_energy, b$parsimony$free_energy)
  # A separately evaluated model can have MORE cases than the joint comparison.
  stopifnot(lambda_model(d, ~ x3, items2, 1/12, na.action = "complete")$n > c$model_b$n)
})

test("matrix missingness joins the same common complete-case selection", {
  d <- d2
  d$y_item_1[2] <- NA
  x <- as.matrix(d[c("x1", "x2")])
  x[8, 2] <- NA
  c <- lambda_compare(d, x, ~x3, items2, 1/12, na.action = "complete")
  equal(c$cases$excluded, c(2L, 8L))
  equal(c$model_a$parsimony$free_energy,
    lambda_model(d[-c(2,8), ], x[-c(2,8), ], items2, 1/12)$parsimony$free_energy)
})

test("lm case omissions are accepted only when they match the common rows", {
  d <- d2
  d$x1[4] <- NA
  fit <- lm(y_mean ~ x1, d, na.action = na.omit)
  c <- lambda_compare(d, fit, ~x2, items2, 1/12, na.action = "complete")
  equal(c$cases$excluded, 4L)
  d$x2[7] <- NA
  fails(lambda_compare(d, fit, ~x2, items2, 1/12, na.action = "complete"), "cases/order")
})

test("two explicit models and a common matrix are required", {
  fails(lambda_compare(d2, ~x1, indicators = items2, T_q = 1/12), "model_b.*missing")
  fails(lambda_compare(model_a = ~x1, model_b = ~x2, indicators = items2, T_q = 1/12), "data.*missing")
  fails(lambda_compare(d2, NULL, ~x2, items2, 1/12), "model must be")
  fails(lambda_compare(list(d2, d2), ~x1, ~x2, items2, 1/12), "common data frame")
})

test("comparison retains each formula coding and common numerical settings", {
  ctrl <- lambda_control(rank_tol = 1e-9)
  c <- lambda_compare(d2, ~group, ~group, items2, 1/12,
    contrasts_a = list(group = "contr.treatment"),
    contrasts_b = list(group = "contr.helmert"), numerical_control = ctrl, capacity = "smooth")
  a <- lambda_model(d2, ~group, items2, 1/12, contrasts = list(group = "contr.treatment"),
    numerical_control = ctrl, capacity = "smooth")
  equal(c$model_a$parsimony, a$parsimony)
  equal(c$model_a$lambda_h, c$model_b$lambda_h)
  equal(c$model_b$settings$numerical_control, ctrl)
  stopifnot(abs(c$delta_F) > .01)
})

test("transform-created missing and nonfinite values never alter cases silently", {
  d <- d2
  d$bad <- d$x1
  d$bad[1] <- 0
  fails(lambda_model(d, ~ I(1/bad), items2, 1/12, na.action = "complete"), "finite")
  fails(suppressWarnings(lambda_compare(d, ~ log(-abs(x1)), ~x2, items2,
    1/12, na.action = "complete")), "transformation created missing")
  d$y_item_1 <- NA_real_
  fails(lambda_compare(d, ~x1, ~x2, items2, 1/12, na.action = "complete"), "Fewer than two")
})
