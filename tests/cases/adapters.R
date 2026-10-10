# Simulation fixture plus deterministic coding and missingness perturbations.
d2 <- read.csv(fixture("ideal_costly_model3.csv"))
items2 <- paste0("y_item_", 1:6)
d2$group <- factor(rep(c("a", "b", "c"), length.out = nrow(d2)))
d2$binary <- factor(rep(c("a", "b"), length.out = nrow(d2)))
fit2 <- lm(y_mean ~ x1 + x2 + x3, data = d2)
formula2 <- lambda_model(d2, y_mean ~ x1 + x2 + x3, items2, 1/12)

test("measurement API preserves scale and case provenance", {
  a <- lambda_measure(d2, items2)
  b <- estimate_lambda_star_items(as.matrix(d2[items2]))
  equal(a$lambda_star, b$lambda_star)
  equal(a$temperature, b$temperature)
  equal(a$input_rows, seq_len(nrow(d2)))
  equal(lambda_measure(as.matrix(d2[items2]))$lambda_star, a$lambda_star)
  d <- d2
  d$y_item_1[3] <- NA
  fails(lambda_measure(d, items2), "Missing values")
  c <- lambda_measure(d, items2, na.action = "complete")
  equal(c$cases$excluded, 3L)
  equal(c$input_rows, setdiff(seq_len(nrow(d)), 3L))
  equal(c$row_ids, rownames(d)[-3])
  equal(c$lambda_star, lambda_measure(d[-3, ], items2)$lambda_star)
})

test("matrix, formula, lm and H1 paths agree", {
  x <- as.matrix(d2[c("x1", "x2", "x3")])
  a <- lambda_model(d2, x, items2, 1/12)
  b <- lambda_model(d2, fit2, items2, 1/12)
  c <- estimate_lambda_partition_items(as.matrix(d2[items2]), x)
  f <- thermodynamic_model_state(d2$y_mean, x, c$temperature_mean, 1/12)
  for (result in list(a, b, formula2)) {
    equal(result$lambda_h, c$lambda_h)
    equal(result$kappa, c$kappa)
    equal(result$parsimony$free_energy, f$free_energy)
    equal(result$parsimony$rss, sum(residuals(fit2)^2))
    equal(unname(result$basis$X), unname(center_columns(x)))
    equal(rownames(result$basis$X), rownames(d2))
    equal(result$model$rank, 3L)
    equal(result$model$columns, colnames(x))
  }
  equal(lambda_model(d2, ~ x1 + x2 + x3, items2, 1/12)$lambda_h, formula2$lambda_h)
})

test("two-group adapter agrees with pooled-variance t-test", {
  a <- lambda_model(d2, ~ binary, items2, 1/12)
  fit <- lm(y_mean ~ binary, d2)
  tt <- t.test(y_mean ~ binary, d2, var.equal = TRUE)
  model_f <- a$captured_observed / (a$residual_observed / a$df_perp)
  equal(unname(tt$statistic^2), model_f)
  equal(a$parsimony$rss, sum(residuals(fit)^2))
  equal(a$df_h, 1L)
  equal(a$model$intercept_removed, TRUE)
})

test("ANOVA, ANCOVA and interaction use the full encoded model", {
  for (f in list(y_mean ~ group, y_mean ~ group + x1,
                 y_mean ~ group * x1, y_mean ~ poly(x1, 2))) {
    fit <- lm(f, d2)
    a <- lambda_model(d2, f, items2, 1/12)
    b <- lambda_model(d2, fit, items2, 1/12)
    x <- model.matrix(fit)[, -1, drop = FALSE]
    c <- lambda_model(d2, x, items2, 1/12)
    equal(a$df_h, fit$rank - 1L)
    equal(a$parsimony$rss, sum(residuals(fit)^2))
    equal(a$parsimony$free_energy, b$parsimony$free_energy)
    equal(a$parsimony$free_energy, c$parsimony$free_energy)
    equal(a$model$term_assign, attr(model.matrix(fit), "assign")[-1])
  }
})

test("recoding preserves projection but can change geometry and F_H", {
  # Deterministic 3-level treatment versus orthogonal Helmert coding.
  a <- lambda_model(d2, ~ group, items2, 1/12,
                    contrasts = list(group = "contr.treatment"))
  b <- lambda_model(d2, ~ group, items2, 1/12,
                    contrasts = list(group = "contr.helmert"))
  equal(a$lambda_h, b$lambda_h)
  equal(a$parsimony$rss, b$parsimony$rss)
  stopifnot(abs(a$parsimony$log_volume - b$parsimony$log_volume) > .1,
            abs(a$parsimony$free_energy - b$parsimony$free_energy) > .01)
  equal(b$model$contrasts, list(group = "contr.helmert"))
  equal(b$model$factor_levels$group, levels(d2$group))
})

test("lm stored contrasts survive changes to global contrast options", {
  fit <- lm(y_mean ~ group + x1, d2, contrasts = list(group = "contr.sum"))
  local({
    old <- options(contrasts = c("contr.helmert", "contr.poly"))
    on.exit(options(old))
    a <- lambda_model(d2, fit, items2, 1/12)
    b <- lambda_model(d2, ~ group + x1, items2, 1/12,
                      contrasts = list(group = "contr.sum"))
    equal(a$basis$X, b$basis$X)
    equal(a$parsimony$free_energy, b$parsimony$free_energy)
    equal(a$model$contrasts, list(group = "contr.sum"))
  })
})

test("response mismatches and external formula variables are rejected", {
  d <- d2
  d$y_mean[1] <- d$y_mean[1] + .1
  fails(lambda_model(d, y_mean ~ x1, items2, 1/12), "response must equal")
  fails(lambda_model(d2, x1 ~ x2, items2, 1/12), "response must equal")
  fails(lambda_model(d2, cbind(y_mean, y_mean) ~ x1, items2, 1/12), "response must equal")
  external_x <- d2$x1
  fails(lambda_model(d2, ~ external_x, items2, 1/12), "must be columns of data")
  # Explicit tolerance is recorded, not silently applied to estimates.
  d$y_mean <- d2$y_mean + 1e-7
  fails(lambda_model(d, y_mean ~ x1, items2, 1/12, response_tol = 0), "response must equal")
  a <- lambda_model(d, y_mean ~ x1, items2, 1/12, response_tol = 1e-5)
  equal(unname(a$composite_response), unname(rowMeans(d[items2])))
  equal(a$settings$response_tol, 1e-5)
})

test("lm validates saved cases, design and response against common data", {
  fails(lambda_model(d2[nrow(d2):1, ], fit2, items2, 1/12), "cases/order")
  fails(lambda_model(d2, lm(y_mean ~ x1, d2, subset = x1 > 0), items2, 1/12), "cases/order")
  d <- d2
  d$x1[1] <- d$x1[1] + 1
  fails(lambda_model(d, fit2, items2, 1/12), "coding does not match")
  d <- d2
  d$y_mean <- d$y_mean + 1
  wrong_fit <- lm(y_mean ~ x1 + x2 + x3, d)
  fails(lambda_model(d2, wrong_fit, items2, 1/12), "coding does not match")
  fails(lambda_model(d2, lm(y_mean ~ x1, d2, model = FALSE), items2, 1/12), "saved model frame")
  fails(lambda_model(d2, lm(y_mean ~ 0 + x1, d2), items2, 1/12), "include an intercept")
  fails(lambda_model(d2, fit2, items2, 1/12, contrasts = list(group = "contr.sum")), "already specifies")
})

test("weights, offsets and other model families are rejected", {
  fails(lambda_model(d2, lm(y_mean ~ x1, d2, weights = rep(1, nrow(d2))), items2, 1/12), "Weighted")
  fails(lambda_model(d2, lm(y_mean ~ x1, d2, offset = x2), items2, 1/12), "Offsets")
  fails(lambda_model(d2, y_mean ~ x1 + offset(x2), items2, 1/12), "Offsets")
  fails(lambda_model(d2, glm(y_mean ~ x1, data = d2), items2, 1/12), "no GLM")
  fails(lambda_model(d2, lm(cbind(y_mean, x1) ~ x2, d2), items2, 1/12), "univariate")
  fails(lambda_model(d2, ~ x1 + (1 | group), items2, 1/12), "Random-effect")
})

test("rank deficiency, lost levels and invalid explicit row identity fail", {
  fails(lambda_model(d2, ~ x1 + I(2*x1), items2, 1/12), "rank deficient")
  fails(lambda_model(d2, ~ 0 + group, items2, 1/12), "rank deficient")
  fails(lambda_model(d2[d2$group != "c", ], ~ group, items2, 1/12), "constant or redundant")
  x <- as.matrix(d2[c("x1", "x2")])
  rownames(x) <- rownames(d2)
  fails(lambda_model(d2, x[-1, ], items2, 1/12), "same input rows")
  fails(lambda_model(d2, x[nrow(x):1, ], items2, 1/12), "row identities")
  fails(lambda_model(d2, cbind(intercept = 1, x), items2, 1/12), "constant")
  fails(lambda_model(d2, ~ group, items2, 1/12,
                     contrasts = list(typo = "contr.sum")), "factor predictors")
})

test("mean response and F_H temperature stay coupled under item scaling", {
  a <- lambda_model(d2, ~ x1 + x2 + x3, items2, 1/12, measurement_scale = "item")
  equal(a$lambda_h * length(items2), formula2$lambda_h)
  equal(a$lambda_star * length(items2), formula2$lambda_star)
  equal(a$kappa, formula2$kappa)
  equal(a$parsimony, formula2$parsimony)
  equal(a$parsimony$temperature, a$temperature_mean)
  equal(a$temperature, a$sigma2_item)
})

test("null model and temperature singularity remain explicit", {
  a <- lambda_model(d2, ~ 1, items2, 1/12)
  b <- lambda_model(d2, matrix(numeric(), nrow(d2), 0), items2, 1/12)
  equal(a$df_h, 0L)
  equal(a$parsimony$volume, 1)
  equal(a$parsimony$free_energy, sum((d2$y_mean - mean(d2$y_mean))^2))
  equal(a$parsimony$free_energy, b$parsimony$free_energy)
  d <- data.frame(i1 = 1:8, i2 = 1:8, x = 1:8)
  c <- lambda_compare(d, ~1, ~x, c("i1", "i2"), 1/12)
  equal(c$status, "singular_temperature")
  stopifnot(is.na(c$delta_F), is.na(c$model_a$lambda_star),
            is.na(c$model_b$lambda_h), is.na(c$model_b$kappa),
            is.na(c$model_b$parsimony$free_energy))
})

test("required data, indicators, resolution and settings are validated", {
  fails(lambda_model(d2, ~ x1, items2), "T_q.*missing")
  fails(lambda_model(d2, ~ x1, items2, 0), "positive")
  fails(lambda_model(d2, ~ x1, "x1", 1), "at least two")
  fails(lambda_measure(d2, c("group", "x1")), "numeric vectors")
  fails(lambda_model(d2, ~ x1, items2, 1, response_tol = -1), "nonnegative")
  fails(lambda_model(d2, list(formula = ~x1), items2, 1), "model must be")
  fails(lambda_model(d2, ~ x1, items2, 1, capacity = "smooth", eta_lookup = list()), "lookup")
  fails(lambda_model(d2, ~ x1, items2, 1, na.action = "omit"), "arg")
  fails(lambda_measure(as.matrix(d2[items2])[1, , drop = FALSE]), "at least two rows")
})

test("S3 summaries and printing expose estimates without full matrices", {
  a <- summary(formula2)
  equal(a$estimates$estimate, c(formula2$lambda_star, formula2$lambda_h,
                              formula2$lambda_perp, formula2$kappa))
  equal(a$parsimony$free_energy, formula2$parsimony$free_energy)
  stopifnot(is.null(a$basis))
  out <- capture.output(returned <- print(formula2))
  stopifnot(identical(returned, formula2), any(grepl("point estimates", out)))
  capture.output(print(lambda_measure(d2, items2)))
})
