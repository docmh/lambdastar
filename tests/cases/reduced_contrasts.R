test("reduced factor contrasts agree across formula, matrix, lm and comparison bootstrap", {
  set.seed(721)
  d <- data.frame(g = factor(rep(letters[1:3], each = 30)))
  contrast <- matrix(c(-1, 0, 1), ncol = 1,
    dimnames = list(letters[1:3], "linear"))
  x <- contrast[as.integer(d$g), , drop = FALSE]
  rownames(x) <- rownames(d)
  signal <- 1.3*x[, 1] + .5*c(1, -2, 1)[as.integer(d$g)] + rnorm(90, sd = .4)
  items <- paste0("i", 1:6)
  d[items] <- replicate(6, signal + rnorm(90))
  d$y <- rowMeans(d[items])
  cfg <- lambda_bootstrap_control(B = 37, conf.level = .8, seed = 721,
    min_valid = 20, keep_draws = TRUE, keep_indices = TRUE)
  reference <- lambda_model(d, x, items, 1/12, bootstrap = cfg)
  f <- y ~ factor(g)
  fitted <- lm(f, d, contrasts = list("factor(g)" = contrast))
  a <- lambda_model(d, f, items, 1/12,
    contrasts = list("factor(g)" = contrast), bootstrap = cfg)
  b <- lambda_model(d, ~g, items, 1/12,
    contrasts = list(g = contrast), bootstrap = cfg)
  c <- lambda_model(d, fitted, items, 1/12, bootstrap = cfg)
  equal(a$model$contrasts[["factor(g)"]], contrast)
  equal(a$model$columns, "factor(g)linear")
  for (fit in list(a, b, c)) {
    equal(fit$basis$rank, 1L)
    equal(fit$estimates, reference$estimates)
    equal(fit$parsimony$free_energy, reference$parsimony$free_energy)
    equal(fit$bootstrap$indices, reference$bootstrap$indices)
    equal(fit$bootstrap$draws, reference$bootstrap$draws)
    equal(fit$bootstrap$ci, reference$bootstrap$ci)
    stopifnot(all(fit$bootstrap$ci$status == "ok"),
              nrow(fit$bootstrap$failures) == 0L)
  }
  cmp <- lambda_compare(d, f, fitted, items, 1/12,
    contrasts_a = list("factor(g)" = contrast), bootstrap = cfg)
  equal(cmp$delta_F, 0)
  equal(cmp$bootstrap$indices, reference$bootstrap$indices)
  equal(unname(cmp$bootstrap$draws[, 3:5]), unname(reference$bootstrap$draws[, 3:5]))
  equal(unname(cmp$bootstrap$draws[, 6:8]), unname(reference$bootstrap$draws[, 3:5]))
})

test("two reduced contrasts retain their directions with rebuilt transformations and interactions", {
  set.seed(722)
  d <- data.frame(g = factor(rep(letters[1:4], each = 24)), x = rnorm(96))
  contrast <- cbind(linear = c(-3, -1, 1, 3), quadratic = c(1, -1, -1, 1))
  rownames(contrast) <- letters[1:4]
  coding <- list("factor(g)" = contrast)
  f <- ~scale(x)*factor(g)
  x <- model.matrix(f, d, contrasts.arg = coding)
  signal <- drop(x %*% c(0, .8, .5, .4, .2, -.1)) + rnorm(96, sd = .3)
  items <- paste0("i", 1:6)
  d[items] <- replicate(6, signal + rnorm(96))
  cfg <- lambda_bootstrap_control(B = 19, conf.level = .8, seed = 722,
    min_valid = 10, keep_draws = TRUE, keep_indices = TRUE)
  fit <- lambda_model(d, f, items, 1/12, contrasts = coding, bootstrap = cfg)
  equal(fit$basis$rank, 5L)
  equal(fit$model$contrasts, coding)
  equal(fit$model$columns, colnames(x)[-1])
  stopifnot(all(fit$bootstrap$ci$status == "ok"), nrow(fit$bootstrap$failures) == 0L)
  for (b in seq_len(cfg$B)) {
    sample <- d[fit$bootstrap$indices[b, ], ]
    design <- model.matrix(f, sample, contrasts.arg = coding)
    y <- rowMeans(sample[items])
    k <- length(items)
    pairs <- combn(k, 2)
    temperature <- mean(apply(pairs, 2, function(j)
      var(sample[[items[j[1]]]] - sample[[items[j[2]]]]))) / (2*k)
    total <- max(0, (length(y)-1)*(var(y)-temperature)/temperature)
    ols <- lm.fit(design, y)
    equal(ols$rank, 6L)
    captured <- sum((ols$fitted.values-mean(y))^2)
    h <- min(total, max(0, captured/temperature-(ols$rank-1)))
    equal(unname(fit$bootstrap$draws[b, ]),
      c(temperature, total, h, total-h, h/total), tol = 1e-9)
  }
})
