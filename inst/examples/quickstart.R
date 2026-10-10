# Run after installation:
# source(system.file("examples", "quickstart.R", package = "lambdastar"))
library(lambdastar)
set.seed(17)
d <- data.frame(x = rnorm(120), z = rnorm(120),
                group = factor(rep(c("a", "b", "c"), 40)))
stable_signal <- d$x + .5*d$z + rnorm(120, sd=.5)
items <- paste0("item", 1:4)
d[items] <- replicate(4, stable_signal + rnorm(120))
d$score <- rowMeans(d[items])
T_q <- 1/12 # Example resolution assumption on the MEAN scale; choose substantively.

print(lambda_measure(d, items))
print(lambda_model(d, ~group, items, T_q)) # Whole ANOVA model, not a partial effect.
fit <- lambda_model(d, score~x+z, items, T_q,
  bootstrap=lambda_bootstrap_control(B=499, seed=17, keep_draws=TRUE))
print(fit)
print(fit$estimates) # Paper estimator by default; also used in every bootstrap draw.
print(lambda_model(d, ~x+z, items, T_q, estimator="population")$estimates)
print(confint(fit))
print(confint(fit, parm=c("lambda_star", "lambda_h", "kappa"), level=.90))

comparison <- lambda_compare(d, ~x, ~x+z, items, T_q)
print(comparison)
stopifnot(all.equal(comparison$model_b$parsimony$free_energy,fit$parsimony$free_energy))

# Actual application default: bootstrap=TRUE -> B=5000, 95%, seed=1.
# B=499 above is an explicit short-example override.
# Inspect fit$bootstrap$ci and fit$bootstrap$diagnostics for undefined draws.
# Optional conditional_quantiles=TRUE reports valid-only quantiles separately.

# Conditional model with group effects treated as nuisance.
print(lambda_model(d, ~x+z, items, T_q, nuisance=~group))
# Schools are resampled whole, within site; type7 matches the paper's percentiles.
d$school <- rep(seq_len(20), each=6)
d$site <- factor(rep(c("north","south"),each=60))
cluster_fit <- lambda_model(d, ~x+z, items, T_q, nuisance=~site,
  bootstrap=lambda_bootstrap_control(B=99,min_valid=50,cluster="school",
    resampling="cluster",strata="site",percentile_interpolation="type7"))
print(confint(cluster_fit,parm="kappa"))

# Continue with an explicit X_H/reference/basis workflow and interval diagnostics:
# system.file("doc", "MODEL_SPACE_WORKFLOW.md", package="lambdastar")
# source(system.file("examples", "intervals_and_targets.R", package="lambdastar"))
# intervals <- lambdastar_intervals_example()  # B=5000; B=199 for practice
