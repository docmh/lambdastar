#' Estimator selection and ordinary percentile bootstrap
#'
#' estimator="paper" is the default. With nuisance design B (including an
#' intercept), let d=N-rank(B), p=rank(M_B X_H), y=M_B Y_mean and Q_H its
#' captured energy. On the mean scale, raw estimates are lambda_star=||y||^2/T-d,
#' lambda_h=Q_H/T-p and lambda_perp=RSS/T-(d-p). T is the parallel-item
#' disagreement divided by d*(K-1)*K. Item scale changes the denominator,
#' retaining the composite-mean bath correction. Estimates constrain total
#' signal to be nonnegative and captured signal to lie between zero and total.
#'
#' The selected points are in estimates and the top-level lambda/kappa fields.
#' rho_star is lambda_star/(lambda_star+d) for the selected estimator.
#' paper retains the manuscript estimates and raw values. inference and settings
#' identify the estimator independently of whether bootstrap is requested.
#'
#' estimator="population" explicitly selects the alternative random-design
#' estimator. Its targets at retained N are lambda_star=(N-1)*Var(S)/T,
#' lambda_h=(N-1)*Var(Pi_H S)/T and kappa=Var(Pi_H S)/Var(S), where Pi_H is the
#' linear population projection. Its raw components are V_S=var(Y_mean)-T_mean
#' and V_H=var(Y_mean)-RSS/(N-1-p). This correction assumes a linear conditional
#' mean and homogeneous residual variance. It is different from the manuscript
#' fixed-design correction when uncaptured stable signal remains. Set
#' S_hat=max(0,V_S), H_hat=min(S_hat,max(0,V_H)); estimates multiply these
#' components by (N-1)/T. A non-null model requires positive residual df.
#' population stores these components and estimates only when explicitly selected.
#' This option currently supports intercept adjustment and case resampling;
#' additional nuisance or cluster combinations are rejected explicitly.
#'
#' With population inference and fixed case strata, variance uses observed
#' proportions w_g, sample variances s_g^2 and means m_g:
#' `sum w_g [1-(1-w_g)/n_g] s_g^2 + sum w_g (m_g-sum w_g m_g)^2`.
#' Each such stratum needs two persons. Paper estimation uses projected energies;
#' specifying resampling strata does not change its point estimator.
#'
#' Ordinary case bootstrap draws N persons with replacement; all measurements,
#' predictors and nuisance variables travel together. Optional strata retain
#' observed case counts. Cluster bootstrap draws the observed number of whole
#' clusters within each stratum; sample size can vary. Cluster labels are scoped
#' within strata. Nuisance effects are refitted on each draw and d is recomputed.
#' Cluster names and strata are sorted for reproducible sampling; case-stratum
#' ordering follows first occurrence. Complete-case selection is shared by all
#' models and precedes sampling. Factor coding is held fixed.
#'
#' Every draw uses the selected point estimator, including its constraints.
#' Kappa is computed inside each draw. Direct equal-tail percentiles form
#' intervals for T and lambda/kappa. percentile_interpolation="boot" uses
#' boot::boot.ci(type="perc") and normal-quantile interpolation; "type7" uses
#' stats::quantile(type=7), reproducing the manuscript's interpolation. These
#' choices change interpolation only, never the estimator or resampling design.
#' A constant distribution returns identical endpoints. Quantile messages are
#' retained; unavailable quantiles produce NA with quantile_unavailable.
#' Changing confidence level reuses stored draws.
#'
#' Zero lambda_star is a valid value in its own distribution. With zero estimated
#' total signal, kappa is undefined and stored as NA with reason undefined_kappa.
#' At positive signal, kappa=0 and kappa=1 are valid observations. Original T=0
#' is singular and skips inference. A resampled T=0 remains a temperature value;
#' the dependent lambda/kappa values in that draw are undefined.
#'
#' Invalid draws are counted per parameter with reasons; none are replaced.
#' A nonfinite original estimate yields status undefined_estimate. Otherwise,
#' any undefined kappa draw yields undefined_draws; other invalid draws yield
#' failed_draws. The affected regular interval has NA limits. Other parameters
#' retain their own intervals. min_valid applies to quantile calculation.
#' Optional conditional_quantiles=TRUE returns quantiles of finite draws in
#' a separate conditional_quantiles table; these are conditional descriptive
#' quantiles, not regular confidence intervals. n_undefined counts zero-signal
#' kappa draws; n_failed counts all nonfinite parameter values.
#'
#' Stored draws contain actual parameter estimates including zeros and NA.
#' Optional original row indices and draw_reasons make each draw inspectable.
#' Factor coding is fixed; lost levels and rank loss have explicit diagnostics.
#' Data-dependent formula transforms are rebuilt; formulas must be deterministic.
#' The caller's RNG kind and state are restored. Intervals are pointwise;
#' F_H and delta_F remain point estimates. Interval availability, width and
#' coverage are distinct properties; kappa coverage has no defined target at
#' true zero signal.
#' A runnable regular/weak-signal example is available at
#' `system.file("examples", "intervals_and_targets.R", package="lambdastar")`.
#' Source that file and call lambdastar_intervals_example(); its retained synthetic
#' datasets have known population targets, while the estimator receives only
#' the observed measurement and predictor columns. The installed
#' `doc/MODEL_SPACE_WORKFLOW.md` explains why a fixed-profile projection, a
#' finite-N population target and the expectation of the paper statistic differ.
#'
#' Fixed-reference designs from [lambda_reference()] use the same paper statistic
#' on departures from a known origin, without an implicit intercept. With no
#' nuisance, d=N and ~1 retains the mean direction. Reference values are subset
#' with the original persons on every draw. [lambda_paired()] forms matching
#' after-before item differences before estimating T; its T_q is supplied on
#' the mean-difference scale. This path does not extend estimator="population".
#' The derivation is installed as `doc/REFERENCE_MEANS.md`.
#' @references
#'   \url{https://stat.ethz.ch/R-manual/R-devel/library/boot/html/boot.html}
#'   \url{https://stat.ethz.ch/R-manual/R-devel/library/boot/html/boot.ci.html}
#' @seealso [lambda_bootstrap_control()], [confint.lambda_measurement()]
#' @name lambda_inference
NULL
