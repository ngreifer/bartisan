# bartisan tasks

Implements Linero (2025), "Generalized Bayesian Additive Regression Trees Models: Beyond Conditional Conjugacy", with SoftBart-style soft decision rules from Linero and Yang (2018).

This file is organized by subject, not by session. Each entry states the problem, what was done, what it measured, and — where it applies — what was tried and abandoned. **Status** and **To Do** are the only sections that describe the present; everything under **Log** is a record of how it got here and should be read as history.

## Status

463 tests and 2312 expectations passing, 0 failures, 0 warnings and 0 skips when the `skip_on_cran()` fits are run locally with `NOT_CRAN=true`. `R CMD check` on the built tarball, with examples, tests and all ten vignette rebuilds, reports `Status: OK` with no notes and no warnings; inside it the suite skips 67 and passes 1560, which is the configuration a CRAN machine runs. `--as-cran` has not been made to complete on this machine, for a reason that looks environmental rather than packaged; see the pre-submission report at the end of this file.

**What exists.** A C++ engine (`utils`, `slice`, `hypers`, `family`, `polyagamma`, `node`, `mcmc`, `model`) and an R interface following `glm()`: `bartisan()`, `bartisan_control()`, `predict()`, `print()`, `summary()`, family normalization, parallel chains with convergence diagnostics, and `custom_family()` for a likelihood written in R. Families: Gaussian, binomial (logit/probit/cloglog/any link from R), Poisson, negative binomial, gamma, ordinal (logit/probit/cloglog), multinomial (symmetric or reference-coded), multinomial probit with a drawn latent covariance, three AFT variants, location-scale, zero-inflated Poisson and negative binomial, ordered beta, and a Dirichlet process mixture for the error distribution. Missing predictors handled natively by MIA and kept by default. Data augmentations, on by default, for the binomial, ordinal, multinomial and zero-inflated families, and for the negative binomial under hard rules. `marginaleffects` support, so counterfactual estimands come with posterior intervals. Group-level random intercepts through lme4's `(1 | group)` notation, on every additive predictor. Posterior predictive draws for every family that has a sampler, and with them the interfaces to `loo`, `bayesplot`, `performance` and `posterior`. `bartisan_control()` organized into modeling decisions, advanced settings and validation toggles, with a per-forest `num_trees` vector, one `gate` argument covering hard and soft rules, and a `sparsity` argument standing in for the four DART hyperparameters. Documentation, `README.Rmd`, ten vignettes with a shared `references.bib`, and `_dev/benchmark.Rmd`. No `NEWS.md`: nothing has been released, so there is no previous version for a user to have seen.

Benchmark, Friedman function, n = 1000, p = 10, 50 trees, 1000 warmup plus 1000 saved, best of 3, scored against the true regression function on a held-out thousand. Reproducible with `_dev/benchmark.Rmd`. **The table below is not comparable with versions of it from before 2026-09-02**: until then bartisan's Gaussian rows omitted `family = gaussian()` and so timed a `dpm()` fit, and the competing packages were passed test data inside the timed call while bartisan's `predict()` ran outside it. Both are fixed; see the Log entry.

| Task | Package and call | Seconds | ESS | RMSE |
|---|---|---|---|---|
| gaussian | dbarts bart() | 0.281 | 24.0 | 0.225 |
|  | bartisan hard rules | 0.905 | 21.8 | 0.193 |
|  | stochtree bart(num_gfr = 0) | 1.189 | 19.5 | 0.245 |
|  | stochtree bart(num_gfr = 5) | 1.248 | 22.5 | 0.228 |
|  | BART wbart() | 1.932 | 23.1 | 0.231 |
|  | bartMachine bartMachine() | 2.082 | — | 0.234 |
|  | bartisan soft, smoothstep gate (default) | 2.887 | 50.1 | 0.130 |
|  | bartisan soft, smootherstep gate | 2.979 | 56.5 | 0.126 |
|  | bartisan soft, logistic gate | 4.206 | 54.4 | 0.130 |
| probit | dbarts bart() | 0.360 | 51.4 | 0.113 |
|  | bartisan hard, augment = TRUE | 1.094 | 35.8 | 0.107 |
|  | stochtree bart(binary/probit) | 1.876 | 34.0 | 0.122 |
|  | BART pbart() | 2.055 | 39.4 | 0.117 |
|  | bartisan soft, augment = TRUE | 3.104 | 64.2 | 0.086 |
|  | bartisan soft, augment = FALSE | 42.082 | 72.9 | 0.085 |
| logit | bartisan soft, augment = TRUE | 3.973 | 178.8 | 0.082 |
|  | bartisan soft, augment = FALSE | 25.036 | 253.1 | 0.088 |
|  | BART lbart() | 48.884 | 68.6 | 0.108 |
| ordinal | bartisan hard, logit, augment = TRUE | 1.872 | 44.2 | — |
|  | bartisan hard, probit, augment = TRUE | 2.405 | 43.6 | — |
|  | bartisan soft, logit, augment = TRUE | 3.919 | 51.7 | — |
|  | bartisan soft, probit, augment = TRUE | 5.152 | 99.7 | — |
|  | bartisan hard, cloglog, augment = TRUE | 6.421 | 35.9 | — |
|  | stochtree bart(ordinal/cloglog) | 8.102 | 27.0 | — |
|  | bartisan hard, logit, augment = FALSE | 32.143 | 46.8 | — |
|  | bartisan hard, cloglog, augment = FALSE | 38.243 | 32.9 | — |
|  | bartisan hard, probit, augment = FALSE | 51.737 | 55.2 | — |
| gamma (log link) | bartisan hard rules | 10.486 | 26.1 | 0.192 |
|  | bartisan hard, no shortcut | 17.996 | 29.0 | 0.209 |
| log-logistic AFT | bartisan hard rules | 2.094 | 54.9 | 0.305 |
| negative binomial | bartisan hard, augment = TRUE | 10.012 | 28.3 | 0.290 |
|  | bartisan hard, augment = FALSE | 17.774 | 45.8 | 0.303 |
| poisson | bartisan hard rules | 5.112 | 21.0 | 0.189 |
|  | bartisan hard, no shortcut | 9.864 | 22.5 | 0.196 |

Four things this says.

**Hard rules are within 1.8x of dbarts** on the two tasks dbarts supports, at the same mixing and slightly better accuracy: 0.426 s against 0.241 s on the Gaussian task and 0.482 s against 0.272 s on probit, best of six runs each. That ratio has come down 3.1x → 2.4x → 2.03x → 1.85x → 1.77x as each of the entries below landed.

**Note on measurement.** The table above is `_dev/benchmark.Rmd` at two replicates, which is noisy at the ten to thirty percent level — the Gaussian hard-rule cell has read 0.441 and 0.593 on consecutive runs of the same build, with a best-of-five standalone measurement of 0.426 either side of it. Any claim about a ratio near two needs more replicates than the document's default, and the two ratios quoted above are from a longer run for that reason.

**bartisan is now the faster of the two on the ordinal complementary log-log model**, which is the one task stochtree supports and dbarts does not: 3.03 s against 4.01 s.

**bartisan is faster and more accurate than every other package here** on both tasks — stochtree, bartMachine and BART included — once hard rules are used.

**Soft rules are the accuracy argument, not a tax.** They cost three to five times the hard-rule time and cut held-out error by 35–40%, which makes the default configuration the most accurate fit in the table, dbarts included. With a bounded gate the cost falls to 3.2x.

**The ordinal RMSE column is not comparable across links.** The three links put the additive predictor on three different scales, so those numbers compare scales rather than fits; seconds and ESS are comparable. It became *readable* for the probit rows once the identification changed, because the reported predictor is now centered and so is the generating one.

## To Do

Split by whether it stands between the package and a submission. `_dev/SHIP.md`
holds the release assessment, and the pre-submission report at the end of this
file is the current read on what is left.

Three items came off this list rather than being carried forward, and the
reasoning for each is in the Log rather than here: the grow-from-root warm start
and the soft random tree features, both struck because the burn-in transient
they would shorten is 34 to 70 sweeps against a default warmup of 200; and the
joint update for correlated nuisance parameters, which `tweedie()` made testable
and which turned out not to be needed. Two more came off as already done.
`predict(type = "density")` no longer returns `NaN` silently, since it warns and
names the draws responsible, which is what `_dev/SHIP.md` records and this list
had not caught up with. And the rename: what collided with the archived CRAN
package `genBart` was the *old* name `genbart`, the rename to `bartisan` landed
in `cdc3278`, and the item survived the rename it describes. It was carried into
a pre-submission report as a blocker before anyone checked the two names against
each other, which is the sort of thing to check rather than read.

### Before submission

- [ ] **Give the package a release version.** `--as-cran` notes that `0.0.0.9000` "contains large components", which is its way of saying a development version should not be submitted. `0.1.0` is the obvious choice.

- [ ] `vignette("bartisan")` does not cover the bounded gates or either ordinal augmentation; `vignette("families")` covers the augmentations but not the gates. The first runs on a reduced chain (20 trees, 300 draws, n = 400) and builds in about 85 seconds.

### After submission

Ordered by expected value. None of these is load-bearing for a workflow the
package claims to support, which is what puts them here.

- [ ] A **formal variable-selection test** rather than a threshold on `prop_used`, which `dbarts`, `SoftBart` and `bartMachine` all have and this does not. `variable_importance(draws = TRUE)` now returns the per-draw counts, so a caller can compute the posterior probability that one predictor takes more rules than another and does not have to read a difference off the table; what is missing is a calibrated null to test against, which is what the permutation approach in those packages supplies. The partial-dependence half of this item is closed: *marginaleffects*' `plot_predictions()` draws it, and `variable_importance(plot = TRUE)` covers the usage ranking.
- [ ] **Correlated random effects across additive predictors**, which is the one part of the random-effects feature that is not there. The obstacle is the absence of mixed second derivatives in the `Family` interface, and the alternative needs a prior mean threaded through 27 places; see the assessment.
- [ ] A Bayesian-bootstrap dispersion draw for `Gamma()` and `negbin()`, from Pearson residuals rather than the assumed likelihood. Cheapest available improvement to interval calibration; see the quasi-likelihood entry.
- [ ] A joint tridiagonal update for the ordinal cutpoints, which is what makes inference on the thresholds usable when there are many of them. The obstacle is the ordering constraint, not the algebra.
- [ ] Consider the `draw_prior` move from SoftBart, which proposes a whole fresh tree and helps escape local modes. It needs an `L`-dimensional Laplace proposal for the new leaves, so it is real work, not a port.
- [ ] `custom_family()` has no posterior predictive distribution, since a log density supplies no way to draw from it. An optional `rng` argument alongside the density would give it one, and would make `simulate()`, `pp_check()` and `r2()` work for a user-written likelihood. Small work; the design question is whether to also ask for a mean.
- [ ] Missing data, further work: nothing forces the three missing-value rules to be equally likely, and a variable with a handful of missing values probably does not want a third of its rules spent on splitting by missingness. A prior weight on the third rule is a one-line change and an open question.
- [ ] A `quasi()` family parameterized by link, variance function and dispersion update rule. Needs a documented weakening of the exactness claim.
- [ ] A lighter-tailed prior on the leaf scale, or an upper bound, would remove the separation pathology at the cost of changing the default prior. Not done unilaterally; the warning is the interim measure.
- [ ] **Relative survival on top of `ph()`**, per Basak et al. (2024): the excess-hazard model needs one extra Bernoulli draw per sweep, `d_i ~ Bernoulli(lambda_E / (lambda_E + lambda_P))`, with the population hazard supplied as one number per subject from a life table. Cheap now that `ph()` exists -- a nuisance draw and a data column. Narrow audience (cancer registries), so worth doing only on request.

## Speeding up the survival models

Four changes, all of them from the same observation: the leaf-level target's *shape* is what costs, not the non-conjugacy. Measured against the direct likelihood on the same data and seed.

**The rate generalization.** `Family::exp_sign()` returned +1 or -1, so the sampler recognized `exp(eta)` and `exp(-eta)` and nothing else. It is now `exp_rate()`, returning the rate in `c + a mu + b exp(r mu)`. The arithmetic generalizes cleanly once written out: `b exp(r ref) = curve / r^2` and `a = slope - curve / r`, which at r = ±1 reduces to the sign arithmetic it replaced, since there `1/r == r` and `r^2 == 1`. `exponential_mode()`'s Newton step needed the same treatment (`score = a + r ex - mu prec`, `info = prec - r^2 ex`). Verified as a bitwise no-op for the Poisson and the gamma, whose rates are ±1. The member is now spelled `rate_` rather than `sign_`, which is what it is.

**`location_scale()`'s scale forest.** Its log density is `const - eta1 - (y - eta0)^2 exp(-2 eta1) / 2`, which is the exponential form at rate -2 -- previously unreachable, so the forest ran on the general path. Declaring `TARGET_EXP_DOWN` with `exp_rate() = -2` required one other change: its `info_unit(h = 1)` returned the *expected* information, a constant 2, because that is quieter, and the exponential form reads its coefficients off the curvature the target actually has. Switching to the observed `2 r^2` makes the extraction exact. The two paths then agree to 2e-15 on the mean predictor and 5e-6 on the log-scale one (the Newton tolerance), and the fit went from **7.7 seconds to 3.6** under hard rules at 50 trees, with the same accuracy. Soft rules cannot use the form, so `location_scale()` at the default gate is unchanged.

**`weibull_aft()` needed no augmentation at all.** The plan had been to impute censored times from a truncated Gumbel. Writing out the likelihood showed the imputation was unnecessary:

    delta (r - log sigma) - exp(r) = c - (delta / sigma) eta - exp(y / sigma) exp(-eta / sigma)

with `r = (y - eta) / sigma`. That is already the exponential form, at rate `-1 / sigma`, and *censoring does not break it*: delta = 0 drops the linear term and leaves the shape intact. The rate is the same for every observation, which is the form's requirement, and `info_unit` was already the observed curvature. So it is two overrides and no new sampler: 6.7 seconds to **4.2** under hard rules, reproducing the general path to 1.3e-15. Under soft rules it stays general, and no augmentation would change that, because the truncated-Gumbel route lands in the exponential form too -- which soft rules cannot use. This makes the Weibull the slowest of the three AFT families at the default gate.

**The two augmentations that do pay.** The log-normal and log-logistic have no exponential form, and for them the imputation is the whole point. Right-censoring is what makes the direct likelihood expensive: a failure contributes a density in the predictor and a censored observation a survival function, and the two have different shapes, so the target has no exploitable form and every trial leaf value costs its own pass. Imputing the failure time above its censoring time replaces the survival term with a density and makes every contribution quadratic. Unlike the exponential form, **this route survives soft rules**, which is why it is the one that helps at the default gate.

The log-normal imputes from a normal truncated below; the log-logistic imputes from a logistic truncated below (by inverting the CDF through the upper tail probability, so that a censoring time far above the predictor does not collapse onto the endpoint) and then draws a Pólya-Gamma precision, the same device the ordinal logit uses. `sigma` is drawn from the *observed*-data likelihood with the imputations integrated out, and only then are they redrawn -- partially collapsed, and free here, because it is the same sum over observations either way.

| Family | Rules | Speed | ESS ratio | ESS per second |
|---|---|---|---|---|
| `lognormal_aft()` | soft | 11.2 to 19.7x | 0.83 to 1.27x | **14 to 16x** |
| `lognormal_aft()` | hard | 19.7 to 29.3x | 0.74 to 1.01x | **20 to 22x** |
| `loglogistic_aft()` | soft | 8.4 to 9.3x | 0.58 to 0.84x | **4.9 to 7.7x** |
| `loglogistic_aft()` | hard | 11.3 to 12.3x | 0.87 to 0.89x | **9.8 to 11x** |

Both beat the estimates in the old To Do entries (10x and 7x on raw speed), and the log-normal's hard-rule number puts it beside `ordinal("probit")` as one of the two largest gains in the package.

**How the augmentations were checked.** Posterior agreement, not just speed, and on the right scale. Comparing posterior means directly is misleading: two different samplers agree only to Monte Carlo error, so the question is whether the between-sampler gap is bigger than that error. Running two seeds of each sampler and measuring every gap in units of the pointwise posterior standard deviation, on the log-logistic at 50% censoring:

| Comparison | Mean gap | Max gap |
|---|---|---|
| augmented against augmented | 0.13 | 0.54 |
| direct against direct | 0.15 | 0.48 |
| **augmented against direct** | **0.09** | **0.30** |

The between-sampler gap is *smaller* than each sampler's own Monte Carlo error. Across censoring from 10% to 70%, `sigma`'s posterior mean agreed to within 0.002 with a posterior standard deviation of 0.02 to 0.03, and its posterior standard deviation to within 0.002.

## Is `Gamma()` worth keeping?

Asked directly, and the answer is a narrow yes. Measured on 800 training and 800 test observations, 50 trees, 500 draws after 500 warmup, four replicates, over four shapes of positive error around the same mean function, with Jacobian-corrected log scores.

| Errors | `Gamma("log")` | `Gamma("inverse")` | `dpm()` | `gaussian()` on `log(y)` | `ordinal("probit")`, 25 bins |
|----|----|----|----|----|----|
| gamma, RMSE | 0.401 | 0.535 | 0.856 | **0.388** | 0.414 |
| gamma, log score | **-2024** | -2028 | -2102 | -2041 | -- |
| lognormal, RMSE | 0.234 | 0.265 | 0.400 | **0.225** | 0.232 |
| lognormal, log score | -1525 | -1527 | -1600 | **-1521** | -- |
| heavy tail, RMSE | 2.454 | 3.257 | 2.215 | **2.110** | 2.269 |
| heavy tail, log score | -2359 | *NaN* | **-1947** | -2438 | -- |
| mixed, RMSE | 1.606 | 1.624 | 1.808 | **1.419** | 1.610 |
| mixed, log score | -1935 | -1938 | **-1913** | -1987 | -- |
| seconds | 6.7 | 11.1 | 1.3 | **1.1** | 1.8 |

**`Gamma("log")` is the only family that gets the best predictive density when the gamma is the truth.** That is what a correctly specified family should do, and it is the reason to keep it. But it never has the best RMSE for the conditional mean, and it is the slowest of the five by four to six times, because its exponential form pays only under hard rules and the default gate is soft.

**`Gamma("inverse")` is worse on every measure and can return NaN.** The link message was already justified on speed; this is the stronger reason. The additive predictor is unconstrained and the inverse link needs it positive. Over eight replicates on the heavy-tailed setting, five had draws that went non-positive -- 0.001% to 0.04% of them -- and each produced NaN densities for two to five of 800 test points. Recorded as a To Do, because `predict()` says nothing about it.

**The gamma's tails are its weak point.** A first pass reported the heavy-tailed log score as -Inf, which was one replicate poisoning a mean rather than a uniform failure. Per replicate over eight: median -2359, worst -557785, with two replicates containing a point whose density underflowed. `dpm()` had a median of -1947 and a worst of -15602 and never underflowed. A fixed shape cannot supply density in a tail it does not have.

**`ordinal("probit")` on 25 bins is the most consistent of the five for the conditional mean** -- within 3% of the best in all four settings, including the heavy-tailed one where the gamma is 8% behind it, at under a third of the gamma's time, with no assumption about the error's shape. It cannot give a density on the original scale, which is the price.

**One qualification to an earlier recommendation.** `dpm()` is roughly twice the RMSE of everything else for the conditional mean under gamma or lognormal errors. The advice to prefer `dpm()` over `gaussian()` was established on symmetric errors and does not carry over to the raw scale of a skewed positive outcome, where its flexible error absorbs signal that belongs in the mean. On `log(y)` it is fine. The vignette now says so.

## The gamma link, decided rather than reported

The link is now replaced rather than composed: any `Gamma()` whose link is not `log` is fitted on `log` with a message, and only an explicit `Gamma("log")` (or the string `"Gamma"`, which is this package's own spelling) is silent. `stats::Gamma()` is untouched, so attaching the package still cannot change what `glm()` does.

The reason base R defaults to the inverse link is that it is the **canonical** link for the gamma. In exponential-family form the gamma's natural parameter is `-1/mu`, and `glm()` follows McCullagh and Nelder in defaulting every family to its canonical link -- Gaussian to identity, binomial to logit, Poisson to log, gamma to inverse, inverse Gaussian to `1/mu^2`. Canonical links earn that status because the observed and expected information coincide there, which makes IRLS exactly Newton-Raphson, and because `X'y` is then sufficient. It is a theoretical convention, not a practical recommendation, and the gamma is the case where the two come apart: the canonical link does not keep the mean positive, which is a known wart in ordinary GLM practice too and the reason most applied gamma regression uses the log link anyway.

That wart is worse here than in a GLM, which is what settled it. A GLM's IRLS can step out of the parameter space and be caught; a sampler's additive predictor is unconstrained by construction and *will* visit negative means. Measured: over eight replicates on heavy-tailed data, five had draws where the predictor went non-positive -- 0.001% to 0.04% of them -- and each produced NaN densities for two to five of 800 test points. Plus 11.1 seconds against 6.7, and worse RMSE in all four settings.

The cost of the change: composed links are no longer reachable for the gamma, so `Gamma("sqrt")` -- whose inverse at least stays non-negative -- is overruled along with the rest. `Gamma` was removed from `native_links` to say so. A caller who genuinely wants one now needs `custom_family()`. Judged worth it: the link was a documented feature almost nobody wants, and leaving it available meant leaving the inverse link available too.

## `Beta()`, and what it is actually worth

Added because `ordbeta()` without it was an odd gap: the package could model a proportion that reaches its bounds but not one that does not. It is the interior of `OrdBetaFamily` with the endpoint machinery removed -- about eighty lines of C++ -- and it is capitalized for the reason `stats::Gamma()` is, since `base::beta()` is the beta function.

Measured against the alternatives, two replicates each, 50 trees, 300 draws after 300 warmup, RMSE for the conditional mean on held-out data:

| | `Beta()` | `ordbeta()` | `gaussian()` on `logit(y)` |
|---|---|---|---|
| n = 200, interior only | **0.0223** | 0.0251 | 0.0224 |
| n = 500, interior only | 0.0206 | **0.0204** | 0.0257 |
| n = 500, 1% at a boundary | 0.0222 (clamped) | **0.0210** | 0.0299 |

The honest reading, which is weaker than the case for adding it: **`Beta()` is 11% better than `ordbeta()` at n = 200 and a tie at n = 500.** The theory predicts exactly that -- `ordbeta()`'s two cutpoints have nothing to identify them without boundary observations, and unidentified nuisance parameters cost more when there is less data -- but the cost is small and vanishes with sample size. It is also only 13% faster.

So the argument for the family is not accuracy. It is that `ordbeta()` on interior data leaves two parameters wandering: the fitted cutpoints reached magnitudes of 18 to 19 against slice bounds of ±30, which is a diagnostic that looks like a pathology and is really just an unidentified parameter doing what unidentified parameters do. `Beta()` reports one precision, which means something. Choose between the two on whether the response *can* reach a boundary.

The row that does earn its keep unambiguously is the third: with boundary observations, `ordbeta()` beats a clamped `Beta()` by 5% and `gaussian()` on the logit by 30%, with coverage of 0.98 against 0.85. And `gaussian()` on `logit(y)` is worth knowing about: 100 times faster, competitive at n = 200, clearly worse by n = 500 (25% worse RMSE, coverage 0.89 against 0.96).

## The families vignette, cut in half

Rewritten from 9,414 words to 4,917 -- 48% shorter -- on the brief that the document exists to orient a user and support a practical choice, not to survey the method. What came out, as a rule: anything a user cannot act on. The augmentation ESS tables (they live in `?bartisan_control`), the `dpm()` centering derivation, the normal-inverse-chi-square baseline and its `nu`/`q`/`k_s` calibration, the Murray parameterization's leaf-prior scaling, the trace constraint and the inverse-Wishart prior, the multinomial-probit correlation sweep table, the ordinal binning table, the 200-observation continuous table, the `polr()` chart-matching demonstration, and a mea culpa about an earlier version of the documentation that belonged in this file rather than in a vignette.

What went in: the four-way continuous comparison with `location_scale()` added, the positive-outcome comparison re-run without `Gamma("inverse")`, `ordinal()` as a fourth count family, and the corrected defaults. Citations fell from 14 to 4, which is the honest consequence of cutting the method exposition; `references.bib` keeps all 19 verified entries, since pandoc emits only the cited ones and the rest are work worth keeping.

**Adding `location_scale()` to the continuous comparison is what made that table worth having.** With four families over five error shapes it now says one thing per row instead of needing two tables and a paragraph of hedging:

| Errors | `gaussian()` | `dpm()` | `location_scale()` | `ordinal("probit")` |
|---|---|---|---|---|
| normal | 0.135 / -1439 | 0.140 / -1440 | 0.137 / -1442 | **0.129** |
| $t_3$ | 0.094 / -1412 | **0.080 / -1253** | 0.101 / -1400 | 0.105 |
| skewed | 0.078 / -1106 | **0.053 / -956** | 0.078 / -1115 | 0.089 |
| bimodal | 0.149 / -1657 | **0.060 / -1240** | 0.136 / -1656 | 0.166 |
| heteroskedastic | 0.168 / -1684 | 0.165 / -1653 | **0.136 / -1541** | 0.147 |
| seconds | **1.6** | 2.1 | 17.4 | 2.7 |

`dpm()` for a badly *shaped* error, `location_scale()` for a *varying* one, everything level when the errors are normal, and `ordinal()` never far off in either direction. `location_scale()` beating `dpm()` by 112 log points on the heteroskedastic row is a cleaner statement of the division than the two separate tables managed.

**Three claims in the draft did not survive checking against the numbers**, all of them inherited from the longer version:

- "bimodality worth 417 log points **with a third the RMSE**" -- the ratio is 0.060/0.149 = 0.40, so 40%, not a third.
- "`c(50, 10)` was **three times faster** than `c(50, 50)`" -- 14.5 against 5.9 seconds is 2.5x.
- "`ordinal("probit")` is **within 3% of the best RMSE in every row**" of the positive-outcome table -- it is up to 13.5% off the best column. What is true is that it tracks `Gamma("log")` within about 3% in every row and beats it by 8% on heavy tails. The original claim was wrong when it was first written, not broken by the re-run.

**The heavy-tail log score was dropped from the positive-outcome table rather than reported.** The 4-replicate re-run put `Gamma("log")` ahead (-2250 against `dpm()`'s -2392); the earlier 8-replicate run put it behind (-2359 against -1947). Both are medians, and they disagree on the ordering, because the score is dominated by the single worst test point. The stable statement is the range: over eight replicates `Gamma("log")` ran from -1745 to -557785 and `dpm()` from -1661 to -15602. Medians close, tail risk not, and the gamma carries much more of it. That is what the vignette now says, in place of a number that would have looked authoritative and been noise.

## Making the beta families fast

They were the slowest families in the package -- `Beta()` at 18.4 seconds and `ordbeta()` at 20.4 where the slowest *other* general-target family, `binomial("cloglog")` with the augmentation off, took 3.2 (500 observations, 50 trees, 300 draws). Now 4.3 and 4.9. On a longer run (3 seeds, 500 warmup and 1000 draws) `Beta()` went from 91.5 seconds to 20.0 and `ordbeta()` from 104.1 to 24.9, so **4.6x and 4.2x**, with effective samples per second up 5.0x and 3.7x.

**First, what the cost was not.** Instrumenting the Fisher-scoring loop settled it: `Beta()` used 38,089 score passes and `binomial("cloglog")` 37,697 -- within 1% -- for 6x the time. So the iteration count was never the problem, and neither was the general target as such. It was the four `digamma`/`trigamma` calls per observation per pass.

**The structural observation that fixes it.** The two special functions appear only in the combinations

    psi(mu phi) - psi((1 - mu) phi)      and      psi'(mu phi) + psi'((1 - mu) phi),

and the two shapes always sum to `phi`. `phi` moves only in `update_aux`, once per sweep, so *within a sweep both combinations are functions of the single scalar `mu`*. They are now tabulated on a grid of 2049 points in the additive predictor over (-8, 8), rebuilt once per sweep and linearly interpolated, with the exact functions used outside the grid. Four special-function calls become two loads and a multiply. The per-sweep rebuild costs 2049 evaluations against the millions it replaces.

**Why an approximate derivative is legitimate here**, which is the part worth stating: the score and information shape the *Laplace proposal only*. `logdens_unit` stays exact, and the Metropolis step corrects. What the sampler requires of the fit is that it be a **deterministic** function of the current state, so the birth and death moves rebuild the same proposal -- the comment above `SCORE_TOL` says exactly this -- and a fixed table satisfies it. Interpolation error costs acceptance rate, not correctness.

Verified as such. Two seeds of each version, 1000 warmup and 2500 draws, differences measured in units of the pointwise posterior standard deviation of the fitted mean:

| Comparison | `Beta()` | `ordbeta()` |
|---|---|---|
| exact against exact (Monte Carlo error alone) | 0.06 mean, 0.32 max | 0.11 mean, 0.43 max |
| table against table (Monte Carlo error alone) | 0.09 mean, 0.37 max | 0.09 mean, 0.27 max |
| **exact against table** | **0.05 mean, 0.21 max** | **0.03 mean, 0.16 max** |

The between-version gap is smaller than either version's own Monte Carlo error, and `phi`'s posterior mean agrees to 0.01 and 0.02 posterior standard deviations. The tabulated score also matches a central difference of the exact log density to 9e-7 relative, which is what `test-derivatives.R` now checks.

**A crude approximation was tried first and is worse.** The leading term of the expansion in `1/phi` is elegant -- `psi(a) - psi(b)` goes to `logit(mu)`, which is the predictor itself, and `s^2 (psi'(a) + psi'(b))` goes to `s = phi mu (1 - mu)`, so the whole thing collapses to a weighted least squares of `logit(y)` on the predictor with weight `phi mu (1 - mu)`, the ordinary quasi-likelihood form, with no special functions at all. It runs at the same speed as the table (19.3 seconds) but the proposal is far enough off that effective sample size fell from 241 to **91**, so ESS per second was 4.7 against the table's 13.0. Worth recording because the derivation looks compelling on paper: the approximation is good in the bulk and useless in the tails, since `a = mu phi` falls to 1.7e-4 at a predictor of -8, where the asymptotic series diverges. A plain series is not an option for the same reason.

**A second, exact win in `ordbeta()`.** Its cutpoint updates were rebuilding the whole likelihood -- three log-gammas per observation -- on every slice evaluation, for a beta density that does not depend on the cutpoints at all. The likelihood splits: the endpoint-and-middle part depends on the cutpoints and not on `phi`, the interior beta density on `phi` and not on the cutpoints. Since slice sampling only needs the log density up to an additive constant, each update now evaluates its own half and drops the other. That is 30.6 seconds to 24.9. It is exactly invariant -- shifting the log density by a constant shifts the slice level by the same constant, leaving the slice set and the random-number stream untouched -- and the draws came out bit-identical, which is the check that it is a reorganization and not a change. `lgamma(phi)` was also hoisted out of both families' precision loops, being one of the three log-gammas `lbeta` computes.

**A claim in the earlier version of this entry was wrong.** It said the information is the expected one, so "two of the four calls buy only variance reduction and the observed version may be affordable". Working out the observed information shows it is `s^2 (psi'(a) + psi'(b)) - s (1 - 2 mu) C`, which needs *both* trigammas and then some. There was no saving on that route.

**No exact augmentation exists.** Asked directly, and the obstruction is structural. Writing the log density in the predictor,

    phi mu(eta) logit(y) + phi log(1 - y) - lgamma(mu phi) - lgamma((1 - mu) phi) + lgamma(phi),

`eta` enters through `mu = expit(eta)` in two places: inside a linear exponent, giving `exp(c expit(eta))`, and inside two log-gamma normalizers. The Polya-Gamma identity needs `e^{a eta} / (1 + e^{eta})^b`, with `eta` appearing *linearly* in the exponent; `exp(c expit(eta))` is not of that form, which is the same obstruction recorded for `ordbeta()` earlier. The gamma-ratio representation `Y = G1 / (G1 + G2)` does not help either: conditioning on the two gammas leaves `phi [mu log g1 + (1 - mu) log g2] - lgamma(mu phi) - lgamma((1 - mu) phi)`, with the normalizers untouched, because the difficulty is the eta-dependent normalizing constant rather than the shape of the kernel. The order-statistic representation needs integer shapes, and `mu phi` is not one.

## Nuisance parameters for a custom likelihood: the pinned-forest route

Recorded because the measurement overturned the first answer. The question was how to let `custom_family()` draw a nuisance parameter, and the initial recommendation was a declared interface -- `aux_start`, bounds, a log prior -- with the package slice-sampling each coordinate. The alternative raised was to give the nuisance a *forest pinned at depth zero*, which is a scalar with a normal prior drawn by machinery that already exists.

Two objections were offered against it and both were weak. The first, that the user would have to transform for positivity, is no objection at all: they already write `exp(eta[, 2])` for a real second predictor. The second, that pinning a forest and controlling its prior scale would need new structure, was wrong on inspection -- there is already one `Hypers` per forest, `sigma_mu` is already a per-forest vector, and `update_sigma_mu` is already applied as `hypers[h]->update_sigma_mu`. Only the R-side plumbing feeds them scalars.

**It works today with no changes**, on a hand-written Gaussian with a second predictor for `log(sigma)`, 800 observations, truth 0.400:

| | sigma | posterior sd | note |
|---|---|---|---|
| `num_trees = c(50, 1)`, read off one observation | 0.402 | 0.0188 | |
| the same, averaged over observations | 0.402 | 0.0103 | |
| splitting denied (`gamma` at 1e-8) | 0.404 | **0.0102** | |
| `gaussian()`, drawing sigma properly | 0.399 | 0.0097 | |
| theory, `sigma / sqrt(2n)` | -- | **0.0100** | |

**But the catch is the whole story: without pinning it is not a scalar.** The one tree splits, and the nuisance predictor took a single value in only **8% of draws** -- median 232 distinct values, maximum 800. So an unpinned "nuisance forest" is really an unadvertised dispersion surface fitting noise, which is why reading one observation gave a posterior 1.9 times too wide while averaging over observations landed on the theoretical width. Deny the splits and the posterior sd is 0.0102 against a theoretical 0.0100, with exactly one distinct value per draw.

That reframes the enabling change from "nice to have" to "the thing that makes it correct", and it is small: `gamma` as a vector, matching `num_trees`.

The residual cost is reporting rather than sampling. A pinned forest is a constant column of `fit$eta`, so `fit$aux` stays empty, `fit$rhat` has no entry for it, and the user extracts it with `exp(predict(fit, type = "link", draws = TRUE)[[2]][, 1])`. The attractive resolution is to route a pinned forest into `fit$aux` under a supplied name, which makes "pinned at depth zero" and "nuisance parameter" the same object and gives the declared interface's ergonomics with none of its sampler.

## Nuisance parameters for a custom likelihood, implemented

Built the pinned-forest route. `custom_family(aux_names =, aux_start =)` declares them, `logdens` gains a third argument, and they come back in `fit$aux` -- so the interface separates predictors from parameters while the mechanism is a forest pinned at depth zero for each one.

**What the engine does.** `Family::num_pinned()` is a new virtual, zero for everything but `RFamily`; the engine pins the trailing `num_pinned()` forests -- one tree, `gamma = 0` so no birth is possible, `update_sigma_mu = false` -- and keeps them out of `eta`, `sigma_mu`, the variable counts and the encoded forests, so `num_forest` is what the caller asked for. Almost no new structure was needed: `Hypers` was already per forest and `sigma_mu` already a per-forest vector. Per-forest `gamma` and `update_sigma_mu` stayed internal, as the caller has no reason to set them.

**Three bugs, each found by measurement rather than by reading.**

The first was mine and specific to this work. `RFamily::unpack` read the nuisance value once per block, from the first row, on the reasoning that a pinned forest has one leaf holding every observation and so shifts them all together. True of most calls -- but `log_f_pair_at` stacks *two* values of one component in a single block of `2n` rows to evaluate both ends of a Metropolis move in one pass. Both halves therefore got the first half's value, so the likelihood did not respond to the proposal at all: the reported change in the target was exactly the leaf prior's, `-(0.043^2) / (2 * 1.5^2) = -0.00041` against an observed `-0.000409`, which is what identified it. The nuisance column is constant only in *runs*, and the block methods now find them. For a family with no nuisance parameters, or any call that is not a paired one, there is exactly one run and nothing changes.

The other two were in the shared sampler, exposed rather than caused by this feature: a pinned forest's single leaf carries a whole parameter, so it can start far from its conditional mode, where an ordinary leaf in a fifty-tree forest never is.

**Newton's method can be thrown past the mode and then crawl.** On `f(s) = -n s - A exp(-2 s) / 2` the step is `s + 1/2 - n / (2u)` with `u = A exp(-2 s)`. From a start above the mode `u` is small, the step overshoots to about `-27`, and from there each iteration gains exactly `1/2` -- so `MAX_SCORE_STEPS = 50` is reached still far away. A trust region of four standard errors on the step converges in about six iterations instead.

**A Laplace proposal is an independence proposal, and cannot be accepted from far out.** This is the more interesting one. Once the fit was converging properly the proposal was still rejected every sweep: the target improved by 222 log points and the proposal cost 449, because the reverse density of a tight Gaussian at the mode, evaluated at a current value twenty standard errors away, is astronomically small. Damping the proposal's *location* toward the current value makes the step local; reversibility is kept by building the reverse proposal the same way from the proposed value, which is free because both come from the one shared fit, and when the cap does not bind both reduce to the fit itself, so the ratio is exactly the undamped one.

Neither cap binds for a chain near its mode: the full suite passes unchanged and the timings of every family are the same to within noise (`Beta()` 4.3 against 4.4 seconds, `location_scale()` 1.6 against 1.8, `gaussian()` 0.2 either way). Both are fixed constants, so the fit remains a deterministic function of the state, which is what the birth and death moves need.

**What it delivers.** A hand-written Gaussian with its scale drawn, started at 3.0 against a truth of 0.40, on 800 observations:

| | mean | posterior sd | ESS |
|---|---|---|---|
| `custom_family()` with a drawn `log_sigma` | 0.3988 | 0.0099 | 821 |
| `gaussian()`, conjugate step | 0.3990 | 0.0102 | 834 |
| theory, `sigma / sqrt(2n)` | -- | 0.0100 | -- |

Matching a conjugate Gibbs sampler on the mean, the spread *and* the effective sample size, from a start seven times off. Two parameters at once work as well: a t likelihood with a drawn scale and drawn degrees of freedom recovers `sigma` at 0.386 and puts the degrees of freedom at 49 on normal data.

**A note on what was not built.** The first design proposed for this was a declared interface -- `aux_start`, bounds, a log prior -- with the package slice-sampling each parameter, and the pinned forest was dismissed on two grounds that did not survive contact: that the caller would have to transform for positivity, which they already do for a real predictor, and that pinning would need new structure, which it did not. The measurement that settled it was that an *unpinned* nuisance forest is not a scalar at all: its one tree split, and the predictor took a single value in only 8% of draws. Pinning was the thing that made it correct, not a refinement.

## Survival, part two: proportional hazards and the discrete-time route

Two papers read and assessed. `ph()` was built; the other method was measured and left as a documented recipe, because it needs no code.

### Basak, Linero, Maringe and Rubio (2024), and the claim it retracts

The paper is *relative* survival -- excess hazard against population life tables, for cancer registries -- but its engine is a piecewise-exponential proportional hazards BART, and that is the part worth having. The observation that made it cheap: after their augmentation the log-likelihood in the predictor is

    delta_i * eta - Lambda_0(y_i) * exp(eta),

which is `a eta + b exp(eta)` at rate exactly +1 -- **the same exponential form `poisson()` already uses**. So a leaf update is one pass over the node, one row per subject, no data expansion. The baseline is a nuisance vector with an exact gamma conditional and an O(N + B) update, and the trailing-bin recursion in the paper is what keeps it O(N) rather than O(NB).

**This retracts something the documentation asserted in three places.** It said proportional hazards was "deliberately absent" because Cox's partial likelihood couples observations through risk sets and so does not decompose over the observations reaching a leaf. The first half is true and is why the *partial* likelihood cannot be used. The conclusion was wrong: the *full* likelihood of the piecewise-exponential model does decompose, approaches the partial likelihood as the bins shrink (Sinha, Ibrahim and Chen 2003), and is of a form the sampler already had a fast path for. Proportional hazards was reachable all along.

**Was it worth it?** Yes, and the measurement is two-sided, which is the useful part. 800 training and 800 test observations, 50 trees, 500 draws after 500 warmup, against a proportional-hazards truth under two baselines:

| Baseline | Family | Held-out log score | `S(t|x)` RMSE (worst t) |
|---|---|---|---|
| turns over | **`ph()`** | **-618** | **0.050 (0.063)** |
| | `loglogistic_aft()` | -866 | 0.081 (0.136) |
| | `weibull_aft()` | -876 | 0.100 (0.169) |
| Weibull | `ph()` | -870 | 0.046 (0.066) |
| | **`weibull_aft()`** | **-774** | **0.044 (0.059)** |

About 250 log points when the baseline turns over; about 96 back when the parametric baseline is right. Note the asymmetry in *where* the cost falls: on the Weibull truth the two are level on `S(t|x)` (0.046 against 0.044), so the flexible baseline costs the density at a point rather than the survival curve.

**What it does *not* buy, which is worth knowing before recommending it.** The risk *ordering*. An oracle PH model given the true baseline recovered `r(x)` at rmse 0.162 while the accelerated failure time families, after a best-case rescaling, managed 0.145 to 0.177 -- the log-logistic *beat* the oracle. A monotone reparameterization of time barely disturbs the ordering, so if all that is wanted is who is at higher risk, the parametric families were already enough. The case for `ph()` is the hazard's shape, the survival curve, and a hazard-ratio reading.

Two traps hit while prototyping, both mine: the starting intercept is folded into the offset at fit time, so predicting with a zero offset returns the forest *without* the level -- which does not matter for a centered comparison of `r(x)` and matters entirely for `S(t|x)`, where it first showed up as the oracle looking catastrophically worse than everything else. And `prepare_surv()` returned only `log_time`, so `a$time` was silently NULL, which surfaced as bin edges of `c(0, NA)`.

**Not built: the relative-survival part.** It is the paper's actual novelty and it is a narrow one -- it needs population life tables as an input, which is a data interface rather than a family. Worth knowing that it is *cheap* on top of what now exists: one extra Bernoulli draw per sweep, `d_i ~ Bernoulli(lambda_E / (lambda_E + lambda_P))`, with the population hazard entering as one number per subject. That is a nuisance draw and a data column, not a new sampler.

### Sparapani, Logan, McCulloch and Laud (2016), and why nothing was built

Not the same as our accelerated failure time families, so the question was live. It is a discrete-time hazard: expand to one row per subject per grid time up to their own, then probit BART on the binary event indicator **with time as a covariate**. Nothing is assumed about proportionality, so it handles crossing survival curves, which no other family here can.

But it needs no new code -- it is `binomial("probit")` on expanded data, and that is the *fastest* family in the package. So the deliverable is a recipe in the vignette, not a family. Measured against the other two, 600 training and 600 test observations:

| Family | PH truth | crossing truth |
|---|---|---|
| `ph()` | **0.035 (0.048)** | 0.065 (0.083) |
| `weibull_aft()` | 0.102 (0.183) | 0.064 (0.079) |
| discrete-time probit | 0.046 (0.072) | **0.036 (0.048)** |

A clean division. `ph()` is best under proportional hazards and the discrete-time model is a close second there, because it nests proportional hazards and pays only a little for the freedom. When hazards cross it is the only one that copes, and by a factor of 1.8 over both others -- which are then equally wrong.

**Possible follow-up, not done:** the ergonomics. The expansion is six lines and the survival curve is a running product, both now written out in the vignette, but a `survival_expand()` helper plus an `S(t | x)` output would make the route usable without the reader assembling it. The same missing piece applies to `ph()`, whose `S(t | x)` is also a documented recipe rather than a `predict()` type. That is the one interface gap this work leaves.

## The survival predict gap, and a marginaleffects bug behind it

Asked whether `predict()` for the survival families was *missing* something or whether `S(t | x)` was extra. It was missing, and the argument that settles it is the package's own consistency: every family whose response is discrete already returns its full predictive distribution through `type = "prob"` -- binomial, ordinal, multinomial -- while the survival families returned only a point summary, the median. The survival analogue of `prob` is `S(t | x)` over `t`, and it was not there. The corroborating smell was the recipe written into the vignette a moment earlier, which reached into `fit$family_opts$edges` and grepped `^lambda[0-9]+$`: asking a caller to depend on internal layout to get the primary output of a survival analysis is the sign the computation belongs inside.

`predict(type = "survival", times = ...)` now covers all four survival families. One column per time, or a draws by rows by times array. Checked against a closed-form exponential truth: RMSE 0.018 for `ph()` and 0.027 to 0.054 for the three accelerated failure time families, monotone in `t` for every row, inside (0, 1) everywhere, and the averaged matrix equal to the mean of the draws array.

**Looking for the payoff turned up a bug that had nothing to do with `times`.** The reason a survival curve matters is the estimand -- a contrast in t-year survival is what survival analysis is usually asked for -- so the natural test was `avg_comparisons(fit, type = "survival", times = 1)`. It failed. So did `avg_comparisons(fit, variables = "trt")` with no `type` at all, for `weibull_aft()` as well as `ph()`, while the same call on a binomial fit was fine. **marginaleffects has never worked for any survival family in this package.**

The cause: a survival response is a two-column matrix, a model frame keeps it as a single matrix column, and `data.table::setDT` reads its 2n cells as a column of length 2n -- "Column 2 ['trt'] is length 300 but column 1 is length 600; malformed data.table". An error naming data.table, from a package two steps away, for a reason in neither. `get_data.bartisan()` now splits matrix columns into ordinary ones; the estimands need the predictors, not the response.

With both in place the estimand recovers the truth. Contrast in survival between treated and untreated, truth `exp(-t e^0.8) - exp(-t)`:

| t | estimate | 95% interval | truth |
|---|---|---|---|
| 0.5 | -0.261 | -0.311 to -0.211 | -0.278 |
| 1 | -0.238 | -0.283 to -0.191 | -0.260 |
| 2 | -0.119 | -0.153 to -0.088 | -0.124 |

Every interval covers; the point estimates are shrunk toward zero by 0.03 to 0.05, which is the prior doing what it does.

**One limitation that is not ours to fix.** marginaleffects checks the dots against a whitelist hardcoded per model class inside `sanity_dots()` -- no option, no generic -- so it warns that it does not recognize `times`, while passing it through, which is what the warning says. Documented as expected rather than worked around, and the tests suppress it deliberately.

## CoxBART, and whether it differs from `ph()`

Asked whether the CoxBART of Linero, Basak, Li and Sinha (2022, sec. 2.3) is the model `ph()` already fits. **Not identical, but the same leaf-level target, and empirically indistinguishable.**

Their construction: the Cox partial likelihood is an *integrated* likelihood (Sinha, Ibrahim and Chen 2003). Give the discrete-time model `S(t | x) = exp{-e^{g(x)} sum_{t_l <= t} phi_l}` one jump per observed event time with the improper data-dependent prior `pi(phi) ∝ prod [delta_i phi_i^{-1} + (1 - delta_i) delta_0(phi_i)]`, and integrating the jumps out returns the partial likelihood exactly. Conditional on the jumps the likelihood is

    prod_i phi_i^{delta_i} exp[ delta_i g(X_i) - e^{g(X_i)} sum_{j : Y_j <= Y_i} phi_j delta_j ],

per observation `delta_i eta - C_i exp(eta)` -- **the identical target `ph()` has**, at rate +1. The whole difference is `C_i`:

| | `ph()` | `coxph()` |
|---|---|---|
| `C_i` | `Lambda_0(Y_i)`, piecewise linear over B bins | sum of the jumps at or before `Y_i`, a step process |
| baseline resolution | about `n^{1/3}` bins | one jump per event time |
| baseline prior | `Gam(a, b)` per bin, `b` drawn | improper, `Gam(0, 0)`-like per jump |
| marginal | the full likelihood | the Cox *partial* likelihood |
| leaf prior in the paper | -- | log-Gamma, for a conjugate integrated likelihood |

So `coxph()` is the `B -> infinity` limit of `ph()` under an improper prior -- which is exactly the remark Basak et al. (2024) make in passing and which I had recorded but not connected. The leaf prior is a further difference: the paper uses log-Gamma for conjugacy, this implementation keeps the package's normal leaf prior and corrects with the Metropolis step, so the two are slightly different models with the same intent.

`update_aux` is the whole of the new code: draw `phi_i ~ Gam(1, sum over the risk set of exp(eta))` for each event -- an exponential -- then accumulate. Both passes are O(N) over the pre-sorted times, with tie groups added in full before being read, which is the Breslow convention.

**Measured.** Four paired replicates on a baseline that turns over, 700 train and 700 test:

| | `S(t|x)` RMSE | `r(x)` RMSE |
|---|---|---|
| `coxph()` | 0.0360 | 0.128 |
| `ph()` | 0.0348 | 0.121 |
| paired difference | +0.0012 (sd 0.0024) | +0.0068 (sd 0.0111) |

Both within a standard error of zero. Across three truths, one replicate each, RMSE of `S(t | x)`:

| Truth | `coxph()` | `ph()` | `weibull_aft()` | `loglogistic_aft()` |
|---|---|---|---|---|
| baseline turns over | **0.043** | 0.044 | 0.093 | 0.087 |
| Weibull baseline | 0.061 | 0.062 | **0.052** | 0.063 |
| crossing, not proportional | 0.071 | **0.070** | 0.075 | 0.071 |

The two proportional-hazards families agree to the third decimal everywhere and win or lose together, which is what the shared target predicts. Same cost too: 5.7 against 5.6 seconds.

**External validation, which was the check worth doing.** On a linear truth, `coxph()`'s predictor correlates 0.982 with `survival::coxph()`'s linear predictor, with a spread 3% smaller -- the leaf prior shrinking toward zero. That is the evidence that the augmentation really is the partial likelihood rather than something adjacent to it.

**Recommendation, recorded in the docs:** `ph()`, for practical reasons rather than statistical ones -- a proper prior, a baseline reportable at full resolution, a genuine likelihood so the information criteria mean what they usually do, and prior weights, which `coxph()` refuses because the partial likelihood is derived without them. `coxph()` earns its place as the published method exactly, and as a benchmark.

## Can the bin choice be removed? Yes, and it turns out not to be worth it

Asked whether `ph()` could use one bin per event, or `coxph()` could gain weights and the likelihood-based tools -- combining the best of both. Both turned out to be possible, and the measurement then inverted the premise.

**`ph(num_bins = Inf)` was tried and does not work.** Built it, measured it, removed it. Over four paired replicates the survival function got *worse* by +0.017 (sd 0.006) on a turning-over baseline and +0.010 (sd 0.006) on a Weibull one -- both several standard errors, and 26% to 44% worse in relative terms. A smaller `lambda_shape` made it far worse still (0.175 against 0.056), which killed my first explanation: I had guessed the `Gam(a = 1, b)` prior was adding one pseudo-event per bin, but shrinking the shape should then have helped.

The real mechanism is that **`ph()` puts its prior on hazard *rates* with a common rate parameter.** A bin holding one event over a tiny width has a genuinely enormous hazard, of order one over its exposure; but `lambda_b ~ Gam(1 + A_b, b_lambda + B_b)` with a tiny `B_b` caps it at about `1 / b_lambda`. So the narrow bins are massively over-shrunk, and lowering the shape only removes what little data signal was there. A prior on rates cannot be spacing-agnostic.

**`coxph()` escapes that because it works in jumps.** A jump in the cumulative hazard is of order one over the risk set however narrow the gap, so the same prior is well-scaled at any resolution. That is the Gamma process prior of Kalbfleisch (1978), and adding it to `coxph()` is a small change: `jump_b ~ Gam(precision * hazard * gap_b, precision * gap_b)`, whose posterior is `Gam(prior + A_b, prior_rate + risk_sum_b)`. At `precision = 0` it is the improper prior the paper uses and returns the partial likelihood; positive precision gives a full likelihood and admits weights.

So `coxph(precision > 0)` does deliver everything asked for: no grid anywhere, prior weights, and an exact full likelihood -- `predict(type = "density")` now reproduces `fit$loglik` to 2e-13. And it is indistinguishable from `ph()` on both `r(x)` and `S(t | x)` at every precision tried, paired differences ranging from -0.0024 to +0.0024 with standard errors of the same size.

**And yet `ph()` is still the answer, for a reason that only showed up in the diagnostics.** `loo()` on `coxph(precision = 1)` gave an elpd 670 points worse than `ph()`'s despite a *better* in-sample likelihood. That is not a bug -- the density path was verified exact -- it is the fine baseline doing what a fine baseline does:

| | effective parameters (`p_loo`) | Pareto-k above 0.7 |
|---|---|---|
| `coxph(precision = 1)` | 674 | 56% |
| `ph()` | 17 | none |

674 effective parameters for 500 observations is leave-one-out reporting that it cannot do its job: with one jump per event time, an observation's own density is inflated by a parameter that only it informs. **The bin count is not a nuisance to be eliminated; it is the regularization that makes the model's own likelihood usable for comparison.** That reframes the original question rather than answering it, which is the useful outcome here.

Two bugs found along the way, both mine, both in the density path for `coxph()`: `logdens_unit` omitted the hazard factor altogether (it lived only in `reported_loglik`, where `ph()` carries it in `compute_eta_free`), and `set_aux` restored the baseline without rebuilding the cumulative hazard from it, so a restored draw evaluated at the constructor's prior mean. Together those made `loo()` meaningless in a way that looked like a modelling result. The `fit$loglik` against sum-of-log-density check is what caught them and is now a test.

Also fixed while here: the `bin_of` construction was an O(N * B) linear scan in both families, which is fine on a coarse grid and quadratic on a fine one; it is a binary search now. And `summary()` shows the ends of a long nuisance block with a count of what it omitted, rather than several hundred rows.

## The Henderson accelerated failure time model, implemented as `dpm_aft()`

`log T = m(x) + W` with `W` a mean-constrained Dirichlet process mixture and censored log-times imputed -- Henderson, Louis, Rosner and Varadhan (2020). The prediction was that both halves already existed and only the join was missing, and that held: `DPMAFTFamily` inherits from `DPMFamily` and adds about sixty lines. Everything about the mixture -- the Polya urn, the atom draws, the concentration, the centering, `error_density()` -- is inherited untouched. The CRTP static dispatch works through the inheritance because the derived class does not override `score_info_unit`, so `Concrete<DPMFamily>`'s qualified call still reaches the right one.

**What was added.** The imputation draws a censored log-time from the component it currently sits in, truncated below at its censoring time -- conditioning on the label is what makes it an ordinary Gibbs step, and the label is redrawn immediately afterwards given the value drawn. And an observed-data `reported_loglik` that credits a censoring with the mixture's *survival* rather than its density, with a matching `dpm_survival()` on the R side for the density and survival-curve routes.

**Measured**, 700 training and 700 test observations, held-out log score and RMSE of `S(t | x)`:

| Errors | `dpm_aft()` | `lognormal_aft()` | `loglogistic_aft()` | `weibull_aft()` |
|---|---|---|---|---|
| bimodal | **-607 / 0.029** | -818 / 0.098 | -847 / 0.100 | -879 / 0.115 |
| log-normal | -438.1 / 0.0264 | **-438.0 / 0.0266** | -444 / 0.031 | -467 / 0.057 |
| heavy tailed | **-509 / 0.036** | -556 / 0.065 | -516 / 0.040 | -566 / 0.072 |

210 log points and a third of the error in `S(t | x)` on a two-component error, and **within 0.1 log points of `lognormal_aft()` when a single normal is right** -- the same "costs nothing when the simpler assumption holds" property `dpm()` has against `gaussian()`, which is what makes it reachable-for rather than specialist. The price is speed: 5 to 15 seconds against about 1.

The correctness check that matters: `rowSums(predict(type = "density", log = TRUE))` reproduces `fit$loglik` to 2.8e-13, draw by draw. Those are two independent implementations of the observed-data likelihood -- one in C++ over the mixture atoms, one in R over the reported components -- agreeing including the censored contributions.

**A trap found while comparing, and worth knowing about.** `predict(type = "density")` is not on the same measure for every survival family: the accelerated failure time families report the density of `log T`, `ph()` and `coxph()` the density of `T`. Each is self-consistent -- the density of that family's own response -- but the caller supplies `(time, status)` to both and is not told. In the first head-to-head this showed as `ph()` scoring about 1000 log points worse than every accelerated failure time family; the correction is `sum(log t)`, which came to 1042, and on a common scale `ph()` was 168 points *better* rather than 1000 worse. Documented in `?predict.bartisan` and the vignette; `S(t | x)` is the metric that is comparable throughout.

**A note on the working tree.** `coxph()` in the tree is more developed than the version described in the entry above: it gained a `precision` argument putting a Gamma process prior on the jumps, with `precision = 0` the improper limit that reproduces the published partial likelihood. Two claims in the older NEWS entry -- that `coxph()` reports the partial likelihood and refuses weights -- are true only at zero precision, and have been corrected in place rather than left to contradict the newer entry above them.

## One proportional hazards family, and the bin count is not a knob

`coxph()` was built, measured, and removed. The question it existed to raise -- whether the bin count in `ph()` is a choice the caller is being made to make -- turned out to be answerable directly, and the answer removed the reason for a second family.

**Swept `num_bins` from 4 to 250, three replicates at 700 observations, against a baseline hazard that turns over and against a Weibull one:**

| Bins | `S(t|x)` RMSE | `r(x)` RMSE | `p_loo` | Pareto-k above 0.7 |
|---|---|---|---|---|
| 4 | 0.047 / 0.044 | 0.176 / 0.163 | 19 / 17 | none |
| **9** (default at this n) | **0.041 / 0.040** | **0.155 / 0.157** | 23 / 21 | none |
| 20 | 0.047 / 0.036 | 0.184 / 0.139 | 33 / 33 | none |
| 50 | 0.048 / 0.035 | 0.185 / 0.135 | 63 / 60 | none |
| 100 | 0.042 / 0.037 | 0.156 / 0.139 | 103 / 100 | none |
| 250 | 0.047 / 0.040 | 0.147 / 0.139 | 206 / 204 | 0.6% |

Turning-over baseline first, Weibull second. **The estimates are flat over a sixtyfold range** -- no trend in either column, every difference inside the replicate spread. And `loo()` works throughout: the Pareto-k diagnostics stay clean to 100 bins and are barely troubled at 250.

So the earlier framing was too pessimistic. The loo failure is not a consequence of "many bins": it appears only at *one bin per event*, where `coxph(precision = 1)` had 674 effective parameters for 500 observations and 56% bad k. Between nine bins and two hundred and fifty there is no penalty at all. The bin count is a regularization dial with a wide flat optimum, and the default sits in the middle of it.

**Which settles the design.** One family, `ph()`. `num_bins` is demoted to an advanced argument documented as being for checking the insensitivity rather than for tuning -- the caller never chooses it, and it provably does not matter. `coxph()` is gone: it offered no accuracy, and its one distinguishing feature, a grid-free baseline, is exactly what breaks `loo()`. The Gamma-process work is recorded above in case a grid-free version is ever wanted for its own sake.

A test now pins the insensitivity, so the claim cannot rot.

**Re-measured later, on the truths used in `vignette("survival")`** (700 observations, three replicates, default 9 bins), the plateau is narrower than this table suggests: flat from 4 to 100 bins, then a consistent ~20% rise at 250 in both truths and in both `S(t|x)` and `r(x)`, with `bad_k` reaching 1.3%. The replicate spread there is 0.009 and 0.005, so the 250-bin degradation is about 1.3 to 1.6 standard errors -- marginal individually, but it appears in both truths and both metrics, which is what makes it real rather than noise. The design conclusion does not move: the default sits near the bottom of a twenty-five-fold plateau. But "flat over a sixtyfold range" was too strong, and the vignette says "flat from 4 to 100, with over-parameterization visible at 250" instead.

## A survival vignette, and the comparison behind it

The survival section of `vignettes/families.Rmd` had grown to 151 lines -- five families, four embedded results tables, the discrete-time route, and three cautions -- inside a document whose stated purpose is to help a reader make a practical choice rather than to survey the modeling space. It was split out.

`vignettes/survival.Rmd` now holds the long version and `families.Rmd` keeps 43 lines: the table of five families with what each one's predictor means, one fit, `type = "survival"`, a four-line decision rule, and a pointer. Every results table moved.

**The comparison is new work, not a transcription.** Six data-generating truths crossed with the five families plus the discrete-time route, five replicates, 700 training and 700 held-out observations, 50 trees and 500 draws after 500 warmup -- and a censoring sweep from none to 70% on the truth where the families disagree most. Reproducible from `_dev/survival-sim.R`, `_dev/survival-bins.R` and `_dev/survival-results.R`; the last writes `vignettes/survival-results.rds`, which the vignette reads so that it builds without refitting. Figures use ggplot2, which was already in Suggests but had not been used in a vignette.

Two things the comparison established that were not known before.

**The log score can be made comparable across all five families.** The earlier vignette said it could not: the accelerated failure time families report the density of $\log T$ and `ph()` the density of $T$, so the two differ by $\sum \log t$. But censored observations contribute $S(t)$, which carries no measure at all, so the correction is `-sum(status * log(time))` -- events only. Applied, `weibull_aft()` and `ph()` score within a few points of each other on a Weibull truth instead of hundreds apart, which is the check that the correction is right rather than merely plausible. The vignette now gives the correction as a function rather than telling the reader to avoid the comparison.

**The families disagree far more about the density than about the ordering.** Across every proportional-hazards and accelerated failure time truth, all five recover the ordering of subjects by survival nearly perfectly whether or not they have the shape right; the differences are concentrated in $S(t \mid x)$ and the log score. The one exception is the crossing-hazards truth, where the ordering itself moves with time. So the choice of family matters for an absolute probability at a horizon and barely matters for a comparative question -- which is worth telling a reader before they agonize over it.

The rest confirmed what was already recorded: `ph()` wins when the baseline turns over, `dpm_aft()` wins by a wide margin on a bimodal error and costs nothing when a single normal is right, and the discrete-time route is the only option that copes with crossing curves while being a slower and slightly worse `ph()` under proportionality.

**One trap re-encountered.** A verification fit written as `Surv(time, status) ~ .` put `time` itself in as a predictor, and `dpm_aft()`'s reported `error_sd` came back at 0.021 against a truth of 0.7 -- which read as a scale bug in a family written by another agent. It was the formula. With the predictors named, `error_sd` recovers 0.716 uncensored and 0.684 at 27% censoring. This is the third time the `~ .` trap has cost a measurement in this project; the first draft of the new vignette had the same formula in four chunks and it was fixed there too.

## Renamed to bartisan

`genbart` became `bartisan` throughout: package, function, C++ namespace, documentation topics, file names, the repository directory. The fitted object's class is `bartisan_fit` rather than `bartisan`, which is the one part of the rename that is not a substitution.

**Order mattered, and getting it wrong would have been quiet.** `"genbart"` appears as a string in two different roles: as the class (`inherits(x, "genbart")`, `class(out) <- "genbart"`, `expect_s3_class(fit, "genbart")`) and as the package (`vignette("genbart")`, `asNamespace("genbart")`, `test_check("genbart")`, `future.packages`). Those go to different targets. So the class pass ran first, with context-anchored patterns, and only then the global `genbart` -> `bartisan` sweep, by which point no class usage was left to catch.

The S3 method suffixes needed the same care. The pattern is `\.genbart\b`, where the word boundary does not match before an underscore, so `predict.genbart` is renamed while the Rcpp entry points `.genbart_fit`, `.genbart_predict` and the rest are left for the global pass. A naive `.genbart` -> `.bartisan_fit` would have produced `.bartisan_fit_fit`.

Totals: 9 files in the class pass, then 92 files and 1504 occurrences in the global pass, then eight file renames through `git mv` so the history follows.

**Alignment had to be repaired, and needed three attempts.** `genbart` is seven characters and `bartisan` is eight, so every continuation line aligned to a paren on a renamed call sat one space short. The first attempt aligned to the *last* open paren on the line rather than the outermost unclosed one, which fixed inner arguments and left outer ones wrong. The second fixed the outer ones and thereby broke the inner ones again, because each shift moves the next nesting level. Iterating to a fixed point converged after one further sweep, 356 lines in total. Only lines exactly one space short were touched, which is the fingerprint of the rename rather than of hand layout, and a scan afterwards found none left.

**The directory move needed the sandbox off**, since renaming a folder writes to its parent, and `.Rprofile` needed it too, being a protected startup file. Note that after the move every write to the tree needs it, because the sandbox's writable root is pinned to the path the session started in. `.git` moved with the tree, the branch and history are intact, and its extended attributes are unchanged: it still carries no `com.dropbox.ignored`, which is correct and was left alone.

A snapshot of the tree went to the session scratchpad before any of this, since 82 files of uncommitted work were at stake and nothing here was committed.

Verified after the move: 1358 tests pass, `R CMD check` reports `Status: OK`, all nine vignettes rebuild, `bartisan:::.bartisan_optimized()` is `TRUE`, and no file outside `.git` contains the string `genbart`.

## The examples moved from birthwt to the RHC data

`MASS::birthwt` was 189 rows with one continuous outcome. The worked examples now use the SUPPORT right heart catheterization data: a binary `death` and the time `days` it took, so the same event serves a binary analysis and a right-censored one. The package ships a random 1500 of the 5735 patients, drawn under a fixed seed in `data-raw/rhc.R` so the shipped file is reproducible. Added as `data/rhc.rda`, built by `data-raw/rhc.R` from <https://hbiostat.org/data/repo/rhc.csv>, with the thirteen covariates from the maintainer's own worked example.

Verified the reconstruction against the reduced file the maintainer supplied: every shared covariate agrees to within 5e-06, which is CSV rounding, and `RHC` matches `swang1` exactly. One discrepancy worth recording: the published description calls `death` "died at 60 days", but the variable is death during follow-up. Its mean is 0.649, which matches the raw `death` column; 60-day death would be 0.404. The variable was used as it is.

**The change is not cosmetic, and three of the vignettes got better for it.** The causal vignette has real confounding in a known direction, since sicker patients were both more likely to be catheterized and more likely to die; the raw mortality gap is 9.2 points and adjustment brings it to 5.9 without changing sign. And the outcome being binary forced two genuine improvements: `pp_check()` is a weak check when there are only two values to get right, so `diagnostics` gained a decile calibration plot, and `comparison` now compares links rather than families, which came out as a clean negative (logit, probit and cloglog within 1.5 elpd, standard errors around 2).

**Two defects found on the way**, both since fixed rather than worked around; see "The subsample, and the two defects it stopped hiding" below.

**Build time was the binding constraint** at the full 5735 rows: a four-chain binomial fit took about 53 seconds and the survival fit about 63, and `R CMD check` rebuilding all nine vignettes ran to roughly half an hour. The 1500-row subsample is what brought that down. `comparison` still uses the default single chain, since leave-one-out needs draws rather than chains.

## The `. + external` formula warning, and the ACIC evidence for the propensity score

**The warning.** `bartisan(bwt ~ . + ps, data = d)` with `ps` living in the calling environment emits `'varlist' has changed (from nvar=9) to new 10 after EncodeVars() -- should no longer happen!`. It comes from base R's `terms.formula`, not from this package: `glm()` emits it identically on the same formula.

Isolated the trigger, which is narrower than it first looks. It is specifically `.` combined with a term that is not in `data`:

| formula | data | result |
|---|---|---|
| `y ~ . + ps` | `ps` in the environment | **warns** |
| `y ~ .` | `ps` a column of `data` | no warning |
| `y ~ a + b + ps` | `ps` in the environment | no warning |
| `y ~ .` | no extra term | no warning |

So `.` expands against `data`, the extra term is appended afterwards, the variable count changes, and R's C-level `EncodeVars` notices. Writing the formula out avoids it, and so does putting the variable in the data frame, which is what the vignette now does: `d$ps <- ps` and then `bwt ~ .`. One line, keeps the `.`, no warning, and the score is visibly part of the analysis dataset.

Deliberately not fixed inside the package. bartisan calls `terms()` the way `glm()` does, and matching `glm()`'s formula handling is a feature; pre-expanding the dot ourselves would duplicate base behaviour and risk diverging from it, to suppress a cosmetic warning that base R also emits for every other modelling function.

**The regularization-induced confounding section now leads with evidence rather than with absent functionality.** It previously named Bayesian causal forests and pointed at `bcf` and `stochtree`, which reads as an apology. It now keeps @hahn2020 for the mechanism, since that is his contribution, and reports what the 2016 Atlantic Causal Inference Competition found: @dorie2019 for the competition, and @carnegie2019, who decomposed which features of a BART fit mattered.

Read the Carnegie comment rather than citing it from the abstract, which changed what I wrote. Bias was small for every BART variant and including the propensity score cut average absolute bias by about a tenth, so on bias it is a refinement rather than a rescue. The differences showed in interval coverage: base BART covered at 83.4%, and the best combination without targeted learning, ten chains plus symmetric intervals plus the propensity score, reached 91.9%. Running several chains contributed alongside the score, which is a second reason for the `chains = 4` the vignettes already use. Carnegie also notes that ignorability and overlap held by construction in every competition dataset, so the value of modelling treatment assignment is plausibly larger when they are strained, which is exactly the case the section is about; that caveat is carried over.

**Two smaller fixes to the maintainer's edits.** An inline expression used `ac <- .Last.value` to pull the interval bounds into the prose; `.Last.value` is not set inside a knitr chunk, so the object was not the estimand and `abs()` failed on it. The result is now assigned in the chunk and referenced, which achieves the intent, and `conf.hi` was corrected to `conf.high`. Separately, removing the per-fit seeds shifted every number, and one prose claim landed on a knife edge: the third subgroup interval now ends at exactly 0.000, which made "two of the three intervals exclude zero" both awkward and fragile. That sentence no longer depends on the pattern.

## The causal vignette, revised

Six changes on the maintainer's instruction, most of them correcting things I had got wrong rather than matters of taste.

**The propensity score is now fitted with `bartisan()`** rather than `glm()`, which is the obvious thing given what the vignette is about, and it changes the picture: the BART score spans 0.13 to 0.78 where the logistic one spanned 0.03 to 0.97. The extreme logistic values are the artefacts a saturated parametric model produces when a covariate pattern happens to be perfectly predictive, and they make a positivity assessment look worse than it is. Overlap is plotted with `cobalt::bal.plot()`, added to Suggests. The working call passes the score as `distance = data.frame(prop.score = ps)`; `var.name = "prop.score"` with the score supplied separately is rejected.

**The weighted analysis is removed.** Two reasons, both mine to have caught. Propensity weights are not frequency weights, and a bartisan fit treats prior weights as though they were, so it is not established that a posterior interval from a weighted fit has the coverage a weighted estimator needs. Until that is settled the intervals should not be reported as they stand. Separately, and independently of the package, **g-computation for the ATE should not pass weights to `avg_comparisons()`**: the averaging that turns individual contrasts into an average effect is over the target population's covariate distribution, which for the ATE is the sample as observed, so weighting it again applies the reweighting twice. My original version did exactly that. Both points are now stated in a short section rather than demonstrated.

**The outcome family is `dpm()`, not `gaussian()`.** It estimates the error's shape instead of assuming it, costs almost nothing when a normal would have done, and is the default for a numeric outcome anyway, so using `gaussian()` here was a step backwards from what the package does on its own.

**Added: the potential outcomes** through `avg_predictions(fit, variables = "smoke")`, which reports $E[Y(0)]$ and $E[Y(1)]$ at 3072 and 2773 grams. Their difference is the ATE reported in the next section, and showing both is more informative than the difference alone.

**Added: a moderation analysis.** `by = "race"` gives three subgroup effects (-342, -306, -245) of which two intervals exclude zero and one does not, which is exactly the pattern people misread as moderation. `hypothesis = ~pairwise` gives the three differences, all with intervals covering zero comfortably. The vignette says plainly that comparing whether one interval excludes zero and another does not is not a test, and that these should be read as three noisy estimates of one effect. Note that `hypothesis = "pairwise"` as a string is rejected by the current marginaleffects; the formula form is required.

## The debug build did not link, and -O2 had been hiding it

Reported from RStudio: after `pkgbuild::clean_dll()`, loading the package failed with `symbol not found in flat namespace '__ZN7bartisan27OrdinalLogitAugmentedFamily9OMEGA_MINE'`, plus an unused-variable warning in `model.cpp`.

**Three things had to line up, which is why it stayed hidden.**

1. `OrdinalLogitAugmentedFamily::OMEGA_MIN` and the matching member of `LoglogisticAFTAugmentedFamily` are `static constexpr double`. Such a member is implicitly `inline`, and so needs no out-of-line definition, only from C++17 onward.
2. The package declared no `CXX_STD` and `R CMD config CXX` emits no `-std=` flag, so the compiler's own default applied. On this toolchain that is C++14: `clang++ -dM -E -x c++ /dev/null` reports `__cplusplus 201402L`.
3. `std::max()` takes both arguments by const reference, so `std::max(rpg(...), OMEGA_MIN)` binds the member to a reference. That is an ODR-use and needs the symbol.

`R CMD INSTALL` builds at `-O2`, where the constant is folded into the instruction stream and no reference is emitted, so the missing definition never mattered. `pkgbuild::compile_dll()` uses `-UNDEBUG -Wall -pedantic -g -O0`, where it does.

Reproduced in isolation before changing anything, which is what pinned all three factors at once:

| standard | -O0 | -O2 |
|---|---|---|
| gnu++11 | **1 undefined** | 0 |
| gnu++14 | **1 undefined** | 0 |
| gnu++17 | 0 | 0 |

**Fixed twice over.** `src/Makevars` now sets `CXX_STD = CXX17`, which is the honest declaration: the code already relies on C++17. And both call sites were rewritten as `drawn < OMEGA_MIN ? OMEGA_MIN : drawn`, which is what `std::max` is defined to compute, so the results are unchanged while nothing binds a reference. The second change means the package links under any standard, not just the one now declared.

Checked afterwards that no other `static constexpr` member is at risk: `WEIGHT_TOL`, `PROBIT_DIRECT`, `TAB_N` and `TAB_L` all appear only in arithmetic and comparison, which are lvalue-to-rvalue conversions rather than ODR-uses. That matches the single missing symbol in the report.

`model.cpp`'s `total_trees` was genuinely dead, left by the refactor that split reported forests from engine forests; the identically named variable 300 lines later is a different function and is used. A clean `-Wall -pedantic -O0` build now emits zero warnings.

**Incidental correction.** With this fixed, `R CMD check` outside the agent sandbox reports `Status: OK`: no warnings, no notes. The `OMP: Warning #179` line reported in every earlier entry, and an `nm` cache-file NOTE, are both artifacts of the sandbox denying writes to `TMPDIR`, not properties of the package.

## The remaining five workflow vignettes, and what writing them found

`effects`, `importance`, `diagnostics`, `comparison` and `causal` written, completing the series. All use the same birth weight data as `workflow`, so a reader is never learning a new dataset and a new idea at once, and each ends by pointing at the next.

**One bug, found the same way as the last one: by running the documented workflow on real data.** `marginaleffects::predictions(fit, newdata = )` returned the wrong number of rows -- covered in the previous entry. Nothing new surfaced in the five vignettes themselves beyond the two items below, which were investigated and left alone deliberately.

**`avg_slopes()` looked like a bug and is not.** On the birth weight fit it reports several hundred grams per pound where a one-unit comparison gives about 2.5. The cause is `x_transform = "quantile"`, the default: the transform is the predictor's empirical distribution function, so the fitted function is piecewise constant in the original scale and a difference quotient is either zero or a whole step over a tiny denominator. Measured directly, perturbing `lwt` by 0.0001, 0.01 or 0.1 changes the prediction by exactly nothing; only at 1, the integer spacing, does it move. With `x_transform = "range"` the derivative is stable across step sizes (3.51, 3.51, 3.48) and `avg_slopes()` agrees with `avg_comparisons()` to two decimals.

I built a guard that raised an error from `get_predict()` and then removed it, for two reasons. The behaviour is already documented at `?bartisan-marginaleffects` and pinned by a test that asserts exactly this inflation, so it was a deliberate decision rather than an oversight, and overriding it unilaterally was wrong. And the guard could not be softened into a warning: marginaleffects swallows warnings raised inside `get_predict()`, so only an error is deliverable, which is too blunt for documented behaviour. The caution now lives in `vignette("effects")`, where a reader meets the function.

Two notes from building the guard, in case it is ever revisited. The call can be identified reliably: `avg_slopes()` puts `slopes` in the call stack and passes `internal_call` in the dots, while `comparisons()` and `predictions()` do neither. And match the function name exactly rather than searching the stack for the word, because a helper named `stop_if_slopes` matches its own frame, which cost a debugging cycle.

**Content worth keeping across the series.** The correlated-predictor demonstration in `importance` is the sharpest: two nearly identical columns, only one in the truth, and the forest gives the copy `prop_used` 1.00 and 75 splits while the real variable gets 0.096 and 0.1 splits. `avg_comparisons()` follows the usage, attributing 1.42 to the copy and 0.00 to the original. Moving both together recovers the true 1.5. So the model has the relationship and cannot say which column owns it, which is the honest statement of what importance can and cannot support.

The `comparison` vignette's family comparison came out as a negative result and is reported as one: `gaussian()`, `location_scale()` and `dpm()` on the same predictors land within 1.5 elpd of each other with standard errors around 1, and loo flags the differences as indistinguishable. The plain family is adequate here, which is the outcome that makes the flexible ones safe to try.

## The bartisan vignette as the theory document, and two measurements behind it

`vignette("bartisan")` rewritten. It had been a family-by-family tour that `families` now covers, opening with a first-model walkthrough that `workflow` now does better, and it contained no mathematical statement of the model at all. It is now the reference document: the sum-of-trees model and its generalization to an arbitrary density, the three parts of the prior, soft rules, sparsity, then backfitting, the reversible-jump tree moves, the Laplace approximation, the three target shapes, augmentation, nuisance parameters and the Dirichlet process mixture. References added to `references.bib` from Zotero: Hill (2011), Hill, Linero and Murray (2020), Hahn, Murray and Carvalho (2020), Murray (2021).

**The sparsity claim I made earlier was wrong.** `sparsity = TRUE` is and was the default; the comparison that produced the claim had changed `num_trees` at the same time and I attributed the difference to the wrong argument. `vignette("workflow")` has been corrected.

**Measured properly, on the Friedman function with five real predictors and five noise ones, three replicates at n = 500, mean `prop_used` for the noise:**

| Trees | `sparsity = FALSE` | `sparsity = TRUE` |
|---|---|---|
| 10 | 0.28 | 0.09 |
| 20 | 0.50 | 0.08 |
| 50 | 0.95 | 0.09 |
| 100 | 1.00 | 0.14 |

The real predictors sit at 1.00 in every cell. So the prior works, and it works at every tree count. The "use fewer trees for variable importance" advice comes from @chipman2010 itself, verified in the paper: counting splits "is less effective when m is large because the redundancy offered by so many trees tends to mix many irrelevant predictors in with the relevant ones". The `sparsity = FALSE` column reproduces that exactly. The DART prior addresses the same problem directly, so the advice is a workaround for its absence rather than a general property of variable selection. Their own conclusion anticipates this: "Prior specifications for variable selection via BART are part of our ongoing research."

On `birthwt` nothing separates even with sparsity on, and five added pure-noise columns interleave with the real predictors at 0.76 to 0.92. That is n = 189 with a weak signal, not a failure of the prior: there is nothing for it to concentrate on.

**The `eta` row of `fit$rhat` is not a convergence problem, and more draws makes it worse.** Sweeping the draw count on the Friedman function at n = 400, four chains: 500 draws gives eta rhat 1.29, 2000 gives 1.33, 5000 gives 1.37, with ESS around 10 throughout, while `sigma` sits at 1.03 in all three. The per-observation distribution is elevated as a whole (median 1.09, 79% above 1.05), so it is not merely the maximum being an extreme-value statistic.

Checked against `dbarts` on identical data, which is **worse**: median per-observation rhat 1.22, 90th percentile 1.46, maximum 1.87, and 99% above 1.05, against bartisan's median 1.09 and maximum 1.29. So slow mixing of the fitted values is a property of BART samplers rather than of this one, and bartisan is the better of the two here. The vignette says so without naming the comparison; the measurement is recorded here.

Worth considering for the `rhat` table: reporting a quantile of the per-observation values, or the proportion above a threshold, alongside the maximum. A single worst-case number over hundreds of observations reads as alarming and is not actionable.

## The getting-started vignette, and the marginaleffects bug it uncovered

`vignette("workflow")` written as the "Getting started" entry point: one complete analysis of `MASS::birthwt`, every step a single call at the defaults. Deliberately shallow, with a pointer at the end of each section to the vignette that goes deeper.

**Writing it found a real bug, which is the argument for writing vignettes against real data.** `marginaleffects::predictions(fit, newdata = )` returned the wrong number of rows whenever `newdata` was not the whole training frame: 1 row in gave 25 out, 3 gave 5, 10 gave 12. `lm` on the same data was correct, so it was ours.

Tracing what `get_predict.bartisan()` actually received showed marginaleffects passing five rows for a three-row request, the first two carrying `rowid = -1`. Those are marginaleffects' own scratch rows, and it drops them again by that marker after the predictions come back. `get_predict.bartisan()` was rebuilding the column as `seq_len(ncol(draws))`, which overwrote the markers with 1 and 2, so nothing was dropped and the scratch rows landed in the output. Carrying `newdata$rowid` through when it exists is the whole fix.

The effect estimates in the vignette moved by about 10% once this was corrected (smoking went from -277 to -292 grams), so every `avg_comparisons()` number produced before this was mildly wrong. It never showed up in the test suite because the existing marginaleffects tests all used the default `newdata`, which is the one case that worked.

**Two smaller findings, recorded rather than acted on:**

- `fit$rhat` is `NULL` with the default `chains = 1`, so a caller who never sets `chains` gets no convergence diagnostics at all and no indication that any exist. The vignette uses `chains = 4` and explains why. Worth considering whether the default should be higher, as it is in most Bayesian packages.
- `variable_importance()` on a default fit is close to uninformative: every predictor sits between 0.92 and 1.00 on `prop_used`, because without `sparsity = TRUE` nothing is ever excluded. The separation only appears with sparsity on. Worth considering whether `sparsity` should default to `TRUE`, or whether `variable_importance()` should say so when it is off.

**Overlap created.** `vignette("bartisan")` still opens with "What this package is for" and "A first model", which `workflow` now does better, and its middle is a family-by-family tour that `families` now covers. Its genuinely distinctive material is the mechanics: soft rules, missing values, the conditional density, chains and speed. It has no mathematical statement of the model at all. Reorganizing it into the theory document is the natural next step and has not been done.

## `weibull_aft()`: 15% for free, and why the rest needs a decision

Profiled the same way as `dpm_aft()`, and the answer was the opposite: 99.6% of the fit is in `.Call`, so this one is genuinely the sampler.

Sampling the process put `bartisan::Family::score_info_unit` -- the *base class* fallback -- at the top of the family functions. Two findings behind it:

1. `AFTFamily` was declared `: Family`, not `: Concrete<AFTFamily>`, so it ran the generic accumulate loops at four non-inlinable virtual calls per observation. It was the only remaining survival family not using the CRTP path: the two augmented AFT families and `PHFamily` all do.
2. It did not override `score_info_unit`, so the score and the information came from two separate virtual calls, each forming `r` and the exponential again.

Fixed both -- `Concrete<AFTFamily>` plus a fused `score_info_unit` writing the same expressions in the same order, so the results are bit-identical. Checked that nothing derives from `AFTFamily` first, since making a class `Concrete` statically binds the accumulate loops to it (the `DPMAFTFamily : DPMFamily` case is the one place that pattern already appears and is safe only because it overrides no unit function).

**9.68 -> 8.19 seconds.** Only 15%, so the virtual dispatch was not the story.

The rest is the exponential-form Laplace machinery itself, and that is inherent rather than a defect: `poisson()`, which takes the same `TARGET_EXP_UP` path on comparable data, costs 3.6 seconds against `gaussian()`'s 1.17 -- a 3x penalty for the iterative mode-finding over the closed-form quadratic draw. Cost depends only weakly on `sigma` (8.87 at 0.25, 8.18 at 1.0), so the rate being `-1/sigma` rather than `+-1` is not a lever either.

**What would actually fix it, and why it was not done.** The Weibull AFT error is exactly `log(E)` with `E ~ Exp(1)`, and the standard route to conditional Gaussianity for that is the Frühwirth-Schnatter finite mixture-of-normals approximation to the log-exponential density -- the same device used for Poisson regression. It would plausibly bring the family to roughly `lognormal_aft()`'s ~1.1 seconds. But it is an **approximation**, and every other augmentation in this package is exact; adopting it would weaken a claim the package currently makes cleanly. That is a design decision for the maintainer, not a bug fix, and it is no longer urgent now that `weibull_aft()` is not the default. Recorded rather than implemented.

## The `Surv` default moved to `dpm_aft()`

Once the speedup below landed, `dpm_aft()` was the second-cheapest survival family and the most accurate on average, so it replaced `weibull_aft()` as the family inferred from a `Surv` response.

Evidence: best or tied-best on four of six truths, never worse than third, and level with the correctly specified family on the two truths where one existed; 2.8 seconds against `weibull_aft()`'s 9.7.

**The cost of the change is interpretive, not statistical.** The reported predictor changes meaning -- `weibull_aft()`'s is a log time ratio, `dpm_aft()`'s is $E[\log T \mid x]$, which is a time ratio only if the error is symmetric, which is the assumption the family exists to avoid. `type = "survival"` and `type = "response"` are unaffected, being on the same scale for every survival family.

The maintainer was asked whether to flag that estimand change at fit time and chose not to, after the trade-off was put to them. The standard inferred-family message was kept -- it is uniform across every response type and names the family chosen, so removing it for `Surv` alone would have made survival the one response type that changes model silently. No *extra* estimand warning was added. The distinction is documented in `?bartisan`, both vignettes, and NEWS instead.

One consequence needed handling: `dpm_aft()` cannot take prior weights, so `default_family()`'s existing weights guard, which covered `dpm()`, was extended to cover it, with a message naming the three `*_aft()` families and `ph()` as the weighted alternatives.

Fixed in passing: `README.Rmd` still described the numeric default as `dpm()` "with ten or more distinct values" and `gaussian()` below that. That threshold was removed some time ago -- every numeric response gets `dpm()` -- so the README had been wrong about it independently of this change.

## `dpm_aft()` was ten times slower than it needed to be, and the sampler was not why

Reported as too slow to be a default: 13.3 seconds against `lognormal_aft()`'s 1.2 on 700 observations, 50 trees, 500 draws after 500 warmup.

**The localization mattered more than the fix.** Four measurements, each cheap, narrowed it without guessing:

1. `dpm()` against `gaussian()` on the same log-times: 1.29 against 0.93. The Dirichlet process mixture costs 0.36 seconds. So the mixture machinery is not it.
2. `dpm_aft()` against `dpm()` on *identical uncensored data*, fitting the same model: 12.3 against 1.24, with near-identical fitted mixtures (28.4 against 25.7 clusters). A 9.9x gap that nothing about the model explains.
3. Censoring swept 0 to 60%: 12.4, 8.1, 7.5, 7.2 seconds. **More censoring is faster**, so the truncated-normal imputation is not it either -- and that sampler is inversion-based, not rejection-based, so it was never a candidate.
4. Tree count swept 1, 5, 50: 7.65, 6.71, 7.22 seconds -- **flat**, where `dpm()` scaled 0.17 to 0.75. A fixed per-sweep cost, not the tree loop.

`sample(1)` on the running process then put 40% of the main thread in `Rf_pnorm5`, reached through `math3_2` in `stats.so` and `R_doDotCall` -- that is R's *vectorized* `pnorm`, called from R code, not from `bartisan.so` at all. `Rprof` finished it: `stats::pnorm` 66% self time, `dpm_survival` 83% total, `.Call` 11%.

**The defect.** `response_scale()`'s `dpm_aft` branch finds the median survival time by bisecting the mixture's survival function, 60 steps per draw. It initialized `lo <- rep.int(-30, ncol(e))` and `hi <- rep.int(30, ncol(e))` -- one entry per observation -- and called `dpm_survival(object, s, mid)` on that vector. But `dpm_survival()` is a function of the error value alone: the median is a property of the error distribution that every observation shares, and only the shift `e[s, ]` differs between them. Every entry of that length-700 vector computed the same number. 500 draws x 60 steps x 28 components x 700 observations is 588 million `pnorm` evaluations to produce 500 scalars.

Bisecting a scalar and adding it to the row: 13.3 -> 2.27 seconds. Hoisting `mixture_at()` out of the bisection, which was rebuilding the component matrix on all sixty steps of each draw: 2.27 -> 2.17. `dpm_survival()` took an optional `components` argument for the second.

**Both are pure redundancy removal.** The fitted values are bit-identical -- checked by running the old algorithm verbatim against the new `predict()` output, max absolute difference exactly 0.

**What this says about the model.** The sampler was never slow. `DPMAFTFamily` inherits `DPMFamily`'s `TARGET_QUADRATIC`, because conditional on the component labels and the imputed log-times the model is exactly Gaussian -- the same fast path `lognormal_aft()` reaches by augmentation. Its C++ time is 0.9 seconds against `lognormal_aft()`'s ~1.1. **So there was nothing to augment and nothing a warm start would have helped**: burn-in length was not the constraint, and the two obvious "make the model cheaper" routes would both have been wasted work. Standing order that this vindicates: measure where the time is before designing a speedup.

Final standing, 700 observations at the package defaults: `lognormal_aft()` 1.10, `loglogistic_aft()` 1.23, `dpm_aft()` 2.79, discrete-time probit 6.26, `ph()` 6.53, `weibull_aft()` 9.68. The most flexible family is the second cheapest, and the *default* for a `Surv` response is the most expensive.

## A note on working alongside concurrent edits

Partway through this the tree turned out to contain a `dpm_aft()` family and a roxygen note about the density measure differing across the survival families, neither of which came from this session. The maintainer had implemented the Henderson et al. (2020) model recorded as a To Do above while this work was in progress.

Nothing was lost -- `DPMAFTFamily` sits well after where `CoxPHFamily` was, so removing the latter by locating its banner comments did not touch it, and it still resolves and fits. But the near miss is the lesson: a removal that finds its target by index between two markers is only as safe as the assumption that nothing has moved. Checking for unexpected identifiers before a sweeping edit, rather than after, is the cheap version of that check.

## Design decisions

The whole family-specific surface of the algorithm reduces to three per-observation quantities, as a function of the additive predictor `eta`:

- `logdens(i, eta)` — the log density or mass of observation `i`
- `dlogdens(i, eta, h)` — first derivative with respect to component `h`
- `info(i, eta, h)` — minus the second derivative, or the Fisher information

Everything else — tree structure, birth/death/change moves, the Laplace proposal, the sparsity prior — is shared. New families are therefore additions to `src/family.cpp` alone. Where analytic derivatives are awkward the base class falls back to central differences, as Linero's own Weibull and generalized-gamma code does, but with a larger step for the second derivative: his `1e-6` amplifies round-off by `1e12` in the second difference.

Soft rules drop into the same framework by the chain rule. A leaf contributes `w_i * mu` to observation `i` for a membership weight `w_i` in `[0, 1]`, so the gradient picks up a factor `w_i` and the information a factor `w_i^2`. Hard rules are the case `w_i` in `{0, 1}`, which recovers Linero's expressions exactly, so there is a single code path — and because multiplying by 1.0 is exact, a hard tree stores no weights at all and the two paths still agree bit for bit.

Multi-predictor families (multinomial, location-scale) carry `H` independent forests, each with its own sparsity prior and leaf scale, and the family exposes partial derivatives with respect to each. `H = 1` is the common case.

**Two shapes of target are exploited where the family declares them.** A leaf value enters the predictor linearly, so the log target over a leaf inherits the shape of the log density — and where that shape is known, one pass over a node determines the whole function and everything after it is arithmetic. `Family::target_form()` returns `TARGET_GENERAL`, `TARGET_QUADRATIC`, `TARGET_EXP_UP` or `TARGET_EXP_DOWN`. See the two log entries on the conjugate shortcut and the exponential form.

### Departures from the reference implementation

Four deliberate differences from Linero's `FlexBart`, each verified:

1.  **Soft rules with a non-conjugate likelihood.** Linero (2025) names this an open extension. Soft rules make the two child leaves of a split dependent, so the reversible-jump move gets a *bivariate* Laplace proposal, with the off-diagonal information `sum_i w_iL * w_iR * info_i`. That term is identically zero for hard rules, so the bivariate proposal collapses to Linero's independent pair and one implementation serves both.
2.  **Deterministic Laplace fits.** `FlexBart` starts Fisher scoring from the node's previously stored mode and stops at `|score| < sqrt(info)/10`, which makes the resulting proposal depend on the sampler's history rather than only on the current state — a quiet violation of detailed balance. Here every fit starts at zero and converges to a fixed tolerance, so the birth and death moves provably build the same proposal, which Linero (2025) calls essential.
3.  **Corrected death-move transition ratio.** The published `R_DEATH` has the primes on `|L|` and `|NOG|` swapped relative to its own derivation, and `FlexBart` evaluates the reverse birth probability on the pre-collapse tree, so it uses 0.5 even when collapsing the root, where the correct value is 1. Both are fixed.
4.  **Bounded slice sampler.** `FlexBart`'s interval expansion is an unbounded `while (true)` that spins forever if the log density returns a non-finite value, which is reachable when a nuisance parameter wanders into a region where the likelihood underflows. All three loops are capped.

## `diagnose()`, and a drift statistic that had to be abandoned

**One call for convergence and mixing**, in `R/diagnose.R`. Split-R-hat, bulk and
tail effective sample size for every scalar the sampler draws, for the fitted
function over its worst 5% of observations, and for the total number of splitting
rules in the forest -- then the checks that failed and what to change. All from
the stored draws; no other package.

Three things it does that the plain `fit$rhat` table does not:

- **Works with one chain**, by folding it into halves, rather than reporting
  nothing. The one-chain warning is still the first line of advice.
- **Keys the checks to the share of a row's components that failed**, not to the
  worst one. The worst of a thousand per-observation R-hats is extreme even when
  every chain has converged, so a threshold on a maximum condemns every fit; the
  table shows the worst-5% boundary for reading and the checks use the share.
- **Repeats R-hat on the second half of the draws alone.** Discarding the early
  retained draws is exactly what more `num_burn` would have done, so if that
  fixes R-hat then warmup was the problem and the advice says only "raise
  `num_burn`"; if it does not, the chains have settled in different places and
  the advice is the other list. That decomposition is what makes the output
  advice rather than a menu.

**The drift statistic was built, calibrated, and thrown away.** A Geweke-style z
on each chain's first half against its own second half would answer the warmup
question directly. Three variance estimators were calibrated against stationary
AR(1) series, where by construction there is nothing to find:

| estimator | fires at rho = 0.995 (want 5%) | misses a 6-sd trend |
|---|---|---|
| each half's own effective sample size | 31% | no |
| the whole chain's effective sample size | under 3% | 73% of the time |
| batch means, sqrt(n) batches | 92% | no |

The first is anti-conservative because a window shorter than the autocorrelation
time overstates its effective sample size; the second is blind to trends because
the trend inflates the autocorrelation estimate that sets its own error bar; the
third because batches of length sqrt(n) are still correlated. A forest is sticky
enough to sit where all three fail, and there is no threshold that separates slow
mixing from non-convergence in one chain -- which is why the literature uses
multiple chains. Re-running an already-calibrated statistic on a subset avoids
the problem. The measurements are in the roxygen block so the next person does
not repeat them.

**A feature comparison against the other packages** is now in `README.Rmd`,
checked against the installed versions rather than from memory. Two things it
corrected: *flexBART* 2.0.3 is no longer the minimal package an earlier session
described -- it has a formula interface, a `family` argument, heteroskedastic
regression, VCBART, the DART prior and four chains by default -- and *stochtree*
supports ordinal cloglog, which was not obvious from its exports. The table names
what this package lacks: no threads inside a chain, no cross-validation over
hyperparameters, no grow-from-root warm start, no JSON serialization, no
recurrent-event or competing-risks survival, and no formal variable-selection
test. Partial dependence, interaction detection and estimands go through
*marginaleffects* and are marked as a helper's rather than as this package's.

Also fixed: the vignette rename in 3b5eb75 left three dangling
`vignette("workflow")` references, in `effects.Rmd`, `implementation.Rmd` and
`diagnostics.Rmd`. `workflow.Rmd` is now `bartisan.Rmd`.

## The benchmark measured the wrong thing twice, and a progress bar

**Two errors in `_dev/benchmark.Rmd`, in opposite directions.** The reported gap
to dbarts, "within a factor of 1.8", is really about three.

The first: the Gaussian section called `bartisan(y ~ ., data = dtr)` with no
`family`. A numeric response with no family reaches `default_family()`, which
returns `dpm()` -- a Dirichlet process mixture for the error distribution. So
every row labelled "Gaussian" timed a DPM against three packages fitting a
Gaussian, and its RMSE column described a different model. Worth 1.35x.

The second, and larger: the competitors were called as `bart(xtr, ytr, xte, ...)`,
which evaluates a thousand test points at every draw *inside* the timed call,
while bartisan's `predict()` ran after the timer stopped. dbarts is 0.281 s
without the test matrix and 0.467 s with it, so the comparison charged dbarts for
two thirds again as much work. Fixed by giving nobody test data in the timed call
and predicting afterwards, which needed `keeptrees = TRUE` for dbarts and a
post-fit `predict()` for BART and stochtree.

Two settings that looked unfair and are not, measured in
`_dev/parity-benchmark.R`: the Dirichlet sparsity prior and the drawn leaf scale,
neither of which dbarts has, are both free to within noise (0.95x and 0.99x).
`keeptrees = TRUE` costs dbarts 1.03x. Matched, the ratio is 2.87.

**What dbarts does that this package cannot.** Reading its source: nodes hold a
pointer into one shared index array, so a split is an in-place partition with no
allocation and no copy; each node caches a mean and an effective count computed
in the same recursion as that partition, so the marginal likelihood reads two
doubles; and predictors are pre-discretized to `misc_xint_t` cutpoint codes, so
the split test is an integer compare over one or two bytes and the partition is
hand-vectorized with runtime dispatch to `partition_avx2`, `_sse4_1`, `_sse2` or
`_neon`. BART has none of the three -- `getsuff()` loops over all n on every move
-- which is most of why it is six times slower than dbarts here.

All three are unavailable under a soft rule, where an observation reaches both
children with a weight: there is no partition, the target is not a function of two
summaries, and the gate needs the covariate value rather than a side. They are
available under `gate = "hard"`, but only by splitting a code path that soft and
hard share on purpose (`Node` caches `idx` plus `wt`, and `split_support()`
already skips the weights when rules are hard). Not attempted: it is a second
sampler for the configuration that loses on accuracy.

**A C++ audit found no correctness bugs.** Checked: the birth and death move
probabilities, including the two boundary cases where the root is a stump or is
being collapsed, and their evaluation order relative to the tree mutation -- all
correct, and `p_birth_after_death()` documents a bug in Linero's reference code
that this fixes. The node pool recycles only leaves and deletes anything else,
so no subtree leaks. `ForestGuard` is RAII and already covers the interrupt and
R-error paths. All five `-Wfloat-equal` sites are exact `-Inf` sentinels or
bit-identical cache invalidation. Every unguarded-looking `c(cat - 1)` and
`cuts(k - 1)` is in fact guarded upstream. Of 1819 warnings under
`-Wall -Wextra -Wconversion -Wshadow -Wold-style-cast`, 1650 are `int`-into-
`arma::uword` noise, 123 are unused virtual-override parameters, and the two real
ones were shadowed names, now renamed.

**A progress bar, through progressr.** The sampler is one C++ call, so it reports
by calling back into R: `control$progress` is a function of no arguments and
`control$progress_ticks` says how many times to call it. `R/progress.R` builds
the reporter in the calling session, which is what lets it work under
future.apply -- the closure is captured by the engine, sent to each worker, and
the conditions it signals are relayed back. Sized for every chain at once, so
four chains fill one bar once. Capped at 50 reports per chain: the callback is
nothing next to a sweep but a handler that redraws a bar is not.

Verified at 50/50, 150/150, 150/150 and 150/150 for one sequential chain, three
sequential, three under multisession and three under multicore. The multisession
case appears to fail under `devtools::load_all()` and does not: the worker runs
`library(bartisan)` and gets the installed build, so the test is only meaningful
after `R CMD INSTALL`. Progress consumes no randomness -- the same seed gives the
same draws whether or not anything is listening, which is a test.

## Log: how the sampler got fast

The whole arc, in one sentence: the complaint was a probit fit to `lalonde` taking **82.3 s** where dbarts took 0.20 s, and it is now **0.72 s** with hard rules — but almost none of that came from making the code tighter. It came from noticing that the leaf-level target has an exploitable *shape* far more often than it looks like it does.

### Where the time goes now

Instrumented per move type, n = 1000, p = 10, 50 trees, 500 + 500 draws.

| | hard | soft |
|---|---|---|
| total | 0.454 s | 2.157 s |
| birth | 20.8% | 13.3% |
| death | 17.1% | 8.3% |
| change | 36.2% | 17.6% |
| leaf refresh | 18.1% | 9.2% |
| bandwidth move | — | 46.5% |
| &nbsp;&nbsp;of which `rebuild_support` | — | 34.1% |
| &nbsp;&nbsp;of which the likelihood difference | — | 5.5% |

For **hard** rules the four targets inside the moves — the family evaluations — are what is left, and they are irreducible without changing the model. For **soft** rules it is the bandwidth move, and within that the gate evaluations that rebuilding every membership weight requires. Both remaining costs are the algorithm rather than the implementation.

### Three traps that corrupted measurements before any of this

Worth stating first, because every number below would have been wrong without them, and two of the three cost real time more than once.

**`-O0` builds, three times.** `devtools`, `pkgload` and `roxygen2::roxygenize()` all compile at `-O0` by default (pkgbuild's "extra flags"), which is a five- to twentyfold slowdown that looks like nothing at all: the fits are correct, the tests pass, only the clock is wrong. The first occurrence overstated costs by up to 14x. The third was a user reporting timings 5–20x mine on comparable hardware; reproducing it took one build and matched on all three cases to within 2%:

| | reported | `-O0` build | `-O2` build |
|---|---|---|---|
| probit soft, direct | 68.9 s | 70.1 s | 13.6 s |
| probit soft, augmented | 37.6 s | 37.5 s | 2.40 s |
| probit hard, augmented | 16.6 s | 16.5 s | 0.75 s |

It is now *detected* rather than remembered. `__OPTIMIZE__` is defined by the compiler when it is optimizing, so `bartisan:::.bartisan_optimized()` reports the truth about the loaded library, `bartisan()` warns once per session when the answer is no, and `_dev/benchmark.Rmd` refuses to run at all. The project `.Rprofile` sets `options(pkg.build_extra_flags = FALSE)` so `load_all()` produces an optimized library in the first place.

**Stale object files.** `roxygenize()` leaves `-O0` objects behind and a subsequent `R CMD INSTALL` sees them as up to date. Any timing claim must come from a clean build; `devtools::test()` is enough to contaminate them again.

**Nothing beyond `-O2` is worth setting.** Measured on this package, `-O3` is within 1%, `-mcpu=apple-m1` makes no difference, and `-flto` is slightly *worse*. The optimization level is not the variable; being at `-O0` is.

### Why it was slow: the decomposition that made the rest possible

Guessing would have been hopeless, so the sampler was instrumented to count how many times it visits an observation, against the number a conjugate sampler needs. That separates *how much work* from *how expensive each unit of work is*, and the two had completely different explanations.

| | Visits | vs conjugate minimum | Time | ns per visit |
|---|---|---|---|---|
| gaussian, hard | 108.9M | 4.4x | 0.66 s | 6.1 |
| gaussian, soft | 192.9M | 7.9x | 1.19 s | 6.2 |
| logit, hard | 124.0M | 5.0x | 1.53 s | 12.3 |
| probit, hard | 130.5M | 5.3x | 8.66 s | **66** |
| probit, soft | 226.7M | 9.2x | 15.75 s | 69 |

Fisher scoring was converging in 1.9 to 2.8 steps, not the 50 it was allowed, so the visit count was only 4–9x a conjugate sampler's. The whole probit anomaly was **66 nanoseconds per visit against the Gaussian's 6**. Multiplying out — 5.3x the visits, 4x the per-visit machinery, 11x for probit's special functions — gives 233x against the 216x observed. The decomposition is complete and none of the three factors is the one that would have been guessed.

**The general lesson, which drove everything after it.** The expensive part of a non-conjugate sampler is not the non-conjugacy of the *likelihood*; it is that the leaf-level target is not quadratic. Wherever a data augmentation exists that makes it quadratic, the whole Laplace apparatus collapses into a conjugate draw and the cost falls by most of an order of magnitude.

### The conjugate shortcut: where the target is quadratic, the Laplace fit is exact

Not approximately exact — exact. **A leaf value enters the additive predictor linearly, so if the family's log density is quadratic in the predictor then the log target is quadratic in the leaf value; and the leaf prior is Gaussian, so it is quadratic too.** One pass over the node determines the entire function. Writing the target's value, score and information at any reference point as `(f0, d1, d2)`,

```
log f(mu) = f0 + d1 (mu - ref) - d2 (mu - ref)^2 / 2
```

is an identity, not a Taylor approximation, and the Laplace fit is `mean = ref + d1/d2`, `sd = d2^{-1/2}` — the conditional posterior itself. Fisher scoring reaches the mode in one step from anywhere, the fitted normal *is* the conditional posterior, and the leaf refresh is a Gibbs step with acceptance one. The two-child case is the same statement with a gradient and a 2x2 curvature.

Implemented as `Target1` and `Target2`, which all four moves use uniformly. Where the family is quadratic they do one fused pass in the constructor and answer everything from arithmetic; where it is not, they forward to exactly what the sampler did before. That uniformity is the point: there is no second copy of the move logic to drift.

Two things fell out of getting it right.

**`Context::log_f_score_info*()`, one pass for three quantities.** The sampler previously took separate passes for the log target and for its derivatives. They are wanted at the same point, so they are computed together.

**The general path was stopping short for quadratic targets.** Fisher scoring broke out when the score fell below `SCORE_TOL` standard errors — sensible when convergence is asymptotic, wrong when one step from anywhere is exact. It was returning a mode off by up to a hundredth of a standard error, and a worse proposal for it. Fixed by skipping the tolerance check when `quadratic` is set, which also made the two paths comparable: they agree to 1e-15 over twenty iterations across soft and hard rules and four families. Before the fix they differed by 3e-5, which took tracking down and was the tolerance, not an error in the closed forms.

`exact_quadratic = FALSE` forces the general path, which is what makes the agreement testable. Measured on lalonde: probit hard 1.50 → 0.94 s, gaussian hard 1.42 → 0.92 s.

### The exponential form: a second shape, for the count models

Hill et al. (2020, sec. 3.1.5) observe that the count models share one shape. For a single additive predictor it is `c + a*eta + b*exp(s*eta)`, and it sits alongside the quadratic form as a second case the sampler can exploit: three numbers from one pass over a node determine the whole function, so the Laplace fit iterates on scalars rather than on the data.

The coefficients need no new family interface. From one evaluation of the target's value, score and information at any point, with the prior's contribution subtracted, `b*exp(s*mu0) = L''(mu0)`, `a = L'(mu0) - s*L''(mu0)` and `c = L(mu0) - a*mu0 - L''(mu0)` recover all three — expressed reference-free, so the fit does not depend on where the expansion was taken, which is what the reverse move needs in order to rebuild the same proposal.

Covers **Poisson** (`s = +1`) and **gamma with a log link** (`s = -1`). **Hard rules only**: a soft rule gives each observation its own exponent `exp(s*w*mu)`, and a sum of those is not a function of three numbers, so `Target1` and `Target2` check and decline. Worth **1.86x on Poisson and 1.89x on gamma**.

Two things had to be got right. **The gamma family was reporting the expected information, not the observed one.** Most families that do so have a reason — the observed version can go negative — but the gamma's response is strictly positive, so its observed curvature `shape * y / mu` cannot. It now reports the true curvature, which is what the recovery reads its coefficients off, and which incidentally makes its Laplace fit a genuine second-order match. And **with hard rules the two children separate exactly**: every observation reaches one child or the other, so nothing contributes to the cross curvature and `Target2` becomes two independent one-dimensional problems.

Verified by checking the premise rather than the outcome: fit the three coefficients from three evaluations and predict elsewhere. Exact to 5.7e-14 for Poisson and gamma, off by 1.16 for binomial — so the check can fail.

### Data augmentations: what each is worth

Every augmentation rewrites the likelihood as the margin of a Gaussian one (or, for the negative binomial, a Poisson one), which makes the target a shape the sampler can exploit. Every one of them also costs mixing, because the chain now has to move the latent variables and the predictor in alternation — so the ratio to judge is effective samples per second, and it differs a lot by family.

| Family | Rules | Speed | ESS | ESS/s | Verdict |
|---|---|---|---|---|---|
| `ordinal("probit")` | hard | 26–30x | 0.73–0.90x | **22–24x** | default |
| `ordinal("logit")` | hard | 14.5–15.3x | 0.87–1.02x | **13–15x** | default |
| `ordinal("probit")` | soft | 14x | 0.79–0.95x | **11–13x** | default |
| `ordinal("logit")` | soft | 7.0–7.2x | 0.70–0.82x | **5.0–5.9x** | default |
| `binomial("probit")` | either | 5.7x | 0.66x | **3.8x** | default |
| `binomial("logit")` | either | 2.6x | 0.81x | **2.1x** | default |
| `negbin()` | hard | 1.7–1.9x | 0.61–1.14x | 1.2–2.0x | default |
| `negbin()` | soft | 1.1x | 0.71x | 0.8x | off |
| `multinomial()` | either | 4.2x | 0.37x | 1.6x | by name |

**Probit** is Albert and Chib (1993): a latent normal with mean `eta` and unit variance, truncated to the sign of `y`, has `P(z > 0) = Phi(eta)`. This is the dbarts insight, and it is why dbarts fits a probit model as fast as a Gaussian one.

**Logit, negative binomial and multinomial** rest on one identity (Polson, Scott and Windle 2013): a likelihood proportional to `exp(kappa psi) / (1 + exp(psi))^b` becomes, after introducing `omega ~ PG(b, psi)`, proportional to `exp(kappa psi - omega psi^2 / 2)`, which is Gaussian in psi. Only what psi, kappa and b are differs:

| | psi | kappa | b |
|---|---|---|---|
| binomial logit | `eta` | `s - w/2` | `w` |
| negative binomial | `eta - log theta` | `(y - theta)/2` | `y + theta` |
| multinomial | `eta_j - log C_j` | `y_j - n/2` | `n` |

The multinomial case is the one worth explaining. Conditional on the other categories, the likelihood of category j is *exactly* binomial-logistic in `eta_j - log C_j`, where `C_j` is the sum of `exp(eta)` over the others. So no stick-breaking decomposition is needed — the structure the sampler already has, updating one forest at a time, is the structure the augmentation wants. It does mean the latent variables are redrawn per forest rather than per sweep, since `C_j` moves during a sweep; hence `Family::before_forest()`.

**The Polya-Gamma sampler.** Devroye's alternating-series method for integer `b`, which is exact; for non-integer `b`, the sum-of-gammas representation at 20 terms with the tail replaced by a gamma matched on **two** moments. Matching the variance as well as the mean is what buys the shorter sum — the first version used 200 terms with the tail matched on the mean alone, so this is 10x cheaper. What is left out is a third-moment discrepancy in a component carrying 0.3% of the standard deviation. Validated against the closed-form mean and variance of PG(b, c) across 35 combinations of `b` in {0.5, 1, 2, 3.7, 5, 20, 40} and `c` in {0, 0.5, 2, 6, 20}: every mean within 2.2 standard errors, every variance within 2%.

**The negative binomial abandoned Polya-Gamma entirely.** Its `b = y + theta` is not an integer, so every draw went through the series: roughly 200 gamma draws per observation per sweep against the 150 observation visits the forest update itself cost. The augmentation was more expensive than the likelihood it replaced, and measured **0.5x** — a net loss. The Poisson-gamma mixture is the better route and the one the review points at: a negative binomial count is a Poisson count whose rate is drawn from a gamma, and introducing that rate leaves the predictor's contribution as `-theta*eta - lambda*theta*exp(-eta)` — the exponential form, with **one** gamma draw per observation as the entire cost. `theta` is still drawn from its collapsed conditional on the true negative binomial likelihood, with the rate redrawn after it, which is a valid partially collapsed Gibbs step in that order. That turned a 0.5x loss into a 1.2–2.0x gain under hard rules. The range is two problems of different shape and is quoted as a range rather than at the flattering end: this is the one augmentation whose value genuinely depends on the data, and it is worth turning off if its diagnostics look poor.

**The ordinal probit** is the largest gain in the package, because the target it replaces is the most expensive: `log_prob` is a difference of two cumulative normals plus a log, and with no exploitable shape every trial value of a leaf parameter costs its own pass. Profiling said the work per sweep was identical to a Gaussian fit's — 1000 observation-visits per tree, exactly the sample — so it was never bookkeeping; the leaf update alone was 43% of the fit. Conditional on the latent normal all of it collapses to one pass of arithmetic:

| | direct | augmented | speed | ESS | correlation with truth |
|---|---|---|---|---|---|
| K = 3, soft | 50.2 s | 3.56 s | 14.1x | 63 → 50 | 0.9843 → 0.9854 |
| K = 3, hard | 23.0 s | 0.76 s | 30.2x | 51 → 37 | 0.9766 → 0.9735 |
| K = 5, soft | 52.5 s | 3.79 s | 13.9x | 66 → 63 | 0.9889 → 0.9899 |
| K = 5, hard | 25.1 s | 0.98 s | 25.6x | 41 → 37 | 0.9800 → 0.9785 |

**The cutpoint update is the part to be careful about.** Their full conditional given the latent variables is uniform between the largest latent value in the category below and the smallest in the category above. That is Albert and Chib's own step and it is a trap: the interval has width O(1/n), so the cutpoints barely move and the sampler gets *worse* the more data there is. Instead they are drawn from the ordinal likelihood with the latent variables integrated out — the existing slice sampler, reused unchanged — and the latent variables are redrawn immediately afterwards. That is a partially collapsed Gibbs sampler (Van Dyk and Park 2008) and the standard remedy (Cowles 1996). The existing cutpoint code was factored into `update_ordinal_cuts()` so the direct and augmented families share one implementation rather than two that can drift.

**The ordinal logit needed no Kolmogorov-Smirnov sampler**, which is the route the literature takes (Holmes and Held 2006) and which would have needed an alternating-series sampler of its own — the same shape of risk as the saddlepoint sampler declined below. Polson, Scott and Windle's Theorem 1 at `a = 1`, `b = 2` reads

```
e^x / (1 + e^x)^2  =  (1/4) E[exp(-w x^2 / 2)],   w ~ PG(2, 0)
```

and the left-hand side *is* the standard logistic density. So a logistic residual is a mean-zero normal with precision `w`; and since PG(b, c) is PG(b, 0) tilted by `exp(-c^2 w / 2)`, the conditional of the precision given a residual `r` is exactly **PG(2, |r|)** — an integer-parameter draw, which the exact Devroye sampler already in `polyagamma.cpp` covers. Nothing approximate enters and no new sampler was needed.

Checked against the density before building on it: Monte Carlo of `(1/4) E[exp(-w x^2/2)]` matched `dlogis(x)` to within 0.07% out to x = 5, and the tilting property separately. Both are in the test suite, because the whole family rests on that one identity. Measured 15.3x (hard) and 7.0x (soft), with cutpoints matching the direct fit to two decimals; the gain is about half the probit's because the precision is an extra draw per observation that the probit's unit-variance latent does not need.

**Poisson is not helped at all**, which was worth checking since it looks like it should be. Polya-Gamma applies to likelihoods proportional to `exp(kappa psi)/(1 + exp(psi))^b`; the Poisson's `y eta - exp(eta)` is not of that form and no substitution makes it so. Nor are gamma, the AFT families, ordered beta or location-scale. Fruhwirth-Schnatter et al. (2009) give a Poisson augmentation through inter-arrival times and a finite normal mixture, but the mixture is an approximation to a log-gamma density rather than an identity, so it would trade exactness for speed in a way nothing else here does.

**Correctness was checked on the reported likelihood, not on timings.** For each augmented family the reported log likelihood is compared against the original model's, computed by hand in R — `dbinom` with trial counts, `dnbinom`, the multinomial from its own fitted probabilities, the ordinal from its cutpoints — and agrees to 1e-8 or better. That is the sharp test: the sampler's target is the augmented density, so anything wrong in the rewriting shows up as a likelihood on the wrong scale.

One diagnostic trap: comparing the augmented and direct multinomial samplers on `eta` gave a correlation of 0.898, which looks like a bug. It is not. The symmetric coding is unidentified, so `eta` includes a direction pinned only by the prior, and two chains need not agree on it. Comparing identified quantities — fitted probabilities — gives 0.995 at 1500 draws and 0.997 at 5000. **Always compare identified quantities.**

### DPMBART: a Dirichlet process mixture for the error distribution

BART assumes i.i.d. normal errors, and that assumption does most of the work in
its uncertainty quantification. `dpm()` drops it, following George, Laud, Logan,
McCulloch and Sparapani (2019): each observation gets its own error mean and
variance drawn from a Dirichlet process, so the error distribution is whatever
mixture of normals the data ask for.

**It is cheap because the target does not change.** Conditional on the mixture the
log density is `-(y_i - mu_i - eta)^2 / (2 sigma_i^2)` -- still exactly quadratic
in the predictor, so `TARGET_QUADRATIC`, the closed-form leaf draw, and no
Laplace approximation. The only cost over a Gaussian fit is the mixture update:
1.7 s against 1.3 s at n = 1000, 50 trees, 500 + 500 draws. The paper reports the
total roughly doubling; here it is about 1.3x, because the Gaussian baseline this
is measured against is itself faster than theirs.

The baseline `G_0` is the conjugate normal-inverse-chi-square rather than this
package's half-Cauchy, and it has to be: the Escobar and West draws that make the
mixture update a few lines *are* the closed forms conjugacy provides. Everything
in it is calibrated off a linear fit the way BART calibrates its own scale prior
-- `nu = 10` and `q = 0.95` (tighter than BART's 3 and 0.90, because the mixture
covers small errors with extra components rather than with one component's left
tail), and `k_0` set so the marginal of a component mean reaches the largest
residual at `k_s = 10` of its own scale units. The concentration gets Rossi's
tapered prior, and is drawn on a grid: `P(I = k | alpha)` is proportional to
`alpha^k Gamma(alpha) / Gamma(alpha + n)` times a Stirling number that does not
involve `alpha` and so cancels, which makes the grid update exact rather than
approximate.

**The mixture is reported the way the trees are.** Its component count changes
every draw, so it goes into a flat vector with one offset per draw --
`mixture_flat` and `mixture_start`, next to `forest_flat` and `tree_start` --
rather than into a matrix. `combine_chains()` grew a shared helper for
concatenating that shape, which the forests now use too.

**The likelihood is the mixture's own predictive**, occupied components weighted
by their sizes plus the Dirichlet process's chance of opening a new one, whose
kernel is the baseline's marginal `t`. It is computed twice -- in C++ for
`fit$loglik`, and in R from the stored components for
`predict(type = "density")` -- and the two agree to 4.5e-13, which is the check
that the stored components really are the ones the sampler used.

#### Only the sum of the fit and the error mean is identified

Nothing forces the mixture to be centred, so `f` and `E[e]` are individually
unidentified. Measured on a heavy-tailed example at n = 500: the level of the
predictor and the error mean had standard deviations of 0.750 and 0.751 across
draws with a correlation of **-0.987**, while their sum had a standard deviation
of **0.122** -- six times smaller. The bias of `type = "response"` (the sum) was
-0.003 against 0.144 for `type = "link"` (the trees alone).

So `type = "response"` is the conditional mean and the thing to compare against a
truth or across families, and `type = "link"` carries the drift. Documented, and
the benchmark below uses the response scale for every family so the comparison is
of like with like.

#### Measured against gaussian and heteroskedastic gaussian

The paper's three error distributions plus a heteroskedastic one, which is the
case `location_scale()` exists for. n = 1000 train and test, `f(x) = 10x^3`,
50 trees, 500 + 500 draws, two replicates. RMSE, coverage and width are for the
regression function on the training data; the score is the held-out predictive
log density.

| errors | family | seconds | RMSE | coverage | width | score |
|---|---|---|---|---|---|---|
| normal | `gaussian()` | 1.3 | 0.263 | 0.949 | 0.988 | **-2097** |
| | `location_scale()` | 12.2 | 0.255 | 0.958 | 0.973 | -2100 |
| | `dpm()` | 1.7 | 0.258 | 0.946 | 0.957 | **-2097** |
| t3 | `gaussian()` | 1.3 | 0.264 | 1.000 | 1.301 | -2642 |
| | `location_scale()` | 11.8 | 0.281 | 0.997 | 1.319 | -2607 |
| | `dpm()` | 1.9 | **0.251** | 0.999 | **1.206** | **-2467** |
| skewed | `gaussian()` | 1.3 | 0.277 | 0.983 | 1.587 | -2680 |
| | `location_scale()` | 12.1 | 0.250 | 1.000 | 1.564 | -2680 |
| | `dpm()` | 1.6 | **0.190** | 1.000 | **0.993** | **-2456** |
| heteroskedastic | `gaussian()` | 1.3 | 0.284 | 0.921 | 0.831 | -1886 |
| | `location_scale()` | 14.3 | **0.243** | 0.963 | 0.786 | **-1734** |
| | `dpm()` | 1.7 | 0.303 | 0.925 | 0.779 | -1853 |

Four readings.

**It costs nothing when the errors are normal.** Identical score to `gaussian()`
to four figures, and the same RMSE and coverage. That is the property the paper
spends its prior specification on -- "the strengths of the standard BART approach
is not lost when the errors are close to normal" -- and it holds.

**It wins by a wide margin when the errors are heavy tailed or skewed.** 175 and
224 log points of held-out score over `gaussian()`. On the skewed case it also
cut the interval width for the regression function from 1.59 to 0.99 *while
holding coverage*, and cut RMSE by a third. That is the paper's real point: BART's
intervals under non-normal errors are not merely wide, they are the wrong shape.

**It is not the family for heteroskedasticity, and the measurement says so
plainly.** `location_scale()` wins that row by 119 log points and `dpm()` is
worse than `gaussian()` on RMSE there. A mixture makes the error distribution
flexible but keeps it the *same* distribution at every `x`. The two families
answer different questions and the documentation now says which is which.

**`location_scale()` costs about nine times a Gaussian fit** -- 12 to 14 seconds
against 1.3 -- where `dpm()` costs 1.3x. That is two forests and a target that is
not quadratic in the second one, against one forest and a mixture update. Worth
knowing when choosing between them on a problem where either would do.

### The family is inferred from the response when none is named

`family` defaults to `NULL` and is read off the response: `Surv` or a two-column
matrix of times and events gives `weibull_aft()`, an ordered factor `ordinal()`,
a logical or a two-level factor or numeric zeros and ones `binomial()`, a factor
with more than two levels `multinomial()`, and anything else `gaussian()`. A
message says which was chosen, and naming `family` is what silences it -- the
same action that changes it, so there is no second argument for the message.

The family is now resolved *after* the model frame is built rather than before,
because the response is not available until then. Nothing else moved.

**Two rules are deliberately not what a reader might guess, and both are
refusals to guess.** A count is not read as `poisson()`: a non-negative integer
response is often Poisson and often not, and the Poisson variance assumption is
strong enough that making it silently would be a modelling decision taken on the
caller's behalf. Gaussian is the weaker guess and the one whose failure is easy
to see. And a numeric response with two values that are not zero and one --
`c(1, 2)` -- is Gaussian rather than binomial, because which value is the success
is not something to infer.

One thing this broke and how: `is_binary()` first read "at most two values, all
in {0, 1}", which made a constant response of all ones a binomial with no
failures. That turned an existing test's clear "no variation" error into a silent
fit. It now reads "exactly two values", so a constant numeric response goes to
Gaussian and complains, which is what a degenerate response should do whatever
the family.

### The zero-inflated families: two latent variables, 4 to 10x in ESS/s

This was on the To Do list as "an exponential-form route for the zero-inflated
families", and the framing was wrong in a useful way: what blocks them is a
*mixture*, not a link, and the fix is one augmentation before the exponential
form is even reachable.

The zero contributes `log[pi + (1 - pi) P_0(mu)]` -- a log-sum-exp of the two
components -- so neither predictor has a shape and `dlogdens_unit` for a zero
costs several transcendentals on top of that. Introducing `z_i`, the indicator of
whether observation i is a structural zero, separates them completely:

- the **count** forest sees `prod over {z = 0}` of the count likelihood, which is
  a plain Poisson (exponential form) or negative binomial;
- the **inflation** forest sees a Bernoulli logistic likelihood in `z`, which is
  the Polya-Gamma case already in the package, so its target is quadratic.

`z` is drawn from its exact conditional: zero whenever `y > 0`, and for `y = 0`
one with probability `pi / (pi + (1 - pi) P_0)`. It is redrawn before each forest
rather than once a sweep, so each forest moves under the indicator the other has
just been fitted with.

For `zi_negbin()` a second augmentation goes on top: the non-structural
observations get the Poisson-gamma rate, which turns their target from a general
one into the exponential form. `z` is drawn with that rate integrated out and the
rate redrawn immediately afterwards -- a partially collapsed Gibbs step in that
order (Van Dyk and Park 2008), and better than conditioning `z` on a stale rate.
`theta` is still drawn from the true zero-inflated likelihood. The rate of a
structural zero is never drawn, because its observation's contribution to the
count target is multiplied by `1 - z`.

Measured at n = 500, p = 5, 50 trees, 500 + 500 draws, with ESS the median over
observations of the effective sample size of the fitted mean:

| | rules | speed | ESS | ESS/s | RMSE |
|---|---|---|---|---|---|
| `zi_poisson()` | hard | 7.1x | 1.42x | **10.1x** | 0.447 vs 0.501 |
| `zi_poisson()` | soft | 4.6x | 0.85x | **3.9x** | 0.394 vs 0.387 |
| `zi_negbin()` | hard | 9.4x | 0.84x | **7.9x** | 0.570 vs 0.648 |
| `zi_negbin()` | soft | 5.9x | 0.98x | **5.8x** | 0.462 vs 0.562 |

The gain is much larger than the "smaller than the multinomial's" this was
predicted to be, and the reason is the direct family rather than the
augmentation: the log-sum-exp made it the second most expensive target in the
package after the ordinal probit's. Note that it pays under soft rules too, where
the count forest gets no shape shortcut at all -- half the benefit is simply that
a per-unit kernel of arithmetic replaced one of transcendentals.

On by default for both kinds of rule.

### The multinomial-Poisson transformation: implemented, measured, rejected

The other half of that To Do item was Murray's (2021) route to the multinomial: a
multinomial with total `n_i` and probabilities `softmax(eta_i)` is the conditional
law of independent Poissons with rates `lambda_i exp(eta_ij)` given their sum, so
introducing `lambda_i` with the scale-invariant prior `p(lambda) ∝ 1/lambda`
makes the categories independent Poissons and every forest gets the exponential
form. Integrating `lambda` back out returns the multinomial likelihood exactly,
so it is an augmentation and not a reparameterization.

**Two things recorded here were wrong.** First, this was described as "a larger
change than the negative binomial's was" because it "changes the identification"
-- that the extra Poisson total would loosen the symmetric and reference codings.
It does not: `lambda` is integrated out exactly, so the marginal posterior of the
predictors is the one the direct family targets and neither coding is disturbed.
`lambda` is precisely the per-observation level the softmax already leaves free,
made explicit. The family was about 60 lines and touched nothing else.

Second, and the reason it is not shipped: **the prize was supposed to be the
mixing, and the mixing is not what decides it.** Head to head under hard rules on
the same data and seed, against the Polya-Gamma route the package already had:

| | seconds | ESS | ESS/s |
|---|---|---|---|
| K = 3, multinomial-Poisson | 4.62 | 170.6 | 36.9 |
| K = 3, Polya-Gamma | 1.19 | 134.9 | **113.4** |
| K = 6, multinomial-Poisson | 5.41 | 107.5 | 19.9 |
| K = 6, Polya-Gamma | 1.59 | 156.3 | **98.6** |

The prediction about mixing was right at K = 3 -- one scalar latent per
observation does couple to the predictor less tightly than one Polya-Gamma
variable per category, 170.6 against 134.9 -- and it is worth nothing, because
Polya-Gamma is 3.9x faster. At K = 6 the Poisson route loses on both counts.

**The mechanism, and the general lesson: the quadratic form beats the exponential
form, and it is not close.** A quadratic target is reached exactly in one Fisher
scoring step, the fitted normal *is* the conditional posterior, and -- the part
that dominates -- the per-unit kernel is pure arithmetic. The exponential form
still iterates Newton on scalars and carries one `exp()` per observation per
evaluation. That is the same finding as "what the bounded gates save is the
`exp()`", arrived at from a different direction. Reaching for the exponential
form because it is the shape a Poisson has, when a quadratic rewriting of the
same likelihood exists, is backwards.

The family was deleted rather than kept behind a flag: two routes to one
likelihood with no way to choose between them is worse than one, and the
measurement is here.

**What did come out of it.** The head-to-head forced a re-measurement of the
Polya-Gamma multinomial against the direct family, and it is far better than the
recorded row said: 9.6x in ESS/s under hard rules and 10.1x under soft, against
the 1.6x that had kept it off by default. It is now a default. The old row said
4.2x speed and 0.37x ESS; the new one says 14.5x and 0.66x (hard). The two cannot
be reconciled -- different problem, and ESS measured over the fitted
probabilities here rather than whatever it was then -- so the old row is replaced
rather than averaged with, and the new one is stated with its configuration.

### Multinomial probit: correlated categories, and the sampler that fits the trees in the normalized space

`multinomial()` is a logit, and a multinomial logit cannot express dependence
between the categories at all. The probit version can, and that is the whole
reason for it: the outcome is the largest of `C + 1` latent utilities, and
differencing against a reference category leaves
`W_i ~ MVN(eta_i, Sigma)` with one forest per component and a covariance matrix
to draw.

**The target is exactly quadratic, which is why it is cheap.** Conditional on
`W` and `Sigma` the log density is `-(1/2)(W_i - eta_i)' P (W_i - eta_i)` with
`P = Sigma^{-1}`, so the score in component h is `sum_k P_hk (W_k - eta_k)` and
the information is `P_hh`: `TARGET_QUADRATIC`, the closed-form leaf draw, no
Laplace approximation and no Metropolis ratio. A fit costs about twice a
multinomial logit's on the same data (7.2 s against 2.9 s at n = 800, 50 trees,
500 + 500 draws), which buys the covariance.

**Which sampler.** Xu et al. (2025) give two proposals and compare them against
Kindo, Wang and Pena's (2016) original. Their Algorithms [P1] and [P2] measured
indistinguishable on every figure in the paper -- accuracy, the covariance
estimates, the autocorrelation of tree depth, and the application table to two
decimals -- so [P2] was implemented, being the one with no expansion parameter at
all: draw `W` by a Gibbs sweep of truncated normals, fit the forests, draw the
unnormalized covariance from its inverse Wishart conditional, and rescale to
`trace(Sigma) = C`. What separates both from [KD] is that [KD] fits the trees to
the *unnormalized* utilities, which keeps being rescaled underneath the
stochastic search; their tree depths come out at 6 and 9 against about 2.

Identification is the trace constraint of Burgette and Nordheim (2012) rather
than pinning a diagonal element, which keeps it symmetric in the categories and
makes a two-category fit *exactly* binary probit -- one latent variable, trace
one, `P(S = 1) = Phi(eta)`. Measured, the two agree to a correlation of 0.997 on
the predictor.

**Validated against an independent reimplementation.** The behavior of the
covariance draw looked wrong at first -- see the next paragraph -- so Algorithm
[P2] was written out again from the paper as forty lines of plain R for a
*linear* multinomial probit, with no bartisan involved. It reproduces the C++ to
the second decimal: with the mean fixed at the truth and a true correlation of
0.7, the R reference gives 0.758; with the mean fitted, 0.419, against bartisan's
0.38 on a comparable problem. That is what established the implementation is
right and the surprise is the model's.

**Two things a user has to be told, both found by measurement.**

*The correlation is attenuated, so read its sign and not its magnitude.* A
nonparametric mean absorbs part of the dependence between categories. On a linear
truth at n = 900, true correlations of 0, 0.5 and 0.8 came back as about 0.28,
0.55 and 0.83; the R reference shows the same thing happening as soon as the mean
is estimated rather than known (0.758 -> 0.419 with a *linear* mean of the
correct functional form). The paper's own comparison of samplers turns on the
sign of `sigma_12` for exactly this reason.

*The inverse Wishart degrees of freedom are not exposed, and that is deliberate.*
The obvious knob would be `nu`, and it does the opposite of what it looks like.
With the scale matrix held at the identity, raising `nu` does not pull the
correlations towards zero: `Psi = I` contributes about 1 against a residual
scatter of order `N`, so it is swamped, and all a large `nu` does is make the
draw concentrate on that scatter -- whose correlation the truncation in the
latent draw inflates. Measured against a true correlation of 0.7: `nu = 3` gives
0.38, `nu = 10` gives 0.79, `nu = 50` gives 0.97, `nu = 300` gives 0.996. The R
reference reproduces it (0.997 at `nu = 300`), so it is the prior
parameterization and not the code. Shipping the knob would have shipped a way to
get a confident wrong answer, so `nu = C + 1` -- Imai and van Dyk's (2005) choice,
and the paper's -- is fixed.

**The likelihood has no closed form**, since a category probability is a
`C`-dimensional Gaussian orthant probability. Two different simulators, for two
different jobs. The reported `loglik` uses a *fixed* set of standard normal draws
held by the family, so it is a deterministic function of `eta` and `Sigma` and
the chain sees no Monte Carlo noise -- only a bias that is the same at every
iteration, which is what a convergence diagnostic needs. Predictions use fresh
draws, which is unbiased, and their error is per posterior draw and averages down
over them: a few hundred replicates give an accurate posterior mean from a chain
of a few hundred draws, which is why `replicates` defaults to 200 rather than to
something that looks more careful.

`augment` does not apply: the latent variables are the model rather than a
rewriting of it, so there is nothing to turn off, and `augment = "mnp"` is
refused as an unknown name.

**It is a link on `multinomial()`, not a family of its own.** The two are
different enough inside the engine to be separate `Family` classes -- the probit
carries a covariance matrix, and under the same coding has one fewer forest --
but they are one model to the caller, and the package's convention is that the
link is an argument (`binomial()`, `ordinal()`). So `multinomial("probit")` it
is, and `family_label()` reports it as a multinomial with a probit link rather
than by the engine's name for it.

### The multinomial probit is reference-sensitive, and the symmetric fix is not worth it

Reference coding is what Murray's symmetric multinomial removed for the logit, so
the question is whether the probit needs the same treatment. It is
reference-sensitive, measurably, and about half as much as the reference-coded
logit.

Measured at n = 700, three categories, 50 trees, 500 + 500 draws, on the fitted
probabilities, which are the identified quantity. The spread across the three
possible reference categories only means something next to the spread across
three *seeds* at a fixed reference, since anything the two share is Monte Carlo
noise:

| | across references (max / mean) | across seeds (max / mean) | ratio of means |
|---|---|---|---|
| `multinomial("probit")` | 0.141 / 0.0263 | 0.092 / 0.0140 | **1.9x** |
| `multinomial(reference = )` | 0.156 / 0.0257 | 0.045 / 0.0063 | **4.1x** |

So both exceed noise and the probit is the better behaved of the two. Held-out
error barely moves either way -- RMSE against the truth was 0.042, 0.047, 0.046
across the probit's three references, and the symmetric logit's 0.043 sits inside
the reference-coded logit's 0.040 to 0.047 -- which is consistent with Xu et al.
reporting their Table 2 accuracies as stable to two decimals across reference
levels. Accuracy summaries are stable; individual predicted probabilities move by
up to 0.14.

**A symmetric parameterization is possible and was not built.** It would fit
`C + 1` forests for the raw utilities with no differencing and no zero threshold
-- the latent draw becomes a Gibbs sweep truncated only by `Z_l <= Z_winner`,
which is fully symmetric in the categories -- and Murray's argument carries over:
a proper leaf prior gives a proper posterior and the identified quantities are
recovered. The reason not to is what it does to the covariance. The probabilities
depend on `Omega` only through `var(Z_k - Z_l) = Omega_kk + Omega_ll -
2 Omega_kl`, which is invariant to `Omega -> Omega + a 1' + 1 a'` for any vector
`a`, so `Omega` would carry `C + 1` unidentified directions plus the overall
scale, against the *one* unidentified direction the symmetric logit has. The
reported covariance -- the thing the probit link exists to give you -- would be
uninterpretable without differencing it back, and its posterior would wander
freely in those directions. Buying that to halve a sensitivity already at 1.9x
noise is the wrong trade. Recorded rather than done, and the measurement is here
if it ever looks worth revisiting.

### Soft rules: the gate, and the bandwidth move

Soft rules are charged for in two places, and it took instrumentation to see which.

**Support inflation.** A logistic gate never saturates: dropping a weight needs it below `WEIGHT_TOL` = 1e-10, which needs the observation about 23 bandwidths from the cutpoint — 2.3 units on a predictor scaled to [0, 1], so it never happens.

| | obs per leaf (of 1000) | work per tree | seconds |
|---|---|---|---|
| hard rules | 379 | 1000 | 0.23 |
| soft, bandwidth 0.02 | 768 | 1983 | 0.60 |
| soft, bandwidth 0.10 | 1000 | 2656 | 0.57 |
| soft, bandwidth 0.30 | 1000 | 2707 | 0.55 |

Even at a bandwidth of 0.02, three quarters of the sample reaches every leaf. With the bandwidth held fixed the soft/hard time ratio is 2.5x, exactly the work ratio: nothing else is going on.

**The bandwidth move.** It is a Metropolis step per tree per sweep, and each attempt rebuilds every membership weight in the tree. Instrumented at 50,000 attempts with 58% rejected, it was **48% of a soft-rule fit**, of which `rebuild_support` was 31% on the accepted path and another 17% rolling back a rejection.

**The bounded gates.** `gate = "smoothstep"` is the Beta(2, 2) CDF, `t^2 (3 - 2t)`, on a bounded interval; `"smootherstep"` is Beta(3, 3). Past the interval the answer is exactly zero or one, so the observation takes one side outright, the subtree on the other side is never visited, and the gate is a polynomial rather than an `exp()`. The half-width is `pi * sqrt((2a + 1) / 3)` times `bandwidth` for the Beta(a, a) gate — 4.06 and 4.80 — which equates the gates' standard deviations, so `bandwidth` means the same amount of smoothing whichever is chosen.

Five candidates were considered, written as Beta(a, a) CDFs so the trade is explicit — more derivatives means a wider kernel for the same smoothing, and the kernel width is what is paid for:

| gate | kernel | derivatives | support, in sd |
|---|---|---|---|
| linear | Beta(1, 1) | 0 | 1.73 |
| smoothstep | Beta(2, 2) | 1 | 2.24 |
| raised cosine | — | 1 | 2.30 |
| smootherstep | Beta(3, 3) | 2 | 2.65 |
| logistic | — | all | unbounded |

The raised cosine is **strictly dominated** and was rejected without being built: smoothstep's smoothness class, a wider kernel, and a `cos()` to pay for it. The other three were built and measured over three test functions and six replicates each — all at 0.68 to 0.71 times the logistic's time, with accuracy indistinguishable from it and from each other (paired t of −0.63, +0.30, −0.09).

**That the narrowest gate is no faster than the widest was the surprise, and it falsifies the reasoning that motivated the bounded gates in the first place.** Holding the bandwidth fixed so support width is the only thing varying:

| gate | bw = 0.03 | bw = 0.10 | bw = 0.30 |
|---|---|---|---|
| logistic | 0.92 | 0.81 | 0.81 |
| linear | 0.53 | 0.61 | 0.54 |
| smoothstep | 0.60 | 0.68 | 0.56 |
| smootherstep | 0.59 | 0.62 | 0.56 |

At a bandwidth of 0.30 a Beta(2, 2) gate has half-width 1.22, wider than the whole unit interval, so it truncates **nothing** — and it is still 1.45x faster than the logistic. **What a bounded gate saves is the `exp()`, not the work on the far side of the cutpoint.** The support-inflation measurement above is still correct; it simply was not the thing costing the time.

So no gate improves on smoothstep for a reason worth having. `"smootherstep"` was kept because it is free and gives a twice-differentiable fit, which matters if the fit is going to be differentiated. `"linear"` was **dropped**: fastest by about a twentieth, within noise, and it gives up a differentiable fit, which is the thing soft rules exist to provide.

**Fixing the bandwidth is faster and more accurate, and is still not the default.** Holding it at its prior mean is 2.4x faster again and, on the first function tried with three seeds, more accurate. Six seeds on three functions said otherwise:

| | Friedman | smooth | three-step |
|---|---|---|---|
| logistic, drawn | 0.2081 | 0.1772 | 0.1874 |
| logistic, fixed | 0.1930 | 0.1635 | **0.4186** |
| smoothstep, drawn | 0.2067 | 0.1697 | 0.1901 |
| smoothstep, fixed | 0.1852 | 0.1568 | **0.4624** |

On a function with jumps, fixing the bandwidth more than doubles the error — which is the point of the move: it is what lets the rules sharpen towards hard ones where the truth is sharp. Paired over all 18 fits, fixed is worse (t = +2.4 for both gates). The default stays, and `bandwidth_every` exposes the trade as a stride rather than an on/off switch. This is the clearest case in the file of a speed win that would have cost real accuracy, and of three seeds on one data-generating process pointing the wrong way.

**Two exact improvements to the move itself.**

*A tree with no splits has no gate at all*, so the bandwidth does not enter the likelihood and its full conditional is exactly the prior. Drawing from it directly is an exact Gibbs step and free, where the Metropolis version rebuilt every membership weight in the tree and evaluated the whole likelihood to decide a move that cannot change the fit. Verified where the answer is known in closed form: with the branching probability set to zero so every tree stays a single leaf, 60,000 draws gave mean 0.10020 and sd 0.10004 against an Exp(0.1) prior, KS p = 0.80, and quantiles matching to four decimals. It also mixes better — independent draws rather than a random walk.

*A rejected proposal is rolled back from a snapshot* rather than by evaluating every gate a second time. `SupportStore` holds the whole subtree's memberships in one flat buffer, in pre-order, reused across attempts. Verified bit-identical to the rebuild it replaces, in isolation, across all 20 configurations of the equivalence harness. This took the rollback from 17% of a soft fit to 0.7%.

### Tree bookkeeping: three items, and the one that mattered was not the documented one

The previous round's analysis said the gap to dbarts was the *number of passes* a move makes — six for a birth against about two — and named two fusions. Both were pursued; a third thing found on the way was worth four times as much as either.

**What the passes actually cost, measured rather than argued.** Calling `make_base`/`make_base_children` twice per move instead of once costs 6.1% of a hard fit and 4.2% of a soft one, across all three move types. Since birth is a third of those calls, dropping `make_base` from the birth move — the documented fusion — is worth about **2%**, upper bound. And it would trade a sequential read of `base[k]` for a gather of `eta[idx[k]]`, giving some of that back. **Not done**, and the measurement is the reason rather than an argument.

**Folding `make_child_weights` into `split_support`** — the other documented fusion — *was* done, because it removes a gate evaluation per observation rather than a buffer pass, and the gate is the expensive primitive for soft rules. Dividing a node's support and recording what each side got are now one pass where they were two. Bit-identical.

**The node pool, which was not on the list.** Profiling `birth_leaves` at 10% of a hard fit and 8% of a soft one made it obvious: it allocated two `Node` objects per birth attempt, and 64% of births are rejected, so most of that was `new` immediately followed by `delete`. Nodes are now recycled through a free list on the `Tree`, and a recycled node's index and weight vectors keep their capacity, so a support does not have to grow them again either. Every field a child starts with goes through one `Node::init_as_child()`, so a recycled node and a fresh one cannot drift apart. `change_rule` was also allocating two local vectors per proposal; those became reusable `Context` buffers.

Together, bit-identical, and worth **5% on both hard and soft**.

**And then the thing that dwarfed all of it, which had been considered and rejected on bad reasoning.** `split_support` filled the children's index and weight vectors with `push_back`. The earlier note in this file argued against changing it because `resize()` value-initializes the unused tail and that memset might cost more than the capacity check it saves. That was wrong. Sizing the vectors to the parent's support, filling by index and trimming is **20.6% faster on a soft fit** and neutral on a hard one — four vectors' worth of per-element capacity tests, in the innermost loop of the whole sampler, against one memset. More than the two structural fusions put together.

**The membership weights and the fused sums, from the round before, were both smaller than predicted and are worth recording as such.** A hard tree no longer stores weights (every one is exactly 1.0) and the leaf sums no longer materialize a node's predictors into a buffer before reading them back. Both bit-identical; together about 10%, against a predicted factor of several. The prediction assumed memory bandwidth was the constraint. It is not: a node holds a few thousand observations at most, so the buffer is 3–8 KB and never leaves L1. Isolating the fusion with `block_eval = TRUE`, which forces the buffer back, gives 0.477 vs 0.459 hard and 3.284 vs 3.257 soft — about 4%.

**Devirtualizing the leaf sums** (`Concrete<Derived>`, the curiously recurring template pattern, so the loops resolve each family's `logdens_unit` and `score_info_unit` statically) was worth 1.05x to 1.20x, well short of the 4x the per-visit analysis implied. That analysis was right when made and stale by then: it predated the conjugate shortcut, when family evaluations dominated. `generic_accumulate = TRUE` forces the dispatched path and the two are bitwise identical across eight configurations, which is what made a refactor of every hot family safe.

### The bit-identity harness, which is what made all of this safe

Twenty configurations — five families, hard and soft rules, both gates, the MIA path, and each of the `block_eval`, `generic_accumulate` and `exact_quadratic` diagnostic flags — run from a fixed seed and compared element by element against a stored baseline. A change that is supposed to be a pure optimization has to come out **bit-identical**, not close.

That standard earned its keep repeatedly. It caught the membership-weight change producing 1e-15 differences on soft rules, which turned out to be my own ternaries sitting inside additions (`base[k] + (wt ? wt[k]*mu : mu)`) blocking the fused multiply-add the original allowed; hoisting the branches out of the loops fixed it and was faster anyway. It isolated the split-less-tree Gibbs draw as the *only* source of difference when two bandwidth changes landed together, by rebuilding a variant with that one line disabled and confirming all 20 still matched. And it is what allowed the `split_support` rewrite — a change to the innermost loop of the sampler — to be adopted in one step.

### Earlier performance work, for the record

Before any of the shape-exploiting work above, a first pass of ordinary optimization gave **1.6x** aggregated over five families:

1.  **Score and information in one pass.** Fisher scoring always wants both at the same point, so computing them separately doubled the loop and the dispatch. For families that difference the log density it also lets three evaluations do the work of five. About 33% for the numeric families.
2.  **Analytic derivatives for `ordinal()` and `ordbeta()`**, removing numeric differentiation from the two slowest families — `ordinal()` 2.1x. For the ordered beta the two digamma contributions from `log Beta(a, b)` collapse to their difference, because the shapes move in opposite directions.
3.  **Hoisting the eta-free part of each log density** out of the hot path. Terms that do not involve the predictor cancel from every acceptance ratio, so the sampler never needs them, while a reported density adds them back. Removes the factorial term from Poisson, all three log-gamma terms from the negative binomial, the shape normalizer and `log(y)` from gamma, `lgamma(phi)` from ordered beta and the log scale from Gaussian. Poisson 1.6x.
4.  **The binomial family was computing everything twice or more.** `score_info_unit` was not overridden, so the base class called `dlogdens_unit` and `info_unit` separately — four `pnorm` and two `dnorm` calls for one Fisher-scoring evaluation, where the two quantities share every term. Then: only the smaller probit tail needs `pnorm`, the other follows from `log1mexp`; `dnorm` on the log scale is two multiplications, not a function call; a binary response multiplies one of the two log-density terms by zero; and `log(exp(eta))` in the cloglog path is `eta`. Six special-function calls became one. **Probit 8.66 s → 2.63 s.**
5.  **Fisher-scoring tolerance from 1e-3 to 1e-2.** The effect on the acceptance rate is of order the square of the tolerance. What matters for reversibility is that the tolerance is *fixed*, not that it is tight.
6.  **Adaptive bandwidth proposal.** The multiplicative random walk was accepting 55% against the 44% optimal in one dimension, i.e. its steps were too small. Robbins-Monro tuning of the log step with a vanishing gain, run during warmup and frozen before the retained draws, moved acceptance to 0.39 and raised the bandwidth's per-tree effective sample size from a median of 85 to 120 out of 1000, at no cost in time.
7.  **Reusable buffers on `Context`**, so a sweep no longer allocates three vectors per proposal.

`test-derivatives.R` checks every family's analytic score against a central difference of its own log density, which is the test that would catch an error in any of this.

`ordbeta()` resisted all of it and remains the most expensive family. Varying the tree count isolated the nuisance-parameter updates at only 5% of its runtime, so the cost is the per-tree work, and within that the two log-gamma calls every evaluation of its density needs.

### Parallel chains

The chain is the only parallel axis this sampler has. A sweep conditions on the one before it, and the per-move work is far too small to justify synchronizing a within-chain split.

Implemented through `future.apply::future_lapply()` rather than by choosing a backend, so `plan(multisession)`, `plan(multicore)`, a cluster and mirai's `plan(mirai_multisession)` all work and the package takes a Suggests dependency rather than an Imports. `future.seed = TRUE` gives each chain an L'Ecuyer stream, which is why one `set.seed()` reproduces the whole run whatever the backend — verified: sequential and four-worker multicore runs are bitwise identical.

Measured on `lalonde`, 4 chains of 500 draws: 10.14 s sequential, 3.50 s on four workers, 2.9x. Sublinear because the problem is small and the workers have to be started; the fixed cost is per-run, not per-draw.

The fiddly part was the stored forests. `tree_start` indexes into a flat vector of tree records at a position running iteration, then forest, then tree, so pooling chains means shifting every chain's offsets past the total length of those before it. `expect_predictor_invariant()` over a three-chain fit is the test that catches getting that wrong.

## Log: statistical behavior

### Why credible intervals miss nominal coverage

Measured on the Friedman function, n = 250, nominal 95%, averaged over 12 replicate datasets, reporting the mean absolute bias of the posterior mean alongside the mean posterior standard deviation:

| Setting | Coverage | abs bias | post sd | ratio |
|---|---|---|---|---|
| default (k = 2, 50 trees) | 0.945 | 0.261 | 0.338 | 0.77 |
| 4 chains | 0.950 | 0.258 | 0.339 | 0.76 |
| 4x longer chain | 0.951 | 0.262 | 0.343 | 0.76 |
| k = 1 (wider leaf prior) | 0.943 | 0.260 | 0.334 | 0.78 |
| k = 0.5 | 0.949 | 0.264 | 0.344 | 0.77 |
| 200 trees | 0.939 | 0.264 | 0.334 | 0.79 |
| hard rules | 0.966 | 0.290 | 0.414 | 0.70 |

The mechanism is **bias, not mixing**. The posterior mean sits about 0.77 posterior standard deviations from the truth on average, and an interval centered on a biased point estimate loses coverage in proportion to that ratio. Four chains and a fourfold longer chain move coverage by 0.005, which rules out mixing in this regime — unlike He and Hahn's setting (n = 10000, p = 30), where warm-starting moved coverage from 0.74 to 0.96 and mixing clearly was the binding constraint. Undersmoothing through the leaf prior does not help either, because it inflates the posterior spread and the bias in step. Hard rules over-cover by carrying a wider posterior, not by being less biased.

The one lever that did work is getting the mean-variance relation right; see the quasi-likelihood entry.

### Quasi-likelihood (Linero 2026): a clean extension, with a caveat

Linero, "Bayesian Nonparametric Quasi Likelihood", JASA 2026, replaces the log density with Wedderburn's quasi-deviance, so the only distributional assumption is a mean-variance relation `var = phi * V(mu) / w`. The existing `Family` interface covers it: the score is `(w / phi) * (y - mu) * mu'(eta) / V(mu)` and the information `(w / phi) * mu'(eta)^2 / V(mu)`, which is the GLM working weight and slots straight into the three-quantity contract. Prior weights are already Linero's `omega`, and `1 / phi` is a uniform temper of all three quantities.

Two things do change.

- **`compute_eta_free()` must stay empty for such a family.** That machinery holds the part of a genuine log density that does not involve the predictor — exactly the normalizing constant a quasi-likelihood does not have. Relatedly, `phi` cannot be drawn by slice-sampling the objective, because the objective carries no information about it; it has to come from Pearson residuals.
- **The exactness claim weakens.** Those dispersion updates are incompatible with the tree conditional, so the chain has a well-defined stationary distribution rather than being a sampler from the posterior.

The advantage that matters here is calibration, and it is the one lever found that improves coverage. Linero's own comparator in his section 5.2 is precisely this package's `Gamma("log")` with a drawn shape, and the following reproduces his finding. Inverse-gamma data, which shares the gamma's mean-variance relation but has far heavier tails, n = 250, 10 replicates:

| True dispersion | Drawn | Coverage |
|---|---|---|
| 0.50 | 0.27 (54%) | 0.900 |
| 1.00 | 0.38 (38%) | 0.898 |
| 2.00 | 0.45 (22%) | 0.884 |

Drawing the dispersion from the assumed likelihood underestimates it badly, by more as the true dispersion grows, which narrows the intervals and costs coverage. A Bayesian-bootstrap dispersion draw would be a few lines inside the existing `update_aux` and is the cheapest available route to better calibration.

The advantage is *not* robustness to a misspecified variance function; the paper is explicit that it has none, and its Figure 3 shows coverage collapsing to 0.45 when the variance function is wrong.

### Bootstrapping: one chain per replicate

Dirichlet weights work directly, because prior weights act as exact frequency weights, so `weights = n * rexp(n) / sum(rexp(n))` is a Bayesian bootstrap draw with no changes to the package. Nominal 95%, n = 250, 40 replicates, averaged over 8 datasets:

| Construction | Coverage | Width |
|---|---|---|
| ordinary posterior, one chain | 0.927 | 1.24 |
| Bayesian bootstrap, posterior mean per replicate | 0.891 | 1.22 |
| Bayesian bootstrap, one draw per replicate | 0.985 | 1.82 |

Neither bootstrap construction is the model's posterior, and the choice of summary decides which object you get. Summarizing each replicate by its posterior mean gives the weighted likelihood bootstrap — the sampling distribution of the *estimator* — and it covers **worse** than the posterior, because the estimator is biased and bootstrapping does not remove bias. Taking a single draw per replicate pools posterior and bootstrap variability and over-covers with intervals 47% wider. So a single chain per replicate is computationally fine and statistically coherent only under the second reading, as a nonparametric posterior over the population functional in the sense of Lyddon, Holmes and Walker, and it is conservative rather than calibrated.

### The number of trees, and why nobody puts a prior on it

Binary response, n = 400 train and 400 test, 8 replicates, everything else at the defaults:

| Trees | Held-out log score | Test RMSE | Posterior mean sigma_mu | sqrt(m) sigma_mu | Seconds |
|---|---|---|---|---|---|
| 5 | -226.5 | 0.542 | 0.870 | 1.95 | 1.5 |
| 10 | -227.2 | 0.530 | 0.661 | 2.09 | 2.2 |
| 25 | -226.6 | 0.521 | 0.410 | 2.05 | 5.3 |
| 50 | -227.6 | 0.537 | 0.289 | 2.05 | 10.7 |
| 100 | -227.5 | 0.537 | 0.196 | 1.96 | 21.8 |
| 200 | -227.3 | 0.520 | 0.133 | 1.89 | 54.3 |

The whole range of the held-out log score is 1.05 on a total of 227. Paired within replicate, 200 trees beats 50 by 0.31 with a standard deviation of 1.27 across replicates, and 10 trees beats 50 by 0.39 with 2.08. Neither is a signal.

The mechanism is in the last two columns. The default leaf scale is proportional to `1 / sqrt(m)`, so the ensemble's prior standard deviation at any point, `sqrt(m) sigma_mu`, is `3 / k` for every `m` by construction — and the drawn scale keeps the *posterior* version of that quantity between 1.89 and 2.09 while `sigma_mu` itself falls by a factor of 6.5. The number of trees is very nearly not identified: it changes the richness of the approximation, not the model the prior induces on `f`.

That is the answer to why no BART package does it. A prior on `m` would be integrating over a direction the likelihood barely distinguishes, at up to 40 times the cost, and it would need reversible jump over whole ensembles — the acceptance ratio wants the marginal likelihood of an added tree, which is available in closed form only in the conjugate case this package exists to escape. Meanwhile the thing worth adapting *is* already adapted: the half-Cauchy draw of `sigma_mu` adjusts the ensemble's amplitude, and the Dirichlet sparsity prior adjusts which predictors it spends on.

Caveat: the test function here is smooth and five-dimensional. A rougher target should reward more trees.

### Separation runs the leaf scale away

`sigma_mu` is drawn under a half-Cauchy prior, which has no upper bound. Where the predictors nearly separate a binary response the likelihood rewards an unbounded predictor, and that prior is not enough to hold the scale down. On fully separated data (`y = 1(x1 > 0.5)`, n = 200, 20 trees) the drawn scale averaged 2.79, 9.29, 8.44 and 4.22 over the four quarters of a 1600-draw chain — wandering, not settling, which is what a barely proper posterior looks like — and the additive predictor reached 110. Pinning the scale with `update_sigma_mu = FALSE` brought the maximum predictor to 6.9 at the default value and 4.1 at 0.15. On the same predictors with a non-separable response the drawn scale sat at 0.32 to 0.42 and the predictor at 4.5.

`bartisan()` now warns when the posterior mean of the leaf scale settles more than five times above its prior median. The threshold has room: a genuinely strong signal needs about twice the default scale, and five times corresponds to an ensemble prior standard deviation of 7.5 on the log-odds scale, which is essentially never a real signal.

### Firth-type penalization: wrong direction for the usual bias

Firth's penalty is `+ (1/2) log det I(beta)`, the Jeffreys prior, and it removes the `O(1/n)` bias of the logistic MLE — a bias *away* from zero, which is why Firth's estimator is a shrinkage estimator and the standard remedy for separation.

BART's bias points the other way. On 12 replicates of a binary problem (n = 400, 50 trees, defaults), regressing the error of the posterior mean predictor on the truth gives a slope of **−0.235** with intercept 0.218: the fitted predictor is pulled toward the null fit, by +0.24 in the bottom quintile of the truth and −0.39 in the top. That is shrinkage bias, and adding a shrinkage penalty makes it worse.

There is a genuine information-dependent component — adding `1 / [p(1-p)]` to that regression gives a coefficient of −0.066 with `t = −15.4` — but it lifts `R^2` from 0.213 to only 0.250 against the 0.213 the truth alone explains, and its sign again adds to the shrinkage. So the Firth-shaped part of the bias is both small and pointing the wrong way.

The separation problem Firth was built for *does* appear here, but not in the leaf values: it appears in the drawn leaf scale, which the prior on the leaf values cannot see. The remedy is on that parameter, not a penalty on the leaves. Two further obstacles: `log det I` is over all leaves of all trees at once, so it does not decompose into the leafwise term the Laplace proposal needs; and it would change the target rather than the proposal, so the exactness claim would have to be restated.

### Symmetric multinomial, per Murray (2021)

The reference category was doing real damage. On a three-category problem with a distinct nonlinear surface per category, the mean absolute error of the fitted probabilities was 0.043, 0.058 and 0.048 depending on which category was made the reference: a 35% spread produced by a choice the analyst has no basis for making.

Murray's fix is to stop making it. Fit one forest per category, leave the model unidentified — adding any function of x to every category's forest leaves the probabilities alone — and rely on the proper leaf prior to keep the posterior proper.

**This is a change of parameterization, not of method.** Murray needs a conditionally conjugate leaf prior, so he builds one: a mixture of generalized inverse Gaussians, reached through a gamma data augmentation with one latent variable per covariate value. None of that is needed here, because the Laplace proposal never integrates the leaf out. What transfers is only the redundant parameterization and its prior calibration:

- `H` goes from `num_cat - 1` to `num_cat`; the softmax runs over all of them instead of carrying an implicit zero. The score is `1(y = j) - pi_j` and the information `pi_j (1 - pi_j)` in both cases.
- The leaf scale is divided by `sqrt(2)`. Each log-odds contrast is now a difference of two forests, so its prior variance doubles; Murray's section 4.3 makes the same adjustment, which is where his default `a_0 = 3.5 / sqrt(2)` comes from.
- The intercept becomes the centered log category proportions rather than log odds against the reference.

Worth recording from his supplement (S.1.1): for two categories the two parameterizations induce *exactly* the same prior on the identified function, because the difference of two independent ensembles of `m` trees with a symmetric leaf prior is an ensemble of `2m` trees. For more than two they differ: the symmetric prior makes contrasts against a common category correlated at 1/2 and exchangeable over categories, where reference coding makes them independent and asymmetric.

**Measured**, 12 replicates, n = 400 train and test, 50 trees, 800 iterations, reference coding averaged over every choice of reference:

| | 3 categories | | 5 categories | |
|---|---|---|---|---|
| | symmetric | reference | symmetric | reference |
| RMSE of fitted probabilities | **0.0779** | 0.0829 | **0.0562** | 0.0604 |
| log loss | **0.9074** | 0.9140 | **1.4430** | 1.4500 |
| classification error | 0.3921 | 0.3944 | 0.6650 | 0.6695 |
| coverage of 90% intervals | **0.832** | 0.790 | **0.840** | 0.791 |
| interval width | 0.2096 | 0.2081 | 0.1588 | 0.1569 |
| effective sample size (of 400) | **267** | 206 | **284** | 260 |

Paired within replicate, the RMSE difference is −0.0051 (sd 0.0042, t = −4.2) at three categories and −0.0042 (sd 0.0028, t = −5.2) at five. Against the *best* single choice of reference the symmetric coding is ahead by 0.0014 at three and behind by 0.0028 at five — so it is not that it beats every reference, it is that it lands near the best one without having to know which that is. Within a replicate the best and worst reference differ by 10% of RMSE at three categories and **37% at five**. The point is to delete that dial.

Two unexpected results. The symmetric coding **mixes better** (ESS 267 against 206), presumably because each forest's level is anchored by its prior rather than every forest having to move together with a fixed reference. And its intervals are better calibrated by 4–5 points, though both under-cover.

Cost is exactly the extra forest: 1.509x at three categories against a predicted `c/(c-1)` of 1.500, 1.235 against 1.250 at five, 1.104 against 1.143 at eight. So 50% at three categories, falling away as categories grow. Reference coding is kept as `multinomial(reference = )`, which is the right choice when a particular contrast is the estimand.

### Ordinal models with many thresholds

`rms::orm()` fits proportional-odds models with thousands of intercepts by exploiting sparsity in the information matrix, and the same structure applies here. **Cutpoint k enters the likelihood only through categories k and k + 1**, being the upper limit of one and the lower limit of the other. Summing a cutpoint's slice-sampler target over just those two groups makes a sweep over every cutpoint cost O(n) in total, because the group sizes add to n twice over, instead of O(n * num_cat).

Implemented, with observation indices grouped by category once at construction. Measured at K = n on a continuous response ranked into n categories, which is what `orm()` does: the cutpoint block fell to **0.4% of runtime at n = K = 1600** (0.07 s of 16.1 s), and the whole fit runs to n = K = 3200 in under two minutes with correlation 0.99 against the truth.

Three cautions found along the way.

- **An earlier reading of this was wrong.** Comparing K = 4 against K = n at the same n suggested the cutpoints were 77% of runtime. They were not: a finely graded response supports deeper trees, so that comparison confounds cutpoint cost with tree cost. Direct instrumentation settled it. The superlinear growth in n at K = n is tree work.
- **The computational barrier is gone; the mixing one is not.** Adjacent cutpoints are tightly coupled, and updating them one at a time gives a median effective sample size of only 7 to 21 per 100 draws. The regression function itself mixes fine, so this matters for inference on the thresholds — which for a continuous response *are* the baseline distribution function. The natural fix is again `orm()`'s: the cutpoints' information matrix is **tridiagonal**, so a joint Newton step and a joint Gaussian proposal both cost O(K) by the Thomas algorithm. The obstacle is the ordering constraint, which a Gaussian proposal does not respect; the usual reparameterization to log-gaps destroys the banded structure.
- Storage grows as `num_draws * K`. At n = K = 3200 with 1000 draws the cutpoint matrix alone is 26 MB.

## Log: features and interfaces

### Ordinal: a third link, and the chart the cutpoints are reported in

**The complementary log-log link** was added to match what `stochtree` offers, and
it turned out to bring its own augmentation. The cumulative cloglog model *is* the
discrete proportional hazards model, which is not just an analogy: writing the
survivor as

```
P(Y > k) = exp(-exp(c_k - eta)) = exp(-L_k * exp(-eta)),   L_k = exp(c_k)
```

says exactly that a latent waiting time with exponential rate `exp(-eta)` has
passed `L_k`. So the category is which interval between the transformed cutpoints
the time lands in, and conditional on it the log density in eta is
`-eta - T*exp(-eta)` — the **exponential** form, the same shape as the gamma
family, rather than the quadratic one the probit and logit augmentations reach.

Measured on n = 1500, 50 trees, 400 + 400 draws:

| | direct | augmented | speed | ESS | ESS/s |
|---|---|---|---|---|---|
| hard | 18.75 s | 3.65 s | 5.1x | 55 → 58 | 5.2x |
| soft | 35.27 s | 12.67 s | 2.8x | 93 → 53 | 1.6x |

Hard rules are where it pays, because that is where the exponential form applies;
under soft rules what is left is one `exp()` per observation instead of a
difference of two extreme-value distribution functions, which is a real but small
gain bought with worse mixing. On by default for both, since both are positive.

Correctness was checked by generating from the cloglog link and confirming that
the cloglog fit beats the logit and probit ones on the log likelihood
(−1817 against −1846 and −1840) while recovering the cutpoint gaps — 1.15 and 1.96
against a truth of 1.1 and 1.9. That is the sharp test: a wrong density gives a
plausible fit but not the best one on data from its own link.

**The identification changed.** Only the differences `c_k - eta_i` are identified,
so one location has to be pinned, and the sampler pins the first cutpoint at zero.
That is a fine chart to *work* in and a poor one to *report* in: the cutpoints are
then not comparable with `polr()`'s, and the reader has to remember that the
predictor carries an offset. Draws are now recorded in the chart where the
predictor has mean zero over the fitted sample and every cutpoint is free.

The mechanism matters for why this is safe. `Family::report_shift()` returns the
amount, and it is applied to three things together: the recorded predictor, the
recorded cutpoints, and **the recorded leaf values**, at 1/num_trees per tree. A
tree's membership weights sum to one for every observation, so subtracting the
same amount from all of a tree's leaves moves that tree's contribution by exactly
that amount — which means the stored forest replays to the reported predictor and
the predictor invariant still holds (verified at 6.6e-14). The sampler is
untouched, and every identified quantity with it.

Two things found along the way.

- **`polr()` does not actually center its predictor.** It drops the intercept
  column from the design matrix, which leaves the linear predictor's mean at
  whatever it happens to be — measured at −0.486 on the test problem. So "like
  polr" is the *free cutpoints* half of the convention; mean zero is the natural
  way to pin the location when there is no intercept column to drop. The
  cutpoints recovered in that chart match the truth translated into the same
  chart: −1.291, −0.227, 0.565 against −1.277, −0.177, 0.623.

  The relationship to `polr()` is exact, and there is a test for it. On a linear
  truth with n = 3000, where `polr()` is correctly specified:

  | | cutpoints |
  |---|---|
  | `polr`, predictors as given | −1.041, 0.327, 1.864 |
  | `polr`, predictors centered | −0.479, 0.888, 2.426 |
  | `polr$zeta - mean(polr$lp)` | −0.479, 0.888, 2.426 |
  | bartisan | −0.479, 0.904, 2.444 |

  So **bartisan's chart is `polr()`'s chart with the predictors centered**, and
  either centering the predictors or subtracting the mean of the linear predictor
  puts the two side by side. No chart makes them agree automatically, because
  `polr()`'s convention is a property of its design matrix and a forest has no
  columns; centering is the convention that both can be put in.
- **Two categories are exempt, deliberately.** There the model *is* binary
  regression, the single boundary is conventionally folded into the intercept, and
  reporting it as a free cutpoint against a centered predictor would put the same
  fit on a different scale from `binomial()`. There is a test for that
  correspondence and it should keep passing.

### Ordinal predictions on the latent and mean scales

Two prediction types taken from `WeightIt::predict.ordinal_weightit()`, which is
where the conventions come from.

**`type = "mean"`** reports the probabilities weighted by the category labels read
as numbers, so a response with levels `"1"`, `"2"`, `"4"` has a mean between one
and four. The model stays ordinal; only the reporting treats the categories as
numbers, which is the point — it gives a single summary without assuming a
numeric response at the modelling stage. `values` overrides the labels, and labels
that cannot be read as numbers are an error naming the offenders rather than a
guess. Checked against WeightIt: correlation 0.997 on a linear truth, and the mean
of the predictions matched the observed mean of the response read as numbers to
three decimals (2.598 against 2.599). For a binomial response with levels `0` and
`1` it is exactly `type = "response"`, which the tests require.

**`type = "stdlv"`** divides the predictor by the standard deviation of the latent
variable it indexes, `sqrt(var(eta) + var(e))`, where `var(eta)` is over the
fitted sample so the divisor is a property of the model rather than of whatever is
being predicted. `var(e)` comes from the link: 1 for probit, `pi^2/3` for logit,
`pi^2/6` for cloglog.

**The location convention took some settling, and WeightIt is right where it first
looked wrong.** Their `mu` is `-digamma(1)` for cloglog, which is `+gamma`,
whereas the latent error of a cloglog model — a smallest extreme value variate —
has mean `-gamma`. Simulation confirmed the sign: their `stdlv` differs from
`(lp - gamma)/sd` by exactly `2 * gamma / sd`. But it is not an error. They are
shifting the latent so that *its error* is mean zero, which moves the error's mean
into the index and gives `+gamma`; the mirror-image constant for `loglog`
(`digamma(1)`, i.e. `-gamma`) is the same convention applied to the largest
extreme value. That reading is confirmed by their cloglog fit recovering the truth
on data generated from the smallest extreme value latent, so the orientation of
their model matches this one.

Against WeightIt on a linear truth, for all three links: correlation 0.994 to
0.997, standard deviations agreeing to within 0.8%, and the difference a constant
equal to `-mean(lp)/sd` to three decimals — which is the same chart difference as
for the cutpoints, since bartisan centers the predictor and WeightIt drops the
intercept column. A standardized quantity is used for differences, which that
constant leaves alone.

### Posterior predictive draws, and the packages they unlock

`predict()` gave the mean, the link and the density but never drew Y, which was
the largest functional gap against other BART packages. There is now one sampler
covering every family that has one, and with it the standard interfaces:
`posterior_predict()`, `posterior_epred()`, `posterior_linpred()` and `log_lik()`
on the *rstantools* generics, `simulate()` on base R's, `loo()` and `waic()`,
`bayesplot::pp_check()`, `posterior::as_draws()`, and
`performance::model_performance()` / `r2()`.

**The mean is not derived twice.** For the families that have one, the sampler
takes it from `response_scale()` — the function `type = "response"` already uses,
and the one that already knows about a link the engine does not carry natively.
Only the families whose predictive distribution is not "noise around the mean"
are written out separately: the AFT families (the response scale reports a
*median*), the zero-inflated pair (the mean mixes the two components, and the
sampler needs them apart), the ordered beta (two point masses and a density), and
the two categorical families.

**Validated against the C++ log density, which is an independent statement about
the same distribution.** `iterations` accepts repeats, so `rep(1L, R)` gives R
independent draws at one fixed parameter value; the empirical distribution of
those is then compared against `type = "density"` evaluated at the same value.
For a discrete family that is a frequency against a probability at every point of
a grid; for a continuous one it is quantiles against the numerically integrated
density. Every one of 20 configurations agreed:

| Family | Check | Worst discrepancy |
|---|---|---|
| gaussian, `Gamma`, `location_scale` | quantiles vs integrated density | 0.002 of the range |
| poisson, negbin | frequency vs probability | 0.0035 |
| binomial logit/probit/cloglog, and with 5 trials | frequency vs probability | 0.0043 |
| ordinal logit/probit/cloglog | frequency vs probability | 0.0072 |
| multinomial | frequency vs probability | 0.0011 |
| `zi_poisson`, `zi_negbin` | frequency vs probability | 0.0052 |
| `ordbeta` | both point masses and the interior mass | 0.003 |
| AFT weibull/lognormal/loglogistic | log-time quantiles | 0.014 |

**Scale conventions, which are where this kind of thing goes wrong quietly.** A
binomial replicate is a *proportion*, because that is the scale the likelihood was
written on — so binary data come back as 0 and 1 and trial data as a fraction. A
categorical replicate is an integer category index, because a matrix cannot hold a
factor; `simulate()` returns factors instead, since its result is a data frame and
can. An AFT replicate is a time, not a log time, and it is an *event* time: the
predictive distribution does not know about censoring, so `pp_check()` warns that
the comparison is not like for like. `custom_family()` has no sampler at all — a
log density supplies no way to draw from it — and says so.

**One footgun closed by making it an error.** For a binomial response the trials
live in the prior weights, which are not a function of the predictors, so they
cannot be reconstructed for `newdata`. Defaulting to one trial would answer a
question about counts with a plausible 0/1, so `posterior_predict()` refuses and
asks for `weights`, which is what `predict()` already does about an offset.

**Leave-one-out is documented as strained rather than offered as a number.** PSIS
importance sampling needs finite-variance weights, and a forest is flexible enough
that a single observation can dominate the leaves it lands in; high Pareto k is
common rather than exceptional here. The documentation says so and points at the
held-out log score, which the package can compute exactly.

**`loo` is given the chain structure.** The draws are stacked chain by chain, so
`chain_id` is that block structure; without it `relative_eff()` would treat
dependent draws as independent and understate the standard errors.

### marginaleffects reaches the mean and standardized-latent scales

`get_predict()` was mapping three type names and rejecting the rest, so
`type = "mean"` and `type = "stdlv"` — the two ordinal scales added just before —
were unreachable through *marginaleffects* even though `predict()` had them.
Both are one number per observation, so they need nothing but the name. `"probs"`,
`"lp"` and `"lv"` are accepted as aliases, since those are the names the same
quantities go by for the WeightIt classes in *marginaleffects* itself.
`"class"` and `"density"` stay out, and are refused by name: neither is a number
per observation that an average or a contrast could be taken of.

The dots cannot simply be forwarded to `predict()`, because *marginaleffects*
puts arguments of its own there (`mfx`, and whatever the caller passed to the
estimand function) and a stray name would match one of `predict()`'s arguments
partially. Only `values`, `iterations`, `offset`, `weights` and `log` are taken,
by exact name. *marginaleffects* warns that it does not recognize `values`,
which is correct — it is this package's argument — and the value is used anyway.

**A convention worth knowing, found by a test failure that was the test's
fault.** *marginaleffects* centers a posterior at its **median**;
`predict()` reports its **mean**. On the same draws for the same fit those read
2.402 and 2.393, and the difference is the skewness of the posterior rather than a
disagreement. Documented, and the tests now compare like with like.

### The exponential form under soft rules: measured, and the ceiling is 5–10%

This sat in the To Do list on the reasoning that the hard-rule gain was 1.86x for
Poisson and 1.89x for gamma, so recovering it for the default configuration should
be worth something like that. It is not.

The obstruction is real: a soft rule gives observation i the exponent
`exp(s * w_i * mu)`, and a sum of terms with different exponents is not a function
of three numbers, so `exponential_usable()` declines. The proposed fix was to
bucket the weights, since a saturating gate puts most of them at or near 0 and 1 —
`2 + K` numbers for K distinct weights.

The ceiling was measured before writing any of it, by forcing
`exponential_usable()` to return `true` regardless of `soft` in a scratch build.
That gives statistically wrong answers, which is the point: it costs what a
*perfect* bucketing would cost — one exponential per Newton step instead of K, and
no bucket-building pass — so whatever it saves is strictly more than the real
thing could. On the Friedman function at n = 1000, p = 10, 50 trees, 200 + 200
draws, two runs each:

| | general path (shipped) | forced exponential form | gain |
|---|---|---|---|
| poisson soft | 2.77 s | 2.63 s | **1.05x** |
| gamma soft | 4.28 s | 3.87 s | **1.10x** |
| poisson hard | 0.69 s | 0.69 s | (already on) |
| gamma hard | 1.21 s | 1.22 s | (already on) |

So the whole optimization is worth at most 5% on Poisson and 10% on gamma, and
the achievable version is worth less. The reason the hard-rule case gained 1.86x
and this does not is that the shape shortcut removes the *repeated* passes of the
Newton iteration, not the first one; under soft rules every observation is in
every node's support, so that first pass is O(n) whatever happens and the repeats
turn out to be a small share of a soft fit. The decomposition in "Where the time
goes now" says the same thing from the other side: under soft rules the bandwidth
move is 46.5% of the fit and the leaf refresh is 9.2%, and the exponential form
does not touch the first of those at all.

Dropped from the To Do list. The premise ("few distinct weights") was never the
weak part; the expected value was.

### The zero-inflated and multinomial exponential route: still open, and what it is for

Murray (2021) reaches the multinomial through the multinomial-Poisson
transformation: a multinomial with total `n` and probabilities `softmax(eta)` is
the conditional law of independent Poissons with rates `exp(eta_j)` given their
sum, so introducing a gamma latent for that sum makes the categories independent
Poissons — and a Poisson forest has the exponential form.

**The prize is the mixing, not the speed.** The Polya-Gamma route that ships is
4.2x faster and mixes at 0.37x, which is why `multinomial()` is augmented only
when named rather than by default. Murray's route replaces one Polya-Gamma draw
per observation *per category per forest* with one gamma draw per observation, and
a single scalar latent should couple to the predictor less tightly than one latent
per category does — so the plausible gain is on the 0.37x, which is the number
that is actually holding the augmentation back.

**Why it is a larger change than the negative binomial's was.** That one was
local: same single additive predictor, same output, one gamma draw and a different
`logdens`. This one changes the identification. Independent Poissons with rates
`exp(eta_j)` carry one parameter more than the multinomial does — the total — and
the gamma latent absorbs it, so the per-category predictors are no longer
constrained the way the symmetric and reference codings constrain them now. That
reaches `category_probs()`, the recorded output, the coding option, and the
`before_forest()` machinery.

**The zero-inflated case needs a second augmentation before it needs this one.**
Its `y > 0` term is already the exponential form in the count predictor. The
`y = 0` term is not: it is `log_sum_exp(log_expit(eta_2), log1m_expit(eta_2) +
log_p0(eta_1, theta))`, a mixture of the two components rather than either one. A
Bernoulli indicator for "structural zero" splits that, after which the count part
is a clean Poisson (exponential form) and the inflation part is a binomial
logistic (Polya-Gamma, worth 2.1x in ESS/s). The gain would be real but smaller
than the multinomial's, because a zero-inflated fit divides its time across two
forests and only one of them gains.

### The standardized latent variable, extended to binary responses

`type = "stdlv"` was ordinal-only, on the reasoning that the latent variable is
what an ordinal model cuts. A binary response is the same construction with one
threshold, so it now works there too, for the probit, logit and complementary
log-log links -- the ones whose error distribution has a name and therefore a
variance to divide by. A `cauchit` fit is refused, correctly: its error has no
variance.

**The complementary log-log error enters the two families with opposite signs,
and getting that wrong would have been invisible.** A normal or logistic error is
symmetric, so it does not matter whether `e` or `-e` is the thing added to the
index. A smallest extreme value error is not symmetric. The ordinal model here
writes `P(Y <= k) = G(c_k - eta)`, which is `P(eta + e <= c_k)`, so its additive
error has mean `-gamma`. The binomial model writes `P(Y = 1) = G(eta)`, which is
`P(e <= eta)`, so its latent is `eta - e` and the additive error has mean
`+gamma`. Reusing the ordinal constant would have put the binary answer off by
`2*gamma/sd` -- about 0.80 on the fit this was measured on, against a quantity
whose own standard deviation was 0.42.

Checked by deriving the constant outside the code and comparing:
`mean(stdlv)` came out `-0.24686` against a predicted `(mean(eta) - gamma)/sd` of
`-0.24674`, where the ordinal convention would have given `+0.55676`. The scale
was checked separately for all three links against `sd(eta)/sqrt(var(eta) +
var(e))`: agreement to 0.3%.

**A corollary worth recording:** a two-category ordinal complementary log-log fit
is *not* the binomial complementary log-log model. For a symmetric error,
`1 - F(c - eta) = F(eta - c)` and the two coincide; `G` is not symmetric, so they
do not. Measured, the two fits correlate at 0.99 but differ in scale by 18%,
where the probit and logit pairs agree up to Monte Carlo error. The
two-category collapse the ordinal identification chart promises is a statement
about the logit and probit links.

### The residual-scale prior was ignoring the weights

`residual_scale()` anchors the prior on the Gaussian `sigma` with the residual
standard deviation of a linear fit, which is a much better anchor than the
marginal spread. It was fitting that line unweighted, so a weighted analysis
anchored its prior on the wrong observations -- and the failure is quiet, because
a prior scale that is too large produces a fit that merely looks under-confident.

Now weighted, with one decision in it: only the relative sizes of the weights can
matter, since a residual variance is per observation, so the weights are
normalized to average one *over the rows they keep*. Normalizing over all rows
instead would make a zero weight shrink every scale rather than drop a row --
measured, that put the scale a factor of `sqrt(2)` out on a half-zeroed sample.
With the normalization over kept rows, zeroing out a noisy half reproduces the
clean half's scale exactly, and a constant weight reproduces the unweighted
answer exactly, which is what keeps every existing fit unchanged.

### Missing values are kept by default, and the old default never worked

`na.action` now defaults to `na.pass`, so missing predictors reach the splitting
rules instead of taking their rows with them. MIA has been implemented for a
while; what changed is that it is on unless asked otherwise, which is the right
way round — a tree can do something better with a missing value than either
imputing it or discarding the row, and that is most of the reason to use one.

**A bug surfaced by making the change, and it is the more interesting half.** The
default in the signature had never been used. `model.frame()` is called through a
call rebuilt from `match.call()`, which records only what the caller actually
wrote, so an argument left at its default is *absent* from the reconstructed call
and `model.frame()` falls back on `getOption("na.action")` — usually `na.omit`.
Setting the formal's default to `na.pass` therefore changed nothing at all until
the value was injected into the call explicitly. Any future default for that
argument would have been swallowed the same way, silently. `NULL` is still passed
through as "whatever the session option says", which is what `lm()` and `glm()` do
with it.

### marginaleffects: the draws are the interface

A fitted forest has no coefficient vector and no variance-covariance matrix, so
the delta method `marginaleffects` uses for a frequentist model has nothing to
work with. It has something better: every estimand is computed by pushing all the
draws through the same transformation and summarizing at the end, so an interval
is a posterior quantile and a nonlinear estimand needs no approximation. That is
the same path the package takes for `brms` and `rstanarm`.

What the integration needed, in the order the obstacles appeared:

1. **The base generics.** `formula()`, `terms()`, `model.frame()`, `nobs()` and
   `family()`, which is what `insight` falls back on for a class it does not know.
   The model frame is now retained in the fit, as `glm()` retains it, because a
   counterfactual grid is built from the data the model saw.
2. **`get_predict()` with the draws attached**, as an observations-by-draws matrix
   in `attr(, "posterior_draws")`. For a categorical family the rows are stacked
   category by category and the draws matrix has to be stacked the same way; the
   tests check that the point estimate equals the row mean of the draws, because a
   transposed or mis-stacked attribute produces plausible nonsense rather than an
   error.
3. **Stubs for the frequentist path**: `get_coef()` returning an empty vector,
   `set_coef()` returning the model, `get_vcov()` returning `NULL`. Returning
   rather than erroring is what routes it to the draws.
4. **Registering the class.** `marginaleffects` validates a fit's class against a
   list of the ones it knows *before* any of those methods is reached, and exposes
   `options(marginaleffects_model_classes = )` as the way an outside package adds
   its own. That is done in `.onLoad()`. Without it every call fails with "Models
   of class bartisan are not supported", which is what it did at first.

**A trap worth the documentation it got.** Slopes are numerical derivatives, and
the default `x_transform = "quantile"` maps each predictor through its empirical
distribution function — a *step function*. So the fit is a step function of the
original predictor whatever the decision rules are, and the difference quotient
diverges as the step shrinks. On a surface whose average derivative is zero:

| step | `x_transform = "quantile"` | `x_transform = "range"` |
|---|---|---|
| 1e-4 | −4.79 | −0.28 |
| 1e-2 | −0.40 | −0.30 |
| 5e-2 | −0.29 | −0.25 |

`"range"` is stable across a 500-fold change in the step; the default is not. This
is documented in `?bartisan-marginaleffects` rather than fixed, because the
quantile transform is the default for a good reason — it makes the cutpoint prior
invariant to monotone reparameterization — and predictions and comparisons, which
evaluate the fit at two points a substantive distance apart, are unaffected.

A family with several additive predictors is refused on the `"link"` scale, where
there is no single quantity to be talking about, and works on the response and
probability scales, where there is.

### Random intercepts: implemented, and what the leaf machinery gave for free

`(1 | group)` in the formula, in lme4's notation, parsed with `reformulas`
(`findbars`, `nobars`, `subbars`). Several grouping factors are allowed and
`(1 | a/b)` expands to nesting before this code sees it.

**The implementation is almost entirely reuse, which was the argument for doing
it.** A random intercept is a scalar with a Gaussian prior entering the predictor
with weight one for the observations in its level — which is a leaf with the gate
removed. So each level is a `Node` carrying that level's observations, and it is
updated by `update_scalar()`, the function the leaf refresh calls. The only
change to the sampler was factoring that function out of `update_leaf_params()`.
Everything follows:

- the quadratic closed form for a Gaussian or augmented family,
- the exponential form for Poisson and gamma,
- the general Laplace-and-Metropolis path otherwise,
- and `half_cauchy_update_precision_mh()`, which already took a vector of values
  with a common half-Cauchy prior, for the scale `tau`.

One detail worth recording: each term carries its own `Hypers` with
`soft = false`, purely so `exponential_usable()` sees it. A random intercept's
weights really are one whatever the decision rules are, so the exponential form
is available even in a soft-rule fit — and it would not be if the level nodes
pointed at the forest's own hyperparameters.

Verified by recovery, which is the test that matters for a latent quantity:

| | correlation with the truth | tau | truth |
|---|---|---|---|
| gaussian, 60 groups of 10 | 0.94 | 1.07 | 1.2 |
| two grouping factors | 0.97, 0.95 | 1.04, 0.47 | 1.0, 0.6 |
| probit | 0.93 | 1.02 | 1.0 |
| poisson | 0.98 | 0.52 | 0.5 |
| location-scale, mean forest | 0.97 | 0.76 | 0.8 |
| location-scale, log-sd forest | 0.95 | 0.53 | 0.5 |
| zero-inflated Poisson, count | 0.98 | 0.72 | 0.8 |
| zero-inflated Poisson, inflation | 0.75 | 0.41 | 0.5 |

**Every additive predictor gets its own set**, which is what makes the
multi-forest cases above work: a zero-inflated count model has a group effect on
the count part and another on the inflation part, with separate scales. The tests
check that the two sets are matched to the right forests rather than swapped, by
generating different effects for each and requiring each to correlate better with
its own.

**Only random intercepts.** A slope is refused with a message that says what to
do instead — put the variable in the fixed part, where a tree can split on the
group and the variable together and get an interaction of any shape, which is
strictly more general than a linear slope varying by group. The reason for the
restriction is structural rather than incidental: a slope is a scalar multiplying
a covariate, so it is not a leaf with the gate removed and the machinery above
does not apply to it.

Two smaller things the implementation needed:

- **The model frame is built from `subbars(formula)`** so the grouping variables
  are present and get the same missing-value handling as everything else, while
  the design matrix comes from `nobars(formula)` so they are not also predictors.
  That split introduced a bug worth remembering: `.` in a formula expands against
  the columns of whatever is passed as `data`, and the model frame has two extra
  columns — `(weights)` and `(offset)` — so `y ~ .` picked them up as predictors
  and `predict()` then went looking for them in `newdata`. The frame is now
  handed over without them.
- **A level absent at fitting time** gets the prior mean of zero, with a warning,
  which is what `lme4` does with `allow.new.levels`. There is a test that two
  copies of the same row differing only in the group differ by exactly that
  group's intercept, which is what says the unseen one really got zero.

### Correlated random effects across forests: not implemented, and why

The natural next step is to let a group's intercepts on different additive
predictors correlate — for a zero-inflated count model, to let a group that tends
to produce zeros also tend to produce low counts. It is not here, and the
obstacle is specific rather than a matter of effort.

**A joint update is unavailable.** Correlating `b_g^{(1)}, ..., b_g^{(H)}` under a
multivariate normal prior wants a joint H-dimensional Laplace step, and that needs
the family's *mixed* second derivatives `d2 logdens / d eta_j d eta_k`. The
`Family` interface has no such thing: `info(i, eta, h)` is per-predictor, and
nothing in the package has ever needed a cross term. Adding it means a derivation
per multi-forest family — location-scale, both zero-inflated families,
multinomial, ordered beta — plus a two-dimensional central difference for the
numeric fallback. The two-child `Target2` does not help: its cross term is between
two leaves of one forest at one `h`, not between predictors.

**The tractable route is component-wise, and its cost is where the prior lives.**
Updating `b^{(h)}` conditional on the others needs only the h-th likelihood, with
the multivariate normal conditional as its prior — which is Gaussian with a
*nonzero mean* and a reduced variance. That would work, and it needs a prior mean
threaded through the leaf update. The prior appears at **27 places** in
`mcmc.cpp`, and one of them is the exponential-form coefficient recovery, where
the prior is subtracted off to isolate the likelihood's own value, slope and
curvature and then added back. A mistake there is silent and would affect every
fit, not only correlated ones. Defaulting the mean to zero would keep the
bit-identity harness as a safety net, so this is doable — it is simply a larger
and riskier change than the feature justifies on its own.

**There is also no syntax for it.** lme4's `(1 | g)` describes one response; it
has nothing to say about which of several additive predictors a term belongs to,
let alone which pairs correlate. `glmmTMB` and `brms` reach for multivariate
formula syntax, which is a much larger interface change. A control flag would do,
but a flag that silently decides an H-by-H covariance structure is not a good
interface.

**And the identification is weak for H above two.** An H-by-H covariance is
H(H+1)/2 parameters estimated from as many levels as the grouping factor has. For
a five-category multinomial that is 15 parameters, and the intercepts themselves
are already only weakly identified against the forests' own levels.

What is available instead covers most of the practical benefit: independent
intercepts on every predictor, so a group effect on the count part *and* on the
inflation part, each with its own scale. What is missing is only the correlation
between them.

### Random effects: the assessment this was built from

This is the measurement that decided the feature was worth building, kept because
it is also the answer to "should I use this or put the group in as a predictor".

`dbarts::rbart_vi()` and `stan4bart` both add a group-level random intercept:
`eta_i = f(x_i) + b_{g(i)}` with `b_g ~ N(0, tau^2)`. The question is whether this
package needs it, given that a group can already be handed in as a factor
predictor and split on like anything else.

Measured: a smooth fixed part plus a group intercept with `tau = 1`, unit residual
noise, four replicates per cell, reporting RMSE against the true conditional mean.

| groups | per group | bartisan, group as a factor | group ignored | `rbart_vi` |
|---|---|---|---|---|
| 5 | 100 | **0.190** | 0.844 | 0.279 |
| 25 | 20 | **0.306** | 0.924 | 0.341 |
| 100 | 5 | 0.541 | 0.944 | **0.468** |
| 250 | 4 | 0.708 | 0.960 | **0.493** |

So the answer is regime-dependent and the crossover is around a hundred groups.
With **few large groups** the factor predictor is *better* than a random intercept
— the group means are well determined without pooling, and a factor predictor can
additionally interact the group with the covariates, which an intercept cannot.
With **many small groups** the random intercept wins and the margin grows: at 250
groups of four it cuts the error by 30%, which is the partial-pooling regime doing
exactly what it is for.

**One measurement error worth recording**, because it inverted the conclusion.
`rbart_vi`'s `yhat.train` is the *fixed* part only; the random effects are in
`ranef`, and `fitted()` combines them. Scoring `yhat.train` against the truth made
`rbart_vi` look worse than ignoring the group entirely (0.95 against 0.92), which
should have been obviously impossible and was the tell.

**The design, if it gets built.** The reuse is unusually high, and that is the
argument for doing it:

- **Drawing `b_g` is drawing a leaf.** The target is a Gaussian prior times the
  likelihood over a set of observations with unit weight — which is exactly what
  `Target1` handles. Point it at the group's index set with `nullptr` weights and
  `tau` in place of `sigma_mu`, and the whole apparatus follows: the quadratic
  closed form for a Gaussian or augmented family, the exponential form for the
  Poisson and gamma, the general Laplace-and-Metropolis path otherwise. No new
  sampler.
- **Drawing `tau` is drawing `sigma_mu`.** `half_cauchy_update_precision_mh()`
  already takes a vector of values with a common half-Cauchy prior, which is
  exactly the `b` vector.
- **Cost is about one extra tree per sweep.** Each `b_g` update touches its own
  group's observations, so the whole set of them is one pass over the sample.
- What is actually new: a `group` argument (`rbart_vi` calls it `group.by`; the
  lme4 `(1 | g)` formula syntax is a larger job than it looks), storage for the
  `b` and `tau` draws, prediction for a group not seen at fitting time — either
  zero or a draw from `N(0, tau^2)`, and the choice should be explicit — and the
  predictor invariant extended to include the random part.

Built, as described above; the design sketched here is what was implemented, and
the reuse was as high as predicted. The measurement remains the guidance for *when
to use it*: with few large groups the factor route is better, and the crossover is
around a hundred groups.

### Missing predictor data: implemented as MIA

Built with `na.action = na.pass` as the switch: a standard R idiom that reads correctly ("pass the missing values through"), needs no new argument, and leaves `na.omit` behaving exactly as before.

**All three of Twala's rules, not two.** The `na_left` bit alone gives variants A and B. Variant C — split on whether the value is there — is what makes a variable usable when only its *absence* carries signal, and it is not reachable in practice from the other two: isolating the missing values needs a cutpoint below the node's minimum observed value, which the uniform cutpoint prior hits with probability of order 1/n. So `na_rule` is a three-way choice, drawn uniformly, and only for variables that actually have missing values — which keeps the prior, and the sequence of random numbers, identical on complete data.

**Nothing in the sampler changed.** Birth and change both propose a rule from its prior, so the extra factor appears in the prior and the proposal alike and cancels. That claim is tested, not asserted: with the likelihood flattened the target is the tree prior, which says nothing about missing values, so the distribution of tree sizes has to be the same as with complete data. It is, with missing values in every predictor.

**A missing value takes a hard path even through a soft tree.** Its gate is 0 or 1 by the rule. That is correct — there is nothing about being absent to smooth over — and it keeps the two children's weights summing to the parent's, which everything else rests on.

**One bug found, and it was mine.** `expect_predictor_invariant()` failed immediately: the replayed forests disagreed with the recorded predictor on exactly the observations in one leaf. The change move saved the old rule field by field and put it back on rejection, and that list did not know about `na_rule`. So a rejected change left the node with the new missing-value rule and the old variable and cutpoint, and `split_support()` built a partition the stored record did not describe. Fixed, and then the field list was removed entirely: `Node::rule()` and `Node::set_rule()` carry the whole rule in one object, so the next field added cannot be forgotten the same way. The first run of the comparison below was done on the broken build and reported MIA at eight times the error of every alternative, which is what that bug looks like from the outside.

**Measured**, 8 replicates, n = 400 train and test, 30% of `x1` missing, held out.

| Mechanism | Method | RMSE, all rows | RMSE, observed rows | Test rows scored |
|---|---|---|---|---|
| MCAR | MIA | 0.356 | 0.105 | 400 |
| | fill + indicator | 0.354 | 0.100 | 400 |
| | median fill | 0.359 | 0.129 | 400 |
| | complete cases | — | 0.098 | 284 |
| MAR | MIA | 0.365 | 0.099 | 400 |
| | fill + indicator | 0.366 | 0.102 | 400 |
| | median fill | 0.373 | 0.128 | 400 |
| | complete cases | — | 0.094 | 282 |
| Informative | MIA | 0.356 | 0.100 | 400 |
| | fill + indicator | 0.353 | 0.099 | 400 |
| | median fill | 0.362 | 0.131 | 400 |
| | complete cases | — | 0.094 | 284 |

MIA beats plain median imputation everywhere (paired, +0.025 to +0.031 with paired sds of 0.020 to 0.029) and **ties fill-plus-indicator** (within 0.005 in every cell, and 4 of 8 replicates each way on a harder design). Complete-case analysis is marginally best on the rows it can reach, and cannot predict 29% of the test set at all.

That last point is worth stating plainly rather than spinning: on these designs MIA's advantage over the crude alternative is not accuracy. It is that no fill value has to be chosen, the design does not double in width, missingness at prediction time needs no matching indicator columns built by hand, and it is one argument rather than a preprocessing pipeline. **The earlier claim in this file that fill-plus-indicator is "MIA done badly" was wrong, and the measurement is what says so.**

Where MIA is not merely convenient is the informative case in isolation, which the table understates because `fill + indicator` is handed the same information by construction. With the response depending on whether `x1` was recorded and on nothing else about it, the fitted means come out at 1.98 and −0.01 against a truth of 2 and 0, with 29.6 of 30 splitting rules landing on `x1`. Median imputation alone cannot represent that at all.

The caveat to state plainly: MIA targets `E[Y | X observed, pattern]`. That is right for prediction and wrong if the estimand is a regression or causal effect defined on complete data, where MIA's answer depends on the missingness mechanism. `mice` is the documented route for that.

**Deliberately not done.** `predict()` refuses a missing value in a column that had none at fitting time, rather than defaulting it somewhere: no rule in the forest carries an answer, so every one of them would send it the same arbitrary way and the prediction would be quietly meaningless. A column that is constant where observed but sometimes missing is kept rather than dropped, since its missingness is still something a rule can split on.

**Two routes rejected.** *Surrogate splits* (CART's device, `rpart`'s default): a surrogate has to be computed from the association between predictors *within the node*, which means recomputing it for every proposal, and it is a deterministic construction with no place in the parameter vector, so it does not fit a sampler whose whole job is to average over rules. It also assumes missingness is uninformative, which MIA does not need. *Bayesian imputation of X inside the sampler*: the principled version, and a great deal of machinery — it needs a joint model for the predictors — for a gain that only materializes when that model is right.

### Arbitrary links from R

The compiled families accept only links for which the additive predictor is unconstrained, which left out cauchit among others. Two fixes were available and one is much better.

The rejected one: write each new link in C++. Cheap per link, unbounded in total, and no help to anyone who wants a link nobody anticipated.

The chosen one: **compose the caller's link onto the scale the compiled family already works on.** For the five families with a single mean and a conventional link, the compiled code reads the predictor on a known scale — the mean for Gaussian, the log mean for the counts, the log odds for binomial — so a link `g` is honored by mapping `theta = t(eta)` with `t = native_link . g^{-1}` and applying the chain rule. One decorator, no per-family code, and every link `stats::make.link()` knows works immediately, as does a `link-glm` object written by hand.

Two decisions inside that deserve recording.

**The information drops the term in `t''`.** The exact second derivative of the composite is `l''(theta) t'^2 + l'(theta) t''`, and the second term can be negative, which would break the Laplace proposal's curvature. It is the score times `t''`, and the score has expectation zero, so dropping it leaves exactly the *expected* Fisher information of the composite whenever the wrapped family reports the expected information. Since these numbers only build a proposal and the acceptance ratio uses the exact log density, this costs a little efficiency and no correctness. Verified against the closed form for cauchit binomial to 1e-10, through both evaluation paths.

**The nuisance parameters stay with the wrapped family.** `gaussian("log")` must still draw `sigma`, and `sigma` lives on the wrapped family's scale, so `update_aux` transforms the predictors before delegating, and `aux_names`, `aux_values`, `set_aux`, `log_norm_const` and the eta-free part all pass through.

Not offered for `ordinal()`, `multinomial()`, the AFT families, `location_scale()`, the zero-inflated families or `ordbeta()`. In those the link is not a map from a single predictor to a single mean — an ordinal model's cutpoints sit *inside* the cdf, so composing on eta alone cannot change it — and pretending otherwise would be wrong rather than merely unsupported. `custom_family()` is the route for those.

### custom_family(): a likelihood from R

The interface a family has to satisfy is a log density and its first two derivatives with respect to each predictor. Nothing about that requires compiled code, and `score_info_numeric` already produces both derivatives from three evaluations of the log density. So `custom_family()` takes the log density as an R function, with `num_predictors` forests, an optional analytic `derivatives`, and a `start` value in place of the intercept-only fit it cannot compute.

Deliberately excluded: drawing a nuisance parameter (a dispersion has to be fixed inside the closure), a non-numeric response, and a response-scale prediction, since the package cannot know what the mean of an arbitrary density is. `predict(type = "response")` returns the predictors instead, which is stated rather than silent.

Checked against the compiled Poisson: identical log density up to the dropped `lgamma(y + 1)`, score agreeing to 1e-6 and information to 1e-4, which is the accuracy of the central difference at step 1e-4. A two-predictor location-scale density written by hand recovers both surfaces (correlations 0.98 and 0.83).

### The blocked evaluation path, and why it is opt-in

Both R-backed routes pay a fixed cost per call, so calling per observation is hopeless: a 500-observation, 50-tree, 2000-iteration fit makes order 10^9 per-observation family evaluations. The sampler always evaluates a whole leaf at a time, so the fix is to hand over the leaf: one call per leaf per Fisher-scoring step, order 10^6 for the same fit.

The first attempt made that the *only* path. It cost **30% of the runtime for a Gaussian response** (2.58 s against 1.99 s), 9% for Poisson and nothing measurable for binomial. A Gaussian log density is a handful of operations, so materializing a vector of n of them costs more than the density it computes. Specializing the block fill for a single predictor, where the copy is provably redundant, recovered none of it: the cost is the vector, not the copy.

So both paths exist and `Context::blocked` chooses; `Family::wants_block()` is true only for the two R-backed families. Two code paths for one calculation is exactly the kind of duplication that drifts into a bug, so it is pinned by a test rather than by care: `block_eval` forces the blocked path for a compiled family and `test-flexible.R` runs four families both ways from the same seed. The chains are kept short on purpose — over hundreds of iterations a rounding difference will eventually flip an accept/reject and diverge for a reason that is not a defect.

### Later families

Three families beyond Linero's own set, all fitting the existing framework without changes to the sampler:

- **Zero-inflated Poisson and negative binomial** use the multi-forest machinery: the log count mean and the logit structural-zero probability each get a forest, so the excess-zero mechanism is nonparametric rather than a constant. The observed information of a mixture can go negative, so the proposal uses the complete-data expected information instead, which is positive by construction; the Metropolis step absorbs the approximation. Named `zi_poisson()` rather than `zip()` because `zip()` would mask `utils::zip()`.
- **Ordered beta** (Kubinec 2023) needed only a log density: the derivatives come from the inherited central differences. Unlike `ordinal()`, both cutpoints are drawn, because the predictor also enters the beta mean and is therefore identified without pinning one.
- **Conditional density** (`predict(type = "density")`) calls the family's own C++ log density through a new entry point, so each distribution is defined once. This required a `log_norm_const()` hook for the terms the sampler is free to drop but a reported density is not — in practice the binomial coefficient.

### Rank-normalized diagnostics

`fit$rhat` is a data frame with a row per quantity and three columns: rank-normalized folded split R-hat, bulk effective sample size, and tail effective sample size (Vehtari, Gelman, Simpson, Carpenter and Buerkner 2021).

Why each piece is there. **Rank-normalization** replaces the draws by the normal scores of their pooled ranks, which guarantees the finite variance the formulas assume whatever the posterior looks like, and makes the diagnostic invariant to any monotone reparameterization. **Folding** — the same computation applied to the distance from the median — catches chains that agree about the middle and disagree about the spread; the reported R-hat is the larger of the two. **The tail effective sample size** answers a different question from the bulk one: a posterior mean is an average over every draw and converges quickly, while an interval endpoint depends on the few draws out in the tail.

Validation, since a diagnostic that is quietly wrong is worse than none. Split R-hat matches the textbook formula computed by hand to 1e-10 on four cases, two of them cases the diagnostic is supposed to flag. ESS matches the `posterior` package — the authors' own implementation — to within 1.6% on five cases spanning iid, heavy-tailed and autocorrelated chains. And ESS matches the closed form `MN(1-rho)/(1+rho)` for an AR(1) chain to within 10%, on the conservative side, which is the intended behavior of Geyer's initial positive sequence with the monotonicity correction.

**One bug this surfaced.** An ordinal model's first cutpoint is pinned at zero for identifiability, so there is nothing to diagnose — and the code reported an R-hat of `-Inf` with a warning and a fabricated effective sample size of about 6. The variance guard did not catch it because the sample autocovariance of a constant is a rounding error rather than exactly zero. Both now return `NA`, silently, with a test.

One trap worth recording: comparing against `posterior` initially looked like a disagreement on R-hat, including `posterior` reporting 1.00 for four identical monotonically drifting chains. It was the harness — `posterior` was not splitting the input it was handed. Hand computation settled it. **Compare implementations on inputs where you can also work out the answer yourself.**

## Assessed and not adopted

### Multithreading below the chain level

The profile argues against it. The work is spread across many small per-leaf loops rather than a few large ones — at n = 1000 a leaf's support is a few thousand elements, the same order as thread synchronization overhead — and the two structurally parallel axes are unavailable: backfitting is sequential across trees by construction, and the forests of a multi-predictor family each condition on the others' current predictor.

There is also a portability obstacle specific to this package. Every family calls into R's math library (`lgammafn`, `digamma`, `dnorm4`), which is not documented as thread-safe; XBART sidesteps this by using `std::random` and its own numerics and touching no R API inside a thread. Doing the same here would mean replacing R's special functions, and OpenMP on macOS additionally needs `libomp`.

### XBART (He and Hahn 2021)

Their grow-from-root sweep replaces the reversible-jump tree moves with a recursive pass that samples a cutpoint proportional to the marginal likelihood of the resulting split, and gets 20–28x over BART MCMC. Three findings decided against porting it:

- **It is a hard-rule technique.** Its speed comes from presorting each predictor once and maintaining sorted index vectors, so all candidate cutpoints for a variable share a single cumulative-sum pass and each child's statistics follow from the parent's by subtraction. Under soft rules there is no partition to sort: every observation reaches every leaf with a weight that itself depends on the candidate cutpoint. All of it dies, and soft rules are the default.
- **It is not a posterior sampler**, and its authors say so: the grow-from-root step is "not a proper full conditional" and the estimator is "a greedy stochastic approximation". Their only theorem establishes that *a* stationary distribution exists for a modified version, not that it is the BART posterior. Their own tables show 95% intervals covering as little as 0.50.
- **The conjugacy substitution is possible but limited.** A one-step Laplace expansion about the parent's mode gives a criterion structurally identical to theirs with the count and residual sum replaced by the information and score sums, which are still additive and still prefix-summable. But iterating Fisher scoring per candidate would cost a factor of the grid size, so the criterion would have to stay one-step, and it would still only serve the hard-rule path.

What *is* worth taking is their warm start, which does not touch the transition kernel. It is in the To Do list.

### Ultimate Polya-Gamma samplers (Zens, Fruhwirth-Schnatter and Wagner 2024)

The paper's contribution is in two parts, and the package's position on each is
different.

**The representation is already here, where it applies.** Their equation (5)
writes the logistic density as a Polya-Gamma normal scale mixture,
`f(e) = (1/4) E[exp(-w e^2 / 2)]` with `w ~ PG(2, 0)`, so that `w | e ~ PG(2,|e|)`
is a tilted Polya-Gamma draw. That is exactly the identity the ordinal logit
augmentation in this package rests on, reached independently from Polson, Scott
and Windle's Theorem 1 at `a = 1`, `b = 2`. Independent corroboration of a
derivation that had been arrived at here from scratch, and nothing to do.

**The boosting is the paper's real contribution, and it does not port.** They add
two working parameters to the latent utility equation -- a location `gamma` and a
scale `delta` -- and alternate: draw `gamma` from the working prior, shift the
utilities, redraw `gamma` from its conditional *with the coefficients integrated
out*, shift back, then draw the coefficients. That middle step is what makes the
shift free, and it needs the coefficient vector marginalized. From their own
replication code (`Simulations_Logit/algorithms/LOGIT_V2.R`), every iteration
does

    Bn = chol2inv(chol(A0.inv + t(X * omega) %*% X))
    beta.draw = sqrt(delta.star / delta) * bn + t(chol(Bn)) %*% rnorm(P)

-- an explicit `P x P` inverse and a joint block draw of every coefficient. BART
has neither. Its "coefficients" are the leaf values of every tree, a set whose
dimension changes every iteration, and backfitting exists precisely so that this
matrix is never formed. There is a cheap special case -- marginalize only the
*level* of the predictor, which is one direction and so a scalar -- and it is
worth knowing that it exists, but see below for why it is not worth building.

A second obstacle: iMDA needs a latent utility with a threshold, so that the
observed outcomes restrict `gamma` to `[max z_i(y=0), min z_i(y=1)]`. The paper
says this itself about the original Polya-Gamma sampler, and it applies to this
package: the binomial *logit* augmentation here is the marginal Polya-Gamma form,
which has a Gaussian pseudo-response but no utility and no threshold. Only the
probit and ordinal probit augmentations, which are Albert and Chib latent
normals, could carry the move at all.

**And the pathology it fixes is not present.** This is the part that decided it.
UPG is aimed at a level that has to travel to its posterior region and then random
walk there in tiny steps, which is what happens when the intercept carries a
near-flat prior -- theirs is `N(0, 100)`. Measured here at n = 2000, p = 5, 50
trees, 1000 + 1000 draws, with ESS taken on the level of the predictor
(`rowMeans(eta)` per draw):

| positives | logit augmented | probit augmented | logit direct | probit direct |
|---|---|---|---|---|
| 1019 (51%) | 862 | 354 | 1000 | 1000 |
| 113 (5.7%) | 143 | 100 | 245 | — |
| 27 (1.4%) | 75 | 24 | 28 | 17 |
| 4 (0.2%) | 7.6 | 10 | — | — |

So mixing does collapse with imbalance. Three further measurements say it is not
UPG's problem:

- **Displacing the anchor changes nothing.** The predictor is anchored at the
  intercept-only fit, and a user offset is the only way to move that anchor. At
  27 positives, offsets of 0, -3 and +3 gave ESS(level) of 23.8, 24.5 and 27.3
  (probit) and 75.3, 54.8, 62.1 (logit). If the level had to travel, displacing
  it by three units on the probit scale would have shown up. It does not.
- **From a cold start the level arrives in about 25 draws**, displaced or not:
  -0.47, -1.95, -2.31 at draws 1, 10, 25, then flat, with `offset = +3` and with
  `offset = 0` alike. There is no slow approach to shorten.
- **The level is not a separately stuck coordinate.** At 27 positives ESS(level)
  is 23.8 against ESS of the *centered* predictor of 50.2 -- the whole fit mixes
  at that rate. A block move on the level cannot fix a shape that is equally slow.

The cross-check that settles it: **dbarts, an independent implementation of the
same augmentation with the same anchoring, reproduces the number exactly.** On the
same data, `ESS(level) 24.4`, `acf1 0.952`, `sd(level) 0.104`, against bartisan's
`24.8`, `0.948`, `0.101` under hard rules. Two implementations agreeing to three
digits is a property of BART with 27 events, not a defect in either.

The mechanism is the prior. The leaf prior here is proper and tight --
`sigma_mu = 3 / (k sqrt(num_trees))`, about 0.21 at the defaults -- and the
predictor is anchored at the null fit, so the level's conditional is sharp and it
starts where it belongs. UPG's near-flat intercept prior is the regime where the
random walk is slow, and this sampler never enters it.

**Finally, even in the paper it is a trade rather than a free win.** Figure 3's
lower panels plot inefficiency against the true intercept: UPG's curve is flat
where the Polya-Gamma sampler's is U-shaped, and the Polya-Gamma curve dips
*below* UPG's near a balanced intercept. It buys robustness to imbalance at a
cost when balanced.

What a user with 27 events actually needs is more draws, and to know that the
number is the information in the data. That is now said in `?bartisan_control`.

### Windle, Polson and Scott's saddlepoint Polya-Gamma sampler

It was on the list as the prerequisite for the negative binomial augmentation being worth anything. The Poisson-gamma route serves that purpose better, is exact, and needs no new sampler. The saddlepoint method is also itself an approximation, whose envelope could not be validated against the paper from here — and shipping a delicate approximate sampler on the strength of a half-remembered derivation would be the wrong trade for a package whose selling point is exactness. The same reasoning applied later to the Kolmogorov-Smirnov sampler for the ordinal logit, and there an exact route was found.

### Cox proportional hazards

Deliberately not supported: the partial likelihood couples observations through risk sets, so it does not decompose into a sum of per-observation terms over the observations reaching a leaf, which is what the leafwise Laplace approximation requires. The AFT families cover the same ground.

## Validation

Recorded here because it is the evidence for believing the sampler is correct.

- **Detailed balance.** Shrinking the prior weights to nothing makes the likelihood flat, so the target collapses to the tree prior. Observed 2.5096 leaves per tree against 2.5087 expected by backward recursion — 0.04% error — and the forest scale matches `sqrt(num_trees) * sigma_mu` to 0.6%. This is the test that would catch any wrong acceptance ratio, and it is in `test-recovery.R`. The flat-likelihood premise was itself verified: the Gaussian residual scale is still drawn, and a near-zero weighted sum of squares could in principle send the precision to infinity, but the half-Cauchy prior rejects the blow-up and the likelihood contributes 2.3e-9 relative to the prior curvature.
- **The whole prior distribution of tree sizes**, not just its mean. Against an exact reference computed by recursion, a single-tree chain reproduces the distribution with a chi-square p of 0.98 once scaled to the effective sample size, and 0.48 on a thinned subsample.
- **All three acceptance ratios**, re-derived from scratch: the prior ratio, the move-type probabilities and eligible-node counts in each direction, and the requirement that the two directions build the *same* proposal, including the boundary cases of a single-leaf tree and a collapsing root.
- **Predictor/forest invariant.** Replaying the stored trees must reproduce the additive predictor recorded while sampling. Agreement is ~1e-9 across every family, and drift is bounded at 7e-10 after 2000 draws. This is the single most productive test in the suite: it caught the double-counted parent contribution in the birth and change moves, the missing `na_rule` in the change move's rollback, and the chain-pooling offset arithmetic.
- **Signal recovery.** On the Friedman function, correlations with the truth of 0.94 (Gaussian), 0.94 (binomial logit and probit), 0.97 (Poisson), 0.94 (negative binomial), 0.99 (gamma), 0.90 (ordinal), 0.94–0.97 (AFT). Nuisance parameters recover: gamma shape 3.96 against 4, AFT scale 0.78 against 0.8 under 30% censoring, ordinal cutpoint gaps 1.39 and 2.91 against 1.5 and 3.
- **Weights are frequency weights.** Fitting with `weights = 3` and fitting on data with every row replicated three times produce *bit-identical* chains under the same seed. So they are correct for frequency weights and binomial trials; for survey or balancing weights the point estimates are right but intervals reflect the sum of the weights rather than the actual information.
- **Densities match R.** `predict(type = "density")` agrees with `dnorm`, `dpois`, `dbinom`, `dgamma` and `dnbinom` to ~1e-16.
- **A saved fit round-trips** through `saveRDS` and predicts identically in a fresh session — a known failing of several BART packages.
- **Every family's analytic score** against a central difference of its own log density, in `test-derivatives.R`.
- `Laplace2` by hand: solving `L' d = z` gives covariance equal to the inverse precision, and the density constant `-log(2pi) + 0.5 log|P|` is right.

### Bugs found by auditing rather than by tests failing

1.  `BinomialFamily::log_norm_const` called `Rf_lchoose(w, w * y)`. The success count is recovered as trials times proportion, and **that product is not reliably an integer in floating point** — 1354 of the ~20,000 (successes, trials) pairs with trials <= 200 give a non-integer. `lchoose` then rounds and warns once per observation; a Bayesian bootstrap on a binomial model emitted 50+ warnings. Rewritten with log-gamma, which is the natural continuous extension, returns exactly zero for binary data at any weight, and still matches `dbinom` to 1.5e-15.
2.  `total_loglik()` omitted `log_norm_const`, so the reported `loglik` was off by the sum of log binomial coefficients — 297.9 in a 60-observation example — while the documentation called it "the log likelihood".
3.  **`sigma_mu_ramp = 1` silently disabled the leaf-scale update for the entire sampling phase.** The restore was conditioned on reaching iteration `num_ramp`, which never happens when the ramp spans all of warmup.
4.  The multinomial probability computation in `category_probs()` exponentiated the predictor directly and would overflow. Replaced with a shifted softmax.

### Claims that were wrong and were corrected by measurement

Kept together because the pattern is the lesson.

- "Fill-plus-indicator is MIA done badly." It ties MIA on accuracy in every cell measured.
- "The cutpoint block is 77% of an ordinal fit at K = n." It is 0.4%; the comparison confounded cutpoint cost with tree cost.
- "The O(n) bandwidth update dominates." It was 8% at the time — though after everything else got faster it became 46%, which is a different statement.
- "The remaining gap to dbarts is family dispatch." It was, before the conjugate shortcut; devirtualizing bought 1.05–1.20x, not the 4x implied.
- "Removing the materialized buffers will be worth 15–25%." About 4%: the buffers never leave L1.
- "`resize()`'s zero-fill will cost more than the capacity check it saves." The opposite, by a factor of five: it is worth 20% of a soft fit.
- "A bounded gate is faster because it truncates the far side." It is faster because it has no `exp()`; at a bandwidth where it truncates nothing it is still 1.45x faster.
- "Fixing the bandwidth is faster and more accurate." Only on smooth functions; on a step function it more than doubles the error.
- Coverage was documented as "about 90% to 94%". Measured at the defaults it is 0.95 (Gaussian), 0.91 (binomial), 0.96 (Poisson), 0.96 (gamma).
- The timing table could not be reproduced to the precision it was stated at: re-measurement on a clean build came out 12–25% higher in every cell. Restated as an anchor plus ratios.
- "Raising the inverse Wishart degrees of freedom shrinks the latent correlations towards zero." It pushes them towards one: 0.38, 0.79, 0.97, 0.996 at `nu` of 3, 10, 50, 300 against a truth of 0.7. `Psi = I` is swamped by a residual scatter of order N, so all a large `nu` does is concentrate the draw on that scatter. The knob was written, measured, and removed.
- "The exponential form for soft rules would pay off on the default configuration the way it did on hard rules." A forced-on scratch build put the ceiling at 1.05x for Poisson and 1.10x for gamma, against 1.86x and 1.89x under hard rules. Dropped.
- "The multinomial-Poisson transformation changes the identification, so it is a larger change than the negative binomial's." It changes nothing: the gamma latent is integrated out exactly, so neither coding is disturbed. The family was 60 lines.
- "The prize for the multinomial-Poisson route is the mixing." The mixing prediction was right and irrelevant: Polya-Gamma is 3.9x faster and wins on ESS/s by 3.1x at K = 3 and 5.0x at K = 6. The exponential form loses to the quadratic form whenever both are available.
- "The multinomial augmentation is a modest gain bought with a severe loss of mixing, worth 1.6x." Re-measured, 9.6x and 10.1x in ESS/s with mixing at 0.66x and 1.09x. It is now a default.
- "The zero-inflated gain would be smaller than the multinomial's, because only one of its two forests gains." 3.9x to 10.1x in ESS/s, and it pays under soft rules where neither forest gets the exponential form -- the direct target's log-sum-exp was the expensive part, not the missing shape.
- "The Gaussian hard-rule fit regressed by 35%." It had not: two consecutive benchmark runs of the same build read 0.441 s and 0.593 s, and a best-of-five standalone measurement read 0.426 s both times. `_dev/benchmark.Rmd` defaults to two replicates, which is not enough to support a claim about a factor near two.

## Log: what limits the parallel pass is memory bandwidth (a correction)

`_dev/diagnose-timing.R` run on ten cores, which is the measurement the entry
below could not take. `_dev/diagnose-timing.rds` holds it.

| n | draws | 1 worker | 2 | 4 | 8 | share of the fit |
|---|---|---|---|---|---|---|
| 500 | 6.4 MB | 0.59 s | 0.36 s | 0.24 s | 0.22 s | 53% |
| 2000 | 25.6 MB | 2.40 s | 1.35 s | 0.85 s | 0.81 s | 79% |
| 8000 | 102.4 MB | 9.37 s | 5.10 s | 3.19 s | 2.80 s | 78% |

**Two things it says.** Deferring the pass was right: sequentially it is 78% of a
four-chain fit's own time at any size worth parallelizing, so it was most of what
a `diagnose()`-calling workflow spent. And the speedup **capped at about three**
however many workers were given: 2.72x, 2.98x, 3.35x at eight workers, with
eight barely beating four (2.80 s against 3.19 s at n = 8000).

**The cap had a findable cause.** Fitting `T(p) = O + W/p` gives R-squared of .985
to .995 with an overhead `O` of 0.15 s, 0.46 s and 1.48 s against data of 6.4,
25.6 and 102.4 MB, so `O` tracks the *size of the draws* rather than anything
computational. And observed `T(8)` was consistently worse than that model
predicts (2.80 s against 2.45 s at n = 8000), which says `O` grows with the
worker count too: the marginal 4-to-8 step should have saved 0.97 s of work and
saved 0.39 s, so 0.58 s went somewhere, about 0.145 s per added worker at 102 MB
-- which is what one 102 MB serialization costs.

Measured directly with `future::getGlobalsAndPackages()`: the closure referred to
`wide`, so `wide` was a *global*, and a global goes to **every** worker.
**102.4 MB each, 819 MB of serialization at eight workers.**

**The fix is to cut the blocks in the calling session** and map over those, so the
mapped element is the worker's own columns and `wide` is never a global. Globals
per worker went from 102.4 MB to 0.03 MB, and the total crossing the boundary is
now one copy of the draws whatever the worker count. It costs holding the blocks
alongside the draws for the length of the pass, which is one extra copy in the
main session.

**It also made the pass bit-identical across worker counts**, which it was not
before: the same synthetic case that had 587 of 3600 entries differing in the
last bits now differs in none at 2, 3, 4 and 5 workers. Every difference before
was in the effective sample sizes and never in either R-hat column, so it was in
the FFT path, but the mechanism was never established and is not claimed here.
The test keeps its tolerance rather than asserting exactness, on the grounds that
an unexplained agreement is not something to depend on.

### The correction, after re-running on ten cores

**The fix did nothing, and the diagnosis above is wrong.** Everything from "The
cap had a findable cause" onward is left standing as the record of what was
believed; this is what measurement says instead.

Re-run with the slicing in place, against the same run above at eight workers:

| n | before | after |
|---|---|---|
| 500 | 2.72x | 2.58x |
| 2000 | 2.98x | 2.77x |
| 8000 | 3.35x | 3.52x |

Unchanged, within run-to-run variance. Head to head at n = 8000, the two forms
were the same speed at every worker count (the slicing +0.9% at eight, +5.7% at
four), even though `getGlobalsAndPackages()` reports 102.4 MB of globals per
worker for the index form against 0.03 MB for the slicing. **The 819 MB was a
size, and a time cost was inferred from it that does not exist.** Exporting the
same object to a persistent `multisession` worker is not what its size suggests,
and the honest lesson is that a size measured with `getGlobalsAndPackages()` is
not a measurement of anything until it is timed.

**What the ceiling actually is.** Profiled at n = 8000 on eight workers: the
whole pass computes in 9.34 s sequentially, one eighth of the columns takes
1.17 s with the machine to itself, and the parallel run takes 2.36 s. Slicing
the blocks costs 0.016 s and `cbind`-ing the results 0.000 s, so there is no
serial section to speak of. Giving each of `p` workers the *same* 1000 columns
and timing them concurrently:

| workers | per-worker compute | against one |
|---|---|---|
| 1 | 1.17 s | 1.00x |
| 2 | 1.26 s | 1.08x |
| 4 | 1.32 s | 1.13x |
| 8 | 1.97 s | 1.69x |

The work per worker is identical, so a pass with cores to spare would be flat.
1.17 x 1.69 + 0.37 of transfer accounts for the 2.36 s measured, so `O` in the
`T(p) = O + W/p` fit is not a serial section that could be removed; it is
contention, which that model can only represent as one.

**And the contention is the core topology, not memory bandwidth.** Calling it
bandwidth was the same mistake as the 819 MB: a cause assigned without a test
that could rule it out. The machine is an M4, **four performance cores and six
efficiency ones**, not ten equal ones. Running p copies of a compute-bound loop
over 8 KB of data, which fits in L1 and so cannot be bandwidth-limited, against
p copies of a loop that streams 100 MB and reuses nothing:

| workers | compute-bound, 8 KB | bandwidth-bound, 100 MB |
|---|---|---|
| 2 | 1.08x | 1.04x |
| 4 | 1.16x | 1.14x |
| 6 | 1.43x | 1.51x |
| 8 | 1.65x | 1.59x |

The same curve, with the knee just after four. So the ceiling is four fast cores
plus a fractional contribution from the slow ones, which is what 3.5x on "eight
workers" means, and it is a fact about this machine rather than about the pass.
More than four equal cores should scale further; the script is worth re-running
if that ever gets measured.

**Reverted the slicing.** No speed difference, and it holds a second 102 MB copy
of the draws for the length of the pass. Both forms are bit-identical to the
sequential result at 2, 3, 4, 5 and 8 workers, so the claim above that the
slicing is what made the pass reproducible across worker counts is also wrong;
whatever the 587 differing entries were, it was not this.

**Also corrected: the `columns < 400` floor.** Its comment said the hand-off
costs more than the work below that size. Measured on four workers, the split
pays at every size tried: 3.09x at 2000 columns, 2.01x at 400, 3.03x at 200,
2.07x at 50. Lowered to 100, and the comment now says what the floor is really
for, which is that the first parallel call in a session has to start the workers.

**The one claim that survives.** Deferring the pass out of `bartisan()` was
right: sequentially it is 53%, 77% and 77% of a four-chain fit's own time at
n = 500, 2000 and 8000.

## Log: the diagnostics moved to diagnose(), and the pass got parallel

Follows the profiling entry below, and settles the two questions it left open.

### `fit$rhat` is gone; `diagnose()` is the only route

`diagnose()` never read `fit$rhat` -- `diagnosis_table()` recomputes everything
from the draws, including its own per-observation loop -- so a multi-chain fit
followed by `diagnose()` paid the cost twice. It now runs in one place, when it
is asked for. `chain_diagnostics()` was deleted; `scalar_draws()` moved to
`R/diagnose.R`, which is its only remaining caller.

Measured, four chains at n = 2000: **the fit went from 7.9 s to 3.34 s**, and
`diagnose()` is 2.36 s when wanted.

**This made the previous entry's progress fix redundant, and would have broken
the bar if left alone.** The convergence pass had just been given 50 steps of the
fit's progressor; with the pass gone those steps would never be spent and the bar
would have stalled at 80% forever. `progress_reporter()` is back to sizing the
sampling alone, which is now the whole of what a fit spends time on, and the pass
carries its own `diagnosis_reporter()` inside `diagnose()`.

Touched, because they all named a thing that no longer exists: `bartisan()`'s
`@returns`, prose in `?diagnose`, `?bartisan_control`, `?bartisan-interop` and
`?bartisan-families`, five tests in `test-chains.R` and one each in
`test-interop.R` and `test-random.R`, and four vignettes. Two of those were only
found by the test suite: the first sweep grepped `test-chains.R` alone.

**`vignette("diagnostics")` was rebuilt around `diagnose()`**, since it had been
built around printing `fit$rhat`. Its two worked examples now show
`diagnose(fit)$table`, and both of their reading paragraphs had to be rewritten,
because the numbers they described were no longer the numbers on screen: the old
text said "`rhat` above 2 ... effective sample sizes under 10" where the table
now reads 1.26 to 1.86 and 6 to 30, and "the log likelihood and the residual
standard deviation have converged" where the log likelihood is 1.14. The
replacement leans on `rhat_late`, which the old table did not carry and which is
what separates a short warmup from chains that settled apart.

### The pass runs over a `future` plan

`diagnosis_columns()` splits the columns into one chunk per worker when there are
at least 400 of them and a plan with more than one worker, and is the same
`vapply()` otherwise -- the choice `run_chains()` already makes for the chains.
No random numbers are drawn and the columns are independent, so chunking is safe;
`cut()` gives contiguous ascending chunks and `cbind()` restores the order.

**It is not bit-identical across worker counts, and the reason is worth
recording.** On a real fit the checks and the advice come out identical and the
table agrees to ten significant digits, with a maximum relative difference of
2.1e-16. Chased to its source: chunking *sequentially* is bit-identical, the
draws arrive in the worker bit-intact, `diagnosis_stats()` on an exported matrix
is bit-identical, and every difference observed was in `ess_bulk`/`ess_tail` and
never in either R-hat column. Those two are the statistics that use the FFT, so
this is R's FFT differing in the last bit or two between processes, not anything
in the chunking. Documented in the test rather than papered over.

**The speedup is unmeasured and deliberately so.** `availableCores()` reports 1
in the sandbox this was written in, so four workers contend for one core and the
1.16x observed there means nothing. `_dev/diagnose-timing.R` sweeps sizes and
worker counts, drops any worker count above `availableCores()` so a reported
speedup cannot be workers fighting each other, and discards the first call under
each plan because starting the workers is a cost of the plan rather than of the
pass. Run it on real hardware before any speedup goes in the documentation.

## Log: why the diagnostics pass cost more than the fit

Profiled rather than guessed at. At n = 2000, four chains, 400 draws each:

| | before | after |
|---|---|---|
| `chain_diagnostics()` | 4.85 s | **1.84 s** |
| `diagnose()` | 5.31 s | **2.34 s** |
| sampling, four chains | 4.0 s | 4.0 s |

**Where the time went, and it was not the arithmetic.** The scalar rows cost
0.006 s for two quantities; the whole 4.85 s was the per-observation loop, and
within it ESS was 90% (2.1 ms of 2.33 ms per observation) against R-hat's 10%.
Two causes, both R-level overhead:

1. `stats::acf()` was called once per split half-chain per observation, so 32,000
   times for a 2000-observation fit, each building an `acf` object it then threw
   away. The profile was `colnames`, `mode`, `outer`, `deparse1`, `cbind` -- the
   bookkeeping, not the covariance.
2. It asked for `lag.max = draws - 1`, every lag, on the assumption that Geyer's
   initial positive sequence stops early. On a *forest's* fitted values it does
   not: measured, the sequence runs to **98 lags out of 98 available**, median
   and max, because `eta` mixes slowly (its R-hat is 1.2 and its ESS about 15 on
   the same fit). So `rho()` was a closure called ~98 times per observation and
   `kept <- c(kept, ...)` reallocated the vector 49 times.

**Two fixes, and the second is worth less than it looks.** Replacing `acf` with a
zero-padded FFT autocovariance (`autocovariance()`, Wiener-Khinchin, one
`mvfft()` pair for every chain and lag at once) took the pass from 4.85 s to
2.19 s. Preallocating `kept` and indexing a precomputed `rho` vector instead of
calling a closure took it to 1.84 s. What is left is `rank()` and `sort.int()` at
about 35% of the remainder, from the three or four `rank_normalize()` calls each
observation needs, which is irreducible in R without sharing them across the
three statistics or moving the loop to C++.

**On correctness, one claim had to be withdrawn.** The FFT change was first
reported as bit-identical over 1800 statistics on real draws; that test shimmed
`autocovariance` into a closure's environment and the shim silently did not take
effect, so it compared the function against itself. Measured properly, the FFT
route agrees with `acf` to 1.8e-15 on the covariance and the ESS values to
2.7e-12 absolute on quantities of order 10 to 1500, which is floating-point noise
against thresholds of 400 and 1.01, and is the same estimator Stan computes this
way. The Geyer refactor *is* exact: the `kept` vector and `extra` are **bitwise
identical** across 6000 randomized cases spanning 4 to 150 draws, 1 to 8 chains,
and autocorrelations from -0.98 to 0.995.

### The redundancy that is still there

`diagnose()` does not read `fit$rhat`; `diagnosis_table()` recomputes everything
from the draws, including its own per-observation loop. So a four-chain fit
followed by `diagnose()` pays the per-observation cost **twice**, 1.84 s and then
2.34 s, for one answer. `print()` and `summary()` never touch `fit$rhat`, and no
vignette chunk outside `diagnostics.Rmd` and one chunk of `implementation.Rmd`
does either.

Deferring it to `diagnose()` would take it out of every fit. Not done here,
because it removes a documented element of the return value and
`vignette("diagnostics")` is built around displaying `fit$rhat` -- it prints the
table, has a section called "The Rows of the `rhat` Table", and reads
`too_short$rhat` -- so it is a design decision about the package's diagnostic
surface rather than a mechanical change, and it is being put to the author.

### On parallelizing it

The export is cheap: 25.6 MB of `eta` draws reach four `multisession` workers in
0.16 s, so a `future_lapply` over chunks of observations would not be dominated
by serialization. What could not be measured here is the speedup itself, because
`parallelly::availableCores()` reports **1** in this sandbox, so four workers
contend for one core and the observed 1.16x says nothing about a real machine.
Left unimplemented rather than shipped on an unmeasurable benefit.

## Log: the progress bar finished before the fit did

Reported as *progressr* working incorrectly with several chains under a
multisession plan, looking like a bar for the first chain only.

**The relaying was never the problem, and the report's diagnosis was wrong while
the observation was right.** Instrumented with `handler_debug(uuid = TRUE)`: under
`multisession` with four chains, all 200 progression conditions arrive in the
calling session, from **one** progressor uuid and **one** owner session, against
a `max_steps` of 200. Sizing and relaying were both already correct.

**What was actually wrong is the bar's denominator.** It covered the sampling and
nothing after it, and the thing after it is `chain_diagnostics()`, which runs
*only when `chains > 1`*. Measured at n = 2000, four chains: sampling 4.0 s,
`chain_diagnostics()` 4.8 s, `fitted_from_eta()` 0.00 s. So the whole
post-sampling cost is the convergence pass, it exists only in the multi-chain
case, and a `future` plan shortens the sampling while leaving it untouched. The
fraction of the run the bar covered therefore fell as chains and workers were
added:

| plan | bar reached 100% at | of a run lasting | share covered |
|---|---|---|---|
| sequential, 1 chain | 0.7 s | 1.0 s | 73% |
| sequential, 4 chains | 3.0 s | 8.3 s | 36% |
| multisession, 2 workers, 4 chains | 1.7 s | 6.7 s | 25% |
| multisession, 4 workers, 4 chains | 1.1 s | 6.0 s | **18%** |

A bar that races to 100% in the first fifth of the run and then sits there is
exactly what "a progress bar for the first chain" looks like from outside.

**The fix** gives the convergence pass its own share of the bar:
`PROGRESS_DIAG_TICKS` of 50, added to the progressor's steps when the caller says
the pass will run, and spent by a `progress_stepper()` that turns the columns of
the per-observation and per-level loops into at most that many reports. After it,
all three multi-chain configurations reach 100% at 100% of the run.

The two phases are charged *fixed* shares rather than shares proportional to what
they will cost, which is not knowable in advance, so the bar no longer advances
uniformly in time: it moves quickly through the sampling and then slowly through
the diagnostics. That is documented, and it is the right trade against a bar that
lies about being finished.

**Still not covered**, and left alone deliberately: the one-chain case reaches
100% at about 73% of a one-second run, the remainder being model-frame setup
before sampling and object assembly after. That is a small fixed overhead rather
than something that grows with the fit, so charging it to the bar would add
machinery for no benefit.

**Two testing notes worth keeping.** *progressr* reports nothing in a
non-interactive session unless `progressr.enable` is set, which made the first
three attempts at measuring this show zero conditions and look like a much worse
bug than it was; always validate the instrument on a textbook `progressor()` loop
first. And `handlers(global = TRUE)` cannot be measured with a
`withCallingHandlers()` wrapper, because the global mechanism needs an empty
handler stack, so that path stays untested here.

**The existing tests encoded the old contract** and both failed on the fix, which
is what they were for: `3L * PROGRESS_TICKS + 2L` became
`3L * PROGRESS_TICKS + PROGRESS_DIAG_TICKS + 2L`. A new test asserts the
difference a second chain makes is one chain's ticks *plus* the pass's share, and
that `diag_ticks` is zero for one chain before any fitting happens.

## Log: *BART*'s multinomial support is binary fits, verified

The feature table's checkmark for *BART* on the multinomial row was misleading.
Read from the installed source of *BART* 2.9.10 rather than from the help page,
which says only that "P(Y=y | x) = F(f(x))" and hides the mechanism:

- **`mbart2()`** is `for (h in 1:K) gbart(x.train, (y.train == h) * 1, type =
  "pbart"/"lbart")`, so $K$ independent one-vs-rest binary fits on the full data,
  followed by `prob = exp(yhat_h) / sum_h exp(yhat_h)`. The normalization is
  applied to the outputs afterward; no forest is fitted to a multinomial
  likelihood at any point. This is exactly "binomial applied separately to each
  category".
- **`mbart()`** is the same loop with `cond <- which(y.train >= cats[h])` and
  `gbart(x.train[, cond], (y.train[cond] == h) * 1, ...)`, so $K - 1$ binary fits
  on *nested subsets*, combined as a continuation-ratio product. That is an exact
  factorization of the multinomial mass function and so a coherent model, but
  each conditional gets its own independent forest and prior, and the
  factorization runs over sorted categories, so its prior on a probability vector
  is not exchangeable in them.

Noted in the table as `✓ **separate binary fits**` with a paragraph giving the
mechanism, in the convention the table already uses for qualified checkmarks.

**A test that did not work, recorded so it is not repeated.** The obvious check
for `mbart()`'s asymmetry is to relabel the categories and compare fits, but
relabelling changes the order data reach the sampler, so the RNG stream diverges
and Monte Carlo noise swamps the effect: *BART* moved by a max of .061 and
*bartisan*'s symmetric `multinomial()` by .057 on the same test, which
discriminates nothing. The asymmetry is a property of the factorization and is
established by reading it, not by simulating it. The claim in the vignette is
therefore about the model's definition and says nothing about the size of any
practical consequence.

## Log: what the AFT predictor means, and empty forests everywhere

### The vignette was wrong about `dpm_aft()` and the time-ratio reading

The claim was that `dpm_aft()`'s predictor "is only a time ratio to the extent
that the error density is symmetric, which is exactly what `dpm_aft()` declines
to assume." That is false, and it conflated a contrast with a level.

**A contrast is a log time ratio for every AFT family, whatever the error.** With
$\log T = \eta(x) + W$ and $W$ independent of $x$, every quantile of $T$, the
mean of $T$, and the geometric mean of $T$ all scale by $e^{\Delta\eta}$. Checked
on a truth with an error of skewness $-1.00$ and a log time ratio of exactly .8:
the ratios of the 10th, 25th, 50th, 75th and 90th percentiles came out .794,
.811, .782, .757 and .792 in logs, the mean ratio .781, the geometric-mean ratio
.787. Symmetry has nothing to do with it.

**What the error's shape does govern is what $e^\eta$ is on its own**, and there
the odd family is `weibull_aft()`, not `dpm_aft()`. All three parametric families
are written $r = (y - \eta)/\sigma$, so $\log T = \eta + \sigma G$, and the
question is where each $G$ sits:

| family | $G$ | $E[G]$ | $e^\eta$ is |
|---|---|---|---|
| `lognormal_aft()` | standard normal | 0, symmetric | the median of $T$, and its geometric mean |
| `loglogistic_aft()` | standard logistic | 0, symmetric | the same |
| `dpm_aft()` | centered DP mixture | 0, asymmetric | the geometric mean, $\exp E[\log T]$, not the median |
| `weibull_aft()` | standard Gumbel-min | $-\gamma$ | the Weibull scale, the 63.2nd percentile |

Verified by fitting each on its own truth at $n = 3000$: `weibull_aft()` recovers
the generating *location* to within .005 but sits **+.286** from $E[\log T]$,
which is $\gamma\sigma = .577 \times .5 = .289$; `lognormal_aft()` and
`dpm_aft()` recover $E[\log T]$ to within .001. And directly:
$e^\eta$ falls at the 63.0 percentile of $T$ against a theoretical 63.2, where
the median is 2.27 and $e^\eta$ is 2.72.

### Why the Weibull carries a hazard ratio too

The premise that the structural part is shared is right, and that is exactly why
all four have the AFT reading. The PH reading is an extra property of the
Gumbel-min error, which is the only one making an AFT model proportional hazards
as well: $h(t \mid x) = \sigma^{-1} t^{1/\sigma - 1} e^{-\eta/\sigma}$, whose
ratio between two covariate values is free of $t$. Checked: the hazard ratio is
0.201897 at $t$ of .5, 1, 2 and 5, and equals $e^{-\Delta\eta/\sigma}$ exactly.
So one $\eta$ gives a log time ratio of $\Delta\eta$ and a log hazard ratio of
$-\Delta\eta/\sigma = -k\,\Delta\eta$.

Fixed in three places that each carried the error independently: the estimand
table and the "Three Estimands, Not One" section of `vignette("survival")` (now
"Two Estimands, Not Three", with a table for where the level sits), the survival
reference table in `vignette("families")`, and `?bartisan-families`.

### `~ 1` already generalized, and now says so

It works for every family taking more than one formula, because the change was
in family-agnostic code (`resolve_vc()`, `resolve_split_matrix()`, and the
per-forest `gamma`). Confirmed on `gaussian_ls()`, `Gamma_ls()`, `zi_poisson()`,
`zi_negbin()` and a two-predictor `custom_family()`: zero splits in the second
forest and zero variance in its predictor in all five.

The multinomial families are the one exclusion and it is the right one: their
forests are the levels of one vector-valued parameter, which is what
`joint_forests()` marks, so they refuse a formula list outright.

Statistically checked rather than assumed: `zi_poisson()` with `~ 1` on its
inflation part recovers a constant structural-zero probability of .313 against a
truth of .300, with the count surface at .075 RMSE on the log mean; and
`gaussian_ls()` with `~ 1` gives sigma .7012 against `gaussian()`'s .6987 and a
truth of .70, mean surfaces correlating at .998. The scalar is drawn under the
leaf prior rather than the built-in family's own prior on its nuisance parameter,
so these agree closely rather than exactly, which the documentation now says.

### Method in rows, condition in columns

Applied to every table that measures methods under conditions: both comparison
tables in `vignette("families")`, the two `pivot()` tables in
`vignette("survival")` (one helper feeds both), the sparsity table in
`vignette("implementation")`, and the tree-count, sparsity-strength and
categorical-rule tables in `?bartisan_control`. Transposition was done by a
script that reparses the markdown so cell contents, including bolding, survive
verbatim rather than being retyped.

Left alone, deliberately: the long-format tables that already have the method in
rows (the benchmark table, the augmentation table, the relative-cost table), the
decision tables ("if X, use Y"), and the package feature matrix in
`vignette("implementation")`, whose 35 features would become 35 columns and whose
reader compares packages on one feature along a row.

The transposition moved unit labels out of the corner cell, so the column headers
carry them now ("10 predictors", "10 per level", "5 trees"), and five prose
references to rows became references to columns.

## Log: gaussian_ls, Gamma_ls, and a review of custom_family

### `location_scale()` became `gaussian_ls()`

Renamed everywhere, including the engine's family string, which is user-visible:
it is what `print()` reports and what `augment` takes as a name. 28 files, and
the C++ class became `GaussianLSFamily`. Entries in this file from before the
rename still say `location_scale()`; that is history and was left alone.

### `Gamma_ls()`: a gamma whose dispersion is a forest

`Gamma("log")` draws one shape for the whole sample, which asserts a constant
coefficient of variation. `Gamma_ls()` puts a second forest on the log
*dispersion*, so the shape is `exp(-eta1)` per observation.

**The second predictor is a dispersion rather than a shape** so that both
location-scale families read the same way, with a larger second predictor
meaning more spread. The two differ by a sign, so nothing is lost.

**The target forms are split**, which is what made this cheap enough to be worth
having. In the mean, the log density is `-s eta0 - s y exp(-eta0)`, the same
exponential form at rate -1 that `Gamma("log")` uses, with both coefficients
free of `eta0` even though `s` now varies by observation. In the log dispersion
there is nothing to exploit, because `lgamma(exp(-eta1))` has no form, so that
forest takes the general path. Measured, the fit costs about 5.4x
`Gamma("log")` at 50 and 20 trees, essentially all of it in the second forest.

**The information for the log dispersion is the expected one**,
`s^2 trigamma(s) - s`, not the observed one. The observed version differs from
it by the score, so they agree at the mode, and the expected one is guaranteed
positive because `trigamma(s) > 1/s` for every `s > 0`. The observed one can go
negative away from the mode, which would give the Laplace proposal a negative
variance. The mean forest keeps the observed curvature, as `Gamma("log")` does,
since a strictly positive response cannot make it negative.

**Three checks, and the second is the one that matters.** Both components' scores
match central differences to 1.7e-9, and the mean component's information to
2.4e-8. The log density matches
`dgamma(y, shape = 1/phi, rate = shape/mu, log = TRUE)` to **7.1e-15**, which is
what establishes that the parameterization is the one claimed rather than merely
self-consistent. And a dispersion that varies with a predictor is recovered at a
correlation of .973, worth 62.9 log points over `Gamma("log")` on the same data.

### An intercept-only forest, which `~ 1` now means

Checking `Gamma_ls()` against `Gamma()` needs the scale forest held constant, and
`bartisan(list(y ~ x, ~ 1), ...)` was an error: "no predictor left to split on".

The mechanism to allow it already existed. A nuisance parameter carried as a
trailing forest is pinned by giving it a **branching probability of zero**
(`src/model.cpp`), which makes every tree in it a stump, so the forest is one
drawn scalar. `gamma` is already a per-forest control. So a forest whose formula
names no predictor now sets `gamma = 0` for that forest instead of erroring, and
`resolve_split_matrix()` gives it a uniform placeholder column that is never
read. The error survives for the case it was actually protecting against, a
formula whose terms `split_prior` has all zeroed, and its message now says so.

**One trap in this.** `resolve_vc()` has an early return for the common case of
no `vc()` terms, so the first version of the change had no effect at all on any
formula without a varying coefficient, which is every formula this feature is
for. The forest kept splitting; `fit$counts` showed 714 rules in a forest that
was supposed to have none. Both return paths now carry `pinned`.

Verified: with `~ 1` the scale forest has **0** splitting rules and its predictor
has zero variance across observations, and against `Gamma("log")` on the same
data the implied shapes are 4.09 and 4.13 against a truth of 4, the mean surfaces
correlate at .993, and the log scores differ by 3 points in 1300.

### `custom_family()`: the changes are sound, with two holes now closed

Reviewed as asked. The change from `aux_start = 0` to `aux_start = NULL` is a
genuine fix, not just a tidy-up: the old code decided whether nuisance
parameters existed with `missing(aux_start)`, so
`do.call(custom_family, list(logdens = f, aux_start = 0))` manufactured a
parameter and then errored, where the same call written literally did not.
Reproduced against the old version. A `NULL` default makes the two
indistinguishable by value. Naming a parameter through `aux_start = c(shape = 1)`
is also new and useful.

Two holes in the new branch, both fixed:

- **A partially named `aux_start` produced an empty parameter name.**
  `c(a = 1, 2)` has names `c("a", "")`, no duplicates, so it passed the guard and
  `aux_names` became `c("a", "")`. That empty name would have labeled a column of
  `fit$aux` and a row of `summary()`. `aux_names = c("a", "")` was already an
  error, so the same input was accepted or refused depending on which argument
  carried it.
- **Duplicate names were silently discarded.** `c(a = 1, a = 2)` fell back to
  `aux1`, `aux2` with no message, where `aux_names = c("a", "a")` errors.

Both came from the check living in the other branch. `aux_names` is now derived
first and checked once, whichever argument it came from, with the message naming
that argument.

### `error_density(plot = TRUE)`

Returns a `ggplot` of the posterior mean density with the pointwise interval as
a ribbon; `plot = FALSE` still returns the data frame. *ggplot2* is a soft
dependency, so `TRUE` checks for it and says what to do instead.
`vignette("survival")` keeps the hand-drawn version, because there the point is
to overlay the normal a `lognormal_aft()` fit would have assumed, which is what
the values are for.

### The multinomial probit correlations need no reparameterization

Asked whether to constrain them through a bounded transform, on the evidence
that `vignette("families")` reported an interval reaching 1.07.

**The vignette was reporting interval *widths*, and said so ambiguously.** The
sentence read "95% intervals from .38 to 1.07 wide", which scans as endpoints.
Rewritten to name the widths explicitly.

**No transform is needed, and an element-wise one would be wrong.** `Sigma` is
drawn from an inverse Wishart and rescaled to `trace(Sigma) = C`, so every draw
is positive definite by construction and every correlation it implies is
strictly inside (-1, 1); measured over 400 draws, the range was [.058, .622] and
the trace held at 2.000000. Squashing each correlation through a scaled normal or
logistic distribution function would also be the wrong repair past three
categories, where the correlations are entries of a matrix that has to be
positive definite *jointly*, which a per-entry transform cannot enforce. Weak
identification here shows up as a wide posterior, which is the honest signal.

### The two-category `cloglog` footnote, and what it says

`ordinal()`'s documentation says a two-category response is "exactly binary
regression on the same scale". That is true for `logit` and `probit`, whose
errors are symmetric, and **false for `cloglog`**, whose smallest extreme value
error is not. The cumulative-link form puts the error in with the opposite sign
to the binomial's, so with the single cutpoint pinned at zero the ordinal model
reports `P(Y = 2) = exp(-exp(-eta))`, the log-log link, against the binomial's
`1 - exp(-exp(eta))`.

Checked per draw rather than asserted: the ordinal fit's probabilities match
`exp(-exp(-eta))` to **1.1e-16** and differ from `1 - exp(-exp(eta))` by .264.
Reversing the levels and negating the predictor turns either into the other.

### The positive-outcome comparison, re-measured against the right rivals

`_dev/positive-sim.R`, which is now reproducible rather than pasted from a lost
script. The old table compared `Gamma("log")` against `gaussian()` and `dpm()`
fitted to `log(y)`, and both of its columns were meaningless as a result: the
RMSE compared an estimate of $\log E[Y]$ against an estimate of
$E[\log Y]$, and the log score compared densities taken with respect to
different measures. `gaussian("log")` is the like-for-like Gaussian, since the
composed link puts the forest on the log mean while leaving the error additive.

`dpm()` has no link argument and appears only on the raw scale. Its additive
predictor is *defined* to be the conditional mean, which is what identifies the
centered mixture, so a log link would be a different construction rather than a
composition. Worth knowing before anyone reaches for `dpm("log")`: it silently
passes `"log"` to `nu` and errors there.

800 train and test, 50 trees, 500 draws after 500 warmup, four replicates,
medians. Cells are RMSE / log score:

| Errors | `Gamma("log")` | `Gamma_ls()` | `gaussian("log")` | `dpm()` | `ordinal("probit")` |
|---|---|---|---|---|---|
| gamma, constant dispersion | 0.395 / **-1938** | **0.378** / -1940 | 0.433 / -2096 | 0.793 / -2005 | 0.434 |
| gamma, varying dispersion | 0.620 / -2218 | **0.461 / -2156** | 1.265 / -2466 | 1.442 / -2251 | 0.858 |
| lognormal | 0.548 / -2074 | 0.537 / **-2073** | 0.574 / -2366 | 0.913 / -2109 | **0.397** |
| heavy tail | 0.465 / -2097 | 0.440 / -2132 | 1.229 / -2261 | 0.628 / **-2042** | **0.432** |
| seconds | 7.6 | 48.7 | 16.8 | **2.4** | 3.1 |

**`Gamma_ls()` has the property that makes a flexible family safe to default
to.** It cuts RMSE 26% and gains 62 log points where the dispersion varies, and
where it does not it is two log points from `Gamma("log")`, well inside the
replicate spread. That is `dpm()`'s relationship to `gaussian()`, one level up.

**`gaussian("log")` is the useful negative result.** It targets the right
estimand and is still behind on log score in every row by 150 to 310 points,
because the error stays additive and of constant variance where the data are
multiplicative and skewed. It is also *slower* than `Gamma("log")`, at 16.8
seconds against 7.6, because a link the engine does not compile is applied from
R once per leaf visited. Worth remembering as the cost of the composition path.

**The heavy-tail log score is stable this time** and is reported, where the old
table had to withhold it. The old truth put the score at the mercy of one test
point; a contaminated gamma whose two components both have mean one does not,
and the across-replicate range is now about 50 points rather than five orders of
magnitude.

### A trap worth recording: `roxygenise(load_code = "installed")`

It documents whichever build is on `.libPaths()`, not the source tree. Run
against the previously installed version it silently **dropped
`export(Gamma_ls)` from NAMESPACE** and omitted the function from the `\usage`
block, while leaving every prose mention of it in place, so the page looked
right. `R_LIBS` has to point at the build under test.

## Log: the documentation voice pass

### Applying the `/r-doc-style` skill across the package

The skill at `~/.config/agents/skills/r-doc-style/` is derived from the
documentation of the published packages with agent-written lines excluded by `git
blame`. Applied here to the roxygen blocks, the README and seven of the nine
vignettes; `causal.Rmd` and `bartisan.Rmd` were held out of scope. Working notes
and the operative rule list are in `_dev/DOCSTYLE.md`.

**The diagnostic was the useful part, and it was mechanical.** Measured against
the skill's twenty-four deltas before any edit, the package was uniformly in a
different register than the corpus rather than unevenly so, which is what made
the pass tractable: 0 uses of `Default is` against 6 of `Defaults to`; 0 uses of
`Allowable options include`; 0 type markers of the form `` `logical`; ``; **0
glosses** against the corpus's 363 `i.e.`/`e.g.`; 1 `Note that` against 54; 0
`@note` blocks against 19; 0 `@inheritParams` against 38; **230 dashes as
subclause delimiters** (75 in roxygen, 155 in markdown) against 7 in twelve
thousand corpus lines; and in the vignettes **73 instances of "you" against 0 of
"we"**, which is the exact inversion of the corpus's 157 to 35.

What was already right: every one of the 124 `@param`s started lowercase,
`@returns` was used rather than `@return` on every exported page, titles were
sentence case without periods, bold lead-ins were already the device for parallel
definitional material, and none of the skill's banned words appeared anywhere.

**Two factual errors surfaced from reading rather than from the style rules.**
`bartisan()`'s `@returns` named the returned class as `bartisan`, where the code
assigns `bartisan_fit`. And `split_prior`'s `@param` read "Weights must be finite
and not than negative", a mangled sentence that was presumably meant to say
nonnegative, which is what the validation actually enforces.

**The README was the structural decision.** At 8,881 words it was ten times the
house length of 450 to 1050, and it was functioning as a second manual: it
duplicated the control page, the families vignette and the effects vignette, and
carried four sections that existed nowhere else. It is now 795 words in the house
shape (Overview, Installation, Examples, Citing, Questions and Bug Reports), with
one worked example on `rhc` that runs. Every duplicated section was checked
against the roxygen and the vignettes before being cut; the four that had no
other home were **moved rather than deleted**, into `vignette("implementation")`,
which is the "how it works" reference document and the right place for them:

- the feature comparison against *dbarts*, *BART*, *flexBART*, *SoftBart*,
  *bartMachine* and *stochtree*;
- the benchmark table, the two measurement errors it corrects, and the
  phase-by-phase account of where the time goes;
- the three correctness checks and the interval-coverage numbers;
- the notes and limitations, including the ordinal identification chart.

`_dev/README-moved.md` holds the pre-move text verbatim in case any of it wants a
different home. `implementation.Rmd` went from 409 to 908 lines and still renders
in 3.4 seconds.

**Two things worth remembering for the next pass.** The vignette headings were
the least mechanical part: the file had essayistic headings ("A trap: the
densities must be on the same scale", "Three things importance is not", "Why
slopes are unreliable here") where the corpus names the thing documented and puts
the triggering argument in parentheses. Converting those is a judgment call per
heading and cannot be scripted. And the `#` level-one headings inside `@details`
and inside `implementation.Rmd` were wrong in a way that renders without
complaint: in roxygen a bare `#` becomes a top-level help-page section rather than
a subsection of Details, and in a vignette `#` collides with the YAML title. Both
became `##`.

## Notes

### Candidate names

**Settled: the package is `bartisan`, and the rename landed in `cdc3278`.** What collided was the *old* name `genbart`, case-insensitively against the archived CRAN package `genBart`; `bartisan` collides with nothing. The table below is the shortlist that was drawn up at the time and is kept as a record of what was considered.ed package names, and are free. Also note `flexBART`, `SoftBart`, `dbarts`, `bartMachine`, `bartCause` and `stochtree` exist, and that `gbart` is the main *function* in the `BART` package, so it should be avoided even though the name is free.

| Candidate | Reading |
|---|---|
| `anybart` | BART for any likelihood — the actual contribution, short and unambiguous |
| `glmbart` | signals the `glm()` interface and the family system; understates AFT, location-scale and custom families |
| `omnibart` | same idea as `anybart`, slightly more formal |
| `laplacebart` | names the mechanism; good for a methods audience, longer to type |
| `bartlap`, `lapbart` | shorter forms of the same |
| `bartfam`, `familybart` | the family system |
| `beyondbart` | echoes Linero's subtitle, "Beyond Conditional Conjugacy" |
| `bartleby` | memorable, says nothing; Melville's scrivener was a copyist |

`flexBART` is the closest competitor on interface — formula-based, heteroskedastic mean-and-variance ensembles, varying-coefficient models — so a name that reads as a variant of it is worth avoiding.

### Traps and housekeeping

- **`.git/` does not carry the `com.dropbox.ignored` xattr, and now has history in it.** It carries `com.dropbox.attrs`, so Dropbox is tracking it; the ignore xattr was never applied, and it is now 696 KB with one commit rather than the 14-file skeleton it was. Applying the xattr retroactively is the action `~/.config/agents/AGENTS.md` § "Dropbox Sync Exclusions" warns against, and rebuilding the directory in the canonical order was declined earlier, so it is still untouched. Normal commits are append-mostly and Dropbox handles them; the latent problem is a later `git gc` or `filter-branch`, which rewrites paths Dropbox has already indexed. Worth deciding before the repository grows.
- `README.md` is generated from `README.Rmd` with `knitr::knit()`, and needs re-knitting whenever the Rmd changes. `figure/` holds its one plot and is committed so GitHub can render it.
- **`_dev/` is mostly not committed.** `.gitignore` keeps `_dev/benchmark.Rmd`, which README.md and this file both point at, and excludes the rest. `_dev/Reproduce/` is Linero's JASA replication package — seven third-party GPL-2 packages — which is here as a reference and is not ours to redistribute; publishing it under this repository's name is a decision for the maintainer, not a side effect of committing. `_dev/benchmark.html` is regenerable output.
- **Error-message regexes in tests must not span a line break in the message's source string.** testthat pins `cli.condition_width` when it runs a package's tests, which stops cli reflowing a condition message, so the source string's own indentation survives into the message — under `R CMD check` only. A regex crossing one of those breaks passes when the tests are run any other way and fails under check, which is how three of them got through. The scratch test runner now sets the same option so the two agree.
- **Never reuse a seed for the predictors and for the response.** `sim_x(seed = k)` followed by `set.seed(k)` makes the noise a deterministic function of the predictors, because both draw from the same restarted stream. For a continuous response the linear correlation is only about 0.008, so it hides; for `rbinom(n, 1, p)` it is catastrophic — with `p = plogis(2 * x1 - 1)` on the same stream that produced `x1`, every draw came out zero, which is what surfaced it. Twenty-two tests written across several sessions had the pattern and were changed to offset the response seed. It costs nothing to avoid and a recovery test built on coupled noise is not testing what it claims.
- **`na.action` defaults to `na.pass`** now, so a test that expects rows to be dropped has to ask for `na.omit` explicitly. And note what the fix to that default exposed: `model.frame()` is called through a call rebuilt from `match.call()`, so *any* argument of `bartisan()` that is forwarded to `model.frame()` and left at its default is absent from that call and picks up `model.frame()`'s default instead. Adding a default to `subset`, `weights` or `offset` would be swallowed the same way.
- The package is installed only in a scratch library, because the sandbox cannot write to the system R library. Reinstall outside the sandbox to use it from a normal session.
- `R CMD check --as-cran` reports three CRAN-incoming issues that are not code defects: the name collision above, a development version component, and a GitHub URL that 404s because nothing has been pushed. Plain `R CMD check` is `Status: OK`.
- The MCMC engine is adapted from Linero's `FlexBart` (GPL-2) in `_dev/`. The package is GPL (>= 2), which is compatible; Linero is credited in `DESCRIPTION` as contributor and copyright holder.
- Theory caveat: Linero's posterior-concentration theorem assumes a bounded variance function, which excludes Poisson with the log link and the gamma family. The sampler is still exact; only the asymptotic guarantee is unproven for those.

### Benchmarking against bartMachine and stochtree

`bartMachine` was silently absent from a whole benchmark run while reporting as installed. Two causes, both now fixed in `_dev/benchmark.Rmd`:

- It needs `options(java.parameters = c("-Xmx8g", "--add-modules=jdk.incubator.vector", "-XX:+UseZGC"))` set **before** the JVM starts, i.e. before anything touches the package. It is now set in the first chunk.
- It needs `JAVA_HOME` set. On this machine `/usr/bin/java` reports "Unable to locate a Java Runtime" despite two Temurin JDKs being installed, and every fit fails. `JAVA_HOME=/Library/Java/JavaVirtualMachines/temurin-21.jdk/Contents/Home` in `~/.Renviron` fixes it for every R session.

The deeper problem was that the document could not tell the difference. Backend detection called `set_bart_machine_num_cores()`, which succeeds even when fitting does not, and `record()` dropped any fit whose timer returned NA — so a broken backend produced a table with no row for it, indistinguishable from an uninstalled one. Detection now fits something tiny per backend and reports the error text, and failed fits are listed rather than dropped.

`bartMachine` reports no posterior draws of the fitted value through `bart_predict_for_test_data()`, so its ESS column is genuinely unavailable, not missing by accident.

`stochtree`'s grow-from-root warm start (`num_gfr`) is a different algorithm from the MCMC every other package here runs, so the Gaussian task times it both at 0 and at its default of 5. It supports continuous/identity, binary/probit, binary/cloglog and ordinal/cloglog.

### The family documentation moved to a vignette, and two claims in it were wrong

The `bartisan-families` help page had grown to roughly 380 lines of Details, most of it exposition rather than reference. It is now `vignette("families")`, with a section per family, a section on the inferred default, a section on choosing a family, per-family guidance on choosing among the links (and among the zero-inflated and AFT models, which differ by model rather than by link), and a section on `custom_family()`. The help page keeps the family table, the identification facts that output could be misread without, and a pointer to the vignette. Both vignettes now carry citations in `vignettes/references.bib`; every entry was checked against CrossRef, which caught a `\references` block where the Murray (2021) entry had been split in half by the George et al. (2019) entry pasted into the middle of it.

**The multinomial probit correlation claim was wrong, and wrong in an interesting way.** The documentation said the correlation is "attenuated" because a nonparametric mean absorbs part of the dependence, and then printed 0 → 0.28, 0.5 → 0.55, 0.8 → 0.83 as evidence — three numbers all *larger* in magnitude than the truth. The test carried the same wrong comment.

Re-measured at 1000 draws after 1000 warmup, sweeping the true correlation, with the sum of trees as the mean:

| True | n = 900 | 95% interval | n = 3000 | 95% interval |
|---|---|---|---|---|
| -0.6 | -0.748 | [-0.883, -0.504] | -0.588 | [-0.787, -0.407] |
| -0.3 | -0.347 | [-0.793, 0.144] | -0.043 | [-0.288, 0.492] |
| 0.0 | -0.567 | [-0.855, 0.047] | 0.083 | [-0.176, 0.245] |
| 0.3 | 0.110 | [-0.488, 0.584] | 0.424 | [0.243, 0.563] |
| 0.6 | 0.456 | [0.117, 0.736] | 0.654 | [0.490, 0.768] |
| 0.8 | 0.648 | [0.392, 0.844] | 0.600 | [0.464, 0.729] |

Eight draws of the data at a true correlation of zero and n = 900 gave posterior means of -0.214, -0.189, 0.330, -0.024, 0.095, -0.600, -0.205 and -0.102.

So the parameter is **weakly identified rather than attenuated**. At 900 observations the sweep is not even monotone, intervals run 0.38 to 1.07 wide on a parameter confined to (-1, 1), and a true zero can come back at -0.6. By 3000 it behaves: right sign everywhere, within about 0.2 of the truth, intervals 0.26 to 0.78. The problem is variance, not bias, and the documentation now says to read the fitted probabilities rather than the covariance. The test was rewritten to assert the deterministic constraints (trace, positivity, unit ball) plus separation of the two extremes, which is all that holds at a size a test can afford; the old `abs(estimate at zero) < 0.4` assertion passed only because 400 draws had not yet reached where 1000 draws go.

### avg_comparisons() returning exactly zero: neither package's bug

Reported as `avg_comparisons(fit, variables = "treat", newdata = subset(treat == 1))` giving `Estimate 0, 2.5% 0, 97.5% 2154` on `MatchIt::lalonde`. Two checks settled it.

The plumbing is exact. The draws `marginaleffects` receives match a hand computation from `predict()` — build the treat = 0 and treat = 1 frames, take `rowMeans(dh - dl)` per draw — to 5.9e-12, elementwise.

The zero is a **posterior atom**. In any draw where no tree splits on the contrasted variable, the fit does not depend on it, so the two counterfactual predictions are identical and the difference is exactly zero. The Dirichlet sparsity prior (`update_s`) is what makes those draws common. `marginaleffects` centers a posterior at its median, so once the atom holds more than half the mass the reported estimate is exactly zero however large the rest is. Reproduced at the defaults: `treat` was in none of the 50 trees in 64% of draws and the contrast was exactly zero in 65%, median 0, mean 197.

Also visible in the same measurement, and worth its own note: **the variable-selection state mixes slowly.** Four chains at the defaults put `treat` in 82%, 46%, 78% and 100% of draws and gave average contrasts of 649, 426, 847 and 1039. With `update_s = FALSE` the same four chains gave 1327, 1345, 1417 and 1356. A predictor whose splitting proportion has gone small is rarely proposed and so is hard to get back in, which is the known stickiness of the DART prior rather than a defect in this implementation — `update_s_param()` and `update_alpha_param()` were both re-derived against Linero's Dirichlet conditional and are correct, and `alpha_scale` defaults to the group count as it should. But it means a single chain can look much more settled than the posterior is, and the `marginaleffects` help page now says to run several.

Not changed: the defaults. `num_trees = 50` with the sparsity prior on is a deliberate configuration and the benchmark above is built on it. Whether 50 trees is too few for the DART prior to mix at is an open question worth measuring, and is the one item this round added to the To Do list.

### The control surface, reorganized

`bartisan_control()` had 29 arguments in no particular order, several of which existed only so a test could check that two code paths agree. It now has the same settings in three declared groups, stated in the description and marked in each `@param`: modeling decisions (`num_trees`, `gate`, `sparsity`, `k`, `bandwidth`, the chain lengths, `augment`, `x_transform`), advanced settings (everything from `gamma` to `num_print`), and three toggles that exist for internal validation (`block_eval`, `exact_quadratic`, `generic_accumulate`).

Three substantive changes came with it.

**`num_trees` takes a vector, one value per additive predictor.** A scalar is recycled, so the common case is unchanged. The engine stored forests as a rectangle -- `(iter * num_forest + h) * num_trees + t` -- and now stores them back to back with a per-forest offset, which is the only indexing change; the same substitution applies to the bandwidth matrix's columns and to `bartisan_predict`. The leaf scale divides by the square root of each forest's *own* tree count, so shrinking one forest leaves the prior on the sum it forms unchanged. `print()` and `summary()` say "2 forests of 50 and 10 trees" when the counts differ.

**`soft` is gone and `gate` decides both questions.** `gate = "hard"` (or `"step"`) gives the step functions of standard BART; `"smoothstep"`, `"smootherstep"` and `"logistic"` give soft rules and name the gate's shape. They were one decision pretending to be two: a hard rule has no gate shape to pick, and the old pair allowed `soft = FALSE, gate = "logistic"`, which had a test asserting that the second argument was ignored. `soft` survives as an internal field derived from `gate`, because `predict()` and the engine both need it.

**`sparsity` replaces four hyperparameters for the common case.** `TRUE` (the default) is Linero's (2018) DART prior, `FALSE` is a uniform prior over predictors, and `"none"`, `"weak"`, `"moderate"`, `"strong"` name four strengths. It sets `update_s`, `update_alpha`, `alpha_shape_1` and `alpha_shape_2` together; any of those supplied directly wins. The four knobs are a poor interface for what a caller wants to say, because `alpha / (alpha + P)` is Beta(a1, a2) and moving the selection pressure means moving two numbers in opposite directions at once.

### How many trees, and does the default depend on anything

Friedman function, n = 1000 train and 1000 test, p = 10, four chains of 500 draws after 500 warmup, held-out RMSE against the true regression function:

| Trees | Soft (smoothstep) | Hard | Soft seconds | Hard seconds |
|---|---|---|---|---|
| 5 | 0.286 | 1.149 | 5.4 | 4.4 |
| 10 | 0.281 | 0.682 | 5.1 | 4.1 |
| 20 | **0.270** | 0.558 | 6.0 | 4.3 |
| 50 | 0.284 | **0.521** | 9.3 | 5.6 |
| 100 | 0.289 | 0.531 | 15.1 | 8.3 |
| 200 | 0.319 | 0.510 | 28.2 | 11.0 |

**Soft rules need far fewer trees than hard ones**, and 200 -- the default in most BART packages -- is actively worse for them than 20. **The two want different counts**: hard rules are still improving at 200 where soft rules peaked at 20 and then degraded 12%.

That argues for a gate-dependent default, and the answer is still no, because point accuracy is not the only thing a tree count buys. On `MatchIt::lalonde`, four chains, average contrast on `treat`:

| Trees | `sparsity = TRUE` spread | `sparsity = FALSE` spread | P(contrast exactly 0), sparsity on |
|---|---|---|---|
| 10 | 45% | 93% | 0.39 |
| 20 | 122% | 35% | 0.32 |
| 50 | 127% | **9%** | 0.20 |
| 100 | 127% | 10% | 0.23 |
| 200 | 113% | 7% | 0.18 |

So 20 soft trees costs 35% between-chain disagreement where 50 costs 9%, against a 5% gain in Friedman RMSE. **50 stays the default for both gates**, and the two curves above are documented so that someone optimizing for prediction can drop to 20 and someone using hard rules can raise towards 200.

**Where the default should arguably vary is the number of forests, and there the answer is to document rather than to default.** `location_scale()`, n = 1000, smooth mean and log-linear standard deviation:

| `num_trees` | Seconds | Mean RMSE | Log-SD RMSE | Log score |
|---|---|---|---|---|
| `c(50, 50)` | 14.5 | 0.092 | 0.050 | -1188 |
| `c(50, 20)` | 8.3 | 0.094 | 0.047 | -1188 |
| `c(50, 10)` | 5.9 | 0.093 | 0.051 | -1188 |
| `c(50, 5)` | 4.9 | 0.094 | 0.046 | -1187 |
| `c(20, 5)` | **2.9** | **0.084** | **0.041** | **-1184** |
| `c(50, 1)` | 3.9 | 0.096 | 0.048 | -1188 |

A Gaussian fit on the same data is 1.4 s, so `c(50, 50)` is 10x a Gaussian fit and `c(50, 5)` is 3.5x, at the same accuracy to three decimals. Even one scale tree holds up here -- but that is because this truth's log standard deviation is linear in one predictor, and how many trees a variance surface needs depends on how complicated it is. Under-parameterizing it silently would show up as intervals that are wrong, which is the thing `location_scale()` exists to get right. So it is documented in `?bartisan_control` with the table, not made the default. Making it one is a one-line change to `resolve_num_trees()` if that judgment is ever revisited.

**A correction this produced.** The `marginaleffects` help page said a larger `num_trees` removes the atom at zero "almost entirely". It does not: with the sparsity prior on, the contrast was exactly zero in 20% of draws at 50 trees and 18% at 200. The earlier single-chain measurement that suggested otherwise was one lucky chain, which is exactly the failure mode the same page warns about. Corrected.

### Infinite BART (Battiston and Luo 2025): assessed, not implemented

arXiv 2511.20087 proposes `Y_i = sum_k W_ik g(X_i; T_k, mu_k) + eps_i` with `W` an n-by-infinity binary matrix under a three-parameter Indian buffet process prior. Two claimed features: the number of trees is learned, and each observation uses only a subset of the trees, which induces soft clustering with heterogeneous regression functions.

**It is implementable here, and the way in is neat.** Every accumulation and every `eta` commit in this engine goes through a node's index list and its membership weights, and the weights are already fractional because that is what a soft rule is. So `W_ik` is a multiplicative factor on the *root* membership weight of tree k: a masked observation then contributes nothing to any leaf's sums and receives nothing from the tree, with no change to the accumulators, the structure moves or the leaf draws. The row update needs `logdens_unit(i, eta_i)` with and without tree k's contribution, which the `Family` interface already exposes, so it would generalize to every family rather than just the Gaussian one, and it is O(nK) arithmetic per sweep. The dynamic column count is the only real work, and it can be avoided with the standard finite Beta-Bernoulli truncation at `K_max`, exactly as `dpm()` already truncates a Dirichlet process.

**It was not implemented, because the measurement does not support it.** Two reasons.

First, **the clustering cannot reach a new observation.** `W`'s prior does not depend on the predictors -- the paper says so in its discussion -- so at a new point the predictive mean is `sum_k E[W_k] g_k(x)`, a fixed re-weighted sum of trees. There is no per-observation tree selection at prediction time, by construction. Whatever the model gains has to come from `W` acting as an allocation device *during fitting*, not from clustering the test set.

Second, **the paper's headline gain is mostly a weak baseline.** Its clearest win is the clustered Friedman example of section 4.2: five groups of 40 observations, each with a Friedman regression function on a different window of five of nine predictors, group label never observed. Reported: classic BART 38.14, infinite BART 29.80, mean test MSE over ten 4:1 splits. Reproducing that design here, same n, same Beta marginals, same ten splits:

| Fit | Mean test MSE | sd over splits |
|---|---|---|
| paper's classic BART | 38.14 | — |
| paper's infinite BART | 29.80 | — |
| dbarts, 200 trees (its default) | 32.52 | 8.36 |
| dbarts, 10 trees | 33.46 | 8.14 |
| bartisan default (soft, sparsity on, 50 trees) | 31.65 | 7.48 |
| bartisan soft, `sparsity = FALSE`, 50 trees | **31.33** | 7.67 |
| bartisan hard, `sparsity = FALSE`, 200 trees | **31.29** | 7.72 |

Two independent modern BART implementations land at 31.3 to 32.5 where the paper's classic BART reads 38.14, so roughly 80% of the gap it reports closes without any of its machinery. What is left, 1.5 units, is a fifth of the between-split standard deviation and is on data drawn from a different seed, so it cannot be resolved without running their code. The sweep also shows `sparsity = FALSE` beating `sparsity = TRUE` at every tree count on this design, which makes sense: all nine predictors matter to some group, so a variable-selection prior is working against the truth.

**What would change the decision.** A covariate-dependent prior for `W` -- which the paper names as future work -- would make the clustering reach new observations and would turn this into a dependent Dirichlet process style conditional density model, which is a different and more interesting proposition. A direct comparison against their implementation on identical draws would settle the residual 1.5 units. Neither is cheap, and the finite-truncation prototype described above is the way to get the second if it is ever wanted.

### Infinite BART, implemented as a prototype and measured: it does not learn the tree count

The assessment above said the finite Beta-Bernoulli truncation was the cheap way to build this and that the root-weight trick was the way in. Both held. The prototype lives in a scratch copy of the package, not here, because the measurements say it should not ship.

**What was built.** `Tree` gained a `mask`, a 0/1 column of the weight matrix `W`. A masked observation is simply absent from the root's index list, so the tree never sees it: no accumulator changes, no weight vector for hard rules, and `reseat_tree()` rebuilds the root and calls the existing `rebuild_support()` when a column changes. `eval_live()` walks the live nodes to get what a tree *would* give an observation it is switched off for, which is what the row update needs. The update precomputes that K-by-n table once, does O(nK) arithmetic, and reseats only the columns that moved, so it costs about the same order as a sweep. `pi_k ~ Beta(a/K, 1)` conjugately, and `a` by slice sampling. Because the row of `W` for a new observation is unknown, the stored trees are scaled by `pi_k` on the way out, which makes `predict()` the plug-in predictive mean; in-sample `eta` uses the realized weights. About 200 lines across `node.h`, `node.cpp`, `mcmc.cpp` and `model.cpp`.

**The implementation is right.** Fixing the concentration at 1e4 drives `W` to all ones -- 49.8 of 50 trees per observation, mean `pi` 0.995 -- and the fit reproduces plain bartisan: test RMSE 0.431 against 0.419, Friedman n = 400, p = 10. That dense limit is the check that matters, because the model is *defined* to reduce to BART there.

**It does not select a small number of trees.** Friedman, n = 300, p = 30, six replicates, RMSE against the true regression function:

| Fit | RMSE (sd) | seconds | active trees |
|---|---|---|---|
| bartisan, 5 trees | 0.662 (0.166) | 0.1 | -- |
| **bartisan, 10 trees** | **0.507 (0.029)** | 0.1 | -- |
| bartisan, 20 trees | 0.529 (0.070) | 0.2 | -- |
| bartisan, 50 trees | 0.551 (0.108) | 0.5 | -- |
| bartisan, 200 trees | 0.610 (0.083) | 2.0 | -- |
| ibp, truncation 50 | 0.884 (0.239) | 1.3 | 46.5 |
| ibp, truncation 200 | 1.025 (0.110) | 5.1 | 175.3 |
| ibp, truncation 200, concentration fixed at 2 | 0.786 (0.088) | 4.4 | 59.9 |

It keeps almost every tree it is given, and it costs a factor of two in RMSE and a factor of fifty in time against the 10-tree fit that wins.

**And the count it reports depends on where the chain starts.** Truncation 200, four chains per row, same data:

| Setting | active trees (range over chains) | concentration | RMSE |
|---|---|---|---|
| start dense, Gamma(1, 1) prior | 161.0 (150--184) | 48.4 | 1.006 |
| start from the prior, Gamma(1, 1) | **79.0 (74--82)** | 17.3 | 0.883 |
| start dense, Gamma(0.05, 0.01) as in the paper | 200.0 (200--200) | 2428 | 0.717 |
| start from the prior, Gamma(0.05, 0.01) | **114.2 (106--128)** | 29.5 | 1.165 |
| start from the prior, concentration fixed at 1.2 | 11.0 (9--12) | 1.2 | 1.175 |

Same prior, same data, different starting point: 161 against 79, and 200 against 114, with ranges over four chains that do not overlap. The number of trees is not being learned; it is being remembered. The one row that does pick a small number is the row where the concentration was fixed by hand -- which is choosing the tree count, one level of indirection away -- and it is the worst fit in the table, 1.175 against 0.507 for a 10-tree bartisan fit. At a matched effective tree count the per-observation subsetting costs a factor of 2.3, because on a homogeneous problem it is noise.

Two other things the table says. With the paper's own weak concentration prior the chain runs to 2428 and every tree is active, which is standard BART -- and that row has the *best* RMSE of the five, which is the model telling you what it wants. And bartisan's leaf-scale warning fired on three of the sparse fits, correctly: the predictor is weakly identified when each observation sees a random subset of the trees.

**The variable-importance claim, which is the paper's headline, is better served by the sparsity prior already here.** Friedman, n = 300, p = 30, six replicates; separation is the smallest importance among the five real predictors minus the largest among the 25 noise ones, so positive means a clean split:

| Fit | Separation | Real predictors in the top five |
|---|---|---|
| bartisan 200 trees, sparsity off | +0.0042 | 4.8 / 5 |
| bartisan 200 trees, DART | +0.0426 | 5.0 / 5 |
| bartisan 50 trees, DART | +0.0433 | 5.0 / 5 |
| bartisan 10 trees, sparsity off | +0.0124 | 4.5 / 5 |
| **bartisan 10 trees, DART** | **+0.0836** | 5.0 / 5 |
| ibp, truncation 200 | **-0.0013** | 4.5 / 5 |
| ibp, truncation 200, DART | +0.0784 | 5.0 / 5 |

The paper's premise checks out: 200 trees with no sparsity prior barely separates the real predictors from the noise. But the Indian buffet process is not the fix. On its own it makes the separation *negative*, worse than plain BART; the two rows where it looks good are the rows where DART is on, and DART reaches the same place at a tenth of the cost without it.

**Verdict.** The prototype answers the question it was built for. Learning the number of trees is not what this model does: it replaces one choice with two -- a truncation and a concentration prior -- and returns an answer that depends on initialization. Kept in the scratch tree in case a covariate-dependent prior for `W` ever makes the clustering reach new observations, which is the change that would make the model a different proposition.

### McCartan and Huang (2026): their ablation replicated here, and where it stops holding

*Seeing the Forest for the Trees: The Gaussian Process Limit of BART* (arXiv 2607.28844) proves that a symmetric-tree BART prior converges weakly to a Gaussian process as the number of trees goes to infinity, derives the kernel, and shows that ridge regression on *random tree features* -- leaf indicators from trees drawn from the prior and never updated -- attains minimax-optimal rates depending only logarithmically on the covariate dimension. The empirical claim underneath it is an ablation: once the number of trees is large, neither Bayesian averaging, nor learning the tree structure, nor asymmetric trees does much for out-of-sample R-squared.

**Nothing in the package was changed on the strength of this.** What follows is the replication and what it suggests for later.

**Setup.** Four datasets -- `airquality` (n = 111, p = 5), `MASS::Boston` (506, 13), a 1200-row sample of `ggplot2::diamonds` (9), and Friedman with n = 500 and p = 30, so 25 irrelevant predictors. Four 75/25 splits each, 400 draws after 400 warmup, predictors mapped through the training ECDF. Ablation (a) is a single final draw against the full posterior mean, which is a cruder version of theirs -- they condition on the final tree structure and integrate the leaves. Ablation (b) is random tree features against full BART, with the tree structures drawn from bartisan's own branching prior and the ridge penalty by leave-one-out. Their ablation (c), symmetrized trees, was not run; their own answer there is "no effect".

Differences in R-squared, averaged over the four datasets, negative meaning the ablated model is worse:

| Trees | Bayes, hard | Learning, hard | Bayes, soft | Learning, soft |
|---|---|---|---|---|
| 5 | -0.026 | -0.419 | +0.003 | -0.319 |
| 20 | -0.037 | -0.283 | -0.012 | -0.158 |
| 75 | -0.044 | -0.158 | -0.058 | -0.094 |
| 200 | -0.060 | -0.062 | -0.022 | -0.032 |
| 500 | -0.050 | -0.040 | -0.054 | -0.022 |

**Both of their findings replicate.** Ablating Bayesian averaging costs a small amount that does not depend much on the tree count. Ablating tree learning costs a great deal at five trees and almost nothing at five hundred. The shape is theirs.

**Their section 6 conjecture about soft varieties is confirmed.** They speculate that adapting the *type* of random feature to the data "may yield improved performance at a smaller number of trees". Random features built from bartisan's smoothstep gate against the same features built from hard splits, mean R-squared over the four datasets:

| Trees | Hard features | Soft features | Difference |
|---|---|---|---|
| 5 | 0.350 | 0.507 | **+0.157** |
| 20 | 0.549 | 0.693 | **+0.144** |
| 75 | 0.697 | 0.778 | +0.080 |
| 200 | 0.794 | 0.840 | +0.045 |
| 500 | 0.812 | 0.845 | +0.033 |

A soft gate is a better random basis, and exactly as they guess, the advantage is largest where the trees are fewest. The same thing shows up in the ablation table: the learning gap closes faster under soft rules at every tree count, because the prior-drawn soft basis is already closer to what learning would have produced.

**Where their conclusion stops holding is the other open question they name.** Their last paragraph asks whether the hierarchical variable-selection prior of Linero (2018) can be approximated by a penalty on random-feature coefficients. Measured on the Friedman design with 25 irrelevant predictors:

| Trees | Soft, no sparsity prior | Soft, DART | Soft random features |
|---|---|---|---|
| 5 | 0.921 | **0.950** | 0.168 |
| 20 | 0.950 | **0.958** | 0.436 |
| 75 | 0.946 | **0.957** | 0.669 |
| 200 | 0.937 | **0.958** | 0.877 |
| 500 | 0.925 | **0.958** | 0.874 |

Three things. DART is **flat in the tree count** -- 0.950 to 0.958 from five trees to five hundred -- where plain soft BART peaks at twenty and then decays. Random features never catch it: the gap is still 0.084 at five hundred trees and has stopped closing. And on the three datasets where most covariates matter, the same gap is 0.02 to 0.04 by five hundred trees, which is their result. So **"tree learning does not matter once T is large" is conditional on the covariates mostly mattering.** When they do not, what is being learned is which variables to split on, and a basis drawn from a uniform prior over predictors cannot represent that however many features it has.

That also explains why this package's defaults do not move. bartisan's configuration is soft rules with DART at fifty trees, and at that point on the curve the learning ablation still costs 0.09 to 0.16, not 0.02.

**What is worth following up, in order.**

- [ ] **Soft random tree features as a fast approximate fit.** The R prototype is about sixty lines and reached 0.840 average R-squared at 200 features against full soft BART's 0.872, in a fraction of the time. Two uses: a `random_features()` estimator for when a fit is needed inside a loop, and, more interestingly, a warm start for the MCMC -- which is the existing grow-from-root To Do item arrived at from a better direction, since these features come from the prior and cost one ridge solve.
- [ ] **Random tree features for the non-Gaussian families.** Their section 5.2 point is that random features slot into any linear predictor. Here that would mean a penalized GLM on the feature matrix, which reaches every family the package has without a sampler. Whether the uncertainty holds up outside the Gaussian case is open; their section 5.3 evidence is Gaussian only.
- [ ] **A sparsity-aware feature draw.** Drawing the splitting variable from the DART proportions of a short pilot run, rather than uniformly, is the obvious way to give random features the one thing the measurement above says they lack. This is their closing question and the table gives it a concrete target: 0.874 to beat 0.958 on Friedman with p = 30.
- [ ] **Reconsider whether `sigma_mu` should be tuned rather than drawn.** Their figure 3 bottom row shows the ablation patterns become much less variable across datasets once the leaf prior variance is tuned by cross-validation, and they flag incorrect tuning of it as the reason several datasets misbehave. bartisan draws it under a half-Cauchy, which is a third option neither of them tested; whether the drawn version lands where the tuned one does is a cheap thing to check and would say something about the leaf-scale warning this package emits.

**Not suggested by any of this:** changing the tree-count default, changing the gate default, or turning the sparsity prior off. The ablation's message is that computation spent on structure learning has diminishing returns at large T, and this package is not at large T -- it is at fifty trees with a basis and a prior that both make the learning worth more, not less.

### `Gamma_shape()` removed, `Gamma()` masked so that the default link is log

`Gamma_shape()` existed for one reason: `stats::Gamma()` has no slot to carry a fixed shape, and the package's convention is that a family which draws a nuisance parameter also lets you fix it (`negbin(theta =)`, `ordbeta(phi =)`, `dpm(alpha =)`). Realistically nobody fitting BART knows a gamma shape, so the convention was not worth a second family function and it is gone. The shape is drawn, as it always was; the engine still supports holding it fixed and nothing exposes that.

Removing it exposed something worse, which is why this entry exists at all. **`stats::Gamma()` defaults to `link = "inverse"`, and the two functions therefore differed in the default link, not only in the argument.** The inverse link is the worst case for this sampler: its inverse sends a negative predictor to a negative mean, whose log is not a number, so the proposal is rejected. `compose_link()` already carried a comment saying a default `Gamma()` fit produces "dozens of them", but nothing surfaced it and three documentation tables listed `Gamma()`'s link as `log`, which was true of `Gamma_shape()` and false of `Gamma()`. Measured on 600 observations and 50 trees, fitted mean against the truth:

| Call | Seconds | RMSE |
|---|---|---|
| `stats::Gamma()` -- inverse link | 7.2 | 0.664 |
| `Gamma("log")` | 3.8 | 0.606 |

So `Gamma()` is now exported from this package with `link = "log"`, which **masks `stats::Gamma()`**. It is otherwise the same function -- it returns `stats::Gamma(link)` unchanged, so a link name, a `link-glm` object and `glm()` all still work. The cost of the mask is that `glm(y ~ x, family = Gamma())` gets the log link while bartisan is attached; that is documented in three places and `stats::Gamma()` still reaches base R's default. `Gamma("log")` and the old `Gamma_shape()` produced bit-identical draws from one seed, which is what confirms the two were the same model.

Separately, and generally rather than for the gamma alone: **`bartisan()` now says when a composed link's inverse does not cover the whole additive predictor.** `warn_restricted_link()` evaluates the caller's `linkinv` on a grid and checks it against the domain the engine's own link needs -- positive for `log`, the unit interval for `logit`. It fires for `Gamma("inverse")`, `Gamma("identity")` and `poisson("identity")`, and stays quiet for compiled links and for composed links that do cover the line, such as `binomial("cauchit")`.

One incidental fix: `test-bartisan.R` matched the inferred-family message on the literal `"using"`, and that message had been sentence-cased to "Using ...", so the regex silently stopped matching -- testthat reports a non-matching regexp with the same wording it uses for no message at all, which is what made it look like the message had disappeared. It now matches on `"family = "`, which no capitalization rule touches.

### gaussian(), dpm() and ordinal() on a numeric response: dpm dominates

The question was whether `dpm()` should be recommended over `gaussian()` in general, and whether `ordinal()` is a legitimate choice for a continuous outcome. Both hold up. 200 training and 200 test observations, 50 trees, 500 draws after 500 warmup, four replicates, errors centered so that every family is estimating the same conditional mean:

| Errors | `gaussian()` | `dpm()` | `ordinal("probit")` | `ordinal("logit")` |
|---|---|---|---|---|
| normal | 0.263 / -280 | **0.252 / -279** | 0.263 | 0.294 |
| t3 | 0.264 / -281 | **0.213 / -261** | 0.307 | 0.244 |
| skewed | 0.183 / -237 | **0.145 / -212** | 0.177 | 0.180 |
| bimodal | 0.356 / -334 | **0.199 / -261** | 0.353 | 0.415 |
| heteroskedastic | 0.281 / -315 | 0.289 / **-312** | **0.274** | 0.282 |

RMSE against the true regression function, and for the two continuous families the held-out predictive log score. At 1000 observations, with the full grid:

| Errors | `gaussian()` | `dpm()` | `ordinal("probit")` | `ordinal("logit")` |
|---|---|---|---|---|
| normal | 0.146 / -1433 | 0.143 / -1434 | **0.139** | 0.152 |
| t3 | 0.135 / -1425 | **0.102 / -1247** | 0.133 | 0.127 |
| skewed | 0.089 / -1216 | **0.073 / -1020** | 0.091 | 0.097 |
| bimodal | 0.165 / -1660 | **0.049 / -1078** | 0.150 | 0.175 |
| heteroskedastic | 0.171 / -1596 | 0.150 / -1563 | 0.145 | **0.141** |

Level on normal errors to within one log point, ahead by 178 and 196 on heavy tails and skewness, and ahead by **582** on bimodal errors with a third the RMSE. One refinement over the smaller sample: heteroskedasticity was a wash at n = 200 and at n = 1000 `dpm()` is ahead of `gaussian()` by 33 log points, with `ordinal()` ahead of both on error. `location_scale()` is still the family that actually finds the pattern.

**`dpm()` does not pay for its flexibility.** On normal errors, where `gaussian()` is exactly right, it came out slightly ahead on both measures. That is what makes it a default rather than a specialist tool, and it confirms the reading that `gaussian()`'s remaining advantages are not statistical. They are: prior weights, which `dpm()` refuses (verified -- `gaussian()`, `ordinal()` and `location_scale()` all take them and `dpm()` errors); an identified additive predictor, since `dpm()` identifies only the sum of the fit and the error mean; one interpretable `sigma`; and 1.4 times the speed at a thousand observations.

**`ordinal()` on a continuous outcome is a real method**, and the cutpoint structure is why. Every distinct value becomes a category, the cutpoints absorb the marginal distribution, and the forest explains only the ordering, so nothing is assumed about the error and the model for the cumulative probability is invariant to a monotone transformation of the response. This is `rms::orm()` with a forest in place of the linear predictor. It never won by much, but it had the lowest error and the only above-nominal coverage on the heteroskedastic row -- the one setting where the other two are misspecified -- which is what a model with no error distribution should do.

**Bin the outcome, and bin it hard.** One cutpoint per distinct value means n cutpoints. t3 errors, n = 1000:

| Cutpoints | RMSE | Coverage | Seconds |
|---|---|---|---|
| 10 bins | 0.118 | 0.96 | 4.8 |
| **25 bins** | **0.115** | 0.97 | **4.6** |
| 50 bins | 0.123 | 0.97 | 5.6 |
| 100 bins | 0.152 | 0.95 | 6.1 |
| 250 bins | 0.161 | 0.96 | 11.8 |
| every value | 0.134 | 0.97 | 73.3 |
| `gaussian()` | 0.133 | 0.97 | 2.7 |
| `dpm()` | **0.096** | 0.97 | 4.2 |

Twenty-five bins is sixteen times faster than no binning *and* slightly more accurate, because a cutpoint vector with a thousand weakly-identified entries is worse conditioned than one with twenty-five. Bimodal errors put the optimum at fifty bins (0.148); anywhere from ten to fifty is fine and the choice inside that range hardly matters.

Timings in these tables were taken with other jobs on the machine, so read the ratios within a table rather than the absolute seconds across tables -- the same unbinned cell read 73 seconds in one run and 185 in another under heavier load.

Documented in `vignette("families")` with a head-to-head section and a new "A continuous outcome as ordinal" section, in the `bartisan-families` help page, in the README, and in `NEWS.md`.

### dpm reports the conditional mean on the predictor, and is the numeric default

Two changes that go together, and the first is what makes the second reasonable.

**The reporting chart.** Nothing in the DPMBART model forces the error mixture to be centered, so the sampler works in a chart where only the *sum* of the predictor and the error mean is identified and each piece alone wanders. That was documented as a limitation and it was the one respect in which `dpm()` was harder to use than `gaussian()`. It is now a reporting question rather than a modelling one: `DPMFamily::report_shift()` returns minus the mixture's mean, so the recorded predictor moves up by it, `mixture_flat()` reports component means with the same amount taken out, and the chart the draw is recorded in has the mixture at mean zero and the whole conditional mean on the predictor. This is exactly the device the ordinal families use for their cutpoints, and `model.cpp` already distributed a shift across the recorded leaf values, so the stored forest still replays to the reported predictor.

Measured on a skewed example, n = 600:

| | before | after |
|---|---|---|
| sd of the predictor's level across draws | 0.750 | **0.021** |
| bias of `type = "link"` against the truth | 0.144 | **-0.021** |
| `type = "link"` vs `type = "response"` | differ by the drift | identical |
| log likelihood rebuilt in R vs the sampler's | 4.5e-13 | **3.4e-13** |

That last row is the check that matters and it caught a real bug on the first attempt. Shifting the atoms is not enough: the predictive density's *new-component* term is the baseline `G_0`, which lives on the raw chart, so after centering the R-side reconstruction was evaluating it in the wrong place and the invariant broke to 0.47. The fix is to report the shift, which `aux` now does as **`center`** in place of `error_mean` -- the raw mixture's mean, relabelled as the bookkeeping quantity it is rather than an estimate of anything, and used by `dpm_predictive()` to place the baseline term and by the posterior predictive sampler to place a fresh component. The error mean in the reported chart is zero by construction, and `error_density()` is centered: integral 1, mean 4e-04, sd matching `error_sd`, skewness 1.34 against a true 1.63 on a gamma error.

**The default.** A numeric response now gets `dpm()` rather than `gaussian()`. The comparison in the previous entry is the argument: `dpm()` matches `gaussian()` when the errors really are normal and beats it, sometimes by a factor of three, when they are not, so there is no error distribution on which the old default was the better choice.

With one boundary. A mixture cannot separate an error distribution from a mean when the response takes a handful of values, so a numeric response with fewer than ten distinct values keeps `gaussian()`. Rounding a continuous response onto k equally spaced levels, n = 600, held-out RMSE:

| Levels | `gaussian()` | `dpm()` | Ratio |
|---|---|---|---|
| 3 | 0.175 | 0.662 | **3.78** |
| 5 | 0.150 | 0.151 | 1.00 |
| 8 | 0.138 | 0.141 | 1.02 |
| 12 to 80 | ~0.14 | ~0.14 | 0.99 to 1.03 |
| continuous | 0.136 | 0.138 | 1.01 |

The break is between three and five, so ten is a conservative line, and in the range where the two tie the tie goes to the family that also takes weights. One case not caught by a distinct-value count: a *clamped* five-level scale, with mass piled on both end values, gave 0.122 against 0.224. Point masses at the boundary are the shape a mixture handles worst, and a distinct-value guard does not see them -- `ordinal()` is the right family there and the vignette says so.

Since `dpm()` refuses prior weights, a weighted fit with no family named is an **error** naming the alternatives rather than a silent substitution: dropping the weights and swapping the family are both defensible and only the caller can say which was meant.

### augment is a sampling setting, and what it would take to extend it

`augment` moved from the modeling group to the advanced group in `?bartisan_control`. It belongs there: a rewriting targets exactly the same posterior as the direct likelihood, so nothing about the model changes and what changes is how fast the sampler gets there. The default is the set of rewritings measured to pay, and nobody should have to think about it. Agreed with, not argued against.

**Which families could gain, and what it would take.** The families whose leaf target is `TARGET_GENERAL` with no augmented counterpart are the three accelerated failure time families, `ordbeta()`, and the second predictor of `location_scale()`. Everything else is already at the best form available -- `gaussian()` and `dpm()` are quadratic, `poisson()` and `Gamma()` are the exponential form, and the binomial, ordinal, multinomial, negative binomial and zero-inflated families all have augmented counterparts already.

The prize is worth stating first, because it is large. Timed on 1000 observations, 50 trees, 300 draws after 300 warmup, 31% censoring, against a Gaussian fit on the same design at 3.2 s:

| Fit | Seconds | Against `gaussian()` |
|---|---|---|
| `gaussian()` | 3.2 | 1.0x |
| `weibull_aft()` | 17.9 | 5.6x |
| `loglogistic_aft()` | 21.9 | 6.8x |
| `lognormal_aft()` | 31.9 | **10.0x** |
| `location_scale()`, `c(50, 50)` | 13.2 | 4.1x |

So a general target costs five to ten times a quadratic one here, which is the same order as the 14x that `ordinal("probit")`'s augmentation was worth. Concretely:

- **`lognormal_aft()` is censored normal regression.** Impute the censored log-times from a normal truncated below at the observed time; conditional on the completed data the likelihood is a plain Gaussian, so the target is **quadratic** and the leaf draw is closed form with acceptance one. `truncated_normal_between()` is already in `utils.h`, the augmented-family pattern is established in four places, and `sigma` draws from its conditional given complete data. This is the largest single performance win available in the package and the cheapest of the three to write.
- **`loglogistic_aft()`** is the same imputation from a truncated logistic, followed by the Pólya-Gamma step that `ordinal("logit")` already uses -- a logistic variate is a normal whose precision is Pólya-Gamma(2, |r|). Two augmentations composed, both already present, and the target lands quadratic.
- **`weibull_aft()`** imputes from a truncated Gumbel, after which the log density is `-(t - eta)/sigma - exp(-(t - eta)/sigma)`: the *exponential* form rather than the quadratic one, and with rate `1/sigma` rather than 1.
- **`ordbeta()`** has no clean route. Its beta component's log density is a constant times `logit^-1(eta)`, and neither the Pólya-Gamma identity nor a latent normal makes a sum of expits quadratic. Not worth pursuing.

**The one change that unlocked two of these was not an augmentation at all.** `TargetForm` had `TARGET_EXP_UP` and `TARGET_EXP_DOWN`, and `Family::exp_sign()` returned +1 or -1 -- the machinery was hardcoded to `exp(±eta)`. Both remaining cases needed a *rate*. Generalizing `exp_sign()` to `exp_rate()` reached both; see the entry below.

A caution on measurement, learned here: **`exact_quadratic = FALSE` is not a proxy for a general target.** On a Gaussian fit it costs 1.16x and on a Poisson fit nothing at all, because the underlying target really is quadratic or exponential and the iteration terminates immediately whichever path it takes. The five-to-ten-times figures above are direct comparisons between families, which is the only way to see it.


## The subsample, and the two defects it stopped hiding

The shipped `rhc` is a random 1500 of the 5735, and `rhc` and `death` are integer 0/1 rather than factors. Neither the roxygen block nor `data-raw/rhc.R` mentions build time as the reason for the subsample; the documented reason is what it is, a random sample, and there is no case for advertising to a reader that the analysis they are about to run is slow.

The 0/1 coding is documented for the reason that motivated it: a contrast between two levels is then a single number rather than one row per level of a factor.

### `~pairwise` worked once binomial stopped being a categorical family

The expected fix was the 0/1 recoding, and it was not the fix. The duplicated rows never came from the outcome's *type*; they came from `me_type()` listing binomial alongside `ordinal`, `multinomial` and `mnp`, which puts the fit on the probability scale and returns one group per outcome category. Removing binomial from that list is the change. Its two probabilities sum to one, so reporting both gave every estimand twice as mirror images where `glm()` gives one row, and `hypothesis = ~pairwise` then paired across outcome levels as well as subgroups. `type = "prob"` still asks for both columns for anyone who wants them.

With that removed, `avg_comparisons(fit, variables = "rhc", by = "card", hypothesis = ~pairwise)` returns the single clean row it should, and the vignettes no longer name `"b4 - b2 = 0"`.

### The `. - days` defect was in the *stored* formula

The vignettes now write model formulas out in full, which is closer to how a reader would write one and removes the need to subset the data. That alone sidesteps the defect, but the defect was also worth fixing, and the diagnosis in the entry above was wrong about where it lived. `model.frame()` keeping `days` is real and shared with `glm`, but it is not what collapsed the table: `insight::find_formula()` reads the formula stored on the fit, and from the *unresolved* `death ~ . - days` it concluded the model's one predictor was `days`, the variable the formula removes.

`bartisan()` now resolves `.` before storing, via `terms(formula, data = data) |> update(. ~ .)`, and stores the resolved formula as well as using it. Left alone when the formula has random-effect terms, where `.` would expand over the grouping variables too.

### Variable importance at n = 1500, and the claim that had to be withdrawn

At 5735 rows the three noise controls sorted below every real predictor. At 1500 they do not: one lands above five of the real predictors and another sits in the middle of them. The vignette says so. The lesson is better than the one it replaced, because it is the lesson the section exists to teach: the ranking alone was never the finding, the gap against a variable known to be noise is, and here that gap says the top of the table is real and the bottom half is not distinguishable from random numbers.

Every other number in the five RHC vignettes was re-read off a fresh render rather than scaled, and the prose corrected against it.

## A one-sided formula reported a missing object named `.`

`bartisan(~ x1 + x2, data = d)` gave `object '.' not found`. The `.` expansion
added for the `- days` defect is what produced it: `update(~ x1 + x2, . ~ .)`
returns `. ~ x1 + x2`, inventing a response named `.`, and the model frame then
went looking for a variable by that name. Before the expansion existed the
message was no better, since a formula with no response fell through to the
family default and came back as "The dpm family requires a numeric response".

Fixed at the top of `bartisan()` rather than inside the expansion:
`arg::arg_formula(formula, one_sided = FALSE)` already says
"`formula` must be a two-sided formula", so the case never reaches the
expansion at all. Tested for both `~ x1 + x2` and `~ .`.

## Varying-coefficient BART works through `custom_family()`, and the interface is open

Full write-up in `_dev/bcf-interfaces.md`; the runnable proof is
`_dev/bcf-proof.R` and the mixing follow-up `_dev/bcf-mixing.R`.

`g(mu_i) = f_0(x_i) + z_i * tau(x_i)` with both surfaces forests, written as a
`custom_family()` with `num_predictors = 2` and the chain rule supplied as the
analytic derivatives. It samples, and it recovers `tau(x)` better than the
alternative available today: correlation with truth 0.986 against 0.969, RMSE
0.150 against 0.233, ATE 0.964 against a truth of 0.981. Binomial the same on
the logit scale.

**The device, and the gap behind it.** A custom likelihood is handed `y` for the
rows reaching a leaf and is not told which rows they are, so a treatment vector
held outside cannot be aligned. The proof carries `z` inside the response and
unpacks it, exactly for a binary outcome and by a scale shift otherwise. The
reason that was necessary is worth recording on its own: **`custom_family()`
cannot see per-observation covariates.** Any hand-written censored or truncated
likelihood needs its indicator and hits the same wall.

**Two architectural findings, both checked rather than assumed.** `LinkFamily`
in `src/family.cpp` is already a decorator that transforms the predictor inward
and applies a chain rule to `d1` and `d2` outward, so a varying-coefficient
decorator is the same object with a different map, and it would work for every
family without per-family code. And because the map is linear in each `eta_j`, a
target quadratic in `mu` stays quadratic in each `eta_j`: `gaussian()`, `dpm()`
and every augmented family keep their closed-form leaf draws. The exponential
form does not survive, because `exp_rate()` returns one scalar per predictor and
a varying coefficient makes the rate vary by observation, so `poisson()` and
`Gamma()` would drop to `TARGET_GENERAL`.

**A mixing claim withdrawn.** The proof's rhat of 1.39 on the prognostic forest
looked like the weak separation of the two surfaces that the BCF literature
warns about. It is not: on the same data, the compiled `gaussian()` with `z` as
an ordinary predictor and no varying coefficient anywhere mixes the same
(`sigma_mu` 1.36, `eta` 1.35), and a hand-written one-forest Gaussian is worse
(1.54). The varying coefficient does not degrade mixing.

## `sigma_mu` mixes poorly across every design tried

Fell out of the above and is not related to it. `_dev/sigma-mu-mixing.R`,
compiled `gaussian()`, n = 1000, four chains, 750 saved each:

| Design | `sigma_mu` rhat | `sigma_mu` ess_bulk | `eta` rhat | `eta` ess_bulk |
|---|---|---|---|---|
| Friedman, high signal | 1.37 | 8.9 | 1.37 | 9.2 |
| Friedman, low signal | 1.22 | 13.3 | 1.06 | 70.5 |
| pure noise | 1.85 | 5.8 | 1.49 | 7.5 |
| one linear predictor | 1.24 | 12.5 | 1.01 | 538.4 |
| Friedman, `update_sigma_mu = FALSE` | -- | -- | 1.29 | 10.5 |

The `eta` column is a maximum over 1000 per-observation values, so it is
inflated by selection and should be read the way `vignette("diagnostics")`
already reads it, as quantiles rather than a max. `sigma_mu` is one scalar per
forest, so no such excuse applies: ess_bulk under 15 out of 3000 draws, on every
design including the package's own benchmark.

Investigated below.

## The ramp left the leaf scale frozen for most of warmup

Found while investigating the above. It is a real bug and it is **not** the
cause; see the next entry.

`sigma_mu_ramp` holds `sigma_mu` at a fraction of its target over the first part
of warmup, per Linero (2025) Remark 2, and switches its update off while doing
so. Nothing switched the update back on when the ramp ended. The ramp block was
the only thing that ever touched the flag, and the restore sits after the warmup
loop, so with the default `sigma_mu_ramp = 0.25` the scale was frozen for the
remaining 75% of warmup and took its first draw at the first *retained*
iteration. On a 500-iteration warmup it jumped from 0.212 to about 0.8 in one
step, and the trees then spent the sampling phase equilibrating to a leaf scale
four times larger than the one they had been fitted under: the mean of the first
50 retained draws was 0.756 against 0.880 for the last 50.

Fixed with an `else if (iter == num_ramp)` branch that hands the scale back for
the rest of warmup. It fires only when the ramp is shorter than warmup, so
`sigma_mu_ramp = 1` still relies on the restore below the loop, which is why
that restore was moved out of the loop in the first place. After the fix the
systematic opening drift is gone. Two regression tests.

## Why `sigma_mu` mixes badly: three hypotheses, two refuted, and an honest stop

`_dev/sigma-mu-cause.R`. The ramp fix changed the mixing numbers not at all
(Friedman `sigma_mu` rhat 1.370 before, 1.371 after; pure noise 1.848 before,
1.888 after), so the cause is elsewhere.

n = 600, Friedman, four chains:

| trees | burn/save | acceptance | `sigma_mu` rhat | `sigma_mu` ess | `sd(eta)` ess | `eta` rhat |
|---|---|---|---|---|---|---|
| 10 | 750 | 0.89 | 1.230 | 12.4 | 2501 | 1.496 |
| 50 | 750 | 0.94 | 1.322 | 10.1 | 2768 | 1.224 |
| 200 | 750 | 0.97 | 2.015 | 5.4 | 1687 | 1.108 |
| 50 | 4000 | 0.95 | 1.418 | 8.1 | 8023 | 1.197 |

**Not the sampler.** The step is an independence Metropolis proposal from
`Gamma(n/2 + 1, 2/sse)`, the conditional posterior of the precision under a flat
prior, corrected by the half-Cauchy term. Read against the definitions the
proposal and the correction are both right, and the acceptance rate is 0.89 to
0.97, so nothing is sticking.

**Not the forest's mixing.** `sd(eta)`, a global aggregate of the fitted surface
and the same kind of quantity, has an ess in the thousands while `sigma_mu` has
single digits. `sigma_mu` is not inheriting anything.

**Not slow mixing either.** Its ess does not grow with the run: 10.1 at 750
saved draws and 8.1 at 4000, a 5.3-fold increase in length. A quantity mixing
slowly but ergodically would have gained roughly in proportion. Each chain
settles somewhere and stays.

**What is actually going on, as far as the evidence reaches.** The per-chain
means tell it. At 50 trees: 1.083, 1.083, 1.189, 1.085. At 200 trees: 0.450,
0.465, 0.569, 0.526. The within-chain distribution is very tight, because the
proposal is the flat-prior posterior over hundreds of leaf parameters and has a
coefficient of variation around `1 / sqrt(n_leaves)`, measured at 0.09. rhat is
the ratio of between-chain to within-chain variance, so a small but persistent
offset between chains against a very tight within-chain spread gives a large
rhat for a disagreement that is practically negligible. Adding trees tightens
the proposal further without shrinking the offset, which is exactly why more
trees makes `sigma_mu` worse while making `eta` better.

**A hypothesis that did not survive.** The obvious explanation for the offsets is
that an additive ensemble can trade leaf magnitude against tree structure
without changing its sum, so `sigma_mu` distinguishes configurations the
likelihood cannot. The test was whether chains with more splits show a smaller
`sigma_mu`. The correlation came out +0.90 at 50 trees and -0.92 at 200, on four
points each. That reverses sign and is four points, so it is noise and the
hypothesis is unsupported. **Why the chains settle at slightly different levels
is not established.**

## The other packages do it too, so the leaf scale is closed

`_dev/sigma-mu-others.R`. Same data, four chains, 750 saved draws each, every
row's rhat and ess computed by this package's own functions:

| Package | Quantity | Prior and sampler | rhat | ess |
|---|---|---|---|---|
| bartisan | `sigma_mu` | half-Cauchy, independence Metropolis | 1.188 | 14.8 |
| dbarts | `k` | chi hyperprior, slice | 1.122 | 21.8 |
| stochtree | `sigma2_leaf` | inverse-gamma, **conjugate Gibbs** | 1.163 | 17.1 |
| BART | `k` | not drawn: fixed at 2 | -- | -- |
| bartisan | `sd(eta)`, for contrast | -- | 1.000 | 2440 |

**stochtree settles it.** Its prior is conjugate, so its leaf scale is drawn
exactly from its full conditional. Nothing can mix better than an exact draw,
and it gets ess 17.1 out of 3000 with rhat 1.163 -- the same range as the other
two, on the same data. The sampler is therefore not the cause anywhere, which is
what the acceptance rate had already suggested for this package. Three
independent implementations with three different samplers show the same
behaviour, and the fourth avoids the question by never drawing the parameter,
which is itself a comment on how much anyone gets out of it.

So this is a property of the leaf scale in an additive tree ensemble, not a bug
in bartisan, and bartisan's numbers are unremarkable among its peers. Closed. Why
the chains settle at slightly different levels is still not established and no
longer worth establishing.

**What remains is presentational, and is a real problem.** `fit$rhat` puts a
`sigma_mu` row next to the `eta` rows at equal status, that row sits above any
conventional threshold on ordinary data including the package's own benchmark
design, and a reader following `vignette("diagnostics")` has no way to know it is
expected. Nothing about the fitted function is affected: on every design where
`eta` has signal to mix on, `eta` is fine while `sigma_mu` is not. Two things
worth deciding, neither done here because both change what the package reports:

- Drop `sigma_mu` from `fit$rhat`, or mark it as a hyperparameter whose rhat is
  not a reason to distrust the fit.
- Say so in `vignette("diagnostics")` either way, now that there is a
  cross-package answer to point at.

## Per-predictor splitting priors: `split_prior`

`bartisan_control(split_prior = c(x1 = 3, x3 = 0.5))`. Every predictor starts at
a weight of 1, the named ones take the value given, and the prior probability of
splitting on a predictor is its weight over the total, so on three predictors
that example gives 3/4.5, 1/4.5 and 0.5/4.5.

One weight per *term*, not per design-matrix column, so a factor is named once
and its levels share the weight. That is the same granularity the sparsity prior
already uses, and `make_group_probs()` already labels its columns by term, so the
weights match by name against `colnames(group_probs)`.

Implemented by fixing `Hypers::s_` at the normalized weights and turning
`update_s` and `update_alpha` off. The C++ change is one optional constructor
argument; an empty vector keeps the old uniform initialization, and a full one
also forces the two flags off inside the constructor, so the guarantee holds
whatever reaches it rather than depending on R having set them.

**It overrides `sparsity`, and warns only when asked for both.** The two answer
different questions: `sparsity` is for not knowing which predictors matter and
wanting the prior to find out, `split_prior` is for knowing something and wanting
it honored. Drawing `s` from a Dirichlet centered on the weights would answer
neither, so the weights are held fixed. Since `sparsity = TRUE` is the default
and therefore not a request, the warning fires only when the caller passed
`sparsity` explicitly, which `missing()` distinguishes.

**Verified against the specification.** On four pure-noise predictors, where
nothing in the data prefers any of them, the realized share of splitting rules
matches the weights: `c(x1 = 8)` gave 0.740 against 8/11 = 0.727, and
`c(x1 = 3, x3 = 0.1)` gave 0.577 / 0.201 / 0.199 / 0.023 against
0.588 / 0.196 / 0.196 / 0.020.

**A documentation claim corrected before it shipped.** The first draft said
`split_prior` leaves `prop_used` at 1 for every predictor, since nothing is
dropped from the forest. The test failed. A predictor can still miss out on a
rule in some draw of a small forest, and at 20 trees `prop_used` ran 0.89 to
1.00. The distinction that survives is better than the one I wrote: under fixed
weights the exclusion is sampling variation at a fixed probability and goes away
as trees are added, reaching 1 for every predictor at 50 and at 200, whereas the
sparsity prior on the same data at 50 trees left the *most heavily weighted*
predictor out of 41% of draws. Documented as measured.

Errors rather than silence on a name the model does not have, since a weight is a
claim about a particular predictor and a typo would otherwise change nothing and
say nothing. Zero and negative weights refused, with the reason: a zero would
forbid splitting rather than discourage it, and dropping the predictor from the
formula is the way to say that.

## Zero is a legitimate `split_prior` weight

Was refused with a message suggesting the predictor be dropped from the formula
instead. Wrong: a weight of zero is a use for the argument, not a mistake in it.
It holds the predictor out of every tree while leaving it in the model frame, so
the formula, `predict()` and `newdata` all stay as they are, and for the
varying-coefficient work it is how a covariate is given to one forest and
withheld from another.

Now allowed. Negative and non-finite still refused. Every predictor at zero is
refused, since nothing would be left to split on, and that check belongs in
`resolve_split_weights()` rather than in `bartisan_control()`: the control
function sees only the names the caller gave, so `c(x2 = 0)` on a three-predictor
model looks all-zero there while the two unnamed predictors default to one. The
first version of the check was in the wrong place and a test caught it.

Verified that a zero weight really does hold the predictor out: `splits` and
`prop_used` both come back exactly 0 for it, the term is still in
`terms(fit)`, and `predict()` still works.

## The atom at zero in the ATE, and what actually removes it

`_dev/ate-atom.R`. The question was whether it is expected that a credible bound
comes out exactly zero whenever the treatment's `prop_used` is below the
confidence level, and whether the treatment should be given a larger splitting
weight.

**It is exact arithmetic, not a coincidence.** Under the DART sparsity prior a
draw in which the treatment is in no tree makes the contrast exactly zero, so the
posterior of the ATE is a mixture: mass `1 - prop_used` at the point zero, the
rest spread over nonzero values. If the nonzero part is one-signed, the 2.5%
quantile is zero as soon as `1 - prop_used > 0.025`, that is as soon as
`prop_used < 0.975`. Measured on the RHC fit: `prop_used` 0.898 and the share of
draws whose average contrast is exactly zero is 0.102, which is `1 - 0.898` to
three places.

RHC, four chains:

| | `prop_used` | atom mass | ATE |
|---|---|---|---|
| `sparsity = TRUE` (default) | 0.898 | 0.102 | 0.0531 [0.0000, 0.1058] |
| `sparsity = FALSE` | 0.999 | 0.001 | 0.0623 [0.0137, 0.1114] |
| `split_prior = c(rhc = 1)` | 0.999 | 0.001 | 0.0623 [0.0137, 0.1114] |
| `split_prior = c(rhc = 5)` | 1.000 | 0.000 | 0.0668 [0.0172, 0.1154] |
| `split_prior = c(rhc = 20)` | 1.000 | 0.000 | 0.0670 [0.0183, 0.1187] |

**A larger weight on the treatment is not the answer, and is barely even an
answer.** What creates the atom is the prior being *able to drop* the treatment,
not the treatment being under-weighted. Fixing the weights at all removes it:
`split_prior = c(rhc = 1)` is uniform weights and reproduces `sparsity = FALSE`
exactly. Going from a weight of 1 to 5 to 20 moves the estimate 0.0623 to 0.0668
to 0.0670, which is inside the noise and flat after the first step.

**And the vignette's headline number is affected.** With the atom gone the
interval no longer reaches zero: 0.0623 [0.0137, 0.1114] rather than
0.0531 [0, 0.1058]. So `vignette("causal")` currently reports an interval whose
lower bound is a property of the prior rather than of the data, under prose
saying the direction is reasonably clear -- which understates it. Not changed
here, because the vignette is being edited.

The general statement, which `?bartisan_control` already makes for contrasts and
which this makes concrete for the ATE: a variable-selection prior on the variable
whose contrast is the estimand is answering a question the analysis did not ask.
For a causal estimand, either `sparsity = FALSE` or any `split_prior` is right,
and `split_prior` is the better of the two when the other predictors are many
enough that treating them all alike is wasteful.

## Should `sparsity = TRUE` stay the default? Yes, and the guidance had to get sharper

`_dev/sparsity-default.R` and `_dev/sparsity-effect-size.R`.

**The first pass found no atom at all**, at any sparsity level, with a treatment
effect of 0.5 against residual sd 1: `prop_used` 1.000 for the treatment at
p = 10 and at p = 50, estimates and intervals identical across the four
settings. So the atom is not a property of the prior on its own. It is what the
prior does to a predictor weak enough to be droppable, which is why RHC shows it
(`prop_used` 0.898, an ATE of 0.06 on a probability scale) and this design does
not. Effect size was the variable to sweep, and sweeping it turned the default
question into a measurable one.

**Prediction: any sparsity beats none, and the strength hardly matters.** RMSE
against the true regression function on held-out data, Friedman, five relevant
predictors, mean of three replicates:

| predictors | none | weak | moderate | strong |
|---|---|---|---|---|
| 10 | 0.446 | 0.385 | 0.400 | 0.374 |
| 50 | 0.465 | 0.346 | 0.372 | 0.362 |

Between 13% and 26% better than `"none"`, and `"weak"` is within noise of
`"strong"`. Worth noting that the gain is no larger at p = 50 than at p = 10,
which is not what the high-dimensional framing of DART would suggest.

**A weak contrast: sparsity is actively harmful, not just cosmetically.** Binary
treatment among 20 predictors, continuous outcome, residual sd 1, n = 800, five
replicates:

| true effect | setting | estimate | bias | atom | covers |
|---|---|---|---|---|---|
| 0.05 | none | 0.031 | -0.019 | 0.08 | 0.80 |
| 0.05 | moderate | 0.000 | -0.050 | 0.89 | 0.40 |
| 0.05 | strong | 0.000 | -0.050 | 0.92 | 0.40 |
| 0.10 | none | 0.131 | +0.031 | 0.03 | 1.00 |
| 0.10 | moderate | 0.029 | -0.071 | 0.69 | 1.00 |
| 0.20 | none | 0.161 | -0.039 | 0.05 | 1.00 |
| 0.20 | moderate | 0.094 | -0.106 | 0.55 | 0.60 |
| 0.20 | strong | 0.085 | -0.115 | 0.51 | 0.60 |
| 0.50 | none | 0.475 | -0.025 | 0.00 | 1.00 |
| 0.50 | moderate | 0.474 | -0.026 | 0.00 | 1.00 |

The estimate is attenuated by half or more and the 95% interval covers at 0.40
to 0.60. That is a wrong answer presented confidently, not a reporting quirk, and
it is the thing the atom was a symptom of. At a true effect of 0.5 every setting
agrees to three places, because the prior never has reason to drop a predictor
that is earning its splits.

**So a moderate default does not split the difference, and that was the useful
finding.** `"moderate"` is as bad as `"strong"` on the contrast side (coverage
0.40 and 0.60, bias -0.05 to -0.12) and `"weak"` is as good as `"strong"` on the
prediction side. The argument behaves close to a switch rather than a dial, so
there is no middle setting to retreat to.

**Default kept at `TRUE`.** Whoever does not set it is more likely predicting
than estimating a treatment effect, the literature expects DART on by default,
and the cost of the wrong choice is asymmetric in the other direction too: with
sparsity off you lose 13% to 26% of predictive accuracy, which is a worse answer
but not a misleading one, whereas with it on you can get a halved effect with a
60% interval. What changed is the documentation, which now carries both tables
and a recommendation indexed by estimand rather than the previous general
warning.

**The principled fix is the varying-coefficient work, not a new sparsity
setting.** What is wanted is a sparsity prior over the covariates that cannot
touch the treatment, and putting the treatment in its own forest is exactly that:
it is why Bayesian causal forests separate the prognostic and treatment surfaces
in the first place. `split_prior` is the interim answer and turns the prior off
altogether, which is coarser.

## `chains` moved, `num_save` became `num_draws`, `sigma_mu` left the rhat table

`chains` is a `bartisan_control()` argument. Backward compatible without doing
anything: `merge_control()` takes its allowed set from
`names(formals(bartisan_control))`, so every existing `bartisan(..., chains = 4)`
call still reaches it through `...`, which is how the vignettes and 28 test call
sites keep working unchanged.

`num_save` renamed to `num_draws` at 373 sites across R, C++, tests, vignettes
and `_dev`.

`num_burn`, `num_draws` and `num_thin` each gained a sentence saying what raising
it does: warmup buys convergence, draws buy precision, and thinning buys neither
and is only for holding down memory.

`sigma_mu` is out of `scalar_draws()` and so out of `fit$rhat`, with the
cross-package numbers as the comment explaining why. `fit$sigma_mu` is unchanged.
`vignette("diagnostics")` lost its `sigma_mu.*` row description and gained a
paragraph saying the quantity is deliberately absent and what the other packages
do. Two tests.

**One test had to change and one had to move.** `test-gate.R` deliberately fits a
single tree to a clean step, which is near-separable, and now that the ramp fix
lets the leaf scale actually equilibrate it climbs high enough to trip
`warn_runaway_scale()`. The test is about gate shape, so it pins the scale with
`update_sigma_mu = FALSE`. And `test-marginaleffects.R` asserted that a survival
credible interval covers the truth at seed 43, where the truth fell 0.004 inside
the accelerated failure time interval against 0.03 to 0.10 at other seeds; any
change to the sampler flipped it, and the ramp fix did. Ten seeds all cover, so
the assertion is sound and the seed was the problem. Moved to seed 3.

## More than one chain no longer needs future.apply

It refused to run without it. Parallelism is how fast the chains are, not
whether the model is fitted, and several chains run one after another is still
what makes the convergence diagnostics available, so the absence of an optional
package should not be an error.

**The seeds are generated here now rather than by `future.seed = TRUE`.** That
was the part worth getting right: the two branches would otherwise draw from
different streams, and the same script would give different answers depending on
whether future.apply happened to be installed, which is a worse failure than
being slow. `parallel_streams()` advances one L'Ecuyer stream per chain from the
session's own seed and both branches consume the same list, so the draws are
identical either way -- verified, and a stronger guarantee than the package made
before. `parallel` moved to Imports for `nextRNGStream()`.

The session's generator is put back afterwards, kind included, since the streams
need L'Ecuyer and the session did not ask for it. Two tests, one of them mocking
`rlang::is_installed()` to take the sequential branch with future.apply present.
The old test asserting the refusal is gone.

## Per-forest arguments, including the formula

The interface change behind varying coefficients, and useful on its own. Any
argument that could mean something different for each forest of a multi-forest
family may now be given once, to apply to all, or one per forest: positionally,
or keyed by the forest names. That covers `formula`, `num_trees`, `k`,
`sigma_mu`, `sparsity`, `split_prior`, `bandwidth`, `gamma`, `beta`, the four
`alpha` arguments, and the three `update_` flags.

```r
bartisan(list(y ~ x1 + x2, ~ x2 + x3), data = d, family = location_scale(),
         num_trees = c(mean = 50, log_sd = 10))
```

**Per-forest predictors cost no new machinery in the engine, because zero
splitting weights already do the job.** The frame is built from the union of
every forest's formula, and each forest is then held to its own subset by zeroing
its splitting weights on the terms its formula leaves out. So a predictor absent
from one forest's formula is present in that forest's data and never split on,
which is the right semantics, and the tree and node code is untouched. This is
what the zero weight in `split_prior` was for. Measured: on
`list(y ~ x1 + x2, ~ x2 + x3)` the mean forest takes 42.3 and 33.4 rules on `x1`
and `x2` and exactly 0 on `x3`, and the scale forest exactly 0 on `x1`.

`split_prior` became a matrix with one column per forest, and the engine reads
column `h`. The other per-forest settings are read through two small lambdas
that index by forest, so R sends one value per forest and the engine never has
to decide what a scalar means. `NULL` is returned when there is nothing to say --
no weights asked for and every forest using every predictor -- so the ordinary
call reaches the engine exactly as before, and a test checks the draws are
identical between `y ~ x` and `list(y ~ x)`.

**A forest a named argument does not mention keeps that argument's default**
rather than borrowing the value chosen for another forest. Naming one forest is
the natural way to say "leave the other alone", and the alternative would make
`k = c(log_sd = 8)` silently set the mean forest's `k` to 8 as well.

**The multinomial families are the exception**, as asked. Their forests are the
levels of one categorical parameter and act together, so every per-forest
argument takes a single value and more than one is an error naming the reason.

**`sparsity` had to be vectorized** rather than just spread, because it is one
argument standing in for four: `resolve_sparsity()` now resolves each element on
its own and carries the names through, and `bartisan()` spreads the four derived
settings. That makes `sparsity = c(mean = TRUE, log_sd = FALSE)` work, which is
the setting the atom-at-zero finding wants and the shape a BCF fit needs.

`bartisan-families` gained a "Several additive predictors" section with the forest
order and names per family, linked from `bartisan()`'s `formula` argument.
`bartisan_control()` gained an "Arguments that vary by forest" section. The names
are the ones already in `fit$eta` and the diagnostics table -- `mean` and
`log_sd` for `location_scale()`, `count` and `zero` for the zero-inflated
families -- rather than `mu` and `sigma`, so that one set of names labels
everything.

### Six bugs the tests found, none of which I would have found by reading

Written down because every one of them was in code that looked right.

- **A named list may carry the response on any element.**
  `list(log_sd = ~ x3, mean = y ~ x1)` is legitimate, and both
  `split_formula_list()` and `bartisan()` required the response on element 1.
  Names say which forest each formula is for, so position carries nothing; an
  unnamed list is still positional and must lead with the response.
- **`union_formula()` read the response off element 1** for the same reason, and
  so built `x3 ~ ...` from the list above. It now finds whichever element carries
  one, since the union is built before the family is known and therefore before
  the list can be put in forest order.
- **`as.integer(num_trees)` dropped the forest names**, so
  `num_trees = c(mu = 8)` reached the resolver unnamed and the bad name went
  unreported.
- **`~ .` on a later formula expanded over the response**, making the outcome a
  predictor of itself. The first formula's left-hand side is put back before the
  terms are taken.
- **`x1:x2` and `x2:x1` became two predictors.** Term labels are compared by
  their set of variables now, so an interaction written either way is one term.
- **`forest_masks()` returned a vector, not a matrix**, with a single predictor
  group, and everything downstream indexes by column.

Two more, smaller: a list of one formula left `formula` as a list and broke
`terms(fit)`; and `split_prior`'s own names are predictor names, so a bare named
vector cannot also be read as keyed by forest -- `c(x3 = 0)` is a weight on `x3`
for every forest, and a *list* is what says per-forest.

`arg::when_not_null()` does not accept a plain closure as its checker, which
three recovery tests found; the optional per-forest checks use plain R instead.
And one existing test matched on the old `num_trees` error message, which now
comes from the per-forest resolver and names the forests the family does have.

## Three papers assessed: VCBART, flexBART, SBT

All three are the VCBART citation neighborhood, which is worth noticing: the
VCBART discussion cites Luo and Pratola's sharding paper as its scaling route and
Deshpande's flexBART for its categorical decision rules.

### VCBART (Deshpande, Bai, Balocchi, Starling, Weiss; Bayesian Analysis 2026)

The varying coefficient model with a BART ensemble per coefficient, over the
*modifiers* rather than the covariates. BCF is named in the paper as a special
case with one covariate.

**Build it.** It is the paper the package is closest to and the gap is one
decorator: the map from forests to the likelihood is
$\mu_i = \eta_{i0} + \sum_j x_{ij}\eta_{ij}$, and by the chain rule the score
and information for forest $j$ scale by $x_{ij}$ and $x_{ij}^2$. That is
`_dev/bcf-interfaces.md`'s finding arrived at independently, and this term's
per-forest work supplies the rest: a formula per forest is exactly the
covariate-versus-modifier split, and per-forest `sparsity` is their per-ensemble
modifier selection.

**Two of their stated future work items are things this package already has.**
Their hard trees cap the recoverable smoothness at Holder $\alpha_j \le 1$ and
they name soft rules as the fix and as ongoing work; bartisan has had soft rules
from the start. They also describe a heteroskedastic VCBART needing a variance
ensemble in the manner of Pratola (2019); that is `location_scale()`.

**One thing they have that we do not**, and it is not the varying coefficient: a
compound-symmetry within-subject correlation with $\rho$ drawn, which makes the
leaf update an intercept-free conjugate linear regression per leaf rather than a
scalar draw. bartisan models repeated measures with a random intercept instead,
which is a different model rather than a worse one, but it is not the same thing
and should not be described as if it were.

**Their diagnostics advice corroborates the `sigma_mu` decision made this term.**
They tell users to track $\sigma$ rather than trees, note that individual trees
are not identified, and report needing 20,000 to 50,000 iterations for
$\hat{R} < 1.1$ while 2,000 gives good predictions and calibrated intervals. That
is the same shape as the finding here: the ensemble mixes slowly on quantities
nobody reports and fast enough on the ones they do.

### flexBART (Deshpande; arXiv 2211.04459)

One-hot encoding a $K$-level factor lets a tree express only $2^K - K$ of the
$B_K$ partitions of its levels, because a rule on a single indicator can only
peel off one level at a time. At $K = 5$ that is 27 of 52; at $K = 10$ it is
1,014 of 115,975. The paper re-implements BART with rules that assign subsets of
levels to each branch, plus a decision-rule prior that produces spatially
contiguous clusters by deleting an edge from a random spanning tree.

**This is a real gap here, and it was checked rather than assumed.**
`build_design()` calls `contrasts(..., contrasts = FALSE)`, so a factor becomes
$K$ full dummy columns, and the engine's rules are thresholds on single columns.
bartisan is exactly in the described regime. What `make_group_probs()` already
does is the *other* half: it makes a factor one unit for the sparsity prior, so
selection treats it as a whole. Partitioning of its levels is untouched by that,
and the two should not be conflated.

The cost is partial pooling, not fit. On five levels with a strong signal
bartisan recovers the level means fine, because four splits isolate five levels.
It bites where levels are many and thin, which is what their baseball and census
tract examples are.

**Worth doing, and it is the deepest of the three**, because a decision rule
would have to become a subset rather than a threshold, which reaches the node
representation, the prediction path and the missing-value handling. One design
question has no obvious answer: a soft gate is a smooth function of a threshold,
and there is no evident soft analogue of "this subset of levels goes left", so
soft rules and categorical subsets would need a decision about how they compose.

### SBT (Luo and Pratola; arXiv 2306.00361)

A sharding tree on an auxiliary uniform variable partitions the data into $B$
shards, a separate BART is fitted to each, and predictions are a weighted
aggregate. The theory gives posterior concentration for the aggregate and shows
the weights should equalize $w_b^{-1}\varepsilon_{b,n}$, so equal shards want
equal weights. Prediction draws a fresh $u_*$ per iteration.

**Do not build it.** It is a scalability device rather than a modeling extension,
and it changes what the model is: the fit becomes a mixture over shardings, so
every family, every estimand and `marginaleffects` would have to account for the
sharding, and prediction stops being a function of the covariates alone. The
payoff is parallelism across shards, and this package's cost is not
sharding-shaped -- the measured wins have been leaf-target form and augmentation,
which cut five to ten times off the general families and leave the interface
alone. Recorded as the route to look at if sample size ever becomes the binding
constraint, which it is not.

## Subset splitting rules for a factor, and what they are actually worth

`_dev/categorical-priors.R`, `_dev/categorical-check.R`,
`_dev/categorical-benchmark.R`. flexBART's contribution, written from the paper
rather than from their code.

### What the package was doing, established rather than assumed

A factor became `K` full indicator columns, the sparsity prior picked the factor
as one *group*, then one indicator *column* uniformly, and the rule was a
threshold on that column -- so it peeled a single level off the rest. What
`make_group_probs()` already did was the other half of the problem, selection,
and that half was right; the partition of the levels was untouched by it.

Simulating the tree prior directly, at K = 10:

| | reachable partitions | mean co-clustering | spread | singleton levels per tree |
|---|---|---|---|---|
| one-hot | 1,014 of 115,975 | 0.770 | 0.003 | 1.18 of 2.18 leaves |
| subsets | all | 0.459 | 0.004 | 0.33 of 2.48 leaves |

A typical one-hot tree put one level alone and the other nine together.

**The shuffle-and-treat-as-ordinal idea was tested and is worse than the status
quo.** A threshold on one fixed order cuts a contiguous block of it, so it
reaches the `2^(K-1)` interval partitions, and `2^(K-1) < 2^K - K` for every
K >= 2: 16 against 27 at K = 5. It is also *rigid*, which is the more serious
objection: co-clustering spread 0.230 against 0.003, so adjacent levels are
pooled almost always and distant ones almost never, under an order chosen at
random. Reshuffling per tree removes the rigidity (spread 0.001) but each tree
still reaches only interval partitions. Recorded and not pursued.

**Two incidental findings, both fixed by the same change.** Between 1% and 7% of
factor rule draws landed on an indicator the path had already used up, giving an
empty child, because `get_limits()` tracked an interval per column and not
whether the column was exhausted. And 81% of cutpoints on a 0/1 column left one
level with a *fractional* membership weight under the default soft gate, which is
meaningless for a category. A categorical rule is now always hard, in a soft tree
too: a gate is a smooth function of a distance and there is no distance between
two levels.

### The implementation

A rule on a categorical group holds a bitmask of the levels that go left.
`Node::mask` is empty for a numeric rule, so the common case allocates nothing
and a node from the pool keeps its capacity. The available levels come from
intersecting over every ancestor that split on the same group, masks being
absolute level sets. The prior is the uniform distribution over the `2^m - 2`
non-degenerate subsets of the `m` available levels, drawn by assigning each level
to the left with probability one half and rejecting the two degenerate draws --
which is drawing from the prior, so the rule still cancels out of the acceptance
ratio the way the variable and the cutpoint do.

The engine gets an integer matrix of level codes alongside X, and a rule's `var`
is then a column of that rather than of X. Testing a level is a shift and an
`and`, which is where the efficiency is: flexBART keeps a `std::set<int>` per
rule. Which groups are categorical is decided by the columns rather than the
terms -- indicators, exactly one set per row -- which admits a factor and an
interaction of factors and correctly excludes a factor crossed with a numeric
predictor.

`categorical = "onehot"` in `bartisan_control()` keeps the old rule, expressed by
telling the engine that no group has levels. Useful for comparison and it is what
the benchmark below uses.

**Verified against the thing it is for.** On a one-tree forest on pure noise at
K = 5, the sampler visits 51 of the 52 partitions, 25 of which one-hot cannot
form at all, and 41% of draws are in one of those 25. Under
`categorical = "onehot"` the count of such draws is exactly zero rather than
merely small.

### Three bugs, two of them mine and one pre-existing

- **`Node::Rule` did not carry the mask.** The change move restores the old rule
  when its proposal is rejected, and restoring `var` and `group` from a
  categorical rule while leaving the mask cleared leaves a node whose rule says
  numeric and whose `var` indexes the codes matrix, or the reverse. That read
  out of bounds on the first fit mixing a factor with a numeric predictor. Found
  by instrumenting `gate()` rather than by reading, after three wrong guesses.
- **`get_limits()` compared `y->var == var` across rule kinds.** A categorical
  rule's `var` is a column of level codes, so it could collide with a numeric
  column index and constrain a cutpoint for no reason.
- **`predict_tree()` and `predict_accumulate()` were dead** -- nothing in the
  package or the tests called them -- and both read a rule as a threshold on a
  column of X. Deleted rather than fixed: a dead path that silently mishandles a
  categorical rule is a trap for whoever calls it next.

### What it is worth, measured

20 levels in 4 clusters of 5 sharing a mean, 50 trees, one chain, five
replicates, RMSE against the true mean function. Hard rules for the first two
rows, which is what flexBART has, so the comparison is of the categorical rule
and not of soft against hard:

| n (per level) | subset, hard | onehot, hard | subset, soft | onehot, soft | flexBART |
|---|---|---|---|---|---|
| 200 (10) | **0.3364** | 0.3635 | 0.3132 | 0.3223 | 0.3398 |
| 500 (25) | **0.2732** | 0.2878 | 0.2467 | 0.2434 | 0.2814 |
| 2000 (100) | 0.1385 | 0.1370 | 0.1059 | 0.0969 | 0.1475 |

**Superseded: this table was five replicates reported as means, with no standard
error, and it did not survive twelve replicates and a paired analysis.** See the
next entry. The corrected reading is that under hard rules subset is better at
every size, and under soft rules the two are indistinguishable.

`subset, hard` matches or beats flexBART at every sample size, which is the check
that the implementation is right rather than merely different. Timing: bartisan
is three times faster at n = 200 and 1.4 times slower at n = 2000, where
flexBART's observation-to-leaf bookkeeping pays off. Subset rules cost bartisan
20% to 27% more time than one-hot above the smallest n.

## Correcting the subset-rule benchmark: five replicates were not enough

The table in the entry above was five replicates, reported as means with no
standard error, and two of the things I concluded from it were not there. Every
method sees the same data within a replicate, so the comparison is paired and the
standard error of the paired difference is what a gap has to beat. At twelve
replicates, with that standard error:

| per level | subset, hard | onehot, hard | subset, soft | onehot, soft | flexBART |
|---|---|---|---|---|---|
| 10 | 0.3705 | 0.3998 | **0.3435** | 0.3445 | 0.3715 |
| 25 | 0.2617 | 0.2725 | 0.2265 | **0.2201** | 0.2709 |
| 100 | 0.1421 | 0.1488 | 0.1098 | **0.1059** | 0.1499 |

Paired standard errors against the best row run 0.0025 to 0.0112.

**What changed.** I had said subset and one-hot were "level at a hundred
observations per level" under hard rules; they are not, subset is better there
too (0.1421 against 0.1488), by about 4% against a standard error around half
that. And I had said the ordering "reverses" under soft rules "by one or two
percent"; it does not reverse in any sense the data support -- under soft rules
the largest gap between the two is 0.006 against a standard error of 0.005, so
**they are indistinguishable**, and the apparent reversal was replicate noise.

**What held.** Subset beats one-hot under hard rules at ten observations per
level, 7.3% against a standard error of about 2%, which is the one clear
categorical result. And `subset` with hard rules matches or beats flexBART at
every size (0.3705 against 0.3715, 0.2617 against 0.2709, 0.1421 against
0.1499), which was the check that the implementation is right rather than merely
different.

**The largest number in the table is not about categorical rules at all.** Soft
rules beat hard rules by 7% to 26% at every size, which is far more than either
categorical choice is worth. Worth remembering when reading any of the rest.

The default stays `"subset"`: it is right about the prior, it wins under hard
rules, and it never loses beyond noise. `?bartisan_control` now carries this
table and says plainly that it will not visibly improve a fit rather than
implying it will.

**The methodological lesson, since this is the second time this session.** Five
replicates and a difference of a few percent is not a measurement. The marginal
RMSE varies far more across seeds than the paired difference does -- the same
design gave 0.31 and 0.39 for one configuration under two seed sets -- so
reporting means without pairing hides the only comparison that is stable. Paired
differences with standard errors from here on.

## A factor as a predictor or as a random intercept

`_dev/factor-random-vs-fixed.R`. Now that a rule can pool levels, the fixed
route can express what a random intercept does, so the comparison is worth
making. Twenty levels, ten replicates, paired, RMSE against the true mean
function. Two truths: level effects as independent normal draws, which is the
random intercept's own prior, and level effects taking four distinct values,
which is a partition and so the subset rule's.

| truth | per level | fixed, subset | fixed, onehot | random intercept | both |
|---|---|---|---|---|---|
| iid | 10 | 0.3880 | 0.3830 | **0.3784** | 0.3941 |
| iid | 50 | 0.1547 | **0.1500** | 0.1548 | 0.1551 |
| clustered | 10 | 0.3883 | 0.3824 | **0.3753** | 0.3773 |
| clustered | 50 | 0.1537 | **0.1509** | 0.1534 | 0.1576 |

Paired standard errors run 0.0014 to 0.0084.

**Everything is within noise of everything else, and that is the finding.** The
largest gap that clears two standard errors is `both` at ten observations per
level on the iid truth, which is *worse* than any single route -- putting the
factor in as a predictor and as a random intercept at once costs something and
buys nothing. Otherwise no route beats another by more than about two standard
errors, in either direction, on either truth.

**Neither prior wins on the truth that matches it**, which is the part I did not
expect. The clustered truth is exactly a partition of the levels and the subset
rule is exactly a prior over partitions, and it comes last there. The iid truth
is exactly the random intercept's prior and the random intercept wins by 0.005
against a standard error of 0.005. With four distinct effects of -3, -1, 1 and 3
and ten or more observations per level, every route estimates each level well
enough on its own that the prior over how levels group has almost nothing left to
do, and a fifty-tree ensemble builds the structure additively whatever any single
tree can express.

So the practical answer is that the choice is not worth agonizing over at this
number of levels, and the reason to prefer a random intercept is what it always
was: it is a statement that the levels are exchangeable draws and that new levels
are expected, which a predictor cannot represent at all. Where the two should
start to separate is many more levels with very few observations each, which this
design does not reach.

**A bug in the first version of this script, worth recording because it announced
itself.** The level effects were drawn inside the simulation from a seed that
included which dataset was being generated, so the training and test sets got
*different* level effects and every method scored an RMSE larger than the
standard deviation of the truth. A fit that cannot beat predicting the mean is
the signature of a target that is not there.

## flexBART's bookkeeping: already here, and what they have that we do not

Their paper attributes part of their speed to caching which observations reach
which leaf and updating it incrementally, rather than looping over the whole
dataset on every tree update.

**bartisan already does that**, and this was checked two ways rather than
asserted. In the code, `Node::idx` is that cache, `split_support()` divides a
parent's cache between its two children and touches nothing else, `save_support`
and `restore_support` snapshot it for rollback instead of recomputing, and the
node pool keeps the vectors' capacity across births. The structural moves are
local: a birth picks a leaf, a death picks a branch whose children are both
leaves, and `change_rule` picks a node from `not_grand_branches`, so no move
touches more than one node's support.

And in the measurement, `_dev/bookkeeping-cost.R`, Friedman with p = 10, hard
rules, 50 trees, best of three:

| n | bartisan | dbarts | flexBART |
|---|---|---|---|
| 500 | 0.24 s (2.78x) | 0.09 s | 0.32 s (3.72x) |
| 2000 | 0.82 s (3.19x) | 0.26 s | 0.62 s (2.40x) |
| 8000 | 2.65 s (2.45x) | 1.08 s | 1.63 s (1.51x) |

**bartisan's ratio to dbarts is flat in n** -- 2.78, 3.19, 2.45, no trend -- which
is the decisive fact. A missing observation-to-leaf cache costs O(n) per tree
update where the cache costs O(support), so its absence would show as a ratio
that grows with n. It does not. Turning off the sparsity draw and fixing the leaf
scale changed nothing either, so the constant factor is the generalized
machinery, which is what it is for.

**flexBART's ratio falls with n** -- 3.72, 2.40, 1.51 -- so its bookkeeping buys
scaling over *dbarts as well*, and it is faster than bartisan above n of about
1000 and slower below it. So there is something there, and it is not the thing
bartisan is missing; it is something neither bartisan nor dbarts does.

**What it would be here, and why it is awkward.** bartisan stores a support
vector on every node, so total storage is O(n x depth) and a structural move
allocates and copies the affected node's support. A single tree-level
observation-to-leaf map is O(n) and a move rewrites only the entries that
actually change. That is a real gain and it grows with n.

The obstacle is soft rules. With a soft gate an observation reaches many leaves
with different weights, so "which leaf does observation i reach" is not a
function and there is no single map to keep -- which is exactly why `Node::wt`
exists alongside `Node::idx`. So flexBART's representation is available only for
hard rules, and adopting it means a second support representation maintained in
parallel with the first, correct under rollback, under the bandwidth move, and
under the categorical rules just added. Given that soft rules are the default and
are worth 7% to 26% of RMSE against hard ones, the configuration this would speed
up is not the one most fits use.

Recorded as available and not taken. The place a soft-rule fit actually spends its
extra time is the bandwidth move's `rebuild_support()`, which is a full O(n x
depth) rebuild every `bandwidth_every` sweeps and has no flexBART analogue,
because flexBART has no bandwidth. That is the better target if soft-rule speed
becomes the goal.

## Log: `bcf()` showed two progress bars

Reported as "progressr makes multiple bars for VC models (noticed in `bcf()`)".
Counting distinct progressor UUIDs under `handler_debug(uuid = TRUE)` located it
somewhere else:

| call | progressors |
|---|---|
| plain fit | 1 |
| fit with `vc()` terms | 1 |
| `bcf(propensity = TRUE)` | **2** |
| `bcf(propensity = FALSE)` | 1 |

So varying coefficients are not involved: `bcf()` fits a propensity model and
then an outcome model, and each `bartisan()` call built a progressor of its own.
The caller asked for one fit and watched two bars, the first of which filled
while the call was nowhere near half done.

Suppressing the propensity bar was rejected: it reintroduces the problem just
fixed for the convergence pass, a bar that covers part of the wall clock and
then sits at 100% while the call keeps working. The fix shares one bar instead.
`the$claimed_progress` holds a reporter a wrapper has claimed;
`progress_reporter()` returns it rather than building a second;
`shared_reporter()` sizes it across the fits that will actually run, and `bcf()`
releases the claim in its own `on.exit()`.

Two details that took measuring rather than reasoning:

- The shared tick count has to be the **smallest** across the fits, not each
  fit's own. One reporter goes to both, and a fit reports at most once per
  sweep, so promising a 10-sweep fit 50 reports leaves the bar 40 short. With
  the minimum every chain of every fit spends exactly `ticks`, and the bar
  fills exactly in all six cases measured, including deliberately mismatched
  ones.
- Sizing must not become a validation site. `progress_spec()` merges the
  control the way `bartisan()` will, which can error on a bad argument name;
  that error belongs to the fit, with the fit's wording, so it is caught and
  the call goes without a bar rather than complaining early.

The claim also covers the adaptive retry, which fits the outcome model a second
time when a drawn coding does not apply -- verified with a Poisson `bcf()`,
which takes that path: one bar, filled exactly. Verified released after a
`bcf()` that errors, and a plain fit afterwards sizes itself normally.

## Log: `control$augment` in a fit is now an answer, not a request

`bartisan_control(augment = )` takes a flag or the names of the engine families
a rewriting may apply to, and the fit stored that request verbatim. Reading
`fit$control$augment` and finding `c("binomial", "ordinal", ...)` tells you
nothing about the fit in hand. It is now a `logical` saying whether a rewriting
was actually applied.

It has to come from C++ rather than be worked out in R. Whether the rewriting
applies depends on the data, not only on the family and link: `augmented_base()`
asks each candidate's `applies()`, and probit's wants a Bernoulli response.
Measured on the same data, three trials per observation:

| family | `control$augment` |
|---|---|
| `binomial("probit")`, Bernoulli | `TRUE` |
| `binomial("probit")`, 3 trials | `FALSE` |
| `binomial("logit")`, 3 trials | `TRUE` |

Polya-Gamma carries any number of trials and Albert-Chib does not, which is a
distinction no R-side reimplementation of the rule would get right for free. So
`.bartisan_fit()` reports `augmented` alongside the draws and `bartisan()`
overwrites the control's element with it. `combine_chains()` needs no change: it
starts from the first chain, and the value depends on the data and the family,
not on the chain.

The request is not lost, and the documentation says where it went:
`attr(control, "supplied")` still holds what the caller asked for.

## Log: what four chains on four cores actually cost

Asked directly: is a four-chain fit on four cores as fast as a one-chain fit?
Nearly. Measured on the M4 (four performance cores), `num_burn = num_draws =
400`, 50 trees, hard rules:

| | n = 2000 | n = 8000 |
|---|---|---|
| 1 chain, sequential | 0.73 s | 2.98 s |
| 4 chains, sequential | 3.01 s (4.11x) | 12.12 s (4.06x) |
| 4 chains, 4 workers | 0.89 s (**1.22x**) | 3.81 s (**1.28x**) |
| 4 chains, 8 workers | 0.94 s | 3.95 s |

So four chains cost about a quarter more wall clock than one, a speedup of 3.4x
and 3.2x over running them in sequence. Sequential four chains costs 4.1x one
chain, so the per-chain work really is the whole cost and the setup that
`bartisan()` hoists out of the chain closure is not a meaningful share of it.

**None of the remaining 25% is overhead that could be engineered away.** Two
measurements settle it. Returning the draws costs 0.02 s at n = 2000 and 0.05 s
at n = 8000 for all four chains together, against fits of 0.89 s and 3.81 s, so
serialization is nothing. And running p concurrent one-chain fits that return
only a scalar reproduces the whole gap: 1.18x per worker at p = 4 for n = 2000
and 1.37x for n = 8000. It is the four-fast-cores ceiling from the entry above.

That proxy is *worse* than the real four-chain fit at n = 8000 (1.37x against
1.28x) because each worker in the proxy rebuilds the model frame, the design
matrix and the quantile transform, which `bartisan()` does once in the calling
session. Which is a small confirmation that hoisting the setup out of `engine()`
was worth doing.

**Eight workers for four chains is slower than four**, at both sizes. There are
only four tasks, so the extra workers do nothing but get started. `workers` is
worth setting to the chain count rather than to `availableCores()`.

**`multicore` saves the startup and nothing else.** At n = 2000, four chains:

| plan | cold plan + first fit | warm fit |
|---|---|---|
| `multisession` | 1.48 s | 0.90 s |
| `multicore` | 0.98 s | 0.95 s |

Setting the plan alone is 0.21 s for `multisession` and 0.00 s for `multicore`,
the rest of the difference being each worker loading the package. Steady state
is a wash. Worth knowing for a single fit in a script, not worth recommending
generally, since forking is unavailable in RStudio and `supportsMulticore()` is
`FALSE` there.

**Where the remaining time is.** Since four chains cost 1.25x one chain, the wall
clock of a fit is one chain's sampling, and chain parallelism has nothing left to
give. Making fits faster means making a single chain faster: the flexBART support
representation recorded as available-and-not-taken below, and
`rebuild_support()`, which is a full O(n x depth) rebuild every
`bandwidth_every` sweeps and is where a soft-rule fit actually spends its extra
time.

## Log: the convergence pass, 35% to 40% faster, and the chunking that does not help

### The suggestion that does not hold up

Asked whether chunking into as many pieces as there are columns beats chunking
into `nbrOfWorkers()` pieces. It does not. At n = 8000, all agreeing bitwise
with the sequential result:

| | 4 workers | 8 workers |
|---|---|---|
| `nbrOfWorkers()` chunks (current) | 3.90 s | 4.53 s |
| one chunk per column, default scheduling | 5.01 s | 4.62 s |
| 4x workers chunks, one future each | 5.74 s | 6.19 s |

Under `future_lapply`'s default scheduling the change is close to a no-op,
because it re-chunks whatever it is given into one future per worker; the extra
cost is the per-element bookkeeping. Forcing one future per chunk, which is what
would actually buy dynamic load balancing on a machine whose cores differ, is
worse still and gets worse the more chunks there are: the dispatch dominates,
and the chunks are equal in cost anyway, so there is no imbalance to correct.
More chunks was monotonically worse at every multiple tried, up to 6.28 s at 16x
workers.

### What did help: doing less per column

Profiled per column, `rank_normalize()` was the cost, and it ran **seven** times
for every column. Three changes, in order of what they were worth:

**The tail effective sample size does not need ranking at all.** The indicator
takes two values, so rank-normalizing it is an affine map, and an effective
sample size is built from ratios of autocovariances and so invariant to one.
Over randomized cases the two agree to 8e-16 or exactly. It is also what the
quantity is defined as. Two of the seven rankings, gone.

**R-hat and the bulk effective sample size start from the same ranking.** Both
began by computing `rank_normalize(x)` on the same draws. Computed once in
`diagnosis_stats()` and handed to both, which is exact.

**The within-chain variances without `apply()`.** `mean(apply(y, 2, var))` splits
the matrix into a list and calls a closure per column, four times per column of
draws. `colSums()` on the centered matrix is the same two passes, 4.3x faster,
agreeing to the last bit.

Then one smaller one: with no ties, average ranks are just the inverse of the
ordering, so `order()` answers what `rank(ties.method = "average")` does at 1.67x
the speed. Draws of a continuous quantity have no ties and those are every
column the pass walks bar a handful of scalar rows; ties and missing values fall
back to `rank()`. Worth 3% to 5% end to end, not the 20% the microbenchmark
implied.

Together, at n = 4000, alternating old and new implementations inside one R
session, six rounds each:

| | min | median |
|---|---|---|
| before | 6.505 s | 6.605 s |
| after | 3.913 s | 4.300 s |

**39.8% on the min, 34.9% on the median, faster in six of six rounds**, with the
table agreeing to 1.6e-15. Rankings per column went from seven to four, which
the shimmed call counts confirm: 28,021 against 16,012 over 4000 columns.

### A measurement trap worth recording

**pueue's default group runs four tasks at once, and this machine has four
performance cores.** Two timings queued together compete for the same cores, and
three measurements taken that way were wrong by up to 20% before the overlap was
noticed. Worse, a `R CMD INSTALL` into the library a queued test run was reading
truncated that run. Timing jobs go one at a time, and the queue gets checked
before each.

Cross-run variance on this machine is around 20% even with the queue empty,
which is larger than several of the effects being chased here. Every claim above
that is smaller than that comes from a **paired** comparison inside one session,
alternating implementations, rather than from comparing two runs.

## Fixed: a restricted forest silently turned the sparsity prior off

Asked whether `bcf()` could take a sparsity prior on the effect forest but not
the control function. It can be *written*, and for a varying-coefficient model
in general it works. In `bcf()` with a propensity score, which is the default,
it silently does nothing.

Localized by spying on the engine call. `bcf(propensity = FALSE)` with
`sparsity = c(FALSE, TRUE)` reaches the engine as `update_s = FALSE,TRUE` and no
split matrix. With a propensity score it reaches the engine as `update_s =
FALSE,TRUE` **and a 9x2 `split_prior` matrix**, and then

- `resolve_split_matrix()` builds that matrix whenever `!all(masks)`, that is
  whenever some forest may not split on some predictor, even though the caller
  gave no `split_prior` at all. A `vc()` term is exactly such a restriction, and
  the propensity score is what makes it a strict one: without it the moderators
  are all the covariates and `all(masks)` holds, which is why the bug needs a
  propensity score to show up.
- `Hypers::Hypers()` (`src/hypers.cpp:58`) then sets `update_s = false` and
  `update_alpha = false` because a column was supplied, on the reasoning that
  "weights the caller supplied are a statement about the predictors, not a
  starting point for one". Which is right for a `split_prior` a caller gave and
  wrong for a mask the package generated.

Measured: with a propensity score, `sparsity = TRUE`, `FALSE`, `c(FALSE, TRUE)`
and `c(TRUE, FALSE)` all give **bitwise identical** fits. Without one they all
differ. No warning is emitted, unlike the explicit `split_prior` case, which
does warn.

The two mechanisms are being sent down one channel. A mask says *which*
predictors a forest may split on; `split_prior` says *hold these probabilities
fixed*. They should compose -- a Dirichlet drawn over the allowed subset, zero
elsewhere -- rather than the second cancelling the first. That needs the engine
to take the mask separately from the weights, so `fixed_s` is true only when the
caller actually fixed them.

### The fix

The two statements now travel down two channels. `resolve_split_matrix()`
returns both a `prior` and a `mask`, and it fills `prior` only when the caller
actually supplied weights; `model.cpp` reads `control$split_mask` alongside
`control$split_prior`; and `Hypers` keeps an `allowed_` index of the groups a
forest may split on. The Dirichlet is drawn over those groups alone, with
concentration `alpha / length(allowed)`, and `s` is held at zero off them, so
`sample_var()` still cannot propose a predictor the formula did not name.
`update_alpha_param()` averages `log_s` over the allowed groups, since the
masked entries are negative infinity. `fixed_s`, and with it the refusal to
draw, now means only what it says: the caller fixed the weights.

Verified. With a propensity score, `sparsity = FALSE`, `TRUE`, `c(FALSE, TRUE)`
and `c(TRUE, FALSE)` now give four different fits where all four used to be
bitwise identical, and a `vc()` model whose moderators are a strict subset does
too. The mask is still absolute: the effect forest in a `vc(z, ~ x1 + x2)` model
took 2,314 splits, none of them on `x3` to `x8`. An explicit `split_prior` still
fixes the weights, `sparsity` making no difference to it.

Nothing else moved. Built the previous commit into a second library and compared
ten configurations bitwise -- `bcf()` with and without a propensity score, plain
fits with the prior on and off, soft rules, a restricted and an unrestricted
`vc()` fit, an empty `~ 1` forest, a factor, and an explicit `split_prior` --
and every one is identical. The only behavior that changes is the one that was
broken.

## Log: sparsity in bcf(), now that it can be asked for

With the mask bug fixed, the original question is answerable. Twenty
covariates, only `x2` moderates, n = 800, five replicates, propensity score
included. `cateRMSE` and `cor` are against the true conditional effect;
`x2/noise` is that predictor's share of the effect forest's splits over the
average of the other nineteen; `atom@0` is the share of draws in which the
effect forest took no splits at all.

| tau | sparsity | ATE | cateRMSE | cor | x2/noise | atom@0 |
|---|---|---|---|---|---|---|
| 1 + 1.5 x2 | `FALSE` (default) | 1.009 | 0.417 | 0.963 | 6 | 0.00 |
| | `TRUE` | 1.011 | 0.303 | 0.987 | 3689 | 0.00 |
| | control off, effect on | 0.994 | **0.254** | 0.986 | 401 | 0.00 |
| | control on, effect off | 0.992 | 0.394 | 0.967 | 6 | 0.00 |
| 1 + 0.3 x2 | `FALSE` (default) | 1.008 | 0.223 | 0.717 | 2 | 0.00 |
| | `TRUE` | 0.998 | **0.174** | 0.866 | 22 | 0.00 |
| | control off, effect on | 1.001 | 0.217 | 0.724 | 8 | 0.00 |
| | control on, effect off | 0.999 | 0.219 | 0.757 | 2 | 0.00 |
| 1 + 0 x2 | `FALSE` (default) | 1.005 | 0.150 | | 1 | 0.00 |
| | `TRUE` | 1.017 | **0.143** | | 1 | 0.00 |
| | control off, effect on | 1.003 | 0.161 | | 1 | 0.00 |
| | control on, effect off | 1.019 | 0.158 | | 1 | 0.00 |

**The failure mode that motivates `sparsity = FALSE` does not appear.** The
atom at zero is exactly zero in every cell, the ATE is within 0.02 of the truth
everywhere, and nothing is attenuated. That is the structural argument holding
up: the prior can only drop a predictor a forest splits on, and the treatment
is not one, so there is no mass to pile at zero. Where the prior does drop
every moderator it leaves a constant effect, and a constant effect estimates
the ATE rather than zero, which is why the null row is the one it costs least
on.

**Sparsity on the effect forest buys conditional-effect accuracy.** Root mean
squared error against the true conditional effect falls from 0.417 to 0.254
under strong moderation, and the effect forest goes from spending six splits on
the real moderator per noise predictor to four hundred.

**Under weak moderation the control forest wants it too**, which the earlier
measurements in `?bartisan_control` predict: for prediction any sparsity beats
none, and the control function is a prediction problem. `TRUE` on both beats
the effect forest alone, 0.174 against 0.217.

Not enough to move the default on: one data-generating process, five
replicates, and the ATE column here comes from recentring the effect forest's
draws rather than from a proper contrast, so treat the coverage column as
indicative only. What it does settle is that the reason the documentation gave
for the default was the wrong reason, which the documentation now says.

## Log: bcf() now uses the package's sparsity default

`bcf()` used to force `sparsity = FALSE` on the outcome model. That is the right
setting for a contrast on a predictor **a forest splits on**, where the
variable-selection prior can drop the predictor whose contrast is the estimand
and leave a point mass at exactly zero. It is not the situation `bcf()` is in:
the treatment is the coefficient, carried by a forest of its own, so no
splitting proportion can drop it.

The override is gone, so `bcf()` inherits `sparsity = TRUE` like everything
else. The measurements in the entry above are what justify it: over three
moderation strengths the point mass is exactly zero at every setting, the
average effect is within 0.02 of the truth everywhere, and the conditional
effect is recovered better with the prior on (root mean squared error 0.303
against 0.417 under strong moderation, 0.174 against 0.223 under weak). Either
forest can still be set on its own, and `sparsity = FALSE` still does what it
always did when a caller asks for it.

`?bcf` said the wrong thing and now says this. `vignette("causal")` listed
`sparsity = FALSE` among the five settings `bcf()` chooses for the caller; it is
four now, and the vignette says why it is not among them. The vignette's own
`bartisan()` fits keep `sparsity = FALSE` and are unaffected, because those are
exactly the single-forest case the setting is for: the treatment is one
predictor among many there. All of its reported numbers come from evaluated
chunks rather than prose, so nothing there goes stale.

## Assessment: the two warm-start items, and what to do instead

Both To Do items aimed at the same thing, shortening burn-in: soft random tree
features as a warm start, and a grow-from-root warm start for hard rules. The
McCartan entry already noticed they are one item, since prior-drawn features
cost one ridge solve and come from the prior. **Neither is worth building.** The
measurement that settles it is not about either technique; it is about how much
burn-in there is to save.

### The transient is 34 to 70 sweeps, against a default of 500

Fitted with `num_burn = 0` so that every sweep is retained, the log likelihood
reaches within two standard deviations of its eventual level at:

| fit | sweep |
|---|---|
| Gaussian, soft | 36 |
| Gaussian, hard | 55 |
| `ordinal()`, hard | 61 |

And on the case that should be worst for warmup, 30 predictors with 25
irrelevant, the share of splits on the five that matter goes 0.20 at sweep 1,
0.59 at 25, 0.84 at 100, and reaches its plateau by sweep 34 under soft rules
and 48 under hard.

**So a warm start can remove at most 34 to 70 sweeps of a 1000-sweep run, or
about 5%.** It cannot even reach that: it would have to be *better* than what
the sampler does in those sweeps, and a basis drawn from a uniform prior over
predictors is what the McCartan replication above measures as worse on exactly
the sparse case where warmup takes longest (0.874 against DART's 0.958). The
sampler would spend its early sweeps undoing the warm start.

**The premise in the XBART entry is false as measured.** It said burn-in length
had become the binding constraint on an ordinal fit once the augmentations
landed. Burn-in length is not binding: it is already five to ten times longer
than the transient. What binds is effective sample size per sweep, which a warm
start does not touch, since it changes the starting point and not the kernel.

### What the random-feature prototype is actually worth

Sixty lines of R, n = 1500, p = 10 Friedman, out-of-sample R-squared:

| trees | features | R-squared | seconds |
|---|---|---|---|
| 5 | 13 | 0.196 | 0.02 |
| 20 | 53 | 0.722 | 0.05 |
| 50 | 132 | 0.876 | 0.25 |
| 200 | 515 | 0.900 | 5.33 |
| 500 | 1261 | 0.908 | 45.04 |

Full soft BART with DART at 50 trees, four chains, reaches **0.975 in 14.0 s**.

Two things worth recording. The residual variance at 50 features is five times
BART's (0.124 against 0.025), so this is a real accuracy sacrifice and not a
free lunch; and the ridge solve is cubic in the feature count, so "a fraction of
the cost" stops being true past a few hundred features (45 s at 500 trees is
three times the full fit). Where it stands up is the middle of the table: 0.876
in a quarter of a second, 56 times faster than the sampler. That is a case for a
`random_features()` estimator for use inside a loop, which is a separate
deliverable from anything about speed, and not a case for a warm start.

### What was done instead

`num_burn`'s default went from 500 to 200. Measured across seven designs and
four families, out-of-sample error at 200 is within a standard error or two of
500 and better on three of them, effective sample size per second improves by
1.4 to 2.0 times, and wall clock falls by about a third. `dpm()` is the one
family with a real cost, 4 standard errors and about 2.5% of its error, which
`?bartisan_control` now says.

Both items are struck. The third McCartan item, a sparsity-aware feature draw,
is untouched and is the only one of the three that addresses the measured
weakness.

## Log: why a propensity model mixes badly, and the setting that fixes it

Reported as the propensity model in `vignette("causal")` never mixing well, even
at 8 chains and 10,000 draws. It is not the chains or the draws. On `rhc` with
its 13 covariates, four chains of 1000 draws after 500 warmup:

| setting | loglik R-hat / ESS | eta R-hat / ESS | mean `sigma_mu` | seconds |
|---|---|---|---|---|
| default | 1.163 / 18 | 1.099 / 28 | 0.420 | 16.6 |
| `sparsity = FALSE` | 1.099 / 29 | 1.057 / 57 | 0.401 | 17.1 |
| `probit` link | 1.188 / 16 | 1.087 / 34 | 0.210 | 19.5 |
| `augment = FALSE` | 1.109 / 28 | 1.058 / 56 | 0.382 | 175.1 |
| 200 trees | 1.114 / 26 | 1.030 / 146 | 0.208 | 95.6 |
| hard rules | 1.077 / 36 | 1.077 / 39 | 0.289 | 10.9 |
| **`update_sigma_mu = FALSE`** | **1.006 / 249** | **1.020 / 235** | 0.212 | 25.2 |

**It is the leaf scale.** This is the failure `warn_runaway_scale()` names, one
notch below the threshold that makes it warn. Predicting who was treated is a
classification problem the covariates do well at, the likelihood rewards an
ever-larger predictor, and a half-Cauchy prior with no upper bound is not enough
to pin the scale down. It does not run away here, it merely fails to settle:
`sigma_mu` averages 0.420 against the 0.21 it takes when it is fixed or when 200
trees identify it. Everything downstream inherits that, which is why R-hat sits
above 1.1 however many draws are taken.

The two other things that help say the same thing. More trees works because each
leaf then carries less and the scale is better identified, and it costs 5.8x the
time to get a third of the benefit. Turning off the sparsity prior is worth about
2x on its own, which is the separate and already-documented fact that the
variable-selection state mixes slowly.

**Held out, on four 75/25 splits:**

| setting | eta ESS | ESS/second | log score | AUC |
|---|---|---|---|---|
| default | 91 | 5.02 | -0.5764 | 0.7390 |
| `sparsity = FALSE` | 80 | 4.40 | -0.5771 | 0.7373 |
| `update_sigma_mu = FALSE` | 301 | 16.61 | -0.5806 | 0.7318 |
| both | **443** | **27.51** | -0.5803 | 0.7324 |
| 200 trees | 185 | 4.06 | -0.5780 | 0.7363 |

**5.5x the effective sample size per second for about 1% of AUC.** For a score
that goes on to be a covariate in the outcome model that is a good trade, and
badly mixed is the worse failure. Not made a default for `bcf()`'s propensity
model on the strength of one dataset; it is written up under
`?bartisan_control`'s `update_sigma_mu` so that it is findable.

## Log: what settings a propensity model wants, and a simulation that could not say

`_dev/propensity-settings.R`, 25 replicates on four designs, eight
configurations of the propensity model scored on the **effect** rather than on
the propensity model. The earlier answer to this question scored the propensity
model's own held-out AUC, which is the wrong quantity: the score is a covariate
in the control function, nothing is reported about it, and a setting can predict
treatment slightly worse while leaving less confounding in the effect.

**The comparison came out empty, and the anchors are what said so.** The design
carries two anchors, an `oracle` handed the true propensity score and a `none`
fitted without one at all. Paired on absolute bias:

| design | oracle minus none | paired SE | t |
|---|---|---|---|
| rhc, linear | +0.0022 | 0.0130 | 0.17 |
| lalonde, linear | -0.0255 | 0.0140 | -1.82 |
| rhc, nonlinear | +0.0126 | 0.0073 | 1.74 |
| lalonde, nonlinear | +0.0189 | 0.0113 | 1.66 |

**Knowing the true propensity score is worth nothing here**, and on two designs
the sign says a fit without one did better. So no setting of the model that
estimates it could matter, and none did: every one of the six is within 1.5
standard errors of the default on every design. The diagnostic that catches this
without reading any of that is the anchor gap against the spread among the
settings under test:

| design | anchor gap | spread among the six |
|---|---|---|
| rhc, linear | 0.0129 | 0.0381 |
| lalonde, linear | 0.0247 | 0.0366 |
| rhc, nonlinear | 0.0171 | 0.0309 |
| lalonde, nonlinear | 0.0241 | 0.0105 |

Three of the four have the settings varying by more than the whole range the
score could possibly be worth, which is the signature of noise rather than
signal.

**Why the score had nothing to do.** Both surfaces were built on the same strong
covariates, so the control function, which is a forest over all of them, absorbs
selection by itself. A propensity score earns its place under
regularization-induced confounding: the selection has to sit in a direction the
outcome model shrinks away. That takes misalignment between the two surfaces,
which this design did not have and which the 2016 competition carries as an
explicit factor (`alignment`, at 0, 0.25 and 0.75), alongside `overlap.trt`.

A `targeted` confounding variant was added for this: the covariates that drive
treatment are ones the outcome barely depends on, among twenty pure-noise
columns for the sparsity prior to shrink. Running.

**Also worth recording: coverage was useless as a metric here**, 0.996 across
all 800 fits at a mean interval width of 1.68. Whatever the effect-forest
posterior for a sample average is doing, it is far too wide to separate anything
at this number of draws, and a design that leans on coverage needs to check that
first.

### The targeted variant did not separate them either

| design | oracle minus none | paired SE | t | anchor gap | spread among six |
|---|---|---|---|---|---|
| rhc, nonlinear, targeted | +0.0006 | 0.0077 | 0.07 | 0.0039 | 0.0199 |
| lalonde, nonlinear, targeted | +0.0097 | 0.0094 | 1.03 | 0.0156 | 0.0159 |

Treatment driven by covariates the outcome barely depends on, among twenty
noise columns, still leaves the true propensity score worth nothing. Every
setting is within two standard errors of the default on both designs;
`undersmooth` is the only one whose sign is consistent across them and it
reaches t = -2.01 on one and -1.53 on the other, which is not a finding.

**So the answer is that no propensity setting matters, in six designs across two
datasets, because the propensity score itself does not.** With a forest over
every covariate in the control function, selection is absorbed without it. What
would make it matter is the regime the score exists for, which is strong
selection, poor overlap, and a prognostic surface the outcome model cannot
represent. Two attempts to build that by hand failed. The 2016 competition
carries `overlap.trt` and `alignment` as designed factors and its response
surfaces are step and exponential rather than smooth, so that is the instrument
for this question and `_dev/acic2016.R` is ready for it.

**No default was changed.** The `update_sigma_mu = FALSE` recommendation from the
mixing entry above still stands on its own terms, which are mixing and speed
rather than bias; this bench says it costs nothing in bias either, which is the
most it can say.

### A bug this turned up in the harness, worth keeping in mind

Coverage came out 0.998 over the first 1200 fits. The estimand was being
averaged from the effect forest's own draws with a single shift applied, and a
varying-coefficient model splits the level between the control function and the
coefficient in a way that moves from draw to draw. That movement inflates the
interval without touching the mean, so bias and error were right and coverage
was void. `coef(draws = TRUE)` applies the recentering that `vc_recenter()`
does. Over 20 replicates: posterior standard deviation 0.468 before and 0.058
after, against a sampling standard deviation of 0.080 for the estimate itself,
and coverage 1.00 before and 0.85 after.

Both benches are fixed. The coverage and width columns of
`_dev/propensity-settings.rds` predate the fix and are void; bias, error and the
conditional-effect error are unaffected.

**0.85 against a nominal 0.95 is itself worth measuring properly.** Twenty
replicates put a standard error of 0.08 on it, so it is suggestive rather than
established, but a `bcf()` interval that is three quarters as wide as the
sampling spread of its own point estimate is the kind of thing the competition
bench reports as a matter of course.

## Fixed: a vc() model never refreshed its family's augmentation

Reported as `bcf()` giving an average effect closer to zero than plain BART on
`rhc`, which is backwards: the literature has BART attenuating a treatment
effect and `bcf()` existing to stop it.

**It was a bug, and not regularization.** `Family::before_forest()` is a no-op by
default and `VaryingCoefficientFamily` did not override it, so for every `vc()`
and `bcf()` model the inner family was never told a sweep had begun. An
augmented family redraws its augmentation there, and never did. For a logit
binomial that means the Polya-Gamma weights stayed at the 1 the constructor set,
for the whole run: a fixed pseudo-likelihood, far too precise, which pins the
predictor near the working response and compresses everything that varies.

**Measured against a known truth.** A log-odds effect of exactly 1, no
confounding, n = 2000, four replicates, so anything but recovery is the fit's
own doing:

| fit | coefficient | recovered |
|---|---|---|
| `vc()`, logit, augmented | 0.171 | **17%** |
| `vc()`, logit, `augment = FALSE` | 0.900 | 90% |
| `vc()`, probit, augmented | 0.548 | 91% |
| plain, logit, augmented | 0.878 | 88% |

**Why probit escaped.** A probit's latent draw is in `update_aux()`, which the
wrapper does forward. The families that refresh in `before_forest()` are the
logit binomial, the negative binomial, the two zero-inflated ones, the
multinomial and the multinomial probit, and all of them were affected under
`vc()`. That is why `vignette("causal")`, which fits a probit, looks sensible
and a logit `bcf()` returns a null effect.

**The tell was the log likelihood.** The broken fits reported a *positive* one,
+37.5 where every other fit on the same data reports about -830. Nothing
Bernoulli is positive. That is the same missing forward seen from the other
side: `reported_loglik()` was not forwarded either, so the base class reported
the augmented density rather than the likelihood the augmentation is a device
for.

**The fix** gives the wrapper `before_forest()` and `reported_loglik()`
overrides that combine the predictor and hand it to the inner family, which is
exactly what its `update_aux()` already did. On `rhc`, `bcf()` on a logit
binomial goes from an average effect of +0.0060 with an interval 0.043 wide to
+0.0641 with an interval 0.096 wide, against plain BART's +0.0588 and 0.112. The
contradiction the report started from is gone.

Regression test in `test-varying.R`: recovery within sight of the truth, and a
log likelihood that is negative.

**What this invalidates.** Any measurement in this file taken on a `vc()` fit of
an augmented family. The propensity-settings bench is Gaussian throughout and
the ACIC bench has Gaussian outcomes with a plain, non-`vc()` propensity model,
so neither is affected. `vignette("causal")` fits a probit and is not affected
in its estimates, though its `bcf()` log likelihood was being misreported.

## Log: the ACIC bench, and what it says about the propensity score

Six settings chosen to vary the two factors a propensity score should be
load-bearing for, four simulations each, five configurations: 120 fits on the
4802 x 58 competition covariates, scored on the sample average effect on the
treated.

| configuration | bias | RMSE | coverage | width | PEHE |
|---|---|---|---|---|---|
| `none` | +0.0096 | **0.0568** | 0.88 | 0.157 | 1.035 |
| `oracle` | +0.0098 | 0.0583 | 0.83 | 0.168 | 1.093 |
| `default` | +0.0135 | 0.0797 | 0.88 | 0.204 | 1.185 |
| `fixed_scale` | +0.0137 | 0.1067 | 0.88 | 0.176 | 1.350 |
| `both` | +0.0012 | 0.1082 | 0.79 | 0.177 | 1.363 |

Paired against the default on absolute bias, over the 24 datasets:

| configuration | difference | paired SE | t |
|---|---|---|---|
| `none` | **-0.0133** | 0.0055 | **-2.42** |
| `oracle` | **-0.0154** | 0.0058 | **-2.68** |
| `fixed_scale` | +0.0109 | 0.0122 | 0.90 |
| `both` | +0.0119 | 0.0134 | 0.88 |

**Two things, and the second is the surprise.** The settings of the propensity
model still do not matter: `fixed_scale` and `both` are within one standard
error of the default, which is now four benches saying the same thing. But
`none` and `oracle` are both **better than the default**, by about two and a
half standard errors, and they are indistinguishable from each other. Knowing
the true propensity score is worth nothing; *estimating* one costs something.

That is a statement about how the score enters rather than about how it is
fitted, and it is sharpest where the score is supposed to earn its place. Root
mean squared error by cell:

| configuration | one-term, align 0.75 | full, align 0.75 | one-term, align 0 |
|---|---|---|---|
| `none` | 0.0279 | **0.0485** | **0.0810** |
| `oracle` | **0.0263** | 0.0654 | 0.0724 |
| `default` | 0.0381 | 0.0719 | 0.1116 |

Poor overlap with no alignment, which is the hardest cell and the one the score
exists for, is where the estimated score hurts most: 0.1116 against 0.0810
without it.

Not enough to change `bcf()`'s default on. Twenty-four datasets, two chains,
400 draws, and the competition's own runs are 100 simulations per setting. What
it does say is that the open question is not which settings estimate the score
but whether putting an estimated score in the control function is helping at
all, which is a different question from the one four benches have now failed to
answer.

**Coverage is short of nominal everywhere**, 0.79 to 0.88 against 0.95, on 24
datasets per configuration so a standard error of about 0.07. Suggestive rather
than settled, and consistent with the 0.85 measured separately on the
semi-synthetic bench. This is the third measurement pointing the same way and
the instrument for settling it is simulation-based calibration, not more of
these.

## Log: simulation-based calibration, and what the interval question turned out to be

The effect intervals had come in short of nominal three times: 0.85 on the
semi-synthetic bench, 0.79 to 0.88 across the ACIC configurations. The right
instrument for that is simulation-based calibration, and the right framing is
that these are **credible** intervals, which were never promised frequentist
coverage.

**SBC passes, cleanly.** 300 replicates, a parameter drawn from the model's own
prior, data simulated from it, and the rank of the truth among 100 thinned
posterior draws:

| bin | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| count | 28 | 34 | 31 | 36 | 29 | 27 | 29 | 23 | 33 | 30 |

Thirty expected in each. Chi-square 4.2 on 9 degrees of freedom, **p = 0.90**.
Half the ranks fall in the middle half, which is what uniform means, and the
mean rank is 49.6 of 100 against 50. **The 95% interval covers at 0.957, with a
standard error of 0.012.**

So the sampler draws from the posterior it defines, and its intervals are right
for the model. Nothing needs improving.

**The two results are not in conflict, they are the definition of the
distinction.** When the truth comes from the prior, coverage is nominal, which
is what SBC just showed. When the truth is a fixed function the prior does not
favor, which is every one of the benches, the prior shrinks the estimate toward
zero -- measured at 88% to 92% recovery across the binomial cells -- and an
interval that is correctly narrow for a Bayesian misses more often than one time
in twenty. That is a credible interval behaving exactly as advertised and not a
confidence interval underperforming.

**What SBC here could and could not cover**, because the model's prior is
empirical the way every BART implementation's is:

- The leaf scale is fixed through `sigma_mu` in `bartisan_control()`, so that
  part of the prior is exact.
- The residual scale for a Gaussian is `residual_scale(y, x)` and nothing fixes
  it, so this used a `binomial()`, which has no such parameter.
- The predictor is centered at `qlogis(mean(y))`, which no prior produced. The
  calibrated quantity is therefore a **contrast** between two observations,
  where the centering is common and cancels. Calibrating the level instead would
  have been measuring the centering, which absorbs it entirely.
- The tree prior is data-free and is replicated exactly: branch with probability
  `gamma * (1 + depth)^-beta`, a variable drawn uniformly, and a cutpoint
  uniform on what the nearest same-variable ancestor left, from [0, 1] down.

`_dev/sbc.R` holds it. What it does not yet cover is a soft gate, a drawn
sparsity prior, or a family with a nuisance parameter, and each of those is a
separate replication of a separate piece of prior.

## Log: stage 1 of the mixing plan -- the predicted n-trend is not here

`_dev/sbc.R` across three sample sizes and both gates, 200 replicates a cell and
600 at the one that needed them. The question is Tan et al.'s: does the sampler's
approximation of the posterior degrade as the sample grows?

| gate | n | reps | chi-square | p | middle | mean rank | coverage |
|---|---|---|---|---|---|---|---|
| hard | 250 | 200 | 7.9 | 0.54 | 52% | 47.7 | 0.960 |
| hard | 1000 | 200 | 9.8 | 0.37 | 46% | 49.1 | 0.925 |
| hard | 4000 | 200 | 2.5 | 0.98 | 48% | 51.5 | 0.960 |
| soft | 250 | 600 | 9.4 | 0.40 | 47% | 50.0 | 0.938 |
| soft | 1000 | 200 | 4.2 | 0.90 | 48% | 50.4 | 0.935 |
| soft | 4000 | 200 | 7.5 | 0.59 | 53% | 52.9 | 0.935 |

**No cell departs from uniform**, the smallest p over six being 0.37, and
coverage runs 0.925 to 0.960 throughout. The pooled trend in mean rank against
log n is null under both gates: +1.36 per log-n under hard rules (t = 1.30) and
+0.95 under soft (t = 1.13).

**The shape is what settles it, not the p-values.** A chain that cannot reach
the high-posterior region leaves the truth in the tails and the middle empty, so
the signature is a U. There is none at any size:

| gate | n | middle | lower tenth | upper tenth |
|---|---|---|---|---|
| hard | 4000 | 48% | 10% | 12% |
| soft | 4000 | 53% | 9% | 14% |

against 50%, 10% and 10% for uniform. At n = 4000 under hard rules, which is
the case the papers actually analyze, the fit is the cleanest of the six.

The intervals meanwhile concentrate the way they should: under hard rules the
width falls 3.60, 2.76, 2.05 across the three sizes while the spread of the truth
stays flat near 2.0.

**One marginal signal, chased and dismissed.** At 200 replicates the soft trend
was +2.51 per log-n with t = 2.40, driven by the n = 250 cell sitting two
standard errors below 50. Six hundred replicates there moved its mean rank from
45.9 to 50.0 and the slope to t = 1.13. It was the small-n end, not the large-n
end where the theory predicts trouble, and it was noise.

**What this does to the plan.** Stages 2 and 3 are speculative for this package
and nothing is built. The soft-gate hypothesis in `_dev/MIXING.md` is *also* not
supported: hard rules are equally clean, so smoothness is not the explanation.
What is left is the real difference from the sampler all three papers analyze,
that a birth here draws its two leaves from a Laplace fit rather than
integrating them out, which changes the acceptance geometry and is the thing
worth understanding if this is ever revisited.

**What would change the verdict**, and neither is cheap: sizes past 4000, since
the bound is asymptotic and n = 4000 may simply be small; and a
data-generating process carrying the deep isolated signal Kim and Ročková
construct, which a forest drawn from the branching prior produces only by
accident. A prior draw is the right generator for SBC and the wrong one for
finding the worst case, and those are different experiments.

## The `diagnose()` progress bar, and the 336 MB closure behind it

The bar under a `multisession` plan started at zero, sat there for the first
third of the pass, jumped to 98 and finished. The arithmetic in
`progress_budget()` was not the problem; it had been fixed already and the ticks
add up exactly. The problem was *when* the ticks could be delivered.

**What was measured.** A `withCallingHandlers()` on `progression` timestamps each
report as it reaches the calling session. On a lalonde fit with four chains and
5000 draws, the ticks fired on the workers at 0.50s, 0.68s, 0.82s and on through
the pass, and the first twelve of them arrived at the master together at 1.31s
of a 2.8s pass. Replacing `future_lapply()` with an explicit loop showed why: the
four futures were created at 0.08s, 0.23s, 0.54s and 0.93s, and the first arrival
was 0.93s. **While the session is sending futures out it is the one thing that
cannot report**, because a worker's progress reaches a bar only when the session
looks for it, and it is busy.

**Two things were being sent that did not need to be.** Serializing the draw
matrix costs 0.065s, so the volume was not obviously the issue and an earlier
measurement had concluded that slicing the chunks "was within noise". That
measurement was wrong, and so was the first re-measurement here, in the same way:
writing `parts[[k]]` inside a future exports the whole of `parts`, so both arms
shipped the same bytes. Binding the slice to a local first separates them, and
then whole-matrix export costs 0.78s of dead bar against 0.34s.

The larger one was invisible to any reading of the code. `diagnosis_reporter()`
takes `envir = parent.frame()` so that the progressor finalizes with the caller,
and that argument stays in the frame the reporter closure was written in. The
closure was four lines long and **serialized to 336 MB**, because `envir` is
`diagnose()`'s frame and `diagnose()`'s frame holds the fit. Both ways of leaking
a frame were present at once: writing the closure where the frame is, and leaving
an argument unforced, since an unforced promise holds the environment it came
from. `reporter_for()` fixes both with `force(p)` and a frame containing nothing
else, and all three reporters now go through it, so a multi-chain fit stops
sending `bartisan()`'s frame to every chain worker as well.

| | first tick | worst gap | pass |
|---|---|---|---|
| before | 0.78s | 0.71s | 2.28s |
| after | 0.21s | 0.25s | 1.94s |

Results are identical to the sequential pass. `diagnosis_block()` moved to top
level for the same reason as the rest: a closure written inside
`diagnosis_columns()` carries `wide` with it however few columns the worker got.
The dispatch loop polls the futures already running after each one is created, so
their reports relay while the rest are still being sent.

`tests/testthat/test-invariants.R` gained a regression test that serializes the
reporter and the stepper and asserts they are under 100 KB with 16 MB sitting in
the frame beside them. `make_unit_map()` in `R/utils.R` had the same shape and
is dealt with in the entry below, where it turned out to be much the larger of
the two.

## R-hat's null is `1 + chains / ess`, and what that does to the thresholds

`diagnose()` was giving confident and wrong advice on fits that fail its
per-observation row: "the chains have each settled somewhere different" when the
true situation was that there were not enough effective draws for R-hat to say
anything at all.

**The calibration.** `_dev/ess-rhat-calibration.R` runs the package's own
`diagnosis_stats()` against a stationary AR(1) with rho = 0.99, where every chain
has the same distribution by construction and there is nothing to find, 400
replicates a cell. R-hat's null mean is `1 + m/S` for `m` chains carrying `S`
effective draws between them, and the fit is exact:

| chains | draws | true ESS | mean R-hat | `1 + m/S` | P(R-hat > 1.01) |
|---|---|---|---|---|---|
| 2 | 8000 | 80 | 1.026 | 1.025 | 86% |
| 4 | 4000 | 80 | 1.050 | 1.050 | 100% |
| 8 | 2000 | 80 | 1.103 | 1.100 | 100% |
| 16 | 1000 | 80 | 1.205 | 1.200 | 100% |

and along the other axis, four chains lengthened:

| true ESS | mean R-hat | 90th pct | P(R-hat > 1.01) |
|---|---|---|---|
| 25 | 1.162 | 1.279 | 100% |
| 100 | 1.041 | 1.068 | 100% |
| 200 | 1.020 | 1.035 | 87% |
| **400** | **1.010** | 1.017 | **44%** |
| 800 | 1.005 | 1.008 | 2% |
| 1600 | 1.002 | 1.004 | 0% |

At 400 effective draws over four chains the null mean *is* 1.010, which is where
the pairing of the two defaults comes from: `rhat_max = 1.01` and `ess_min = 400`
are one threshold stated twice, not two independent tests. Sixteen chains need
1600 for the same 1.01. The estimator also loses a little at short chain lengths,
recovering 80 effective draws as 83 over two chains and 63 over sixteen.

**What changed.** A `rhat readable` check fires when the row that drove the R-hat
failure carries too few effective draws for `rhat_max` to be reachable, compared
against `chains / (rhat_max - 1)` rather than against `ess_min`, so it tightens
correctly as chains are added. Its detail names the null value, so 1.030 over six
chains at 147 effective draws can be read against the 1.027 it would average
anyway. The advice then sends the reader to `num_draws` instead of to a
disagreement, and the warmup note drops its causal claim in the same case.

**A known limitation, not fixed.** `FAIL_SHARE = 0.2` keys the per-observation
R-hat check to the share of columns above `rhat_max`, and its comment justifies
0.2 as a noise floor on the reasoning that about 5% of columns exceed a 95%
critical value under the null. 1.01 is not a 95% critical value; at 400 effective
draws it is close to the median, and 44% of stationary columns exceed it. So the
share test is not calibrated at the boundary and will fire on a fit sitting
exactly at the recommended effective sample size. The readability check covers the
regime where this misleads; making the share test itself calibrated needs a
threshold that moves with the effective sample size, and was not attempted.

## What `diagnose()`'s worst row is actually measuring

The additive predictor now gets two rows, the worst 5% of observations and their
average, because they routinely disagree by two or three orders of magnitude and
the gap is the diagnostic.

`_dev/where-mixing-lives.R`, four chains of 2000 draws so 8000 kept:

| | n | avg R-hat | avg ESS | median unit ESS | worst 5% ESS | % over 1.01 |
|---|---|---|---|---|---|---|
| friedman, gaussian | 800 | 1.000 | 7754 | 60 | 21 | 99% |
| lalonde, gaussian | 614 | 1.000 | 8104 | 166 | 33 | 82% |
| lalonde, dpm | 614 | 1.001 | 7805 | 9 | 6 | 100% |
| lalonde, dpm, `sparsity = FALSE` | 614 | 1.000 | 8154 | 31 | 8 | 100% |
| rhc, binomial | 1500 | 1.001 | 3905 | 554 | 224 | 48% |

The average over observations carries close to one effective draw for every draw
kept, and individual observations carry tens. This is a property of forests
rather than of a hard dataset: the clean simulated Friedman fit, correctly
specified and with a Gaussian likelihood, has 99% of its observations above 1.01
and a worst-5% effective sample size of 21. Read against the null above, its
worst-5% R-hat of 1.135 is *below* the 1.190 that four chains at 21 effective
draws average anyway, so there is no disagreement there to find at all. Only the
`dpm()` row exceeds its null (1.789 against 1.667), and the family is worth
naming separately: on the same data and the same draws it costs a factor of 18 in
per-observation effective sample size against `gaussian()`, on an outcome that is
23% zeros with a long right tail.

## Chains and draws are not interchangeable, and lalonde is the worked example

`_dev/chains-vs-draws.R`, lalonde with `dpm()`, three seeds a cell.

Sixteen thousand kept draws, moved between chains:

| chains | draws | worst-5% R-hat | worst-5% ESS | average ESS |
|---|---|---|---|---|
| 2 | 8000 | 1.403 | 6 | 15882 |
| 4 | 4000 | 1.455 | 10 | 15982 |
| 8 | 2000 | 1.676 | 13 | 15942 |
| 16 | 1000 | 1.815 | 23 | 15986 |

The average's effective sample size does not notice the arrangement, which is
what it should do: effective sample size depends on the total. R-hat rises
monotonically, which is also what it should do, and is the null moving rather
than the fit changing.

Four chains, lengthened, gives an average effective sample size of 3981, 8015,
15719, 31525 and 63143 against totals of 4000, 8000, 16000, 32000 and 64000.

**The worked example.** `_dev/lalonde-reproduce.R` runs the two configurations
that were actually tried on `vignette("causal")`'s `fit_earn`, three seeds each:

| configuration | kept | worst-5% R-hat | worst-5% ESS | verdict |
|---|---|---|---|---|
| 6 chains, 20000 burn, 10000 draws | 60000 | 1.030, 1.010, 1.018 | 147, 473, 282 | warns, 3 of 3 |
| 4 chains, 30000 burn, 30000 draws | 120000 | 1.004, 1.006, 1.009 | 869, 479, 669 | passes, 3 of 3 |

Both changes pushed the same way and both were needed. Doubling the kept draws
roughly doubled the per-observation effective sample size, and dropping from six
chains to four lowered the bar R-hat had to clear from 600 to 400. The average
row reports R-hat 1.000 and an effective sample size equal to the number of draws
kept in all six fits.

**The controlled version**, which separates the two changes: twelve chains of
10000 draws after the same 30000 warmup, so the same 120000 kept draws as the
configuration that passes, arranged differently.

| configuration | kept | seconds | worst-5% ESS | worst-5% R-hat | R-hat's null | verdict |
|---|---|---|---|---|---|---|
| 4 chains x 30000 | 120000 | 134 | 869, 479, 669 | 1.004, 1.006, 1.009 | 1.006 | passes, 3 of 3 |
| 12 chains x 10000 | 120000 | 256 | 674, 490, 693 | 1.016, 1.027, 1.018 | 1.019 | warns, 3 of 3 |

The effective sample sizes are the same to within their own noise, means of 672
and 619, which is what "effective sample size depends on the total" means when it
is measured rather than asserted. Every one of the twelve-chain R-hats is at its
own null to three decimal places, so there is no disagreement there at all, and
the `rhat readable` check fires on all three. Before that check existed the
advice on these fits was to raise `num_burn` and `num_draws` together because
"the chains have each settled somewhere different", which was untrue.

The twelve-chain fits also took twice as long, since twelve chains over four
workers is three waves and each chain pays its own 30000-sweep warmup. More
chains is the cheaper axis only up to the worker count.


## The unit maps were carrying the design matrix, twice

Found by auditing for the shape of the progress-reporter bug rather than by
noticing anything wrong: every function in `R/` that returns a closure was listed
and checked for whether it forces what it captures. `make_unit_map()` was the
only one that did not, and it is called once per predictor column with its result
stored in the fit.

Both hazards were live at once, and the estimate made from reading the code was
an order of magnitude short. The closures were written in `make_unit_map()`'s
frame, so they held `x` and `ux` whether they used them or not; that is the part
reading finds, and it is the smaller part. The larger part is that **the
two-value branch never mentions `type`**, so that argument was never forced, and
an unforced promise holds the frame the call was made from. That frame is
`unit_transform()`'s, which has the design matrix in it under `x` and a second
copy under `out`. So the map for a binary predictor, whose entire content is two
numbers, retained the whole design matrix twice. Adding nothing but `force(type)`
took one from 213 KB to 4.7 KB in isolation.

Every branch now goes through a constructor that forces what it captures and is
written where nothing else is in scope: `constant_map`, `midpoint_map`,
`range_map`, `ecdf_map`. Measured on a 400-row fit with two continuous
predictors, a binary one, a three-level factor and a four-level ordered factor,
so ten columns after contrasts:

| | before | after |
|---|---|---|
| one binary or indicator map | 148.1 KB | 0.5 KB |
| one quantile map of a continuous column | 16.8 KB | 8.6 KB |
| one range map of a continuous column | 10.0 KB | 0.6 KB |
| all ten maps, quantile | 1218.4 KB | 20.8 KB |
| all ten maps, range | 1127.0 KB | 4.9 KB |

The per-map figures overstate what a fit pays, since one `saveRDS()` of the
whole list writes the shared frame once. On the RHC data, 1500 rows and 15
columns after contrasts, the list as the fit stores it goes from 788 KB to 84 KB
under the quantile transform and from 721 KB to 5.6 KB under the range one,
against a design matrix of 186 KB. The 84 KB that remains is fifteen `ecdf`
objects holding the values they were built from, which is what a quantile map is.

`predict()` is bit-identical before and after, with and without `newdata`,
under both transforms, on a fit with numeric, binary, factor and ordered-factor
predictors, as are the stored `eta` draws and the maps evaluated over a grid
running past both ends of the training range. The full suite passes.

`tests/testthat/test-invariants.R` gained the regression test, next to the
reporter one because it is the same invariant: the binary and constant maps must
not scale with the column at all, and the range map must not either. It was
checked against the old implementation, where all four assertions fail.

## Eight error messages were losing their second half

`family = glmmTMB::tweedie()` reported that the supported families are
`"custom"`, `"zip"`, `"mnp"`, `"aft"`, `"beta"` and the rest, which are the
engine's own family strings and several of them are not names a caller can
write: `"zip"` and `"zinb"` are reached through `zi_poisson()` and
`zi_negbin()`, `"beta"` through `Beta()`, `"mnp"` through
`multinomial("probit")`, `"aft"` through the three `*_aft()` functions, and
`"custom"` through `custom_family()`. The list came from `names(valid_links)`;
it is now `bartisan_family_names`, which is what the string branch of the same
function already validated against, so the two paths agree.

Looking for a place to put the hint about `custom_family()` turned up the
larger problem. `arg::err(m, .call, .envir, ...)` builds its message from `m`
alone and passes `...` to `rlang::abort()`, where a named argument becomes a
condition *field* rather than a message bullet. So `arg::err("main", i = "hint")`
compiles, runs, stores the hint on the condition, and prints only the first
half. Eight calls were written that way and every one of them was losing its
second half:

| | what was invisible |
|---|---|
| `R/families.R` `custom_family()` | that `logdens` is called as `logdens(y, eta, aux)` |
| `R/families.R` `default_family()`, twice | which families keep prior weights |
| `R/families.R` `as_bartisan_family()` | the new `custom_family()` hint |
| `R/predict.R` | `times = c(1, 5)`, the example of what to pass |
| `R/response.R`, three times | the alternative family for weights, and `ordbeta()` for a beta response at 0 or 1 |

The vector form, `arg::err(c("main", i = "hint"))`, is what `arg::msg()` is
already called with everywhere and it renders correctly. All eight are converted.
Worth watching for: nothing warns, the argument is accepted, and the message
looks complete unless it is read against the source.

The `family` argument's documentation said that ordinary `family` objects "are
used unchanged, including their links", which is what made an unsupported one
look like it should work. It now says that the object is accepted when the
distribution it names is one the package implements, since the likelihood is the
package's and a `family` object carries a link and a variance function rather
than a density, and points at `custom_family()`.

### Could a custom family be built from a `family` object automatically?

For an exponential dispersion family, yes in principle: the per-observation log
density is `-dev.resids(y, mu, wt) / (2 * phi)` plus a term that does not involve
`mu`, and the sampler needs the density only up to a term free of `eta`. Checked
against `poisson()`, `gaussian()` and `binomial()`, where the offset from
`dpois()`, `dnorm()` and `dbinom()` is constant as `mu` moves. `dev.resids()` is
per-observation, which `aic()` is not, so it is the usable half of the object.

Two things stop it being worth building on its own.

The dispersion cannot be estimated. The dropped term is `c(y, phi)` and it
carries all of the dependence on `phi`, so `-dev / (2 * phi)` alone rises
monotonically toward zero as `phi` grows: on four Gaussian observations it goes
from -3.150 at `phi = 0.5` to -0.016 at `phi = 100`, while the true log
likelihood peaks at -5.251 near `phi = 1`. An adapter would have to take a fixed
known dispersion, and it would have to refuse `loo()`, `waic()` and `pp_check()`,
which need the pointwise density rather than a shifted one.

And it would not cover the case that prompted the question.
`glmmTMB::tweedie()$dev.resids` is a stub that returns `NA` with a warning,
`aic()` returns `NA`, and the object's whole content is a link, a name and
`variance = phi * mu^power`; the density is in glmmTMB's TMB C++. What the
adapter would reach is `inverse.gaussian()`, `quasipoisson()` and hand-written
exponential dispersion families, at fixed dispersion and without the model
comparison tools. Not started.

## `tweedie()`, and why the intractable density did not matter

A family for a non-negative response with a point mass at zero and a continuous
positive part: earnings, spending, rainfall, claims. The compound Poisson-gamma
with `1 < p < 2`, one forest, `mu = exp(eta)` and `Var(y) = phi mu^p`, so the
share of zeros follows from the mean rather than having a predictor of its own.

**The objection that turned out not to hold.** The density has no closed form at
a positive response and normalizing it takes an infinite series, which is the
usual reason a Tweedie is not implemented. Writing it in exponential-dispersion
form separates the two halves:

```
log f(y; mu, phi, p) = (1/phi)[y mu^(1-p)/(1-p) - mu^(2-p)/(2-p)] + log W(y, phi, p)
```

and `log W` contains no `mu`. Verified against the compound Poisson-gamma sum
before any of it was written: `log f` minus the bracket agrees to six decimals at
`mu` of 0.8, 2.0 and 5.0. So the bracket is the whole of `logdens_unit()` and the
series is `compute_eta_free()`, the hook `src/family.h` already had for terms that
"cancel from every acceptance ratio", refreshed when a nuisance parameter moves.
The series is evaluated once per sweep, not once per leaf.

Three more things fall out of the same algebra. `P(y = 0)` is closed form,
`exp(-mu^(2-p)/(phi(2-p)))`, and is exactly what the bracket reduces to at
`y = 0`, so the zero and positive cases need no branch. The score is
`mu^(1-p)(y - mu)/phi`. The observed second derivative is
`((2-p)mu^(2-p) + (p-1) y mu^(1-p))/phi`, whose two terms are both positive for a
power in `(1, 2)` and a non-negative response, so the true curvature is usable as
it stands and there is no need to fall back on the expected information.

**Cost.** Against a Gaussian fit at n = 500 with 50 soft trees, `tweedie()` costs
5.5 times as much and `Gamma("log")` costs 5.2, so the series does not show up in
the total. Its length depends on the response and the dispersion but never on the
mean, so it does not grow as the forest moves: at the dispersion the RHC-scale
earnings data imply it runs to twenty terms, peaking at `j` between 1 and 5, and
at a dispersion eighteen times smaller the worst observation needs 143. The
summand's width scales as the square root of its peak index, so the cheap regime
is wide. `-lgamma(j+1) - lgamma(j alpha)` depends on the power alone and is
tabulated, which is what keeps the inner loop to two multiplies and a lookup.

**`power` is validated on the open interval (1, 2)**, and the whole of it is
usable. The family degenerates at both ends, to a Poisson multiple at 1 and to a
gamma with no mass at zero at 2, but the interior holds up: the series' shape
parameter `(2 - p) / (p - 1)` runs from 999 at `p = 1.001` to 0.001 at
`p = 1.999`, and across that range every density is finite, the zero stays the
closed form, and a fit runs. The bound was widened from `[1.01, 1.99]` after the
family was written, so the concern was that the term selection would overflow at
the edges; measured, it does not, and there is a test over the range now.

That widening also broke a test, in a way worth not repeating.
`expect_error(tweedie(power = 1), "1.01")` asserted `arg`'s *wording* of the
bound rather than the contract, so rewording the bound turned a passing suite
into an `R CMD check` ERROR with nothing wrong in the package. The test now
matches on the argument name, which `arg` always includes, and asserts that
`1.001` and `1.999` are accepted and that `1`, `2`, `0.5` and `3` are not. The
general rule: an `expect_error()` regexp should pin down which argument was
rejected, not how the rejection reads.

**Power fixed by default**, at 1.5, which is the one place this family treats a
nuisance parameter differently from the rest. The power and the dispersion are
identified jointly through the share of zeros, the power weakly so at these
sample sizes, and a badly determined power drags the dispersion with it.
`power = NULL` draws it, on a logit scale over `(1, 2)` with the Jacobian that
keeps the prior uniform there, and costs about a third again because the slice
sampler cannot use the table.

**What it was checked against.** `tests/testthat/test-tweedie.R` sums the
compound Poisson-gamma representation directly, sharing no code with the engine,
and the two agree to 1e-8 across powers of 1.2, 1.433, 1.5 and 1.8, dispersions
of 0.4 to 180, and responses from 0 to 60000. The score and information are
checked against their closed forms and against the engine's own numerical
differencing. A fit recovers a known mean function and dispersion, `simulate()`
puts back a point mass of the right size, and the whole suite passes with
`NOT_CRAN=true`, which is 0 failures and 0 skips.

On lalonde `re78` with four chains and 2000 draws it draws `phi = 105`, fits a
mean zero share of 0.219 against the observed 0.233, and its posterior predictive
reproduces the zeros (0.221), the mean (6748 against 6793) and the spread (8277
against 7471). `diagnose()` warns where it warns for every family on that
response: the worst 5% of observations carry 46 effective draws while their
average carries 4537.

**References are wanted and were not written.** The compound Poisson-gamma is
Jorgensen's exponential dispersion models and the series is Dunn and Smyth, but
neither paper is in the Zotero library, and `PAPERS.md` says to leave a citation
out rather than reconstruct one. The `@references` block on `?bartisan-families`
and `references.bib` for `vignette("families")` both want an entry.

## A stronger `variable_importance()`, and what the audit around it found

The table was a bare data frame with four numbers per predictor and no way to
ask anything the summary did not already answer. Four changes, in the order they
matter to a reader.

**`prop_splits`**, each predictor's share of the forest's splitting rules,
computed within a draw before averaging so the shares add to one. This is the
column to compare across fits: `splits` counts rules and so scales with
`num_trees`, and a fit of 40 trees reports several times the count of a fit of 10
without being several times as informative. Measured, it is not invariant either.
The leading predictor's share moved from .51 at 10 trees to .64 at 40 on the same
data, because the forest size also changes how the rules get allocated. So the
claim in the documentation is the weaker true one: it removes the dependence on
`num_trees`, and the ranking is the part that carries over.

**`draws = TRUE`**, which returns the per-draw counts instead of a summary of
them, the way `coef()` and `predict()` already do. This is what makes a
comparison the table cannot express available at all: the posterior probability
that the forest spends more rules on one predictor than another is a proportion
of draws, and two predictors adjacent in the ranking are usually not
distinguishable.

**`plot = TRUE`**, following the `plot` argument on `error_density()`, and a
`plot()` method beside it; see the entry below.

**A `<bartisan_importance>` class with a `print()` method**, which drops the two
interval columns from the displayed table and keeps them in the object, and which
says when the fit had `sparsity = FALSE`. That last one is the point of the class:
`prop_used` is only readable as a selection rule under the sparsity prior, and
without it every predictor keeps a share of the rules and the column sits near 1
throughout. A reader looking at the table had no way to see which fit they had.

### What the audit turned up

The static scan found nothing: no `1:n`, no `sapply()`, no missing
`drop = FALSE`, no float equality, no deprecated *ggplot2* calls in the code as
it stood. Profiling the user-facing calls on a 1500-row fit with 500 draws found
nothing either. Everything is a hundredth of a second except
`predict(newdata = )`, which is 1.37s for 3000 rows and 97% of that is inside
`.Call`, so there is no R-level overhead to remove.

Four real things came out of writing the tests rather than out of reading:

`geom_errorbarh()` is deprecated in *ggplot2* 4.0.0, which the new plot used and
the test surfaced as a warning. Replaced with `geom_errorbar(orientation = "y")`.

`expect_predictor_invariant()` in the test helper claims an invariant that does
not hold for two shapes of fit, and said nothing about either. A `vc()` fit
combines its forests into one predictor on the way out while `eta` keeps them
apart, so the two differ by 3.3 rather than by rounding; a fit with an offset
needs the offset passed to `predict()`. Nothing in the suite had exercised it on
those, so the trap was live for the next person adding a family, which is exactly
how it was found. The helper now says so.

`variable_importance()` and `summary()` compute the same usage numbers from the
same counts in two places. Left as it is, and noted: they agree today, the
summary's matrix is a deliberately different presentation, and a test asserts the
`prop_used` column matches between them, which is what would catch a drift.

The joint update for correlated nuisance parameters, on the To Do list against
the day a two-nuisance family arrived, is now testable because `tweedie()` is
one. Drawing the dispersion and the power in separate slice steps at n = 2000
over four chains of 2000 draws gives R-hat 1.000 and 4735 effective draws for the
dispersion and 1.001 and 4823 for the power, against a posterior correlation
between them of -0.51. Fixing the power raises the dispersion's effective sample
size to 7854, so the separate updates cost about 40% of it and nothing else.
That is not the pathology Linero reports for `sigma` and the generalized gamma's
shape, so the item comes off the list.

## Pre-submission report

Where the package stands as of this pass. Everything below was run rather than
recalled.

### What passes

`R CMD check` on the built tarball, with examples, tests and all nine vignette
rebuilds, reports **Status: OK** with no notes and no warnings. Inside the check,
which is the configuration a CRAN machine runs, the suite is
**0 failures, 0 warnings, 67 skipped, 1560 passing**; the skips are
`skip_on_cran()` on the slow recovery fits, and running them locally with
`NOT_CRAN=true` gives **388 tests, 1933 passing expectations, 0 failures and 0
skips**.

### Blocking

**The version.** `--as-cran` notes that `0.0.0.9000` "contains large components",
which is how it says a development version should not be submitted. `0.1.0`.

### Verification still owed

**`--as-cran` did not complete here, and the reason is probably not the package.**
It hangs at `checking use of S3 registration`, with the subprocess alive at 0% CPU
holding a unix socket and a `com.apple.netsrc` handle, which is a network wait.
The same tarball without `--as-cran` is Status: OK and passes that step, and
`tools:::.check_S3_methods_needing_delayed_registration("bartisan")` returns in
seconds when called directly, so the step's own work is not what blocks. This
needs one run on a machine with unrestricted network before submission, because
the step is where CRAN would find a conditionally registered method that should
not be. The one note `--as-cran` did produce before hanging is the version string
above.

**A second platform.** win-builder and macbuilder have never been run, and
neither has R-devel. The package compiles C++17 and links RcppArmadillo, which is
where a platform difference would show up.

**The spell check.** `--as-cran` runs one and did not get that far.

### Worth doing, not blocking

**Two references for `tweedie()`.** The compound Poisson-gamma is Jorgensen's and
the series is Dunn and Smyth's, and neither paper is in the Zotero library, so
per `PAPERS.md` the prose was written without a citation rather than with a
reconstructed one. `?bartisan-families` and `vignette("families")`'s
`references.bib` both want an entry.

**`vignette("bartisan")` gaps.** It does not cover the bounded gates or either
ordinal augmentation.

**`custom_family()` has no posterior predictive draws**, so `simulate()`,
`pp_check()` and `r2()` are unavailable for a user-written likelihood.
`_dev/SHIP.md` lists this as blocking and this file does not; the case for the
lighter reading is that `custom_family()` is the escape hatch rather than a
supported path, and every compiled family has a sampler. Worth an explicit
decision either way before submission rather than an implicit one.

**The DART inclusion probability is still not stored.** `_dev/SHIP.md` asked for
it on the importance accessor. `prop_splits` is the empirical share of the rules
and answers the same question from the draws, but the posterior of the Dirichlet
`s` itself would be the better column and needs the sampler to record it, which
is a C++ change and a new fit component.

### What is not on the list

The two struck warm-start items, the joint nuisance update, and the silent `NaN`
from `predict(type = "density")`, all for reasons given above the To Do list.

## `plot()` methods for the two functions that had a `plot` argument

`variable_importance()` and `error_density()` were the only two, and both now
have a method as well: `plot.bartisan_importance()` and
`plot.bartisan_error_density()`. `error_density()` needed a class to dispatch on
and returns `<bartisan_error_density>`, which is the same data frame it always
returned with a class on top, so the use `vignette("survival")` makes of the
values is unaffected.

The argument was kept rather than replaced, and **it calls the method** rather
than sitting beside it. That is the whole design decision here: two entry points
to one drawing cannot drift, and a test asserts they produce the same `data` and
the same `labels` rather than trusting that they will. The ggplot objects
themselves carry environments and do not compare equal, which is why the test
compares the pieces.

What the method buys over the argument is that the table can be subset first,
which is what makes a wide model readable: `plot(head(imp, 10))` and
`plot(subset(imp, prop_used > .9))` both work, because `[` on a data frame
subclass keeps the class. The vignette says so and the tests cover both.

The message about *ggplot2* being needed moved to `require_ggplot2()` in
`R/utils.R`. Four places would otherwise have said it in four wordings: two
arguments and two methods.


## Shared forests: what was read, what was built, and what was measured

Linero, Sinha and Lipsitz (2020), *Semiparametric mixed-scale models using
shared Bayesian forests*, Biometrics 76(1) 131-144. Read from the paper and its
supplement rather than from the authors' code, which is MIT licensed and so
cannot be drawn on here.

### The model in the paper

M model components `h_1(x), ..., h_M(x)`, each a sum of T trees, sharing **the
tree structures** and nothing else: the basis functions `psi_t^l(x)`, which say
which leaf an observation falls in, are identical across components, and each
component has its own leaf values in those leaves. So the partitions are the
same and the values are not. Their marginal likelihood factors over the leaves
into a product of one integrated likelihood per component, which is what makes
the Metropolis-Hastings acceptance computable: the supplement writes it as
`L_theta(t, l) * L_u(t, l)` for their hurdle models.

Two things in the paper decide what is worth building. First, the leaf prior may
be multivariate and correlated across components, but need not be: they say
plainly that "the rule-sharing interpretation of our approach still applies even
if Sigma_ij = 0" and that "substantial gains are possible even with
Sigma_ij = 0". So the gain does not come from a correlated leaf prior. Second,
they attribute the gain to variable selection: the informative component "can do
a much better job of selecting the relevant predictors" and the uninformative one
inherits that.

It is worth being clear how this differs from what the package already has.
`bcf()` shares an entire *function* between components, which the paper notes is
a stronger assumption than sharing topology. Sharing topology is stronger than
sharing only a prior over which predictors get used, which is what was built.

### Is there a prize? (`_dev/shared-forests-sim.R`)

Their Section 4 design, rebuilt from the description: a binary Z and a
continuous Y driven by the same Friedman function, conditionally independent
given x, Y informative and Z weak, n = 250, scored by the cross-entropy between
the true and fitted `Pr(Z = 1 | x)` on a held-out thousand. Three arms, of which
the third is the point: `separate` is a probit fit to Z alone with all P
predictors, which is what the package does today and what the paper compares
against, and `oracle` is the same fit given only the five predictors that matter,
which is the ceiling a perfect transfer of variable selection would reach.

| P | oracle | separate | ratio |
|---|---|---|---|
| 5 | 0.0442 | 0.0445 | 1.00x |
| 20 | 0.0488 | 0.0656 | 1.34x |
| 50 | 0.0463 | 0.0755 | 1.63x |
| 100 | 0.0429 | 0.0794 | 1.85x |
| 250 | 0.0447 | 0.1046 | 2.34x |

Twenty replicates. This reproduces the shape of their Figure 2 from scratch: the
oracle is flat in P while the separate fit degrades steadily, so the prize is
entirely a variable-selection prize and it grows with the number of irrelevant
predictors. At P = 5 there is nothing to win, which is the other half of the
result and the reason the feature is off by default.

### What was built: `share_sparsity`

Not the paper's model. `bartisan_control(share_sparsity = TRUE)` pools the
splitting counts of a family's forests and draws **one** Dirichlet from the
total, so every forest reaches for predictors in the same proportions. The
forests keep their own trees, their own cut points and their own leaf scales.
This is the variable-selection content of shared forests without the shared
topology, which the measurement above says is where the gain lives.

The engine change is small, which is the argument for trying it before the large
one. `update_forest()` gained a flag to skip its own `s` draw; `update_shared_s()`
in `mcmc.cpp` pools the counts after every forest has moved and hands one draw to
the rest through `Hypers::copy_s_from()`. Three things it has to get right, all
of them tested:

Only the forests that asked for a drawn Dirichlet take part, since `sparsity` is
a per-forest setting and a forest given the uniform prior must keep it.

The forests that share must be able to split on the same predictors, because a
pooled Dirichlet over different supports is not one distribution. A `vc()` term
holds its coefficient forest to the modifiers, so that case is real and is an
error rather than a silent pooling. It is also narrower than it looks: a forest
allowed only one group does not draw its proportions at all, so a single-modifier
`vc()` term is a no-op rather than an error.

One forest is nothing to share, so it is ignored rather than an error, since a
caller may set this globally and change family.

### The effect, visible in the splitting counts

A homoskedastic response fitted with `gaussian_ls()` at P = 25: the scale forest
has no signal of its own, so what it concentrates on is either noise or whatever
the mean forest found. Its top three by `prop_splits` go from `x7 x13 x23` under
the default to `x4 x1 x2` when the prior is shared, against a truth of `x1`, `x2`
and `x4`. That chunk is in `vignette("families")`.

### Does it help? (`_dev/share-sparsity-sim.R`)

`gaussian_ls()`, n = 400, 20 replicates, each component scored by root mean
squared error against the truth on a held-out thousand. The mean and the log
standard deviation are driven either by the same five predictors, where the
assumption behind sharing holds, or by disjoint sets of five, where it does not.
Ratios are default error over shared error, so above 1 means sharing helped.

| agree | P | `log_sd` | `mean` |
|---|---|---|---|
| yes | 5 | 1.02x | 0.99x |
| yes | 25 | 1.13x | 1.03x |
| yes | 100 | **1.23x** | **1.12x** |
| no | 25 | 0.97x | 0.96x |
| no | 100 | 0.98x | 0.96x |

And `zi_poisson()`, the same design, with the count log mean and the logit of
the zero-inflation probability standing in for the two components:

| agree | P | `zero` | `count` |
|---|---|---|---|
| yes | 5 | 1.04x | 0.97x |
| yes | 25 | 1.30x | 1.03x |
| yes | 100 | **1.40x** | 1.02x |
| no | 25 | 0.95x | 0.91x |
| no | 100 | 1.00x | 0.92x |

The shape is the one the oracle experiment predicted. Nothing at P = 5, because
there is no selection problem to transfer. A quarter off the log standard deviation's error at
P = 100 and two fifths off the zero-inflation probability's, which is the
component with the least signal of its own in each family and so the one with the
most to inherit. The stronger component gains a little too, up to a tenth for
`gaussian_ls()`, which was not expected: a scale forest that has stopped
splitting on noise leaves less for the mean forest to explain.

The cost when the assumption is false is 2% to 8%, several times smaller than the
gain when it holds. That asymmetry is the argument for offering the setting at
all; it is not an argument for making it the default, because a caller who has
not thought about whether the components share predictors has not made the
statement the setting makes.

### The paper's own model, built: `share_forests`

The plan above was carried out, and the note in it about the acceptance ratio was
wrong in a way worth recording, because the same mistake is the natural one to
make and it is invisible in some families.

`bartisan_control(share_forests = TRUE)` gives every forest of a multi-forest
family the same trees. The forests are still one `Tree` object per component, so
`Node`, the flat forest storage, `encode_tree()`, `predict()`'s replay,
`variable_importance()` and the reported per-forest bandwidths all see exactly
what they saw before; what changes is that the t-th tree of every component
holds the same partition, and a move is accepted against the sum of what it does
to all of them. Three new tree moves in `mcmc.cpp` (`shared_birth`,
`shared_death`, `shared_change`), a shared bandwidth move, and
`update_shared_forests()` to drive them. `Node::mirror_birth()` and
`Node::mirror_rule()` take a rule and the division of the support from another
tree of the same shape rather than drawing and computing them again;
`twin_of()` finds the corresponding node by walking the path from the root, so
nothing depends on the order any traversal happens to use.

Which forests take part: every reported one. A varying coefficient's forest as
much as a second parameter's, since both are forests over the same predictors
and what was asked for is one partition rather than one per kind of forest. A
family's pinned nuisance forests are left out, because a forest held at a single
leaf has no topology to share and its branching probability is zero, which would
make a shared birth impossible for the whole group.

#### The likelihood does not separate across components, and the first version assumed it did

The plan said the ratio was the sum over components of
`log_f_after - log_f_before + log_g_reverse - log_g_forward`. It is not, and the
first implementation did exactly that and was wrong.

`Context::log_f()` gives the log density of a node's observations as a function
of *one* component, with the others read from `eta`. Summing those differences
over components, each measured against the values the others held before the
move, is the joint difference only if the log density is additively separable in
the components. A `gaussian_ls()` density is
`-log(sigma) - (y - mu)^2 / (2 sigma^2)`, in which the mean and the log standard
deviation appear together, so it is not.

The fix is to visit the components one at a time and write each into the
predictor before weighing the next, so that component k's difference is measured
against the values components 0..k-1 have just taken. The differences then
telescope to the joint difference exactly. The price is that the predictor is
written during the move, so a rejection restores it from a saved copy of the
node's own entries rather than by never having written it.

That makes the visiting order part of the proposal, which is the second thing
the plan missed. Each component's leaf values are proposed from a Laplace fit
conditioning on the other components as they stand at that moment, so which fit
the reverse move would use depends on the order it visits them in. Ascending for
a birth conditions component k on `{< k split, > k merged}`; descending for a
death conditions it on the same set, which is what makes the pair reversible
with the fits each already computes. A change move is its own reverse, so a
fixed order would need a fit it never computes: it draws the direction, and the
reverse draws the opposite with the same probability, so the two cancel.

**How it was checked.** A temporary build recomputed `total_loglik()` before and
after every shared move and compared it with the accumulated delta, with the
leaf-prior terms that `log_f()` also carries subtracted off. Over 400 moves in
each of seven configurations (all three gates, bandwidth on and off,
`gaussian_ls()`, `zi_poisson()`, a `vc()` fit, and the Dirichlet sparsity prior
on) the largest gap was 8.6e-13. Reinstating the simultaneous commit in
`shared_birth` alone took `gaussian_ls()`'s largest gap to **4.2 log-likelihood
units**, so the check has teeth.

It also would not have caught the bug on `zi_poisson()`, where the gap stayed at
6e-14 with the bug in place: that family's two components are close enough to
separable that the wrong ratio is nearly the right one. A family-by-family
smoke test would have passed. This is the argument for checking the identity
rather than the fit.

The instrumentation is not in the tree. What is: `expect_predictor_invariant()`
over all four gates and both families in `test-share-forests.R`, and the
structural invariants below, which are exact rather than statistical.

#### Decided: not shipping it, and not offering a way to reach it

**2026-09-08.** The model stays in the tree; the argument does not.
`bartisan_control()` no longer has a `share_forests` parameter, there is no
`@param` for it, and `?bartisan_control` and the vignettes say nothing about it.
The engine still implements it in full, and the control list still carries the
field -- it is now read from `the$share_forests`, the package's internal flag
environment, so the only way to engage it is a reference to that environment and
an assignment into it. `tests/testthat/test-share-forests.R` goes through
`with_shared_forests()` in the test helper, and the `_dev/shared-topology-*.R`
scripts through a local `share_forests()` of their own.

The reasoning is in the comment beside the field in `R/control.R` and comes down
to the measurements above: the assumption pays 1.2x to 1.4x on the median
quantity when it is true, costs up to 1.8x when it is false, and saves 1.00x to
1.12x of the run time either way. An asymmetric bet on an assumption a caller
has probably not examined is not a setting.

Two notes for anyone reviving it. The C++ guards still name `` `share_forests` ``
in their error messages, which is the right name to see if the flag is ever
flipped, and a wrong one for a caller to be shown -- so if it is ever offered
again under a different name, those strings move with it. And the flag is read
when `bartisan_control()` builds its list, not when the fit runs, so it has to
be set *before* the control object is constructed; the timing script got that
backwards once and measured two unshared arms without noticing.

#### A scoring bug worth remembering: `set.seed()` keeps the RNG *kind*

The first accuracy run produced a table in which sharing changed nothing
anywhere, with a bias-squared of 2.36 on a mean whose truth had a standard
deviation of 1.08. A fit that misses by more than the target's own spread is not
a fit, so the number was checked rather than reported, and the check was
`cor(fitted, truth)`, which came back **-0.028** with the marginal distribution
of the fitted values matching the truth's almost exactly (sd 1.09 against 1.08,
mean 3.20 against 3.15, range 0.47-5.67 against 0.20-5.60).

Right marginals and no correlation means a permutation, which means the fits
were being scored against a *different draw of the test design*. The cause:

    f <- function(seed) { set.seed(seed); runif(3) }
    f(10000)                      # 0.4419 0.4739 0.3327   Mersenne-Twister
    RNGkind("L'Ecuyer-CMRG"); f(10000)   # 0.3939 0.7027 0.1354

`set.seed(seed)` with no `kind` argument keeps whatever kind is current.
`future.seed = TRUE` puts every worker on L'Ecuyer-CMRG so the streams are
independent, while the main process stays on Mersenne-Twister. The script built
the test set inside the worker and then rebuilt "the same" test set in the main
process to score against, and the two seeds gave two different datasets.

The fix is not to pass `kind =` but to stop regenerating: the test sets are
built once in the main process and handed to the workers, so there is only ever
one draw of them. Training sets are still generated per replicate in the worker,
which is harmless -- their identity does not matter as long as both arms of a
replicate see the same one, and both do.

`_dev/shared-topology-results.R` now reports the fit-to-truth correlation
alongside the decomposition and stops outright below 0.3. The correlation is not
interesting as a result; it is there because it is the one number that catches
this, and no amount of reading a mean squared error does.

The lesson generalizes past this script: **anything regenerated from a seed on
both sides of a `future` boundary is suspect.** Timings and split counts from
that run were unaffected, since neither is scored against a truth.

#### What makes it testable

If the t-th tree of every forest is the same tree then every forest takes the
same rules, so the per-draw split counts and the per-tree bandwidths must agree
to the last bit, and must not agree when the forests are separate. That is
`expect_identical(fit$counts[[1]], fit$counts[[2]])` and the same on the
bandwidth block, and it is what `test-share-forests.R` mostly consists of. Two
forests can agree on totals by accident; they cannot agree on a whole matrix of
counts by predictor by draw.

## `bandwidth_every = 10` is faster and mixes worse, and the default stays at 1

`_dev/bandwidth-every.R` found the candidate at one sample size, one predictor
count and one family: 1.6x faster, effective draws per second slightly up on
smooth means and up by half on step means. That is a range of one, which is the
design I had just faulted Souto & Louzada for, so `_dev/bandwidth-confirm.R`
varied it: twelve cells over n in {250, 500, 1000}, p in {10, 50}, `gaussian()`,
`binomial()` and `gaussian_ls()`, smooth and three-step mean functions, ten
paired replicates each, 240 fits, 1h 22m.

| quantity | `every = 10` against `every = 1` |
|---|---|
| seconds per fit | 1.12x to 1.64x faster; the low end is `gaussian_ls()` |
| worst-quantity ESS | 0.79x |
| ESS per second | 1.26x geometric mean, up in 10 of 12 cells |
| RMSE | 1.2% worse on average, up in 7 of 12 cells |
| 95% coverage | .9713 against .9718 |

The two cells where ESS per second falls are not distinguishable from a wash
(p = .46 and p = .81); three of the ten gains are significant. So the mixing
loss and the speed gain very nearly cancel, and a quarter more effective draws
per second is what is on offer.

What decides it is where the RMSE penalty lands. The two cells with a
significant accuracy difference are `gaussian()` step at n = 250, p = 50 (+4.8%)
and `gaussian_ls()` step (+6.6%), both step functions, which is exactly the case
the bandwidth update exists to handle: the rules have to sharpen toward hard
ones and drawing the bandwidth is how they do it. A default should be safe on
the hard case rather than fast on the easy one, so the default stays at 1 and
`R/control.R` now records the measurement and says when raising it is
reasonable (a mean function known to be smooth, and a compute-bound fit).

## The variable-selection gap was half a gap

`vignette("implementation")` listed "formal variable-selection test" as absent,
filled by a helper package, and named two things other packages do. Reading the
source settles the first half without fitting anything.
`SoftBart::posterior_probs()` is, in its entirety,

```r
varimp <- colMeans(fit$var_counts)
post_probs <- colMeans(fit$var_counts > 0)
median_probability_model <- which(post_probs > 0.5)
```

which is the `prop_used` column of `variable_importance()` and a cut at .5. Both
packages default to the sparse Dirichlet prior that makes the column readable,
so this is the same quantity computed the same way from the same object, and the
only thing missing here was the name: the posterior inclusion probability, and
@barbieri2004's median probability model at the .5 cut. Named now, in
`variable_importance()`'s details.

The permutation test is the real gap. `bartMachine::var_selection_by_permute()`
refits a 20-tree model on 100 response permutations, builds the null
distribution of the splitting shares, and offers three thresholds: pointwise at
the per-predictor 1 - alpha quantile, simultaneous at the 1 - alpha quantile of
the row maxima, and simultaneous at mean + c SD with c bisected for 1 - alpha
coverage. `BART::mc.wbart.gse()` offers the third one only, as a standalone
function taking `x.train`/`y.train`. Neither can read a fit from another
package, so neither integrates. `bartMachine::cov_importance_test()` is a
separate thing again, permuting named covariates rather than the response and
comparing pseudo-R-squared to the null.

Whether it is worth writing is a measurement rather than a judgment, and
`_dev/varsel-check.R` is it: 20 replicates at n = 300, p = 25 with three real
predictors, and a null where none of them matter, crossed with
`sparsity = TRUE`/`FALSE`, comparing the median probability model against all
three permutation thresholds on the same fitted counts. The null arm is the
point: a selector's size is not visible on a design where something is true.

## A per-arm residual variance in `bcf()`: no, and not as an unbuilt family either

The BCF modification that gives each treatment arm its own residual variance
needs no new feature here: `gaussian_ls()` with `log_sd = ~ z` puts the
treatment in the scale forest's formula. What is *not* available is a DPM error
whose scale varies by arm, since `dpm()` carries one additive predictor and one
drawn scale. `_dev/hetero-dpm.R` asks whether either is worth having, with
`dpm()` fitted separately per arm as the stand-in for the family that does not
exist, which bounds above what building one could buy.

Ten paired replicates per cell, $n = 600$, $t_3$ errors, a 3:1 arm SD ratio in
the heteroscedastic cell, truth $\tau = 1 + 1.5 x_1$, 60 fits, 40m 36s.

**Heteroscedastic truth, which is the case the modification is for.** Nothing
separates the three arms. ATE coverage is 9/10, 10/10 and 10/10, which at ten
replicates carries an MCSE of .095 and therefore no information. CATE RMSE is
.206 (`dpm()`), .199 (`gaussian_ls(z)`) and .217 (per-arm `dpm()`), and none of
the three paired comparisons approaches significance (p = .66, .30, .58).
`gaussian_ls(z)` gets its coverage with intervals 31% wider than the single
DPM's, per-arm DPM with 4% wider. The flexible error absorbs the
heteroscedasticity, so paying for the scale forest buys nothing measurable.

**Homoscedastic truth, which is what the default would also have to survive.**
Here there is separation, and it runs against both scale-varying models. Single
`dpm()` beats `gaussian_ls(z)` on CATE RMSE by .037 (7 of 10, p = .014) and
beats per-arm `dpm()` by .068 (**10 of 10**, p < .001), with `gaussian_ls(z)`
still 23% wider. Splitting is the worse of the two because each fit sees half
the data and loses the pooling of the prognostic surface, which is the reason
`bcf()` pools in the first place.

So: no upside in the case the feature exists for, a measurable and in one form
unanimous downside otherwise. Not a default. And since per-arm `dpm()` bounds
above what a heteroscedastic DPM family could deliver, and that bound is the
*worst* CATE performer in the homoscedastic cell, there is no case for building
the family either. A user who wants it can still write
`list(mean = y ~ ... + vc(z), log_sd = ~ z)` with `gaussian_ls()`; the
measurement says what that costs.

One process note. At 6 of 10 replicates the heteroscedastic cell looked like a
clean win for `gaussian_ls(z)` on CATE RMSE, 6 of 6 paired with a mean
difference of .031. Four more replicates took that to 8 of 10 and .007, p = .66.
The sign held and the magnitude did not, which is the usual shape of a partial
result read too early.

## Interval coverage, remeasured: conservative, not deficient

`_dev/coverage-calibration.R`, 40 replicates at $n = 500$ with 10 predictors, at
the package defaults, scoring the 95% interval for the additive predictor
pointwise against the linear predictor the data were generated from.

| family | coverage | MCSE | \|bias\| / posterior SD |
|---|---|---|---|
| `gaussian()` | .964 | .006 | .719 |
| `binomial()` | .961 | .008 | .739 |
| `poisson()` | .973 | .005 | .690 |
| `Gamma("log")` | .970 | .005 | .700 |
| `binomial()`, 4 chains and 5x the draws | .965 | .010 | .717 |

The vignette said .95, .91, .96, .96, blamed the binomial on a binary response
carrying the least information, and put the bias-to-SD ratio near .8. None of
that survives. Coverage is .96 to .97 against a nominal .95, so the intervals
are mildly conservative; the binomial sits .008 below the other three, which is
not a difference at 40 replicates (p = .35); and the ratio runs .69 to .74. The
one claim that does hold is the one about mixing: four chains and five times the
draws move neither coverage nor the ratio (p = .20, p = .18).

The honest caveat is that this is one data-generating process. The mean function
has moderate spread and is mapped to each family's link, which keeps the
binomial's success probabilities informative. A response whose probabilities sit
near 0 or 1 would be a different measurement, and that, rather than anything
here, is the reason to treat pointwise intervals as approximate. The vignette
now says so.

## The permutation test is worth having, and cannot be stacked on the prior

`_dev/varsel-check.R`, 20 replicates at $n = 300$, $p = 25$ with 3 real
predictors, 20 trees, 40 permutations, alpha = .05, crossed with `sparsity` and
with a null design where nothing matters.

Under the global null, family-wise error:

| selector | `sparsity = TRUE` | `sparsity = FALSE` |
|---|---|---|
| .5 cut (median probability model) | **1.00** (7.1 of 25 selected) | **1.00** (25 of 25) |
| pointwise permutation | .90 | .80 |
| simultaneous max | .10 | .10 |
| simultaneous SE | .40 | .10 |

With signal, power over the 3 real predictors:

| selector | `sparsity = TRUE` | `sparsity = FALSE` |
|---|---|---|
| .5 cut | .78 (FWER .35) | 1.00 (FWER 1.00, 13 false) |
| pointwise | .73 (.15) | .87 (.95) |
| simultaneous max | **.12** (.00) | .75 (.30) |
| simultaneous SE | **.37** (.00) | .82 (.45) |

Three things follow. First, the .5 cut is not a test and must not be documented
as one: it fires on every replicate of a pure null. That is not a defect in the
rule, which @barbieri2004 posed as predictive model choice, but it is a defect
in how the package described it, now fixed in `variable_importance()`.

Second, the simultaneous permutation thresholds do hold their size, .10 against
a nominal .05 at 20 replicates, and keep useful power. So the gap is real and
worth closing.

Third, and this is the part that would be easy to get wrong in an
implementation: **the two corrections do not compose.** With the sparsity prior
on, the max threshold's power falls to .12, because the prior shrinks the null
distribution by the very mechanism that shrinks the observed shares, and
thresholding one against the other subtracts the effect twice. Anything built
here has to fit both the observed and the permuted forests with
`sparsity = FALSE`, and say why.

Cost is not the obstacle: 80 replicates of 41 fits each ran in 6m 03s, so a
hundred fits of a 20-tree model is seconds for a cheap family. It scales with
the family, though, and `ordbeta()` at 22.6x a Gaussian fit would not be cheap.

### Recommendation, not built

A `variable_selection()` taking a `<bartisan_fit>`, refitting with
`sparsity = FALSE` at a small `num_trees`, and reporting the three thresholds
with the null distribution attached. Not written: the assessment was the task,
and the design constraint above is the thing worth deciding on before code.

## `bcf()` returns `<bcf_fit>`, and the plan for what dispatches on it

`_dev/plot-methods.md` is the design document: a `causal_effect()` engine
computing estimands on the response scale by per-draw g-computation, a
`summary()` on a BCF fit that calls it with defaults, and plot methods on the
result objects rather than a `what =` switch on the fit.

The class is in. `bcf()` now returns `c("bcf_fit", "bartisan_fit")`, prepended
rather than replacing, so `predict()`, `pp_check()`, the *marginaleffects*
methods and the `arg::arg_is()` guards all keep dispatching. The audit for
exact-class comparisons (`class(x) == "bartisan_fit"`, `identical(class(x), ...)`,
`class(x)[1]`) found none anywhere in `R/` or `tests/`, and `test-bcf.R:208`
passes unchanged because `expect_s3_class()` uses `inherits()`.

One of the plan's four checks is already answered rather than assumed:
intervening on the treatment does not perturb the stored propensity score.
`bcf_newdata_score()` returns identical scores under `transform(d, z = 0)` and
`transform(d, z = 1)`, because the score is a model of $z$ on $X$ and the
intervention touches only $z$. That is what makes g-computation on a BCF fit
correct without special handling, and it wants a test rather than a second
round of reasoning.

The reason the design puts everything on the response scale is worth keeping
here too. `coef(fit)[, "z"]` is a link-scale contrast; on a logit fit the average
of those is the average conditional odds ratio, which is not the marginal odds
ratio. A plot method drawing the coefficient forest would be the easy thing to
draw and wrong in a way that looks right.

## What data augmentation buys, at fifteen replicates instead of three

The vignette's augmentation table was rebuilt twice. The first rerun used three
replicates and had to be thrown away: the speed column's median coefficient of
variation was 3%, but the worst-quantity ESS ratio's was 66%, and within a
single cell that ratio swung by a median factor of **4.5** across the three
replicates, with `ordinal("probit")` soft running .38 and 6.61 on two of them.
Publishing "augmentation triples the effective sample size" off three draws of
a statistic that varies fourfold is exactly the mistake the numbers were being
rechecked to avoid.

The rerun that counts: 15 replicates, `n = 400`, `p = 8`, 50 trees, 500 warmup
and 1000 draws over 2 chains, 480 fits, 2h 38m. Three changes to the design,
each aimed at the variance rather than at the mean:

- a **median-quantity** ESS alongside the worst-quantity one, since a minimum
  over many quantities is inherently high-variance and the median says whether
  the whole chain moved or only its worst corner;
- **ratios of means** rather than means of ratios, because a per-replicate ratio
  of two noisy ESS estimates has a heavy right tail and averaging those tails is
  what produced the 3.81 in the discarded run;
- an **80% bootstrap interval** over replicates, so the table shows its own
  precision.

The verdict, and it changes a claim rather than a default. Augmentation is 1.8x
to 31x faster, median 9.2x. Its mixing cost is **not** universal: the median
worst-quantity ESS ratio is 1.046, and of sixteen cells seven have an 80%
interval entirely below 1, two entirely above, and seven include it. Effective
draws per second favor augmentation in **16 of 16**. So the vignette's "every
one of them trades speed for mixing" was wrong, and it now says what was
measured instead.

Two smaller things. The median-quantity ratios run .52 to 1.15 where the
worst-quantity ones run .37 to 1.58, so the extremes in that column are mostly
the minimum being a volatile statistic rather than the chain as a whole moving
that far. And `negbin()` is the marginal default at 1.1x effective draws per
second, which the vignette now says out loud; the old table had it at 0.5x,
which would have made the default indefensible rather than merely close.

`R/control.R`'s own augmentation table needed no change. It was measured over two
problems and reported as ranges, and its `negbin()` hard row (1.7 to 1.9x speed,
0.61 to 1.14x ESS, 1.2 to 2.0x ESS per second) contains the new point estimates
(1.80x, 0.61x, 1.10x). Reporting a range over problems rather than a point over
one turns out to have been the more durable choice, which is worth remembering
the next time a single-problem number goes into a vignette.

## Three things the full suite caught that a filtered run could not

`estimate_effect()` and the plot methods were built and tested against filtered
runs (`test_local(filter = ...)`), which passed. The full suite then found two
failures and an error, and the interesting part is why the filtered runs missed
them.

**A filtered testthat run is a different run, not a weaker one.** `fit_effect()`
in `test-estimate-effect.R` sets no seed, so its fit depends on whatever RNG
state the preceding files left behind. Run alone the file passed; run after
thirty-five others it did not. Reproduced deterministically with
`set.seed(99); invisible(runif(17))` before the file, which is worth keeping as
the trick for this class of thing.

**The bug it exposed was real and constant, not intermittent.** `effect_focal()`
was edited so that `ATC`'s default focal is the *control* level rather than the
treated one, while `keep` still read `ATC` as the complement of focal. Composed,
those two made `estimand = "ATC"` average over the treated: `z != 0`. ATT and
ATC returned the same units every time, and only whether an assertion could see
it depended on the RNG. The fix follows the edit's intent, which is the cleaner
semantics anyway: **`focal` names the group averaged over in both cases**, and
the estimands differ only in which level it defaults to, the second for `ATT`
and the first for `ATC`. `keep` collapses to one expression, and the identity
that catches this if it ever comes back is already in the tests: the ATE is the
size-weighted average of the ATT and the ATC.

One consequence to note: with three or more treatment levels `ATT` and `ATC` now
differ only in the level named, so `ATC` is `ATT` with the control as focal.
That falls out of the semantics rather than being chosen.

**Two stale test regexes were left by the `R/` reorganization**, both of the same
shape: a hand-written `arg::err()` message replaced by a standardized one, with
the test still matching the old text. `test-methods.R` expected
"must be a fit" where `arg::arg_is()` now says "must inherit from class", and
`test-families.R` expected "not supported by" where the message is now "is not a
supported `family`". Neither is a behavior change; both are the cost of matching
on message text, which is worth paying only where the message *is* the contract.

Final state: 36 files, 427 tests, **2088 passing, 0 failures**.

## Tables moved out of `?bartisan_control`

`@details` on `bartisan_control()` had grown to 637 lines carrying eight tables
of simulation results, and none of them is something a user reads a help page
for. They are the evidence for the defaults, which makes them a developer's
record, so they live here now and the help page keeps only what follows from
them. What a reader needs is the conclusion and the condition under which to
depart from it; what the tables answer is "how do you know", which is this
file's job.

The one exception is the augmentation table, which had a better version in
`vignette("implementation")` by the time it was cut: 15 replicates with a
bootstrap interval against this one's single pass. The help page points there
rather than repeating either.

Reproduced verbatim below, each with the sentence that set it up.

| design | 100 | 200 | paired SE |
| --- | --- | --- | --- |
| Gaussian, soft, n = 1500, p = 10 | -0.001 | +0.001 | 0.002 |
| Gaussian, soft, n = 500, p = 30 | +0.003 | -0.002 | 0.003 |
| Gaussian, hard, n = 1500, p = 10 | -0.001 | -0.002 | 0.003 |
| `ordinal()`, hard | -0.021 | -0.007 | 0.011 |
| `negbin()` | -0.002 | +0.001 | 0.003 |
| `binomial()` | -0.019 | -0.016 | 0.009 |
| `dpm()` | **+0.011** | **+0.006** | 0.002 |

Every design but the last is within a standard error or two of the longer
warmup, and several are better with the shorter one. Effective sample size per
second improves everywhere, by 1.4 to 2.0 times, because the sweeps saved were
producing nothing. `num_draws` went from 500 to 800 at the same time, which
spends some of what warmup gave back on draws that do count towards an
effective sample size; the two together still run in less time than the old
pair did.

| Rules | 5 trees | 10 | 20 | 50 | 100 | 200 |
|---|---|---|---|---|---|---|
| Soft rules | 0.286 | 0.281 | **0.270** | 0.284 | 0.289 | 0.319 |
| Hard rules | 1.149 | 0.682 | 0.558 | **0.521** | 0.531 | 0.510 |

Two things to read off it. **Soft rules need far fewer trees than hard ones**,
which is what makes 200 (the default in most BART packages) actively worse here
than 20. And **the two kinds of rule want different counts**, since hard rules
are still improving at 200 where soft rules peaked at 20.

| `num_trees` | Seconds | Mean RMSE | Log-SD RMSE | Log score |
|---|---|---|---|---|
| `c(50, 50)` | 14.5 | 0.092 | 0.050 | -1188 |
| `c(50, 20)` | 8.3 | 0.094 | 0.047 | -1188 |
| `c(50, 10)` | 5.9 | 0.093 | 0.051 | -1188 |
| `c(50, 5)` | **4.9** | 0.094 | 0.046 | -1187 |
| `c(20, 5)` | **2.9** | 0.084 | 0.041 | -1184 |

A Gaussian fit on the same data takes 1.4 seconds, so `c(50, 50)` costs 10
times a Gaussian fit and `c(50, 5)` costs 3.5 times, at the same accuracy to
three decimal places. That is not the default, because how many trees a
variance surface needs depends on how complicated it is, and silently
under-parameterizing it would show up as intervals that are wrong, which is the
thing `gaussian_ls()` exists to get right. It is worth setting by hand.

| Sparsity | 10 predictors | 50 predictors |
|---|---|---|
| `"none"` | 0.446 | 0.465 |
| `"weak"` | 0.385 | 0.346 |
| `"moderate"` | 0.400 | 0.372 |
| `"strong"` | 0.374 | 0.362 |

For a contrast on a predictor whose signal is weak, any sparsity is actively
harmful, and not only in the atom-at-zero sense above. A binary treatment
among 20 predictors, continuous outcome, residual standard deviation 1,
n = 800, five replicates, with `covers` the share of replicates whose 95%
interval contains the truth:

| true effect | setting | estimate | atom | covers |
| --- | --- | --- | --- | --- |
| 0.05 | `FALSE` | 0.031 | 0.08 | 0.80 |
| 0.05 | `TRUE` | 0.000 | 0.89 | 0.40 |
| 0.10 | `FALSE` | 0.131 | 0.03 | 1.00 |
| 0.10 | `TRUE` | 0.029 | 0.69 | 1.00 |
| 0.20 | `FALSE` | 0.161 | 0.05 | 1.00 |
| 0.20 | `TRUE` | 0.094 | 0.55 | 0.60 |
| 0.50 | `FALSE` | 0.475 | 0.00 | 1.00 |
| 0.50 | `TRUE` | 0.474 | 0.00 | 1.00 |

The prior attenuates a weak effect by half or more and its interval covers
well below its nominal rate. A strong effect is untouched, because the prior
never has reason to drop a predictor that is earning its splits, so this is a
weak-signal failure rather than a general one.

| | predictors | weaker component | stronger component |
|---|---|---|---|
| `gaussian_ls()`, same 5 | 5 | 1.02x | 0.99x |
| `gaussian_ls()`, same 5 | 25 | 1.13x | 1.03x |
| `gaussian_ls()`, same 5 | 100 | **1.23x** | 1.12x |
| `gaussian_ls()`, disjoint | 100 | 0.98x | 0.96x |
| `zi_poisson()`, same 5 | 5 | 1.04x | 0.97x |
| `zi_poisson()`, same 5 | 25 | 1.30x | 1.03x |
| `zi_poisson()`, same 5 | 100 | **1.40x** | 1.02x |
| `zi_poisson()`, disjoint | 100 | 1.00x | 0.92x |

Ratios above 1 are reductions in root mean squared error against the default.
The weaker component is the log standard deviation and the zero-inflation
probability respectively, which are the ones with less signal to find the
relevant predictors from on their own. Sharing buys nothing at five
predictors, because there is no selection problem to transfer, and takes a
quarter to two fifths off the weaker component's error at a hundred. When the
assumption is false it costs 2% to 8%, since the pooled prior pulls each
forest towards the other's variables. Turning it on is therefore a statement
about the data, and one whose downside is a good deal smaller than its
upside. `variable_importance()` is where to check it: with a shared prior the
forests report similar `prop_splits`, and a fit that wants them different will
show that under the default.

| Rule | 10 per level | 25 per level | 100 per level |
|---|---|---|---|
| subset, hard | 0.3705 | 0.2617 | 0.1421 |
| onehot, hard | 0.3998 | 0.2725 | 0.1488 |
| subset, soft | 0.3435 | 0.2265 | 0.1098 |
| onehot, soft | 0.3445 | 0.2201 | 0.1059 |

Under hard rules `"subset"` is better at every size, clearly so at ten
observations per level, where the gap is 7% against a standard error of 2%,
and by 4% at the two larger sizes, where the standard error is around half
the gap. Under soft rules, which is the default, **the two are
indistinguishable**: the largest gap is 0.006 against a standard error of
0.005. Soft rules are worth far more than either choice, which is the biggest
number in the table and the one to act on.

| Family | Rules | Speed | Effective sample size | ESS per second |
|---|---|---|---|---|
| `binomial("probit")` | either | 5.7 to 7.2x | 0.66 to 0.75x | **3.8 to 5.3x** |
| `binomial("logit")` | either | 2.6 to 3.7x | 0.81 to 1.8x | **2.1 to 6.7x** |
| `ordinal("probit")` | soft | 14x | 0.79 to 0.95x | **11 to 13x** |
| `ordinal("probit")` | hard | 26 to 30x | 0.73 to 0.90x | **22 to 24x** |
| `ordinal("logit")` | soft | 7.0 to 7.2x | 0.70 to 0.82x | **5.0 to 5.9x** |
| `ordinal("logit")` | hard | 14.5 to 15.3x | 0.87 to 1.02x | **13 to 15x** |
| `ordinal("cloglog")` | soft | 2.8x | 0.57x | 1.6x |
| `ordinal("cloglog")` | hard | 5.1x | 1.05x | **5.2x** |
| `negbin()` | hard | 1.7 to 1.9x | 0.61 to 1.14x | 1.2 to 2.0x |
| `negbin()` | soft | 1.1x | 0.71x | 0.8x |
| `multinomial()` | soft | 9.3x | 1.09x | **10.1x** |
| `multinomial()` | hard | 14.5x | 0.66x | **9.6x** |
| `zi_poisson()` | soft | 4.6x | 0.85x | **3.9x** |
| `zi_poisson()` | hard | 7.1x | 1.42x | **10.1x** |
| `zi_negbin()` | soft | 5.9x | 0.98x | **5.8x** |
| `zi_negbin()` | hard | 9.4x | 0.84x | **7.9x** |
| `lognormal_aft()` | soft | 11.2 to 19.7x | 0.83 to 1.27x | **14 to 16x** |
| `lognormal_aft()` | hard | 19.7 to 29.3x | 0.74 to 1.01x | **20 to 22x** |
| `loglogistic_aft()` | soft | 8.4 to 9.3x | 0.58 to 0.84x | **4.9 to 7.7x** |
| `loglogistic_aft()` | hard | 11.3 to 12.3x | 0.87 to 0.89x | **9.8 to 11x** |

The ranges are two problems of different size and shape, which is a fair
picture of how much this varies: what an augmentation costs in mixing depends
on the data, not only on the family. The negative binomial is the marginal
case (a clear gain on one problem and a slight one on the other) and is worth
turning off if its diagnostics look poor.
## The negative binomial's soft-rule exclusion did not replicate, and a run that measured nothing

`augment = TRUE` excluded `negbin()` under soft rules, on a single measurement
of 1.1x speed, 0.71x effective sample size and **0.8x** effective draws per
second: a net loss. That was the least-supported default in `bartisan_control()`
after the fifteen-replicate augmentation rerun, which had covered the negative
binomial under hard rules only. Rechecked at the same fifteen replicates:

| cell | speed | ESS, worst | 80% interval | ESS, median | ESS/sec |
|---|---|---|---|---|---|
| soft rules | 1.18x | 1.110x | [0.68, 1.64] | 0.857x | **1.31x** |
| hard rules | 1.80x | 0.611x | [0.47, 0.79] | 0.835x | 1.10x |

The hard cell reproduced the earlier run to three digits, which is what says the
rerun measured the same thing. The soft cell is a modest gain rather than a
loss, and it arrives differently: hard rules buy more time (1.80x) at a real
cost in mixing (per-replicate ESS ratio median 0.52, paired p = .024), while
soft rules buy less time (1.18x) at no measurable cost (p = .53). Both come out
ahead on the ratio that matters, so `augment = TRUE` now means every family with
a rewriting, whatever the rules, and `resolve_augment()` no longer takes `soft`.

### The first attempt measured a model against itself

The run before this one reported 1.0x speed and a 2.06x ESS ratio for the soft
cell, and the speed was the tell: a rewriting that changes nothing takes exactly
the same time. It changed nothing, because the benchmark's "on" arm passed
`augment = TRUE`, which under soft rules resolved to a family list *without*
`negbin` -- the very default being tested. Both arms fitted the same model, and
the 2.06x was two identical configurations differing only in their RNG stream.

Two things worth keeping from that. A benchmark that toggles a default cannot
use the default's own spelling to turn the thing on; the "on" arm has to name
what it wants, which is why the cell now carries `on = "negbin"` and the script
says why. And a ratio of exactly 1.00 on a quantity that should have moved is
worth more suspicion than a ratio that looks wrong: a wrong number invites
checking, while a number that lands on the null reads as a finding.

## A documentation audit, and the duplication list it left behind

Six read-only agents audited the roxygen and the vignettes in parallel: the
`bartisan_control()` `@param` blocks, its `@details`, every functional claim in
every vignette against the source, the causal vignette's structure, what the new
API obliges the other vignettes to say, and duplication between help pages and
vignettes. What follows is what was acted on and what was not.

### `?bartisan_control`, shortened again and stripped of measurements

Details went from 637 lines this morning to 185, with **no simulation figure
left in it**, and the `@param` block from 191 to 158. The arguments were
reordered so the ones adjusted most come first (`chains`, the three chain
lengths, then `num_trees`, `gate`, `sparsity`), in the signature as well as the
documentation, which is safe because every internal call is by name. Two
sections were cut almost entirely because a vignette does the same job better:
"Gaussian Rewritings", which `vignette("implementation")` covers with the
seventeen-row augmentation table and the derivations, survives as a short note
saying which families have no rewriting at all, that being the one fact the
vignette does not carry. Cutting it orphaned the Albert and Chib and Cowles
references, which are removed from the page and live on in `references.bib`.

The rule the trimming followed, in the maintainer's words: give the qualitative
result that informs the choice, not the number from one simulation, and point at
a vignette when the short version is incomplete.

### Twenty-five inaccurate or stale claims, three of them written this session

The claims pass verified every functional statement against `R/` and `src/`
rather than against other documentation. What it found, with the three that were
mine at the top:

- `vignette("causal")` said the predictions are differenced within each draw and
  then averaged. `effect_marginal()` does the opposite: it averages the potential
  outcomes over units first and contrasts them afterward, which is what makes a
  ratio the marginal one. The roxygen said so correctly and the vignette I had
  just written contradicted it.
- The same vignette said `summary()` on a fit reports the potential outcomes.
  Only `summary.bcf_fit` does; `summary.bartisan_fit` reports the forests.
- It also merged two different `focal` behaviors: required with more than two
  levels, guessed with a message at two uninformative ones.
- `@param chains` claimed split-R-hat "is reported in the `rhat` element". No
  such element exists on a fit; `diagnose()` computes it. The same mistake
  appeared twice more in the same file.
- `predict()`'s `type = "response"` restricted the median-survival reading to the
  accelerated failure time families. `ph()` returns one too, by exact inversion
  of the piecewise-linear cumulative baseline. Checked, not assumed: survival at
  the predicted value comes back at .509.
- Three vignettes quoted "500 draws after 500 warmup" as the defaults, which
  stopped being true when they became 200 and 800.
- `ordered_beta()` does not exist; the function is `ordbeta()`.
- `Beta()` was missing from the list of families that accept a composed link,
  though `native_links` includes it.
- A two-column numeric matrix was documented as always binomial; times and 0/1
  events reach `dpm_aft()` first, which `default_family()` confirms.
- `binomial("probit")` was called the fastest family in the package; the
  relative-cost table puts `gaussian()` ahead of it.
- `variable_importance()` was said to return three columns; it returns six.
- `eta.eta` was described as one row of the diagnostics table; there are two, an
  average and a worst-5%, and which one binds depends on the estimand.
- The README's copy of the comparison table had drifted from the vignette's in
  two cells and a footnote.

### The duplication list, mostly not acted on

Forty-one overlaps between roxygen and vignettes, eight of them verbatim. The
five cheapest and clearest were fixed: the lalonde atom-at-zero figures, which
were the same sentences in `R/control.R` and `R/marginaleffects.R`; the
within-chain drift calibrations in `R/diagnose.R`; and three stale pointers that
promised measurements the help pages no longer carry. The rest is a real task
and wants its own pass, in rough order of severity:

- **the atom at zero** has six copies, two of them roxygen-to-roxygen;
  `vignette("effects")` should own it and the help pages should point
- **per-forest arguments** appear in four roxygen blocks, with a byte-identical
  example in `R/bartisan.R` and `R/varying.R`
- **the family table** is in `R/families.R` and `vignette("families")` and will
  drift; worth keeping both, worth a note saying which is canonical
- **the default-family lookup** has three copies of the table and three of the
  same aside
- **the accelerated failure time estimand derivation** has three copies, one of
  them a table rendered as prose
- **`num_bins`** is covered fully in both `R/families.R` and
  `vignette("survival")`

None of these is wrong, which is why none of them was urgent; they are the cost
of a package whose help pages were written before its vignettes.

## The duplication pass, and three changes to how an effect prints

### What the effect object prints

`print()` on a `<bartisan_effect>` now shows the average potential outcomes
beneath the contrast, because a difference of a few points means one thing
against a baseline of .6 and another against .05 and a reader should not have to
reach for `attr()` to see which. `potential_outcomes = FALSE` turns it off.

That made `summary.bcf_fit` redundant, so it is gone. `summary()` on a
`<bcf_fit>` is now the same summary of the forests it is on any other fit, with
a line at the end naming `estimate_effect()`, in the manner of *adrftools*'
`print.effect_curve()`. The point is that the same call means the same thing
whether the model came from `bcf()` or from `bartisan()` with a `vc()` term;
before this, one printed forests and the other printed an effect.

`print()` on a `CATE` object was unusable and is fixed in passing: it dumped one
row per unit, which is 1500 rows on `rhc`. It now reports the quartiles of the
per-unit estimates, which is what `summary.bcf_fit` used to show, and says the
rows are still in the object.

### The contrast label names the quantity

`1 / 0` cannot distinguish a ratio from a log odds ratio, and with numeric levels
it reads as arithmetic on the numbers themselves. Four options were laid out:

- **the comparison in a second column** relocates the ambiguity rather than
  removing it, since `1 / 0` still appears, and spends a column on something
  constant within a call;
- **a formula in the level names**, `log(O(1) / O(0))`, is self-describing and
  needs no legend, but `log(1 / 0)` reads as the log of one over zero, which is
  worse than the status quo in the commonest case;
- **`E[y|z=1]` notation** needs no legend at all and runs to 35 characters,
  repeating the treatment's name in every cell;
- **a formula in symbolic means**, `log(O(Y[1]) / O(Y[0]))`, cannot be misread
  because the bracket marks the level as an index, scales to any level names, and
  costs one legend line per call however many rows the table has.

The last one is in, and it is what *lmw* does. `print()` prints the legend, with
the odds clause only where an odds appears.

### The forest plot

The marginal effect was a band across the panel, which put it underneath every
conditional interval: the one quantity a reader wants to locate was the hardest
to see. It is now a single interval past the right edge in its own color, with a
rule separating it and an axis label. The marks also scale with the unit count,
since at a thousand units a fixed point size merges into a solid block and loses
the spread the plot exists to show.

### The duplication list

Acted on, in severity order: the lalonde atom-at-zero figures that were the same
sentences in `R/control.R` and `R/marginaleffects.R`; the within-chain drift
calibrations in `R/diagnose.R`; the `vc()` naming example that was byte-identical
in `R/bartisan.R` and `R/varying.R`, which now defers to `vc()`; the
default-family lookup table, which `bartisan()` owns since `family` is its
argument and `R/families.R` now points at, keeping only the two deliberate
non-inferences; the accelerated failure time level derivation, compressed to the
claim and the per-family answer with the table left to `vignette("survival")`;
`num_bins`, whose sweep the vignette owns; and the per-forest recycling rule,
which `bartisan_control()` owns while `R/families.R` claims the table of forest
names as canonical.

What remains is the family capability table, in `R/families.R` and
`vignette("families")` in two different orders. Both are wanted, the help page as
a reference and the vignette as an opening, so the roxygen copy is now marked as
the canonical one rather than being cut.

## The binary-outcome check, and a comparison that was described but never run

### `pp_check()` was feeding probabilities' worth of questions to zeros and ones

A binned residual plot and a calibration plot both bin their second argument and
read the outcome within each bin. Replicate outcomes are zeros and ones, so every
bin held one value and every binned mean came back as exactly 0 or exactly 1:
`ppc_error_binned` drew two points per facet and `ppc_calibration` drew a flat
line at the outcome's own mean. Neither errored, and both looked like plots.

The fix is the one *rstanarm* makes for the same two: pass the mean of the
predictive distribution rather than a draw from it. `ppc_calibration()` names
that argument `prep`, and its own family disagrees about whether `prep` or `yrep`
comes second, so it is passed by name. `ppc_loo_calibration()` is not in the set
— it takes replicates and forms the leave-one-out probabilities itself, which is
why it was the one calibration check that had always worked.

The test asserts the binned means are strictly inside the unit interval, which is
what separates probabilities from what used to arrive; before the fix there were
exactly two distinct values and they were 0 and 1.

### `vignette("bartisan")` disclaimed its own check

The fit vignette ran `pp_check()` on a binary outcome with `eval = FALSE` and
then explained that the check says nothing there. It now runs
`type = "loo_calibration"`, which is the check that does say something and holds
each patient out of the probability it judges them against.
`vignette("diagnostics")` owns the reading of it; the intro shows it and defers.
The hand-rolled ten-line calibration plot in the diagnostics vignette is gone
with it, since the package now draws the same thing honestly in one call.

### The logistic regression comparison

`vignette("comparison")` had a section saying a comparison against a Bayesian
logistic regression was worth doing and describing how, without doing it. It now
fits one with `rstanarm::stan_glm()` — precompiled, so it costs about three
seconds and no toolchain — and the result is the useful negative: 16 coefficients
predict this outcome as well as the forest, which is the honest report that the
log-odds are close to linear here. `p_loo` carries the other half, about 17
against the forest's 33.

`loo_compare()` warns that the responses differ. They do not: the hash is
*rstanarm*'s own convention, computed with `digest::sha1()` on its stored
response, and ours is stored as double where theirs is integer, so the hashes
would disagree for identical data. Attaching one would make the warning fire
falsely rather than stop firing, so the vignette explains the warning instead.

### Every print method was writing half its output to stderr

Found by looking at the knitted `vignette("bartisan")` to check that the
potential outcomes had appeared: they had, and `diagnose()` two sections above
showed a bold **What to do** heading with nothing whatever under it.

`cli::cli_bullets()` writes to stderr. `cli::cat_line()`, which `cli_cat()`
wraps, writes to stdout. Every print method in the package mixed the two, so the
tables and headings reached a knitted document and the bullets did not. What was
being lost: every line of `diagnose()`'s advice, the legend naming the levels in
`estimate_effect()`'s contrast labels, the note that a CATE ratio is a
conditional one, `print.bcf_fit()`'s pointer to `estimate_effect()`, and
`variable_importance()`'s warning that `prop_used` cannot be read as a selection
rule under `sparsity = FALSE`. All of it invisible in `capture.output()`, in
every vignette, and in anything else that redirects.

The rule was already written down, in `R/methods.R` above `print_header()`:
print methods use cli's `cat_*` functions because `cli_text()` emits on stderr
and would be invisible to `capture.output()` and to knitr. The `cli_bullets()`
call sites never applied it.

`cli_bullets_cat()` in `R/utils.R` renders the same call through `cli::cli_fmt()`
and cats the result. Rendering is byte-identical — bullets, glyphs, wrapping,
`{.val }` styling — and only the stream changes, which was checked both ways
before replacing the eleven call sites. Indentation inside the message strings
moved with the shorter function name and does not matter, since cli collapses
whitespace; that was checked too rather than assumed.

`test-invariants.R` now asserts that eight print methods write nothing at all to
stderr, and that the two passages that were lost are present on stdout.

### Eight section titles that counted instead of naming

`vignette("survival")`, `vignette("effects")` and `vignette("importance")` had
titles that announced a count rather than a subject: "Two Estimands, Not Three",
"A Comparison Across Six Truths", "Three Traps", "Three Questions". A reader
scanning the table of contents learns nothing from a number, and the surrounding
titles in all three vignettes are descriptive noun phrases, often carrying the
function or argument in parentheses. Two more began with an interrogative, which
reads as a rhetorical setup rather than a heading.

| was | is |
| --- | --- |
| Two Estimands, Not Three | Time Ratios and Hazard Ratios |
| Where the Level Sits | The Level Each Family Reports |
| A Comparison Across Six Truths | The Families Compared by Simulation |
| Three Traps | Details That Are Easy to Get Wrong |
| Three Questions | Predictions, Comparisons, and Slopes |
| The Step for a Numeric Predictor (`variables`) | Choosing the Step (`variables`) |
| What These Estimates Are Not | Descriptions of the Fitted Model |
| What Importance Does Not Measure | The Limits of a Usage Ranking |

Each new title is the sentence the section already opens with: "Predictions,
Comparisons, and Slopes" over a section whose three paragraphs bold exactly those
words, "Descriptions of the Fitted Model" over "Everything here is a description
of the fitted model." The one prose cross-reference, a comment in the setup chunk
of `vignette("survival")` naming the simulation section, moved with them. The
records above this line keep the old titles because they are records of when
those titles were current.

## Three sections of `vignette("comparison")` that named a problem without solving it

### The held-out score was a stub because its punchline was false

The section showed `sum(predict(type = "density", log = TRUE))` on a held-out
split and said it "is the same quantity `elpd_loo` approximates, computed
directly". The two numbers on screen were -847.8 and -172.6, which is not what
the same quantity looks like, and the section stopped there.

They are the same quantity. `type = "density"` averages the draws before taking
the log, which is the form of one pointwise `elpd_loo` contribution, and
`log(colMeans(draws))` equals the `log = TRUE` output exactly. The totals differ
because one sums 1500 terms and the other 300. Per observation: -0.5652 against
-0.5753, the residual gap being leave-one-out's 1499 training rows against the
split's 1200. The vignette prints that table rather than asserting the identity.

What was missing was the use. A log score means nothing alone, so the section now
puts two candidates through the same split and forms the difference with a
standard error from the spread of the per-observation differences, 10.82 +/- 4.94,
which is `loo_compare()`'s arithmetic by hand. It also prices the split: `loo()`
on all 1500 finds the same difference at six standard errors and the split finds
it at two.

### One of the two scale problems has a fix and the other does not

The section named the accelerated failure time against `ph()` case, said the
comparison "can reverse the ordering", and deferred the correction to
`vignette("survival")`. Measured here, uncorrected, the accelerated failure time
model leads by 3390.9; subtracting `log t` from its pointwise contributions, on
events only, puts `ph()` ahead by 131.0. The reversal is shown now instead of
promised. The shortcut, shifting the total by `sum(jacobian)`, agrees to the
digit, because a constant cancels from the importance ratio and leaves the
weights untouched.

Gaussian against `tweedie()` looks like the same problem and is not, which is
worth stating because treating them alike would have produced a false
instruction. Tweedie leads by 1247.4 on `lalonde`. Split at the zeros: at the
positive outcomes the two are within 2 points (-4843 against -4841), and the
whole gap is the 143 exact zeros, where tweedie reports a probability of .2325
and gaussian a density of 3.96e-05 per dollar. No constant relates a probability
to a density, so there is nothing to correct, and the restricted comparison is
the valid one. The question that actually decides it is predictive, not an elpd:
gaussian replicates 0.0000 exact zeros against an observed .2329.

### Tuning and variable selection said what not to do and stopped

Split into two sections, each showing the method rather than only the warning.
Tuning is a grid fixed in advance, chosen with `loo()` on the training half and
assessed on the held-out half, reusing the split from the first section. The
result earns its place: the two orderings disagree, which is the warning's
content rather than a contradiction.

One claim did not survive being checked. "Dropping one predictor from a forest
usually moves the predictive density very little even when the predictor is real"
is too strong: dropping `surv2m` costs 40.5 +/- 9.4 and is plainly detectable.
It is `rhc` that disappears, at 1.79 against a standard error of 2.54, while
`vignette("causal")` puts its effect at six percentage points. The vignette shows
that contrast now, which makes the same point honestly and makes it sharper.

Build cost: `comparison.Rmd` goes from 47s to 105s.

## The atom correction that the previous commit said did not exist

### `loo()` learns the survival measure, and is not allowed to apply it alone

`loo(x, scale = )` and `waic(x, scale = )` take `"time"` or `"log_time"` and put
a survival fit's pointwise densities on the measure named, subtracting or adding
`log t` over the events. A fit already on that scale is returned untouched, so
one `scale` named for every model in a comparison is enough and it reads the same
from either side.

It is deliberately not automatic. Correcting `ph()` on its own initiative would
make `loo()` stop reporting the model's own predictive density: it would no
longer agree with `log_lik()`, and a comparison against a proportional hazards
fit from another package would quietly acquire the error the correction exists to
remove. That trades a visible trap for an invisible one. The argument is the
compromise, with the Details section saying why.

Verified: the argument reproduces the manual `sweep()` to the digit in both
directions, is a no-op on the side that is already right, and leaves the
comparison unchanged at 131 whichever of the two scales is named, since a
constant per observation cancels from a difference.

### Gaussian against tweedie is repairable, and the last commit said it was not

`0a25d3d` claimed "no constant relates the two, so there is nothing to correct".
That is wrong, and the way to see it is that the raw comparison is not invariant
to the units the outcome is recorded in. Refitting `lalonde` earnings in
thousands moves the gap from 1247.4 to 260.1, a shift of 987.3 against a
predicted `143 * log(1000) = 987.8`. A density carries units of one over the
outcome, so all 614 gaussian densities rescale together; the 143 tweedie atoms do
not, because a probability has no units.

The repair follows from the diagnosis. A probability becomes a density when it is
spread over the width the outcome is recorded to, so the atom contributes
`log P(Y = 0) - log(delta)`. Corrected, the gap is 1247 and 1248 in the two unit
systems. In dollars the correction is zero because `log(1) = 0`, which is why the
raw comparison looked reasonable: right by coincidence, not by construction.

`delta` is a fact about the data rather than a tuning knob, and where the zeros
are exact rather than rounded there is no `delta` and the comparison stays ill
posed. That is the case the restricted comparison over the positive observations
still covers.

### What `p_worse` is

`pnorm(0, elpd_diff, se_diff)`, `NA` on the reference row, and at least .5 by
construction because `loo_compare()` sorts before computing it. So .5 does not
mean even odds after weighing evidence; it means the ranking is arbitrary. It
inherits the normal approximation behind `se_diff`, which is what `diag_diff` and
`diag_elpd` flag. Now stated in the vignette.

### Why the convention is a total

A log score is additive, so the sum is a log predictive likelihood and a
difference of two is a log likelihood ratio, which is evidence in nats and ought
to grow with the sample. It is also the scale AIC, WAIC and DIC are quoted on.

The choice costs nothing either way: the mean difference and the total differ by
a factor of n, the standard error differs by the same factor, and the ratio that
decides whether a difference is real is identical, which was checked rather than
assumed. The average is better only when the two halves do not sum over the same
observations, which is the leave-one-out against held-out case and never arises
inside one `loo_compare()` call.

### Structure

All comparison moved out of "When the Approximation Fails", which now stops at
the single-model score and the per-observation table against `loo()`. The
held-out comparison became "The Same Comparison Through the Split" under
"Comparing Two Models", where the price of the split (two standard errors against
six) sits beside the `loo()` result it is being read against.

## Diagnosing the estimand, and an accuracy pass over the remaining vignettes

### The reported quantity mixes worse than anything the fit's table shows

`eta.eta` averaged over observations: R-hat 1.004, 1610 effective draws. The ATE
computed from the same fit: R-hat 1.04, 117. One chain had sat at exactly zero
for 663 consecutive draws, because the splitting prior gave `rhc` no rule and the
contrast of two identical predictions is exactly zero. The fitted function never
stopped moving, the other predictors carrying it, so no parameter the sampler
draws looks stuck and the fit's diagnosis reports nothing. With
`sparsity = FALSE` the same estimand gets 1223 effective draws and R-hat 1.01.

So `diagnose()` is generic now, with a method for `<bartisan_effect>` that
reuses `diagnosis_row()`, `diagnosis_columns()`, the thresholds and the print
method, folds `CATE` over units the way the fit's table folds observations, and
adds one check the fit cannot have: the share of draws sitting at the atom.
`estimate_effect()` stores `chains` and `control` so the method can fold the
draws and name the settings its advice would change.

This complements the fit's diagnosis and does not replace it, since chains that
have settled on different fitted functions can still agree about an average over
them. Both are run in `vignette("diagnostics")`.

The general case is `posterior::summarise_draws()` on anything that can be
arranged as draws by chains, which the diagnostics vignette now shows for a
*marginaleffects* `avg_comparisons()` result. It returns the same four numbers
`diagnose()` gives for the same estimand, which is the point of including it.

Two claims fell to this. `vignette("bartisan")` said "the average row governs an
average effect"; it does not, and no row does. The advice string in
`diagnose()` said an estimand "carries far more effective draws than the table's
worst row does", which invites the same inference.

### `by = ~ x3 > 0` grouped by `x3`

`all.vars()` reduced the formula to the names it mentions, so the expression was
discarded and the grouping ran over a continuous predictor, one group per
distinct value, 250 of them, with no complaint. The right-hand side is evaluated
now, so an expression groups by what it says and a wrong-length one errors.

### The rest of the vignettes, read as a reviewer

`vignette("faq")` still documented `summary.bcf_fit()`, deleted three commits
earlier. The code block ran and produced something entirely unlike what the text
described, which is the worst shape this kind of staleness takes.

The Gamma links were the reverse of the expected direction: `vignette("families")`
was right that every link but `log` is ignored, and `?bartisan-families` was
wrong three times over -- the table row offering `inverse` and `identity`, the
paragraph listing `Gamma()` among the families that compose any link, and the
restricted-range warning naming two links Gamma never reaches. Measured:
`Gamma("inverse")`, `Gamma("identity")` and `Gamma("sqrt")` all come back fitted
on `log`, while `poisson("identity")`, `poisson("sqrt")`, `binomial("cauchit")`
and `Beta("cauchit")` keep theirs. `?predict.bartisan_fit` used the same bad
example for an undefined density and now uses `poisson("identity")`.

Both family tables called the Weibull AFT error "standard Gumbel" where
`src/family.cpp` implements `exp(z - exp(z))`, the smallest extreme value
density; plain Gumbel is conventionally the maximum, and the survival vignette's
own prose already said "smallest extreme value".

A `###` heading in `vignette("families")` had swallowed a paragraph, so the whole
overdispersion discussion rendered as a heading. A comment in `R/families.R` said
a two-valued numeric response falls to "the Gaussian default" when the numeric
default is `dpm()`.

Verified and left alone: the capability table's every link option against every
constructor, the prior-weights rules, the default-family table, all four of
`bcf()`'s settings (`num_trees = c(50, 25)` for the last), `power = NULL`, the
tweedie mean, variance and zero-probability formulas, `variable_importance()`'s
six columns and four printed, the claim that the quantile transform is a step
(it is `stats::ecdf()`), and `ppc_km_overlay`'s ggfortify requirement. Two sweeps
came back clean: every backticked call named in a vignette resolves to a bartisan
export or a Suggests package, and every `vignette()` cross-reference resolves.

## `diagnose()` died on multiple workers, and the error had been dismissed once

Reproduced as `could not find function "diagnosis_block"`. It needs three things
at once: a `future` plan with more than one worker, a pass wide enough to clear
the 100-column floor, and the package loaded with `pkgload::load_all()`.

The convergence pass splits its columns with
`future::future(diagnosis_block(part, chains, step), packages = "bartisan")`.
Written as a bare call, *future* reads `diagnosis_block` as belonging to this
package and drops it from the globals it ships, naming the package in `packages`
instead. Checked directly: `getGlobalsAndPackages()` returns no globals and
infers `packages: bartisan`. The worker then attaches the package, which puts
only its *exports* on the search path, and this function is not one.
`load_all(export_all = TRUE)` attaches the internals to the calling session,
which makes that misreading certain; installed, the lookup resolved anyway. So
the code was relying on a heuristic being right and it was right by luck.

Fixed by naming it in `globals`. Verified under `load_all()` and installed, with
and without a *progressr* handler, on the fit and the `CATE` path, and the
parallel table is `all.equal` to the sequential one. The regression test runs two
real workers and skips where a second cannot be started.

`future.apply::future_lapply()` in `R/bartisan.R` was never exposed to this: it
ships `engine`, a local closure, which future serializes by value.

### This error was seen before and written off

An earlier entry recorded these as "a `pkgload::load_all()` + `future` worker
artifact, not defects -- confirmed by installing and re-knitting (0 errors)".
That was wrong. Installing did make it go away, which is what misled me, but
`load_all()` is the ordinary development workflow and the error there was real
the whole time. A clean re-knit confirmed the vignettes build, not that the code
was sound.

### Whether the manual chunking should be there at all

Asked why this splits chunks by hand rather than calling `future_lapply()`, which
would chunk automatically. The simplification argument holds: it would remove the
`cut()`, the per-chunk steppers, the `resolved()` poll and the explicit `globals`,
and it would have made this bug impossible, since future.apply handles the applied
function's globals -- which is exactly why the chains path never had it. The data
send is answerable by iterating over pre-sliced parts rather than closing over
`wide`.

The progress argument does not hold. The stepper runs on the *worker*, proven by
an instrumentation attempt failing with "could not find function stamp" from
inside `value.Future -> signalConditions`; that stack also shows progressr's
conditions being signalled at collection. So both designs relay in chunk-sized
bursts and the poll only covers the dispatch window. Not measured either way:
three harnesses failed to make worker-side progress relay at all, zero update
events for both designs, so the comments' measured claims about the bar are
unconfirmed and so is any claim that a rewrite would preserve it. That check
wants eyes on a terminal.

## Progress in chunks, and three questions about mixing

### The bar jumped because the workers were in step

Measured separately, which is what separated the causes. Worker ticks are
*generated* spread through the pass (17/46/61/76/100 percentiles), and
`future_lapply` *relays* them progressively -- 22/41/61/80/100 against an
explicit `resolved()` poll's 21/40/60/80/99, so the rewrite cost nothing here.

What jumped was the arithmetic. `PROGRESS_DIAG_TICKS` was 50, divided among the
workers, so on four workers each held about twelve and fired at its own
thresholds; since the workers run in step they crossed those thresholds at about
the same moment and the bar advanced four ticks at a time, an eighth of its
width. The longer each column takes the further apart the jumps, which is why a
fit with many draws showed it most. At 200 ticks the largest gap between reports
falls from 2.8% of the run to 1.3% and the pass costs the same.

Note on method: every "zero progress events" measurement in the previous entry
was worthless, because *progressr* handlers are disabled in non-interactive
R. Those harnesses were measuring a switched-off system. Second time in this
session a conclusion came from a broken harness rather than from checking the
harness first.

### ess_tail above ess_bulk is the normal shape here

Not an error, and systematic: on `rhc` every row has it, `loglik` 23 against
234, `splits.eta` 264 against 731, `eta.eta` averaged 1610 against 2278, its
worst 5% 80 against 305. `ess_bulk` is computed on rank-normalized draws and so
reflects how fast the chain crosses the *centre*; `ess_tail` is the smaller of
the two tail-membership indicators' effective sizes. A slowly drifting level --
what a forest's fitted level does -- is strongly autocorrelated in the middle
while indicators that are mostly zero and flip only on excursions are cheap. So
the pattern says the centre drifts, not that the tails are good.

The ATE inverts it, bulk 117 against tail 64, and splitting the tails says why:
the 5% indicator has ESS 64 and the 95% has 1793. The atom at zero *is* the
lower tail, so `ess_tail` is reporting the stuck atom. Two phenomena, one
column.

### `bcf()` mixing: the propensity score made it worse, not better

At 4 chains and 2000 draws on `rhc`, `propensity = TRUE` gave the ATE bulk 127
and the control forest R-hat 1.32 at ESS 10; `propensity = FALSE` gave bulk 630
and R-hat 1.09 at ESS 34. Five times better without the score, which is the
opposite of the expected direction and worth following up: the score is itself a
fitted BART function of the same covariates, so adding it to the control function
hands that forest two nearly equivalent ways to build one surface, which is the
"fewer ways to represent the same fit" problem the advice already names for
`num_trees`. Seen in one configuration; not concluded.

Neither arm has an atom (0.000 either way), so the `bartisan()` ATE's mechanism
is absent here. Both forests mix far worse than the ATE does, ESS 10 to 48
against 127 to 630, which is the level trade between `f0` and `f1`: only jointly
identified on the treated, so the sampler shifts level between them while the sum
stays put.

### Why ten thousand draws yields a small ESS

Scored on prefixes of one 4-chain, 10,000-draw run, `propensity = FALSE`:

| per chain | bulk ESS | bulk per draw |
| --- | --- | --- |
| 1,000 | 274 | 6.9% |
| 2,000 | 541 | 6.8% |
| 5,000 | 711 | 3.6% |
| 10,000 | 1,200 | 3.0% |

Per-draw efficiency *halves* as the chain lengthens. A short chain cannot observe
autocorrelation at long lags, so its ESS is biased upward; the longer chain sees
the slow mode and reports the more honest number. Absolute ESS still grows, so
the draws are not wasted. At 3% and a single chain, 10,000 draws gives about 300
effective, which is the reported figure. R-hat is 1.003 there, so the chain has
converged and is merely autocorrelated: "mixes slowly", not "fails to mix", and
more draws is the whole remedy.

## Headings underline, and a worked example of fixing mixing with draws alone

### `{.strong}` reads as an accident

cli renders it bold, which in many terminals and in the fonts an editor pane
uses is close enough to the body text that a heading does not look like one. cli
has no `{.underline}` inline class, and an ANSI string embedded in a `{}` slot
does not survive `format_inline()`, which strips it. So `cli_head()` formats
first and styles the result, giving `ESC[4m` where `{.strong}` gave `ESC[1m` and
plain text where the terminal has no ANSI to give it. Ten headings across
`diagnose()`, `estimate_effect()`, `summary()`, `variable_importance()` and
`partial_dependence()`.

### Getting to clean mixing on draws alone

The vignette's list of remedies puts more draws first and notes that the other
two change the model. What it lacked was the demonstration, so `rhc` was measured
until a configuration was found where draws alone genuinely suffice. Three did
not, and the failures were the informative part:

- The full 14-predictor model plateaus. R-hat sits at 1.019 at 8,000 draws and
  1.019 at 12,000, at 269 seconds, so the disagreement there is real and no
  number of draws removes it.
- An `n = 800` subset mixes *worse* per draw, not better: less information, a
  flatter posterior, and the sampler wanders further in it.
- Cutting `num_trees` helps but is a change to the model, which is the thing
  being held fixed.

What works, on the five-predictor model over the whole sample:

| draws | max R-hat | min bulk ESS | warnings |
| --- | --- | --- | --- |
| 200/800 (defaults) | 1.036 | 101 | rhat, rhat readable, bulk ESS, tail ESS |
| 1,000/4,000 | 1.010 | 479 | rhat |
| 2,000/8,000 | 1.008 | 744 | none |

The intermediate row is in the prose rather than the vignette, because it carries
the lesson: the effective sample sizes clear first and R-hat clears last, a hair
at a time, so the final stretch costs the most.

Two cautions went in with it. Effective sample size is non-monotone in draws --
280 at 4,000, 199 at 8,000, 490 at 12,000 on the fourteen-predictor model --
because a short chain cannot see long-lag autocorrelation and reports an
efficiency the chain does not have, so the factor is what to read and not the
figure. And the example says what it does not show: draws sufficed here because
the flagged R-hat was the unreadable kind resting on too few effective draws, and
where R-hat stays put instead, more draws will not help.

Cost: `diagnostics.Rmd` goes from 110s to 194s, which leaves it beside
`bartisan.Rmd` at 198s rather than making it the slowest, so the precompute route
`vignette("survival")` uses was not needed.

## Parallelism where it pays, and the limit nobody had written down

### `partial_dependence()` is the only other site worth it

Every loop in `R/` was read. One qualifies: the grid loop in
`partial_dependence()`, which predicts once per grid point over the whole
sample. Measured on 1500 observations with four workers:

| grid | sequential | 4 workers |
| --- | --- | --- |
| 10 | 93.3s | 32.5s |
| 25 | 242.8s | 72.4s |
| 50 | 350.5s | 129.4s |

A floor was planned, on the model of `diagnosis_columns()`'s hundred columns,
and the measurement says not to have one: even ten points pays 2.9x, because a
point is a whole prediction and costs about nine seconds here. Streams come from
`parallel_streams()` before the branch, so the two paths agree exactly.

Rejected, with reasons, so the next reader does not re-litigate them:

- `interop.R:369`, the Dirichlet process replicate loop, is per-draw and calls
  `sample.int()`, so the RNG order *is* the result and parallelising it would
  change `posterior_predict()` output. Cheap bodies besides.
- The loops over forests and random-effect terms (`bartisan.R:829`, `:838`,
  `diagnose.R:429`, `:450`, `methods.R:127`, `:274`, `interop.R:899`,
  `predict.R:550`) run over one to three components assembling lists.
- `varying.R` and `predict.R` carry most of the package's loops and they build
  design matrices; the expensive part sits under them in C++.

### A fit is 70 MB, and `future` refuses 500

Measured at 1500 observations and 3200 draws: 71.8 MB serialized, of which
`eta` is 38.4 and `forest_flat` 31.0 -- 97% between them. `eta` is draws by
observations, so at 8000 draws and 5000 observations that term alone is 320 MB.

`future.globals.maxSize` defaults to 500 MB and refuses a single export above
it, which now affects three paths: the chains, the convergence pass, and the two
prediction loops. It was documented nowhere. `?bartisan_control` now names all
three axes, the measured speedup, and the limit with its remedy. A user meets
this exactly when the parallelism starts to matter, which is the wrong moment to
meet an undocumented error.

## Weights were one way to flatten a likelihood, and I had mistaken them for the only one

`prior_only = TRUE` zeroes every observation's weight, and the weight multiplies
the log density, the gradient and the curvature at `src/family.h:68`. That
reaches the forests and the leaves. It does not reach an update written against
the response directly, so five families were refused outright. Challenged on
whether that was a mathematical limit or a technical one, and it was mostly
technical.

### The permutation test

The intercept anchor takes the mean and the sd of the response, both of which a
permutation preserves, so under a genuine prior-only fit permuting the response
must leave the draws where they were. Fit twice, once permuted, and read the
largest difference in the `eta` draws:

| family | max abs diff |
| --- | --- |
| `gaussian()` | 7.1e-15 |
| `binomial()`, probit | 0 |
| `gaussian_ls()` | 1.8e-15 |
| `ordinal()` | 0 |
| `multinomial()`, both links | 0 |
| `dpm()` | 4.17 |

Three of the five refusals were wrong. `gaussian_ls()` and `Gamma_ls()` never
leaked; the refusal was a judgment that wide replicates are useless, written as
if it were a correctness claim. Wide replicates are the check working. Measured
at 10 trees, the prior replicates came back at sd 8.1 against the response's
0.84, entirely finite. `ordinal()` does not leak either: its `eta` draws are
identical under permutation and only the cutpoints misbehave.

### What the fix turned out to be

`draw_atom(0, 0, 0, ...)` already reduces to a draw from the
normal-inverse-chi-squared base measure, so the DPM needed no new machinery,
only a question: `bool informative = w(i) > 0.0`. At zero weight the label comes
from the Chinese restaurant prior, the atom from the base measure, and a
weightless observation contributes nothing to the sufficient statistics of the
atom it sits in.

Preferred over threading a `prior_only` flag into the engine, for three reasons.
No plumbing through every family. No second source of truth that can disagree
with the weights. And it fixes a bug that had nothing to do with `prior_only`: a
user passing `weights` containing zeros was having those rows contaminate the
mixture atoms. With unit weights the path is bit-identical to before.

`multinomial(link = "probit")` was never on the refused list and should have
been looked at: `draw_latent()` computed a variance of `1 / (w * prec)`, which
at zero weight is infinite, and the covariance drawn from those utilities was no
longer symmetric. 240 `inv_sympd()` warnings per fit, now none.

### The residue is the RNG stream, not the data

After the fix `dpm()` is not bitwise identical under permutation. It is not a
leak:

| comparison | max abs diff |
| --- | --- |
| two responses, same seed | 0.12 to 0.45 |
| two seeds, same response | 5.3 to 126 |

KS on the two `eta` samples gives p = 0.98. The mixture consumes random numbers
as it goes and the two streams drift apart by a hair. The test asserts the drift
is under a quarter of the prior's own standard deviation, which it clears by an
order of magnitude and which the pre-fix code failed outright.

### `ordinal()` is the one that stays refused

Its cutpoint target is `sum_i w(i) * ordinal_log_prob(...)` at
`src/family.cpp:549`, with no prior term. At zero weight the target is not
flattened but empty, and the slice sampler walks a flat improper density toward
the 1e4 clamp the code puts on the top cutpoint: `cut2` reached 311.6 in 100
draws, and every replicate lands in one category.

That is a gap in the model rather than in the mechanism. It closes the moment a
prior over ordered cutpoints is specified, which would change every ordinal
posterior and so was not slipped in. `multinomial()` is the substitute offered
for an ordered outcome with few enough categories.

## `Matrix` is one call and stays, because the sparsity is the point

Asked whether the dependency earns its place, `Matrix` being in Imports for a
single call: `Matrix::sparseMatrix()` in `make_group_probs()`, which maps
design-matrix columns to formula terms. It earns it twice.

The matrix is columns by terms with exactly one nonzero per row, each design
column belonging to exactly one term, so dense storage is quadratic in the
predictor count:

| p | sparse | dense |
| --- | --- | --- |
| 10 | 2.4 KB | 1.9 KB |
| 100 | 9.3 KB | 84.8 KB |
| 1000 | 79.6 KB | 7.7 MB |

Dense wins only up to about ten predictors, which is to say only on models
small enough not to care.

The second reason is the hotter one. `sample_class_col()` in `src/utils.cpp`
picks a design column within a chosen term by walking that column's nonzeros
with an `sp_mat::const_col_iterator`, which is one step for a scalar predictor.
Dense would scan every row of the column instead, turning a constant into
O(p) on every split proposal.

Removing it would mean changing `bartisan_fit()`'s signature from
`arma::sp_mat` to `arma::mat`, recompiling, and accepting both regressions, in
order to drop a package that has `Priority: recommended` and therefore ships
with every R installation. Not a saving.

## A documentation audit: what ten clean-building vignettes were still getting wrong

The vignettes all built, and the build was not the check. A pass over the ten
against the source found 106 defects that `R CMD check` cannot see, because
nothing validates the prose in a markdown file against the code it describes.
The mechanical checks were already clean: no Rd markup in any `.Rmd`, no ` -- `
anywhere, every `vignette()` cross-reference resolving, all 50 citation keys
resolving in `references.bib`. The defects were in what the sentences claimed.

**Five vignettes could not build without a `Suggests` package.** The worst is
`families.Rmd`, where the chunk guarded on *loo* was the only place `n` and `d`
were created and six later unguarded chunks used them, so without *loo* the
vignette died at the next chunk with `object 'n' not found`. `causal.Rmd`'s
`bal.tab()` chunk carried no chunk options at all. `diagnostics.Rmd` called
`posterior::rhat()` under `eval = run`, where `run` is hard-coded `TRUE` and
`has_post` was sitting unused two lines above. `comparison.Rmd` called
`rstantools::posterior_predict()` guarded on *cobalt*. `survival.Rmd` defined
`has_surv`, spent it on three chunks, and called `library(survival)` unguarded
before any of them. The pattern is the same each time: the flag exists, so the
hazard was noticed, and the guard landed on the wrong chunk.

**An inline expression had been silently deleting a sentence.**
`survival.Rmd` read `res$meta$num_draws` where the field is `num_save`. `$` does
not partially match on a list, so `sprintf()` received `NULL`, returned
`character(0)`, and knitr's inline hook rendered that as the empty string. The
whole opening sentence of "The Families Compared by Simulation" was absent from
the built vignette and had been for as long as the field was misnamed. Nothing
warns. This is the argument for reading the rendered HTML rather than the source
when checking a vignette.

**Three claims about the package were wrong.** `bartisan.Rmd` said out-of-range
predictions are "shrunk toward the overall mean"; `range_map()` clamps with
`pmin(pmax(., 0), 1)` and `ecdf_map()` saturates, so they are held flat at the
boundary instead. `families.Rmd` listed `negbin()` among the families taking a
composed link, where `negbin()` validates with `arg::match_arg(link, "log")` and
the roxygen names it as one of the two exceptions. `survival.Rmd` sent a reader
wanting interval censoring to `custom_family()`, which cannot take a censored
response at all: `logdens` gets `y` as a length-`n` vector, and a `Surv` response
reaches the engine flattened, where it dies with `Mat::operator(): index out of
bounds`. That crash was a real gap and is fixed below.

**The largest class was prose contradicting the chunk printed directly above
it.** In `importance.Rmd` the straggler noise predictor is `x9` at .22, not `x8`
at .4, and `x8` is in fact the lowest of the ten; both `x1` and `x1_copy` sit at
`prop_used` 1.000, where the text said `x1` was used in under 10% of draws; and
`avg_comparisons()` gives the original .347 where the text said "none". In
`bartisan.Rmd` the partial dependence curve runs .870 to .480, not .88 to .45.
In `families.Rmd` the `loo_compare()` output puts `negbin()` last, not
`zi_negbin()`. Every one of these is a seeded chunk whose output moved when
something upstream changed, with the sentence beneath it left behind.

The lesson is not to check the numbers more often. It is that **a sentence should
not commit to a digit a seeded chunk recomputes** unless the digit is the point.
The rewrites read the tables qualitatively ("four of the five noise predictors
near zero. The fifth sits above the other noise variables and still well below
the real ones") and cannot go stale the same way.

Also done in the pass: `README.Rmd` rewritten from scratch, from about 1700 words
to 644, dropping the feature-comparison table that duplicates
`vignette("implementation")` and replacing the syntax tour with one analysis of
`rhc` carried from a binomial fit through `estimate_effect()` to the same event
as a censored time. A `reference:` index in `_pkgdown.yml`, five sections by use,
all 18 non-internal topics, `check_pkgdown()` clean. `bartisan.Rmd` moved from
second person to the first person plural the other eight tutorials use, and its
headings and every vignette title to Title Case. Three spellings of
`marginaleffects_posterior_center` normalized to the one `?bartisan-marginaleffects`
documents. `?bartisan-package`'s `@seealso` had been sending readers to
`vignette("bartisan")` "for how the sampler works", which is
`vignette("implementation")`.

Verified by rebuilding all ten: 13m 12s, ten clean, zero warnings or messages
captured. The duplicate `diagnose(fit)` chunk removed from `diagnostics.Rmd` is
most of the 3m 34s the build lost.

## A `Surv` is a numeric matrix, which is how four families and the engine each failed differently

`check_numeric_response()` dropped a one-column matrix and then called
`as.numeric()`. A `Surv` object is a numeric matrix with two columns, so it
passed the `is.numeric()` check and came out flattened to length `2n`, with the
event indicators appended to the times. Every family that takes one value per
observation funnels through that function, so the same mistake surfaced four
ways and none of them named it:

| Family | What the caller saw |
|---|---|
| `gaussian()` | `'x' and 'w' must have the same length`, from `weighted.mean()` |
| `Gamma("log")` | "requires a strictly positive response", flagging the zeros in the status column |
| `Beta()` | "requires a response strictly between 0 and 1", same cause |
| `custom_family()` | `Mat::operator(): index out of bounds`, from the C++ engine |

The last is the one that mattered, because `vignette("survival")` had been
recommending `custom_family()` as the route to an interval-censored likelihood.
It is not a route: `logdens` receives `y` as a length-`n` vector and has nowhere
to put a censoring indicator. The vignette now says so.

The fix is six lines in `check_numeric_response()`, refusing a response with
more than one column before the coercion rather than after it, with the hint
keyed on whether the response is a `Surv` (name the five survival families) or a
plain two-column matrix (name `binomial()`). Fixing it at the funnel rather than
in `custom_family()` is what makes all four cases report the same thing.

Tested in `test-families.R`, over all five families and both hint branches, plus
the one-column matrix that must still come back as a vector.

## `collapse` is in Suggests for *marginaleffects*, not for us

Nothing in `R/`, `tests/` or `vignettes/` references it, so a dependency sweep
reads it as dead weight and proposes dropping it. It is not: *marginaleffects*
uses *collapse* for its Bayesian calculations, which is the path every
`<bartisan_fit>` takes through `avg_comparisons()` and its relatives. Removing
it would degrade the estimands three vignettes lead with. Recorded here because
the sweep will run again.

## Blocking grid points in `partial_dependence()`: measured, and there is nothing to amortize

`partial_dependence()` calls `predict()` once per grid point, over the whole
sample each time. The obvious optimization is to stack several grid points into
one data frame and make fewer, larger calls. Measured, it does not pay, and at
the sizes where a partial dependence plot is actually slow it costs time rather
than saving it.

**The per-call cost of `predict()` is zero.** Timing one call against row count,
minimum of five runs, 50 trees and 800 draws at 10 predictors:

| rows | elapsed | per 1000 rows |
|---|---|---|
| 1 | 0.0010 s | 1.00 s |
| 10 | 0.0070 s | 0.70 s |
| 100 | 0.0680 s | 0.68 s |
| 500 | 0.3490 s | 0.70 s |
| 2000 | 1.4090 s | 0.70 s |

Regressing elapsed on rows gives a fixed cost of **-1.3 ms**, which is zero
within noise, and a marginal cost of 0.705 ms per row. There is no per-call
setup to spread over a larger batch, because the forest is walked per row and
nothing is rebuilt per call.

**That is where the time goes.** `Rprof` on a 25-point plot at n = 1500 puts
**98.7%** of the run inside `.bartisan_predict()`, the C++ forest walk. The R
side that blocking would actually remove (one `model.frame()`, one
`model.matrix()`, one data frame copy per grid point) is under 1% of the total.
Stacking B grid points does the same forest work on the same number of rows.

**End to end it goes the wrong way at scale.** 25 grid points, 50 trees, 500
draws, stacking all 25 against looping over them:

| n | 25 calls | 1 stacked call | ratio |
|---|---|---|---|
| 200 | 2.35 s | 1.79 s | 1.31x |
| 1000 | 11.30 s | 9.84 s | 1.15x |
| 5000 | 54.54 s | 63.74 s | **0.86x** |

The gains at small n are fractions of a second on runs that were already fast,
and they are within the run-to-run variation of a single-shot measurement. The
loss at n = 5000 is nine seconds, and the reason is the draws matrix: stacking
is `n * grid * num_draws` doubles, 500 MB at that size, against 20 MB unstacked.

Two further reasons not to, beyond the timings. The result is discarded almost
immediately: `at_grid_point()` takes `rowMeans()` of the draws and keeps a
length-`num_draws` vector, so a bigger matrix buys nothing downstream. And each
grid point currently gets its own RNG stream, deliberately, so that a family
whose prediction simulates gives the same draws whether or not a `future` plan
is set; a block would share one stream across its points and change those
results.

**What actually helps is skipping trees**, and that is the entry below. The
predictor is a sum over trees, and a tree that never splits on the plotted
variable contributes the same value at every grid point, so it is evaluated once
rather than `grid` times. The ceiling is large: on an `rhc` fit with 50 trees, an
upper bound on the share of trees using each predictor is 37% for `age`, 31% for
`surv2m`, and under 15% for the other twelve.

Benchmarks run on one core (`parallelly::availableCores()` reports 1 on this
machine), so the parallel path was not measured. It would not change the
conclusion: `future_lapply()` chunks tasks per worker, so the fit is exported
once per worker rather than once per grid point, and blocking removes no
serialization.

## Tree skipping in `partial_dependence()`: 2.6x to 12x, and exact

The entry above measured blocking and found nothing to amortize, because
`predict()` has no per-call cost worth spreading. This is the version of the
idea that works: cut the number of *trees* evaluated rather than the number of
calls.

The additive predictor is a sum over trees. A tree that never splits on a
plotted column returns the same value however that column is set, so it
contributes the same amount at every grid point. Those trees are evaluated once,
on the data as it stands, and only the rest are re-evaluated per grid point and
added to them. Measured on the `rhc` fit at 50 trees, 800 draws and 25 grid
points, against the path it replaces:

| variable | trees that move | before | after | speedup |
|---|---|---|---|---|
| `surv2m` | 32.6% | 24.95 s | 9.72 s | 2.57x |
| `age` | 22.1% | 25.77 s | 7.06 s | 3.65x |
| `paco2` | 15.3% | 26.51 s | 5.61 s | 4.72x |
| `aps` | 11.1% | 27.03 s | 4.48 s | 6.04x |
| `resp` | 4.9% | 28.23 s | 2.87 s | 9.84x |
| `meanbp` | 3.0% | 27.84 s | 2.33 s | 11.96x |

The saving tracks the share of trees that move, which is what it should do. The
curves are identical to the unoptimized ones, not close to them: the largest
difference across those six is `1.11e-16`, and four of the six are exactly zero.

**The pieces.** `.bartisan_tree_uses()` walks each stored tree's records with
the cursor `eval_tree()` uses, taking both subtrees unconditionally since this
is a property of the tree rather than of an observation, and reports whether any
internal node splits on a given set of columns. Rules carrying a level mask index
the level-code matrix and the rest index the design matrix, so the two kinds of
column are looked up in separate sets. `bartisan_predict()` gained a `tree_mask`;
a zero-length mask evaluates everything, so the ordinary path is untouched.

On the R side, `eta_to_type()` was split out of `predict.bartisan_fit()`, which
is everything that happens once the predictors are in hand. Partial dependence
calls it rather than a copy of it, so a grid point's number is the one
`predict()` would have returned. `predict_eta()` gained `tree_mask` and a
`constants` flag: the intercept, the offset and the random effects belong to the
predictor once rather than to any subset of trees, so the base carries them and
the moving part does not.

**Where it declines, and why each case is the way it is.** `pd_tree_mask()`
returns `NULL` and the ordinary path runs when a plotted variable is a `(1 | g)`
grouping factor, because those intercepts are drawn parameters rather than trees
and no choice of trees holds them fixed. Everything else uncertain is marked
*moving* instead, which is only slower and never wrong: a term label that will
not parse, and a `bcf()` propensity score, which is rebuilt from the data at
every grid point and so moves whether or not the plotted variable feeds it.

Two cases needed care. A `vc()` covariate multiplies its own forest, so the part
that does not move still reaches the prediction scaled by something the grid
does; combining after the sum rather than before is what keeps that right. And
when no tree splits on the variable at all, the whole forest lands in the base,
the grid costs one evaluation in total, and the curve comes out exactly flat
rather than flat to rounding.

A term is matched to a variable through `get_varnames()` on the label rather
than by string equality, so `log(x1 + 1)` is correctly marked as moving when the
grid is over `x1`. Keying on the label alone would have held that curve flat with
nothing to show it had.

**What `...` reaches.** `partial_dependence()` forwards `...` to `predict()`, and
the fast path assembles the predictor itself rather than calling `predict()`, so
it reads what it needs by name: `offset`, `iterations`, `weights`, `values`,
`log` and `times`. Any other name, or a positional argument there is no name to
match, hands the whole grid back to `predict()`. The alternative was to accept
the argument and not apply it, which would have been wrong in silence; this is
slower on those calls and cannot be.

## `x_transform`: neither default is dominant, and `range` has the worse tail

Asked whether `range` would do as a default, since it is the option that leaves
a soft-rule fit differentiable and so makes `avg_slopes()` mean something.
Measured rather than argued: `_dev/xtransform-sim.R` fits both transforms to
eight data-generating processes at n = 500 and 2000, 20 replicates each, and
scores the posterior mean of the regression function against the **true** function
on 1000 held-out points. `_dev/xtransform-results.R` prints the tables. RMSE
below is a fraction of the signal's own standard deviation, so it is comparable
across scenarios.

| scenario | n | quantile | range | ratio | range better in |
|---|---|---|---|---|---|
| uniform | 500 | .0593 | .0544 | .92 | 17/20 |
| uniform | 2000 | .0322 | .0310 | .96 | 13/20 |
| lognormal | 500 | .0627 | .0573 | .91 | 16/20 |
| lognormal | 2000 | .0329 | .0314 | .96 | 13/20 |
| outliers | 500 | .0695 | .0520 | .75 | 20/20 |
| outliers | 2000 | .0365 | .0310 | .85 | 18/20 |
| pareto | 500 | .0652 | .0567 | .87 | 19/20 |
| pareto | 2000 | .0340 | .0327 | .96 | 12/20 |
| **outliers_fine** | 500 | .1788 | **.7116** | **3.98** | 2/20 |
| **outliers_fine** | 2000 | .1423 | **.8753** | **6.15** | 0/20 |
| sparse_tail | 500 | .0723 | .0551 | .76 | 20/20 |
| sparse_tail | 2000 | .0343 | .0304 | .89 | 14/20 |
| bimodal_linear | 500 | .0922 | .0331 | .36 | 20/20 |
| bimodal_linear | 2000 | .0496 | .0187 | .38 | 20/20 |
| mixed | 500 | .2091 | .1303 | .62 | 20/20 |
| mixed | 2000 | .1028 | .0639 | .62 | 20/20 |

**`range` is better in seven of the eight**, by 4% to 64%, and the two scenarios
built to punish `quantile` behave as predicted: with two tight clusters and a
truth linear in $x$, the empirical distribution function jumps .5 across the gap,
linear in $x$ becomes a step in $u$, and `range` is nearly three times more
accurate.

**And in the eighth it fails outright.** On `outliers_fine` it is four times
worse at n = 500 and six times worse at n = 2000, with an RMSE near 1.0, which
against a signal standardized to 1 means it is fitting nothing at all. Coverage
of the 95% interval goes .392 and then **.217** against a nominal .95, so the
failure is not visible as extra uncertainty; the intervals are nearly three times
wider than `quantile`'s and still miss.

**It is the cutpoint prior, not the bandwidth.** The first guess was that the
gate width, being a fixed fraction of the transformed scale, over-smooths the
bulk. It is not: fixing the bandwidth, and setting it to .01, leave the failure
untouched (1.006, 1.005, 1.005). What does it is `Node::draw_rule()`, which draws
a cutpoint uniformly on the node's live range. With 1% of the data at 300 times
the scale, the central 98% of $x$ spans **1.34%** of the range at n = 500 and
**.42%** at n = 2000, so a proposal almost never lands where the structure is and
the forest cannot resolve `sin(3x)` there.

That also explains the direction with $n$, which is the tell. The maximum of $n$
draws from a wide component grows with $n$, so the bulk occupies an ever smaller
share of the range and `range` gets **worse** as the sample grows. Every other
scenario improves with $n$ for both transforms.

**Outliers alone are not the condition.** `outliers` and `pareto` are heavy-tailed
too and `range` wins both. The difference is what the truth asks for: those use a
saturating logistic that needs one or two cutpoints across the whole bulk, where
`outliers_fine` oscillates and needs many. The failure needs a compressed bulk
**and** fine structure inside it.

**Conclusion: keep `quantile`.** Not because it is more accurate, since it
usually is not, but because the two have different worst cases. Across these
eight, `quantile` is at worst 2.8 times off the better option with coverage
intact; `range` is at worst 6 times off with coverage at .217, and that case
deteriorates with more data rather than improving. A default is a bet on the
unseen dataset, and bounded downside is the right bet. It also keeps the package
on the same footing as `SoftBart:::trank()`, which is the same rank transform.

`range` remains the right choice when the predictors are known to be free of
extreme outliers, and it is the only one of the two that supports a derivative.
Worth saying in `?bartisan_control` that the choice has a failure mode in each
direction rather than being a matter of taste.

## A smoothed CDF beats both current `x_transform` options on the tail

The entry above left `quantile` as the default because `range` fails badly on
one scenario. The obvious next question is whether the two failures can be
avoided at once, since they have opposite causes: `range` puts the cutpoint
prior in the wrong place, `quantile` puts the coordinate on a step. Those are
separable, and separating them is what the candidates below do.

`_dev/xtransform-candidates.R` runs five arms over the same eight processes, 12
replicates, n = 500 and 2000. `robust` is a logistic squash on a median/MAD
scale, `smoothcdf` a kernel-smoothed empirical distribution function with the
bandwidth at Azzalini's (1981) $n^{-1/3}$ rate for a distribution function
rather than the $n^{-1/5}$ rate for a density, and `winsor` caps at the 1st and
99th percentiles. All three are monotone maps to a bounded interval, which is
what `make_unit_map()` already returns, so each was prototyped by transforming
the column and fitting with `x_transform = "range"`: `range` on a monotone image
is an affine rescale and a cutpoint uniform on a node's live range is affine
equivariant, so the prototype is the model rather than an approximation of it.

The decision-relevant table is not the per-scenario accuracy but how bad each
arm gets, as a ratio to whichever arm won that cell:

| arm | median ratio | worst ratio | mean coverage | where it is worst |
|---|---|---|---|---|
| quantile | 1.18 | 3.01 | .951 | bimodal_linear, n = 500 |
| range | **1.00** | **6.65** | .873 | outliers_fine, n = 2000 |
| robust | 1.13 | 4.23 | .953 | mixed, n = 2000 |
| **smoothcdf** | 1.12 | **1.42** | **.960** | mixed, n = 500 |
| winsor | 1.05 | 4.04 | .956 | mixed, n = 2000 |

**`smoothcdf` dominates the current default on every summary**: better typical
accuracy (1.12 against 1.18), a worst case less than half as bad (1.42 against
3.01), and better coverage (.960 against .951), with no cell below .922. Against
`range` it trades 12% of typical accuracy for a worst case that is 4.7 times
less bad and coverage that never collapses.

**Capping and saturating are the wrong fix, and the reason is instructive.**
Both `robust` and `winsor` repair `outliers_fine`, and both then break on
`mixed`, at .30 against `range`'s .07. `mixed` has a truth containing
$0.1 x_2$ where $x_2$ carries the outlier component, so the extreme values hold
real signal; winsorizing discards it and a logistic squash saturates it away.
`robust` additionally fails on `pareto` (.061 against .032), because a
median/MAD logistic is symmetric and a long one-sided tail is not. A strictly
increasing smooth map compresses the tail without destroying what is in it,
which is the property that separates `smoothcdf` from the other two.

**It also reaches what a cutpoint-side change would reach.** Drawing a cutpoint
uniformly in a coordinate $v = F(x)$ and mapping back gives $c = F^{-1}(u)$, a
quantile-spaced cutpoint in $x$. So a smooth invertible CDF coordinate places
cutpoints where "uniform on the quantile scale, coordinate on the range scale"
would, and differs only in the metric the soft gate measures width in. That
version needs `Node::draw_rule()` rewritten in C++; this one is about fifteen
lines in `make_unit_map()`.

**Cost is not an obstacle.** Naive evaluation is $O(n)$ per point, which at
n = 10000 is .27 s per predictor per call. Tabulating on a 512-point grid and
interpolating, which is how `ecdf_map()` already works, costs 3 to 44 ms once at
fit time and is then free, with an interpolation error of about `1e-3` in the
unit coordinate. A finer grid buys more if that matters.

Not implemented. Adding a third option to a public argument, and any move of the
default, is a decision rather than a measurement. Two things to weigh first: the
arm is never the best in a single cell, only the best worst case, which is the
right property for a default but worth saying out loud; and while a smoothed CDF
restores a usable derivative, the chain rule then multiplies by the density
estimate, so `range` remains the better transform when a derivative is the point.

## `smoothcdf` implemented and made the default

`make_unit_map()` gained a third type and `bartisan_control()`'s `x_transform`
now defaults to it, on the evidence in the two entries above: it dominates the
old default on typical accuracy, worst case and coverage at once, and bounds the
failure that kept `"range"` from being the default instead.

`smooth_cdf_map()` convolves the empirical distribution function with a normal at
the bandwidth rate Azzalini (1981) gives for a distribution function, which is
`n^(-1/3)` rather than the `n^(-1/5)` that is right for a density. The scale is
the smaller of the standard deviation and the interquartile range over 1.349, so
a long tail sets the bandwidth from the bulk rather than from itself. The result
is tabulated on a 1024-point grid over the observed range and read back through
`approxfun()`: evaluating the sum directly costs `O(n)` a point, which is .27 s
per predictor per call at n = 10000, where the table costs 3 to 44 ms once and is
then free. Linear interpolation of an increasing function is increasing, so the
one property everything rests on survives the tabulation.

Normalized on the observed range rather than on the grid, which was the first
version: the kernel puts mass below the smallest observation and above the
largest, so normalizing on the smoothed values left about 3% of the coordinate
where no cutpoint could ever be useful. With outliers at 300 times the scale the
central 98% of a predictor now occupies 98% of the coordinate, against 1.3%
under `"range"`, which is the whole of the fix.

**What a slope costs under it.** Writing the fit as `f(T(x))`, a slope is
`f'(T(x))` times `T'(x)`, and the three transforms differ in what is estimated
there. Measured on a surface whose average slope over the sample is .5485, as the
numerical step shrinks from 1e-1 to 1e-5:

| step | smoothcdf | quantile | range |
|---|---|---|---|
| 1e-1 | .553 | .577 | .544 |
| 1e-2 | .559 | .882 | .543 |
| 1e-3 | .560 | 3.73 | .543 |
| 1e-4 | .561 | 31.8 | .543 |
| 1e-5 | .561 | 311. | .543 |

`"quantile"` diverges because there is no derivative to find: a step function's
difference quotient grows like one over the step. The other two settle, and
`"range"` settles closer, because an affine map has a known constant derivative
where a smoothed distribution function contributes an estimated density. That
density is off by a median of 6% and by as much as 40% in the sparse upper tail
of the fit measured, and it is undersmoothed for the purpose besides, the
bandwidth having been chosen for the distribution function rather than for its
derivative. So `"range"` stays the recommendation when a slope is the quantity
being reported, and that is now what `?bartisan-marginaleffects` says.

**One test had to change, and not because anything broke.** `test-control.R`'s
ramp test read `|log(late/early)|` off a single chain against a threshold of
`log(1.5)`. Across eight seeds that ratio runs from .02 to .76 under either
transform, exceeding the threshold in two of eight both ways, and `smoothcdf`'s
median is lower than `quantile`'s (.12 against .30). The threshold was tighter
than the noise and had been passing on its seed by luck. It now averages four
seeds, which the noise allows and which still leaves the factor of four the bug
produced far above the bar.
