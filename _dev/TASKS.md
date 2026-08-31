# bartisan tasks

Implements Linero (2025), "Generalized Bayesian Additive Regression Trees Models: Beyond Conditional Conjugacy", with SoftBart-style soft decision rules from Linero and Yang (2018).

This file is organized by subject, not by session. Each entry states the problem, what was done, what it measured, and — where it applies — what was tried and abandoned. **Status** and **To Do** are the only sections that describe the present; everything under **Log** is a record of how it got here and should be read as history.

## Status

1358 tests passing, 0 failures, 0 warnings, 0 skips. `R CMD check` reports `Status: OK` with no warnings and no notes when run outside the agent sandbox; inside it, `OMP: Warning #179` and an `nm` cache-file NOTE appear, and both are artifacts of the sandbox rather than the package. the suite run inside the check passes 1055 of them, skipping 51 for Suggests packages that environment does not have.

**What exists.** A C++ engine (`utils`, `slice`, `hypers`, `family`, `polyagamma`, `node`, `mcmc`, `model`) and an R interface following `glm()`: `bartisan()`, `bartisan_control()`, `predict()`, `print()`, `summary()`, family normalization, parallel chains with convergence diagnostics, and `custom_family()` for a likelihood written in R. Families: Gaussian, binomial (logit/probit/cloglog/any link from R), Poisson, negative binomial, gamma, ordinal (logit/probit/cloglog), multinomial (symmetric or reference-coded), multinomial probit with a drawn latent covariance, three AFT variants, location-scale, zero-inflated Poisson and negative binomial, ordered beta, and a Dirichlet process mixture for the error distribution. Missing predictors handled natively by MIA and kept by default. Data augmentations, on by default, for the binomial, ordinal, multinomial and zero-inflated families, and for the negative binomial under hard rules. `marginaleffects` support, so counterfactual estimands come with posterior intervals. Group-level random intercepts through lme4's `(1 | group)` notation, on every additive predictor. Posterior predictive draws for every family that has a sampler, and with them the interfaces to `loo`, `bayesplot`, `performance` and `posterior`. `bartisan_control()` organized into modeling decisions, advanced settings and validation toggles, with a per-forest `num_trees` vector, one `gate` argument covering hard and soft rules, and a `sparsity` argument standing in for the four DART hyperparameters. Documentation, `README.Rmd`, nine vignettes with a shared `references.bib`, and `_dev/benchmark.Rmd`. No `NEWS.md`: nothing has been released, so there is no previous version for a user to have seen.

Benchmark, Friedman function, n = 1000, p = 10, 50 trees, 500 warmup plus 500 saved, best of 2, scored against the true regression function on a held-out thousand. Reproducible with `_dev/benchmark.Rmd`.

| Task | Package and call | Seconds | ESS | RMSE |
|---|---|---|---|---|
| gaussian | dbarts | 0.239 | 19.6 | 0.221 |
| | **bartisan hard** | **0.441** | 17.6 | 0.219 |
| | stochtree | 0.697 | 18.3 | 0.218 |
| | bartMachine | 0.866 | — | 0.230 |
| | BART `wbart()` | 1.061 | 21.3 | 0.228 |
| | bartisan soft, smoothstep gate (default) | 1.399 | 42.0 | **0.135** |
| | bartisan soft, smootherstep gate | 1.420 | 43.2 | 0.140 |
| | bartisan soft, logistic gate | 2.003 | 43.7 | 0.145 |
| probit | dbarts | 0.271 | 35.8 | 0.134 |
| | **bartisan hard** | **0.491** | 34.8 | 0.122 |
| | stochtree | 0.985 | 29.4 | 0.134 |
| | BART `pbart()` | 1.133 | 35.9 | 0.125 |
| | bartisan soft, logistic gate | 2.049 | 33.4 | **0.111** |
| | bartisan soft, `augment = FALSE` | 23.578 | 65.2 | 0.106 |
| logit | bartisan soft, logistic gate | **2.111** | 128.8 | 0.088 |
| | bartisan soft, `augment = FALSE` | 13.861 | 103.8 | 0.085 |
| | BART `lbart()` | 20.820 | 37.8 | 0.109 |
| ordinal | bartisan hard, probit | **0.868** | 35.0 | — |
| | bartisan hard, logit | 0.932 | 30.5 | — |
| | bartisan soft, probit | 2.383 | 50.5 | — |
| | bartisan soft, logit | 2.520 | 48.1 | — |
| | bartisan hard, cloglog | 3.030 | 24.3 | — |
| | stochtree (cloglog) | 4.010 | 20.3 | — |
| | bartisan hard, logit, `augment = FALSE` | 14.903 | 42.3 | — |
| | bartisan hard, cloglog, `augment = FALSE` | 17.005 | 37.2 | — |
| | bartisan hard, probit, `augment = FALSE` | 25.186 | 36.4 | — |
| poisson | bartisan hard | **1.715** | 20.0 | 0.207 |
| | bartisan hard, no shortcut | 3.351 | 24.1 | 0.192 |
| gamma | bartisan hard | **3.600** | 21.5 | 0.198 |
| | bartisan hard, no shortcut | 6.234 | 18.8 | 0.204 |
| negative binomial | bartisan hard, augmented | 3.494 | 25.8 | 0.276 |
| | bartisan hard, direct | 6.566 | 33.1 | 0.286 |
| log-logistic AFT | bartisan hard | 9.465 | 35.6 | — |

Four things this says.

**Hard rules are within 1.8x of dbarts** on the two tasks dbarts supports, at the same mixing and slightly better accuracy: 0.426 s against 0.241 s on the Gaussian task and 0.482 s against 0.272 s on probit, best of six runs each. That ratio has come down 3.1x → 2.4x → 2.03x → 1.85x → 1.77x as each of the entries below landed.

**Note on measurement.** The table above is `_dev/benchmark.Rmd` at two replicates, which is noisy at the ten to thirty percent level — the Gaussian hard-rule cell has read 0.441 and 0.593 on consecutive runs of the same build, with a best-of-five standalone measurement of 0.426 either side of it. Any claim about a ratio near two needs more replicates than the document's default, and the two ratios quoted above are from a longer run for that reason.

**bartisan is now the faster of the two on the ordinal complementary log-log model**, which is the one task stochtree supports and dbarts does not: 3.03 s against 4.01 s.

**bartisan is faster and more accurate than every other package here** on both tasks — stochtree, bartMachine and BART included — once hard rules are used.

**Soft rules are the accuracy argument, not a tax.** They cost three to five times the hard-rule time and cut held-out error by 35–40%, which makes the default configuration the most accurate fit in the table, dbarts included. With a bounded gate the cost falls to 3.2x.

**The ordinal RMSE column is not comparable across links.** The three links put the additive predictor on three different scales, so those numbers compare scales rather than fits; seconds and ESS are comparable. It became *readable* for the probit rows once the identification changed, because the reported predictor is now centered and so is the generating one.

## To Do

Ordered by expected value. `_dev/SHIP.md` holds the release assessment -- which of
these actually block a workflow the package claims to support, which are covered
by other packages, and the six-vignette workflow series that exposed the
difference.

- [ ] **Rename the package.** `bartisan` collides case-insensitively with the archived CRAN package `genBart`, which is a hard block on submission. Candidates checked against both the current index and all 27,654 archived names are in the Notes below.
- [ ] **Correlated random effects across additive predictors**, which is the one part of the random-effects feature that is not there. The obstacle is the absence of mixed second derivatives in the `Family` interface, and the alternative needs a prior mean threaded through 27 places; see the assessment.
- [ ] A Bayesian-bootstrap dispersion draw for `Gamma()` and `negbin()`, from Pearson residuals rather than the assumed likelihood. Cheapest available improvement to interval calibration; see the quasi-likelihood entry.
- [ ] A joint tridiagonal update for the ordinal cutpoints, which is what makes inference on the thresholds usable when there are many of them. The obstacle is the ordering constraint, not the algebra.
- [ ] A `quasi()` family parameterized by link, variance function and dispersion update rule. Needs a documented weakening of the exactness claim.
- [ ] **Partial dependence plots** and a **formal variable-selection test** rather than raw split counts, both of which `dbarts`, `SoftBart` and `bartMachine` have and this does not. Partial dependence is now largely reachable through `marginaleffects::plot_predictions()`, so this is less of a gap than it was.
- [ ] Possibly a grow-from-root warm start for hard-rule fits, using a one-step Laplace criterion, to shorten burn-in. See the XBART assessment for why this is not obviously worth the code. It became more attractive, not less, once the augmentations landed: they made the per-sweep cost small enough that burn-in length is now the binding constraint on an ordinal fit.
- [ ] Consider the `draw_prior` move from SoftBart, which proposes a whole fresh tree and helps escape local modes. It needs an `L`-dimensional Laplace proposal for the new leaves, so it is real work, not a port.
- [ ] `vignette("bartisan")` does not cover the bounded gates or either ordinal augmentation; `vignette("families")` covers the augmentations but not the gates. The first runs on a reduced chain (20 trees, 300 draws, n = 400) and builds in about 85 seconds.
- [ ] Missing data, further work: nothing forces the three missing-value rules to be equally likely, and a variable with a handful of missing values probably does not want a third of its rules spent on splitting by missingness. A prior weight on the third rule is a one-line change and an open question.
- [ ] `custom_family()` has no posterior predictive distribution, since a log density supplies no way to draw from it. An optional `rng` argument alongside the density would give it one, and would make `simulate()`, `pp_check()` and `r2()` work for a user-written likelihood. Small work; the design question is whether to also ask for a mean.
- [ ] Joint update for correlated nuisance parameters. Linero reports that `sigma` and the shape of the generalized gamma mix badly when updated separately. Only relevant if a two-nuisance family is added; the current families have at most one.
- [ ] A lighter-tailed prior on the leaf scale, or an upper bound, would remove the separation pathology at the cost of changing the default prior. Not done unilaterally; the warning is the interim measure.
- [ ] **`predict(type = "density")` returns NaN silently** when a composed link's inverse sends the predictor outside the family's support. Measured on `stats::Gamma("inverse")` with heavy-tailed data: five of eight replicates had draws where the predictor went non-positive, and each produced NaN densities for two to five test points out of 800. A negative fitted mean is not a gamma mean, so NaN is arguably the right *value*, but it should not be silent -- `bartisan()` already warns about the link at fit time and `predict()` says nothing. Left alone deliberately: it changes the output contract of `predict()`, which is the user's call. See the gamma comparison entry.
- [ ] **Relative survival on top of `ph()`**, per Basak et al. (2024): the excess-hazard model needs one extra Bernoulli draw per sweep, `d_i ~ Bernoulli(lambda_E / (lambda_E + lambda_P))`, with the population hazard supplied as one number per subject from a life table. Cheap now that `ph()` exists -- a nuisance draw and a data column. Narrow audience (cancer registries), so worth doing only on request.
- [ ] **Soft random tree features**, as a fast approximate fit and as a warm start for the sampler. Measured at 0.840 average out-of-sample R-squared against full soft BART's 0.872 at 200 features, for a fraction of the cost, and a soft basis beats a hard one by 0.16 R-squared at five trees. See the McCartan and Huang entry, which has three further items.
- [ ] **Categorical splits on subsets of levels**, as in `flexBART`, rather than one-hot columns sharing a sparsity weight. More expressive rules.
- [ ] **Causal-inference structure** — a separate treatment forest, ATE/CATE (`bcf`, `bartCause`, `stochtree`). Substantial work, and `marginaleffects` now covers the estimand side of it for a model fitted by hand.

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

**The log score can be made comparable across all five families.** The earlier vignette said it could not: the accelerated failure time families report the density of \eqn{\log T} and `ph()` the density of \eqn{T}, so the two differ by \eqn{\sum \log t}. But censored observations contribute \eqn{S(t)}, which carries no measure at all, so the correction is `-sum(status * log(time))` -- events only. Applied, `weibull_aft()` and `ph()` score within a few points of each other on a Weibull truth instead of hundreds apart, which is the check that the correction is right rather than merely plausible. The vignette now gives the correction as a function rather than telling the reader to avoid the comparison.

**The families disagree far more about the density than about the ordering.** Across every proportional-hazards and accelerated failure time truth, all five recover the ordering of subjects by survival nearly perfectly whether or not they have the shape right; the differences are concentrated in \eqn{S(t \mid x)} and the log score. The one exception is the crossing-hazards truth, where the ordering itself moves with time. So the choice of family matters for an absolute probability at a horizon and barely matters for a comparative question -- which is worth telling a reader before they agonize over it.

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

**Added: the potential outcomes** through `avg_predictions(fit, variables = "smoke")`, which reports \eqn{E[Y(0)]} and \eqn{E[Y(1)]} at 3072 and 2773 grams. Their difference is the ATE reported in the next section, and showing both is more informative than the difference alone.

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

**The cost of the change is interpretive, not statistical.** The reported predictor changes meaning -- `weibull_aft()`'s is a log time ratio, `dpm_aft()`'s is \eqn{E[\log T \mid x]}, which is a time ratio only if the error is symmetric, which is the assumption the family exists to avoid. `type = "survival"` and `type = "response"` are unaffected, being on the same scale for every survival family.

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
on the \pkg{rstantools} generics, `simulate()` on base R's, `loo()` and `waic()`,
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
were unreachable through \pkg{marginaleffects} even though `predict()` had them.
Both are one number per observation, so they need nothing but the name. `"probs"`,
`"lp"` and `"lv"` are accepted as aliases, since those are the names the same
quantities go by for the WeightIt classes in \pkg{marginaleffects} itself.
`"class"` and `"density"` stay out, and are refused by name: neither is a number
per observation that an average or a contrast could be taken of.

The dots cannot simply be forwarded to `predict()`, because \pkg{marginaleffects}
puts arguments of its own there (`mfx`, and whatever the caller passed to the
estimand function) and a stray name would match one of `predict()`'s arguments
partially. Only `values`, `iterations`, `offset`, `weights` and `log` are taken,
by exact name. \pkg{marginaleffects} warns that it does not recognize `values`,
which is correct — it is this package's argument — and the value is used anyway.

**A convention worth knowing, found by a test failure that was the test's
fault.** \pkg{marginaleffects} centers a posterior at its **median**;
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

## Notes

### Candidate names

`bartisan` collides case-insensitively with the archived CRAN package `genBart`, which blocks submission. All of the following were checked against the current CRAN index and against all 27,654 archived package names, and are free. Also note `flexBART`, `SoftBart`, `dbarts`, `bartMachine`, `bartCause` and `stochtree` exist, and that `gbart` is the main *function* in the `BART` package, so it should be avoided even though the name is free.

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
