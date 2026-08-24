# genbart tasks

Implements Linero (2025), "Generalized Bayesian Additive Regression Trees Models: Beyond Conditional Conjugacy", with SoftBart-style soft decision rules from Linero and Yang (2018).

This file is organized by subject, not by session. Each entry states the problem, what was done, what it measured, and — where it applies — what was tried and abandoned. **Status** and **To Do** are the only sections that describe the present; everything under **Log** is a record of how it got here and should be read as history.

## Status

850 tests passing, 0 failures, 0 warnings. `R CMD check` Status: OK.

**What exists.** A C++ engine (`utils`, `slice`, `hypers`, `family`, `polyagamma`, `node`, `mcmc`, `model`) and an R interface following `glm()`: `genbart()`, `genbart_control()`, `predict()`, `print()`, `summary()`, family normalization, parallel chains with convergence diagnostics, and `custom_family()` for a likelihood written in R. Families: Gaussian, binomial (logit/probit/cloglog/any link from R), Poisson, negative binomial, gamma, ordinal (logit/probit/cloglog), multinomial (symmetric or reference-coded), three AFT variants, location-scale, zero-inflated Poisson and negative binomial, ordered beta. Missing predictors handled natively by MIA and kept by default. `marginaleffects` support, so counterfactual estimands come with posterior intervals. Group-level random intercepts through lme4's `(1 | group)` notation, on every additive predictor. Documentation, `README.Rmd`, `NEWS.md`, a vignette, and `_dev/benchmark.Rmd`.

Benchmark, Friedman function, n = 1000, p = 10, 50 trees, 500 warmup plus 500 saved, best of 2, scored against the true regression function on a held-out thousand. Reproducible with `_dev/benchmark.Rmd`.

| Task | Package and call | Seconds | ESS | RMSE |
|---|---|---|---|---|
| gaussian | dbarts | 0.239 | 19.6 | 0.221 |
| | **genbart hard** | **0.441** | 17.6 | 0.219 |
| | stochtree | 0.697 | 18.3 | 0.218 |
| | bartMachine | 0.866 | — | 0.230 |
| | BART `wbart()` | 1.061 | 21.3 | 0.228 |
| | genbart soft, smoothstep gate | 1.399 | 42.0 | **0.135** |
| | genbart soft, smootherstep gate | 1.420 | 43.2 | 0.140 |
| | genbart soft (default) | 2.003 | 43.7 | 0.145 |
| probit | dbarts | 0.271 | 35.8 | 0.134 |
| | **genbart hard** | **0.491** | 34.8 | 0.122 |
| | stochtree | 0.985 | 29.4 | 0.134 |
| | BART `pbart()` | 1.133 | 35.9 | 0.125 |
| | genbart soft (default) | 2.049 | 33.4 | **0.111** |
| | genbart soft, `augment = FALSE` | 23.578 | 65.2 | 0.106 |
| logit | genbart soft (default) | **2.111** | 128.8 | 0.088 |
| | genbart soft, `augment = FALSE` | 13.861 | 103.8 | 0.085 |
| | BART `lbart()` | 20.820 | 37.8 | 0.109 |
| ordinal | genbart hard, probit | **0.868** | 35.0 | — |
| | genbart hard, logit | 0.932 | 30.5 | — |
| | genbart soft, probit | 2.383 | 50.5 | — |
| | genbart soft, logit | 2.520 | 48.1 | — |
| | genbart hard, cloglog | 3.030 | 24.3 | — |
| | stochtree (cloglog) | 4.010 | 20.3 | — |
| | genbart hard, logit, `augment = FALSE` | 14.903 | 42.3 | — |
| | genbart hard, cloglog, `augment = FALSE` | 17.005 | 37.2 | — |
| | genbart hard, probit, `augment = FALSE` | 25.186 | 36.4 | — |
| poisson | genbart hard | **1.715** | 20.0 | 0.207 |
| | genbart hard, no shortcut | 3.351 | 24.1 | 0.192 |
| gamma | genbart hard | **3.600** | 21.5 | 0.198 |
| | genbart hard, no shortcut | 6.234 | 18.8 | 0.204 |
| negative binomial | genbart hard, augmented | 3.494 | 25.8 | 0.276 |
| | genbart hard, direct | 6.566 | 33.1 | 0.286 |
| log-logistic AFT | genbart hard | 9.465 | 35.6 | — |

Four things this says.

**Hard rules are within 1.8x of dbarts** on the two tasks dbarts supports, at the same mixing and slightly better accuracy: 0.426 s against 0.241 s on the Gaussian task and 0.482 s against 0.272 s on probit, best of six runs each. That ratio has come down 3.1x → 2.4x → 2.03x → 1.85x → 1.77x as each of the entries below landed.

**Note on measurement.** The table above is `_dev/benchmark.Rmd` at two replicates, which is noisy at the ten to thirty percent level — the Gaussian hard-rule cell has read 0.441 and 0.593 on consecutive runs of the same build, with a best-of-five standalone measurement of 0.426 either side of it. Any claim about a ratio near two needs more replicates than the document's default, and the two ratios quoted above are from a longer run for that reason.

**genbart is now the faster of the two on the ordinal complementary log-log model**, which is the one task stochtree supports and dbarts does not: 3.03 s against 4.01 s.

**genbart is faster and more accurate than every other package here** on both tasks — stochtree, bartMachine and BART included — once hard rules are used.

**Soft rules are the accuracy argument, not a tax.** They cost three to five times the hard-rule time and cut held-out error by 35–40%, which makes the default configuration the most accurate fit in the table, dbarts included. With a bounded gate the cost falls to 3.2x.

**The ordinal RMSE column is not comparable across links.** The three links put the additive predictor on three different scales, so those numbers compare scales rather than fits; seconds and ESS are comparable. It became *readable* for the probit rows once the identification changed, because the reported predictor is now centered and so is the generating one.

## To Do

Ordered by expected value.

- [ ] **Rename the package.** `genbart` collides case-insensitively with the archived CRAN package `genBart`, which is a hard block on submission. Candidates checked against both the current index and all 27,654 archived names are in the Notes below.
- [ ] **Correlated random effects across additive predictors**, which is the one part of the random-effects feature that is not there. The obstacle is the absence of mixed second derivatives in the `Family` interface, and the alternative needs a prior mean threaded through 27 places; see the assessment.
- [ ] **Posterior predictive draws of a new outcome.** `predict()` returns the mean and the link but never draws Y, so predictive intervals for a new observation are unavailable. Needs a sampler per family; small work, high value, and the largest functional gap against other packages.
- [ ] The exponential form for **soft** rules. A soft rule gives each observation its own exponent, so three numbers no longer suffice — but a small fixed number of them might, since the membership weights take few distinct values in practice. The payoff would be on the default configuration rather than only on hard rules.
- [ ] An exponential-form route for the **zero-inflated** families and the **multinomial**, which Murray reaches through a further gamma augmentation. A larger change than the negative binomial's was.
- [ ] A Bayesian-bootstrap dispersion draw for `Gamma()` and `negbin()`, from Pearson residuals rather than the assumed likelihood. Cheapest available improvement to interval calibration; see the quasi-likelihood entry.
- [ ] A joint tridiagonal update for the ordinal cutpoints, which is what makes inference on the thresholds usable when there are many of them. The obstacle is the ordering constraint, not the algebra.
- [ ] A `quasi()` family parameterized by link, variance function and dispersion update rule. Needs a documented weakening of the exactness claim.
- [ ] **Partial dependence plots** and a **formal variable-selection test** rather than raw split counts, both of which `dbarts`, `SoftBart` and `bartMachine` have and this does not. Partial dependence is now largely reachable through `marginaleffects::plot_predictions()`, so this is less of a gap than it was.
- [ ] Possibly a grow-from-root warm start for hard-rule fits, using a one-step Laplace criterion, to shorten burn-in. See the XBART assessment for why this is not obviously worth the code. It became more attractive, not less, once the augmentations landed: they made the per-sweep cost small enough that burn-in length is now the binding constraint on an ordinal fit.
- [ ] Consider the `draw_prior` move from SoftBart, which proposes a whole fresh tree and helps escape local modes. It needs an `L`-dimensional Laplace proposal for the new leaves, so it is real work, not a port.
- [ ] The vignette does not cover the bounded gates or either ordinal augmentation. It runs on a reduced chain (20 trees, 300 draws, n = 400) and builds in about 85 seconds.
- [ ] Missing data, further work: nothing forces the three missing-value rules to be equally likely, and a variable with a handful of missing values probably does not want a third of its rules spent on splitting by missingness. A prior weight on the third rule is a one-line change and an open question.
- [ ] `custom_family()` cannot draw a nuisance parameter. The clean way is a user-supplied slice-sampling target, since `slice.cpp` already has the sampler; the awkward part is that the parameter has to reach the log density, so the closure can no longer be opaque.
- [ ] Joint update for correlated nuisance parameters. Linero reports that `sigma` and the shape of the generalized gamma mix badly when updated separately. Only relevant if a two-nuisance family is added; the current families have at most one.
- [ ] A lighter-tailed prior on the leaf scale, or an upper bound, would remove the separation pathology at the cost of changing the default prior. Not done unilaterally; the warning is the interim measure.
- [ ] **Categorical splits on subsets of levels**, as in `flexBART`, rather than one-hot columns sharing a sparsity weight. More expressive rules.
- [ ] **Causal-inference structure** — a separate treatment forest, ATE/CATE (`bcf`, `bartCause`, `stochtree`). Substantial work, and `marginaleffects` now covers the estimand side of it for a model fitted by hand.

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

It is now *detected* rather than remembered. `__OPTIMIZE__` is defined by the compiler when it is optimizing, so `genbart:::.genbart_optimized()` reports the truth about the loaded library, `genbart()` warns once per session when the answer is no, and `_dev/benchmark.Rmd` refuses to run at all. The project `.Rprofile` sets `options(pkg.build_extra_flags = FALSE)` so `load_all()` produces an optimized library in the first place.

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

`genbart()` now warns when the posterior mean of the leaf scale settles more than five times above its prior median. The threshold has room: a genuinely strong signal needs about twice the default scale, and five times corresponds to an ensemble prior standard deviation of 7.5 on the log-odds scale, which is essentially never a real signal.

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
- Storage grows as `num_save * K`. At n = K = 3200 with 1000 draws the cutpoint matrix alone is 26 MB.

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
  | genbart | −0.479, 0.904, 2.444 |

  So **genbart's chart is `polr()`'s chart with the predictors centered**, and
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
for the cutpoints, since genbart centers the predictor and WeightIt drops the
intercept column. A standardized quantity is used for differences, which that
constant leaves alone.

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
   of class genbart are not supported", which is what it did at first.

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
is documented in `?genbart-marginaleffects` rather than fixed, because the
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

| groups | per group | genbart, group as a factor | group ignored | `rbart_vi` |
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
- "The Gaussian hard-rule fit regressed by 35%." It had not: two consecutive benchmark runs of the same build read 0.441 s and 0.593 s, and a best-of-five standalone measurement read 0.426 s both times. `_dev/benchmark.Rmd` defaults to two replicates, which is not enough to support a claim about a factor near two.

## Notes

### Candidate names

`genbart` collides case-insensitively with the archived CRAN package `genBart`, which blocks submission. All of the following were checked against the current CRAN index and against all 27,654 archived package names, and are free. Also note `flexBART`, `SoftBart`, `dbarts`, `bartMachine`, `bartCause` and `stochtree` exist, and that `gbart` is the main *function* in the `BART` package, so it should be avoided even though the name is free.

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
- **`na.action` defaults to `na.pass`** now, so a test that expects rows to be dropped has to ask for `na.omit` explicitly. And note what the fix to that default exposed: `model.frame()` is called through a call rebuilt from `match.call()`, so *any* argument of `genbart()` that is forwarded to `model.frame()` and left at its default is absent from that call and picks up `model.frame()`'s default instead. Adding a default to `subset`, `weights` or `offset` would be swallowed the same way.
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
