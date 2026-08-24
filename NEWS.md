# genbart 0.0.0.9000

* Initial version. `genbart()` fits Bayesian additive regression trees for
  response distributions outside the conditionally conjugate case, using the
  Laplace-approximation reversible-jump sampler of Linero (2025), with optional
  soft decision rules following Linero and Yang (2018).

* Sampler and prior settings may be passed to `genbart()` directly, as in
  `glm()`: any argument of `genbart_control()` given in `...` is merged into
  `control`, overriding a value set there. A name that is not an argument of
  `genbart_control()` is an error rather than being silently ignored.

* `ordinal("probit")` is 14 times faster with soft rules and 30 times faster
  with hard ones, at the same accuracy, through the latent normal of Albert and
  Chib (1993). Conditional on it the leaf target is quadratic, so the sampler
  takes the closed form instead of Fisher scoring plus a Metropolis ratio, and
  the two cumulative-normal evaluations per observation per pass disappear. On
  by default via `genbart_control(augment = )`. The cutpoints are drawn from the
  ordinal likelihood with the latent variables integrated out and the latent
  variables redrawn immediately afterwards -- a partially collapsed Gibbs
  sampler -- because their conditional given the latent variables mixes worse the
  larger the sample. `ordinal("logit")` has no such representation and is
  unchanged.

* `genbart_control(gate = "smoothstep")` offers a soft decision rule whose gate
  has bounded support: past `4.06 * bandwidth` from the cutpoint an observation
  takes one side outright, so the other subtree is never visited and the gate
  costs no `exp()`. Measured 1.4 times faster than the logistic gate with no loss
  of accuracy. The constant makes `bandwidth` mean the same amount of smoothing
  under either gate.

* `genbart_control(bandwidth_every = )` draws each tree's bandwidth every *k*-th
  sweep rather than every sweep. The bandwidth move rebuilds every membership
  weight in a tree and is about half the cost of a soft-rule fit, but drawing it
  is what lets soft rules sharpen towards hard ones, so this is a trade and the
  default is unchanged.

* `ordinal("logit")` is now augmented too, and is 15 times faster with hard
  rules and 7 times faster with soft ones, at unchanged accuracy. The cumulative
  logit has no latent normal, but a logistic variate is a normal whose precision
  is random, and the mixing conditional turns out to be an exact Polya-Gamma(2,
  |r|) draw rather than the Kolmogorov-Smirnov distribution the literature
  usually reaches for -- so it reuses the Devroye sampler already in the package
  and nothing approximate enters. Both ordinal links are now covered by
  `augment = TRUE`.

* `genbart_control(gate = "smootherstep")` adds a twice-differentiable bounded
  gate (the Beta(3, 3) CDF) alongside `"smoothstep"`. It is the same speed, so
  the choice between the two is about smoothness. Measurement corrected an
  assumption in the process: what a bounded gate saves is the `exp()`, not the
  work on the far side of the cutpoint -- at a bandwidth wide enough to truncate
  nothing it is still 1.45 times faster than the logistic gate.

* A hard tree no longer stores membership weights. Every one of them is exactly
  one, so the loops that dominate the sampler were loading a vector of ones and
  multiplying by it. Draws are bit-identical, since multiplying by 1.0 is exact.

* The leaf-level sums no longer materialize a node's predictors into a buffer
  before reading them back; the family reads them from the predictor matrix
  directly. Also bit-identical. Worth about 4% -- less than expected, because
  the buffer was small enough to stay in cache.

* Two prediction types for ordinal models, following
  \pkg{WeightIt}. `predict(type = "mean")` reports the probabilities weighted by
  the category labels read as numbers, so a response with levels `"1"`, `"2"`,
  `"4"` gets a mean between one and four; `values` says what the categories are
  worth when the labels are not numbers. `predict(type = "stdlv")` reports the
  additive predictor divided by the standard deviation of the latent variable it
  indexes, which is what puts fits with different links, or different amounts of
  signal, on one scale.

* **Group-level random intercepts**, written in \pkg{lme4}'s notation:
  `y ~ x1 + x2 + (1 | school)`. Several grouping factors are allowed and
  `(1 | a/b)` expands to nesting. Each additive predictor gets its own set, so a
  zero-inflated count model has a group effect on the count part and another on
  the inflation part, each with its own standard deviation drawn under the same
  half-Cauchy prior the leaf scale uses. `fit$ranef` and `fit$tau` hold the draws;
  `summary()` reports the scales.

  A random *slope* is refused rather than ignored, and the message says what to do
  instead: put the variable in the fixed part, where a tree can split on the group
  and the variable together and get an interaction of any shape. The restriction
  is structural — a random intercept is a scalar entering the predictor with
  weight one over a set of observations, which is what a leaf is once its gate is
  removed, so the sampler's leaf machinery handles it exactly.

  Whether to use this or to put the grouping factor in as an ordinary predictor
  depends on the shape of the data, and `?genbart` says which way round: with few
  large groups the predictor route measures better, and the random intercept wins
  once there are many small groups.

* Ordinal cutpoints are now documented against [MASS::polr()] with a worked
  example, and there is a test: the chart genbart reports in is the one `polr()`
  reports in when its predictors are centered, and the two agree to Monte Carlo
  error on a linear truth.

* **`na.action` now defaults to `na.pass`**, so missing predictor values are kept
  and the splitting rules decide where they go, rather than the rows being
  dropped. A tree can do something better with a missing value than either
  imputing it or discarding the row, and that ability was previously off unless
  asked for. Rows with a missing response, weight or offset are still dropped,
  with a warning. Pass `na.action = na.omit` for the old behavior.

* Bug fix, found while making that change: the `na.action` default was never
  used at all. The `model.frame()` call is rebuilt from `match.call()`, which
  records only what the caller actually wrote, so an argument left at its default
  was absent and `model.frame()` fell back on `getOption("na.action")`. Any
  future default for this argument would have been silently ignored the same way.

* `ordinal("cloglog")` joins the logit and probit links. The complementary
  log-log model is the discrete proportional hazards model, which gives it a
  latent variable of its own: an exponential waiting time, conditional on which
  the target is the exponential form the sampler already exploits for the Poisson
  and gamma families. That is 5.1 times faster under hard rules, with mixing
  slightly better rather than worse, and is on by default.

* **Ordinal models with three or more categories now report a centered
  predictor and free cutpoints**, as `MASS::polr()` does, rather than pinning the
  first cutpoint at zero. The cutpoints are then readable as category boundaries
  without having to remember that the predictor carries an offset. This is a
  change of chart and not of model -- the shift is applied to the recorded
  predictor, cutpoints and leaf values together, so every fitted probability is
  unchanged and the stored forest still replays to the reported predictor -- but
  `fit$aux[, "cut1"]` is no longer a column of zeros. Two-category responses keep
  the old chart, which is the one that makes them identical to binary regression.

* Support for the \pkg{marginaleffects} package, so predictions, comparisons,
  slopes and hypothesis tests on any of them can be computed directly from a fit.
  Uncertainty comes from the posterior draws rather than from a delta-method
  approximation. See `?genbart-marginaleffects`, which also documents why slopes
  want `x_transform = "range"`: the default quantile transform is an empirical
  distribution function, so the fit is a step function of the original predictor
  and a derivative of it is not well behaved. The model frame is now retained in
  the fit, as `glm()` retains it, because that is what counterfactual grids are
  built from.

* Faster again: soft-rule fits are about 35% quicker and hard-rule fits about 5%.
  The largest single item was replacing `push_back` with a sized fill in the
  innermost loop that divides a node's support -- four vectors' worth of
  per-element capacity tests, worth 20% of a soft fit on its own. Beyond that,
  tree nodes are recycled rather than allocated on every birth proposal (about
  two thirds of which are rejected), the child weights a two-child target needs
  now come out of the split that produces them rather than from a second pass
  over the same gates, and the change move no longer allocates two vectors per
  proposal. All of those are bit-identical.

* The bandwidth move got two exact improvements. A tree with no splits has no
  gate, so the bandwidth does not enter the likelihood and its full conditional
  is exactly the prior; it is now drawn from the prior directly rather than by a
  Metropolis step that rebuilt every membership weight in the tree to decide a
  move that cannot change the fit. And a rejected proposal -- 58% of them -- is
  rolled back from a snapshot rather than by evaluating every gate a second
  time, which took the rollback from 17% of a soft fit to under 1%.

* Bug fix: a quantity the sampler holds fixed -- an ordinal model's first
  cutpoint, which is pinned at zero for identifiability -- reported an R-hat of
  `-Inf` with a warning and a fabricated effective sample size of about 6, rather
  than `NA`. The autocovariance of a constant is a rounding error rather than
  exactly zero, so it passed the variance guard.

* Faster throughout: the bandwidth move no longer evaluates the whole likelihood
  twice for two numbers it subtracts, the whole-sample log density is statically
  dispatched rather than one virtual call per observation, and the birth, death
  and change moves no longer allocate a vector of node pointers each.

* Supported families: Gaussian, binomial (logit, probit, complementary
  log-log), Poisson, negative binomial, gamma, ordinal, multinomial,
  zero-inflated Poisson and negative binomial, ordered beta, accelerated
  failure time models for right-censored data (Weibull, log-logistic,
  log-normal), and Gaussian location-scale regression.

* `multinomial()` now fits one forest per category and leaves the model
  unidentified, following Murray (2021), rather than singling out a reference
  category whose choice changes the fit. The leaf prior scale is divided by
  `sqrt(2)` so the prior on the identified log odds is unchanged. Pass
  `reference` for the old parameterization.

* Links the package does not compile are accepted for `gaussian()`,
  `binomial()`, `poisson()`, `negbin()` and `Gamma()`, and applied from R by
  composing them onto the scale the compiled family works on. So
  `binomial("cauchit")` works, as does any `link-glm` object, including one
  written by hand.

* `custom_family()` takes a log density as an R function and fits the model that
  goes with it, with any number of additive predictors and optional analytic
  derivatives. It cannot draw a nuisance parameter.

* **`augment` now defaults to `TRUE`.** It rewrites a likelihood as the margin of
  a Gaussian or a Poisson one where that has been measured to pay: the binomial
  family always, and the negative binomial under hard rules. The target is the
  same either way -- only the sampler changes -- but it does change, so a fit
  from an earlier version will not reproduce draw for draw. `augment = FALSE`
  restores the old behavior.

* **The exponential form.** Hill et al. (2020, sec. 3.1.5) observe that the count
  models share the shape `c + a*eta + b*exp(s*eta)`. Like the quadratic form,
  three numbers from one pass over a node determine it everywhere, so the Laplace
  fit is iterated on scalars rather than on the data. It applies to the Poisson
  and gamma families, and under hard rules only, since a soft rule gives each
  observation its own exponent. Worth 1.5x on a Poisson fit and 1.45x on a gamma
  one.

* **The negative binomial is now written as a Poisson-gamma mixture** rather than
  through Polya-Gamma augmentation, which puts it in the exponential form above
  and costs one gamma draw per observation instead of a two-hundred-term series.
  Measured in effective samples per second under hard rules that turns a 0.5x
  loss into a **2.0x gain**, and it improves the mixing rather than degrading it.
  The Polya-Gamma sampler for real-valued parameters is 10x cheaper as well --
  20 terms with a tail matched on two moments rather than 200 with one -- but it
  is no longer on any default path.

* **Leaf sums are accumulated against each family's concrete type**, so the
  compiler can inline the arithmetic instead of dispatching four virtual calls
  per observation. Worth 1.05x to 1.20x; less than it would have been before the
  conjugate shortcut, because the family evaluations are no longer the bulk of
  the work. `generic_accumulate = TRUE` forces the dispatched version, and the
  two agree bitwise.

* **Rank-normalized folded split R-hat and the bulk and tail effective sample
  sizes** (Vehtari et al. 2021) replace the classic split R-hat. `fit$rhat` is now
  a data frame with a row per quantity. The tail column is the one that governs
  an interval endpoint, and it is reported separately because a run can be
  adequate for a posterior mean and nowhere near adequate for a 2.5% quantile.
  Validated against the `posterior` package to within 2% and against the closed
  form for an AR(1) chain.

* **Speed.** On a binary probit fit to the 614-row `lalonde` data with 50 trees
  and 1000 warmup plus 1000 saved draws, the default configuration went from
  82.3 seconds to 4.8, and hard rules from 45.2 to 1.5. Three changes, in order
  of size:

  - `augment` in `genbart_control()` rewrites a likelihood as the margin of a
    Gaussian one, which makes the target quadratic in the additive predictor.
    `TRUE` covers the binomial family: Albert and Chib (1993) for a probit link
    and Polya-Gamma augmentation (Polson, Scott and Windle 2013) for a logit
    one. Measured in effective samples per second, which is the ratio that
    matters since every augmentation trades speed for mixing: probit 3.8x,
    logit 2.1x. `"multinomial"` is 1.6x and available by name, because its
    mixing penalty is severe. `"negbin"` is a measured net loss at 0.5x and is
    there only for completeness -- its Polya-Gamma parameter is not an integer,
    so each draw needs a truncated series rather than the exact method the
    others use. Polya-Gamma does not apply to the Poisson, gamma, survival,
    ordered beta or location-scale families at all.
  - The binomial family shared its score and information, computed the probit
    tails from one call to `pnorm` instead of two, and stopped evaluating the
    term that a binary response multiplies by zero. Probit alone: 3.3x.
  - Where the log density is quadratic in the predictor the Laplace
    approximation is exact, so Fisher scoring is stopped after one step and the
    leaf refresh becomes a Gibbs step with acceptance one. Gaussian fits do a
    third less work.
  - **The conjugate shortcut.** A leaf value enters the predictor linearly, so a
    quadratic log density makes the whole log target quadratic in the leaf value:
    one pass over a node determines it everywhere, and the Laplace fit and the
    log target at any value are then arithmetic rather than another pass. A birth
    move needed six passes over the node and now needs two. Verified against the
    general path to rounding error rather than to a tolerance. This is what a
    conjugate sampler does, and it applies to a Gaussian response and to
    everything `augment` rewrites into one.

* The project's `.Rprofile` sets `options(pkg.build_extra_flags = FALSE)`, so
  `devtools::load_all()` and `devtools::install()` compile with optimization
  instead of `-O0`. Measured on this package, `-O3`, `-mcpu=native` and `-flto`
  are all within one percent of the `-O2` default, so there is nothing further to
  gain from flags -- only from not being at `-O0`, which costs a factor of five
  to twenty.

* `_dev/benchmark.Rmd` times genbart against dbarts, BART and bartMachine on the
  same data, reporting effective samples per second and held-out accuracy
  alongside the clock. It refuses to run on an unoptimized build.

* `chains` runs independent chains and pools them, under whatever backend the
  caller has planned with `future::plan()` -- `multisession`, `multicore`, a
  cluster, or mirai's `mirai_multisession`. One `set.seed()` reproduces the whole
  run whatever the backend, and split-R-hat for the log likelihood, the leaf
  scales, the nuisance parameters and the additive predictor is reported in the
  `rhat` element.

* Missing predictor values are handled natively with `na.action = na.pass`,
  which keeps rows whose predictors are missing instead of dropping them and
  lets each splitting rule say where a missing value goes. This is missingness
  incorporated in attributes (Twala et al. 2008; Kapelner and Bleich 2015),
  including the rule that splits on missingness itself, so a variable whose
  absence carries the signal is usable. The choice is drawn from its prior, so
  no acceptance ratio changes and complete data behaves exactly as before.

* `predict()` supports `type = "density"`, the conditional density of the
  outcome given the predictors evaluated at the outcome's own value. Summed on
  the log scale this gives a held-out log score.

* The bandwidth of a soft decision rule is tuned during warmup towards the
  acceptance rate that is optimal for a one-dimensional random walk, and frozen
  before the retained draws, which roughly halves its autocorrelation.
