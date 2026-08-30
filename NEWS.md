# bartisan 0.0.0.9000

* **The package is now called bartisan.** It was `genbart`, which collided
  case-insensitively with the archived CRAN package `genBart` and could not have
  been submitted under that name. The fitting function is `bartisan()`, the
  settings are `bartisan_control()`, and a fitted model has class `bartisan_fit`.
  Every documentation topic, the C++ namespace and the repository directory
  follow. There is no deprecation path, since nothing was released as `genbart`.

* **Fixed: the package failed to link under a debug build.** Loading the compiled
  code after `pkgbuild::clean_dll()`, which is what RStudio's Install button and
  `pkgload::load_all()` do, failed with
  `symbol not found in flat namespace '__ZN7bartisan27OrdinalLogitAugmentedFamily9OMEGA_MINE'`.

  Three things combined. `OMEGA_MIN` is a `static constexpr` member, which is
  implicitly inline only from C++17; the package declared no `CXX_STD`, so the
  compiler default applied, and that is C++14 on the current Apple toolchain;
  and `std::max()` takes its arguments by const reference, so passing the member
  to it is an ODR-use that needs a definition the older standard does not
  generate. `R CMD INSTALL` never showed it because `-O2` folds the constant away
  and emits no reference. `pkgbuild` compiles with `-O0`, which does not.

  Fixed twice over: `src/Makevars` now declares `CXX_STD = CXX17`, which the code
  has always required, and the two call sites compare rather than call
  `std::max()`, so they link under any standard. `std::max(a, b)` is defined as
  `a < b ? b : a`, so the results are unchanged.

* **Fixed: an unused variable warning** in `src/model.cpp`, left behind by an
  earlier refactor. A clean build with `-Wall -pedantic` is now silent, and
  `R CMD check` reports `Status: OK` with no warnings and no notes.

* **Five more vignettes**, completing the workflow series. Each takes one step of
  a regression analysis and does it with BART, using the same birth weight data
  throughout so the reader is never learning a new dataset at the same time as a
  new idea.

  `vignette("effects")` replaces the coefficient table: average comparisons,
  effects for subgroups, and testing whether they differ, which is how an
  interaction is tested in a model that never had an interaction term.
  `vignette("importance")` covers predictor usage, calibrating it with added
  noise variables, and the correlated-predictor failure, where the forest
  attributes everything to a near-copy of the variable that matters and
  `prop_used` says to keep the wrong one. `vignette("diagnostics")` separates
  convergence from fit and explains which row of `fit$rhat` to act on.
  `vignette("comparison")` covers `loo()`, what a difference in `elpd` has to
  exceed to mean anything, and the Jacobian trap when outcomes are on different
  scales. `vignette("causal")` covers the potential outcomes, the ATE and the ATT,
  moderation, positivity, and regularization-induced confounding, including what
  to do about it without a separate treatment forest. The propensity score is
  fitted with `bartisan()` rather than `glm()` and its overlap is plotted with
  *cobalt*, now in Suggests.

* **`vignette("bartisan")` is now the theory document.** It had grown into a
  family-by-family tour that `vignette("families")` covers and a first-model
  walkthrough that `vignette("workflow")` now does better, and it stated the
  model nowhere. It now gives the sum-of-trees model and its generalization to an
  arbitrary density, the three parts of the prior, soft rules and the sparsity
  prior, and then how the thing is fitted: backfitting, the reversible-jump tree
  moves, the Laplace approximation that removes the conjugacy requirement, the
  three shapes of leaf-level target, data augmentation, nuisance parameters and
  the Dirichlet process mixture. Missing values, several chains, choosing the
  settings and where the time goes are kept.

* **Corrected in `vignette("workflow")`:** the claim that `sparsity = TRUE` has
  to be set for `variable_importance()` to separate predictors. It is the
  default, and always was. Measured on the Friedman function with five real and
  five noise predictors, the noise sits at `prop_used` around 0.09 with the
  default prior and climbs to 0.95 at 50 trees without it, so the prior is doing
  its job at every tree count. What the `birthwt` table in that vignette shows is
  a sample too small to rule anything out, which the text now says.

* **New vignette, `vignette("workflow")`**, which is the "Getting started"
  document. It runs one complete analysis of the `MASS::birthwt` birth weight
  data end to end: fitting, convergence and fit checks, predictor usage, effects
  through *marginaleffects*, a fitted curve, prediction for a new
  observation with both kinds of interval, and model comparison through
  *loo*. Every step is a single call at the defaults, which is the point it
  is making. Each section ends with a pointer to the vignette that goes deeper.

* **Fixed: `predictions()` returned the wrong number of rows** for any `newdata`
  that was not the full training data, and every estimand built on it was
  affected. *marginaleffects* prepends rows of its own, marked
  `rowid = -1`, and drops them again by that marker once the predictions come
  back; `get_predict.bartisan()` regenerated the column as `seq_len(n)`, which
  destroyed the marker and let those rows through.
  `predictions(fit, newdata = d[1:3, ])` returned five rows, two of them
  *marginaleffects*' own, and `avg_comparisons()` estimates moved by around
  10% once the contamination was removed. The column is now carried through when
  it is present. A regression test pins the row count for several sizes of
  `newdata`, for `datagrid()`, and against `predict()`.

* **New function `variable_importance()`.** Returns the predictor-usage table
  that `summary()` prints, as a data frame: mean splitting rules per draw with an
  interval, and `prop_used`, the proportion of draws in which the predictor got
  any rule at all. `prop_used` is the one to read first, and with
  `sparsity = TRUE` it separates signal from noise sharply enough to be used as a
  selection rule. The documentation is explicit that this is not an effect size,
  is not stable under correlated predictors, and is not causal -- three things
  the name "variable importance" invites people to assume.

* **`as_draws()` now carries the additive predictor.** It previously returned
  only the scalars -- `loglik`, the leaf scales and the nuisance parameters --
  which left out the quantity whose convergence actually matters. `fit$rhat`
  already reported `eta`, so the diagnostics knew it mattered and then did not
  hand it to *bayesplot*. A representative spread of ten `eta` columns is
  included by default, chosen across the range of the fitted function rather than
  at random; `eta = FALSE` restores the old behavior and `eta = c(1, 5)` takes
  named observations.

* **`predict(type = "density")` warns instead of returning `NaN` silently** when
  a draw implies a parameter outside the family's support. The value is still
  `NaN`, deliberately: zero would assert that the outcome is impossible, which is
  a stronger claim than the model makes and turns into `-Inf` in a log score,
  where it reads as a legitimately terrible fit rather than an undefined one.

  The warning reports the amplification, which is the part that was surprising.
  Densities are averaged over draws before the log is taken, so **one undefined
  draw makes that whole observation `NaN`**: measured on a `poisson("identity")`
  fit extrapolated past its training range, 0.4% of draws (42 of 12,000) turned
  35% of the returned values (14 of 40) into `NaN`. The condition needs
  extrapolation to arise at all -- `bartisan()` rejects out-of-support proposals
  while sampling, so an in-sample fit is unaffected.

* **`?bartisan-package` is now a map** rather than boilerplate: a table of what to
  reach for -- and which package it lives in -- for predicting, intervals,
  effects, partial dependence, importance, convergence, fit checks and model
  comparison, plus a short orientation for readers new to Bayesian modeling and a
  paragraph on what a flexible fit does *not* buy you.

* **`weibull_aft()` is about 15% faster** and the result is bit-identical.
  `AFTFamily` was the one survival family still declared against the plain
  virtual `Family` base rather than `Concrete<AFTFamily>`, so its leaf
  accumulation ran four non-inlinable virtual calls per observation; and it had
  no `score_info_unit()`, so the score and the information came from two separate
  calls, each forming the residual and its exponential again. Both fixed: 9.7 to
  8.4 seconds on 700 observations.

  The remainder is inherent to the exponential-form Laplace fit rather than a
  defect -- `poisson()` takes the same path and pays a comparable 3x over the
  closed-form quadratic draw. Making `weibull_aft()` genuinely fast would need a
  mixture-of-normals augmentation of its Gumbel error, which is an approximation
  where every augmentation in the package is currently exact; that is a design
  decision and is recorded rather than taken.

* **A `Surv` response now defaults to `dpm_aft()`** rather than `weibull_aft()`.
  The change follows the comparison in `vignette("survival")`, where `dpm_aft()`
  was best or tied-best on four of six data-generating truths and never worse
  than third, matching the correctly specified family to the third decimal on the
  two truths where one existed; and it follows the speedup below, which makes it
  the second-cheapest survival family. An unnamed survival fit is now both more
  accurate on average and about three and a half times quicker.

  **This changes what the reported predictor means.** `weibull_aft()`'s is a log
  time ratio; `dpm_aft()`'s is the conditional mean of log T, which is a
  time ratio only insofar as the error density is symmetric -- and symmetry is
  exactly what `dpm_aft()` declines to assume. Name `weibull_aft()` to get the
  old behavior, and name it when a hazard-ratio reading is wanted, since it is
  the one family that is both an accelerated failure time and a proportional
  hazards model. `predict(type = "survival")` and `type = "response"` are
  unaffected: both are on the same scale for every survival family, which is why
  they are the recommended way to report one.

  `dpm_aft()` cannot take prior weights, so a weighted survival fit with no
  family named is now an error naming the alternatives rather than a silent
  substitution -- the same rule `dpm()` has carried for numeric responses.

* **`dpm_aft()` is about six times faster**, from 13.3 seconds to 2.2 on 700
  observations with 50 trees and 500 draws after 500 warmup. None of this came
  from the sampler, which was never the bottleneck: profiling put `.Call` at 11%
  of the fit and `stats::pnorm()` at 66%.

  The cost was in `response_scale()`, which finds the median survival time by
  bisecting the mixture's survival function. The median being solved for is a
  property of the *error* distribution, which every observation shares -- only
  the additive shift differs between them -- but the bisection was carried out on
  a vector with one entry per observation, every entry holding the same number.
  At 700 observations, 60 bisection steps and 500 draws that is 588 million
  `pnorm()` evaluations to produce 500 scalars. Bisecting the scalar and adding
  it to the row recovers a factor of six; a second change hoists the mixture
  matrix out of the bisection, since it was rebuilt on all sixty steps of each
  draw. `dpm_survival()` gained an optional `components` argument for that.

  Both changes are pure redundancy removal and the fitted values are bit-identical
  to the old ones, which is checked against the previous algorithm run verbatim.
  `dpm_aft()` is now the second-fastest survival family, ahead of `ph()` at 6.5
  seconds and `weibull_aft()` at 9.7 -- the latter being the current default for a
  `Surv` response.

* **New vignette, `vignette("survival")`.** The survival material had outgrown
  its section in `vignette("families")`, which is meant to be read for a
  practical choice rather than as a survey. The new vignette takes the five
  censored-response families on their own: the three estimands they report and
  why a contrast in `type = "survival"` is usually the one to want, which shapes
  of hazard each family can and cannot represent, a comparison across six
  data-generating truths at 50 trees and 500 draws, a censoring sweep from none
  to 70%, the discrete-time route for non-proportional hazards, and three traps
  -- the log-score measure mismatch and its Jacobian correction, `ph()`'s
  predictor being a contrast rather than a level, and the bin-count
  insensitivity. Figures are drawn with ggplot2, now used in a vignette rather
  than only in `_dev/`. The simulations are reproducible from
  `_dev/survival-sim.R`, `_dev/survival-bins.R` and `_dev/survival-results.R`;
  the vignette reads their saved output so that it builds without refitting.
  `vignette("families")` keeps a short survival section pointing at it.

* **`ph()` is the one proportional hazards family, and `num_bins` is not a knob.**
  A `coxph()` family was built -- CoxBART from Linero, Basak, Li and Sinha (2022),
  which reaches the same model by augmenting one baseline jump per event time and
  integrating them out -- and then removed, because measuring it settled the
  question it was there to raise.

  Sweeping `num_bins` over a sixty-fold range, from 4 to 250, against a baseline
  hazard that turns over and against a Weibull one, moved the error in the
  survival function between 0.035 and 0.048 and the error in the log hazard ratio
  between 0.135 and 0.185, with no trend in either and every difference inside the
  replicate-to-replicate spread. **The estimates do not depend on the bin count.**
  What does depend on it is the effective number of parameters, which grows from 17
  at four bins to 206 at 250 -- and that is what the default is protecting. At one
  piece per event time the effective count exceeds the sample size and more than
  half the Pareto-k diagnostics go above 0.7, so `loo()` and `waic()` stop
  working. The bin count is the regularization that keeps the model's own
  likelihood usable, not an arbitrary imposition, and `num_bins` is now documented
  as an advanced argument for checking that rather than for tuning.

  A second sweep, on the truths used in `vignette("survival")`, refines this
  slightly: the plateau runs from 4 bins to 100 rather than all the way to 250.
  At 250 the error rises by about 20% in both truths and both metrics, and the
  Pareto-k diagnostics begin to be troubled. The conclusion is unchanged and the
  default is nowhere near that edge -- the shape is a wide flat plateau with a
  cliff far past the default, not a peak to be found.

* **New family: `dpm_aft()`**, the fully nonparametric accelerated failure time
  model of Henderson, Louis, Rosner and Varadhan (2020):
  log T = m(x) + W with W a mean-constrained Dirichlet process
  mixture of normals, and right-censored log-times imputed each sweep. It is
  `dpm()`'s error model joined to the survival families' censoring, and both
  halves already existed -- the Polya urn, the atom draws, the concentration, the
  centering that makes the predictor the conditional mean of log T, and
  `error_density()` are all inherited unchanged. What is added is the imputation,
  which draws a censored log-time from the component it currently sits in
  truncated below at its censoring time, and an observed-data likelihood that
  credits a censoring with the mixture's survival rather than its density.

  Measured on 700 training and 700 test observations against three error shapes,
  as a held-out log score and the RMSE of S(t | x):

  | Errors | `dpm_aft()` | best fixed-error family |
  |---|---|---|
  | bimodal | **-607, 0.029** | -818, 0.098 (`lognormal_aft()`) |
  | log-normal | -438.1, 0.0264 | **-438.0, 0.0266** (`lognormal_aft()`) |
  | heavy tailed | **-509, 0.036** | -516, 0.040 (`loglogistic_aft()`) |

  So it gains a great deal when the error is badly shaped -- 210 log points and a
  third of the error in S(t | x) on a two-component error -- and costs
  nothing when a single normal is right, where it and `lognormal_aft()` are within
  0.1 log points of each other. The same property `dpm()` has against
  `gaussian()`. It refuses prior weights, as `dpm()` does.

  This error model is the more flexible of the two: the paper's is a location
  mixture with one common scale, and this is a location-scale mixture.

* **Note: `predict(type = "density")` is not on the same scale for every survival
  family.** The accelerated failure time families, `dpm_aft()` included, report
  the density of log T; `ph()` and `coxph()` report the density of
  T. The two differ by the sum of log t, so held-out log scores are
  comparable within each group and not across them -- on one comparison the
  correction was 1042 log points and reversed which family looked better.
  S(t | x) from `predict(type = "survival")` is comparable throughout
  and is the safer cross-family metric.

* `summary()` shows the ends of a long block of nuisance parameters and says how
  many it omitted, rather than printing hundreds of rows.

* **`predict(type = "survival", times = ...)`** returns the survival function
  S(t | x) for the accelerated failure time families and for `ph()`.
  This closes a gap rather than adding a convenience: every family whose response
  is discrete already returned its full predictive distribution through
  `type = "prob"`, while the survival families returned only a point summary --
  the median. The curve is what survival analysis is usually asked for, and
  computing it from the draws by hand meant depending on how the baseline is
  stored. One column per time, or a draws by rows by times array with
  `draws = TRUE`.

* **Fixed: *marginaleffects* did not work for *any* survival family.** A
  survival response is a two-column matrix and a model frame keeps it as a single
  matrix column, which the data.table conversion inside
  *marginaleffects* cannot hold -- so `avg_comparisons()` and its siblings
  failed for every accelerated failure time fit, whatever `type` was asked for,
  with an error naming data.table rather than the cause. `get_data()` now splits
  such a column into ordinary ones.

  Together with the above, the usual survival estimand is now reachable:
  `avg_comparisons(fit, variables = "trt", type = "survival", times = 1)` is the
  difference in one-year survival. One time per call, and
  *marginaleffects* warns that it does not recognize `times` while passing it
  through, since it has no hook for registering an argument.

* **New family: `ph()`**, proportional hazards with a piecewise-constant
  baseline. The predictor is a log *hazard* ratio and the baseline is free to
  take any shape, where `weibull_aft()`'s hazard is monotone and the other two
  survival families' are unimodal. The bin hazards are drawn from their exact
  gamma conditionals and reported as `lambda1`, `lambda2`, ... in `fit$aux`;
  `num_bins` defaults to about the cube root of the sample size. Following Basak,
  Linero, Maringe and Rubio (2024).

  Measured against a proportional-hazards truth, 800 training and 800 test
  observations, it cuts both ways and the sizes are worth knowing: where the
  baseline hazard turns over, `ph()` was worth about **250 held-out log points**
  over the best accelerated failure time family and halved the worst-case error
  in `S(t | x)`; where the baseline really was Weibull, `weibull_aft()` was worth
  about **96 log points** back, though the two were level on `S(t | x)`. What
  `ph()` does not buy is the risk *ordering*, which the parametric families
  already recover about as well as a model given the true baseline.

  This retracts a claim the documentation used to make. Cox's *partial*
  likelihood genuinely cannot be used -- it couples observations through risk sets
  and so does not decompose over the observations reaching a leaf -- but the full
  likelihood of the piecewise-exponential model does decompose, and it is of the
  same exponential form the Poisson family already uses, so a leaf update is one
  pass over the node. Proportional hazards was reachable all along.

* `prepare_surv()` now also returns the times themselves, not only their logs.

* **`custom_family()` can draw nuisance parameters.** Name them with
  `aux_names`, give them a starting order of magnitude with `aux_start`, and
  `logdens` gains a third argument holding their current values. They come back
  in `fit$aux` under those names and are covered by `summary()` and `fit$rhat`
  like any other family's. So a Gaussian written out by hand no longer has to fix
  its scale inside the closure:

  ```r
  custom_family(
    logdens   = function(y, eta, aux) dnorm(y, eta[, 1], exp(aux[1]), log = TRUE),
    aux_names = "log_sigma", aux_start = 0)
  ```

  Started seven times away from the truth, that recovers the scale to within
  Monte Carlo error of what `gaussian()`'s own conjugate step gives -- matching
  it on the posterior mean, the posterior standard deviation and the effective
  sample size alike.

  There is no prior or bounds argument, because a nuisance parameter is carried
  as an additive predictor whose forest is pinned at depth zero -- one tree that
  can never split, so the forest is one scalar -- and drawn by the same
  Laplace-plus-Metropolis step as any leaf. A restricted range is handled as it
  would be for a real predictor, by writing the transform into `logdens`.

* **The Fisher-scoring fit now carries a trust region**, which makes the sampler
  robust to a leaf that starts far from its conditional mode. Two things went
  wrong there and both are fixed: Newton's method could be thrown far past the
  mode on a sharply asymmetric target and then crawl back too slowly to arrive
  within its step limit; and the Laplace proposal, being an independence proposal
  centred on the mode, could not be accepted from far out however much the target
  improved, because the reverse density at the current value was astronomically
  small. Either one left the leaf frozen at its starting value for the whole run
  -- a permanent trap rather than slow mixing. Both caps are fixed constants, so
  the fit stays a deterministic function of the state and the moves stay
  reversible, and neither binds for a chain sitting near its mode, so no existing
  family's timing or draws changed measurably.

* Fixed: the error for prior weights with `dpm()` was built from a malformed
  message and failed to print.

* **`Beta()` and `ordbeta()` are four times faster.** They were the slowest
  families in the package, and not because of their target: they used the same
  number of Fisher-scoring passes as any other general-target family and six
  times the time, all of it in four `digamma`/`trigamma` calls per observation
  per pass. Both combinations those functions appear in are functions of the mean
  alone once the precision is fixed, which it is for a whole sweep, so they are
  now tabulated per sweep and interpolated. On 500 observations with 50 trees,
  `Beta()` went from 91.5 seconds to 20.0 and `ordbeta()` from 104.1 to 24.9,
  with effective samples per second up five- and nearly four-fold. The same
  posterior: over two seeds each the two versions differ by less than either
  differs from itself, and the precision's posterior mean agrees to a hundredth
  of a posterior standard deviation.

* `ordbeta()`'s cutpoint update no longer rebuilds the beta density on every
  slice evaluation, which it does not depend on. Exactly the same draws --
  slice sampling is invariant to an additive constant in the log density -- for
  a fifth off the remaining time.

* **Fixed: `Beta("probit")` and `Beta("cloglog")` silently fitted the logit
  model.** Only the logit link is compiled, but all three were listed as native,
  so the engine ignored the choice and `predict()` then back-transformed as
  though it had not, giving wrong fitted means. The two are now reached by
  composition, as an uncompiled link is for any other family, and any other link
  `stats::make.link()` understands works too.

* **New family: `Beta()`**, beta regression for a response strictly inside the
  unit interval. A forest on the link of the mean and a drawn precision, with the
  `logit`, `probit` and `cloglog` links. It is capitalized for the same reason
  `stats::Gamma()` is -- `base::beta()` is the beta function, and masking it
  would be worse than a capital letter. `ordbeta()` remains the family for a
  response that can sit exactly at 0 or 1; choose between them on whether the
  response *can* reach a boundary, not on whether it happens to in the sample.

* **The gamma family now overrules any link but the log one**, with a message
  saying so. Base R defaults `stats::Gamma()` to the inverse link, because that
  is the canonical link for the gamma; but the canonical link needs a positive
  mean, and this sampler's additive predictor is unconstrained, so a draw that
  wanders non-positive has no gamma density at all -- it is rejected during
  fitting, and `predict()` returns `NaN` there. Over eight replicates on
  heavy-tailed data, five had such draws. The log link was also faster and more
  accurate in every setting measured. Writing `Gamma("log")` silences the
  message; `stats::Gamma()` itself is untouched, so attaching *bartisan* still
  cannot change what `glm()` does. One consequence: composed links are no longer
  available for the gamma, so a genuinely wanted `sqrt` link now needs
  `custom_family()`.

* **The survival models are much faster.** `lognormal_aft()` and
  `loglogistic_aft()` now impute each censored failure time above its censoring
  time, which makes the leaf-level target quadratic. Right-censoring was the
  expensive part: an observed failure contributes a density in the additive
  predictor and a censored observation a survival function, and the two have
  different shapes, so the target had no form the sampler could exploit and every
  trial value of a leaf parameter cost its own pass over the data. Worth 8 to 30
  times the raw speed and **5 to 22 times the effective sample size per second**,
  the largest gains in the package alongside `ordinal("probit")`. Unlike the
  exponential form, this route works under soft rules, so the gain is there at the
  default gate. The same posterior either way: over censoring from 10% to 70% the
  scale's posterior mean agreed to within a tenth of its posterior standard
  deviation, and the predictor agreed more closely between the two samplers than
  either agreed with a second seed of itself.

* **`weibull_aft()` is faster under hard rules**, by a different route and
  without any augmentation: its likelihood is already of the form
  `c + a * eta + b * exp(-eta / sigma)`, censored observations included, so the
  sampler can collapse a leaf update to a single pass. 6.7 seconds to 4.2 at 50
  trees. It reproduces the general path to 1e-15.

* **`location_scale()`'s scale forest is faster under hard rules**, 7.7 seconds
  to 3.6, for the same reason: its target is of that form at rate -2. Reaching it
  meant switching the log-scale predictor's reported information from the expected
  value to the observed one, which is what the closed form reads its coefficients
  off; both are exact, since the Laplace fit is a proposal.

* The `augment` argument accepts `"aft"`, and `TRUE` includes it.

* **`dpm()`'s additive predictor is now the conditional mean.** Nothing in the
  model forces the error mixture to be centered, so the sampler works in a chart
  where only the sum of the predictor and the error mean is identified. Reporting
  is now done in the chart where the mixture has mean zero, which puts the whole
  conditional mean on the predictor: `type = "link"` and `type = "response"` agree
  exactly, and `dpm()` is directly comparable with `gaussian()`. On a skewed
  example the level of the predictor had a standard deviation across draws of
  0.750 before and 0.021 after, and its bias against the truth fell from 0.144 to
  -0.021. The shift is applied to the recorded predictor, leaf values and mixture
  together, so the density is untouched -- the log likelihood rebuilt in R from
  the stored pieces still matches the sampler's own to 3e-13.

* `fit$aux` for `dpm()` reports **`center`** in place of `error_mean`: the shift
  that was taken out, which is bookkeeping rather than an estimate, and which
  says how far from symmetric the fitted error came out. The error mean in the
  reported chart is zero by construction. `error_density()` is centered to match.

* **A numeric response now defaults to `dpm()`** rather than `gaussian()`, when
  it has ten or more distinct values. `dpm()` matches `gaussian()` when the errors
  really are normal and beats it otherwise, so there is no error distribution
  where the old default was better. A coarser numeric response still gets
  `gaussian()`: a mixture cannot separate an error distribution from a mean when
  the response takes a handful of values, and at three distinct values `dpm()` was
  3.8 times worse on held-out error. Since `dpm()` refuses prior weights, a
  weighted fit with no family named is now an **error** naming the alternatives,
  rather than a silent substitution.

* `augment` is documented as an **advanced setting** rather than a modeling
  decision. A rewriting targets exactly the same posterior as the direct
  likelihood; what it changes is how fast the sampler gets there. Nobody should
  have to think about it, which is what the default is for.

* **`dpm()` is now documented as the default worth reaching for on a numeric
  response**, with the comparison behind it. It does not pay for its
  flexibility: on normal errors, where `gaussian()` is exactly right, it came out
  slightly ahead on both held-out error and log score at the same time to one
  decimal place. On heavy-tailed, skewed and bimodal errors it is ahead by a
  great deal -- 0.050 against 0.154 in held-out RMSE on bimodal errors at a
  thousand observations, a factor of three. The reasons to prefer `gaussian()`
  are not statistical: prior weights, which `dpm()` refuses; an identified
  additive predictor, which `dpm()` does not have; one interpretable `sigma`; and
  1.4 times the speed.

* **Using `ordinal()` on a continuous outcome is documented as a method rather
  than a fallback**, with `type = "mean"` for prediction. The cutpoints absorb
  the marginal distribution and the forest explains only the ordering, so nothing
  is assumed about the error distribution and the model for the cumulative
  probability is invariant to any monotone transformation of the response. It had
  the lowest error and the only above-nominal coverage in the one setting where
  the other two families are misspecified, heteroskedasticity. **Bin the outcome
  onto twenty-odd quantiles first**: one cutpoint per distinct value cost 73
  seconds against 2.7 for `gaussian()` at a thousand observations, where
  twenty-five bins was sixteen times faster *and* slightly more accurate.

* **Write `Gamma("log")`, not `Gamma()`.** Base R defaults the gamma family to
  `link = "inverse"`, which is the worst link for this sampler: its inverse sends
  a negative predictor to a negative mean, whose log is not a number, so the
  proposal is rejected -- dozens of times per fit. On 600 observations with 50
  trees, `stats::Gamma()` took 7.2 seconds against 3.8 for the log link and
  fitted the mean slightly worse (RMSE 0.664 against 0.606). `stats::Gamma()` is
  left exactly as base R defines it, so that attaching *bartisan* cannot change
  what `glm()` does; the string form, `family = "Gamma"`, is this package's own
  spelling and resolves to the log link. `bartisan()` reports any link whose
  inverse does not cover the additive predictor when the fit starts, which is the
  case that catches this one.

* **`Gamma_shape()` is removed.** Its only reason to exist was a `shape`
  argument for fixing the gamma shape, and a caller who knows a gamma shape is
  rare. The shape is drawn, as it was before.

* `bartisan()` **says so when a link's inverse does not cover the whole additive
  predictor** -- `Gamma("inverse")`, `Gamma("identity")`, `poisson("identity")`.
  Those proposals are rejected rather than breaking the chain, so the fit is
  valid, but it is slower and less accurate and the reason was previously
  invisible. Compiled links and composed links whose inverse does cover the line,
  such as `binomial("cauchit")`, say nothing.

* `bartisan_control()` is organized into three groups, stated up front and marked
  on each argument: **modeling decisions** (`num_trees`, `gate`, `sparsity`, `k`,
  `bandwidth`, the chain lengths, `augment`, `x_transform`), **advanced
  settings** (`gamma` through `num_print`), and three toggles that exist **for
  internal validation** (`block_eval`, `exact_quadratic`, `generic_accumulate`).

* **`num_trees` takes one value per additive predictor**, not just one per fit. A
  scalar is recycled, so nothing changes for a single-forest family. This is
  worth using: `location_scale()` spends about 90% of its time on the scale
  forest, and a variance surface needs less capacity than a mean surface, so
  `num_trees = c(50, 10)` ran in 5.9 seconds against 14.5 for `c(50, 50)` with
  the same accuracy on both surfaces and the same held-out log score. The leaf
  prior divides by the square root of each forest's own tree count, so shrinking
  one forest leaves the prior on the sum it forms unchanged. `print()` now says
  "2 forests of 50 and 10 trees" when the counts differ.

* **`soft` is gone; `gate` decides both questions.** `gate = "hard"` (or
  `"step"`) gives the step functions of standard BART, and `"smoothstep"`,
  `"smootherstep"` and `"logistic"` give soft rules and name the gate's shape.
  They were always one decision: a hard rule has no gate shape to choose.

* **`sparsity` stands in for four hyperparameters.** `TRUE`, the default, is the
  Dirichlet variable-selection prior of Linero (2018); `FALSE` gives every
  predictor the same splitting probability, which is classic BART; and `"none"`,
  `"weak"`, `"moderate"` and `"strong"` name four strengths. It sets `update_s`,
  `update_alpha`, `alpha_shape_1` and `alpha_shape_2` together, and any of those
  supplied directly still wins.

* `?bartisan_control` gains two sections of measurements: how many trees to use,
  including that soft rules peak at about 20 and get worse past that while hard
  rules keep improving to 200, and what the sparsity prior costs. **A correction
  there:** a larger `num_trees` does *not* remove the atom at zero in a
  contrast's posterior. With the prior on, the contrast was exactly zero in 20%
  of draws at 50 trees and 18% at 200; `sparsity = FALSE` is what removes it. The
  earlier claim came from a single chain that happened to look settled.

* A new vignette, `vignette("families")`, documents every response family at
  length: how to choose one, how to choose among the links a family offers, how
  to choose among the zero-inflated and accelerated failure time models, and what
  `custom_family()` can and cannot do. The help page for the families keeps the
  table, the identification facts that output could be misread without, and a
  pointer to the vignette. Both vignettes now carry their citations in
  `vignettes/references.bib` rather than a hand-written list.

* Corrected the documentation of the multinomial probit's latent correlations.
  It described them as attenuated towards zero, which was wrong, and the numbers
  printed alongside the claim contradicted it. Measured properly, the correlation
  is *weakly identified*: at 900 observations a true value of zero came back
  anywhere from -0.60 to 0.33 across draws of the data, posterior intervals ran
  up to 1.07 wide on a parameter confined to (-1, 1), and the sweep over the truth
  was not monotone. By 3000 observations it tracks the truth to within about 0.2.
  The advice is now to read the fitted probabilities and not the covariance, and
  the corresponding test asserts only what holds.

* Documented why `marginaleffects::avg_comparisons()` can report a contrast of
  exactly zero. The posterior of a contrast has an atom at exactly zero, because
  in a draw where no tree splits on the variable the two counterfactual
  predictions are identical; the Dirichlet sparsity prior makes such draws common
  for a weak predictor; and *marginaleffects* centers a posterior at its
  median, which then lands on the atom. Neither package is computing anything
  wrong, and the help page now says what to look at instead.

* Initial version. `bartisan()` fits Bayesian additive regression trees for
  response distributions outside the conditionally conjugate case, using the
  Laplace-approximation reversible-jump sampler of Linero (2025), with optional
  soft decision rules following Linero and Yang (2018).

* Sampler and prior settings may be passed to `bartisan()` directly, as in
  `glm()`: any argument of `bartisan_control()` given in `...` is merged into
  `control`, overriding a value set there. A name that is not an argument of
  `bartisan_control()` is an error rather than being silently ignored.

* `ordinal("probit")` is 14 times faster with soft rules and 30 times faster
  with hard ones, at the same accuracy, through the latent normal of Albert and
  Chib (1993). Conditional on it the leaf target is quadratic, so the sampler
  takes the closed form instead of Fisher scoring plus a Metropolis ratio, and
  the two cumulative-normal evaluations per observation per pass disappear. On
  by default via `bartisan_control(augment = )`. The cutpoints are drawn from the
  ordinal likelihood with the latent variables integrated out and the latent
  variables redrawn immediately afterwards -- a partially collapsed Gibbs
  sampler -- because their conditional given the latent variables mixes worse the
  larger the sample. `ordinal("logit")` has no such representation and is
  unchanged.

* `bartisan_control(gate = "smoothstep")` offers a soft decision rule whose gate
  has bounded support: past `4.06 * bandwidth` from the cutpoint an observation
  takes one side outright, so the other subtree is never visited and the gate
  costs no `exp()`. Measured 1.4 times faster than the logistic gate with no loss
  of accuracy. The constant makes `bandwidth` mean the same amount of smoothing
  under either gate.

* `bartisan_control(bandwidth_every = )` draws each tree's bandwidth every *k*-th
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

* `bartisan_control(gate = "smootherstep")` adds a twice-differentiable bounded
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

* **`dpm()`**, which is DPMBART (George, Laud, Logan, McCulloch and Sparapani
  2019): a numeric response with the sum of trees for its mean, as `gaussian()`
  has, and a Dirichlet process mixture of normals for its errors instead of one
  normal. Each observation gets its own error mean and variance, drawn from a
  Dirichlet process, so the error distribution comes out as whatever mixture the
  data ask for -- heavy tailed, skewed, bimodal -- rather than the one normal
  BART is committed to. `error_density()` reports the estimated density with a
  pointwise interval, which is the object the method exists to produce.

  Conditional on the mixture the target is still exactly quadratic in the
  predictor, so the leaf draw is the closed form and the only cost over a
  Gaussian fit is the mixture update: measured at 1.7 s against 1.3 s at
  n = 1000, 50 trees, 500 + 500 draws, where the paper reports the total roughly
  doubling.

  Measured against `gaussian()` and `location_scale()` on the paper's three error
  distributions plus a heteroskedastic one, held-out predictive log score:

  | errors | `gaussian()` | `location_scale()` | `dpm()` |
  |---|---|---|---|
  | normal | **-2097** | -2100 | **-2097** |
  | heavy tailed (t3) | -2642 | -2607 | **-2467** |
  | skewed (log gamma) | -2680 | -2680 | **-2456** |
  | heteroskedastic | -1886 | **-1734** | -1853 |

  So it costs nothing when the errors are normal, wins by a wide margin when they
  are heavy tailed or skewed, and is *not* the family for heteroskedasticity --
  a mixture makes the error distribution flexible but keeps it the same
  distribution at every `x`, which is what `location_scale()` is for. On the
  skewed case it also cut the interval width for the regression function from
  1.59 to 0.99 while holding coverage, which is the paper's point about BART
  intervals being wrong rather than merely wide.

  One thing to know: nothing forces the mixture to be centered, so only the sum of
  the fit and the error mean is identified. `predict(type = "response")` is that
  sum and is the quantity to compare against a truth; `type = "link"` is the sum
  of trees alone and carries the drift. Measured, the two move against each other
  with a correlation of -0.987 and their sum has a sixth of either one's standard
  deviation.

* **`family` may be left alone.** It defaults to `NULL`, and the family is then
  read off the response: `weibull_aft()` for a `Surv` object, `ordinal()` for an
  ordered factor, `binomial()` for a logical, a two-level factor or numeric zeros
  and ones, `multinomial()` for a factor with more than two levels, and
  `gaussian()` for anything else. A message says which was chosen; naming
  `family` silences it. A count is deliberately *not* read as `poisson()`, and a
  numeric response whose two values are not zero and one is not read as
  binomial -- in both cases the guess would be a modeling decision made on the
  caller's behalf rather than a reading of the response's type.

* **`multinomial("probit")`**, a categorical model whose latent utilities are
  *correlated*, which a multinomial logit cannot express at all. It is a link on
  the existing `multinomial()` rather than a separate family, as `binomial()` and
  `ordinal()` do it. The outcome is
  the largest of several latent utilities; differencing against a reference
  category leaves one forest per remaining category and a covariance matrix
  Sigma, drawn from its inverse Wishart conditional and normalized by the
  trace constraint of Burgette and Nordheim (2012). With two categories that
  constraint leaves the fit identical to binary probit regression.

  Conditional on the latent utilities the target is *exactly quadratic in every
  component*, so the leaf draw is the closed form rather than a Laplace
  approximation, and a fit costs about twice a multinomial logit's. The sampler
  is Algorithm P2 of Xu et al. (2025), which fits the forests to the normalized
  utilities; the earlier sampler of Kindo, Wang and Pena (2016) fits them to the
  unnormalized ones and grows much deeper trees as a result.

  Two things follow from the multivariate normal. The likelihood has no closed
  form -- a category probability is a Gaussian orthant probability -- so it is
  simulated, with `replicates` draws; the Monte Carlo error is per posterior draw
  and averages down over them. And the estimated correlation is *attenuated*,
  because a nonparametric mean absorbs part of the dependence: read its sign and
  its ordering rather than its magnitude. The inverse Wishart degrees of freedom
  are deliberately not exposed, since raising them with the scale matrix held at
  the identity pushes the correlation towards one rather than towards zero.

* **The zero-inflated families are 4 to 10 times faster in effective samples per
  second**, through a rewriting with two latent variables. What made them slow
  was the zero: `log[pi + (1 - pi) P_0(mu)]` is a mixture of the two components,
  so neither the count predictor nor the inflation predictor had an exploitable
  shape and every trial value of a leaf parameter cost its own pass. Introducing
  the indicator of *which* component produced the observation separates them --
  the count forest then sees a plain Poisson or negative binomial and the
  inflation forest a Bernoulli logistic likelihood, which Polya-Gamma handles.
  `zi_negbin()` puts the Poisson-gamma rate on top of that, so its count forest
  gets the exponential form too. On by default; `augment = FALSE`, or a character
  vector without `"zip"` / `"zinb"`, turns it off. Measured at 3.9x (soft rules)
  and 10.1x (hard) in ESS per second for `zi_poisson()`, and 5.8x and 7.9x for
  `zi_negbin()`, at unchanged accuracy.

* **`multinomial()` is augmented by default.** It was available by name only, on
  a measurement that put it at 1.6x in effective samples per second -- a real but
  modest gain bought with a severe loss of mixing. Re-measured after several
  rounds of sampler work it is 9.6x (hard rules) and 10.1x (soft), with the
  mixing cost at 0.66x and 1.09x rather than the 0.37x recorded before. The old
  number was taken on a different problem and with a different effective-sample
  measure and could not be reproduced; the new one is stated with its
  configuration in `?bartisan_control`.

* `bartisan_control(gate = )` now defaults to `"smoothstep"` rather than
  `"logistic"`. The bounded gate is about 1.4 times faster with no measured loss
  of accuracy, which makes it the better default; `"logistic"` is still there for
  Linero and Yang's (2018) original and for an infinitely differentiable fit.

* **Posterior predictive draws**, and the fit-quality packages they unlock.
  `posterior_predict()` draws replicate outcomes from the fitted model,
  `posterior_epred()` and `posterior_linpred()` give the mean and the additive
  predictor, and `log_lik()` gives the draws-by-observations matrix of
  log-likelihood contributions -- all on the *rstantools* generics that
  *brms* and *rstanarm* use. `simulate()` is the same thing in base R's
  shape, returning factors for a response that was one. On top of them:
  `loo::loo()` and `loo::waic()`, `bayesplot::pp_check()` for any of that
  package's posterior predictive checks, `posterior::as_draws()` for the scalar
  parameters, and `performance::model_performance()` and `performance::r2()` for
  the Bayesian R-squared of Gelman et al. (2019). `fitted()`, `residuals()`,
  `weights()` and `sigma()` round out the accessors. See `?bartisan-interop`,
  which also says why a leave-one-out estimate is strained for a model this
  flexible.

  Every family's sampler was checked against the log density the C++ engine
  reports -- an independent statement about the same distribution -- at 20
  configurations covering every family and link.

* `marginaleffects` now reaches `type = "mean"` and `type = "stdlv"` for ordinal
  fits, which `predict()` had but `get_predict()` was not mapping. `"probs"`,
  `"lp"` and `"lv"` are accepted as aliases for the names those quantities go by
  elsewhere in *marginaleffects*, and `values` is passed through to
  `predict()`.

* Fitting with a link the engine does not carry natively -- the default inverse
  link of `Gamma()` is the common case -- no longer emits a stream of
  "NaNs produced" warnings. A negative predictor there implies a negative mean,
  which the sampler already handles by rejecting the proposal; the only thing the
  warning did was print one line per rejection. Draws are unchanged.

* `predict(type = "stdlv")` now works for a binomial fit as well as an ordinal
  one: both responses are a threshold crossing of a latent variable, so both have
  a latent standard deviation to divide by. The probit, logit and complementary
  log-log links are covered, being the ones whose error distribution is named.
  The complementary log-log error enters the two families with opposite signs --
  an ordinal model is `P(Y <= k) = G(c_k - eta)` and so `y* = eta + e`, while a
  binomial model is `P(Y = 1) = G(eta)` and so `y* = eta - e` -- which is also
  why a two-category ordinal complementary log-log fit is not the same model as
  a binomial one.

* The prior scale for the Gaussian residual standard deviation now reads the
  prior weights. It is estimated from a linear fit, and previously that fit was
  unweighted, so a weighted analysis anchored the prior on the wrong
  observations. Only the relative sizes of the weights matter, so a constant
  weight -- including the implicit one -- gives exactly the old answer.

* Two prediction types for ordinal models, following
  *WeightIt*. `predict(type = "mean")` reports the probabilities weighted by
  the category labels read as numbers, so a response with levels `"1"`, `"2"`,
  `"4"` gets a mean between one and four; `values` says what the categories are
  worth when the labels are not numbers. `predict(type = "stdlv")` reports the
  additive predictor divided by the standard deviation of the latent variable it
  indexes, which is what puts fits with different links, or different amounts of
  signal, on one scale.

* **Group-level random intercepts**, written in *lme4*'s notation:
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
  depends on the shape of the data, and `?bartisan` says which way round: with few
  large groups the predictor route measures better, and the random intercept wins
  once there are many small groups.

* Ordinal cutpoints are now documented against `MASS::polr()` with a worked
  example, and there is a test: the chart bartisan reports in is the one `polr()`
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

* Support for the *marginaleffects* package, so predictions, comparisons,
  slopes and hypothesis tests on any of them can be computed directly from a fit.
  Uncertainty comes from the posterior draws rather than from a delta-method
  approximation. See `?bartisan-marginaleffects`, which also documents why slopes
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

  - `augment` in `bartisan_control()` rewrites a likelihood as the margin of a
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

* `_dev/benchmark.Rmd` times bartisan against dbarts, BART and bartMachine on the
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
