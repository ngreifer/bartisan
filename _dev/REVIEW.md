# Pre-ship review, 2026-09-21

A pass over the package for bugs and for speed, made against the working tree as
it stood that morning (the vignette edits and the `effect_by()` / binary-plot
changes included). Small bugs were fixed in place; everything larger is written
up here with enough of the mechanism that it can be picked up without this
session's context. The standard for a fix was the one `_dev/SHIP.md` uses: a few
lines, no change to what a correct call returns.

Every change to the sampler was checked for **bit-exactness** against the
unmodified build: the same seed must give the same stored draws to the last bit,
which is a far stronger check than the test suite and is what let the
optimizations be accepted or rejected without argument. The one change that is
not exact by construction (the nuisance-update rewrite, item F below) turned out
to be exact in fact: the `Gamma()` and `negbin()` draws matched to the bit on all
three seeds, because a slice sampler's every decision is a threshold comparison
and none of them landed within rounding distance of the threshold. That is a
property of these runs rather than a guarantee, so it is also reported that
every configuration the rewrite does not touch stayed exact.

## Verification record

| Step | Result |
| --- | --- |
| Full suite on the unmodified tree (`NOT_CRAN=true`) | 39 files, 485 tests, 2479 expectations, 0 failures, 0 errors, 6 warnings, 0 skips, 9m 48s |
| Those six warnings | all `warn_runaway_scale()` from the two multinomial offset tests in `test-bartisan.R`; fixture, not offset (see Fixed, 3) |
| Scoped suites after each change (`varying`, `interop`, `partial-dependence`, `estimate-effect`, `methods`, `bartisan`, `marginaleffects`) | 0 failures, 0 warnings; the two parallel-plan tests skip with no second worker |
| Full suite on the final tree (`NOT_CRAN=true`) | 39 files, 485 tests, 2487 expectations, 0 failures, 0 errors, 0 warnings, 0 skips, 15m 05s |
| Bit-exactness, 11 configurations x 3 seeds, unmodified vs. all exact changes | identical `eta`, `sigma_mu`, bandwidth, nuisance and log-likelihood draws on every row |

The benchmark harness is `_dev/review-bench.R`. It installs nothing: each build
goes into its own library with `R CMD INSTALL --library=`, and the script is run
once per library. It writes a hash of every stored draw beside the timing, which
is what the exactness column above is read from. Timings are medians of three
fits at n = 1000, six numeric and two factor predictors, 50 trees, 300 warmup and
300 kept sweeps, one chain. Two runs were made: builds back to back, and builds
interleaved within each configuration. The first run was misleading -- the
configurations that happened to run last lost 5-12% on *every* build, the
untouched `dpm()` among them -- so the table below is from the interleaved run.

## Fixed in this review

1. **`coef(fit, newdata = d)` returned `NaN` on every varying-coefficient
   fit** (`R/methods.R`, `coef.bartisan_fit()`). It handed `predict_eta()`
   `iterations = NULL`, and `predict_eta()` takes the iterations already
   resolved, as `predict()` hands them to it: `as.integer(NULL) - 1L` is a
   zero-length vector, the engine evaluated no draws, and `colMeans()` of a
   0-by-n matrix is `NaN`. Now resolved with `resolve_iterations()` like every
   other caller. Test added to `test-varying.R` (same rows give the same
   coefficients with and without `newdata`; the draws matrix has one row per
   stored draw). A remaining limitation, not fixed: a model fitted with an
   `offset` still needs one passed for `coef(newdata = )`, because
   `predict_eta()` insists on it, although the slope forests do not use it.

2. **`residuals()` and `r2_bayes()` on survival fits** (`R/interop.R`).
   `dpm_aft()` stores the log time and `predict(type = "response")` returns a
   median survival *time*, so `residuals()` subtracted across two scales.
   `ph()` stores the time, so its residual was at least on one scale, but both
   it and `r2_posterior()` treated every censored time as an observed value.
   `residuals()` now refuses `dpm_aft` and `ph` the way it refuses the
   categorical families, pointing at `type = "survival"` and the `km_overlay`
   check; `r2_posterior()` refuses all three survival families with a warning
   and `NULL`, and `model_performance()` skips RMSE for the two `residuals()`
   refuses. `residuals()` for `aft` is left as it was: it deliberately returns
   the log-scale residual against the predictor and says so in a comment.
   Test added to the survival `loo()` test in `test-interop.R`.

3. **Six warnings on every suite run** were `warn_runaway_scale()` firing in
   the two multinomial offset tests, including on the fit *without* an offset:
   ten trees and a hundred warmup sweeps are not enough for the leaf scale to
   settle. Silenced at the fixture with a comment saying why. The suite now
   reports zero warnings.

4. **`vignette("bartisan")` said `sex` and `race` are used "in fewer than half
   the draws"**; for the knitted fit the smallest `prop_used` is .52. One
   sentence changed to "little more than half".

5. **`predictions()` on one row of `newdata` with one survival time failed**
   (`R/marginaleffects.R`, `me_draws()`). A draws-by-rows-by-times array with
   one time was reduced with `drop = TRUE`, which with one row also drops the
   row dimension and hands marginaleffects a vector where it needs a matrix;
   `colMeans()` then failed inside `get_predict()` and marginaleffects reported
   it as being unable to extract the data. Confirmed on a `ph()` fit (one row
   fails, two rows work); now reduced with `drop = FALSE` and the two leading
   dimensions kept. Test added beside the survival estimand test in
   `test-marginaleffects.R`.

6. **Stale objects after a header edit** (`src/Makevars`). R's rules track only
   the source file, so after `node.h` changed, `random.o` -- whose source did
   not -- was reused and linked against a `Tree` constructor that no longer
   existed. `_dev/TASKS.md` records the same trap biting twice before. Two
   lines fix it durably: `$(OBJECTS): $(wildcard *.h)`, and `all: $(SHLIB)`
   *ahead of it*, because make takes the first target of the first rule it
   reads as its default goal and `Makevars` is read before R's `shlib.mk`.
   Without the `all:` line the goal became `RcppExports.o`, one object was
   checked, and the old shared library was installed as it stood -- which is
   exactly what happened to the first two builds of this review before it was
   caught. Verified: `touch src/node.h` now recompiles all eleven objects.

## Sampler changes, and what each measured

Every entry marked *exact* produced identical draws on all eleven
configurations. Timings are the interleaved run; a ratio above 1 is faster.

| Configuration | base | new (A-E) | new2 (+F) | new3 (+G) | new / base | new2 / new | new3 / new2 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `gaussian()`, hard rules | 0.43 | 0.37 | 0.38 | 0.37 | **1.15** | 0.99 | 1.01 |
| `binomial()`, `augment = FALSE` | 8.63 | 8.05 | 7.79 | 7.61 | **1.07** | 1.03 | 1.02 |
| `binomial()` (augmented, default) | 1.30 | 1.21 | 1.25 | 1.24 | **1.07** | 0.97 | 1.00 |
| `dpm()` | 1.25 | 1.19 | 1.18 | 1.16 | **1.05** | 1.00 | 1.02 |
| `Gamma("log")` | 5.27 | 5.10 | 4.92 | 4.89 | 1.03 | **1.04** | 1.01 |
| `poisson()` | 3.64 | 3.56 | 3.64 | 3.57 | 1.02 | 0.98 | 1.02 |
| `Beta()` | 18.16 | 17.95 | 17.91 | 17.61 | 1.01 | 1.00 | 1.02 |
| `negbin()` | 5.31 | 5.30 | 5.24 | 5.08 | 1.00 | 1.01 | 1.03 |
| `multinomial()`, `augment = FALSE` | 30.76 | 30.58 | 30.75 | 30.28 | 1.01 | 0.99 | 1.02 |
| `ordinal()`, `augment = FALSE` | 14.45 | 14.41 | 14.53 | 14.28 | 1.00 | 0.99 | 1.02 |
| `gaussian()`, soft rules (medians) | 0.84 | 0.96 | 1.09 | 1.07 | 0.87 | 0.89 | 1.02 |
| `gaussian()`, soft rules (minima of 9) | 0.825 | 0.763 | 0.773 | 0.738 | **1.08** | 0.99 | 1.05 |
| Total over the 33 fits | 269.7 | 266.1 | 266.2 | 260.9 | 1.014 | 1.000 | 1.021 |

Two readings of that table. First, the exact changes are worth 5-15% on the
paths where the leaf work is cheap -- the Gaussian and augmented families, hard
rules -- and nothing measurable where a leaf costs several Fisher-scoring passes
over transcendental functions, which is the non-augmented families; there the
sampler's time is in the family's arithmetic and not in anything these touched.
Second, the sub-second configurations are bimodal: on this machine a one-second
single-threaded process lands on a performance core or an efficiency core and
runs 30-40% slower on the second, so their medians flip with the run (the soft
Gaussian read +11% in the first run and -13% in the interleaved one) while the
minima -- the performance-core mode -- are stable and agree with the hard-rule
Gaussian beside it. The soft-Gaussian row is therefore given both ways, and
the nine-fit minima are the number to believe. `new2` differs from `new` only
in the two nuisance updates, so its other rows are a second estimate of the
noise: +-3% on the multi-second fits.

**A. Static dispatch for seven more families** (`src/family.cpp`; *exact*).
`Concrete<Derived>` is the CRTP base that inlines a family's `logdens_unit()`
and `score_info_unit()` into the per-leaf sums; without it the generic loop in
`Family` makes two virtual calls per observation per evaluation.
`BinomialFamily`, `NegBinFamily`, `OrdinalFamily`, `MultinomFamily`,
`ZeroInflatedFamily`, `BetaFamily` and `OrdBetaFamily` derived from `Family`
directly. The first of these is the path for a binomial with `augment = FALSE`
or a cloglog link, and `NegBinFamily` is the path for `negbin()` under the
default soft rules, where the augmentation is off. Each is a two-line change:
the base class and the constructor's initializer. Nothing else about them moved.

**B. One gate evaluator per pass** (`src/node.h`, `node.cpp`, `mcmc.cpp`;
*exact*). `Node::gate(i)` read `tree->X`, `tree->bandwidth`,
`tree->hypers->soft` and `tree->hypers->gate` through three pointers, and
indexed `X(i, var)` with a bounds check, once per observation in
`split_support()` and `make_base_children()`. `GateEval` fetches all of that
once per pass and holds the column pointer; the arithmetic and its order are
unchanged, so the weights are the same numbers. `Node::gate()` now goes through
it too, so there is one definition.

**C. A has-missing table for the categorical groups** (`src/node.h`,
`node.cpp`, `model.cpp`; *exact*). `draw_categorical_rule()` scanned the whole
column of `codes` for a negative entry on every proposal of a categorical rule,
to decide whether to draw a missing-value rule. For complete data -- the common
case -- that is a full pass over the sample per proposal, on a node whose own
support may be a hundredth of it. `bartisan_fit()` now computes one flag per
group and every `Tree` reads it. The random-number stream is untouched because
the extra draw is made under the same condition as before. A null table falls
back to the scan, which is what a `RandomTerm`'s private tree passes and never
reaches.

**D. Swap rather than re-split on a rejected change move** (`src/mcmc.cpp`,
`mcmc.h`; *exact*). `change_rule()` restored a rejected proposal with
`set_rule(rule_old); split_support()`, which re-evaluates every gate in the
node to recompute the exact vectors it had a moment earlier. The children's
`idx` and `wt` are now swapped into four `Context` buffers before the new rule
splits the support and swapped back on rejection, both in constant time. The
shared-topology change move (`shared_change()`, near line 2160) still
re-splits its lead and then mirrors O(n) copies to the twins; left alone, since
that path is opt-in and the copies dominate anyway.

**E. Scratch vectors for the categorical rule draw** (`node.h`, `node.cpp`;
*exact*). `available_levels()` and the open-level list were two heap
allocations per categorical proposal; they now live on the `Tree`.

**F. Nuisance updates for `Gamma()` and `negbin()`** (`src/family.cpp`,
`update_aux()` of each; *not exact by construction, exact as measured*). Both draw their dispersion by slice
sampling, and the slice target was a pass over the sample with three
log-gammas, two logs and an exponential per observation -- evaluated four to
eight times per sweep. For the gamma the log likelihood in the shape is a
function of four sums over the sample, so those are taken once and each
evaluation is arithmetic; for the negative binomial `lgamma(theta)` and
`theta log(theta)` come out of the loop as one term times the weight total,
`lgamma(y + 1)` is dropped because a slice sampler only ever differences its
target, and `exp(eta)` is taken once per update, leaving one log-gamma and one
log per observation. The same sum is formed in a different order, so the
target differs from the old one by rounding; the chain has the same stationary
distribution either way, and on these three seeds the draws came out identical
to the bit. Checked by every configuration's hash being identical to build
`new`, and by the suite. Measured at 4% on `Gamma("log")` and 1% on
`negbin()`.

**G. `ARMA_DONT_CHECK_CONFORMANCE`** (`src/Makevars`; *exact*). Armadillo
bounds-checks every `X(i, j)` and `y(i)` unless told not to, and those are in
every hot loop. The often-cited `ARMA_NO_DEBUG` does nothing in Armadillo 12+:
it undefines `ARMA_DEBUG`, which is never on. The working switch also drops the
size checks on whole-matrix operations, so a dimension mistake in the R-facing
entry points would become undefined behavior instead of an error.
Measured at 1.006 to 1.032 across the eleven configurations, the same sign on
every one, about 2% in total. **Not taken**: `src/Makevars` records the switch
and the number, and the check stays on. The safe way to most of the gain is
surgical -- `.at(i, j)` or a held column pointer at the sites that read an
observation by index inside a leaf loop, as `GateEval` (item B) already does
for `X` -- and the family unit functions' `y(i)` and `w(i)` are the large
remaining set, several hundred sites across `src/family.cpp` that a mechanical
edit could cover. Whether 2% is worth that edit is a judgment left open here.

## Larger items, reported and not changed

Ordered roughly by expected payoff. Each is self-contained enough to be picked
up cold.

1. **The non-augmented families are where the time goes.** At n = 1000 and 50
   trees, `binomial()` costs 0.9s with its augmentation and 6.5s without;
   `ordinal()` without it costs 13s, `multinomial()` 28s, `Beta()` 17s. The
   reason is structural, not a hot spot: a non-quadratic leaf refresh runs a
   Fisher-scoring loop from a fixed start (`fit_laplace1_at()`, `SCORE_INIT`,
   `SCORE_TOL = 1e-2`) and then one paired evaluation, so a leaf costs four to
   eight passes where a quadratic one costs one. The fixed start is what keeps
   the proposal a function of the state, so warm-starting is not available.
   Two routes that keep exactness: (i) for `TARGET_EXP_UP`/`EXP_DOWN` families
   the exponential form already collapses this to one pass under hard rules
   (`exponential_usable()` requires `!soft`); extending it to soft rules is
   blocked because each observation then carries its own exponent -- but a
   *quadrature* of the per-observation weights (binning `wt` into a few values)
   would recover a small number of exponentials at the cost of an approximation
   the Metropolis step corrects; (ii) the Beta family's digamma/trigamma tables
   (`BetaFamily::psi_lookup`) are the template for tabulating the ordinal and
   multinomial link functions, which are the transcendental cost in those two.

2. **Prediction walks both children of every hard split.** `eval_tree()` in
   `src/model.cpp` recurses into a zero-weight subtree because the record
   holds no subtree length and the cursor has to reach its end. Adding a
   seventh field to the record -- the record's own length in doubles, written
   by `encode_tree()` after the subtree is known -- lets `eval_tree()` skip a
   child with weight exactly zero, which halves the walk for hard rules and for
   every categorical and missing-value rule (whose gates are always 0 or 1).
   `RECORD_SIZE` is used in exactly two parsers, `eval_tree()` and
   `tree_splits_on()`, and nothing in R reads the record, so the change is
   contained; a stored fit from before it would not load, which before 0.1.0
   is free. `bartisan_predict()` is also called once per draw per tree per
   observation with `out[h](s, i) +=` bounds-checked; item G covers that.

3. **Seven more slice-sampled nuisance updates have the per-observation loop
   of item F.** `src/family.cpp` near lines 1840 (`AFTFamily` sigma), 1947 and
   2073 (the two augmented AFT families), 2655 (`ZeroInflatedFamily` theta),
   2874 (`BetaFamily` phi), 3216-3232 (`OrdBetaFamily` cutpoints and phi),
   3538-3556 (`TweedieFamily` phi and power), 3796 and 4021 (the augmented
   negative binomial and zero-inflated thetas). Same recipe: hoist what depends
   only on the parameter, precompute what depends only on eta, and where the
   likelihood is exponential-family in the parameter reduce to sufficient
   statistics. None of these is exact, so each wants the F-style check.

4. **The Polya-Gamma draw is a fifth to a third of an augmented-binomial
   sweep.** `LogitAugmentedFamily::before_forest()` draws `rpg(1, eta_i)` for
   every observation every sweep, and a `binomial()` sweep at n = 1000 is about
   1.5ms, of which the draws are 0.3-0.5ms. `polyagamma.cpp` is a faithful
   Devroye sampler; the only micro-target is `series_term()`'s
   `std::pow(x, 1.5)` (a multiply and a square root would do). A real change
   would be batching the draws or a table for the two `pnorm()` calls in
   `pinvgauss()` at the fixed truncation point, neither exact.

5. **A subclass of a `Concrete<X>` family cannot override a unit method and be
   heard.** `Concrete`'s loops call `self.Derived::logdens_unit()` with
   `Derived = X`, so an override in a further subclass is bypassed by the hot
   path while `Family::logdens()` still sees it -- two answers for one density.
   `DPMAFTFamily : DPMFamily` is the one such subclass and overrides only
   `update_aux()` and `reported_loglik()`, so it is not live. A `static_assert`
   is not available (the check needs the most-derived type), but a debug-build
   assertion in `Concrete`'s constructor that `typeid(*this) == typeid(Derived)`
   would catch the next one at the first fit rather than at the first wrong
   posterior.

6. **`update_bandwidth()` still saves and restores every support by copy**
   (`save_support()`/`restore_support()`, `src/node.cpp`). The copy is much
   cheaper than the gate evaluations it replaced, and the move runs every
   `bandwidth_every` sweeps, so this is small; the swap trick of item D would
   apply per node if it were ever worth it.

7. **`update_sigma_mu()` allocates a `std::vector` of every leaf value once per
   forest per sweep** (`src/mcmc.cpp`). Trivial cost; a `Context` buffer would
   remove it.

## Smaller observations, for a later pass

- `bartisan()` cannot be called through a `...`-forwarding wrapper when
  `offset`, `weights` or `subset` is among the forwarded arguments:
  `model.frame()` evaluates those in the formula's environment, where the
  forwarded `..5` does not resolve ("..5 used in an incorrect context"). `lm()`
  and `glm()` share the limitation, and `bcf()` avoids it with `do.call()`, so a
  user's wrapper needs the same; evaluating those three arguments before the
  frame call is built would remove it. Found by a test helper in this review,
  which was rewritten inline.
- `print.bartisan_fit()` reports "Draws: N kept after B warmup" without the
  thinning, and `summary()` reports draws without the chain count that
  `print()` gives. Neither is wrong; they could agree.
- `print.bartisan_fit()` validates `digits` only when the fit has nuisance
  parameters.
- `random_predict()` matches a new grouping level by `as.character()` against
  the fitted levels, which is right for factors and integers and would treat
  `1` and `1.0` as one level but `0.1 + 0.2` and `0.3` as two. Not a practical
  concern.
- The vignette's `avg_comparisons(fit)` table costs about 143s of a 4m 37s
  knit, all of it the fifteen-row table; a single-variable call is 8s. This was
  flagged when the section was written and is the user's call.
- `warn_runaway_scale()` fires on small `multinomial()` fits with the symmetric
  coding even without an offset (item 3 above). Whether the symmetric coding's
  unidentified common shift lets the leaf scale drift on short runs, or the
  fixture is simply too short, was not investigated.

## What a future agent should know about this file

The benchmark libraries lived in the session's scratch directory and are gone;
`_dev/review-bench.R` and the numbers above are what survives. To repeat a
comparison: install the two trees into two libraries, run the script once per
library with a label, and `merge()` the two results on `config` and `rep`; the
`hash_eta` and `hash_aux` columns must agree row for row before any timing is
read. Before *any* rebuild after a header edit, the `Makevars` rule now handles
the stale-object trap; before a rebuild after a *flag* edit, `touch` a header or
remove `src/*.o`, since make cannot see a changed flag.
