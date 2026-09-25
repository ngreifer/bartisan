# What is left before shipping, and the vignette series

Two questions, answered together because they turn out to be the same question:
sketching the workflow vignettes is what exposes which gaps are real. A gap that
no vignette needs to work around is not a gap worth holding the release for.

The standard applied throughout: **does this seriously affect a workflow the
package claims to support?** Speedups, polish and features that another package
already covers are excluded, however tempting.

The current state is the 2026-09-24 section, which is a CRAN submission
checklist. The 2026-09-15 and 2026-09-18 sections below it are kept as the
record of how the list got there.

## Where we stand (2026-09-15, after a clean check)

Everything the original list called blocking has landed except one item, and that
one has an argument for reclassifying rather than doing. The vignette series is
written. The counts are in the table below rather than here, so that there is
one place to correct them.

`R CMD check` is **`Status: OK`** on the built tarball, with no ERROR, WARNING or
NOTE. It was 1 ERROR, 1 WARNING and 1 NOTE a round earlier, which was the
finding then: the suite and the knit had not been standing in for a check. All
three are fixed and the fixes are confirmed against a real build, which for one
of them is the only demonstration available, the finding having been a test that
could only fail under `R CMD check`. Part 4 records what they were.

| | |
| --- | --- |
| Blocking, done | the rename, `variable_importance()`, `as_draws(eta = )`, the `predict(type = "density")` warning |
| Blocking, open | `custom_family()` has no posterior predictive draws (item 4); already documented elsewhere, needs one cross-reference to stop being a blocker |
| Landed since, unplanned | `estimate_effect()`, `bcf()`, `diagnose()`, `partial_dependence()`, `kfold()`, `prior_only`, `prior_summary()`, `ranef()` |
| Check findings, fixed and verified | the invariants test asserts against `getNamespaceInfo()` rather than a path only a checkout has; `potential_outcomes` and `digits` documented; `^vignettes/figure$` added to `.Rbuildignore` (Part 4) |
| Declined | `predictive_interval()` and `predictive_error()`, and with them the cross-reference in `vignette("bartisan")` that promised a section `vignette("effects")` never had (Part 3) |
| Still open | `VarCorr()`, a `NEWS.md`, and a version that is not `0.0.0.9000` |
| Verification | suite 39 files / 2300 assertions, ten vignettes knitting with no chunk errors, and `R CMD check` `Status: OK` |

The judgment: **the code is submittable and the paperwork is not.** Part 3 is the
survey that says the first, done by reading an export list rather than guessing
at parity, and Part 4 is the record of what the check found before it went
green. What is left is neither a feature nor a fix: a `NEWS.md`, a version that
is not `0.0.0.9000`, and one cross-reference on the `custom_family()` page.
"What to do next" at the bottom has them in order.

## Re-read 2026-09-18: what is left, item by item

The section above is the read from 2026-09-15. Every open item in it was checked
again against the tree rather than against its own description, which is the
only way this list stays worth anything. Nothing that was closed has reopened.
What follows supersedes the table above where the two disagree.

| Item | State on 2026-09-18 |
| --- | --- |
| Version is `0.0.0.9000` | **Open.** `DESCRIPTION` unchanged. This is the one thing `--as-cran` actually objects to. |
| No `NEWS.md` | **Closed 2026-09-18.** Added, with an initial-submission entry. |
| `custom_family()` caveat not on its own page | **Closed 2026-09-18.** It has its own page now, `?custom_family`, with a section saying what a log density cannot supply and noting that `loo()` and `waic()` are unaffected. |
| `vignette("bartisan")` gaps | **Closed 2026-09-18.** One paragraph on the decision rules, pointing at `vignette("implementation")`. The augmentations are deliberately left out; they are covered in `vignette("families")` and are not a decision the tour asks a reader to make. |
| `--as-cran` never completed | **Open, and the reason to care has not changed.** It hangs at `checking use of S3 registration` on this machine, which looks like a network wait rather than the package. Needs one run somewhere unrestricted. |
| A second platform, and the spell check | **Open.** win-builder, macbuilder and R-devel have still never seen this package. |
| Two `tweedie()` references | **Closed 2026-09-18.** `jorgensen1987` and `dunn2005`, verified by DOI against Crossref and cited in `vignette("families")` and `vignette("implementation")`. |
| A hex logo | **Open, newly listed.** None exists; `_pkgdown.yml` has no `logo` key and `README.Rmd` has no badge row. |
| `VarCorr()` | Post-1.0, unchanged. |

### Three findings that are new since 2026-09-15

**The counts in this file and in `TASKS.md` were both stale, in the same
direction.** The suite is now 477 tests and 2430 expectations over 39 files with
`NOT_CRAN=true`, 0 failures and 0 skips; inside `R CMD check` it is 1898 passing
and **91** skipped. This file said "39 files / 2300 assertions" and `TASKS.md`
said "463 tests and 2312 expectations ... skips 67 and passes 1560". Both are
corrected.

The 91 matter more than the drift does. They are every `skip_on_cran()`, and
they include all four ordinal tests, all four `prior_only` tests and nine of
`test-augment.R` -- precisely the code that changed most recently. `Status: OK`
from a plain check is therefore not evidence that the induced-Dirichlet prior,
the empty-category handling or the `augment` flag work; the `NOT_CRAN=true` run
is the evidence, and it is clean.

**`ordinal("probit")` with `augment = FALSE` stalled once.** *(Set aside
2026-09-18: seen a single time and not reproduced since, so it is recorded
rather than tracked.)* Over 15
replicates at n = 400, fourteen un-augmented hard-rule fits took 30 to 43
seconds and one took **3882**; under soft rules, thirteen took 51 to 65 seconds
and two took about 600. Every augmented fit and both `ordinal("logit")` arms are
tight to a second, so it is specific to the un-augmented cutpoint sampler. That
is a supported, documented setting, which is what makes this the one genuinely
new candidate blocker rather than a curiosity. Under investigation in its own
session; the measurements are in `TASKS.md` under "stalls".

**Nothing is committed.** *(resolved 2026-09-18 in `ff12510`.)* The tree is
clean. The last full check still predates the documentation changes in that
commit, so one more run is owed before a submission tarball is built.

### What that adds up to

The 2026-09-15 judgment holds and gets one qualification: **the code is
submittable, the paperwork is not, and one measurement now wants explaining.**
The paperwork is unchanged in kind -- a version, a `NEWS.md`, one
cross-reference, a logo -- and none of it is work. The qualification is the
probit stall, which should be understood before a release rather than after,
because `augment = FALSE` is a setting the documentation tells people to try
when a fit's diagnostics look poor.

Two pieces of verification are owed that no amount of local checking supplies:
`--as-cran` on an unrestricted network, and one other platform. Both were owed
on 2026-09-15 too.

## Re-read 2026-09-24: the CRAN submission checklist

Features have landed since 2026-09-18, among them `vignette("varying")`, the `.`
in `vc()`, and function-valued `values` in `partial_dependence()`, and none of
them reopened an item. What has moved is that the list was checked against
CRAN's own requirements for a first submission rather than against the
package's, using the `cran-extrachecks` checklist, and that turned up items no
earlier section had. This section supersedes both earlier ones where they
disagree.

The last full check is pueue task 207 on 2026-09-23: `_dev/check.sh
--no-build-vignettes` with `NOT_CRAN=true`, 27 minutes, **`Status: 1 WARNING,
1 NOTE`**, 2566 expectations passing and none failing. Both findings are
paperwork and one is now fixed; see the table.

| Item | State on 2026-09-24 |
| --- | --- |
| WARNING: GNU extension in `src/Makevars` | **Fixed.** `$(OBJECTS): $(wildcard *.h)` (from `6de73cc`) is GNU make only. The dependency is load-bearing in development, since it is what stops a stale object linking against a changed header, and irrelevant on CRAN, which builds from a clean tarball. So the rule stays and the headers are now listed by name. Declaring `SystemRequirements: GNU make` was the alternative, and CRAN's incoming check turns that into a NOTE of its own. Verified by applying R CMD check's own pattern to the file (no hits) and by touching `utils.h`, after which all 11 objects rebuilt. A new header has to be added to the list by hand. |
| NOTE: "No news entries found" | **Open.** `NEWS.md`'s heading is `# bartisan (development version)`, which R cannot parse as a version. Goes with the version bump. |
| Version is `0.0.0.9000` | **Open.** Wants a release number such as `0.1.0`, with `# bartisan 0.1.0` as the NEWS heading. |
| `Authors@R` has no `cph` for the maintainer | **Fixed.** Now `c("aut", "cre", "cph")`. Linero stays `c("ctb", "cph")`, and his comment now reads "Author of FlexBart, the reference implementation of Linero (2025), from which the MCMC engine is adapted", so that it cannot be read as Deshpande's *flexBART*. |
| `Description:` quoted `'glm()'` | **Fixed.** CRAN quotes software names and leaves function names bare. The spell check will still flag Tweedie, SoftBart and Linero, which `cran-comments.md` should explain. |
| Licensing of the adapted engine | **Open, and a question for Linero.** Provenance was checked on 2026-09-24: the engine is adapted from `FlexBart`, the R package in the reproduction materials of Linero (2025) (in `_dev/Reproduce/packages-scripts/packages/`), which is a different package from Deshpande's *flexBART* on CRAN. The evidence is in the code: `src/mcmc.cpp` names two places where it departs from "Linero's code" (the Fisher-scoring tolerance and the birth probability after a root collapse), and the slice sampler keeps FlexBart's structure and variable names. FlexBart's `DESCRIPTION` says `License: GPL 2.0` and it ships no LICENSE file. If that means version 2 only, a derivative cannot be offered under `GPL (>= 2)`, since that would allow GPL-3. Either Linero confirms "or later" in writing, or *bartisan* is licensed `GPL-2`. |
| `cran-comments.md` | **Open.** Does not exist. First submission, check results, and the three spell-check words. |
| README install instructions | **Open.** Says the package is not on CRAN and installs with `pak::pak("ngreifer/bartisan")`. Wants `install.packages("bartisan")`, edited in `README.Rmd` and rebuilt. |
| Vignette build time | **Reduced, still over budget.** Timed on 2026-09-24 with `_dev/vignette-timing.R` (a fresh R process per vignette, against the installed package, with `_R_CHECK_LIMIT_CORES_=TRUE`): **9.2 minutes** for all eleven, every one knitting cleanly. Moving `bartisan`, `causal` and `importance` from four chains to one, with a pointer to `vignette("diagnostics")` for the fit controls, took them from 238, 139 and 22 seconds to 56, 39 and 9. The prediction was under 8 minutes and it was wrong, because the two largest were untouched: `diagnostics` at 140 seconds keeps its four chains, and its `update()` to 2000 + 8000 draws is ten times the default fit's iterations; `comparison` at 128 seconds already used one chain and makes about twenty fits, ten of them inside two five-fold `kfold()` calls. Then `effects` 67, `families` 63, `varying` 47. With compilation, tests and examples on top, the check is well past ten minutes. What is left is a choice between precomputing the heaviest chunks (as `survival` already does), skipping evaluation on CRAN through the `run` switch most of the vignettes already have while shipping HTML built locally, or cutting content. Per-chunk times are in `_dev/vignette-chunk-times.md` (from `_dev/vignette-chunk-timing.R`): 70 of 176 chunks hold 527 of the 547 seconds, and the ten largest about half, led by diagnostics' `mixfixed` (69), families' `zerocompare` (56), comparison's `links` (34), bartisan's `comparisons` (30) and comparison's `survfits` (29). Diagnostics at 2 chains instead of 4 renders in 74 seconds against 143; causal at 2 chains, with the multisession plan restored, in 45 against 39, since the chains run in parallel on CRAN's two cores. Re-timed later on 2026-09-24 after the cuts (families' four-way comparison, comparison's cloglog fit and point-mass section, the faster prior-only fit) and the switch of causal's lalonde fit to `tweedie()`: an estimated **7.6 minutes** on a quiet machine, from a run taken under load and corrected by the slowdown of unchanged chunks (median 1.44). The largest chunks left are diagnostics' `mixfixed` (about 74 s), bartisan's `comparisons` (30), diagnostics' `diagnose_effect` (28), comparison's `survfits` (17), and effects' three two-predictor partial dependence chunks (about 13 each). |
| `--as-cran` never completed | **Open**, unchanged since 2026-09-15. Hangs on a network wait in this sandbox. Needs one run on an unrestricted network; it also covers the URL check. |
| A second platform | **Open**, unchanged. win-builder, the macOS builder and R-devel have never seen the package. All three upload it to an outside service, so they wait on an explicit go-ahead. |
| Tests under CRAN's settings | **Unmeasured.** The 27 minutes are with `NOT_CRAN=true`. Without it the 91 `skip_on_cran()` tests drop out, and how long the rest take is not known. |
| A hex logo | Optional, unchanged. Cheaper before the README and site are final. |
| `inst/CITATION` | Optional, newly listed. |

One side effect of the single chains is worth recording. In `vignette("causal")`, the `lalonde`
fit, `bcf(family = dpm())`, gave average potential outcomes among the treated of 5420 to 5650 for
`Y[0]` in two four-chain builds and 3780 with one chain, while the treated group's observed mean
earnings are 6349. The ATT intervals overlap (183 [-248, 823] and 228 [-222, 971] against 370
[-151, 1090]), so the contrast is not the problem; the level of a `dpm()` fit is, and it evidently
mixes slowly enough that one chain lands somewhere else. Not investigated yet.

Checked and not a problem: every exported help page has `\value` and examples,
there is no `\dontrun{}` and no commented-out example code, there are no `http://`
links and no relative links in the README, the Title is title case and under 65
characters, the Description does not open with a banned phrase and cites its
methods by DOI, the installed package is 4.6 MB and the tarball 2.4 MB, and the
package's own C++ uses no OpenMP.

### What to do next, in order

1. Settle the license with Linero, since the answer may change `DESCRIPTION`.
2. The paperwork: version bump and NEWS heading, `cran-comments.md`, README
   install line.
3. Decide how the vignettes fit CRAN's time budget (see the table), then re-time
   them with `_dev/vignette-timing.R`.
4. A full check, which is owed again after the Makevars change.
5. `--as-cran` somewhere with an unrestricted network, then win-builder,
   the macOS builder and R-devel.
6. Submit.

## Part 1: what is actually missing

### Blocking

**1. The package name.** *(done -- the rename to `bartisan` landed in `cdc3278`)*
What collided, case-insensitively, was the old name `genbart` against the
archived CRAN package `genBart`. `bartisan` collides with nothing on the current
index or among the archived names. The shortlist that was drawn up is kept in
`TASKS.md` as a record.

**2. Variable importance has no accessor.** *(done -- `variable_importance()`)* This is the largest genuine gap
against the stated workflow. The information exists and is correct --
`summary()` prints a "Predictor usage" block, and `fit$counts$eta` is a
draws-by-predictor matrix of split counts -- but nothing is exported. Every user
doing variable selection will reinvent

```r
colMeans(fit$counts$eta > 0)          # posterior probability of being used
colMeans(fit$counts$eta)              # mean splits per draw
```

and the second is the wrong one to reach for first. With `sparsity = TRUE` the
first separates cleanly (measured: 1.00, 1.00 for two real predictors against
0.005, 0.010, 0.095 for three noise ones), which is a genuinely useful selection
procedure that is currently undiscoverable.

*Needs:* one exported accessor returning a tidy data frame -- predictor, mean
splits, posterior probability of use, and the DART inclusion probability when
`sparsity = TRUE` -- plus documentation of which column answers which question.
Small work, high leverage, and vignette 3 cannot be written without it.

**3. `as_draws()` omits the fitted values.** *(done -- `eta` argument)* It exposes six variables: `loglik`,
`sigma_mu.eta`, and the nuisance parameters. Not `eta`. Convergence on `eta` is
*already diagnosed* -- `fit$rhat` carries `eta.eta (worst over observations)` --
so the package knows this is the quantity that matters and then declines to hand
it to `bayesplot`. A user wanting a trace plot of a fitted value has to build the
draws array themselves from `posterior_epred()`.

*Needs:* a `variables` argument on `as_draws.bartisan()`, or including a small
representative set of `eta` columns by default. Small work; vignette 2 needs it.

**4. `custom_family()` has no posterior predictive draws.** *(open, and the only
one still open)* A log density supplies
no way to draw from it, so `simulate()`, `pp_check()`, `r2()` and everything
built on them are unavailable for a user-written likelihood. Already in `TASKS.md`.
Medium work: an optional `rng` argument alongside the density.

`TASKS.md` argues the lighter reading, that `custom_family()` is the escape hatch
rather than a supported path and every compiled family has a sampler. That
reading is the right one *provided the limitation is documented where it is met
rather than discovered*, and it is nearly there already: `?bartisan-interop` and
`vignette("faq")` both state it plainly, while `?bartisan-families`, the page
that documents `custom_family()` itself, does not. One cross-reference settles
it, and the `rng` argument is then a post-1.0 feature rather than a blocker.

**5. `predict(type = "density")` returns NaN silently** *(done -- it warns; the value stays NaN)* when a composed link's
inverse sends the predictor outside the family's support. Already in `TASKS.md`,
deliberately left alone because it changes `predict()`'s output contract. It is a
silent-wrong-answer class of bug and should be settled before release rather
than after -- a warning is enough.

### Not blocking, and worth saying so explicitly

These are real work, but none of them blocks a workflow the package claims:

- **Partial dependence** -- `marginaleffects::plot_predictions()` covers it.
- **Interactions** -- `comparisons(by = )` covers the useful cases; Friedman's H
  is a dozen lines in a vignette.
- **Model selection** -- `loo()` and `waic()` are present and work.
- **Mixing diagnostics** -- `rhat`, `ess_bulk`, `ess_tail` are present and
  already cover `eta`; only the plotting hand-off is missing (item 3).
- **A separate treatment forest (BCF)** -- *this was called the right post-1.0
  feature and then landed anyway, as `bcf()`.* The reasoning held right up to the
  point where writing `vignette("causal")` made the caveat about
  regularization-induced confounding tiresome to keep making, which is the usual
  sign that the feature is cheaper than the paragraph.
- Correlated random effects across predictors, `quasi()`, the joint ordinal
  cutpoint update, grow-from-root warm start, categorical splits on level
  subsets, soft random tree features -- all post-1.0.

### Deliberately not on the list

Further speedups. `weibull_aft()` at ~8 seconds for 700 observations is the
slowest family and could plausibly reach ~1 second with a Frühwirth-Schnatter
mixture-of-normals augmentation of its Gumbel error -- but that augmentation is
an *approximation*, and adopting it would weaken the exactness claim the package
currently makes. That is a design decision, not a bug, and not a release blocker
now that `weibull_aft()` is no longer the default.

### Measured after the fact: what the other BART packages export

Surveyed `dbarts`, `BART`, `bartMachine` and `SoftBart` by their actual export
lists rather than from memory. Every recurring non-fitting feature, and where it
lives here:

| They export | Here |
| --- | --- |
| partial dependence (`pdbart`, `pd_plot`, `pdsoftbart`) | `marginaleffects::plot_predictions()` |
| variable importance (`investigate_var_importance`, `get_var_props_over_chain`, `posterior_probs`) | `variable_importance()` |
| variable selection (`var_selection_by_permute`, `mc.wbart.gse`) | `prop_used` with `sparsity = TRUE` -- a posterior probability rather than a permutation null |
| credible intervals (`calc_credible_intervals`) | `marginaleffects::predictions()` |
| prediction intervals (`calc_prediction_intervals`) | `posterior_predict()` |
| convergence plots (`plot_convergence_diagnostics`) | `fit$rhat`, `as_draws()` + bayesplot |
| fit checks (`plot_y_vs_yhat`, `check_bart_error_assumptions`) | `pp_check()`, `residuals()` |
| interaction detection (`interaction_investigator`) | `marginaleffects::comparisons(by = )`, partially |
| k-fold CV (`k_fold_cv`, `xbart`, `bartMachineCV`) | `kfold()`, added since this table was written |
| model matrix helpers (`makeModelMatrixFromDataFrame`, `dummify_data`, `preprocess_df`, `bartModelMatrix`) | **not needed** -- the formula interface does it |

Two conclusions.

**There is no remaining feature gap.** After `variable_importance()`, everything
those four packages export is either present, covered by a better-maintained
general package, or unnecessary here. The one absence when this was written was
k-fold cross-validation, argued away on the grounds that `xbart` and
`bartMachineCV` exist to tune `k`, `q` and the tree count, which the defaults
here are meant to settle. That argument was half right: tuning is not what a
`kfold()` here is for, but *checking `loo()`* is, since the two estimate the same
quantity and the Pareto diagnostics are the only thing otherwise vouching for
the approximation. `kfold()` exists now and `vignette("comparison")` uses it
that way.

**What is left is discoverability, not function.** A newcomer cannot guess that
intervals come from `marginaleffects::predictions()`, or that a prediction
interval needs `posterior_predict()`. Four of those packages ship one function
per task with an obvious name; this one ships a map. The map has to be good --
hence the rewritten `?bartisan-package` and the vignette series below.

## Part 3: measured against rstanarm

The question that prompted this was whether *bartisan* could stand in for
*rstanarm* as the thing reached for by default. Surveyed the same way as the BART
packages above, by reading `getNamespaceExports("rstanarm")` rather than from
memory, and restricted to what applies to a single-response regression: the
joint longitudinal-survival machinery (`posterior_survfit`, `posterior_traj`,
`ps_check`, `stanjm_*`, `stanmvreg_*`) is a different model class, `plot_nonlinear`
is for GAM terms, `se` and `posterior_interval` are about coefficients a forest
does not have, and `pairs_*` is a geometry diagnostic for a Hamiltonian sampler.

| rstanarm | Here |
| --- | --- |
| `posterior_predict`, `posterior_epred`, `posterior_linpred`, `log_lik` | present, same names |
| `loo`, `waic`, `loo_compare`, `kfold` | present; `kfold()` was the gap and is not now |
| `loo_model_weights` | works already on our `<loo>` objects (checked: 0.889 / 0.111 on two nested fits) |
| `pp_check` | present, and dispatches the whole *bayesplot* set including the `loo_*` checks |
| `bayes_R2` | `performance::r2()` and `r2_posterior()` |
| `prior_summary` | present, and says more than rstanarm's can: the prior here is a function of the response and the tree count |
| `sigma`, `nsamples` | `sigma()`; a draw count is `nrow(fit$sigma_mu)` and wants no accessor |
| `get_y`, `get_x`, `get_z` | `insight::get_data()`, `fit$model` |
| `Surv`, `invlogit`, `logit` | re-exports and one-liners; use the originals |
| **`ranef`, `VarCorr`, `fixef`, `ngrps`** | **absent.** See below |
| `predictive_interval`, `predictive_error` | absent, and deliberately so; both are one line over `posterior_predict()` |
| `loo_predict`, `loo_linpred`, `loo_predictive_interval`, `loo_R2` | absent; the PSIS machinery is there, the four wrappers are not |
| `posterior_vs_prior` | absent, and newly *possible*: `prior_only = TRUE` gives the other half |
| `pp_validate` | absent as an export; the simulation-based calibration it does is in `_dev/sbc.R` and was run |
| `launch_shinystan` | out of scope; wants a shinystan object |

### The one that matters: no `ranef()`

The package advertises `(1 | group)` in the formula, fits it, and stores the
result correctly: `fit$ranef$eta` is a draws-by-levels matrix with columns named
`g:a`, `g:b` and so on, and `fit$tau$eta` is the draws of each term's scale.
`?bartisan` says as much, in as many words ("The intercepts are in `fit$ranef`
and their standard deviations in `fit$tau`"), so this is a weaker version of item
2 than that item was: the information is documented as well as correct. What is
missing is only the accessor. `lme4::ranef()` and `nlme::VarCorr()` both fail
with "no applicable method", and those are the two names a mixed-model user tries
before reading anything.

It is still worth doing before submission, on the same argument that made
`variable_importance()` blocking and for the same reason the map in Part 1's
closing note has to be good: a package that ships a map cannot also ship the one
turning the reader is sure to take without checking. And it is small. Two
methods over matrices that already exist, and `tau` gives the scale for free.

The counter-argument, which is not nothing: a `$`-accessible matrix that the help
page names is a perfectly usable interface, and `ranef()` would be sugar. If the
release is being held for this alone, ship without it.

`fixef` and `ngrps` do not follow. A forest has no fixed effects to return, and a
group count is `length(levels(...))`.

### The rest is convenience, and one is newly cheap

**`predictive_interval()` and `predictive_error()` are not being added**, decided
2026-09-15. They are a quantile and a subtraction over `posterior_predict()`,
which is exported and documented, and nothing here is unreachable without them.
Recorded so the question is not reopened as though it had never been asked.

Looking for somewhere they were being done by hand turned up a dangling
cross-reference rather than a workaround: `vignette("bartisan")` told the reader
that the distinction between an interval on the mean and one on a new observation
"matters a great deal" and sent them to `vignette("effects")` for it, which has
no section on it and never uses the words. With the wrappers declined, the
pointer went instead, which is the half of that choice that costs nothing.

`posterior_vs_prior()` is the interesting one. It was impossible here until
`prior_only` landed and is now two fits and a comparison, which is what
`vignette("diagnostics")` walks through: the prior's fitted risk against the
posterior's over the same quantiles. Whether that deserves a function or stays a
documented recipe is a discoverability judgment rather than a capability one.

### What this does not answer

Parity of *exports* is not parity of *purpose*. Someone reaching for
`stan_glm()` wants a coefficient with an interval on it, and this package will
never return one: that is what the model is. The honest claim is narrower than
"replaces rstanarm" -- it is that for the questions a forest can answer
(prediction, effects, comparison, calibration), there is no longer a workflow
step that rstanarm supports and this does not, once `ranef()` exists.

## Part 2: the vignette series (written)

All six landed, and two more besides. The names moved as they went, so the
sketches below are keyed to the files that came of them:

| Sketched as | Shipped as | Title |
| --- | --- | --- |
| 1 `workflow` | `bartisan.Rmd` | Getting started with bartisan |
| 2 `diagnostics` | `diagnostics.Rmd` | Checking convergence and fit |
| 3 `importance` | `importance.Rmd` | Which variables matter |
| 4 `effects` | `effects.Rmd` | Effects, curves, and interactions |
| 5 `comparison` | `comparison.Rmd` | Choosing between models |
| 6 `causal` | `causal.Rmd` | Causal Inference with BART |
| -- | `implementation.Rmd` | Generalized BART with bartisan |
| -- | `faq.Rmd` | bartisan Frequently Asked Questions |
| reference | `families.Rmd`, `survival.Rmd` | unchanged in role |

The one structural change: `bartisan.Rmd` was the method document and is now the
orientation piece, the method having moved to `implementation.Rmd`. That is the
right way round. A reader who types `vignette("bartisan")` wants to fit
something, not to read about a reversible-jump sampler.

`faq.Rmd` was not planned and came out of the questions that kept recurring
while the other six were written, which is a better way to arrive at an FAQ than
guessing at one.

Everything below is the original sketch, kept as the record of what was intended
against what shipped. The series was organized by **the analyst's workflow** --
the arc a person actually walks from a data frame to a defensible claim.

A design rule for all of them: they must be honest about the places where the
BART answer is *worse* or *harder* than the `lm()` answer, not only where it is
better. A reader who finds out later that BART's variable importance is unstable
under correlated predictors, having not been told, will not trust the rest.

---

### 1. `workflow` (done) -- "A regression workflow with BART"

**The orientation piece, and the one to write first.** One dataset, one question,
start to finish, with every step deliberately shallow and a pointer to the
vignette that goes deep. Someone who reads only this should be able to do a
competent analysis.

The arc: fit → did it converge → does it fit → what matters → what is the effect
→ is this model better than that one → report it.

Deliberately establishes the two habits the rest of the series depends on: **read
counterfactual estimands through `marginaleffects`, never coefficients** (BART
has none to read), and **check convergence on the fitted values, not just on
`sigma`**.

Ends with the comparison to `lm()`/`glm()` stated plainly: what you give up is
the coefficient table and the ability to state the functional form; what you get
is not having to state the functional form. Nothing about a BART fit is a
substitute for knowing what question you are asking.

### 2. `diagnostics` (done) -- "Has it converged, and does it fit?"

Two questions usually run together and worth separating.

*Converged*: what R-hat and ESS mean when the parameter is a function rather than
a number, why `eta`'s worst-over-observations R-hat is the one to watch, running
`chains > 1` and why the default of 1 is not an endorsement, trace and rank plots
through `posterior` and `bayesplot`, and what to do when it has not converged
(more draws, more trees, harder rules).

*Fits*: `pp_check()`, residuals for the families that have meaningful ones, and
the specific failure BART shows when it is under-fitting -- shrinkage toward the
mean at the edges of the predictor space, which looks like bias and is.

Also the honest part: a BART fit can have excellent R-hat on `sigma` and a badly
mixed forest, and the tree structure itself is not identified, so "did the chains
find the same trees" is not a question worth asking. Only functionals of the fit
are.

*Needs item 3 above.*

### 3. `importance` (done) -- "Which variables matter?"

Split counts, `prop_used`, and the DART posterior inclusion probability, with the
distinction between them made sharp: **how often a variable is split on is not
how much it matters**, and neither is a test.

Turning `sparsity = TRUE` into an actual selection procedure, including the
measured separation. Then the three cautions that make this vignette worth
writing rather than a paragraph:

- Under **correlated predictors** the split counts distribute arbitrarily among
  the correlated set; a variable can be genuinely important and rarely split on
  because a collinear partner absorbed it.
- Importance is **not an effect**: a variable can be split on constantly and move
  the prediction very little. `marginaleffects` answers "how much does it move
  the prediction", and that is usually the real question.
- Importance is **not causal**, and the vignette should say so before someone
  reads a ranking as a set of causes.

*Needs item 2 above.*

### 4. `effects` (done) -- "Reading the fitted function"

The vignette that replaces the coefficient table. Predictions, comparisons and
slopes through `marginaleffects`; partial dependence via `plot_predictions()`;
average versus conditional effects and why the average is usually the one to
report.

Then interactions, which is the part with no ready-made tool: `comparisons(by =)`
to get an effect at levels of a moderator, and how to tell a real interaction
from the shrinkage artifact that produces a similar picture -- the honest answer
being that the posterior interval on the *difference of differences* is the test,
not the eyeball.

Nonlinearity: how to show a fitted curve with an interval, and why a soft-rule
fit is the one to show.

### 5. `comparison` (done) -- "Choosing between models"

`loo()` and `waic()`, what they compare and what they cannot. The measure trap
from the survival vignette generalizes: **a log score is only comparable between
models that put their density on the same scale**, and the package has at least
one pair that does not.

What model selection *means* here, which differs from the `lm()` case and is the
reason this is a vignette rather than a paragraph: you are not selecting a mean
structure, because the forest selects that. You are selecting the **likelihood** --
the family, the link, whether the variance is constant, whether the error shape
is assumed. That reframing is the content.

Also: when `loo()` breaks and what it means (the `ph()` bin-count result is the
worked example), and why comparing a BART fit to a linear model by `loo` is a
fair and useful thing to do.

### 6. `causal` (done) -- "Causal inference with BART"

The one where being clear about what the package does *not* do matters most.

G-computation through `avg_comparisons()`: the ATE, the ATT, and effects on the
scale the question is asked in. Overlap and positivity as things to check before
the estimate, not after. Weighting through `WeightIt` when the design calls for
it, and the propensity score as a covariate (Hahn, Murray and Carvalho).

Then the caveat that has to be in this vignette rather than a footnote:
**regularization-induced confounding**. BART's prior shrinks the treatment effect
toward zero along with everything else, and when treatment assignment is strongly
predicted by the covariates that shrinkage biases the effect. BCF's separate
treatment forest is the standard answer and `bartisan` does not have one. Say what
the bias looks like, when it bites, and what to do meanwhile.

Uncertainty: the posterior interval on an average effect is a credible interval
for the estimand under the model, and is not a confidence interval for a causal
effect unless the identification assumptions hold. Worth one paragraph, plainly.

---

### Sequencing (as planned, and what actually happened)

1 first: it is the entry point and the others are its expansions. Then 4, which
is the highest-value single document (most people's real question is "what is the
effect") and needs no new package code. Then 3 and 2, each after its blocking
item lands. Then 5. Then 6, last, because it is the one where being wrong is most
costly and it benefits from the other five being settled.

`families` and `survival` stay as reference documents and get cross-links from 1
and 5. `bartisan` stays as the "how it works" document; it is the only one
organized around the method rather than the workflow, which is correct for it.

That last paragraph is the one the plan got wrong, and it is worth saying why.
Keeping the method document at `vignette("bartisan")` meant the name a reader
guesses first pointed at the hardest document in the set. Renaming it to
`implementation` cost nothing and fixed the entry point.

Otherwise the order held. What it did not predict is how much package code the
*writing* would turn up: `estimate_effect()`, `diagnose()`,
`partial_dependence()`, `bcf()`, `kfold()`, `prior_only` and `prior_summary()`
all came out of a vignette needing something and finding it absent, which was the
premise of doing the two halves together and is the strongest evidence for it.

## Part 4: what `R CMD check` says

Run with `_dev/check.sh` on 2026-09-14, against the built tarball rather than the
source tree. `Status: 1 ERROR, 1 WARNING, 1 NOTE`. None of the three had shown up
in `testthat::test_local()` or in a knit, because none of them can. **All three
are fixed and the check is `Status: OK` as of 2026-09-15**; kept here because
what each one was is more useful than the fact that it is gone, and because the
first is the standing argument for running the check at all.

### ERROR: a test that only passes in a checkout

```
── Error ('test-invariants.R:304:3'): every S3 method on a foreign generic is registered ──
Error in `file(con, "r")`: cannot open the connection
 1. └─base::readLines("../../NAMESPACE", warn = FALSE)
```

The test is there for a good reason, and its own comment gives it: roxygen
attaches `@exportS3Method` positionally, so inserting a helper between the tag
and the function moves the registration onto the helper, and `load_all()` hides
that by registering from source regardless. So it asserts against `NAMESPACE`
itself.

It reads `../../NAMESPACE`, which resolves in `tests/testthat/` and does not
resolve under `R CMD check`, where the tests run from `bartisan.Rcheck/tests/`.
The test written specifically to catch what only breaks in an installed package
is the one that breaks against the installed package.

The fix is not a path. It is to assert against the registry rather than the file:
`getNamespaceInfo("bartisan", "S3methods")` is what actually governs dispatch,
exists in both contexts, and is the thing the test is really about. `NAMESPACE`
was only ever a proxy for it.

### WARNING: two undocumented arguments

```
Undocumented arguments in Rd file 'estimate_effect.Rd'
  'potential_outcomes'
Undocumented arguments in Rd file 'partial_dependence.Rd'
  'digits'
```

Both are real and both are one `@param` each. `digits` is on
`print.bartisan_partial()`, which shares a page with `partial_dependence()` and
never had the tag; `potential_outcomes` is an `estimate_effect()` argument.

### NOTE: a leftover knitr directory

```
The following directory looks like a leftover from 'knitr':
  'figure'
```

This is `vignettes/figure`. `.Rbuildignore` has `^figure$`, which excludes the
root-level one and not this one, and both exist on disk from knitting vignettes
outside a build. `.gitignore` covers both, so neither is committed and the
tarball picked one up anyway, which is the distinction between the two files
worth remembering: `.gitignore` keeps it out of the repository and
`.Rbuildignore` keeps it out of the package. Needs `^vignettes/figure$`.

## What to do next, in order

0. ~~**The three check findings.**~~ Done. Part 4 says what they were.
1. ~~**`ranef()` and `VarCorr()` methods.**~~ `ranef()` is done, registered on
   `nlme::ranef`, which *lme4* re-exports as the same function object so
   both qualifications reach it. `VarCorr()` is not, and is post-1.0: `tau` is
   the only thing it would return and `?bartisan` names it.
2. ~~**Move the `custom_family()` limitation onto its own help page.**~~ Done in
   `ff12510`, and more than a cross-reference in the end: `custom_family()` now
   has its own page rather than a shared one. Original note follows. It is
   already stated on `?bartisan-interop` ("no posterior predictive distribution
   at all, because a log density supplies no way to draw from it") and in
   `vignette("faq")`, so this is one cross-reference rather than new prose. But
   `?bartisan-families` is where `custom_family()` is documented and where a
   reader meets it, and the caveat is not there. That placement is what makes
   item 4 non-blocking rather than merely unfinished.
3. ~~**`NEWS.md` does not exist.**~~ Added. The version is still `0.0.0.9000`,
   which remains item 0 of any submission.
4. ~~**Re-run `_dev/check.sh` to a clean bill.**~~ Done: `Status: OK`. Worth
   re-running before submission itself, since the last three findings all
   arrived from a source the suite and the knit cannot see.
5. ~~**The `ordinal("probit")` stall under `augment = FALSE`.**~~ Set aside:
   seen once, never reproduced.
6. **A hex logo**, which wires into `_pkgdown.yml` and `README.Rmd` and so is
   cheaper before the site and README are final than after.
7. ~~**Commit.**~~ Done in `ff12510`.
8. **Re-run the full check**, which now predates a commit's worth of
   documentation changes. This is item 4 again, and it is the one that keeps
   coming back because the suite and the knit cannot see what it sees.
9. Post-1.0, in no order: `VarCorr()`, the `loo_*` prediction wrappers,
   `posterior_vs_prior()`, an `rng` for `custom_family()`, and the DART
   inclusion probability as a stored quantity (`TASKS.md` has why it is a C++
   change).
