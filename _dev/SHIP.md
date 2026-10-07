# Shipping bartisan 0.1.0

The standard applied throughout: **does this seriously affect a workflow the
package claims to support?** Speedups, polish and features that another package
already covers are excluded, however tempting, and go under "After the release"
below or in `TASKS.md`.

This file holds only what is current. How the list got here (the earlier status
reads, the surveys of what the other BART packages and *rstanarm* export, the
plans for the vignette series, and what the first full check found) is in its
git history before 2026-10-07, and the work behind each item is in `TASKS.md`.

## Where it stands (2026-10-07)

Submitted twice, and rejected twice. The first submission (`dd8f158`) drew a
NOTE for four examples that took more than 10 seconds on Windows and Debian. The
second (`dba4825`) was rejected because the whole check took more than 10
minutes on Windows. Both are fixed, the second in `131ad34`. What is left before
resubmitting is one Windows check, `devtools::check_win_devel()`, to measure
the time there rather than project it.

## The check-time budget

CRAN rejects a first submission whose check takes more than 10 minutes, and the
slowest machine decides. On the second submission, Debian took 23 s on the
examples, 92 on the tests and 470 re-building the vignettes; Windows took 49 s,
216 and 18 minutes. Windows was about 2.3 times Debian at every stage, and the
same code under CRAN's settings here (`NOT_CRAN` unset,
`_R_CHECK_LIMIT_CORES_=TRUE`) was 2.2 times faster than Debian and 5 to 6 times
faster than Windows. So a time measured here is multiplied by about six for
Windows.

At `131ad34`, the vignettes re-build in 13 s in `R CMD check` on a fresh
tarball, and the tests take about 28 s in CRAN mode. That projects to about one
minute and 2.8 minutes on Windows. Its install time is not in the log, so the
total is an estimate, roughly 6 to 8 minutes.

What keeps it there:

- **Vignettes.** Nine of the eleven replay saved results on CRAN
  (`vignettes/saved-results.R`); only implementation and faq run live. A chunk
  edited in any of the nine needs its results refreshed before a submission,
  with `just --justfile _dev/justfile vignettes <name>`, which installs the
  package, re-runs the vignette, and checks that the replay is byte-identical
  to the run. `_dev/check.sh` replays strictly, so a
  chunk edited since its results were saved fails the local check rather than
  warning on CRAN.
- **Tests.** The rule at the top of `tests/testthat/helper-bartisan.R`: a test
  runs on CRAN only if its assertions hold for any seed, on a fit as small as it
  allows. A `bcf()` fit in a test gets `propensity_args = small_propensity()`,
  since the propensity model otherwise runs at `bartisan()`'s defaults and is
  most of the cost.
- **Examples.** Each has to stay under 10 s on CRAN, so under about 1.5 s here.
  A `bcf()` call in an example passes a small `propensity_args`, as `?bcf`'s
  does. Slow parts go in `\donttest{}`, which `--as-cran` runs in a separate
  pass that is not timed.

## What the incoming check will still say

One NOTE on both platforms: a new submission, and "Linero" and "Tweedie" as
possibly misspelled words in `Description`. Both are expected. They go
unexplained, since `cran-comments.md` is not used, by the maintainer's decision
of 2026-10-07.

## Accepted for this release

- **Calibration of the default Poisson fit**, accepted as unsolved by the
  maintainer's decision of 2026-10-07, since every hypothesis tested for it was
  rejected. Simulation-based calibration of a Poisson response under the
  defaults (soft rules, bandwidth drawn per tree) gives ranks whose dispersion
  contrast is near +2 on 600 replicates in every arm that draws the bandwidth,
  across a tenfold range of bandwidth priors, ten times the warmup, four times
  the draws, three response levels and thinning to one draw in 200. The arm that
  fixes the bandwidth is uniform, and so are the soft-rule logit arms; the
  posterior is slightly too narrow. The move itself is verified (its acceptance
  ratio, its likelihood difference and its bookkeeping, the last to 1e-10 after
  72,488 accepted moves in `_dev/bandwidth-predictor.R`), so the cause is not an
  implementation bug in the move and is not yet located. See the `TASKS.md`
  entries of 2026-10-01 and 2026-10-02.
- **A deprecation warning under clang with libstdc++.** rhub's ubuntu-clang
  build warns that `std::get_temporary_buffer` is deprecated, through
  Armadillo's `inv_sympd()` in the multinomial probit code. Clang with libc++,
  which is what CRAN's clang machines use, does not, and rhub's clang22 build
  with libc++ was clean.
- **The `loglik` row of a `ph()` fit mixes slowly.** On the default `ph()` fit
  to `rhc` it reads R-hat 1.18 and bulk ESS 9 over 800 draws, while the survival
  probabilities and the `rhc` contrast from the same draws have bulk ESS 660 to
  810. Not investigated; both submissions went out with it, and `TASKS.md` lists
  it after the release.
- **`custom_family()` has no posterior predictive distribution**, since a log
  density supplies no way to draw from it. `?custom_family` says so, and says
  which methods error because of it. An `rng` argument is listed below.

## Optional before the release

- **A hex logo.** It wires into `_pkgdown.yml` and `README.Rmd`, so it is cheaper
  before the site and README are final than after.
- **`inst/CITATION`.**

## Decided against

Recorded so that these are not reopened as though they had never been asked.

- **`predictive_interval()` and `predictive_error()`** (2026-09-15). They are a
  quantile and a subtraction over `posterior_predict()`, which is exported and
  documented.
- **`fixef()` and `ngrps()`.** A forest has no fixed effects to return, and a
  group count is `length(levels(...))`.
- **A mixture-of-normals augmentation for `weibull_aft()`.** It could plausibly
  take the slowest family from about 8 s to about 1 s at 700 observations, but it
  is an approximation of the Gumbel error, and adopting it would weaken the
  exactness claim the package makes.

## After the release

In no order. These have entries in `TASKS.md`: correlated random effects across
additive predictors; a joint update for the ordinal cutpoints; an `rng` argument
for `custom_family()`; a `quasi()` family; interaction detection, assessed and
not decided; and keeping a row of a matrix response to a custom family whose
other outcomes are observed.

These are recorded only here:

- **`VarCorr()`.** `fit$tau` is the only thing it would return, and `?bartisan`
  names it.
- **The `loo_*` prediction wrappers** that *rstanarm* has (`loo_predict()`,
  `loo_linpred()`, `loo_predictive_interval()`, `loo_R2()`). The PSIS machinery
  is there; the wrappers are not.
- **`posterior_vs_prior()`.** Possible since `prior_only = TRUE`, and walked
  through as a recipe in `vignette("diagnostics")`. Whether it deserves a
  function is a judgment about discoverability rather than capability.
- **The DART inclusion probability as a stored quantity**, which is a C++ change;
  `TASKS.md` has why.
