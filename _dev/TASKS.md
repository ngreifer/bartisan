# bartisan tasks

Implements Linero (2025), "Generalized Bayesian Additive Regression Trees Models: Beyond Conditional Conjugacy", with SoftBart-style soft decision rules from Linero and Yang (2018).

**Status** and **To Do** describe the present. The rest is the record of how the package got here, by topic under seven parts. Each topic opens with one bold sentence giving the problem and where it stands, then one bullet per entry in chronological order. Entries are summaries: the full text of every entry up to 2026-11-05, with its tables, code and working, is in `_dev/TASKS-archive.md` under the same topic, date and title.

**Keeping it this way.** Add a new entry at the end of its topic (or under a new topic with its own bold sentence) as one bullet on one line, `- **DATE: Title.**` and two or three sentences: what was tested, the result with the numbers that decide it, and the decision or what is still open. Longer working belongs in the `_dev/` script's header, not here. A closed To Do item moves to the end of its topic, marked "(closed To Do item)". Whenever an entry changes what a topic's bold sentence says, rewrite that sentence in the same edit. Write each paragraph on one line, with no hard wraps. Undated entries carry the date of the commit that first recorded them, and "by 2026-08-24" marks those older than the repository.

## Status

`_dev/` was cleaned on 2026-09-30 (113 files, 18.8 MB); the rule for what stays is the `.gitignore` comment, and removed scripts are in the git history.

On 2026-09-17, 477 tests (2430 expectations) passed with `NOT_CRAN=true`, and `R CMD check` on the built tarball reported `Status: OK`. A green check alone is not evidence, since CRAN mode skips 91 blocks, including every ordinal and `prior_only` test. `--as-cran` has never completed on this machine, for a reason that looks environmental.

**What exists.** A C++ engine and a `glm()`-style R interface (`bartisan()`, `bartisan_control()`, `predict()`, `summary()`, parallel chains, `diagnose()`, `custom_family()`). Families: Gaussian, binomial (any link), Poisson, negative binomial, gamma, Tweedie, ordinal, multinomial (symmetric or reference), multinomial probit, the AFT families and `ph()`, location-scale families, zero-inflated Poisson and negative binomial, hurdle (two-part) Poisson, negative binomial, gamma and log-normal, beta and ordered beta, and Dirichlet process mixture errors. MIA for missing predictors, data augmentations on by default, `marginaleffects` support, `(1 | group)` random intercepts, posterior predictive draws with the `loo`, `bayesplot`, `performance` and `posterior` interfaces, and `bcf()`. Ten vignettes. 0.1.0 is on CRAN (2026-11-01); the version is 0.1.0.9000, with `NEWS.md` items since.

Benchmark (`_dev/benchmark.Rmd`): Friedman function, n = 1000, p = 10, 50 trees, 1000 warmup and 1000 saved, best of 3, RMSE against the truth on 1000 held out. Not comparable with versions before 2026-09-02 (see "Benchmarks against other packages").

| Task | Package and call | Seconds | ESS | RMSE |
|---|---|---|---|---|
| gaussian | dbarts bart() | 0.281 | 24.0 | 0.225 |
|  | bartisan hard rules | 0.905 | 21.8 | 0.193 |
|  | stochtree bart(num_gfr = 0) | 1.189 | 19.5 | 0.245 |
|  | BART wbart() | 1.932 | 23.1 | 0.231 |
|  | bartMachine bartMachine() | 2.082 | — | 0.234 |
|  | bartisan soft, smoothstep gate (default) | 2.887 | 50.1 | 0.130 |
| probit | dbarts bart() | 0.360 | 51.4 | 0.113 |
|  | bartisan hard, augment = TRUE | 1.094 | 35.8 | 0.107 |
|  | stochtree bart(binary/probit) | 1.876 | 34.0 | 0.122 |
|  | bartisan soft, augment = TRUE | 3.104 | 64.2 | 0.086 |
| logit | bartisan soft, augment = TRUE | 3.973 | 178.8 | 0.082 |
|  | BART lbart() | 48.884 | 68.6 | 0.108 |
| ordinal | bartisan hard, cloglog, augment = TRUE | 6.421 | 35.9 | — |
|  | stochtree bart(ordinal/cloglog) | 8.102 | 27.0 | — |

The full table, with the remaining gates, links and families, is in the archive. Hard rules are within 1.8x of dbarts at the same mixing and slightly better accuracy (best of six: 0.426 s against 0.241 s Gaussian, 0.482 s against 0.272 s probit); the two-replicate table is noisy at 10% to 30%, so ratios near two need longer runs. bartisan beats stochtree on the ordinal cloglog task (3.03 s against 4.01 s) and every package but dbarts on time once hard rules are used. Soft rules cost three to five times the hard-rule time and cut held-out error by 35% to 40%, which makes the default the most accurate fit in the table. Ordinal RMSE is not comparable across links.

## To Do

`_dev/SHIP.md` is the current release assessment. Items struck rather than done: the grow-from-root warm start and soft random tree features (the burn-in transient is 34 to 70 sweeps against a warmup of 200), and the joint update for correlated nuisance parameters (not needed once `tweedie()` made it testable).

### Before submission

- [ ] **Refresh the saved-results vignettes other than `diagnostics`** (2026-11-01). They replay draws from before the single-chain stream change and predate the `add_chains()` advice in `diagnose()`; check their prose numbers against new output, as was done for `diagnostics`.
- [ ] **A hex logo.** None exists, so `_pkgdown.yml` has no `logo` and the README no badge row. `usethis::use_logo()` wires it up once a 1200x1390 image exists; adding it later means rebuilding the site and README.

### After submission

Ordered by expected value. None is load-bearing for a workflow the package claims to support.

- [ ] **The `loglik` row of a `ph()` fit mixes slowly** (moved here 2026-10-07). On the default `ph()` fit to `rhc`: R-hat 1.18, bulk ESS 9 over 800 draws, while the survival probabilities and contrast reach ESS 660 to 810 (`_dev/ph-level-shift.R`). Not looked at; check first whether the log likelihood tracks forest size or the bandwidth.
- [ ] **Test coverage is 91.7%** (2026-10-06, `_dev/coverage.R`, 1h 20m): R 92.9%, C++ 90.8%, all 552 blocks pass under covr. `methods.R` and `interop.R` were covered the same day, which found and fixed two bugs (`prior_summary()` miscounting `ph()` bins, and `kfold()` doubling a formula offset and failing on a matrix one). The covr object is `_dev/coverage-2026-10-06.rds`.
- [ ] **`src/family.cpp` is at 89.5%, and tests cannot take it much further** (2026-10-06). Of the 325 uncovered lines left, 208 are overrides nothing calls, 41 closing braces, 47 numerical guards and 29 guarded on the R side. Deleting the unused overrides, with `Family::dlogdens()` and `Family::info()`, is a code decision, not a test.
- [ ] **Interaction detection, and whether DP-Forests belong here** (assessed 2026-09-30, not decided). Du and Linero (2018, arXiv 1809.08524) show the sparsity prior reports spurious interactions and fix it with per-tree splitting proportions from a DP mixture. First record per draw which predictor pairs share a tree and measure that summary's false positives on their designs; only then consider DP-Forests as a `sparsity` option.
- [ ] **A formal variable-selection test** rather than a threshold on `prop_used`, as dbarts, SoftBart and bartMachine have. `variable_importance(draws = TRUE)` gives the per-draw counts; what is missing is a calibrated null, such as the permutation approach.
- [ ] **A comparison between a family with a point mass and one without.** `loo()` cannot compare `tweedie()` with `gaussian()` or `dpm()` on exact zeros, and `vignette("comparison")` dropped that section on 2026-09-24. Candidate: leave-one-out CRPS or SCRPS through *loo* 2.10's `loo_crps()`, which needs methods for bartisan fits and a comparison function.
- [ ] **Correlated random effects across additive predictors.** Blocked by the lack of mixed second derivatives in the `Family` interface; the alternative needs a prior mean threaded through 27 places.
- [ ] A Bayesian-bootstrap dispersion draw for `Gamma()` and `negbin()` from Pearson residuals, the cheapest improvement to interval calibration.
- [ ] A joint tridiagonal update for the ordinal cutpoints, for inference on many thresholds. The obstacle is the ordering constraint.
- [ ] The `draw_prior` move from SoftBart, proposing a whole fresh tree to escape local modes. Needs an `L`-dimensional Laplace proposal for the new leaves.
- [ ] An optional `rng` argument to `custom_family()`, giving it a posterior predictive distribution and so `simulate()`, `pp_check()` and `r2()`.
- [ ] A prior weight on the missingness rule, so a variable with a few missing values does not spend a third of its rules on them.
- [ ] A `quasi()` family parameterized by link, variance function and dispersion update. Needs a documented weakening of the exactness claim.
- [ ] A lighter-tailed prior on the leaf scale, or an upper bound, to remove the separation pathology. Changes the default prior; the warning is the interim measure.
- [ ] **`fc()` for a fixed coefficient**, a parametric block such as `fc(poly(z, 2))` beside the forest. Against: it needs a prior of its own and answers for `summary()` and `coef()`, while `vc(z, ~ 1)` already fits a constant coefficient. Out of scope for 1.0 unless users ask.
- [ ] **The Gibbs prior of Linero and Du (2023), as `sparsity = "gibbs"`.** Built, measured and removed before 0.1.0 because its benefit is in selection and the selection test above is unbuilt. Restore with `git apply _dev/gibbs-prior.patch`, then `rm src/random.o` before rebuilding, since its object goes stale unnoticed.
- [ ] **Relative survival on top of `ph()`** (Basak et al. 2024): one extra Bernoulli draw per sweep with a life-table population hazard. Cheap; worth doing on request.

Candidates raised on 2026-10-02 while rewriting `vignette("survival")`, none prototyped or checked against a reference:

- [ ] **A smooth M-spline baseline hazard for `ph()`**, keeping the exponential target and exact gamma conditionals. Expect smoothness rather than accuracy; compare `S(t | x)` and `loo()` against the bins. Not a Cox-style baseline, which broke `loo()`.
- [ ] **Counting-process input for `ph()`, `Surv(start, stop, event)`**, for delayed entry, time-varying covariates and recurrent events. The change is in `PHFamily`'s cumulative baseline and `prepare_surv()`; leave-one-out would be grouped by subject.
- [ ] **Competing risks**: combine cause-specific `ph()` fits into cumulative incidence, with a prediction type and an `avg_comparisons()` route.
- [ ] **Cure models**: a cured fraction with its own logit-scale forest, built like the zero-inflated families.
- [ ] **Interval censoring**, `Surv(L, R, type = "interval2")`. Natural for the AFT families through their imputation; harder for `ph()`.

## Calibration and the sampler

### Soft-rule posteriors narrower than calibration allows

**One chain's posterior of a soft-rule contrast reads too narrow in simulation-based calibration for Poisson and the augmented binomial; the structure moves, the bandwidth move, the anchor, the change move, mixing and the family have each been ruled out, bartisan's Poisson leaf update matches Stan's on the same model (its binomial intervals are 1.4% narrower than Stan's), and on the fresh seed block that decides whether block 5000 was chance, the harness is at fault, or the full sampler narrows, the Poisson arm found Stan's ranks uniform and no narrowing in the full sampler, with the binomial arm not yet recorded.**

- **2026-10-01: The soft-rule Poisson SBC deviation is the bandwidth move.** The soft-rule Poisson arm of `_dev/sbc.R` failed SBC (chi-square 28.4 on 9 df, p = 0.001); mixing, a misreplicated prior (`_dev/sbc-prior-check.R`) and the leaf draw were falsified, and a second block (seeds 90000) replicated it on the mean-rank contrast (pooled chi-square 29.3, p = 0.0006). Holding the bandwidth fixed removed it (dispersion -0.24) while soft logit and hard Poisson showed nothing, and `_dev/bandwidth-mixing.R` gave the Poisson bandwidth a median ESS of 68 of 1000 against 190 under logit. Lessons: coverage from 1000 autocorrelated draws carries a finite-chain shortfall (about 0.93) the ranks do not, and editing a script a queued job reads kills the job.

- **2026-10-02: The bandwidth move's likelihood term is correct, and so is the compiled Poisson family.** `.bartisan_loglik_delta()` in `src/model.cpp` showed `Family::loglik_delta()` agreeing with two full `total_loglik()` calls to 1e-8 for five families (test in `test-invariants.R`). A `custom_family()` Poisson matched the compiled one's contrast posterior SD to one percent (`_dev/poisson-paths.R`). What both paths share is the Laplace approximation, and the proposed next test was correlating each replicate's rank with its bandwidth.

- **2026-10-02: Conditioning SBC ranks on the parameter is not a test, and two things it cost.** An `SBC_REPLAY_ONLY` mode put the deviation in the narrowest gates (Spearman -0.094), but conditioning ranks on a component of the generating parameter breaks uniformity even for a correct sampler, so the tercile table measured nothing. The valid version, varying the prior (`SBC_BANDWIDTH`), gave a positive dispersion contrast at prior means 0.03, 0.1 and 0.3 (about +2.2), vanishing only with the bandwidth fixed. Costs recorded: a patch asserting unread text wrote nothing so the unpatched script ran, a filename without the seed block overwrote a result (recovered via `SEED_BASE + r` seeding), and `timeout` does not exist on macOS.

- **2026-10-02: Not the step floor or the thinning, and the harness miscounted its chi-square.** The finding that the deviation grows with the counts (+4.3 at offset 4) was withdrawn on 2026-11-04 as an anchor artifact. `_dev/bandwidth-step-floor.R` found the step floor not binding (acceptance near 0.44) though the bandwidth's ESS fell to 15 of 2000, and thinning to one in 200 left dispersion at +2.09 against a predicted -0.2. The chi-square assumed sixty ranks a bin (fixed to expected counts per value; contrasts unaffected), and a filename without `draws` let task 276 overwrite task 274 (fixed); the default `gate = "smoothstep"`, `update_bandwidth = TRUE` is the deviating configuration.

- **2026-10-02: The bandwidth move's bookkeeping is exact.** `_dev/bandwidth-predictor.R` compared the recorded predictor with `predict()`'s replay of the stored trees over eight soft-rule Poisson fits: the largest gap was 1.0e-10 after 72,488 accepted bandwidth changes, with nothing accumulating. With the acceptance ratio, likelihood difference and bookkeeping all verified, the implementation of the bandwidth move is cleared.

- **2026-10-02: The tree moves are exactly reversible, as written.** Reading `node_birth()`, `node_death()` and `change_rule()` in `src/mcmc.cpp` found the target ratio, tree-shape prior, rule prior, move-choice probabilities and proposal densities all correct; the Laplace fits are deterministic because Fisher scoring starts at `SCORE_INIT`. The `_dev/sbc.R` generator matches the package except `x_transform = "range"`, which rescales predictors by about 0.5%.

- **2026-10-02: A reference sampler for the soft-rule Poisson SBC, stopped.** `_dev/sbc-stan.R` fitted the model with Stan's NUTS (`_dev/sbc-stan.stan`), tree structures given as data, and was stopped after 28 replicates for cost (corrected 2026-11-05: 14.4 minutes a replicate, about 72 core-hours for 300). The 28 ranks are not evidence; Stan's draws were sound (no divergences). A cheaper version needs a smaller model in both samplers.

- **2026-11-04: Is the soft-gate Poisson narrowing the general leaf path or the Poisson family? A binomial split.** A 20-trial binomial (offset -2.7, `anchor = 0`) fitted augmented and direct both narrowed, dispersion z +2.23 and +1.93 (paired z -0.71, tasks 49 and 50), clearing the general leaf path and anything specific to `poisson()`. The prediction that the bandwidth move under a skewed likelihood would then be the suspect was wrong: in `LogitAugmentedFamily` the move sees an exactly Gaussian likelihood, so what separates the narrowing arms from passing soft logit is how much each observation says; a move that does not hold the leaves fixed was next.

- **2026-11-04: Does the soft-gate narrowing reach `negbin()`?** With `negbin(theta = 5)`, the augmented and direct routes read dispersion z +1.75 and +1.28 (paired z -1.00, tasks 51 and 52), so the route makes no difference and soft-rule `negbin()` most likely shares the narrowing at a smaller size, with coverage (.937, .943) near the finite-chain level. Effect sizes order by how informative each observation is (Poisson about +2.6, Bernoulli logit +0.3), which is consistent with the conditioning explanation and not evidence for it.

- **2026-11-04: Does holding the leaves fixed in the bandwidth move cause the narrowing? A collapsed move.** A collapsed bandwidth move behind `the$collapsed_bandwidth` integrates a tree's leaves out exactly and redraws them jointly; it passed its checks (replay to 7.3e-11). On the augmented binomial it read dispersion z +2.78 against +2.23 (paired z +1.24, the wrong direction; task 53), so the claim is false. The move stays in `src/mcmc.cpp`, off by default, uncommitted and without tests; it doubled the contrast's ESS on one fit (65 to 137), unmeasured as an efficiency change.

- **2026-11-04: A single tree, and several chains: two tests of where the narrowing lives.** Prompted by the BART mixing literature, four pooled chains left 20 soft trees narrowed (Poisson +3.20 against +3.79, binomial +2.53 against +2.57), and chains started from the prior (`the$init_from_prior`) changed nothing (+2.95, +2.68), so the chains converge to a too-narrow distribution. A single tree exposed a tie defect (`sum(thin < truth)` broke the atom at zero downward; fixed with a randomized rank in `_dev/sbc.R`), then narrowed (soft +7.76, hard +4.69) and was fixed by four chains (-1.31, -4.37), a mixing problem separate from the 20-tree one. Dropping the change move (`the$no_change_move`) against a same-build baseline at 600 replicates left both families narrowed (Poisson +2.27 against +2.61, binomial +3.41 against +2.23); the earlier paired z of -2.60 belonged to replicates 1 to 200, not the build. The harness now fits replicates in parallel with `prog_parallel()`.

- **2026-11-04: Does the empirical anchor narrow the soft-rule Poisson posterior?** `SBC_NO_CENTER` cancelled the anchor; at offset 0 dispersion read +2.61 against +2.80 (paired z +0.43), so the anchor is not the cause. Reading the code found the Poisson anchor ignores a user offset (`prepare_response()`, `build_offset()`); at offset 4 cancelling it moved +4.64 to +3.08 (paired z +1.67), so "it grows with the counts" was withdrawn. Not predicted: at offset 4 the level is too narrow too (+3.68), on one seed block only.

- **2026-11-05: The tree structures fixed at the generator's: is the narrowing in the structure moves?** `SBC_UNIT_X` and `SBC_FIXED_STRUCTURE` held every tree at the generator's structure (task 143): Poisson read +2.05 fixed against +1.69 free and the binomial +1.46 against +2.38 (paired z +0.53 and -1.30), so the structure moves are not the cause, nor is the 0.5% rescaling. The mean-rank z was positive in all four arms (about +2.6 together), not predicted; fixing the bandwidths too was the next split.

- **2026-11-05: Structure and bandwidths both fixed at the generator's: is the narrowing in the bandwidth move?** With `SBC_KNOWN_BANDWIDTH` only the leaves are drawn (task 145): Poisson +1.64 (paired z -0.65) and the binomial +1.31 (paired z -0.23), so the bandwidths are not the cause. These arms share block 5000's data sets with every arm above; what was left was the harness or slow mixing among coupled leaves, and a Stan GLM reference on the leaves became cheap.

- **2026-11-05: Stan on the leaves-only model: is the narrowing bartisan's, or the harness's?** `SBC_STAN` with `_dev/sbc-stan-glm.stan` (task 148) gave Poisson dispersion +1.56 against bartisan's +1.64, ranks correlated .98 and width ratio 0.999, so bartisan's Poisson leaf update matches Stan and the +1.64 is block 5000's chance, shared by every arm run that day. For the binomial Stan read +0.74 against +1.31 (paired z +1.50), but bartisan's intervals were narrower by a ratio of 0.986 (CI 0.981 to 0.991), read as Polya-Gamma coupling slowing the chain; the prediction that the gap closes at 4000 draws has not been run.

- **2026-11-05: A fresh seed block: chance on block 5000, or the harness, and does the full sampler narrow?** On `SBC_SEED_BASE=300000` (task 149) the Poisson arms read full sampler +0.87 (coverage .953), bartisan leaves-only +1.02 and Stan +0.47, so the harness is sound, block 5000 was chance, and no narrowing is detectable at 600 replicates on the current build; claim (2) was wrong for Poisson. The augmented binomial read full sampler +0.66 (coverage .955), bartisan leaves-only +0.07 and Stan +0.03 (ranks correlated .97, paired z +0.12), so claim (2) was wrong for it too; every earlier binomial reading came from block 5000, so the soft binomial narrowing has no support outside that block. The binomial width gap recurred (bartisan against Stan 0.991, CI 0.986 to 0.996) with uniform ranks, which fits the finite-chain reading. The comparisons between arms on block 5000 stand as comparisons. The October results (+1.62, +1.88, +2.65 on an older build) are the only evidence left for any narrowing and need replication on the current build before anything in the sampler changes.

- **2026-11-05: October's three soft-Poisson readings, rerun on the current build.** *Tested:* whether October's +1.63 (prior mean 0.1, block 90000, `-rep2`), +1.88 (prior mean 0.03) and +2.65 (prior mean 0.3) belong to their data sets or to the build, and whether the full sampler narrows on the current build at any prior. The 0.03 and 0.3 runs used block 5000's seeds, and a different prior mean only rescales the same exponential draws, so they drew block 5000's trees and leaves with wider or narrower gates and share part of its chance. *Measured:* six arms of 600 replicates in October's configuration (raw design, the empirical anchor, the full sampler, one chain): those three settings on their own blocks, and the same three on fresh blocks 400000, 500000 and 600000. Read each dispersion z, and the three fresh arms pooled (sum over the square root of 3). *Meaning:* the same-block arms reproducing October within about 0.6 and the fresh arms pooled under about 1.5 means October's readings were their data sets' chance, the current build has no detectable narrowing at any prior, and the investigation closes; the fresh arms pooled at about 2.5 or more means a small narrowing that depends on drawing the bandwidth remains, to be located with paired arms on fresh blocks; the same-block arms falling to near zero means the build changed something since October, and the fresh arms then describe the current state. *Result* (task 150, 46 minutes): the first branch. The same-block arms reproduce October on the same truths, ranks correlating .97 to .98: +1.50 against +1.63 (block 90000), +1.31 against +1.88 (0.03) and +2.17 against +2.65 (0.3), paired differences -0.32, -1.25 and -1.16, none significant. The two block 5000 arms correlate .68 and .73 in rank with that block's prior-0.1 arm, so they carried its chance. The fresh arms read -1.02 (0.1, block 400000, coverage .945), +1.07 (0.03, block 500000, .943) and +0.16 (0.3, block 600000, .957), pooled +0.12, and their mean-rank z is -0.48 pooled. **No narrowing is detectable on the current build at any of the three priors**, and nothing found over 2026-11-04 and 2026-11-05 calls for a change to the sampler. The internal flags added for the investigation (`collapsed_bandwidth`, `collapsed_check`, `init_from_prior`, `no_change_move`, `fixed_trees`, `fixed_structure`) stay off by default; whether to keep them is open.

### Credible interval coverage

**Whether pointwise credible intervals reach nominal coverage: remeasured on 2026-09-10 at 40 replicates, they are slightly conservative across four families (.96 to .97 at a nominal .95) rather than deficient, with one data-generating process as the caveat the vignette now states.**

- **by 2026-08-24: Why credible intervals miss nominal coverage.** On the Friedman function (n = 250, 12 replicates) default coverage was 0.945 with the posterior mean about 0.77 posterior SDs from the truth. The mechanism is bias, not mixing: four chains or a fourfold longer chain moved coverage by 0.005, leaf-prior changes and 200 trees did not help, and hard rules over-covered (0.966) through a wider posterior. The one lever that worked was the mean-variance relation (see quasi-likelihood).

- **2026-09-05: Log: simulation-based calibration, and what the interval question turned out to be.** SBC on a `binomial()` contrast (300 replicates, `_dev/sbc.R`) passed: chi-square 4.2 (p = 0.90), coverage 0.957 (SE 0.012), so the sampler draws from its posterior. The short coverage on the benches (0.79 to 0.88) is credible intervals under prior shrinkage of a fixed truth, not a defect. A contrast is calibrated because the `qlogis(mean(y))` centering cancels; a soft gate, a drawn sparsity prior and a nuisance-parameter family were not yet covered.

- **2026-09-10: Interval coverage, remeasured: conservative, not deficient.** `_dev/coverage-calibration.R` (40 replicates, n = 500) found coverage .961 to .973 across `gaussian()`, `binomial()`, `poisson()` and `Gamma("log")`, with bias-to-SD ratios .69 to .74; the vignette's .91 binomial figure and ratio near .8 did not survive. Four chains and five times the draws moved nothing. The caveat is one data-generating process, and the vignette now says so.

### Soft rules: the gate and the bandwidth move

**Soft rules buy accuracy through a per-tree gate bandwidth, which is the worst-mixing quantity in a soft-rule fit; the move stays at every sweep (`bandwidth_every = 10` was faster and mixed worse), and `diagnose()` summarizes the bandwidth row over trees because its convergence binds on nothing a fit reports.**

- **by 2026-08-24: Soft rules: the gate, and the bandwidth move.** The bandwidth move was 48% of a soft-rule fit. Bounded gates `"smoothstep"` and `"smootherstep"` (Beta(2, 2) and Beta(3, 3) CDFs) run at 0.68 to 0.71 times the logistic's time with equal accuracy, because they save the `exp()` rather than far-side work; raised cosine was rejected and `"linear"` dropped. Fixing the bandwidth is faster but more than doubles error on a step function (0.4624 against 0.1901), so drawing stays the default with `bandwidth_every` exposed. Two exact improvements: a split-less tree draws its bandwidth from the prior by Gibbs, and rejections roll back from a `SupportStore` snapshot (17% of a soft fit to 0.7%).

- **2026-09-10: `bandwidth_every = 10` is faster and mixes worse, and the default stays at 1.** `_dev/bandwidth-confirm.R` (12 cells, 240 fits) found `every = 10` 1.12x to 1.64x faster with 1.26x ESS per second, but RMSE 1.2% worse on average, with the significant losses on step functions (+4.8%, +6.6%), the case the move exists for. The default stays 1, and `R/control.R` records when raising it is reasonable.

- **2026-10-02: The gate bandwidth is the worst-mixing quantity in a soft-rule fit.** In `_dev/bandwidth-rows.R` the bandwidth's worst-5% row carried the fit's worst R-hat in 17 of 24 fits; at n = 8000 the worst tree's bandwidth reached R-hat 1.38 to 1.52 on 8 to 18 effective draws while the averaged predictor had 1865 to 4021, and worst-tree ESS fell 5 to 20 times from n = 500. It is graded apart because the chains agree on the function; the row is in `R/diagnose.R` as of `dfa1867`, summarized over trees and omitted for hard fits.

### The level a fit shrinks toward: `anchor`

**The level a fit shrinks toward ignored the offset, which loosened the leaf scale and widened Poisson rate intervals by 17% at large exposure gaps; fixed on 2026-11-04 by moving the start with the offset (`anchor_at_offset()`), and `bartisan_control(anchor = )` now sets that level directly.**

- **2026-11-04: Does the offset-blind anchor cost a Poisson rate model anything?** `_dev/offset-anchor.R` (tasks 40 to 42, 100 replicates per arm) found that at a gap of 5.5 the drawn `sigma_mu` was 2.7 times the corrected value, intervals 17% wider and coverage 0.961 against 0.936, with RMSE 1.8% worse (about two SE); interval width depended on exposure units. Fixed the same day: `prepare_response()` records `anchor_by` and `anchor_at_offset()` applies it in `bartisan()` after `expand_for_vc()`, making the Poisson start `log(sum(w * y) / sum(w * exp(offset)))`, with a NEWS item and *patrick* tests in `test-invariants.R` that a constant offset leaves the chain unchanged. The binomial `vc + offset` test became an explicit `expect_lt()` because testthat 3 reads `tolerance` as relative.

- **2026-11-04: `anchor` in `bartisan_control()`: setting the level the fit shrinks toward.** After checking by fitting how *dbarts*, *BART*, *brms* and *rstanarm* separate offset from level, `anchor` sets one link-scale value per additive predictor (`NA` keeps the intercept-only fit), with the offset added to it; `custom_family()` and `ph()` refuse it, and it lives in `apply_anchor()` in `R/response.R`. The draft documentation wrongly said `update_sigma_mu = FALSE` makes the anchor firmer; a probit measurement showed the opposite, and `?bartisan_control` was corrected. `_dev/sbc.R`'s `SBC_NO_CENTER` uses `anchor = 0`; tests are in `test-invariants.R`.

### The leaf scale (`sigma_mu`)

**`sigma_mu` mixes poorly on every design tried and runs away under separation; the slow mixing is shared by dbarts, stochtree and BART and leaves the fitted function alone, so the question is closed, `sigma_mu` left the R-hat table, and separation gets a warning rather than a new prior (a lighter-tailed prior is on To Do).**

- **by 2026-08-24: Separation runs the leaf scale away.** On fully separated binary data the half-Cauchy-drawn `sigma_mu` wandered between 2.79 and 9.29 and the predictor reached 110; pinning it with `update_sigma_mu = FALSE` held the predictor at 6.9. `bartisan()` now warns when the leaf scale's posterior mean settles more than five times above its prior median.

- **by 2026-08-24: Firth-type penalization: wrong direction for the usual bias.** BART's bias is shrinkage toward the null fit (slope -0.235 of error on truth), so Firth's shrinkage penalty would worsen it, and its information-dependent part is small (R-squared 0.213 to 0.250) and of the same sign. Separation shows in the leaf scale, not the leaves; `log det I` also does not decompose leafwise and would change the target.

- **2026-08-31: `sigma_mu` mixes poorly across every design tried.** `_dev/sigma-mu-mixing.R` (`gaussian()`, n = 1000, four chains) found `sigma_mu` bulk ESS under 15 of 3000 draws with R-hat 1.22 to 1.85 on every design, including the package's own benchmark.

- **2026-08-31: The ramp left the leaf scale frozen for most of warmup.** Nothing re-enabled the `sigma_mu` update after `sigma_mu_ramp` ended, so with the default 0.25 it was frozen for 75% of warmup and jumped at the first retained draw. Fixed with an `else if (iter == num_ramp)` branch and two regression tests; it was a real bug but not the cause of the poor mixing.

- **2026-08-31: Why `sigma_mu` mixes badly: three hypotheses, two refuted, and an honest stop.** `_dev/sigma-mu-cause.R` ruled out the sampler (acceptance 0.89 to 0.97), the forest's mixing (`sd(eta)` ESS in the thousands) and slow ergodic mixing (ESS 10.1 at 750 draws, 8.1 at 4000). Chains settle at slightly different levels against a very tight within-chain spread, which inflates R-hat; the leaf-versus-structure trade-off hypothesis was unsupported, and why chains settle apart was not established.

- **2026-08-31: The other packages do it too, so the leaf scale is closed.** In `_dev/sigma-mu-others.R` dbarts (R-hat 1.122, ESS 21.8) and stochtree's exact conjugate Gibbs draw (1.163, 17.1) matched bartisan (1.188, 14.8), so it is a property of the ensemble, not a bug. What remained was presentational: drop or mark `sigma_mu` in `fit$rhat` and say so in `vignette("diagnostics")`, not done in this entry.

### Observation weights

**The residual-scale prior ignored the weights (fixed), and frequency weights on a binary response run the leaf scale away exactly as replicating the rows does, which is the data rather than a weighting bug, so a weighted binary arm has to fix the leaf scale.**

- **2026-08-25: The residual-scale prior was ignoring the weights.** `residual_scale()` now fits its anchoring line weighted, with weights normalized to average one over the rows they keep; normalizing over all rows put the scale `sqrt(2)` out on a half-zeroed sample. Zeroing a noisy half now reproduces the clean half's scale and constant weights reproduce the unweighted answer exactly.

- **2026-10-01: Frequency weights on a binary response run the leaf scale away.** The recovery matrix's "augmented binomial mixes badly" was not the augmentation, the sample size or the omitted grouping factor (`_dev/augmented-mixing.R`, `_dev/contrast-mixing.R`, `_dev/contrast-mixing-cause.R`). In `_dev/weights-runaway.R` weights of 1 to 3 under hard rules ran `sigma_mu` from 0.649 to 1.335 and left the contrast 3 effective draws and an effect of 1.15 (truth 0.8); fixing the scale restored 497 draws and 0.77, and longer warmup did not. A weight on a 0/1 response is separable trials, which `warn_runaway_scale()` caught at 6.3 times but not soft rules at 3.6; weighted and replicated fits agree (`_dev/weights-mixing.R`), so a weighted binary arm must fix the leaf scale.

### How many trees

**Whether the number of trees wants a prior or a data-dependent default: Infinite BART was prototyped and does not learn the count, McCartan and Huang's ablation stops holding at fifty trees, and the default is unchanged.**

- **by 2026-08-24: The number of trees, and why nobody puts a prior on it.** On a binary problem the held-out log score varied by 1.05 of 227 from 5 to 200 trees, while `sqrt(m) sigma_mu` stayed at 1.89 to 2.09, so the tree count is very nearly unidentified. A prior on it would need reversible jump over ensembles at up to 40 times the cost; the caveat is a smooth five-dimensional test function.

- **2026-08-30: How many trees, and does the default depend on anything.** On Friedman, soft rules peaked at 20 trees (RMSE 0.270) while hard rules kept improving to 200, but on `MatchIt::lalonde` 20 trees cost 35% between-chain disagreement against 9% at 50, so 50 stays the default for both gates. `location_scale()` at `c(50, 5)` matched `c(50, 50)` at a third of the time, documented in `?bartisan_control` rather than defaulted. The `marginaleffects` help page's claim that more trees remove the atom at zero was corrected (20% at 50 trees, 18% at 200).

- **2026-08-30: Infinite BART (Battiston and Luo 2025): assessed, not implemented.** It is implementable as a mask on each tree's root membership weight, but the clustering cannot reach a new observation because the IBP prior ignores the predictors. Reproducing the paper's clustered Friedman design, dbarts and bartisan reached 31.3 to 32.5 against the paper's classic BART at 38.14 and infinite BART at 29.80, closing about 80% of the reported gap; a covariate-dependent prior for `W` would change the decision.

- **2026-08-30: Infinite BART, implemented as a prototype and measured: it does not learn the tree count.** A prototype in a scratch copy reproduced plain bartisan in its dense limit, but kept almost every tree (46.5 of 50, 175.3 of 200) at twice the RMSE and fifty times the time of a 10-tree fit. Its active-tree count depended on the start (161 against 79, 200 against 114 under one prior), and alone it gave negative variable-importance separation (-0.0013) where DART separated cleanly. Not shipped.

- **2026-08-30: McCartan and Huang (2026): their ablation replicated here, and where it stops holding.** Both ablation findings replicated, and soft random tree features beat hard ones most at few trees (+0.157 R-squared at 5, +0.033 at 500). Their conclusion holds only when covariates mostly matter: with 25 irrelevant predictors DART stayed at 0.950 to 0.958 while random features reached 0.874 at 500 trees. Nothing changed; four follow-ups were listed (soft random features as a fast fit or warm start, for non-Gaussian families, a sparsity-aware feature draw, tuned versus drawn `sigma_mu`).

### Data augmentation

**Data augmentations (Polya-Gamma and latent-variable rewritings) speed up the binomial, ordinal, multinomial, zero-inflated and negative binomial families; remeasured at fifteen replicates they win on effective draws per second, so since 2026-09-11 `augment = TRUE` applies to every family with a rewriting under either gate, and `control$augment` in a fit reports what was used.**

- **by 2026-08-24: Data augmentations: what each is worth.** Effective draws per second ranged from 22 to 24x for hard-rule `ordinal("probit")` to 0.8x for soft-rule `negbin()` (then off) and 1.6x for `multinomial()` (by name). The negative binomial moved from Polya-Gamma (0.5x) to a Poisson-gamma mixture (1.2 to 2.0x hard), ordinal cutpoints are drawn with latents integrated out (`update_ordinal_cuts()`) to avoid Albert and Chib's O(1/n) trap, and the ordinal logit rests on an exact PG(2, |r|) identity. Poisson cannot be helped; the recorded trap is to compare identified quantities (multinomial `eta` correlated 0.898, probabilities 0.995).

- **by 2026-08-24: Windle, Polson and Scott's saddlepoint Polya-Gamma sampler.** Declined: the exact Poisson-gamma route served the negative binomial better, and the saddlepoint method is an approximation whose envelope could not be validated.

- **2026-08-30: Ultimate Polya-Gamma samplers (Zens, Fruhwirth-Schnatter and Wagner 2024).** The representation is already the ordinal logit's identity, and the boosting does not port because it needs the coefficients marginalized as one block. The imbalance pathology it fixes is absent: displacing the anchor changed nothing, the level arrives in about 25 draws, and dbarts reproduced the level's ESS at 27 events (24.4 against 24.8). `?bartisan_control` now says few events need more draws.

- **2026-08-30: augment is a sampling setting, and what it would take to extend it.** `augment` moved to the advanced group of `?bartisan_control`. A general target cost 5 to 10 times a Gaussian fit (`lognormal_aft()` 10.0x); the plan was censored-normal imputation for `lognormal_aft()`, imputation plus Polya-Gamma for `loglogistic_aft()`, the exponential form for `weibull_aft()`, and nothing for `ordbeta()`, with `exp_sign()` generalized to `exp_rate()`. Trap recorded: `exact_quadratic = FALSE` is not a proxy for a general target.

- **2026-09-03: Log: `control$augment` in a fit is now an answer, not a request.** `fit$control$augment` is now a logical reported from C++ by `.bartisan_fit()`, since whether a rewriting applies depends on the data (`binomial("probit")` with 3 trials is `FALSE`, logit with 3 trials `TRUE`). The request is kept in `attr(control, "supplied")`.

- **2026-09-10: What data augmentation buys, at fifteen replicates instead of three.** A three-replicate rerun was discarded because the worst-quantity ESS ratio swung a median 4.5-fold within a cell. At 15 replicates (480 fits) augmentation was 1.8x to 31x faster (median 9.2x), its median worst-quantity ESS ratio 1.046, and ESS per second favored it in 16 of 16 cells, so the vignette's "every one trades speed for mixing" was corrected. Lesson: ranges over problems, as `R/control.R` reports, are more durable than a single-problem point.

- **2026-09-11: The negative binomial's soft-rule exclusion did not replicate, and a run that measured nothing.** At 15 replicates soft-rule `negbin()` augmentation gave 1.31x ESS per second, so `augment = TRUE` now covers every family whatever the rules and `resolve_augment()` no longer takes `soft`. The first attempt measured the default against itself because its "on" arm passed `augment = TRUE`; lessons recorded: a toggle arm must name what it wants (`on = "negbin"`), and a ratio of exactly 1.00 deserves suspicion.

### Sampler speed

**How the sampler was made fast: the leaf target's shape (quadratic, or exponential at a rate) gets an exact or near-exact shortcut, a bit-identity harness guards every change, and survival and prior-only fits were brought onto the fast paths; the zero-inflated and multinomial exponential route is still open.**

- **by 2026-08-24: Log: how the sampler got fast.** A probit fit to `lalonde` went from 82.3 s to 0.72 s with hard rules (dbarts 0.20 s), mostly from exploiting the shape of the leaf-level target rather than tighter code.

- **by 2026-08-24: Where the time goes now.** Instrumented at n = 1000, a hard fit took 0.454 s and a soft one 2.157 s, of which the bandwidth move was 46.5% (`rebuild_support` 34.1%). Both remaining costs are the algorithm rather than the implementation.

- **by 2026-08-24: Three traps that corrupted measurements before any of this.** `devtools`, `pkgload` and `roxygen2::roxygenize()` build at `-O0`, a 5 to 20 times slowdown; `.bartisan_optimized()` now detects it, `bartisan()` warns, `_dev/benchmark.Rmd` refuses to run, and `.Rprofile` sets `options(pkg.build_extra_flags = FALSE)`. Stale object files require timing from a clean build, and nothing beyond `-O2` helps.

- **by 2026-08-24: Why it was slow: the decomposition that made the rest possible.** Counting observation visits showed only 4 to 9 times a conjugate sampler's work, while probit cost 66 ns per visit against the Gaussian's 6. The lesson that drove the rest: the expensive part is a leaf target that is not quadratic.

- **by 2026-08-24: The conjugate shortcut: where the target is quadratic, the Laplace fit is exact.** `Target1` and `Target2` do one fused pass and a closed-form draw when the family is quadratic and forward to the general path otherwise. The general path had stopped short for quadratic targets (fixed; the paths now agree to 1e-15), and `exact_quadratic = FALSE` forces it; probit hard went from 1.50 to 0.94 s.

- **by 2026-08-24: The exponential form: a second shape, for the count models.** Poisson and log-link gamma targets of the form `c + a*eta + b*exp(s*eta)` are recovered from one pass, hard rules only, worth 1.86x and 1.89x. The gamma family now reports observed rather than expected information; the coefficient recovery was exact to 5.7e-14 and off by 1.16 for binomial, so the check can fail.

- **by 2026-08-24: Tree bookkeeping: three items, and the one that mattered was not the documented one.** The documented `make_base` fusion was measured at about 2% and not done; folding `make_child_weights` into `split_support` and a node pool were bit-identical and worth 5%. Filling `split_support` by index instead of `push_back`, earlier rejected on wrong reasoning, was 20.6% faster on soft fits; weight removal and fused sums gave about 10% (buffers stay in L1), and devirtualization 1.05x to 1.20x.

- **by 2026-08-24: The bit-identity harness, which is what made all of this safe.** Twenty configurations run from a fixed seed must match a stored baseline bit for bit. It caught ternaries blocking fused multiply-adds, isolated the split-less-tree Gibbs draw, and allowed the `split_support` rewrite in one step.

- **by 2026-08-24: Earlier performance work, for the record.** A first pass gave 1.6x over five families: score and information in one pass, analytic derivatives for `ordinal()` and `ordbeta()`, hoisted eta-free terms, a binomial rewrite (probit 8.66 s to 2.63 s), a Fisher-scoring tolerance of 1e-2, Robbins-Monro bandwidth tuning, and reusable buffers. `test-derivatives.R` checks every analytic score; `ordbeta()` remains the most expensive family.

- **by 2026-08-24: The blocked evaluation path, and why it is opt-in.** A blocked-only path cost 30% on a Gaussian fit, so `Family::wants_block()` selects it only for the two R-backed families. `block_eval` forces it and `test-flexible.R` runs four families both ways from one seed.

- **2026-08-25: The exponential form under soft rules: measured, and the ceiling is 5–10%.** A scratch build forcing `exponential_usable()` measured the upper bound at 1.05x for Poisson and 1.10x for gamma under soft rules, since the first O(n) pass dominates. Dropped from the To Do list.

- **2026-08-25: The zero-inflated and multinomial exponential route: still open, and what it is for.** Murray's multinomial-Poisson route would target the Polya-Gamma route's 0.37x mixing rather than speed, but changes identification and reaches `category_probs()`, the coding option and `before_forest()`. The zero-inflated families would first need a structural-zero indicator augmentation.

- **2026-08-30: Speeding up the survival models.** `exp_sign()` became `exp_rate()`, putting `location_scale()`'s scale forest (7.7 to 3.6 s) and `weibull_aft()` (6.7 to 4.2 s, already exponential at rate `-1 / sigma`) on the exponential form under hard rules. Imputation augmentations that survive soft rules gave `lognormal_aft()` 14 to 16x (soft) and 20 to 22x (hard) ESS per second and `loglogistic_aft()` 4.9 to 7.7x and 9.8 to 11x. The between-sampler posterior gap was smaller than each sampler's own Monte Carlo error.

- **2026-09-24: Prior-only fits ran on the slow path.** `prior_only = TRUE` refused every augmentation and evaluated zero-weighted densities; skipping zero-weight observations in `src/family.h` took prior-only logit from 10.89 to 2.66 s (soft) and 6.35 to 1.36 s (hard), bit-identical (`_dev/zero-weight-skip.R`). Hard prior-only fits remain about 1.5 times the posterior, and a +5.2% on `poisson()` was not settled.

### Benchmarks against other packages

**The benchmark measured the wrong thing twice (a `dpm()` fit timed as the Gaussian rows, and the other packages' prediction inside their timed calls while bartisan's ran outside), fixed on 2026-09-02; the corrected table is the one in Status.**

- **by 2026-08-24: Benchmarking against bartMachine and stochtree.** `bartMachine` was silently absent: it needs `options(java.parameters = ...)` before the JVM starts and `JAVA_HOME` set in `~/.Renviron`. `_dev/benchmark.Rmd` now detects backends by fitting something tiny and lists failed fits instead of dropping them; stochtree is timed with `num_gfr` at 0 and 5.

- **2026-09-02: The benchmark measured the wrong thing twice, and a progress bar.** The Gaussian rows timed `dpm()` (worth 1.35x) and competitors predicted inside their timed calls, so the matched gap to dbarts is 2.87, not 1.8; dbarts's partition and cached-summary techniques are unavailable under soft rules and were not attempted. A C++ audit found no correctness bugs (two shadowed names renamed). A *progressr* progress bar (`control$progress`, `R/progress.R`) caps at 50 reports per chain and consumes no randomness; its multisession test is meaningful only after `R CMD INSTALL`.

### Dispersion: quasi-likelihood and the bootstrap

**Routes to better-calibrated dispersion: a Bayesian bootstrap with one chain per replicate is coherent but conservative, and quasi-likelihood is a clean extension with a caveat; neither is built, and the bootstrap dispersion draw and a `quasi()` family are on To Do.**

- **by 2026-08-24: Quasi-likelihood (Linero 2026): a clean extension, with a caveat.** The `Family` interface covers it, with `compute_eta_free()` left empty and `phi` drawn from Pearson residuals, at the cost of the exactness claim. On inverse-gamma data `Gamma("log")` underestimated the dispersion (54% to 22% of truth) with coverage 0.884 to 0.900, making a Bayesian-bootstrap dispersion draw in `update_aux` the cheapest route to calibration; it offers no robustness to a misspecified variance function.

- **by 2026-08-24: Bootstrapping: one chain per replicate.** Dirichlet weights give a Bayesian bootstrap with no package change. Summarizing each replicate by its posterior mean covered 0.891, worse than the posterior's 0.927, while one draw per replicate covered 0.985 with intervals 47% wider: coherent as a nonparametric posterior, but conservative.

### Design decisions

**The family-specific surface of the sampler reduces to three per-observation quantities (log density, its derivative and the information), which is what let the later families (zero-inflated, ordered beta, conditional density) in without changes to the sampler; the departures from Linero's reference implementation are recorded here.**

- **by 2026-08-24: Design decisions.** Families supply `logdens`, `dlogdens` and `info`; soft rules enter by the chain rule (factors `w_i` and `w_i^2`) on one code path, and `Family::target_form()` declares a quadratic or exponential shape. Four departures from Linero's `FlexBart`: a bivariate Laplace proposal for soft rules, deterministic Laplace fits starting at zero, a corrected death-move ratio, and capped slice-sampler loops.

- **by 2026-08-24: Later families.** Zero-inflated Poisson and negative binomial use two forests with complete-data expected information (named `zi_poisson()` to avoid masking `utils::zip()`), ordered beta draws both cutpoints, and `predict(type = "density")` calls the C++ log density through a `log_norm_const()` hook.

## Mixing, diagnostics and parallelism

### Mixing: chains against draws

**Whether many chains or one long chain mix a fit better: the chain-count effect was Monte Carlo error, chains and draws are not interchangeable for what `diagnose()` grades, calibration does not degrade with n and effective sample size falls with n as Tan et al. (2026) predict but without a multimodality signature, and since 2026-11-01 `add_chains()` adds chains to an existing fit.**

- **2026-09-05: Log: stage 1 of the mixing plan -- the predicted n-trend is not here.** `_dev/sbc.R` at n = 250, 1000 and 4000, hard and soft gates, tested whether the posterior approximation degrades with n: no cell departed from uniform (smallest p 0.37), coverage ran 0.925 to 0.960, and no rank histogram showed a U. Stages 2 and 3 were not built and the soft-gate hypothesis in `_dev/MIXING.md` was not supported; n past 4000 and a Kim and Ročková deep-signal generator remain untested.

- **2026-09-05: Chains and draws are not interchangeable, and lalonde is the worked example.** `_dev/chains-vs-draws.R`: average ESS ignores how 16000 draws are split while worst-5% R-hat rises with chain count (1.403 at 2, 1.815 at 16), the null moving. In `_dev/lalonde-reproduce.R`, 6 x 10000 warned and 4 x 30000 passed in 3 of 3 seeds; 12 x 10000 had the same ESS but warned with every R-hat at its null and took 256 s against 134 s, so more chains is cheaper only up to the worker count.

- **2026-09-13: The reported quantity mixes worse than anything the fit's table shows.** The averaged predictor had R-hat 1.004 and 1610 effective draws while the ATE had 1.04 and 117, one chain sitting at exactly zero for 663 draws because the splitting prior gave `rhc` no rule. `diagnose()` became generic with a `<bartisan_effect>` method adding a share-at-the-atom check, `vignette("diagnostics")` shows `posterior::summarise_draws()`, and two claims that the average row governs an estimand were withdrawn.

- **2026-09-14: Why ten thousand draws yields a small ESS.** Bulk ESS per draw fell from 6.9% at 1,000 draws per chain to 3.0% at 10,000, because a short chain cannot see long-lag autocorrelation and overstates ESS; absolute ESS still grew. R-hat was 1.003, so the chain mixes slowly rather than failing to mix, and more draws is the remedy.

- **2026-09-14: Getting to clean mixing on draws alone.** The 14-predictor model plateaued at R-hat 1.019 and an `n = 800` subset mixed worse; the five-predictor model cleared every warning at 2,000/8,000 draws (R-hat 1.008, min bulk ESS 744), ESS clearing before R-hat. The vignette notes ESS is non-monotone in draws and that draws help only when the flagged R-hat is the unreadable kind; `diagnostics.Rmd` went from 110s to 194s.

- **2026-09-19: Many chains against one long chain: a large effect, and it is an artifact.** `_dev/chains-vs-length.R` (Friedman, n = 4000, 25 of 30 predictors irrelevant, equal sweeps) tested whether between-chain diversity widens intervals as in Krantsevich, He and Hahn: coverage rose from 0.949 with one chain to 0.994 with 16 by intervals 36% wider. No arm met `rhat_max = 1.01`, and the next entry showed it is Monte Carlo error, so the table is not a property of chain count.

- **2026-09-19: The chain-count effect is Monte Carlo error, measured.** `_dev/chains-convergence.R` thinned at fixed stored draws: the one-minus-four-chain width gap decayed with exponent −0.44 (sparsity on) and −0.47 (off) against the Monte Carlo −0.5. `sparsity = FALSE` made the gap two to three times larger, so the propensity entry's 2x mixing gain from it does not generalize to designs with irrelevant predictors; 51200 sweeps gave coverage 0.9575 without reaching R-hat 1.01, and `chains = 1` stays.

- **2026-09-19: What this says about the defaults and about `diagnose()`.** The default 800 sweeps overstated interval width by about 5% (coverage 0.938 against 0.958 at 51200) for 1/64 of the compute; the default stays, and `?bartisan_control` should say intervals want more sweeps. `diagnose()` was not failing: R-hat cannot see bias, every arm rightly failed 1.01, and the one-against-four disagreement was the sharpest convergence statement. The claim that `splits.eta` is the row to act on was measured in the next entry not to hold.

- **2026-09-19: Mixing against sample size, against Tan et al. (2026), and the row that is actually worst.** `_dev/mixing-against-n.R` (144 fits, three families, n 500 to 8000) found R-hat rising with n but every cell at or below its null `1 + chains/ess` while ESS fell, so mixing degrades as Tan et al. predict without detectable multimodality; soft rules mixed better in eight of nine cells. The splits row was worst in only 8 of 144 fits, so `?diagnose` was corrected; `aux.cut1` was the worst-mixing quantity (R-hat 1.31 on 23 effective draws), chased in the next entry.

- **2026-11-01: `add_chains()`, and the session stream that was never advanced.** `add_chains(object, chains = 1)` refits through `refit_call()` and pools through `combine_fits()`; `bartisan()` keeps `fit$next_stream`, so added chains equal a larger `chains` (`test-add-chains.R`) and a single chain now runs on stream 1. Bug fixed: `parallel_streams()` restored `.Random.seed` exactly, so consecutive unseeded multi-chain fits were identical. `vignette("diagnostics")` was refreshed (prior fit now four chains), other saved-results vignettes were not, and `NEWS.md` under 0.1.0.9000 carries three items.

### The burn-in transient and warm starts

**Warm starts (grow-from-root, random tree features, XBART-style) were proposed to shorten burn-in, but the transient is 34 to 70 sweeps and does not grow with n, so both items were struck and the default warmup was cut instead.**

- **by 2026-08-24: XBART (He and Hahn 2021), and XBCF (Krantsevich, He and Hahn 2023).** Not adopted: grow-from-root is hard-rule only (soft rules are the default), not a posterior sampler, `Node::draw_rule()` has no cutpoint grid, and `encode_tree()` has no decoder. An informed-proposal birth move would stay exact but is hard-rules-only and not measured; Tan et al. (2026) found warm-starting from a fitted ensemble worsens R-hat, and their tempered tree-acceptance schedule (3 to 1) is the one idea worth costing out, not implemented.

- **2026-09-04: Assessment: the two warm-start items, and what to do instead.** The transient does not grow with n (`_dev/transient-scaling.R`: at sweep 36 hard rules stood at 0.86, 0.81 and 0.84 of their rise across n = 1500 to 50000), so a warm start could save at most about 5%. Trap: the "within two standard deviations" statistic is not scale-free; use the fraction-of-rise table. Random features reached R-squared 0.876 in 0.25 s against BART's 0.975 in 14.0 s, a case for a `random_features()` estimator; both items were struck and `num_burn` went from 500 to 200.

### `diagnose()` and the convergence diagnostics

**`diagnose()` is the single route to convergence diagnostics (rank-normalized R-hat and effective sample size, run over a `future` plan); a drift statistic was abandoned, the pass was made 35% to 40% faster, R-hat's null of `1 + chains / ess` informs the thresholds, and two small errors in its R-hat and effective sample size were fixed against *posterior* on 2026-10-02.**

- **by 2026-08-24: Rank-normalized diagnostics.** `fit$rhat` reported rank-normalized folded split R-hat and bulk and tail ESS (Vehtari et al. 2021), validated by hand, against *posterior* (within 1.6%) and against AR(1). A constant quantity gave R-hat `-Inf` and a fabricated ESS, now `NA`, tested on a two-category ordinal fit. Lesson: compare implementations on inputs whose answer can be worked out by hand; an apparent disagreement with *posterior* was the harness.

- **2026-09-02: `diagnose()`, and a drift statistic that had to be abandoned.** `diagnose()` in `R/diagnose.R` works with one chain by splitting it, keys checks to the share of failing components, and repeats R-hat on the second half to separate short warmup from chains settling apart. A Geweke-style drift statistic was dropped because no variance estimator calibrated on AR(1); the measurements are in the roxygen. `README.Rmd` gained a feature comparison.

- **2026-09-03: Log: the diagnostics moved to diagnose(), and the pass got parallel.** `fit$rhat` and `chain_diagnostics()` were removed so the pass runs only in `diagnose()`, taking a four-chain fit at n = 2000 from 7.9 s to 3.34 s; `vignette("diagnostics")` was rebuilt around `diagnose()`. `diagnosis_columns()` chunks over a `future` plan, differing across worker counts only in the last bits of ESS (FFT); the speedup was unmeasured, `availableCores()` being 1.

- **2026-09-03: Log: why the diagnostics pass cost more than the fit.** The cost was R overhead in the per-observation ESS loop: one `stats::acf()` call per half-chain per observation, and Geyer sequences running all 98 lags. An FFT `autocovariance()` and a preallocated Geyer loop took `chain_diagnostics()` from 4.85 s to 1.84 s. A bit-identical claim was withdrawn because the test shim never took effect; FFT agrees with `acf` to 1.8e-15.

- **2026-09-04: Log: the convergence pass, 35% to 40% faster, and the chunking that does not help.** Chunking finer than `nbrOfWorkers()` was slower at every multiple tried. Rankings per column went from seven to four (tail ESS needs none, R-hat and bulk ESS share one), 39.8% faster on the min at n = 4000. Trap: pueue's default group runs four tasks on four performance cores, skewing concurrent timings by up to 20%, and an `R CMD INSTALL` truncated a test run; time one job at a time with paired in-session comparisons.

- **2026-09-05: R-hat's null is `1 + chains / ess`, and what that does to the thresholds.** `_dev/ess-rhat-calibration.R` found R-hat's null mean is `1 + m/S`, so at 400 effective draws over four chains it is 1.010 and `rhat_max = 1.01` with `ess_min = 400` is one threshold stated twice. A `rhat readable` check against `chains / (rhat_max - 1)` now sends advice to `num_draws`. Not fixed: `FAIL_SHARE = 0.2` is not calibrated at the boundary, where 44% of stationary columns exceed 1.01.

- **2026-09-05: What `diagnose()`'s worst row is actually measuring.** The predictor now gets a worst-5% row and an average row, because `_dev/where-mixing-lives.R` found the average carries about one effective draw per kept draw while observations carry tens (Friedman: 7754 against 21, with worst-5% R-hat below its null). Only `dpm()` exceeded its null (1.789 against 1.667).

- **2026-09-14: ess_tail above ess_bulk is the normal shape here.** Systematic on `rhc` (`loglik` 23 bulk, 234 tail): a drifting level is autocorrelated in the centre while tail indicators flip rarely, so it says the centre drifts, not that the tails are good. The ATE inverts it (bulk 117, tail 64) because the atom at zero is the lower tail.

- **2026-10-02: `diagnose()` computed R-hat and the effective sample size slightly wrong.** A seed-shifted test gave R-hat 1.352 against *posterior*'s 1.366. `rank_normalize()` divided by `n - 1/4` instead of Blom's `n + 1/4` (since `cdc3278`), and `ess_from_split()` counted the last Geyer pair's even term twice; fixed, they agree to about 1e-15 apart from a documented odd-draws convention. `test-interop.R` now compares at 1e-10 and `vignette("diagnostics")` was regenerated.

### Parallel chains, futures and memory

**Chains run in parallel over `future` without future.apply; the parallel pass is limited by this machine's four performance cores (first misattributed to memory bandwidth), the unit maps were carrying the design matrix twice, `partial_dependence()` is the only other place parallelism pays, and the 500 MB `future.globals.maxSize` limit a large fit meets is now documented.**

- **by 2026-08-24: Parallel chains.** Chains run through `future.apply::future_lapply()`, so any `future` backend works, and one `set.seed()` reproduces bitwise across backends. Four chains on `lalonde` took 10.14 s sequentially and 3.50 s on four workers. Pooling shifts each chain's `tree_start` offsets, tested by `expect_predictor_invariant()`.

- **by 2026-08-24: Multithreading below the chain level.** Rejected: per-leaf loops are the size of synchronization overhead, backfitting is sequential across trees, and families call R's math library, which is not documented as thread-safe.

- **2026-08-31: More than one chain no longer needs future.apply.** Chains run sequentially without it; `parallel_streams()` advances one L'Ecuyer stream per chain from the session seed so both branches give identical draws, `parallel` moved to Imports, and the session generator is restored afterwards.

- **2026-09-03: Log: what limits the parallel pass is memory bandwidth (a correction).** `_dev/diagnose-timing.R` showed the pass is about 78% of a four-chain fit's time and caps near 3x. Slicing blocks so `wide` (102.4 MB per worker) was not a global changed nothing: a size from `getGlobalsAndPackages()` is not a time until timed. The ceiling is the M4's four performance cores, not bandwidth; the slicing was reverted and the `columns < 400` floor lowered to 100.

- **2026-09-04: Log: what four chains on four cores actually cost.** Four chains on four workers cost 1.22x (n = 2000) and 1.28x (n = 8000) one chain's wall clock, the four-fast-core ceiling rather than overhead. Eight workers was slower than four, so `workers` should match the chain count; `multicore` saves only startup. Faster fits need a faster single chain, such as `rebuild_support()`.

- **2026-09-05: The unit maps were carrying the design matrix, twice.** `make_unit_map()` held the design matrix twice through its frame and an unforced `type` promise; constructors (`constant_map`, `midpoint_map`, `range_map`, `ecdf_map`) now force what they capture. On RHC the maps went from 788 KB to 84 KB (quantile) and 721 KB to 5.6 KB (range); `predict()` is bit-identical, and `test-invariants.R` has the test.

- **2026-09-14: `diagnose()` died on multiple workers, and the error had been dismissed once.** `could not find function "diagnosis_block"` arose because `future` inferred `packages = "bartisan"` and workers attach only exports; fixed by naming it in `globals`, with a two-worker test. An earlier entry had wrongly dismissed it as a `load_all()` artifact. A `future_lapply()` rewrite would have prevented it; its effect on progress relay was not measured.

- **2026-09-14: Parallelism where it pays, and the limit nobody had written down.** Only `partial_dependence()`'s grid loop qualified (2.9x at 10 points on four workers) and now runs in parallel with no floor; other loops were rejected with reasons. A fit at 1500 observations and 3200 draws is 71.8 MB, and `future.globals.maxSize`'s 500 MB limit, which affects chains, the convergence pass and prediction, is now in `?bartisan_control`.

### Progress bars

**A run of progress-bar defects (finishing before the fit, two bars in `bcf()`, a 336 MB closure behind `diagnose()`'s bar, jumps when workers ran in step) were fixed one at a time, and since 2026-09-17 every `_dev/` simulation script reports through the live-progress widget.**

- **2026-09-03: Log: the progress bar finished before the fit did.** Relaying was correct; the bar did not cover `chain_diagnostics()`, so on four workers it hit 100% at 18% of the run. `PROGRESS_DIAG_TICKS` of 50 now gives the pass a share. Lesson: *progressr* reports nothing non-interactively unless `progressr.enable` is set, so validate the instrument first.

- **2026-09-03: Log: `bcf()` showed two progress bars.** `bcf(propensity = TRUE)` built one progressor per `bartisan()` call. One bar is now shared through `the$claimed_progress` and `shared_reporter()`, sized to the smallest tick count across fits, and the claim covers the adaptive retry and is released on error.

- **2026-09-05: The `diagnose()` progress bar, and the 336 MB closure behind it.** The bar stalled because the session cannot relay while dispatching futures, and `diagnosis_reporter()`'s closure serialized to 336 MB through `envir`, which held the fit. `reporter_for()` forces `p` in an empty frame and dispatch polls running futures (first tick 0.78s to 0.21s). Trap: `parts[[k]]` inside a future exports all of `parts`; `test-invariants.R` caps reporter size.

- **2026-09-14: The bar jumped because the workers were in step.** Jumps came from 50 `PROGRESS_DIAG_TICKS` split among in-step workers; 200 ticks cut the largest gap from 2.8% to 1.3% of the run at no cost. Every earlier "zero progress events" measurement was worthless because *progressr* handlers are disabled in non-interactive R.

- **2026-09-17: Live-progress instrumentation across every simulation script.** 48 of 52 `_dev/` scripts now have `prog_init()`, `prog_tick()`, an `on.exit(prog_end(...))` and a `checkpoint()` inside the loop with a `complete` flag. Pattern-matching audits missed `lapply()` loops and `do.call(bartisan, ...)` fits, and the threshold is 8 seconds, not a sampler. Two edit errors (undefined `n_total`, a deleted `out <- list()`) were caught by re-reading, not by parsing.

## Families

### `dpm()`: the Dirichlet process mixture error

**`dpm()`, a Dirichlet process mixture for the error, is the default for a numeric response and reports the conditional mean; its level in a `bcf()` fit was an unidentified coordinate (fixed 2026-10-01), and since 2026-11-03 the mixture's weights and its baseline's k_0 and lambda are drawn, which brings its intervals level with other packages' mixtures or above them.**

- **2026-08-30: DPMBART: a Dirichlet process mixture for the error distribution.** `dpm()` follows George et al. (2019) with a conjugate baseline calibrated off a linear fit; the target stays `TARGET_QUADRATIC` at about 1.3x a Gaussian fit's cost, and the mixture is stored as `mixture_flat`/`mixture_start`. Only the sum of the fit and the error mean is identified (correlation -0.987), so `type = "response"` is the conditional mean. It tied `gaussian()` on normal errors, won by 175 and 224 log points on t3 and skewed errors, and lost heteroskedasticity to `location_scale()` (119 log points), which costs about nine times a Gaussian fit.

- **2026-08-30: The family is inferred from the response when none is named.** `family = NULL` maps `Surv` to `weibull_aft()`, an ordered factor to `ordinal()`, binary data to `binomial()`, other factors to `multinomial()`, and anything else to `gaussian()`, with a message. Counts are deliberately not read as `poisson()`, and `c(1, 2)` stays Gaussian. `is_binary()` now requires exactly two values, so a constant response errors rather than fitting a binomial.

- **2026-08-30: gaussian(), dpm() and ordinal() on a numeric response: dpm dominates.** At n = 1000 `dpm()` was level with `gaussian()` on normal errors and ahead by 178, 196 and 582 log points on t3, skewed and bimodal errors; `gaussian()` keeps prior weights (which `dpm()` refuses), an identified predictor and speed. `ordinal()` on a continuous outcome is `rms::orm()` with a forest and should be binned: 25 bins was sixteen times faster than none and slightly more accurate, with 10 to 50 bins all fine. Documented in `vignette("families")`, the help page, README and NEWS.

- **2026-08-30: dpm reports the conditional mean on the predictor, and is the numeric default.** `DPMFamily::report_shift()` records draws in the chart where the mixture has mean zero (predictor level sd 0.750 to 0.021); the baseline term lives on the raw chart, so `aux` now reports `center` in place of `error_mean`. A numeric response defaults to `dpm()` unless it has fewer than ten distinct values (3 levels gave an RMSE ratio of 3.78); clamped scales want `ordinal()`, and a weighted fit with no family named is an error.

- **2026-09-03: `error_density(plot = TRUE)`.** Returns a `ggplot` of the posterior mean density with a pointwise ribbon (*ggplot2* is a soft dependency); `plot = FALSE` still returns the data frame, and `vignette("survival")` keeps its hand-drawn overlay.

- **2026-10-01: The level of a `bcf(family = dpm())` fit was the unidentified coordinate.** `VaryingCoefficientFamily` did not forward `report_shift()`, so `bcf()` and `vc()` fits recorded the raw forest (level tracking `center` at -0.997); fixed in `src/family.cpp`, tested in `test-dpm.R` and `test-ordinal.R` (`_dev/dpm-level-shift.R`). On lalonde ATT draws were unchanged and the spread of Y[1] across fits fell from 541 to 44, leaving it about 300 above the observed 6349 as a model gap. `_dev/dpm-diagnose-rows.R` showed slow mixture coordinates reach the error density at ±1 sd (R-hat 1.2); dropping `aux.center` from `diagnose()` and adding density rows was proposed, not done.

- **2026-10-01: The level of a `dpm()` fit mixes slowly (closed To Do item).** The lalonde `bcf(family = dpm())` fit gave treated `Y[0]` between 3780 and 5650 depending on chains against an observed 6349, read as slow mixing of the level; closed the same day as the varying-coefficient wrapper not forwarding `report_shift()`.

- **2026-11-03: The mixture's weights drawn in the engine, `k_s` checked, and the sources.** `DPMFamily::draw_report_weights()` draws the weights ~ Dirichlet(n_1, ..., n_K, alpha) for every reported quantity (Gelfand and Kottas 2002), with the remainder still at the baseline's t; DPM random streams change (NEWS item). The `@param k_s` text had the direction backwards and was corrected against George et al. (2019). In `_dev/survival-skew-coverage.R` bimodal coverage rose (.55 to .70 `dpm_aft()`, .68 to .79 `dpm_aft_ls()`) but stayed short of .90, and skewed did not move; nftbart weights by count shares but draws k0 and b0.

- **2026-11-03: Drawing k_0 rather than fixing it.** Internal options `k0_shape`/`k0_rate` make `DPMFamily::update_baseline_scale()` draw k_0 conjugately; unset, fits match the previous build bit for bit. In `_dev/dpm-k0-draw.R` (160 fits) bimodal coverage under `centered` rose to .81 and .85 from .70 and .79, with k_0's posterior falling to .11 to .19, but `dpm_aft_ls()` lost skewed coverage (-.028, SE .009), which the rule written beforehand puts in the keep-fixed branch; `lio` was worse. Adoption was left to the maintainer.

- **2026-11-03: The mixture alone: drawing k_0 and lambda, with the mean known.** `_dev/dpm-k0-oracle.R` (true mean as offset) added `lambda_shape`: drawing both raised bimodal error-density coverage from .34 to .80 and, against the prediction, skewed from .62 to .89, so the variance prior's floor was a cause; log-normal lost nothing and nu = 3 added little. Skewed S(t | x) coverage is partly an artifact, since 25% of cells have a true S(t | x) of exactly one, so scripts now record `cover_inner` (`_dev/skewed-inner-coverage.R`).

- **2026-11-03: Drawing k_0 and lambda on fresh replicates, and `dpm_aft()` against `lognormal_aft()`.** `_dev/dpm-k0-reps.R` (unseen replicates 6 to 25): `both` raised bimodal coverage to .856 for `dpm_aft()` (from .707) and .875 for `dpm_aft_ls()` (from .801) at one small cost; the earlier skewed cost did not replicate and nu = 3 was rejected. Under the maintainer's rule (`dpm_aft()` stays default while covering at least as well as `lognormal_aft()`), it covered as well or better on five truths and less on log-normal (.928 against .947) and spread varies (.511 against .537); skewed was pending `cover_inner`.

- **2026-11-03: Uncensored: `dpm()` against `gaussian()`, and what drawing k_0 and lambda costs the mean.** `_dev/dpm-k0-numeric.R`: default `dpm()` ties `gaussian()` on the conditional mean, covers the error density and P(Y > c | x) far better on non-normal errors and pays a little on normal ones. Drawing k_0 and lambda raised density coverage (bimodal .603 to .919) but cut conditional-mean coverage beyond two SE in four cells, `dpm_ls()` skewed from .924 to .881, below `gaussian_ls()`'s .945, so by the prior rule `both` could not be default. Untested: lambda drawn with k_0 fixed, or a tighter k_0 prior.

- **2026-11-03: The baseline's k_0 and lambda are drawn by default.** At the maintainer's direction the four DP families gained `update_baseline = TRUE` (gamma priors of shape 1.5 centered on the George et al. values, in `DPMFamily::update_baseline_scale()`); `FALSE` restores the fixed baseline. The known cost, lower conditional-mean coverage (`dpm_ls()` skewed .924 to .881), was weighed against the survival and density gains and is named in `?bartisan-families`. Vignette numbers and precomputed vignettes with a DP family are stale and need rerunning, coordinated with another session.

- **2026-11-03: Other packages' Dirichlet process mixtures, and whether drawing the baseline is what covers.** `_dev/dpm-packages-coverage.R`: BNPmix and *dirichletprocess* covered the bimodal error density with fixed baselines (.91), and `anovaDDP()` covered bimodal S(t | x) at .892 against .721 for fixed-baseline `dpm_aft()` and .947 drawn. So fixing a baseline is not the problem; bartisan's tight k_0 (k_s = 10) and lambda at nu = 10 were the outlier. Chattopadhyay et al.'s under-coverage for *dirichletprocess* reproduced only on the skewed error (.725).

### Location-scale families

**Location-scale families give the spread a forest of its own (`gaussian_ls()`, `Gamma_ls()`, `dpm_ls()`, and the AFT versions); `dpm_ls()` scales the whole error, and whether its scale forest is needed can be decided by cross-validation, which recovers a constant spread more often than PSIS-LOO does.**

- **2026-09-03: `location_scale()` became `gaussian_ls()`.** Renamed everywhere including the user-visible engine family string (28 files; C++ class `GaussianLSFamily`); older entries keep the old name.

- **2026-09-03: `Gamma_ls()`: a gamma whose dispersion is a forest.** A second forest on the log dispersion (shape `exp(-eta1)`), the mean in exponential form and the dispersion on the general path with the expected information `s^2 trigamma(s) - s`, which is always positive. It costs about 5.4x `Gamma("log")`; the log density matches `dgamma()` to 7.1e-15, and a varying dispersion was recovered at correlation .973, worth 62.9 log points.

- **2026-11-01: `dpm_ls()`: the whole error scaled, or each component.** In `_dev/dpm-ls-forms.R` the whole-error form tied or won on four of six scenarios, lost on `component` as predicted and, unpredicted, on `homoskedastic` (RMSE .0900 against `dpm()`'s .0787), because spurious scale-forest variation moves the mean under a skewed W. `_dev/dpm-ls-scale-prior.R` found no scale-forest prior that fixed it; a fixed small leaf scale was not tried. At the maintainer's direction the trade-off is documented rather than pursued.

- **2026-11-01: Heteroskedastic `dpm()` and `dpm_aft()`, which with censoring is NFT BART (closed To Do item).** `dpm_ls()` fits y = f(x) + exp(g(x)) W with W a DPM of mean 0 and variance 1, and `dpm_aft_ls()` does the same on log T with censoring (Sparapani et al. 2023), reported through `report_shift_obs()` and `report_chart()`. The component form was removed because nftbart fits the whole-error model, and `dpm()`/`dpm_aft()` match 797726d bit for bit. `dpm_aft_ls()` was not then measured against nftbart or `dpm_aft()`.

- **2026-11-02: The scale switch inside the sampler.** An experimental reversible-jump switch on the scale forest was calibrated on the prior, but for `dpm_ls()` the approximate, exact blocked and joint-move versions made no switches on data (off moves' log ratio averaged -732), so its 46% recovery was not a posterior average. All of it was removed at the maintainer's direction (`src/mcmc.cpp` back to 797726d; `_dev/scale-switch-sampler.R` untracked). Cross-validated averaging (`_dev/scale-switch.R`) is the route that worked.

- **2026-11-02: The new families in the vignettes, and two reruns.** `_dev/error-shapes-sim.R` with `dpm_ls()` added did not reproduce the published table, so `vignette("families")` § Numeric Responses was rebuilt from it; `dpm_ls()` matched `gaussian_ls()` on heteroskedastic errors. The companion `_dev/survival-sim-ls.R` paired across builds and was superseded the same day by `_dev/survival-sim.R`. The maintainer decided defaults stay single-forest for every response type.

- **2026-11-02: Switching `dpm_ls()`'s scale forest off by cross-validation.** `_dev/scale-switch.R` averaged full and flat (`log_sd = ~ 1`) fits: K-fold stacking recovered 80% of the full fit's constant-spread RMSE loss (.1014 to .0904 against flat .0876) and cost nothing where the spread varies, narrowly failing the prior rule (2.5 SE above flat). PSIS-LOO recovered 56%.

- **2026-11-02: A location-scale log-normal accelerated failure time family (closed To Do item).** `lognormal_aft_ls()` puts a second forest on the log scale of `log T` so curves can cross; `LognormalAFTLSAugmentedFamily` imputes censored times when all weights are one and `augment` is on, and `LognormalAFTLSFamily` uses the direct censored likelihood for prior weights. Tested in `test-aft-ls.R`, `test-derivatives.R` and the invariants matrix; not yet measured in a simulation.

### Positive outcomes: `Gamma()` and its rivals

**Whether `Gamma()` earns its place on positive outcomes: it stays, masked so that its default link is log, `Gamma_shape()` was removed, and the comparison was re-measured on 2026-09-03 against like-for-like rivals after the first one compared quantities on different scales.**

- **2026-08-30: Is `Gamma()` worth keeping?** A narrow yes: `Gamma("log")` alone had the best predictive density when the gamma was the truth, but never the best RMSE and was 4 to 6x slower. `Gamma("inverse")` was worse on every measure and produced NaN densities when draws went non-positive (To Do). `ordinal("probit")` on 25 bins was within 3% of the best RMSE everywhere, and `dpm()` on the raw scale of a skewed positive outcome had about twice the RMSE, so the advice to prefer it does not carry there.

- **2026-08-30: The gamma link, decided rather than reported.** Any `Gamma()` link other than `log` is now replaced by `log` with a message; `stats::Gamma()` is untouched. The inverse link is canonical but lets the unconstrained predictor reach negative means (NaN densities in five of eight heavy-tailed replicates) and is slower; composed gamma links are no longer reachable, `Gamma` left `native_links`, and `custom_family()` is the route for one.

- **2026-08-30: `Gamma_shape()` removed, `Gamma()` masked so that the default link is log.** Removing `Gamma_shape()` exposed that `stats::Gamma()` defaults to the inverse link (RMSE 0.664 at 7.2 s against 0.606 at 3.8 s for log), so the package exports `Gamma()` with `link = "log"`, masking `stats::Gamma()` (documented). `warn_restricted_link()` now warns when a composed link's inverse does not cover the predictor. A test regex on `"using"` silently stopped matching after sentence-casing and now matches `"family = "`.

- **2026-09-03: The positive-outcome comparison, re-measured against the right rivals.** `_dev/positive-sim.R` replaces the old table, which compared log E[Y] with E[log Y]; `gaussian("log")` is the like-for-like Gaussian, and `dpm("log")` silently passes `"log"` to `nu`. `Gamma_ls()` cut RMSE 26% and gained 62 log points where the dispersion varies and was two log points from `Gamma("log")` where it does not. `gaussian("log")` trailed by 150 to 310 log points and was slower than `Gamma("log")` because R-side links run once per leaf.

### The beta families

**The beta families were the slowest in the package; a reworked target made `Beta()` and `ordbeta()` 4.2 to 4.6 times faster, and no augmentation exists for them because the obstacle is the normalizing constant's dependence on eta.**

- **2026-08-30: `Beta()`, and what it is actually worth.** `Beta()` is `OrdBetaFamily`'s interior without the endpoints, 11% better than `ordbeta()` at n = 200 and a tie at n = 500; its case is that `ordbeta()` on interior data leaves two unidentified cutpoints (reaching 18 to 19), so choose by whether the response can reach a boundary. With 1% at a boundary `ordbeta()` beat a clamped `Beta()` by 5% and `gaussian()` on the logit by 30%.

- **2026-08-30: Making the beta families fast.** The cost was the `digamma`/`trigamma` calls, which within a sweep depend only on `mu`, so they are tabulated on a 2049-point grid rebuilt per sweep: 4.6x (`Beta()`) and 4.2x (`ordbeta()`), legitimate because only the Laplace proposal uses them. An asymptotic `logit(y)` approximation cut ESS from 241 to 91 and was rejected; `ordbeta()`'s cutpoint updates now evaluate only their half of the likelihood (bit-identical). No exact augmentation exists, because eta enters the normalizing constant.

### `tweedie()`

**`tweedie()` covers a non-negative response with a point mass at zero, and its intractable density turned out not to matter for sampling; the Jorgensen and Dunn and Smyth references are still wanted.**

- **2026-09-08: `tweedie()`, and why the intractable density did not matter.** In exponential-dispersion form the series term `log W` has no `mu`, so it lives in `compute_eta_free()` once per sweep and the fit costs about what `Gamma("log")` does. `power` is validated on the open (1, 2) and fixed at 1.5 by default (`power = NULL` draws it); `test-tweedie.R` matches a direct compound Poisson-gamma sum to 1e-8. Lesson: an `expect_error()` regexp should pin the rejected argument, not the wording. References are still wanted.

### The zero-inflated families

**The zero-inflated families use two latent variables for a 4 to 10 times gain in effective draws per second, but their estimands still need about 3,300 (`zi_poisson()`) to 13,000 (`zi_negbin()`) draws for 400 effective, which `vignette("families")` now says.**

- **2026-08-30: The zero-inflated families: two latent variables, 4 to 10x in ESS/s.** A structural-zero indicator `z` separates the count forest (Poisson exponential form or negative binomial) from the inflation forest (Polya-Gamma, quadratic), and `zi_negbin()` adds the Poisson-gamma rate in a partially collapsed Gibbs step. ESS/s gains were 3.9x to 10.1x across families and rules, larger than predicted because the direct log-sum-exp target was expensive; on by default for both rules.

- **2026-10-02: The families vignette says how many draws a zero-inflated estimand needs.** `_dev/zi-draws.R`: a treatment contrast gets 122 to 144 effective draws per 1000 for `zi_poisson()` (about 3300 for 400) and 27 to 54 for `zi_negbin()` (about 13000), flat in length; its attenuation (0.68 to 0.75 against 0.8) is the prior. `vignette("families")` now says so.

### The hurdle families

**`hurdle_poisson()`, `hurdle_negbin()`, `hurdle_Gamma()` and `hurdle_lognormal()` give the probability of a zero a forest of its own and fit the positive part separately (a count truncated at zero, or a gamma or log-normal amount); their augmented sampler matches the direct one on coverage, its count intervals match the direct sampler's at four long chains for the Poisson and come within 1.9% for the negative binomial, and with the tree moves taken out it matches the exact posterior of the count predictor and theta to within Monte Carlo error, so the rewriting is right and `augment = TRUE` stays the default.**

- **2026-11-05: Hurdle (two-part) families added.** One `HurdleFamily` covers the four positive parts, and `HurdleAugmentedFamily` makes both forests cheap: Polya-Gamma on the zero forest, and for a count a geometric number K of the zeros the truncation discarded, which leaves `y * eta - (K + 1) exp(eta)` (and the falling exponential for the negative binomial, with the K + 1 Poisson-gamma rates summed into one gamma draw). Tested in `test-hurdle.R` (density against R's d* functions, mean, replicates, augmented against direct), the derivative, quadratic and exponential shortcut tests, and the invariants matrix, where a drawn coding is refused for all but `hurdle_lognormal()`; the gamma row of the exponential test takes 2e-3 because the general path stopped short by up to 9.1e-4 over 12 seeds while agreeing to 1e-15 on 8.
- **2026-11-05: Does the augmented hurdle sampler draw from the right posterior?** `_dev/hurdle-coverage.R` (25 replicates per family, n = 500, default settings, 200 held-out points, each family's own model) found the paired coverage difference on the positive part within 2 SE of zero for every family (Poisson -.003, SE .002; negative binomial -.011, SE .006), coverage .95 to .98 both ways (.89 to .95 for the zero forest, under both samplers alike), and the width ratio within 5% of 1, so `augment = TRUE` stays the default; augmented fits were 2 to 2.7 times faster, timed with seven workers running. The count rewriting's intervals are nonetheless narrower (log width ratio -.026, SE .009, and -.031, SE .016) while the gamma and log-normal controls, which do not use it, sit at 1.00; whether that is the direct sampler's slower mixing overstating width at 800 draws (as recorded for the default run length) or a small error in the rewriting is open, and a long-chain direct fit would tell them apart.
- **2026-11-06: The narrowing is mostly run length, and it is the augmented sampler that is short, not the direct one.** `_dev/hurdle-long-chain.R` (12 data sets per count family, the same process, each fit augmented and direct at 1 x (200 + 800) and at 4 x (500 + 2000); worst R-hat over held-out points 1.03 to 1.05) gave a pooled long-length log width ratio, augmented over direct, of -.007 (SE .006), inside the pre-stated bound of 2 SE and .014, so the rewriting stands and the default is unchanged; within each sampler the augmented width rose with length (+.034, SE .010) and the direct width less (+.014, SE .010), the reverse of the guess in the previous entry. By family, Poisson sits at +.006 (SE .007; +.002 at default length, so the earlier -.026 did not reproduce), while the negative binomial went from -.055 (SE .019) at default length to -.019 (SE .009) long, between the two pre-stated branches, which by the stated rule leads to checking the K and summed-rate conditional with the trees held fixed against an exact sampler for the leaves. Coverage of the count predictor was .94 to .97 in every arm; timings (90 to 1660 s per fit) were taken on a loaded machine and are not comparable.
- **2026-11-06: The negative binomial's augmented steps draw from the exact posterior.** `_dev/hurdle-leaf-check.R` cut the count forest to one tree that never splits (`gamma = 1e-10`, leaf scale fixed at 1, anchor 0), so the count predictor is one number whose posterior, jointly with log theta under its Gamma(0.01, 0.01) prior, is exact on a grid; cells were mu of 0.5 and 3 by theta of 0.5 and 2, theta drawn and held, n = 500, 4 x (1000 + 25000). With theta drawn, the augmented 95% interval width for the predictor over the exact one was .992, 1.019, .988 and .988 (z -0.2, 1.0, -0.2, -2.1 against Monte Carlo error), and .992 to 1.016 for log theta; with theta held, .996 to 1.005; the Poisson controls .994 to 1.004. That is the pre-stated outcome under which the summed-rate draw and the collapsed theta step are right, so the 1.9% left in the forest run belongs to the tree moves' mixing, which this does not exercise, and the default stays. The direct sampler was the one off the grid: at mu = 0.5 with theta drawn its intervals were 5.8% and 3.9% narrow (z -2.9, -2.2), where the zero-truncated likelihood barely separates mu from theta and the predictor's posterior has a long lower tail.

### Multinomial models

**Multinomial models come in symmetric (Murray 2021) and reference codings and as a multinomial probit with a drawn covariance; the multinomial-Poisson transformation was rejected, the probit's sensitivity to the reference category is accepted, and its correlations need no reparameterization.**

- **by 2026-08-24: Symmetric multinomial, per Murray (2021).** One forest per category with no reference and leaf scale divided by `sqrt(2)`; the Laplace proposal needs none of Murray's conjugate machinery. It cut RMSE by .0051 and .0042 (3 and 5 categories) against reference coding averaged over references, landing near the best reference without choosing one (references differ by 37% at five), and mixed and covered better. Cost is the extra forest; `multinomial(reference = )` remains.

- **2026-08-30: The multinomial-Poisson transformation: implemented, measured, rejected.** It is an exact augmentation that leaves identification alone, but the Polya-Gamma route won on ESS/s (113.4 against 36.9 at K = 3, 98.6 against 19.9 at K = 6); lesson: the quadratic form beats the exponential form. Deleted rather than flagged; re-measuring made the Polya-Gamma multinomial the default (9.6x ESS/s hard, 10.1x soft).

- **2026-08-30: Multinomial probit: correlated categories, and the sampler that fits the trees in the normalized space.** `multinomial("probit")` uses Xu et al. (2025) Algorithm [P2] with a trace constraint, a quadratic target at about twice a logit's cost, matched by an independent R reimplementation. The correlation is attenuated by a nonparametric mean (read its sign, not magnitude), and the inverse Wishart `nu` is fixed at `C + 1` because raising it inflates correlations (0.97 at `nu = 50` against a true 0.7). The `loglik` uses fixed normal draws and predictions fresh ones (`replicates = 200`).

- **2026-08-30: The multinomial probit is reference-sensitive, and the symmetric fix is not worth it.** Spread across references was 1.9x seed noise for the probit against 4.1x for the reference-coded logit; held-out accuracy was stable while individual probabilities moved by up to 0.14. A symmetric probit was not built because its covariance would carry `C + 1` unidentified directions.

- **2026-09-03: Log: *BART*'s multinomial support is binary fits, verified.** From *BART* 2.9.10 source, `mbart2()` is K one-vs-rest binary fits normalized afterward and `mbart()` a continuation-ratio factorization whose prior is not exchangeable; the feature table now reads "separate binary fits". Relabelling categories cannot test the asymmetry, since RNG divergence swamps it.

- **2026-09-03: The multinomial probit correlations need no reparameterization.** The vignette's "1.07" was an interval width, now worded explicitly; draws are positive definite by construction (correlations in [.058, .622] over 400 draws), and a per-entry transform could not enforce joint positive definiteness.

### Ordinal models

**Ordinal models gained a cloglog link and a reported chart for the cutpoints; the probit stall without augmentation was contention on the clock rather than the sampler, and the worst-mixing row of an ordinal fit is a threshold the sampler never updates, while inference on many thresholds still wants the joint cutpoint update on To Do.**

- **by 2026-08-24: Ordinal models with many thresholds.** Summing each cutpoint's slice target over its two adjacent categories makes a cutpoint sweep O(n): 0.4% of runtime at n = K = 1600, and fits run to 3200 in under two minutes. One-at-a-time updates mix poorly with many thresholds (7 to 21 effective draws per 100); a tridiagonal joint update is the fix there but conflicts with the ordering constraint, and at K = 20 the slow quantity is the pinned threshold, which it would not reach. Storage grows as `num_draws * K`.

- **by 2026-08-24: Ordinal: a third link, and the chart the cutpoints are reported in.** The cloglog link is the discrete proportional hazards model and brings an exponential-form augmentation (5.2x ESS/s hard, 1.6x soft; on by default). Draws are reported via `Family::report_shift()` in the chart where the predictor has mean zero and every cutpoint is free, with leaf values shifted so the forest replays; this is `polr()`'s chart with the predictors centered. Two categories are exempt to match `binomial()`.

- **by 2026-08-24: Ordinal predictions on the latent and mean scales.** From `WeightIt::predict.ordinal_weightit()`: `type = "mean"` weights probabilities by the numeric labels (`values` overrides; non-numeric labels error), and `type = "stdlv"` divides by `sqrt(var(eta) + var(e))`. WeightIt's `+gamma` cloglog location is a convention, not an error; the two differ by the constant `-mean(lp)/sd`.

- **2026-09-03: The two-category `cloglog` footnote, and what it says.** A two-category `ordinal("cloglog")` gives `P(Y = 2) = exp(-exp(-eta))`, the log-log link rather than the binomial cloglog, so "exactly binary regression" holds only for logit and probit; checked per draw to 1.1e-16.

- **2026-09-17: The error-shapes table was measuring the wrong ordinal model, twice.** The original `_dev/error-shapes-sim.R` was lost and its rebuild is now the design of record; the rebuild's ordinal arm first made every distinct value a category, then mapped bin indices back with `approx()`, doubling RMSE. `bin_response()` in `_dev/positive-sim.R` (bin means) was right (.029 to .043). Lesson: simulations fitting the same family should share the helper that prepares the response.

- **2026-09-17: `ordinal("probit")` with `augment = FALSE` stalls, and the augmentation table hid it.** At 15 replicates one un-augmented probit fit took 3882 s and two soft ones about 600; the table moved from the ratio of means to the median within-replicate ratio (disclosed), giving 15 of 17 rows favoring augmentation and a median worst-quantity ESS ratio of 0.89. `ordsims.sh` rewrote `_dev/positive-results.rds` after it was read; check a result file's mtime against when a number was taken from it.

- **2026-09-17: The `ordinal("probit")` stall was the clock, not the sampler.** Replaying the replicate (`_dev/ordinal-stall.R`) returned identical draws in 29.58 s, and profiling, Newton counters and 60 quiet fits cleared the sampler; all three stalls overlapped other runs, though measured contention (1.4x) does not explain 131x. The clean rerun reproduced all 120 fits' ESS exactly. Lessons: run timings in a `parallel 1` group, and check a suspect timing by replaying and comparing draws; a latent `prob * prob` underflow in `OrdinalFamily::score_info_unit()` was noted, not fixed.

- **2026-09-19: The worst row in an ordinal fit is the threshold the sampler never updates.** `aux.cut1` is `-mean(eta)` in the reported chart, not a cutpoint, and that chart mitigates slow mixing (sampled cutpoints sit at 21 to 49 effective draws); category reversal showed the slow coordinate is the pinned threshold, its R-hat of 1.31 at the `1 + chains/ess` null, reaching the lowest category's probability diluted about eightfold. `diagnose()` was fixed with `zero_if_centered()`, since the centered `eta` average had read as best-mixing; unpinning `cuts(0)` remains unmeasured.

### `custom_family()` and links from R

**`custom_family()` takes a likelihood written in R, with links from R and, since 2026-08-30, nuisance parameters carried as pinned forests; in October a malformed link began to be refused before the fit, `loo()` and `waic()` were fixed for nuisance parameters, and the response is taken as given so `logdens()` can reach per-observation data, while building a family automatically from a `family` object was not started.**

- **by 2026-08-24: Arbitrary links from R.** A caller's link is composed onto the compiled family's native scale (`t = native_link . g^{-1}`), dropping the `t''` term so the information is the expected one (checked for cauchit to 1e-10), with nuisance parameters left to the wrapped family. Not offered for ordinal, multinomial, AFT, location-scale, zero-inflated or `ordbeta()` families, where `custom_family()` is the route.

- **by 2026-08-24: custom_family(): a likelihood from R.** Takes an R log density with `num_predictors` forests, optional analytic `derivatives` (else `score_info_numeric`) and a `start` value; nuisance parameters, non-numeric responses and response-scale prediction were excluded. Matched the compiled Poisson (score to 1e-6, information to 1e-4).

- **2026-08-30: Nuisance parameters for a custom likelihood: the pinned-forest route.** A depth-zero forest carries a nuisance parameter with existing machinery, but only when pinned: unpinned, its tree split and took one value in only 8% of draws, while denying splits gave a posterior sd of 0.0102 against a theoretical 0.0100. The enabling change is `gamma` as a vector, and routing a pinned forest into `fit$aux` was proposed.

- **2026-08-30: Nuisance parameters for a custom likelihood, implemented.** `custom_family(aux_names =, aux_start =)` declares them, pinned through `Family::num_pinned()`, with values in `fit$aux`. Measurement found three bugs: `RFamily::unpack` reading one value per paired block (now per run), Newton overshoot (four-SE trust region), and a Laplace proposal rejected from far out (reversible location damping). A drawn `log_sigma` matched `gaussian()` on mean, sd and ESS from a start seven times off.

- **2026-09-03: `custom_family()`: the changes are sound, with two holes now closed.** `aux_start = NULL` fixes `missing()`-based detection under `do.call()`; partially named and duplicate-named `aux_start` were accepted where `aux_names` refused them, so `aux_names` is now derived first and checked once.

- **2026-09-08: Could a custom family be built from a `family` object automatically?** In principle via `-dev.resids / (2 * phi)` for exponential dispersion families, but the dispersion could not be estimated, `loo()`, `waic()` and `pp_check()` would have to be refused, and `glmmTMB::tweedie()` has only a stub `dev.resids`. Not started.

- **2026-10-06: A malformed link from the caller ran a whole fit and then failed in `colMeans()` (closed To Do item).** `link_functions()` now evaluates `linkinv` and `mu.eta` on a grid and refuses one that does not return a number per value, instead of failing in `fitted_from_eta()` after the fit; test "a link built by the caller is checked before it is used".

- **2026-10-07: `loo()`, `waic()` and `predict(type = "density")` failed on every `custom_family()` fit with nuisance parameters (closed To Do item).** `conditional_density()` in `R/predict.R` now appends each nuisance column as a constant matrix, since the engine counts pinned forests among the predictors; test "a custom fit with nuisance parameters has a pointwise log likelihood", and `?custom_family` notes that terms dropped from `logdens` are missing from what `loo()` reads.

- **2026-10-07: `custom_family()` could not see per-observation data other than `y` (closed To Do item).** The response is now taken as given (`check_custom_response()`) and `custom_callbacks()` subsets it by row with `ss()`, so factors, `cbind()` matrices and `Surv` objects reach `logdens()`; `derivatives` gets the nuisance values as a fourth argument. Tested in `test-flexible.R`, including a bivariate normal recovering its correlation; rows with a missing value in any column are still dropped.

## Survival

### Proportional hazards: `ph()`

**`ph()` is the single proportional hazards family, with a piecewise-constant baseline whose bin count is not a tuning knob (CoxBART adds nothing beyond it); since 2026-10-02 it records draws in the centered chart, so `diagnose()` no longer grades the unidentified level, though its `loglik` row still mixes slowly (on To Do).**

- **by 2026-08-24: Cox proportional hazards.** Deliberately not supported: the partial likelihood couples observations through risk sets, so it does not decompose over a leaf's observations as the leafwise Laplace approximation requires.

- **2026-08-30: Survival, part two: proportional hazards and the discrete-time route.** Built `ph()` from Basak et al. (2024), whose full likelihood has the `poisson()` exponential form, retracting the claim that proportional hazards was out of reach; it won by about 250 log points on a turning-over baseline, lost about 96 to `weibull_aft()` on a Weibull truth, and does not improve risk ordering. The Sparapani et al. (2016) discrete-time probit needs no code (a vignette recipe) and alone copes with crossing hazards; a `survival_expand()` helper and relative survival were not built.

- **2026-08-30: CoxBART, and whether it differs from `ph()`.** CoxBART shares `ph()`'s leaf-level target and is its `B -> infinity` limit under an improper prior; built as `coxph()`, it was within a standard error of `ph()` and correlated 0.982 with `survival::coxph()`'s predictor. `ph()` was recommended for its proper prior, genuine likelihood and prior weights.

- **2026-08-30: Can the bin choice be removed? Yes, and it turns out not to be worth it.** `ph(num_bins = Inf)` was removed (`S(t|x)` worse by +0.017 and +0.010, a common-rate prior over-shrinking narrow bins); a Gamma process `coxph(precision > 0)` matched `ph()` but broke `loo()` (`p_loo` 674 for 500 observations, 56% bad Pareto-k), so the bin count is the regularization that makes `loo()` usable. Two `coxph()` density bugs were fixed and a `fit$loglik` test added.

- **2026-08-30: One proportional hazards family, and the bin count is not a knob.** Estimates were flat over 4 to 250 bins, so `coxph()` was removed and `num_bins` demoted to an advanced argument, with a test pinning the insensitivity. A later re-measurement found a ~20% rise at 250 bins, so the vignette says flat from 4 to 100.

- **2026-10-02: `ph()` records its draws in the centered chart.** `PHFamily::report_shift()` and `aux_values_shifted()` now record draws centered as `ordinal()` does (also under `bcf()`), so the hazard rows went from R-hat 1.6 to 1.00 on `rhc` while identified quantities moved by 7.2e-16 (`_dev/ph-level-shift.R`). The chart-invariant `loglik` row (R-hat 1.18) went to To Do; tests in `test-chains.R` and `test-invariants.R`.

- **2026-10-02: `diagnose()` on a `ph()` fit grades an unidentified coordinate (closed To Do item).** Done 2026-10-02 by centering. The baseline-level rows showed R-hat 1.6 and ESS 3 while the identified `log(lambda_b) + mean(eta)` had ESS 450 to 790.

### The accelerated failure time families and the `Surv` default

**`dpm_aft()` (Henderson et al. 2020) became the default for a `Surv` response once it was made ten times faster, `weibull_aft()` reads as a hazard ratio as well as a time ratio, and since 2026-09-15 a `Surv` response passed to a one-value family is caught at one place instead of failing four different ways.**

- **2026-08-30: The Henderson accelerated failure time model, implemented as `dpm_aft()`.** `DPMAFTFamily` inherits `DPMFamily` with censored-time imputation; it won by 210 log points on a bimodal error, was within 0.1 of `lognormal_aft()` on log-normal errors, and its density reproduced `fit$loglik` to 2.8e-13. Trap: `predict(type = "density")` is on `log T` for AFT families and on `T` for `ph()`, documented in `?predict.bartisan`.

- **2026-08-30: `weibull_aft()`: 15% for free, and why the rest needs a decision.** Making `AFTFamily` `Concrete<AFTFamily>` with a fused `score_info_unit` cut 9.68 to 8.19 seconds. A Frühwirth-Schnatter mixture approximation could go further but would be the package's only inexact augmentation, so it was left as a maintainer decision.

- **2026-08-30: The `Surv` default moved to `dpm_aft()`.** `dpm_aft()` was best or tied-best on four of six truths at 2.8 seconds against `weibull_aft()`'s 9.7, so it became the `Surv` default; the maintainer chose no extra estimand warning. `default_family()`'s weights guard now covers `dpm_aft()`.

- **2026-08-30: `dpm_aft()` was ten times slower than it needed to be, and the sampler was not why.** `response_scale()` bisected the median over a length-n vector though it is shared by every observation; bisecting a scalar and hoisting `mixture_at()` took 13.3 to 2.17 seconds, bit-identical. Lesson recorded: measure where the time is before designing a speedup.

- **2026-09-03: The vignette was wrong about `dpm_aft()` and the time-ratio reading.** A contrast is a log time ratio for every AFT family whatever the error shape; the shape governs only what `e^eta` is alone, which for `weibull_aft()` is the Weibull scale, +.286 from `E[log T]`.

- **2026-09-03: Why the Weibull carries a hazard ratio too.** Only the Gumbel-min error makes an AFT model also proportional hazards, giving a log hazard ratio of minus the contrast over sigma. Fixed in both survival-related vignettes and `?bartisan-families`.

- **2026-09-15: A `Surv` is a numeric matrix, which is how four families and the engine each failed differently.** `check_numeric_response()` flattened a `Surv` to length 2n; it now refuses multi-column responses with a hint naming the survival families or `binomial()`, tested in `test-families.R`. The vignette no longer recommends `custom_family()` for interval censoring.

### Coverage of survival curves on non-normal errors

**`dpm_aft()`'s intervals for S(t | x) under-covered on bimodal and skewed log-time errors; the skewed truth's uncoverable cells were scored away, drawing the baseline helped, and the over-coverage left against `anovaDDP()` comes from the nonparametric mean rather than the mixture.**

- **2026-11-02: A skewed log-time error, and the coverage of S(t | x) intervals.** `_dev/survival-skew-coverage.R` found no skewed-error cost for `dpm_aft_ls()`, but `dpm_aft()`'s 90% intervals covered .55 on bimodal and .67 on skewed errors, the narrowest of any family; causes untested.

- **2026-11-03: Why `dpm_aft()` under-covers S(t | x) on bimodal and skewed errors.** `_dev/dpm-aft-coverage.R` ruled out chains, mean bias and censoring: the mixture's shape posterior is overconfident. Count-share weights carry no spread (Dirichlet-drawn weights raised bimodal coverage .53 to .74 at `k_s = 3`), and `?bartisan-families` may describe `k_s` backwards, to check against George et al. (2019).

- **2026-11-03: The skewed survival truth scored away from its uncoverable cells.** `_dev/skewed-inner-coverage.R`: 25% of cells have true S(t | x) of one, and over the rest coverage is .85 to .88, superseding the sharp-edge reading. Drawing did not help the full model, against prediction; `dpm_aft()` beat `lognormal_aft()` by +.208 (SE .013).

- **2026-11-04: Why `dpm_aft()` over-covers where `anovaDDP()` is at nominal: the mean, not the mixture.** `_dev/linear-mean-coverage.R`: with the true mean as offset `dpm_aft()` covers .928 and .915, and `lognormal_aft()` over-covers equally, so the excess (.958 against .89) comes from the trees. Not tested on a nonlinear mean.

### nftbart, NFT BART and crossing hazards

**nftbart fits `dpm_ls()`'s model rather than the one its paper writes; on censored crossing hazards the default `dpm_aft_ls()` under-covers S(t | x) because it draws the leaf scale, nftbart's better coverage comes from its forests' settings rather than its mixture, and a fixed-scale hard-rule setting closes the gap but is named only as an alternative, since its cost on other truths is unmeasured.**

- **2026-11-02: nftbart fits `dpm_ls()`'s model, not the one its paper writes.** nftbart's code puts DP atoms on the standardized scale, i.e. `dpm_ls()`'s model, unlike the paper's eq. 6. `_dev/nftbart-compare.R`: `dpm_ls()` matched or beat it on score everywhere; the component form was removed 2026-11-02.

- **2026-11-05: `dpm_aft_ls()` against nftbart on censored data.** `_dev/nftbart-survival.R`: `dpm_aft_ls()` was more accurate on all eight truths but covered .838 against .896 on crossing hazards. nftbart failed on the bimodal truth (.583), cause not examined.

- **2026-11-05: Crossing hazards: nftbart's coverage is its forests' settings, not its mixture.** `_dev/nftbart-crossing.R`: the -.070 gap closed through `gate = "hard"` and nftbart's forest settings, not its DPM hyperpriors, at a cost in width and RMSE. The cost on the other seven truths was not run.

- **2026-11-05: Crossing hazards: the step is the fixed leaf scale, and only under hard rules.** `_dev/nftbart-forests.R`: fixing the leaf scale carries the step (+.084), and soft rules undo it. `gate = "hard", k = c(5, 2), update_sigma_mu = FALSE` covered .925 and is named as an alternative only.

### Survival predictions, the vignette and its simulations

**`predict(type = "hazard")` joined the survival predictions on 2026-10-02, the same day the survival vignette was rewritten around what a user does with a fit, stating its simulation results in prose; `_dev/survival-sim.R` was last rerun whole on 2026-11-02 with the location-scale families, and every claim held.**

- **2026-08-30: The survival predict gap, and a marginaleffects bug behind it.** Added `predict(type = "survival", times = ...)` for all survival families. marginaleffects had never worked for any survival family because `data.table::setDT` read the two-column response as length 2n; `get_data.bartisan()` now splits matrix columns, and the survival contrast's intervals covered the truth at all three times.

- **2026-08-30: A survival vignette, and the comparison behind it.** Split out `vignettes/survival.Rmd` with a new comparison from `_dev/survival-sim.R`, `-bins.R` and `-results.R`; the log score is comparable across families with `-sum(status * log(time))`. Trap recorded: `Surv(time, status) ~ .` put `time` in as a predictor, the third time `~ .` cost a measurement.

- **2026-10-01: The survival simulations rerun under `sparsity = FALSE`.** Rerun with the same seeds; every ranking the prose leans on held, though survival RMSE rose 7.7% on average. Three claims were rewritten: `dpm_aft()` "costs little" rather than nothing, it degrades faster than `ph()` under heavy censoring, and the discrete-time curves now cross correctly; the vignette now reads `res$meta$num_draws`.

- **2026-10-02: The survival vignette, rewritten around what a user does with a fit.** Rewritten live on `rhc`, where `ph()` led `dpm_aft()` by 37 elpd (SE 13); old file at `_dev/survival-old.Rmd`. Found, not fixed: no `posterior_predict()` sampler for `ph()` or `dpm_aft()`; `dpm_aft()` mixes slowly on `rhc` (bulk ESS 8).

- **2026-10-02: Swap in the rewritten survival vignette (closed To Do item).** Done 2026-10-02: the old vignette is `_dev/survival-old.Rmd` beside `_dev/survival-results.rds` (allowlisted under Parked), and help pages describing it were reworded (`?bartisan-families` now says 4 to 100 bins).

- **2026-10-02: `predict(type = "hazard")` for the survival families (closed To Do item).** Done 2026-10-02 as `hazard_scale()` in `R/predict.R`, registered with marginaleffects and tested against the derivative of `-log S` to 1e-10 for all five families. It also fixed `dpm_aft()` reading the mixture by row rather than stored draw under `predict(iterations = )`, pinned by a test.

- **2026-11-02: `_dev/survival-sim.R` rerun whole, with the location-scale families.** 376 fits; every vignette claim held and no prose changed. The location-scale families cost little on constant spread and led where the spread varies (.058 and .063 against `dpm_aft()`'s .096).

## Predictors and variable selection

### The sparsity prior: the default, and `bcf()`

**The DART sparsity prior was the default until the atom it puts at zero in an effect's posterior made it the wrong default for causal estimands; since 2026-09-29 `sparsity = FALSE` is the default, `bcf()` follows it, and turning it on for the effect forest alone was measured and rejected.**

- **2026-08-31: Should `sparsity = TRUE` stay the default? Yes, and the guidance had to get sharper.** Superseded 2026-09-29, but the measurements stand (`_dev/sparsity-default.R`, `_dev/sparsity-effect-size.R`). For prediction any sparsity beat none by 13% to 26% RMSE, with `"weak"` within noise of `"strong"`; on a weak treatment contrast (effect 0.05 to 0.2) the prior halved the estimate or worse and 95% intervals covered at 0.40 to 0.60, with `"moderate"` as bad as `"strong"`, so the argument behaves like a switch, not a dial. Default kept `TRUE` then; documentation gained both tables and estimand-indexed advice.

- **2026-09-04: Fixed: a restricted forest silently turned the sparsity prior off.** In `bcf()` with a propensity score, all `sparsity` settings gave bitwise identical fits because `resolve_split_matrix()` built a split matrix from the `vc()` mask and `Hypers::Hypers()` (`src/hypers.cpp:58`) then fixed the weights. Fixed by sending mask and weights down separate channels (`control$split_mask`, `Hypers` `allowed_` index, Dirichlet over allowed groups with concentration `alpha / length(allowed)`); `fixed_s` now means only caller-supplied weights. Ten other configurations compared bitwise against the previous commit were unchanged.

- **2026-09-04: Log: sparsity in bcf(), now that it can be asked for.** Twenty covariates, one moderator `x2`, n = 800, five replicates: the atom at zero was exactly zero in every cell and the ATE within 0.02 of the truth, since the treatment is a coefficient, not a split predictor. Sparsity on the effect forest cut CATE RMSE from 0.417 to 0.254 under strong moderation; under weak moderation `TRUE` on both forests was best (0.174 against 0.217). Not enough to move the default (one DGP, five replicates), but it showed the documented reason for the default was wrong.

- **2026-09-04: Log: bcf() now uses the package's sparsity default.** `bcf()` stopped forcing `sparsity = FALSE` on the outcome model and inherited the package default (then `TRUE`; still inherited after the 2026-09-29 switch to `FALSE`), justified by the entry above (CATE RMSE 0.303 against 0.417 strong, 0.174 against 0.223 weak). `?bcf` and `vignette("causal")` were corrected; the vignette's own single-forest `bartisan()` fits keep `sparsity = FALSE`.

- **2026-09-29: The sparsity prior is off by default.** `bartisan_control(sparsity = FALSE)` became the default on Linero's recommendation, because the prior can drop the predictor whose contrast is the estimand; `sparsity = TRUE` is now framed as for many predictors and selection, not effect estimation. Docs, vignettes and tests moved with it; all eleven vignettes re-knit (7.5 minutes) with numbers and two `vignette("diagnostics")` narratives corrected. Not rerun: `_dev/survival-sim.R`, whose fits used the prior on.

- **2026-09-29: Should `bcf()` turn the sparsity prior on for the effect forest? No.** `_dev/bcf-effect-sparsity.R` (640 fits, 40 replicates): `c(FALSE, TRUE)` helped only with nuisance moderators and real heterogeneity (p = 20, two moderators: d RMSE(CATE) -0.048), did nothing at p = 5, and hurt on a constant effect among many candidates (+0.036) by concentrating splits on noise; ATE bias and coverage were unaffected everywhere. `bcf()` keeps `FALSE` on both forests and `?bcf` states the cost. A Zotero literature read found no paper treating DART as a general default and only Linero and Yang (2018) giving rates, all under sparsity.

### Variable selection: importance, the permutation test and the Gibbs prior

**Variable importance alone cannot support a selection claim; a permutation test was found worth having, the Gibbs prior of Linero and Du (2023) was built, measured and removed before 0.1.0 (kept as `_dev/gibbs-prior.patch`), and a formal selection test is still on To Do.**

- **2026-08-31: Variable importance at n = 1500, and the claim that had to be withdrawn.** At 1500 RHC rows, unlike 5735, the noise controls did not sort below every real predictor (one landed above five), so the vignette now says the gap against known noise is the finding, not the ranking. Every number in the five RHC vignettes was re-read off a fresh render.

- **2026-09-08: A stronger `variable_importance()`, and what the audit around it found.** Added `prop_splits` (within-draw shares; removes `num_trees` dependence but is not invariant, .51 to .64 from 10 to 40 trees), `draws = TRUE`, `plot = TRUE` with a `plot()` method, and a `<bartisan_importance>` print that flags `sparsity = FALSE`. Tests surfaced deprecated `geom_errorbarh()` (now `geom_errorbar(orientation = "y")`) and that `expect_predictor_invariant()` does not hold for `vc()` or offset fits. The `tweedie()` joint-nuisance item came off To Do: separate slice updates cost about 40% of dispersion ESS, no pathology.

- **2026-09-10: The variable-selection gap was half a gap.** `SoftBart::posterior_probs()` is `prop_used` cut at .5, so the posterior inclusion probability and the median probability model were already here and were named in `variable_importance()`'s details. The real gap is a permutation test (`bartMachine::var_selection_by_permute()`, `BART::mc.wbart.gse()`); `_dev/varsel-check.R` was written to measure whether it is worth building, with a null arm to expose selector size.

- **2026-09-10: The permutation test is worth having, and cannot be stacked on the prior.** `_dev/varsel-check.R` (n = 300, p = 25): the .5 cut had FWER 1.00 under the global null, so it is not a test and is no longer documented as one; simultaneous permutation thresholds held size (.10) with useful power. With `sparsity = TRUE` the max threshold's power fell to .12, because the prior shrinks the null and observed shares alike, so any implementation must fit with `sparsity = FALSE`. Recommended a `variable_selection()` on a `<bartisan_fit>`; not built.

- **2026-09-14: Weights were one way to flatten a likelihood, and I had mistaken them for the only one.** `prior_only = TRUE` refused five families; a permutation check showed three refusals were wrong (`gaussian_ls()`, `Gamma_ls()`, `ordinal()`'s `eta`). The DPM fix was `bool informative = w(i) > 0.0` in the mixture update, which also fixed zero user weights contaminating atoms; `multinomial(link = "probit")`'s infinite variance at zero weight (240 `inv_sympd()` warnings per fit) was fixed too. `ordinal()` stays refused: its cutpoint target (`src/family.cpp:549`) has no prior and goes empty at zero weight.

- **2026-09-19: Linero and Du (2023): the Gibbs prior, and the failure it explains.** A prior on the partition of branches by splitting variable, with DART and spike-and-forest as endpoints; it explains `_dev/varsel-check.R`'s median-model false positives as a property of DART (their Table 1: precision 0.98 against DART's 0.72 to 0.79). Implementation is confined to `Hypers::sample_var()` because the prior cancels from the acceptance ratio; cost 1.18 to 1.28 times DART. Gains concentrate at small and moderate P; orthogonality to the Laplace step is unverified outside the conjugate case.

- **2026-09-19: `sparsity = "gibbs"`: built, measured, and removed before 0.1.0.** Never committed; `_dev/gibbs-prior.patch` and this entry are the only record and give the rebuild recipe. Equation (5)'s `!= 0` was inverted by the PDF text layer; read equations off the rendered page. `_dev/gibbs-prior-smoke.R` showed 8.5 against 3.0 false selections on a null; `split_prior`, `share_sparsity` and `share_forests` were refused. Owed: a `num_trees` sweep with held-out error, `zeta` against `p`, and a formal selection test.

- **2026-09-19: The Gibbs prior's recall loss is real, `zeta = 1` is the wrong default here, and the headline was overstated.** `_dev/gibbs-recall.R` (600 fits): at `zeta = 1` recovery fell against DART by 0.47 at p = 20 and 1.82 at p = 50, monotone in `zeta`, and `zeta = 0` restored it (but is anti-sparse). With signal present, the false-positive advantage was modest or absent, so the null-design 8.5 against 3.0 headline only rewarded shrinkage. No `zeta` dominated DART across p; the level was removed before 0.1.0 because its benefit is on a selection criterion the package does not expose.

### Splitting priors: `split_prior`

**`split_prior` sets per-predictor splitting weights, zero included, and is the remedy for the atom at zero in an ATE under the sparsity prior, which is exact arithmetic rather than coincidence.**

- **2026-08-31: Per-predictor splitting priors: `split_prior`.** `bartisan_control(split_prior = )` takes one weight per term, unnamed terms default to 1, and fixes `Hypers::s_` with `update_s` and `update_alpha` forced off in the constructor; it overrides `sparsity`, warning only when `sparsity` was passed explicitly. Realized split shares matched the weights on pure noise. A draft claim that `prop_used` stays at 1 was wrong (0.89 to 1.00 at 20 trees, 1 at 50); unknown names error.

- **2026-08-31: Zero is a legitimate `split_prior` weight.** Zero now holds a predictor out of every tree while keeping it in the model frame, which the varying-coefficient work needs; negative, non-finite and all-zero weights are refused. The all-zero check belongs in `resolve_split_weights()`, not `bartisan_control()`, which cannot see unnamed predictors; a test caught the first placement.

- **2026-08-31: The atom at zero in the ATE, and what actually removes it.** `_dev/ate-atom.R`: under DART a draw without the treatment makes the contrast exactly zero, so the atom mass is `1 - prop_used` (RHC 0.898 and 0.102) and the 2.5% bound is zero once `prop_used < 0.975`. Fixing weights at all removes it (`split_prior = c(rhc = 1)` reproduces `sparsity = FALSE`: 0.0623 [0.0137, 0.1114] against 0.0531 [0, 0.1058]); larger weights barely move it. `vignette("causal")`'s interval was flagged as prior-driven but not changed then.

### Shared forests

**Shared forests (Linero, Sinha and Lipsitz 2020) were read and built twice: `share_sparsity` shipped, while the paper's own model, `share_forests`, stays in the code but not as an argument, after the first version wrongly assumed the likelihood separates across components.**

- **2026-09-08: Shared forests: what was read, what was built, and what was measured.** An oracle experiment (`_dev/shared-forests-sim.R`) showed the prize is variable selection (separate/oracle 1.00x at P = 5 to 2.34x at P = 250). `share_sparsity` (one pooled Dirichlet via `update_shared_s()`) gained up to 1.23x for `gaussian_ls()` log sd and 1.40x for the zero-inflation part at P = 100, costing 2% to 8% when components disagree; not default. `share_forests` was built, but the first version wrongly summed per-component likelihood differences (fixed by sequential visiting; check gap 8.6e-13, 4.2 units with the bug, invisible on `zi_poisson()`), then hidden behind `the$share_forests`. Lesson: `set.seed()` keeps the RNG kind, so anything regenerated from a seed across a `future` boundary is suspect.

### Factor predictors: subset splitting rules

**Subset splitting rules for a factor were implemented and re-benchmarked once five replicates proved too few, a factor was compared as a predictor and as a random intercept, and flexBART's leaf bookkeeping turned out to be here already.**

- **2026-08-31: flexBART (Deshpande; arXiv 2211.04459).** One-hot rules reach only `2^K - K` of the partitions of K levels (1,014 of 115,975 at K = 10), and bartisan's `build_design()` was in that regime; `make_group_probs()` handled selection, not partitioning. Judged worth doing; the open design question was how soft gates compose with subset rules.

- **2026-08-31: Subset splitting rules for a factor, and what they are actually worth.** Implemented as a level bitmask on `Node::mask` with a uniform prior over non-degenerate subsets, categorical rules always hard; `categorical = "onehot"` keeps the old rule. Shuffle-as-ordinal was worse than one-hot. Fixed `Node::Rule` not carrying the mask and `get_limits()` cross-kind collisions; dead `predict_tree()`/`predict_accumulate()` were deleted. The five-replicate benchmark was superseded by the next entry.

- **2026-08-31: Correcting the subset-rule benchmark: five replicates were not enough.** At twelve paired replicates, subset beat one-hot under hard rules at every size (7.3% at ten per level), was indistinguishable under soft rules, and matched or beat flexBART. Soft rules beat hard by 7% to 26%, more than either categorical choice. Default stays `"subset"` and `?bartisan_control` carries the table. Lesson: report paired differences with standard errors, not means of five replicates.

- **2026-08-31: A factor as a predictor or as a random intercept.** `_dev/factor-random-vs-fixed.R`, twenty levels, iid and clustered truths: every route was within noise of every other, except using both at once, which was worse; neither prior won on its matching truth. Prefer a random intercept for exchangeability and new levels; separation is expected only with many thin levels, which was not tested. A first-version bug gave train and test different level effects, signaled by RMSE above the truth's sd.

- **2026-08-31: flexBART's bookkeeping: already here, and what they have that we do not.** `Node::idx` and `split_support()` already cache observation-to-leaf support; `_dev/bookkeeping-cost.R` showed bartisan's ratio to dbarts flat in n (2.78, 3.19, 2.45) while flexBART's falls (3.72 to 1.51). A tree-level map works only for hard rules, so it was recorded as available and not taken; the better soft-rule target is `rebuild_support()` in the bandwidth move.

### `x_transform` and the smoothed CDF

**Neither `x_transform` option dominated and `range` had the worse tail; a smoothed CDF beat both and has been the default since 2026-09-16, with its kernel and bandwidth choices recorded and its documentation moved to `vignette("implementation")`.**

- **2026-09-16: `x_transform`: neither default is dominant, and `range` has the worse tail.** `_dev/xtransform-sim.R`: `range` beat `quantile` in seven of eight DGPs by 4% to 64%, but on `outliers_fine` was 4 to 6 times worse with coverage .217, worsening with n. The cause is the uniform cutpoint prior in `Node::draw_rule()` on a compressed bulk, not the bandwidth. `quantile` kept for its bounded worst case (2.8x, coverage intact); `range` remains right when a derivative is needed.

- **2026-09-16: A smoothed CDF beats both current `x_transform` options on the tail.** `_dev/xtransform-candidates.R`, five arms: `smoothcdf` had the best worst case (1.42 against `quantile`'s 3.01 and `range`'s 6.65) and coverage .960. `robust` and `winsor` fixed `outliers_fine` but broke `mixed` by discarding tail signal. Grid tabulation makes cost negligible. Not implemented in this entry, as a default change is a decision; `range` stays better for slopes.

- **2026-09-16: `smoothcdf` implemented and made the default.** `smooth_cdf_map()` uses a kernel CDF at the `n^(-1/3)` rate, scale `min(sd, IQR/1.349)`, tabulated and normalized on the observed range. Slopes under `"quantile"` diverge (311 at step 1e-5) while `smoothcdf` and `range` settle, `range` closer, so `?bartisan-marginaleffects` recommends `"range"` for slopes. `test-control.R`'s ramp test was tighter than its noise and now averages four seeds.

- **2026-09-16: Epanechnikov and a swept kernel, and why Bergmann and Zaehle's bandwidths are not used.** `.bartisan_smooth_cdf()` computes the map in `O(n + m)` in C++ (exact to `1.2e-15` for Epanechnikov; 4.4 s to 6 ms for ten predictors at n = 10000); Epanechnikov adopted with the plug-in scaled by 2.1. The data-based selectors had worst cases 2.0 to 2.3 times the best arm against the plug-in's 1.07, undersmoothing bimodal predictors back toward a step. Four tests moved for fixture randomness. A text-anchor splice deleted nine functions that `R CMD INSTALL` did not catch; the suite did.

- **2026-09-16: The transform is one `@param` now, and the vignette carries the rest.** The `x_transform` docs went from 48 lines to a 14-line `@param`; failure modes and simulations moved to `vignette("implementation")`. Both DOIs were checked against Crossref metadata and added to `references.bib`.

### Missing predictors

**Missing predictors are handled by MIA with all three of Twala's rules and kept by default, since the old default never worked; whether the missingness rule wants a smaller prior weight is still open (on To Do).**

- **by 2026-08-24: Missing values are kept by default, and the old default never worked.** `na.action` now defaults to `na.pass`. The signature default had never been used: `model.frame()` is called from a call rebuilt from `match.call()`, which omits defaulted arguments, so it fell back on `getOption("na.action")`; the value is now injected explicitly, and `NULL` still means the session option.

- **by 2026-08-24: Missing predictor data: implemented as MIA.** `na_rule` draws among all three Twala variants, only for variables with missing values; the sampler is unchanged and a prior-only test confirms it. A rejected change move did not restore `na_rule`, fixed by `Node::rule()`/`Node::set_rule()`. With 30% missing, MIA beat median fill and tied fill-plus-indicator, withdrawing the claim that the latter is "MIA done badly"; it captured purely informative missingness. MIA targets `E[Y | X observed, pattern]`; `predict()` refuses NAs in columns complete at fit time.

### Random intercepts

**Group-level random intercepts through `(1 | group)` were built on the leaf machinery and beat a factor predictor beyond about a hundred groups; correlated random effects across forests are not implemented (on To Do).**

- **by 2026-08-24: Random intercepts: implemented, and what the leaf machinery gave for free.** `(1 | group)` parsed with `reformulas`; each level is a `Node` updated by `update_scalar()` (factored out of `update_leaf_params()`), with `tau` via `half_cauchy_update_precision_mh()`, and every additive predictor gets its own set. Recovery correlations were 0.93 to 0.98 except zero-inflation (0.75). Slopes are refused. A `y ~ .` bug picked up `(weights)` and `(offset)` from the model frame; unseen levels get zero with a warning.

- **by 2026-08-24: Correlated random effects across forests: not implemented, and why.** A joint update needs mixed second derivatives the `Family` interface lacks; the component-wise route needs a nonzero prior mean threaded through 27 places in `mcmc.cpp`, including the exponential-form recovery, where a mistake is silent. There is also no formula syntax for it and identification is weak above two predictors. Independent intercepts per predictor cover most of the benefit.

- **by 2026-08-24: Random effects: the assessment this was built from.** Against `dbarts::rbart_vi()`, a factor predictor won with few large groups (0.190 against 0.279 at five groups) and the random intercept won beyond about a hundred groups (0.493 against 0.708 at 250 groups of four). Scoring `rbart_vi`'s `yhat.train`, which excludes random effects, inverted the conclusion until caught. The sketched design was built as described.

### Varying coefficients: `vc()`

**Varying coefficients (`vc()`) work through the multi-forest interface; a `vc()` model never refreshed its family's augmentation (fixed 2026-09-05), and the zero-inflated shortfall under `vc()` turned out to be weak identification and too few draws rather than a defect.**

- **2026-08-31: Varying-coefficient BART works through `custom_family()`, and the interface is open.** `_dev/bcf-interfaces.md`, `_dev/bcf-proof.R`, `_dev/bcf-mixing.R`: a two-predictor `custom_family()` recovered `tau(x)` better than the alternative (RMSE 0.150 against 0.233), but `custom_family()` cannot see per-observation covariates. A `LinkFamily`-style decorator would work for every family, keeping closed-form leaf draws but not the exponential form. The rhat of 1.39 was not caused by the varying coefficient; that mixing claim was withdrawn.

- **2026-08-31: VCBART (Deshpande, Bai, Balocchi, Starling, Weiss; Bayesian Analysis 2026).** Judged worth building as one decorator with per-forest formulas and `sparsity`; soft rules and `location_scale()` already cover two of their future-work items. They have a compound-symmetry within-subject correlation that bartisan's random intercept is not equivalent to. Their need for 20,000 to 50,000 iterations for R-hat on trees corroborates the `sigma_mu` decision.

- **2026-09-05: Fixed: a vc() model never refreshed its family's augmentation.** `VaryingCoefficientFamily` did not override `Family::before_forest()` or `reported_loglik()`, so augmented families (logit binomial, negative binomial, both zero-inflated, multinomial, multinomial probit) never redrew their augmentation under `vc()`/`bcf()`: a logit log-odds effect of 1 recovered at 17%, and log likelihood came out positive (+37.5), which was the tell. Probit escaped via `update_aux()`. Fixed with forwarding overrides; RHC logit `bcf()` went from +0.0060 to +0.0641. Regression test in `test-varying.R`; invalidates earlier measurements on `vc()` fits of augmented families.

- **2026-10-01: The zero-inflated shortfall under `vc()`, and what the chase found instead.** `_dev/zi-vc-shortfall.R` falsified the idea that `zero:z` absorbs the count effect: the shortfall was the same with or without a real zero-process effect and shrank with n, so it is weak identification. These fits barely mix (count-effect bulk ESS 6 to 14 of 500 draws), and at n = 1000 `zero:z` reported -0.690 against a truth of zero. Open: nothing says zero-inflated estimands need more draws than the default; ESS against `num_draws` for this family is unmeasured.

## Causal inference and `bcf()`

### `bcf()`: the propensity model

**`bcf()`'s propensity model mixes badly because its leaf scale runs away, and a fixed-scale setting gives 5.5 times the effective draws per second for about 1% of AUC but was not made the default; on the ACIC bench the default score did no better than no score or the oracle, and `gaussian()` stays the family for a continuous treatment, with `propensity_args` to override it.**

- **2026-08-30: The `. + external` formula warning, and the ACIC evidence for the propensity score.** The `'varlist' has changed` warning is base R's `terms.formula`, triggered by `.` plus a term not in `data`; the package deliberately does not work around it. The regularization-induced confounding section now reports the 2016 ACIC evidence (Carnegie 2019: coverage 83.4% to 91.9%).

- **2026-09-04: Log: why a propensity model mixes badly, and the setting that fixes it.** On `rhc` the leaf scale fails to settle (`sigma_mu` 0.420 against 0.21 fixed); `update_sigma_mu = FALSE` lifted loglik ESS from 18 to 249. Documented under `?bartisan_control`'s `update_sigma_mu`, not made a default on one dataset.

- **2026-09-04: Log: what settings a propensity model wants, and a simulation that could not say.** `_dev/propensity-settings.R` scored settings on the effect over six designs: the `oracle` and `none` anchors did not differ, so no setting could matter and none did; no default changed. Harness bug: averaging effect-forest draws with one shift voided coverage, fixed with `coef(draws = TRUE)` (coverage 1.00 to 0.85); the coverage columns of `_dev/propensity-settings.rds` are void.

- **2026-09-05: Log: the ACIC bench, and what it says about the propensity score.** 120 fits: propensity settings again did not matter, but `none` and `oracle` beat the estimated-score default on absolute bias (t = -2.42 and -2.68), most in the poor-overlap, no-alignment cell. Not enough to change the default; open is whether an estimated score helps at all. Coverage 0.79 to 0.88, for simulation-based calibration to settle.

- **2026-09-14: `bcf()` mixing: the propensity score made it worse, not better.** On `rhc`, ATE bulk ESS was 127 with `propensity = TRUE` against 630 without, plausibly two ways to build one surface; seen in one configuration, not concluded.

- **2026-09-25: The family for a continuous treatment's propensity model.** `gaussian()` stays; `dpm()` failed on a point mass at zero (.505 against .902); `propensity_args` overrides it.

### `bcf()`: design and causal estimands

**Design decisions for `bcf()` and its estimands: it returns `<bcf_fit>`, gets no per-arm residual variance, does not force the ATT to match the observed treated mean, runs its fallback refit after `try_fetch()` returns, adopts part of Linero's propensity-score corrections, and leaves difference-in-differences parked.**

- **2026-09-10: A per-arm residual variance in `bcf()`: no, and not as an unbuilt family either.** `_dev/hetero-dpm.R`: on a heteroscedastic truth `gaussian_ls()` with `log_sd = ~ z` and per-arm `dpm()` gained nothing over single `dpm()`; on a homoscedastic one single `dpm()` beat per-arm `dpm()` 10 of 10, so neither a default nor a new family. A 6-of-6 win at six replicates fell to p = .66 at ten.

- **2026-09-10: `bcf()` returns `<bcf_fit>`, and the plan for what dispatches on it.** `bcf()` returns `c("bcf_fit", "bartisan_fit")`, no exact-class comparisons existed; `_dev/plot-methods.md` plans a `causal_effect()` engine by per-draw g-computation on the response scale, since `coef(fit)[, "z"]` is a link-scale contrast.

- **2026-09-22: DiD is parked: no user-facing recommendation until the field settles.** All DiD work moved to `_dev/did/` and removed from the vignettes, Suggests and `R/varying.R`, since none of three designs is a recommendation. Phase 15 (2026-09-30): Wooldridge's extended TWFE needs no package change and removed 90% of linear ETWFE's bias; restart from `_dev/did/DID.md` § "Phase 15", where an ESS at a single chain length means nothing.

- **2026-09-25: The ATT with the observed treated mean (investigated, nothing changed).** A posterior solves no score equations, so forcing the treated mean changes the estimand; imputing only untreated outcomes gives the finite-sample SATT, as *bartCause* does. `_dev/att-observed-mean.R` moved the ATT by at most 65; left as a possible `estimate_effect()` option.

- **2026-09-25: Overlap and the width of the CATE interval.** `_dev/overlap-uncertainty.R` reproduced Li, Ding and Mealli (2023) for per-arm BART, which soft rules do not fix, but bartisan's single-forest designs widen four to nine times in the tails and cover there; `vignette("causal")` was rewritten to say so.

- **2026-09-25: Linero (2024) and Oganisian and Linero (2025): candidates, not decided.** Done: `propensity = TRUE` for a continuous treatment via `gaussian()`, and the score in `vignette("causal")`'s single-forest model. Not started: a score term in the effect forest and a prior check on the confounding bias via `prior_only`.

- **2026-09-28: `try_fetch()` in `bcf()`, and where the fallback fit runs.** Refitting inside the calling handler ran before `run_chains()` restored the RNG stream, changing draws (`_dev/try-fetch-equivalence.R`); the refit now runs after `try_fetch()` returns, identical in all five cases.

## Predictions, comparison and the user interface

### marginaleffects and the prediction scales

**marginaleffects is the interface for counterfactual estimands, reading the draws on the mean and standardized-latent scales; the defects found through it (`avg_comparisons()` returning exactly zero, `~pairwise` for binomial, `by = ~ x3 > 0` grouping by `x3`) are each resolved.**

- **by 2026-08-24: marginaleffects: the draws are the interface.** Estimands push all draws through `get_predict()`'s `attr(, "posterior_draws")`, with frequentist stubs and class registration via `options(marginaleffects_model_classes = )` in `.onLoad()`. Trap documented in `?bartisan-marginaleffects`, not fixed: slopes under the default `x_transform = "quantile"` diverge as the step shrinks (-4.79 at 1e-4 vs -0.28 under `"range"`).

- **2026-08-25: marginaleffects reaches the mean and standardized-latent scales.** `get_predict()` accepts `type = "mean"` and `"stdlv"` (aliases `"probs"`, `"lp"`, `"lv"`), refuses `"class"` and `"density"`, and passes only named `predict()` arguments. Convention: marginaleffects reports the posterior median, `predict()` the mean (2.402 vs 2.393).

- **2026-08-25: The standardized latent variable, extended to binary responses.** `type = "stdlv"` now covers binary probit, logit and cloglog (`cauchit` refused). Trap: the cloglog error enters ordinal and binomial with opposite signs; the ordinal constant would have been off by about 0.80. A two-category ordinal cloglog fit is not binomial cloglog (18% scale difference).

- **2026-08-30: avg_comparisons() returning exactly zero: neither package's bug.** The zero is a posterior atom: `treat` was in no tree in 64% of draws, and marginaleffects reports the median (0; mean 197). Selection under the sparsity prior mixes slowly (chain means 649 to 1039), so the help page says to run several chains; defaults unchanged, and whether 50 trees is too few for DART to mix went to To Do.

- **2026-08-31: `~pairwise` worked once binomial stopped being a categorical family.** The fix was removing binomial from `me_type()`'s categorical list, which had returned mirror-image rows per outcome level; `type = "prob"` still gives both columns.

- **2026-09-13: `by = ~ x3 > 0` grouped by `x3`.** `all.vars()` discarded the expression, silently giving 250 groups; the right-hand side is now evaluated and a wrong-length result errors.

### Posterior predictive draws

**Posterior predictive draws exist for every family with a sampler and open `loo`, `bayesplot`, `performance` and `posterior`; the binned `pp_check()` plots for a binary outcome were fed zeros and ones and now get probabilities.**

- **2026-08-25: Posterior predictive draws, and the packages they unlock.** Added `posterior_predict()`, `posterior_epred()`, `log_lik()`, `simulate()`, `loo()`, `pp_check()`, `as_draws()` and `performance` methods, validated against the C++ log density in all 20 configurations. Binomial replicates are proportions, AFT replicates are event times, `custom_family()` has no sampler, and binomial `newdata` needs `weights`. LOO is documented as strained (high Pareto k); `loo` gets `chain_id`.

- **2026-09-12: `pp_check()` was feeding probabilities' worth of questions to zeros and ones.** `ppc_error_binned` and `ppc_calibration` silently drew degenerate plots from 0/1 replicates; they now get the predictive mean, as *rstanarm* does. The test asserts binned means lie strictly inside (0, 1).

### Model comparison: `loo()` and `kfold()`

**Model comparison through `loo()` and `kfold()`: `kfold()` now refits correctly from inside `lapply()` and on `bcf()` fits, and a correction for comparing a family with a point mass against one without was worked out on 2026-09-13, but the vignette dropped that section on 2026-09-24 and the comparison is back on To Do.**

- **2026-09-12: The logistic regression comparison.** `vignette("comparison")` fits `rstanarm::stan_glm()`, which predicts as well as the forest (`p_loo` about 17 vs 33). `loo_compare()`'s differing-responses warning is false (integer vs double hashes), so the vignette explains it.

- **2026-09-13: Three sections of `vignette("comparison")` that named a problem without solving it.** The held-out score now compares two models (10.82 +/- 4.94; two SEs vs six for `loo()`). Subtracting `log t` reverses AFT vs `ph()` (AFT +3390.9 to `ph()` +131.0); gaussian vs `tweedie()` was judged uncorrectable here. Dropping `rhc` was undetectable (1.79 vs SE 2.54).

- **2026-09-13: The atom correction that the previous commit said did not exist.** Added `loo(x, scale = )` and `waic(x, scale = )` for survival fits, deliberately not automatic. `0a25d3d` was wrong: scoring a zero as `log P(Y = 0) - log(delta)` makes gaussian vs tweedie unit-invariant (1247 and 1248); exact zeros with no `delta` stay ill posed. `p_worse` is at least .5 by construction.

- **2026-09-30: `kfold()` inside `lapply()`, and on `bcf()` fits.** Refits evaluated arguments in the formula's environment; fits now store `fit$call_values` (`call_values()` in `R/interop.R`). `with_bcf_score()` fixed density prediction, so `kfold()` had never worked on `bcf()`. Suite failures from the default change were fixed per test; `R CMD check --as-cran` hung at `checking use of S3 registration`, not from the sandbox.

- **2026-10-06: `kfold()` fails on a `bcf()` fit given its propensity score (closed To Do item).** Each refit now takes its training rows of the score with `ss()` (`bcf_supplied_score()` in `R/bcf.R`); the test "kfold() takes a supplied propensity score fold by fold" failed on the old code.

### `partial_dependence()` speed

**`partial_dependence()` was slow; blocking grid points had nothing to amortize, while skipping the trees that never split on the plotted columns is exact and 2.6 to 12 times faster.**

- **2026-09-16: Blocking grid points in `partial_dependence()`: measured, and there is nothing to amortize.** `predict()` has no per-call fixed cost (-1.3 ms) and 98.7% of time is the C++ forest walk; stacking was 0.86x at n = 5000 and would break per-point RNG streams. The parallel path was not measured.

- **2026-09-16: Tree skipping in `partial_dependence()`: 2.6x to 12x, and exact.** Trees not splitting on plotted columns are evaluated once (`.bartisan_tree_uses()`, `tree_mask`, `eta_to_type()`), with curves identical to 1.11e-16. `pd_tree_mask()` declines for `(1 | g)` factors; unrecognized `...` arguments fall back to `predict()`.

### `bartisan_control()` and the argument surface

**`bartisan_control()` was reorganized into modeling decisions, advanced settings and validation toggles, with per-forest arguments (the formula included, and since 2026-09-18 a per-forest offset at fitting time), `num_draws` in place of `num_save`, and the measurement tables moved out of its help page.**

- **2026-08-30: The control surface, reorganized.** 29 arguments became three groups. `num_trees` takes one value per forest; `soft` is gone and `gate` decides hard vs soft and the shape; `sparsity` (`TRUE` for DART, `FALSE`, or four named strengths) sets the four `update_`/`alpha_shape_` hyperparameters, which win if given directly.

- **2026-08-30: Tables moved out of `?bartisan_control`.** Eight simulation tables (637 lines of `@details`) are archived here as evidence for the defaults, including soft rules peaking at 20 trees vs hard at 50+ and sparsity attenuating weak effects (true 0.05 estimated 0.000, coverage 0.40). The help page keeps the conclusions.

- **2026-08-31: `chains` moved, `num_save` became `num_draws`, `sigma_mu` left the rhat table.** `chains` is a `bartisan_control()` argument still reachable through `...`; `num_save` was renamed at 373 sites; `sigma_mu` left `scalar_draws()` and `fit$rhat`. A seed-fragile test in `test-marginaleffects.R` moved from seed 43 to 3.

- **2026-08-31: Per-forest arguments, including the formula.** `formula`, `num_trees`, `k`, `sparsity`, `split_prior` and other priors take one value or one per forest, per-forest predictors working by zeroed splitting weights; an unnamed forest keeps its default, and multinomial families take single values. Tests found six bugs in code that looked right.

- **2026-09-18: A per-forest offset at fitting time (closed To Do item).** `as.vector(stats::model.offset(mf))` flattened a matrix offset before `build_offset()`; the one-line fix sufficed, and three invariants are now tests.

- **2026-09-18: A test for `offset` beyond one forest (closed To Do item).** Three tests in `test-bartisan.R`; the third pins that under `multinomial()`'s symmetric coding a vector offset is a silent no-op (1.2e-14).

### Formula handling and intercept-only forests

**Formula defects (a one-sided formula reporting a missing object named `.`, `. - days` lost from the stored formula, and `.` in `bcf()`'s `propensity` and `moderators` meaning every column of `data`) were fixed, and `~ 1` now means an intercept-only forest pinned to one drawn scalar.**

- **2026-08-31: The subsample, and the two defects it stopped hiding.** The shipped `rhc` is a random 1500 of 5735 with integer 0/1 `rhc` and `death`, documented as a random sample and with the 0/1 coding's reason.

- **2026-08-31: The `. - days` defect was in the *stored* formula.** `insight::find_formula()` read unresolved `death ~ . - days` as having `days` as its predictor; `bartisan()` now stores the resolved formula, except with random-effect terms.

- **2026-08-31: A one-sided formula reported a missing object named `.`.** The `.` expansion invented a response named `.`; `arg::arg_formula(formula, one_sided = FALSE)` now rejects one-sided formulas first.

- **2026-09-03: `~ 1` already generalized, and now says so.** Confirmed on five multi-formula families; multinomial refuses formula lists. The scalar uses the leaf prior, so it matches built-in families closely, not exactly (`gaussian_ls()` sigma .7012 vs truth .70).

- **2026-09-03: An intercept-only forest, which `~ 1` now means.** A forest naming no predictor gets `gamma = 0`, making it one drawn scalar. Trap: `resolve_vc()`'s early return first skipped the change (714 rules); both paths now carry `pinned`. Against `Gamma("log")`: shapes 4.09 vs 4.13.

- **2026-10-06: `.` in `bcf()`'s `propensity` and `moderators` formulas meant every column of `data` (closed To Do item).** `.` now means the outcome model's covariates (`covariate_terms()` in `R/bcf.R`), tested by "`.` in the propensity and moderator formulas is the outcome's covariates".

### Messages, print and plot methods

**Error messages were losing their second half, print methods wrote half their output to stderr, and headings did not read as headings; each was fixed, `estimate_effect()`'s printed output, contrast label and forest plot were redesigned, and the two functions that had a `plot` argument got `plot()` methods.**

- **2026-09-08: Eight error messages were losing their second half.** `arg::err("main", i = "hint")` stores the hint as a condition field and prints only the first part, with no warning; all eight now use `arg::err(c("main", i = "hint"))`. The family error lists `bartisan_family_names`.

- **2026-09-08: `plot()` methods for the two functions that had a `plot` argument.** Added `plot.bartisan_importance()` and `plot.bartisan_error_density()`; the `plot` argument calls the method, a test compares them, and subsetting first works.

- **2026-09-12: What the effect object prints.** `print()` on `<bartisan_effect>` shows the average potential outcomes; `summary.bcf_fit` was removed, and a `CATE` object prints quartiles.

- **2026-09-12: The contrast label names the quantity.** Labels use symbolic means, `log(O(Y[1]) / O(Y[0]))`, as *lmw* does, with a printed legend.

- **2026-09-12: The forest plot.** The marginal effect moved from a background band to its own interval past the right edge, and marks scale with the unit count.

- **2026-09-12: Every print method was writing half its output to stderr.** `cli::cli_bullets()` writes to stderr, hiding bullets in knitted vignettes; `cli_bullets_cat()` in `R/utils.R` now renders to stdout at eleven call sites, guarded by `test-invariants.R`.

- **2026-09-14: `{.strong}` reads as an accident.** `cli_head()` now underlines ten headings instead of bolding them.

### Interfaces to other packages

**Interfaces to SuperLearner, parsnip and mice were built on 2026-11-02, with the documented limitation that the `dpm()` default can impute outside a variable's support, and coda compatibility is deferred.**

- **2026-09-24: coda compatibility: deferred.** `diagnose()` already reports newer forms of coda's diagnostics. If revisited, build only a thin `as.mcmc.list()` method, with no coda-specific diagnostics.

- **2026-11-02: Interfaces to SuperLearner, tidymodels and mice.** Added `SL.bartisan()`, a `"bartisan"` engine for `parsnip::bart()` with a `"censored regression"` mode, and `mice.impute.bartisan()`. `vignette("workflows")` found BART and logistic regression equal on `rhc` (AUC .736 vs .745). Documented limitation: the `dpm()` default can impute outside a variable's support, remedied by `family` in `blots`; *caret* not done.

## Testing, documentation and release

### The test suite: fragile tests, run time and coverage

**Seed-fragile tests are the recurring failure: they are found by perturbing the sampler and rewritten to check what they mean across seeds (most recently `test-missing.R` on 2026-11-05), and the CRAN-mode run went from 6.6 minutes to 32 seconds by skipping seed-dependent blocks; coverage is 91.7% (on To Do).**

- **2026-09-10: Three things the full suite caught that a filtered run could not.** A filtered run is a different run: unseeded `fit_effect()` exposed under the full suite a real bug where `estimand = "ATC"` averaged over the treated; `focal` now names the group averaged over for both (reproduce with `set.seed(99); invisible(runif(17))`). Two stale message regexes were updated. Final: 2088 passing, 0 failures.

- **2026-09-16: Finding RNG-fragile tests by perturbing the sampler.** Wrapping `.bartisan_fit()` (not `bartisan()`, which breaks `match.call()`) to burn extra uniforms exposed six tests: a coin-flip leading-variable check (now `setequal`, 25/25), correlation floors inside the noise band (longer chains), an unassertable antisymmetry check (rewritten), and two mis-scoped assertions. Four runs then passed 476 tests; the jitter helper is not committed.

- **2026-10-02: Where the test suite's time goes, and what skipping on CRAN would buy.** `_dev/test-timing.R`: skipping seed-dependent blocks was not enough, as 22 seed-free blocks held 4.9 of 6.5 minutes. The 11.9-minute baseline was load-inflated (6.61 idle); `skip_on_cran()` plus `quick_control()` fits brought it to 0.53 (rule in `helper-bartisan.R`). `_dev/test-seeds.R` over offsets 1 to 16 found six fragile blocks, all fixed, two being `diagnose()` bugs.

- **2026-10-02: The tests take 11.9 minutes under CRAN's settings (closed To Do item).** Done: 6.6 minutes on a quiet machine and 32 seconds after the change. The plan skipped the 64 seed-dependent or mixed blocks (saving 5.4 minutes) and shrank the 22 slow seed-free blocks, with shrunk "almost surely" assertions to be run under several seeds.

- **2026-10-06: `test-invariants.R:290` fails under covr, and only there (closed To Do item).** The stepper measured 853 KB under covr against 24 KB otherwise, because covr's tracer forces the source file's `lines` into every srcref. Fixed: `closure_bytes()` passes `serialize()` a `refhook` writing source-file environments as placeholders; both tests pass traced and untraced, and the leaked frame still measures 16 MB.

- **2026-11-03: Two test failures from 797726d's stream change, both fragile tests.** `test-exponential.R`'s Poisson and gamma rows now use tolerance 1e-3 (Fisher scoring stopped up to 1.9e-4 apart); `test-varying.R`'s self-modifying slope now runs 1000 + 500 and is checked one-sided between .5 and 1.35.

- **2026-11-05: `test-missing.R`, "informative missingness is recovered": a fragile check, not a regression.** Bisected to 963a3fe (draws `dpm()`'s weights and baseline), but both builds cleared 1.25 on 18 and 17 of 20 seeds. In `_dev/missing-informative-check.R` over 200 seeds the old form passed 88% with signal and 1% without; four chains against the noise variables' mean passed 200/200 and 0/200, and the test now takes that form.

### Validation and audits

**Audits found bugs that failing tests did not; the `_dev/AUDITING.md` list, including a differential test against dbarts, was worked through on 2026-10-01 and is done apart from the specific items it ends with, and the tail checks of 2026-10-06 fixed the ordinal and beta scores far into the tails.**

- **by 2026-08-24: Validation.** Sampler evidence: detailed balance under a flat likelihood (2.5096 leaves vs 2.5087 expected, `test-recovery.R`), re-derived acceptance ratios, the replay invariant (~1e-9, the most productive test), bit-identical frequency weights. Audits found four bugs (`lchoose` on non-integers, `loglik` missing `log_norm_const`, `sigma_mu_ramp = 1` disabling the leaf-scale update, multinomial overflow), and measurement corrected a list of wrong claims.

- **2026-10-01: Pre-ship checks: decorator audit, identity and chart matrices, recovery across families.** `test-invariants.R` pins each wrapper's inherited hooks. The identity matrix found `BetaFamily::set_aux()` missing `refresh_eta_free()`, and `Beta("probit")`/`"cloglog"` never reaching the sampler; both fixed. The recovery matrix (`_dev/recovery-matrix.R`) found only a defect in its own `ordbeta` generator. Not covered: `mnp()`, most weights and offsets, soft rules.

- **2026-10-01: The audit list, worked through.** Added range checks, offset/weights/`mnp()` arms under both gates, invariant tests, hook counters (`sweep_delivered()`), and calibrated probit and Poisson SBC. The offset arm found the stored intercept in `R/bartisan.R` included the first observation's offset, shifting new-data predictions (0.40 on a Gaussian check); fixed.

- **2026-10-01: The rest of the audit list.** `_dev/differential-dbarts.R` matched dbarts within pre-set bands, with Gaussian intervals 13% narrower because bartisan draws the leaf scale. `_dev/weights-check.R` confirmed weights act as documented. Soft-rule Poisson SBC at 600 replicates was not uniform (chi-square 28.4, p = 0.001, biased low); left open.

- **2026-10-06: Two families' scores fail far into the tails (closed To Do item).** The ordinal logit score underflowed past |z| near 37 and the ordered beta density went `NaN`. Fixed with `logistic_density()` (in `src/utils.h`), `beta_shapes()`, `lgamma_shape()` and `log_diff_logistic()`; the test "the logistic and beta families stay right far into the tails" fails on the old code.

### Vignette content

**The vignettes were written and revised repeatedly for accuracy (the families vignette cut in half, examples moved to the RHC data, `vignette("bartisan")` as the theory document), and the 2026-09-13 read as a reviewer found and fixed the stale claims that remained.**

- **2026-08-30: The families vignette, cut in half.** Cut from 9,414 to 4,917 words, dropping anything a user cannot act on. With `location_scale()` added, the continuous table reads `dpm()` for badly shaped errors, `location_scale()` for varying ones. Three inherited numerical claims were corrected, and the heavy-tail log score dropped because its median ordering flipped between runs.

- **2026-08-30: The causal vignette, revised.** The propensity score is fitted with `bartisan()` (0.13 to 0.78 vs logistic 0.03 to 0.97) and the outcome with `dpm()`. The weighted analysis was removed: prior weights are frequency weights, and g-computation for the ATE should not pass weights to `avg_comparisons()`. Moderation uses `hypothesis = ~pairwise` (the string form is rejected).

- **2026-08-30: The remaining five workflow vignettes, and what writing them found.** Inflated `avg_slopes()` comes from the default `x_transform = "quantile"` step function and is documented in `vignette("effects")`; a guard was built and removed, since marginaleffects swallows warnings from `get_predict()`. The importance demo showed a copied column taking the usage while the joint shift recovers 1.5.

- **2026-08-30: The bartisan vignette as the theory document.** `vignette("bartisan")` became the model reference. An earlier sparsity claim was wrong; with `sparsity = TRUE` noise `prop_used` stays near 0.09 at every tree count. The `eta` R-hat rises with draws (1.29 to 1.37) and dbarts is worse (median 1.22 vs 1.09), so slow fitted-value mixing is a BART property.

- **2026-08-30: The getting-started vignette, and the marginaleffects bug it uncovered.** `get_predict.bartisan()` overwrote marginaleffects' `rowid = -1` scratch markers, so `predictions(newdata = )` returned wrong row counts; fixed by carrying `newdata$rowid` through, moving effects about 10%. Recorded, not acted on: `fit$rhat` is `NULL` at `chains = 1`.

- **2026-08-30: The family documentation moved to a vignette, and two claims in it were wrong.** The help page's 380 lines of Details became `vignette("families")`. The multinomial probit correlation is weakly identified, not attenuated: at n = 900 a true zero can come back at -0.6, at n = 3000 within about 0.2. The test now asserts the deterministic constraints and separation of the extremes.

- **2026-08-31: The examples moved from birthwt to the RHC data.** Worked examples now use a seeded 1500-row subsample of the RHC data (`data-raw/rhc.R`, `data/rhc.rda`). `death` is death during follow-up (mean 0.649), not 60-day death. Build time drove the subsample.

- **2026-09-03: Method in rows, condition in columns.** Every method-by-condition table in the vignettes and `?bartisan_control` was transposed by script; long-format tables, decision tables and the package feature matrix were left alone.

- **2026-09-12: `vignette("bartisan")` disclaimed its own check.** The binary `pp_check()` with `eval = FALSE` was replaced by `type = "loo_calibration"`, and the hand-rolled calibration plot in `vignette("diagnostics")` was removed.

- **2026-09-12: Eight section titles that counted instead of naming.** Count-style and interrogative headings in `vignette("survival")`, `vignette("effects")` and `vignette("importance")` were renamed to descriptive titles (e.g. "Three Traps" to "Details That Are Easy to Get Wrong"); older records keep the old titles.

- **2026-09-13: The rest of the vignettes, read as a reviewer.** `vignette("faq")` documented the deleted `summary.bcf_fit()`, `?bartisan-families` offered Gamma links other than `log` (all fit on `log`), and the Weibull AFT error was mislabeled "standard Gumbel". All fixed.

- **2026-09-18: `vignette("bartisan")` does not cover the bounded gates or either ordinal augmentation (closed To Do item).** Closed with one paragraph on the decision rules pointing at `vignette("implementation")`; the augmentations stay in `vignette("families")`.

### Vignette build time and saved results

**Vignette build time exceeded CRAN's budget, so the vignettes that cost most keep saved results that CRAN replays, refreshed by `just vignettes` (extended to six more after the second submission's check took 26 minutes on Windows); since CRAN's first review (2026-10-18) the save step cannot write into the package directory, and refreshing the saved results other than `diagnostics` is on To Do.**

- **2026-10-02: The heavy vignettes are precomputed, and verified.** Six vignettes became `.Rmd.orig` sources precomputed by `_dev/precompute-vignettes.R` via `rmarkdown::render()`, checked byte-identical to a live render (plain `knitr::knit()` added captions); `R CMD build` took 47 seconds. A quiet-machine re-timing put the cut at about a minute, so only diagnostics, effects and comparison stayed precomputed.

- **2026-10-02: Vignette build time is over CRAN's budget (closed To Do item).** Done 2026-10-02 and changed 2026-10-05: diagnostics, effects and comparison keep saved results in `vignettes/results/`, refreshed by `just vignettes`; the other eight stay live at about 185 seconds. The 2026-09-24 total was 9.2 minutes (`_dev/vignette-timing.R`).

- **2026-10-05: The heavy vignettes keep saved results instead of precomputed sources.** At the maintainer's request `.Rmd.orig` gave way to `vignettes/saved-results.R`, which runs or replays chunks under `BARTISAN_VIGNETTES` and errors on a changed chunk. `just vignettes` (`_dev/refresh-vignettes.R`) renders both modes in place and requires identical bytes and all Suggests.

- **2026-10-05: The saved-results staleness check failed on Windows.** The MD5 was taken from a file written with `writeLines()`, whose CRLF line endings never matched on Windows. `saved.rds` now stores normalized chunk text, compared directly; a mismatch names the line, errors only under explicit `BARTISAN_VIGNETTES=replay`, and warns otherwise. Windows itself was not yet checked.

- **2026-10-07: The second submission's check took 26 minutes on Windows (closed To Do item).** Six more vignettes keep saved results, leaving implementation and faq live; vignettes re-build in 13 s. `test-estimate-effect.R`'s `bcf()` fits use `small_propensity()`, putting the suite at about 28 s; numbers in `_dev/SHIP.md` as of `131ad34`.

- **2026-10-18: CRAN's first review: the vignettes could write into the package directory.** "Run" mode in `vignettes/saved-results.R` wrote into `vignettes/results/`. Writing moved to `_dev/save-vignette-results.R`, injected through the option `bartisan.save_results`; saving was not run on the real vignettes. A `cran-comments.md` was added, reversing the 2026-10-07 decision.

### Documentation audits and style passes

**Documentation passes (voice, audits for stale claims and duplication, headings and enumerations, clefts and argument subjects) were applied across the package; the clefts left in developer-only comments are still open.**

- **2026-09-03: Log: the documentation voice pass.** The `/r-doc-style` skill was applied to roxygen, README and seven vignettes (`_dev/DOCSTYLE.md`). The README went from 8,881 to 795 words, four sections moving into `vignette("implementation")` (original in `_dev/README-moved.md`). Lesson: bare `#` in `@details` or a vignette must be `##`.

- **2026-09-11: A documentation audit, and the duplication list it left behind.** `?bartisan_control` Details went from 637 to 185 lines with no simulation figures. Twenty-five inaccurate claims were fixed, three from that session (e.g. `effect_marginal()` averages before contrasting). Five of 41 roxygen/vignette overlaps were fixed.

- **2026-09-12: The duplication list.** The remaining overlaps were resolved by assigning each an owner (atom at zero, drift calibrations, `vc()` example, default-family table, AFT derivation, `num_bins`, per-forest recycling). The family capability table stays in both places, with the roxygen copy marked canonical.

- **2026-09-15: A documentation audit of ten clean-building vignettes.** 106 prose defects `R CMD check` cannot see: five vignettes failed without a Suggests package, and a misnamed `res$meta$num_draws` silently deleted a sentence. Lessons: read the rendered HTML, and a sentence should not commit to a digit a seeded chunk recomputes. README rewritten (644 words); `_pkgdown.yml` index added.

- **2026-09-16: Nine help pages that should never have existed.** `@keywords internal` still generates an Rd page; `@noRd` suppresses it. All nine `//'` blocks in `src/` now carry `@noRd`, leaving 19 user-facing topics in `man/`.

- **2026-09-16: Bold was the tell, and it is measurable.** `bartisan_control()` had 11 bold spans in about 480 lines, seven times the corpus rate, mostly paragraph openers; two survive. The `@description` went to two sentences and the example gained `set.seed(123)`. The `*Advanced.*` and `*Validation.*` markers are a local convention and stay.

- **2026-09-23: Heading and Enumeration Pass.** Renamed 24 vignette headings and 5 roxygen subheads and rewrote about 50 count-announcing openers, per the r-doc-style skill's rows 25 to 26; committed 2026-09-23, FAQ question headings left as questions.

- **2026-09-25: Cleft and Argument-Subject Pass.** Rewrote about 170 "X is what makes Y" clefts and about 80 argument-subject sentences in user-facing text (rows 29 to 30), uncommitted for review. Open: 103 clefts in `R/*.R` comments, 77 in `src/`, and one in the `R/RcppExports.R` `@noRd` block (fix in `src/`, regenerate with `Rcpp::compileAttributes()`).

### Citations and the bibliography

**Citations were audited by DOI: seven duplicated bibliography entries were merged, a misattributed induced-Dirichlet citation was corrected, and two entries remain uncited.**

- **2026-09-17: The induced-Dirichlet citation was attributed to the wrong document.** `betancourt2019` was built from another paper's reference list; "induced Dirichlet" is named in the 2025 chapter, so the vignette cites both. "Lower bidiagonal" became "bidiagonal". Check citations cheaply with `pandoc --citeproc`.

- **2026-09-17: The bibliography held seven works twice, under two key conventions.** Duplicates were found by DOI; the short-key entries were kept. Three key regexes failed silently; check with `pandoc --citeproc` (a missing key renders as "key?"). `imai2005` needed a new sentence stating the latent covariance prior; `kalbfleisch1978` and `kindo2016` remain uncited.

### Release: pre-ship reviews and CRAN

**0.1.0 (licensed GPL-2 after checking the engine's origin with Linero) went to CRAN after the pre-submission report and the 2026-09-21 pre-ship review, with a categorical `predict()` regression from win-builder and four examples over CRAN's 10 seconds fixed as the checks found them; `_dev/SHIP.md` is the current release assessment.**

- **2026-09-08: Pre-submission report.** `R CMD check` was Status OK. Blocking: the `0.0.0.9000` version. Owed: `--as-cran` hung on a network wait at S3 registration, and no second platform had been run. Not blocking: `tweedie()` references (done 2026-09-18), `vignette("bartisan")` gaps, no posterior predictive draws for `custom_family()` (needs a decision), DART inclusion probability not stored.

- **2026-09-21: Pre-ship review, 2026-09-21.** In `_dev/REVIEW.md`: five bit-exact sampler changes plus dispersion rewrites gave 5-15% on Gaussian, augmented and hard-rule paths; `ARMA_DONT_CHECK_CONFORMANCE` (2%) was not taken. Fixed `coef(newdata = )` `NaN` on varying-coefficient fits, `residuals()` on `dpm_aft()`, `r2_bayes()` with censoring, and one-row survival `predictions()`.

- **2026-09-28: Settle the engine's license with Linero (closed To Do item).** Done: `DESCRIPTION` says `GPL-2`, since the engine adapts `FlexBart`, licensed `GPL 2.0` without "or later".

- **2026-10-05: win-builder found a categorical `predict()` regression.** `1be4e91` let categorical probabilities fall through to `colMeans()` in `response_scale()`, so `predict(type = "response")` and `fitted()` returned one value per category, breaking `bcf()` and `estimate_effect()` with three or more levels; fixed. Seven assertions took the new `match_arg()` message.

- **2026-10-05: Give the package a release version (closed To Do item).** Done: `0.1.0`, in the maintainer's commit.

- **2026-10-05: `cran-comments.md`, and the README's install line (closed To Do item).** `README.Rmd` gives `install.packages("bartisan")` first; `cran-comments.md` was not used by the maintainer's decision of 2026-10-07.

- **2026-10-07: Four examples over CRAN's 10 seconds (closed To Do item).** `partial_dependence`, `estimate_effect`, `print.bcf_fit` and `custom_family` took 10.8 to 13.9 s. The `bcf()` examples now pass small `propensity_args`, and heavy parts went into `\donttest{}`; locally 3.3 s became about 1.1 s. See `_dev/SHIP.md` as of `131ad34`.

### The package name

**The old name `genbart` collided with the archived CRAN package `genBart`, so the package became `bartisan` on 2026-08-30 (`cdc3278`).**

- **by 2026-08-24: Candidate names.** Settled as `bartisan` in `cdc3278`; the shortlist (`anybart`, `glmbart`, `omnibart`, `laplacebart` and others) is kept as a record. Names to avoid included existing packages and `gbart`, the main function of `BART`.

- **2026-08-30: Renamed to bartisan.** The class became `bartisan_fit`; class usages were renamed before the global sweep (1504 occurrences), so Rcpp entry points did not become `.bartisan_fit_fit`. Alignment repair took three attempts (356 lines); the directory move needed the sandbox off. Verified: tests pass, `R CMD check` OK, no `genbart` left.

### Build, dependencies and housekeeping

**Build and dependency traps worth knowing: a debug build that -O2 had been hiding (hence `CXX_STD = CXX17`), `roxygenise(load_code = "installed")`, `Matrix` and `collapse`, which look removable and are not, the rhub libstdc++ deprecation warning through RcppArmadillo's sparse importer (fixed 2026-10-06), and the Dropbox xattr never applied to `.git/`.**

- **by 2026-08-24: Traps and housekeeping.** `.git/` lacks `com.dropbox.ignored` and was left untouched (risk: a later `git gc`). Traps: test regexes must not span a line break in the message source; never reuse a seed for predictors and response (22 tests fixed); defaults on arguments forwarded to `model.frame()` are swallowed.

- **2026-08-30: The debug build did not link, and -O2 had been hiding it.** `static constexpr` `OMEGA_MIN` was ODR-used by `std::max()` under the toolchain's default C++14, which only an `-O0` build exposed. Fixed twice: `src/Makevars` sets `CXX_STD = CXX17` and the call sites use a ternary. Dead `total_trees` removed; the earlier `OMP: Warning #179` and `nm` NOTE were sandbox artifacts.

- **2026-08-30: A note on working alongside concurrent edits.** A maintainer-added `dpm_aft()` appeared mid-session and survived a marker-located removal by luck. Lesson: check for unexpected identifiers before a sweeping edit.

- **2026-09-03: A trap worth recording: `roxygenise(load_code = "installed")`.** It documents the installed build, not the source, and silently dropped `export(Gamma_ls)` from NAMESPACE; `R_LIBS` must point at the build under test.

- **2026-09-15: `Matrix` is one call and stays, because the sparsity is the point.** `Matrix::sparseMatrix()` in `make_group_probs()`: dense storage is 7.7 MB vs 79.6 KB at p = 1000, and `sample_class_col()` relies on sparse column iteration. `Matrix` ships with R, so removal saves nothing.

- **2026-09-15: `collapse` is in Suggests for *marginaleffects*, not for us.** Nothing in the package references it, but marginaleffects uses it on the Bayesian path every fit takes; do not remove it.

- **2026-10-06: The one rhub WARNING: a deprecated libstdc++ function through RcppArmadillo's sparse importer (closed To Do item).** `std::get_temporary_buffer` came via `Exporter<SpMat>`'s `std::stable_sort`. The engine now builds from the `dgCMatrix`'s parts (`csc_parts()`, `sparse_from_parts()` in `src/utils.cpp`); four fits are bit-identical. The other 17 failures were not ours.

### Literature assessed

**Papers read and assessed for adoption (VCBART, flexBART, SBT, and the 2026 sweep: horseshoe leaves, Prado et al., monotone BART, Castillo and Rockova); most were declined or kept as citations, and what they added is on To Do.**

- **2026-08-31: Three papers assessed: VCBART, flexBART, SBT.** All three are the VCBART citation neighborhood.

- **2026-08-31: SBT (Luo and Pratola; arXiv 2306.00361).** Sharding the data into separate BARTs is a scalability device that changes the model into a mixture over shardings. Not to be built; recorded as the route if sample size ever becomes the binding constraint.

- **2026-09-19: Five papers from the 2026 literature sweep, read and assessed.** A sweep returned 328 papers, the citation-graph snowball half being the useful one; five were read in full, XBART derivatives in `_dev/XBART.md`.

- **2026-09-19: Jacobs, van Wieringen and van der Pas (2026): horseshoe leaves.** Horseshoe leaves need the non-conjugate engine bartisan already has, and conditional on scales the leaf prior is Gaussian, so only per-leaf `lambda` and per-forest `tau` draws are new. Whether it fixes the `sigma_mu` separation warning is unmeasured.

- **2026-09-19: Prado et al. (2025): a real fix, for the stump case.** Double-grow and double-prune moves bar main effects in a shared covariate. It maps onto `vc(z, ~ 1)` (the warning at `R/varying.R:371`) but not general `vc()`, and its `beta ~ MVN(b, Sigma_beta)` with `Sigma_beta ~ IW(I, p_1)` supplies the prior `fc()` lacked.

- **2026-09-19: Chipman, George, McCulloch and Shively (2022): monotonicity is per-tree and topological.** The constraint factors over trees as bounds on leaf values, so here it would be a truncation interval on the existing leaf draw. Constraining the additive predictor would give one `monotone =` argument for every family, which no other package can offer.

- **2026-09-19: Castillo and Rockova (2021): a citation, with a caveat.** Adaptive bands with coverage under self-similarity, with functionals the easy case; `_dev/coverage-calibration.R` cannot see such a condition. It covers single-tree CART, not sum-of-trees, so cite for the shape of the result only.

- **2026-09-19: What changed on the To Do list.** The variable-selection test gained two routes (a Gibbs prior making exclusion a posterior event; Ye et al. 2025's no-refit summary), neither built. The leaf-scale prior item was reframed by Jacobs and is still a default change not made unilaterally. `fc()` now has a prior design but stays out of scope for 1.0.

- **2026-09-19: Read and declined.** ASBART is XBART grow-from-root plus a bandwidth pass, with the 10x being 40 against 4000 iterations (`_dev/XBART.md` §2.5). The "Anytime Monte Carlo" concern does not apply: chains stop on iteration count, not wall-clock time.
