> **Status: parked, 2026-09-22.** Nothing in this directory is referenced by
> the package. The DiD section was removed from `vignettes/causal.Rmd` and
> `did` dropped from `Suggests` because *bartisan* is not ready to recommend how
> to do difference-in-differences with BART. The measurements below stand; the
> recommendation is what was withdrawn. See `_dev/TASKS.md` for what was removed
> and `PLAN-categories-B-and-D.md` for the outside design that phase 14 tests.

# Difference-in-differences with bartisan: investigation and vignette plan

The question was whether *bartisan* can fit the DiD-BCF model of Souto and
Louzada Neto (2025), "Forests for Differences: Robust Causal Inference Beyond
Parametric DiD" (arXiv:2505.09706, `soutoForestsDifferencesRobust2025`),
substituting this package's BART and BCF for the XBART and XBCF the paper uses.

**Answer: yes.** The mapping is exact, the recovery is good, every estimand is
reachable from `estimate_effect()`, and the one package change it wanted --
moderating the effect by a variable the control function does not see -- has
since been made. What has to be said out loud is that the specification is not
likelihood-identified and needs enough treated mass to behave.

## The model, and why it is already a bartisan model

Their Eq. 3, after the parallel-trends reparameterization, is

$$Y_{it} = \mu(D_i, t, \mathbf{X}_{it}) + \tau(\mathbf{X}_{it}, k_{it}) \cdot D_{it} + \epsilon_{it}$$

with $D_i$ the ever-treated indicator, $k_{it} = t - G_i$ the event time, and
$D_{it} = \mathbb{I}(G_i \neq \infty)\,\mathbb{I}(k_{it} \geq 0)$ the
treated-now indicator. bartisan's varying-coefficient model is

$$g(\mu_i) = f_0(Z_i) + \sum_j (X_{ij} - c_j) f_j(Z_i)$$

and these are the same model when $c = 0$. That is not a coincidence to be
arranged: `vc()` centers a `0`/`1` covariate at zero by default (`"auto"`
resolves to `"zero"`), so the treatment term is **exactly** zero on every row
where `Dit == 0`, which is the structural zero the paper's reparameterization
is for. Verified: the fitted basis column equals `Dit` to the bit, and rows
with `Dit == 0` contribute nothing to the effect forest, since their score and
information are multiplied by the basis.

So the model is one formula:

```r
bartisan(y ~ ever + t + x1 + ... + vc(Dit, ~ x1 + ... + t),
         data = panel, family = gaussian(),
         control = bartisan_control(num_trees = c(50L, 25L)))
```

`num_trees = c(50L, 25L)` reproduces BCF's asymmetry, fewer trees on the effect
than on the prognostic function. Confirmed working alongside this:
`(1 | unit)` random intercepts, which reach the control function and no
coefficient -- something the paper's model does not have, since its errors are
i.i.d.

### Use `bartisan()` and `vc()`, not `bcf()`

`bcf()` forces `center = "estimate"` on a binary treatment, which draws the
coding rather than fixing it. The contribution at `Dit == 0` is then
$b_0 \tilde\tau(\mathbf{X}, k)$ rather than zero, so the control group's
baseline picks up a function of event time and the structural zero is gone.
For an ordinary BCF that is a harmless reparameterization; for DiD it is the
one property the model is built around. `bcf()` also defaults to estimating a
propensity score, which the paper explicitly declines and which is degenerate
here anyway -- `Dit` is a deterministic function of cohort and period.

## Estimands: `estimate_effect()` reaches all of them

Every estimand in the paper is the varying coefficient averaged over a chosen
set of rows, and `estimate_effect()` is the function that does that averaging.
No user-side programming and no reaching for `coef()`:

| Estimand | Call, with `treated <- subset(d, Dit == 1)` |
| --- | --- |
| ATT | `estimate_effect(fit, treat = "Dit", newdata = treated)` |
| GATT | the same, plus `by = ~ cohort` |
| Event study | the same, plus `by = ~ k` |
| CATT | the same, plus `estimand = "CATE"` |

Checked against a hand-rolled average of `coef(draws = TRUE)`: the ATT and the
CATT agree to 1e-8, so the two routes are the same numbers. `by = ~ k` works
although event time is not a predictor of the model, because `by` is evaluated
in `newdata` rather than against the fit -- a property of the `model.frame()`
rewrite of `effect_by()` made the same day. There is a test for all four in
`test-estimate-effect.R`.

`coef(fit, draws = TRUE)[[1]]` remains the draws-by-rows matrix of
$\tau(\mathbf{X}_{it}, k_{it})$ for anyone wanting an estimand `by` cannot
express.

**The trap that has to be documented,** and it applies to `newdata` as much as
to `coef()`. `coef()` returns the coefficient function on *every* row,
including rows where the treatment is off, and handing the whole panel to
`estimate_effect()` averages over those rows too. Those rows
carry zero weight in the effect forest, so the value there is the forest's
extrapolation and not an estimate of anything. Averaging without first
restricting to `Dit == 1` produces a plausible-looking number that is not the
ATT, and a "pre-treatment effect" that is an artifact. Measured on cell A
below: the unrestricted event-study profile reads about 2.8 at every
pre-treatment event time where the structural truth is zero.

## What the simulation found

`_dev/did-bcf-dgp.R` and `_dev/did-bcf-study.R`. Two of the paper's own DGPs,
20 Monte Carlo replicates each, N = 200 units by T = 8 periods, against static
TWFE with unit and period fixed effects via *fixest*.

Cell A is their DGP 2 / Setting 1: staggered adoption (cohorts at t = 4, 5, 6
plus never-treated), homogeneous $\tau = 3$, fully linear -- TWFE is correctly
specified. Cell B is their DGP 5 / Setting 3: staggered adoption with selection
into cohort on static covariates, effects heterogeneous in $X_3$ and $X_4$,
every covariate entering non-linearly, quadratic time trend.

| | true ATT | bias | RMSE | coverage | CI width |
| --- | --- | --- | --- | --- | --- |
| **A** bartisan | 3.000 | +0.001 | 0.074 | 0.95 | 0.397 |
| **A** TWFE | 3.000 | +0.010 | 0.063 | 1.00 | 0.373 |
| **B** bartisan | 5.461 | +0.007 | 0.142 | 1.00 | 0.609 |
| **B** TWFE | 5.461 | -0.086 | 0.344 | 0.90 | 1.202 |

In cell B the per-observation CATT is recovered at correlation 0.931 with RMSE
0.337 against a true CATT spread of 0.881, so most of the heterogeneity is
found rather than smoothed away.

Three readings. Where the simple model is right, the forest costs essentially
nothing: RMSE 0.074 against TWFE's 0.063, and coverage is nominal at 0.95
rather than TWFE's conservative 1.00. Where it is wrong, the forest is better
by a factor of 2.4 on RMSE with correct coverage where TWFE undercovers. And
the effect forest did not invent dynamics it was free to invent: the truth is
flat in event time and across cohorts, and the estimated event-study and GATT
profiles are flat to within 0.01 in cell A and 0.08 in cell B.

## The two things to change or to write down

### 1. `vc()` modifiers may now name a variable the control function cannot see

**Done 2026-09-21.** `vc(Dit, ~ x + k)` used to error unless `k` was also a
predictor of the model formula, and for DiD that was a real bind: the effect
should vary with event time, and event time has no business in the control
function, which already carries the cohort and the period that determine it.

The restriction turned out to be a typo guard rather than a structural
requirement, and the comment on it said so. A modifier that was not a design
column would have been dropped silently by `intersect(asked, groups)`, fitting
a forest with fewer modifiers than were asked for, and the error existed to
stop that being silent.

A modifier the fixed part does not carry is now given a design column of its
own and reaches the forest that asked for it and no other. The masks already
did the hard part: they are built per parameter from that parameter's own fixed
formula, which does not name these, so every control function excludes them
automatically. `vc_modifiers()` gained a `slope_groups` argument so a
coefficient can reach further than its own control function, and `bartisan()`
appends the terms to the model frame, the design and the stored formula. Only
an explicit `modifiers` formula reaches such a variable; a bare `vc(z)` still
means the predictors of the model, which these are not.

Measured on cell A with `vc(Dit, ~ x1 + x2 + k_fin)`: the effect forest splits
on event time in 15% of draws, the control function splits on it in exactly
0%, and the ATT comes back at 2.921 against a truth of 3.

So the DiD form is now literally the paper's Eq. 3, with event time where the
paper puts it:

```r
bartisan(y ~ ever + t + x1 + ... + vc(Dit, ~ x1 + ... + k),
         data = panel, family = gaussian(),
         control = bartisan_control(num_trees = c(50L, 25L)))
```

**Implications elsewhere.** This is a capability, and like any capability it
can be asserted wrongly. Naming a variable as a modifier only is a substantive
claim: that it moves the effect and not the baseline. Where that is false the
control function cannot absorb it, and the coefficient's forest is the only
place the baseline variation can go, so it lands in the estimand. DiD is the
case where the claim is safe because it is forced -- the control function must
not see the cohort on identification grounds. Elsewhere it should be a
deliberate choice rather than a convenience, and the roxygen says so.

Two smaller consequences, both handled. A typo now reaches `model.frame()`
rather than the modifier check, which would report `object 'x' not found`
without saying which term asked for it; it is caught earlier instead, with the
term named. And the stored formula gains the term, so
`insight::find_formula()` and *marginaleffects* see every variable the model
uses, which is the same reason the frame's formula is what gets stored at all.

### 2. The specification is not identified, and thin panels show it

$\mu$ sees the ever-treated indicator and calendar time, whose interaction is
the treatment, so the likelihood does not separate $\mu$ from $\tau$ and the
split rests on the priors. This is a property of the paper's model, not of
bartisan, and the paper does not discuss it.

It is benign when there is enough treated mass and severe when there is not.
Measured:

- Simulated panel, 200 units and 600 treated observations: coverage 0.95 and
  0.90-1.00, bias under 0.01. Fine.
- `causaldata::castle` (Cheng and Hoekstra castle-doctrine laws, 50 states by
  11 years, 21 ever-treated, 74 treated observations): R-hat for the ATT is
  1.15 to 1.23 and bulk ESS is 13 to 21 **whatever the moderator set** --
  including `vc(Dit, ~ 1)`, where the effect is a single drawn scalar. Turning
  the sparsity prior off moves it to 1.14. Across draws the ATT and the control
  function's fit at treated rows correlate **-0.89**, which is the ridge
  itself.

So the guidance is: run `diagnose()` on the ATT draws, and read a low ESS as
the model telling you the panel cannot separate the two forests rather than as
a sampler that needs longer.

### 3. Under a null, the effect forest absorbs what the control forest underfits

In the hard null cell (DGP 5 / Setting 3 with the effect set to zero) the ATT
comes back biased **+0.18** where the truth is exactly zero. The 95% interval
is wide enough to cover it, so the test stays at its nominal size, but the
mechanism is worth naming because it is the one that would bite on a longer or
more curved panel.

The obvious suspect was the quadratic time trend: $0.2t^2$ spans 12.8 over the
panel, and the treated post-treatment rows all sit at the high end of $t$,
so a control function that underfits the curvature there would leave a residual
for the effect forest to absorb. **Tested, and that is not it.** Twelve paired
replicates, same seeds across specifications:

| Control function | bias | sd | mean interval width |
| --- | --- | --- | --- |
| numeric `t`, 50/25 trees | +0.184 | 0.133 | 0.873 |
| `factor(t)`, 50/25 trees | +0.195 | 0.284 | 0.950 |
| numeric `t`, **150**/25 trees | **+0.064** | 0.109 | 0.665 |

Period dummies do not help at all, so the trend's functional form is not the
problem. Tripling the control forest cuts the bias to about a third and
narrows the interval by a quarter, so what is short is **capacity**: eight
covariates entering through `exp`, squares, absolute values and square roots,
plus a quadratic trend, is more prognostic surface than 50 trees carry, and the
effect forest is where the rest goes.

Capacity can be bought two ways, so the follow-up asked which: more trees, or
deeper ones. The branching probability is $\gamma (1 + d)^{-\beta}$, so a
smaller `beta` on the control forest deepens its trees without adding any. A
period random intercept was tried alongside, as the cheaper way to write a time
effect on panel data. Twelve paired replicates, `_dev/did-bcf-capacity.R`:

| Control function | bias | sd | mean interval width |
| --- | --- | --- | --- |
| 50/25 trees (reference) | +0.184 | 0.133 | 0.873 |
| 150/25 trees | +0.064 | 0.109 | 0.665 |
| 50/25, control `beta = 1` | **+0.271** | 0.220 | 0.969 |
| 50/25, control `beta = 0.5` | +0.175 | 0.221 | 0.842 |
| 50/25, `(1 \| period)` | **+0.051** | 0.105 | **0.458** |

**Depth is not a substitute for trees, and is worse than doing nothing.**
Halving `beta` left the bias where it was and nearly doubled its spread;
`beta = 1` made the bias half again as large. That is the standard BART
intuition holding: the ensemble's capacity comes from the number of weak
learners, and strengthening each one costs stability without buying
flexibility. Whatever a poorly mixing DiD fit needs, it is not deeper trees.

**A period random intercept is the best of the five**, beating even the tripled
forest on bias and beating everything on interval width, which it nearly
halves against the reference. The reason is that the period effect is eight
discrete levels; a random intercept estimates them directly where the forest
has to spend splits approximating a step function in `t`, and those splits are
capacity not spent on the covariates.

**It does not transfer to `mpdta`,** which is the caution to carry. Replacing
numeric `year` with `(1 | year)` there is worse on every reading: interval
width 0.138 against 0.120, ESS 571 against 1146 at 150/25 trees, and 274 at
50/25. Five periods is too few to identify a variance component and there is
little time structure for a forest to struggle with in the first place. So the
random intercept is worth reaching for on a long panel with real time
structure, and is not a default.

Twelve replicates is thin, so the ordering is solid and the magnitudes are not.


### 4. The identification argument does not pick the specification. Measurement does.

Worth recording because the reasoning above is seductive and, taken one step
further, wrong. The control function sees the group and the period, and a
function of those two *is* the treatment, so the model is not likelihood
identified. It is tempting to conclude that the control function should
therefore be told as little as possible -- in particular that it should get a
bare ever-treated indicator rather than the adoption cohort, since only the
cohort lets it form $\mathbb{I}(t \geq G_i)$ exactly.

**Tested on `did::mpdta`, and it is the worst specification of the five.** Four
chains, 500 warmup and 800 kept draws, 150/25 trees, R-hat and bulk ESS on the
ATT draws:

| | ATT | 95% interval | R-hat | ESS |
| --- | --- | --- | --- | --- |
| A. $\mu$: cohort + year; $\tau$: lpop + since | -0.013 | [-0.074, +0.046] | **1.009** | **1146** |
| B. $\mu$: cohort + year; $\tau$: lpop + year + cohort | -0.017 | [-0.083, +0.042] | 1.035 | 160 |
| C. $\mu$: **ever** + year; $\tau$: lpop + since | **+0.093** | **[+0.030, +0.253]** | 1.093 | 28 |
| D. $\mu$: cohort + year; $\tau$: lpop | -0.012 | [-0.071, +0.046] | 1.011 | 1513 |
| E. as A, with 50/25 trees | -0.015 | [-0.084, +0.052] | 1.028 | 248 |

Callaway-Sant'Anna gives -0.042 (se 0.012) and TWFE -0.037 (se 0.013) on the
same data, so C is not merely noisy: it is the wrong sign with an interval
excluding zero. The cohorts differ in their baselines, and a control function
not allowed to say so leaves that difference for the effect forest, which is
the estimand. What the control function is short of on a panel is capacity,
not restraint -- the same conclusion item 3 reached from the null cell, by a
different route.

Two further readings. A against E is the tree count on real data: 248 effective
draws against 1146 for the same number of posterior draws, which is the
clearest evidence for `c(150L, 25L)` anywhere in this file. A against B is the
moderator set: giving $\tau$ the year and cohort as well as event time costs a
factor of seven in ESS for no change in the estimate, because year and cohort
are what the control function is already using, and the ridge widens with every
direction the two forests share.

**So the guidance is A**, and it is what `vignette("causal")` uses: give the
control function the cohort and the period, and give the effect forest the
event time as a modifier the control function does not carry. That
specification is only expressible because of item 1.

### And a pre-trend check is not available from this model

The obvious idea is to fit their Eq. 2 instead, `vc(ever, ~ ... + t + cohort)`,
which lets $\tau$ be non-zero before treatment so that $\hat\tau$ at $k < 0$
is a parallel-trends check. **Tested, and it does not work.** On cell A, where
the pre-treatment truth is exactly zero, that form returns $\hat\tau$ between
0.67 and 0.99 at every pre-treatment event time, and an ATT of 2.754 with a
95% interval of [-0.835, 4.433] against a truth of 3. That is a false positive
for a parallel-trends violation and an interval five times too wide. It is also
exactly what the paper predicts: their whole argument for the reparameterization
is that a forest asked to learn the zero constraint on $k < 0$ does it badly.

The vignette must therefore say that parallel trends is checked outside this
model -- a conventional event-study regression with *fixest*, or the `did`
package -- and not with $\hat\tau$.

## Plan for the vignette section

A new section of `vignette("causal")`, after the BCF material, since it builds
on `vc()` and on the same estimand vocabulary. Eight parts:

1. **What it is for.** Two or three sentences: the estimand is the ATT on panel
   data under parallel trends, and the model is the paper's DiD-BCF, which is
   the DiD design written as a varying coefficient so the effect may vary with
   covariates and with time since treatment. Cite the paper.

2. **The model and the formula.** Their Eq. 3, then the one-line bartisan
   formula, then the sentence that makes it legitimate: a `0`/`1` basis is
   centered at zero, so the treatment term is exactly zero pre-treatment and
   pre-treatment rows cost the effect forest nothing. State that `bcf()` is the
   wrong entry point here and why, in one sentence.

3. **Setting the data up.** The three columns a DiD fit needs beyond the
   covariates -- ever-treated, period, treated-now -- plus cohort for a
   staggered design, in the four lines that build them from an adoption-year
   column.

4. **Reading the estimands off.** The table above, with a worked ATT, a GATT by
   cohort and an event-study plot. The restriction to `Dit == 1` gets its own
   sentence with the reason, not a parenthesis.

5. **Effect heterogeneity**, which is the reason to reach for this rather than
   TWFE: `coef()` per row, and the existing partial-dependence and
   *marginaleffects* machinery applied to it.

6. **What the model needs.** The identification paragraph and the `diagnose()`
   recommendation, with the castle numbers as the concrete case of it failing.

7. **What it cannot do.** Pre-trends, with the measured reason.

8. **Pointers.** The paper; `vignette("effects")`; *fixest*, *did* and *did2s*
   for the conventional estimators and for the pre-trend check.

### Dataset: `did::mpdta`, chosen and verified

`rhc` is not a panel and cannot be used. The section uses **`did::mpdta`**:
500 counties by 5 years, staggered minimum-wage adoption in 2004, 2006 and
2007, 309 never-treated counties and 291 treated observations. It is the
canonical staggered DiD teaching set, it is the application
@soutoForestsDifferencesRobust2025 use themselves, and its one covariate is
log county population, which is exactly the heterogeneity their applied
section reports. `did` is now in `Suggests` and every chunk is guarded by
`has_did`.

Measured rather than assumed: under specification A above it gives R-hat 1.009
and 1146 effective draws on the ATT, which is a fit worth printing, and it
supports a covariate-heterogeneity plot, a GATT by cohort and an event study
by time since treatment.

`causaldata::castle` was tried and rejected: 50 states by 11 years with only 74
treated observations gives R-hat 1.15 to 1.23 and ESS 13 to 21 whatever the
moderator set, including a constant effect. It survives here as the example of
a panel too thin for the method, if a cautionary case is ever wanted.

### Cost

The section is one fit plus a `did::att_gt()` comparison. `mpdta` is 2500 rows
and 175 trees; at four chains it took about 226 s before the varying-
coefficient family was made `Concrete` and had its per-observation allocation
removed, and that change is worth a factor of about 2.4 on this path, so
expect roughly 90 s now. See the follow-up section of `_dev/REVIEW.md`.

---

# Phase 13: the single-fit stacked Callaway–Sant'Anna design

Everything above this line was written for the Souto & Louzada Neto DiD-BCF
specification. That paper is unrefereed and the user has declined to build on
it, so the plan above is superseded. What follows replaces it. The estimand is
now Callaway & Sant'Anna's group-time ATT computed by outcome regression, which
is `did`'s `est_method = "reg"` and is peer-reviewed; `did` is the reference
implementation to check against rather than a competitor to beat.

## The design

For each cohort `g` and each period `t`, take the cohort-`g` units plus the
never-treated, set the outcome to the long difference `Y_it − Y_i,g−1`, and take
the covariates at `g−1`. Label the rows with `cell = g:t`, `cohort = g`,
`since = t − g`. Stack every cell into one data frame and fit one model:

```r
dY ~ cell + X + vc(treat, ~ X + cohort + since_f)
```

`treat` is **"is a cohort-`g` unit"**, not "is treated now". That distinction is
the whole design: it makes every `t` comparable, including `t < g − 1`, where
the same contrast is a placebo.

Three properties follow, and each one fixes a specific failure recorded above.

1. **Unit effects cancel by construction.** The long difference removes them
   exactly, so there is no random effect and no shrinkage bias. The unit-RE
   form in phase 12 returned +3.45 against a truth of 3 with 0/6 coverage
   because the shrunk intercepts under-fit the group baseline and the deficit
   loaded onto τ. Differencing has no such term to get wrong.
2. **The identification ridge is gone.** Within a cell the control function
   sees the cell and `X` but never the unit's cohort, so it cannot reconstruct
   `D_it`. The ridge that gave cor(ATT, control fn) = −0.89 and ESS 13–28 on
   thin panels does not arise.
3. **One fit, every estimand.** `estimate_effect(by = ~ cohort)` is
   `aggte(type = "group")`, `by = ~ since_f` is `type = "dynamic"`, and no
   argument at all is `type = "simple"`. This is the single-model multi-period
   property the user asked for, and the reason to prefer this over per-cell 2×2
   fits.

## Validation against `did` on `mpdta`

One fit, `dY ~ cell + lpop + vc(treat, ~ lpop + cohort + since_f)`, 50/25 trees,
4 chains, 1000 burn / 2000 draws. Post-period cells only. `did::att_gt(xformla =
~ lpop, control_group = "nevertreated", est_method = "reg")`.

| estimand | bartisan | `did` (se) |
|---|---|---|
| overall ATT | −0.0469 | −0.0420 (0.0118) |
| cohort 2004 | −0.0756 [−0.118, −0.037] | −0.0851 (0.0261) |
| cohort 2006 | −0.0432 [−0.078, +0.002] | −0.0204 (0.0182) |
| cohort 2007 | −0.0318 [−0.063, +0.002] | −0.0288 (0.0160) |
| e = 0 | −0.0361 | −0.0211 (0.0124) |
| e = 1 | −0.0560 | −0.0534 (0.0170) |
| e = 2 | −0.0900 | −0.1411 (0.0358) |
| e = 3 | −0.0802 | −0.1075 (0.0334) |

Every cell agrees inside the interval. The two largest gaps, `e = 2` and
`e = 3`, are the two thinnest cells (n = 20 each) and both move toward the
pooled effect: that is the regularization doing what it is for, and bartisan's
intervals there are narrower than `did`'s.

The 2×2 case was checked separately against `DRDID` and matched on both the
point estimate (−0.0277 vs −0.0288) and, once the Bayesian bootstrap over
treated units was used, the width (0.065 vs 0.063).

## The pre-trend check, which is the part that failed before

Emit cells at `t < g − 1` as well — same construction, `since` negative — and
the placebo falls out of the same fit. Compare against `did` with
`base_period = "universal"`, which is the convention the stacking uses; `did`'s
default varying base period is not comparable and will look like disagreement.

| e | bartisan | `did` universal (se) |
|---|---|---|
| −4 | +0.0051 [−0.017, +0.028] | +0.0069 (0.0236) |
| −3 | +0.0039 [−0.015, +0.028] | +0.0276 (0.0177) |
| −2 | +0.0037 [−0.016, +0.028] | +0.0235 (0.0148) |

All three cover 0, and under a DGP where parallel trends holds exactly the
pre-period estimates covered 0 in 4/4 cells in each of 3 replicates. Contrast
the Eq. 2 placebo in phase 8, which returned 0.67–0.99 where the truth was 0:
that check could not fail. This one can.

**Caveat, and it is not small.** Pooling pre-cells into the same fit shrinks the
post estimates toward the near-zero pre-cells, because `since_f` partial-pools
across event time: `e = 2` moves from −0.090 to −0.062 and `e = 0` from −0.036
to −0.006. On the simulation the post-period mean moves from 2.85 (post-only,
6/6 covering) to 3.10–3.47 (with pre-cells). Both cover, but they are not the
same estimator. The honest handling is to fit post-only for the effects and
add pre-cells for the placebo, and to say so rather than presenting one fit
that quietly does both. Do not present the pooled fit as if it were free.

## Two hypotheses tested and refuted, recorded so they are not retried

- **"Unit heterogeneity breaks the level specification."** Added unit effects of
  sd 2 to their DGP 2 / Setting 1. The level form was unaffected in bias and
  coverage (5/5, mean 3.03); only its intervals widened. Persistent
  heterogeneity per se is not what goes wrong.
- **"The `mpdta` gap is regularization at a small standardized effect."** BART
  scales the leaf prior to `sd(y)`, and `effect/sd(y)` is 0.028 on `mpdta` in
  levels against 0.256 in differences, a factor of nine. Scaling their τ down to
  match gave ratios to truth of 0.00–1.89 (level) and −2.29–1.90 (stacked)
  across five replicates: at that effect size their design has no power and the
  estimates are noise, so the comparison cannot discriminate anything.

What does distinguish the two settings is **serial dependence**. Their DGPs have
iid noise of sd 1 and no unit effects, so differencing doubles the error
variance and buys nothing; `mpdta` has level sd 1.51 against innovation sd 0.16,
so differencing removes almost all of it. That is why the two settings rank the
two forms oppositely, and it is a property of real panels rather than of either
estimator. It is also the reason the earlier "the level form is attenuated 60%
on `mpdta`" claim was stated too strongly: it came from one point estimate and
was never checked against its own interval.

## What the vignette section becomes

Replace the DiD-BCF section currently on disk in `vignettes/causal.Rmd`. Keep
the `did` Suggests dependency and the `has_did` guard. The shape:

1. The estimand — ATT(g, t), with a sentence on why the 2×2 is the wrong unit to
   write for a package that can fit all of them at once.
2. Setting up the variables from a bare panel: how to build `cohort` from a
   per-period treatment indicator, which is the step the user asked to have
   spelled out.
3. The stacking, as a short function the reader can copy.
4. The single fit, then `estimate_effect()` three times for simple, group and
   dynamic.
5. The `did` comparison table as the validation.
6. The pre-trend check, with the pooling caveat stated plainly.

`stack_cells()` is a dozen lines and is currently only in the scratchpad. It
should either live in the vignette as printed code or become an exported helper;
that is a design decision, not a detail — an exported `did_stack()` would make
this a supported workflow rather than a recipe.

## Phase 13b: the ridge is in the model after all, and 10 trees answers it

The claim above that "the identification ridge is gone" was wrong, and the way
it was caught is worth recording because the defaults hide it.

`diagnose()` at the package defaults reported the contrast at R-hat 1.013 and
ESS 674, which looks fine. Raising `num_draws` from 800 to 2500 dropped it to
75. ESS cannot fall as a stationary chain lengthens, so this was either drift or
a defect in our diagnostics. Checked against an independent implementation:

| draws | quantity | `diagnose()` | `posterior` |
|---|---|---|---|
| 800 | Y[0] | 1.084 / 32 | 1.084 / 32 |
| 1600 | Y[0] | 1.090 / 28 | 1.090 / 28 |
| 2500 | Y[0] | 1.102 / 25 | 1.102 / 25 |
| 2500 | Y[1] | 1.005 / 5095 | 1.005 / 5073 |

`diagnose()` agrees with `posterior::rhat()` and `posterior::ess_bulk()` to
rounding at every length, so **the diagnostics code is correct**; this is not a
package bug. What the numbers say is that Y[1] behaves normally while Y[0]
drifts, so `mu + tau` is pinned on the treated rows and the split is not.

The mechanism: `cell` is in the fixed part and `cohort + since` are modifiers,
but `cell` **is** the interaction of those two, so on the treated rows the
nuisance can represent anything the effect can. The never-treated rows pin `mu`
only where their covariates lie, leaving slack in the covariate direction. The
design removes the ridge that phase 8 found (the nuisance cannot form `D_it`
from the cell label alone) and introduces a different one.

Four variants at 2500 draws, which is long enough for the drift to show:

| specification | ATT | contrast R-hat / ESS | Y[0] ESS |
|---|---|---|---|
| baseline, `vc(treat, ~ lpop + cohort + since)`, 50/25 trees | −0.0437 | 1.038 / 75 | 25 |
| no `lpop` among modifiers | −0.0437 | 1.025 / 150 | 97 |
| no `cohort`/`since` among modifiers | −0.0455 | 1.014 / 383 | 44 |
| **200/10 trees, all modifiers kept** | **−0.0443** | **1.011 / 858** | **122** |
| `cell` itself among modifiers | −0.0440 | 1.028 / 106 | 41 |

Shrinking the effect forest to 10 trees is the answer, and it keeps every
modifier, so all three aggregations and the heterogeneity section survive. The
ATT does not move (−0.0443 against the baseline's −0.0437 and `did`'s −0.0420)
and the interval width is unchanged, so the mixing is bought with the prior
rather than with bias. Putting `cell` among the modifiers, which is the maximally
collinear case, is near-worst at 106 and confirms the mechanism.

This is BCF's own logic arriving by a different route: the nuisance needs
capacity because it carries a level per cell, and the effect should carry only
deviations from it. The package default of 25 effect trees is too many for this
design.

**Lesson for any future measurement here: read a DiD fit's diagnostics at more
than one chain length.** Every configuration in phase 13 looked acceptable at
800 draws and none of them were.

## Phase 13c: do not pick a configuration from a short-chain ESS

Trying to make the vignette's two fits cheaper produced the trap in its purest
form. Contrast R-hat and bulk ESS for `dY ~ cell + lpop + vc(treat, ~ lpop +
cohort + since)` on `mpdta`, 4 chains:

| trees (nuisance/effect) | burn + draws | wall | contrast R-hat / ESS |
|---|---|---|---|
| 100 / 10 | 400 + 600 | 100 s | 1.005 / **1262** |
| 100 / 10 | 500 + 1000 | 157 s | 1.037 / **85** |
| 50 / 10 | 400 + 600 | 52 s | 1.027 / 147 |
| 200 / 10 | 500 + 1000 | 173 s | 1.016 / 649 |
| 200 / 10 | 500 + 2500 | ~400 s | 1.011 / 858 |

`100 / 10` looks like the best configuration on the table and is the worst one
on it: 1262 effective draws at 600 and 85 at 1000 is not a converged chain, it
is a chain that has not yet drifted. `200 / 10` goes 649 at 1000 to 858 at 2500,
increasing, which is what a stationary chain does.

So **200 nuisance trees is load-bearing and not a comfort setting**. With 100 the
nuisance has too little capacity for a level per cell and the level wanders; the
effect forest at 10 trees cannot fix that, because it is not the part that is
short. The vignette keeps `c(200L, 10L)` and pays 173 s for the post-only fit
and 268 s for the fit with pre-treatment cells.

The general rule this group of scripts has now demonstrated three times: **an ESS
for this model means nothing at a single chain length.** Every wrong conclusion
in phases 13 and 13b came from reading one. Run two lengths and check the ESS
went up.

---

# Phase 14: the Category-B design, coded in bartisan

`did_B_D_and_bart_design.md` (an outside plan, in `~/Downloads`) proposes a
Category-B staggered-DiD model and argues it is the right frame for a Bayesian
package. Its §3.2 model is

    Y_it = alpha_i + lambda_t + f_0(X_i, t)
           + D_it [ tau_{G_i,t} + f_tau(Xdot_{i,G_i}, G_i, K_it) ] + eps_it

with the load-bearing restriction that `f_0` never sees `ever`, `cohort` or
`D_it`, and a free scalar per treated cell so the cell mean cannot leak into the
baseline.

## The formula mapping

Every block has a bartisan equivalent, and the `vc(x, ~ 1)` trick carries the
parametric ones. `vc(z, ~ 1)` gives the coefficient's forest no predictor to
split on, so every tree is a stump and the coefficient is one drawn scalar,
i.e. `z` enters linearly (under the leaf prior, so weakly shrunk rather than
flat).

| plan block | bartisan |
|---|---|
| `alpha_i` | `(1 \| id)` |
| `lambda_t` | `vc(tf, ~ 1)`, `tf` a factor: one scalar per period |
| `f_0(X, t)` | the formula's fixed part, which *is* the control forest |
| `tau_{g,t}` | one `vc(d_g_t, ~ 1)` per treated cell, `d_g_t` a 0/1 column |
| `f_tau` | `vc(Dit, ~ dot_x... + cohort + since)` |

Two mechanics make this work rather than merely typecheck. `vc()`'s `"auto"`
centering is `if (all(x %in% c(0, 1))) 0 else mean(x)`, so every 0/1 column
centers at **zero**: `d_g_t` contributes exactly nothing off its own cell and
`Dit` contributes nothing on an untreated row, which is the structural zero the
design needs. And a `vc()` on a factor creates **one forest per level**, so
`vc(tf, ~ 1)` on eight periods is eight forests; the plan's model on `dgp_A`
has 22 (intercept + 8 periods + 12 cells + `Dit`), and `num_trees` must then be
length 1, length 22, or named. Pinned stump forests want 1 tree, not 50.

The build is `_dev/catB-build.R` (`prep_B()`, `form_B()`, `att_B()`).

## Extracting the estimands: `estimate_effect()` works only on the folded form

Three routes, checked against each other on one fit (`dgp_A`, seed 11):

| route | plan's model | folded model |
|---|---|---|
| `coef()`, summing the `d_g_t` and `Dit` coefficients per treated row | +3.678 [3.522, 3.840] | n/a |
| `predict(draws = TRUE)` g-computation, zeroing every treatment column | +3.678 [3.522, 3.840] | +3.5816 [3.4179, 3.7728] |
| `estimate_effect(treat = "Dit", estimand = "ATT")` | **−0.227, silently wrong** | +3.5816 [3.4179, 3.7728] |

The two hand-rolled routes agree exactly, so either validates the other.
`estimate_effect()` does **not** work on the plan's model, and the way it fails
is the dangerous way: it takes one treatment variable, turns `Dit` off, and
returns the average of `f_tau` alone, because the `d_g_t` columns are separate
variables it never touches. No error, just the wrong number by the whole of
`tau_gt`.

It *does* work, to 1.1e-14 against g-computation, on the **folded** model, where
the whole treatment block hangs off `Dit`:

```r
y ~ X + t + vc(tf, ~ 1) + vc(Dit, ~ dot_x... + cohort + sf) + (1 | id)
estimate_effect(fit, treat = "Dit", estimand = "ATT")            # simple
estimate_effect(fit, treat = "Dit", estimand = "ATT", by = ~ cohort)  # group
estimate_effect(fit, treat = "Dit", estimand = "ATT", by = ~ sf)      # dynamic
```

Folding is legitimate, and cheaper than it looks. `cohort` x `since` **is** the
cell, so `f_tau`'s forest can represent `tau_gt` itself; what is lost is that
the cell effects are now shrunk toward each other under `f_tau`'s prior instead
of free. The plan's §1.5 saturation requirement exists to stop a cell mean
leaking into a baseline that is flexible in the (i, t) direction, and the plan's
own §3.3 argues this baseline cannot form `D_it` at all, so the free `tau_gt` is
a bias-variance choice here, not an identification requirement. Folding
therefore costs shrinkage across cells and buys the whole `estimate_effect()`
interface.

`marginaleffects` was not needed: `predict(draws = TRUE)` g-computation is the
general route, it differences on the response scale, and it is what §3.6
requires for a nonlinear link.

## What the Category-B model recovers, and when

Three DGP cells, 3 replicates each, hard gates, 2 chains, 250 + 500. Truth is
tau = 3. "stacked" is the vignette's long-difference model.

| DGP cell | B bias | B width | B cov | stacked bias | stacked width | stacked cov |
|---|---|---|---|---|---|---|
| `dgp_A` as published (0.75 * ever, no unit persistence) | **+0.510** | 0.391 | **0/3** | −0.292 | 1.029 | 3/3 |
| level shift subtracted, same noise draws | −0.011 | 0.325 | 3/3 | −0.292 | 1.029 | 3/3 |
| level shift kept, unit effect of sd 2 added | −0.007 | 0.420 | 3/3 | −0.286 | 1.070 | 3/3 |

The failure needs **two** conditions at once, which is narrower than it first
appears. `f_0` is forbidden to see `ever` or `cohort`, so a group level shift can
only be carried by `alpha_i`; `dgp_A` redraws its covariates every period and so
has almost no between-unit variance (0.105, implying sigma_alpha ~ 0.32), the
random intercept is shrunk to sigma_alpha ~ 0.13, and the unabsorbed remainder
lands on tau, which is the only other term switched on for treated-post rows.
Remove either condition and the bias goes: subtract the shift and B is exact;
keep the shift but add unit persistence and sigma_alpha becomes large, the RE
stops shrinking, and B is exact again.

Giving `f_0` a pre-treatment mean outcome does **not** rescue it here (3.503
against 3.501, 3.368 against 3.398, still 0/3). On `dgp_A` the pre-period
outcome is dominated by `X beta` noise of sd ~2.1 against a 0.75 signal, so it
carries almost no unit-level information. That would differ on a persistent
panel, but on a persistent panel the RE already works.

The stacked model is **invariant across all three cells** (−0.29 bias
everywhere, to three decimals within a replicate), which is what differencing
buys: a unit-level term of any size cancels exactly. It pays for that with a
standing shrinkage bias of −0.29 and intervals 2.5 times wider, and covers
because it is wide rather than because it is accurate.

## On `mpdta`, B works and still loses

Folded B, `num_trees = c(200, 1 x 5, 10)`, 500 + 1000, 4 chains.

| | B (folded) | stacked | `did` "reg" |
|---|---|---|---|
| overall ATT | −0.0362 | −0.0442 | −0.0420 (se 0.0115) |
| interval width | 0.0625 | 0.0503 | 0.0451 |
| contrast R-hat / ESS | 1.044 / 64 | 1.016 / 649 | — |
| cell estimates outside `did`'s CI | 2 of 7 | 0 of 8 | — |

B agrees with `did` on the overall ATT and is inside its interval, so it is not
wrong. It is simply worse on every other axis, and `sigma_alpha` comes back
**0.077**: `f_0(lpop, year)` absorbs the county variation through `lpop`, so the
random intercept has little left to carry and supplies little of the
within-unit correlation it was introduced for.

The unifying variable is the one phase 13c already identified. Differencing
collapses `mpdta`'s response scale ninefold (level sd 1.51 against innovation sd
0.16), so the same absolute effect is nine times larger relative to the noise
the model must fit, and the differenced model wins. `dgp_A` has iid noise and no
persistence, so differencing doubles the error variance and buys nothing, and B
wins wherever it is unbiased. **B's advantage is real only when the effect is
large relative to level-scale variation**, which real panels rarely satisfy.

Verdict: keep the stacked model in the vignette. B is worth knowing as an
alternative, it must be written in the folded form for `estimate_effect()` to
reach it, and it is not more precise on real panels. Scripts: `_dev/catB-*.R`.


## Citations verified for this work, held here rather than in the package bib

Both were added to `vignettes/references.bib` for the removed section and
reverted with it, so they are recorded here. Each was checked twice, against
the reference implementation's own `citation()` and against Crossref, per
`~/.config/agents/PAPERS.md`.

```bibtex
@article{callaway2021,
  title = {Difference-in-Differences with multiple time periods},
  author = {Callaway, Brantly and Sant'Anna, Pedro H. C.},
  year = {2021},
  journal = {Journal of Econometrics},
  volume = {225},
  number = {2},
  pages = {200--230},
  doi = {10.1016/j.jeconom.2020.12.001}
}

@article{santanna2020,
  title = {Doubly robust difference-in-differences estimators},
  author = {Sant'Anna, Pedro H. C. and Zhao, Jun},
  year = {2020},
  journal = {Journal of Econometrics},
  volume = {219},
  number = {1},
  pages = {101--122},
  doi = {10.1016/j.jeconom.2020.06.003}
}
```
