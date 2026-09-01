# Varying coefficients: what was built, and what the measurements decided

The model is

$$g(\mu_i) \;=\; f_0(Z_i) \;+\; \sum_j (X_{ij} - c_j)\, f_j(Z_i)$$

with a forest for the control function $f_0$ and one for each varying
coefficient. This is VCBART (Deshpande, Bai, Balocchi, Starling and Weiss,
2026); Hahn, Murray and Carvalho (2020) is the case of one binary $X_j$, and
Woody, Carvalho, Hahn and Murray (2020) the case of one continuous one. The
vocabulary here follows Woody et al., because theirs is better: $f_0$ is the
**control function** over **control variables**, $f_j$ is an **exposure
moderating function** over **moderators**.

The interface design is `/Users/NoahGreifer/.claude/plans/cuddly-dancing-axolotl.md`
and the survey of what every other implementation does is in the same place.
This memo records what got built and, more usefully, the four things the
measurements settled differently from how the plan guessed.

## The shape of it

```r
bartisan(y ~ x1 + x2 + vc(z), data = d)          # f0(x1,x2) + z*f1(x1,x2)
bartisan(y ~ x1 + x2 + vc(z, ~ x1), data = d)    # f1 splits on x1 only
bartisan(y ~ x1 + vc(z1) + vc(z2), data = d)     # several
bcf(y ~ x1 + x2, treatment = ~ z, data = d)      # the causal wrapper
```

`vc()` is a formula term, pulled out the way `(1 | g)` is: the covariate reaches
the model frame, so it gets the same missing-value handling as any predictor, and
does not reach the design, so it is not also something the control function
splits on. The parser walks the expression rather than going through `terms()`,
which is what lets it run before `.` is expanded and leave `.`, `. - x1`, `a:b`
and `I(x^2)` spelled as written.

**The engine change is one decorator.** `VaryingCoefficientFamily` in
`src/family.cpp` presents $1 + J$ predictors outward, hands $\mu$ inward, and
scales `d1` by $X_{ij}$ and `d2` by $X_{ij}^2$. The wrapped family never learns
any of this is happening, so **every family gets varying coefficients with no
per-family code**. Because the map is linear in each $\eta_j$, a
`TARGET_QUADRATIC` family stays quadratic and keeps its closed-form leaf draw;
the exponential form does not survive, since `exp_rate()` is one scalar per
predictor and a varying coefficient makes the rate vary by observation, so
`poisson()` and `Gamma()` fall back to `TARGET_GENERAL`.

**Estimands came for free**, as predicted. `get_coef.bartisan_fit()` returns an
empty vector by design, so everything goes through `predict()` on modified
`newdata`; combining the forests inside `predict()` is therefore all it took.
Checked exactly: the contrast through `predict()` and the mean of the
coefficient's draws agree to the last bit (0.9389819584 both).
`avg_comparisons()` differs by 1.1e-4, which is its own grid handling rather
than ours.

## Four things the measurements decided against the plan

### 1. The reference point is not one choice, it is two

The plan said centre at the mean by default, reasoning that it decorrelates the
two forests and helps mixing. On the binary case that is simply wrong.

Recovery of $\tau(x)$ on the `_dev/bcf-proof.R` simulation, and on a continuous
covariate shifted away from the origin:

| covariate | `"mean"` | `"zero"` | `"mid"` |
|---|---|---|---|
| binary 0/1 | cor 0.975, rmse 0.208 | **cor 0.987, rmse 0.144** | -- |
| continuous, near 0 | 0.990 | 0.990 | 0.990 |
| continuous, around 50 | **0.990** | cor 0.422, rmse 1.071 | 0.990 |

Neither wins everywhere, so the default is `"auto"`: zero for a `0`/`1`
covariate, the mean for any other numeric one. For a binary covariate zero is a
value it actually takes and the control function there is the surface among the
untreated -- Hahn et al.'s prognostic score, a quantity with its own meaning. For
a covariate around 50 the control function at zero is an extrapolation, and
recovery collapses.

**And the prior scale, which the plan spent a paragraph on, does nothing.**
`k = 2`, `4` and `8` on the coefficient forest all give correlation 0.975 and
RMSE 0.21 on the same design. The whole of the gap between the interface and the
original proof was the centering.

### 2. A continuous covariate modifying its own coefficient is the point, not an edge case

The plan treated `vc(z, ~ z + x1)` as a curiosity that is "identified and
something a caller might mean". It is more than that: it is how the effect stops
being linear in the covariate. With the truth $y = 2x_1 + z^2$:

| | RMSE of the fitted surface |
|---|---|
| $f_1$ may not split on $z$ | 1.187 |
| $f_1$ may split on $z$ | **0.099** |

Twelve times better, and the fitted coefficient traces the slope of $z^2$, which
is $z$: -1.15, 0.15 and 1.09 at $z$ of -1, 0 and 1. So `vc(z)` on a continuous
covariate is a *linear* dose response with a varying slope, and
`vc(z, ~ z + ...)` is the way out of that assumption. `?vc` says so with these
numbers.

The categorical case is genuinely different and stays a silent removal: a level's
indicator is nonzero only on the rows where that level holds, and the variable is
constant on exactly those rows, so the split separates rows that contribute from
rows that contribute nothing. Wasted rather than unidentified.

### 3. `center = "estimate"` is built, and what it buys is not what the plan said

Hahn et al. section 5.3 replaces the fixed coding with $b_{z_i}$ drawn from
$N(0, 1/2)$ per level. `draw_coding()` in `src/family.cpp` does it in one
conjugate normal draw per level: only the rows at level $k$ carry $b_k$, so the
design is block diagonal and the levels are independent given everything else --
no matrix to build or invert. `refresh_coding_column()` then rewrites the basis,
and `coef()` returns the identified contrasts against the first level.

The plan's arithmetic was right and its emphasis was wrong. It said the coding
buys **invariance** for a numeric binary covariate and **parsimony** for a
factor. Both are true, but the sizes are the opposite of what the plan implied.

**Invariance is real and small.** Fitting the same data with the treatment coded
`0`/`1` and again `1`/`0` and adding the two effects gives zero if the coding
does not matter. Over 12 replicates:

| | `"zero"` | `"mean"` | `"estimate"` |
|---|---|---|---|
| $n = 800$, effect 1.0 | 0.0125 | 0.0074 | **0.0045** |
| $n = 150$, effect 0.2, noise 2 | 0.0830 | 0.0524 | **0.0389** |

The ordering is what Hahn predicts and the drawn coding is best at both signal
strengths, but at $n = 800$ every number is small next to an effect of 1. The
prior difference the whole argument is about is a prior difference, and 800
observations swamp it. It shows up where it should: on weak data the gap widens
sixfold in absolute terms. **The first design measured this at $n=800$ with a
clear effect and found nothing, which was a wrong measurement rather than a
finding.**

And it costs nothing to have: effect recovery is 0.2014 against 0.1991 for a
fixed zero, a paired difference of 0.0023 with a standard error of 0.0140. So
`bcf()` uses it for a binary treatment, which is what Hahn's own software does.

**Parsimony is large, and it cuts both ways.** For a three-level covariate:

| truth | symmetric, a forest per level | `"estimate"`, one shared forest |
|---|---|---|
| every level the same shape | 0.266 | **0.191** (-0.075, se 0.010) |
| each level its own shape | **0.242** | 0.696 (+0.454, se 0.019) |

Twenty-eight percent better when the rank-one restriction holds, and nearly three
times worse when it does not. That is a big enough spread in both directions that
it has to be the caller's choice rather than a default, so the symmetric coding
stays the default and `?vc` states the trade with these numbers.

### 4. `(1 | g)` had to be masked, and the reason is a refusal already on the books

`make_random_effects()` built `terms[h][r]` for every additive predictor, so
under a varying-coefficient model each coefficient would have got its own
group-specific offset -- a random slope, silently, in a package whose
`split_random()` refuses `(x | g)` with "asks for a random slope. Put the
variable in the fixed part of the formula instead". Group intercepts now reach
the control function only. Verified: the coefficient forest's random effects come
back exactly zero.

A group-varying coefficient remains expressible as `vc(z, ~ g + x1)`, which is a
fixed group effect on the coefficient.

## The rest of what is in

- **Categorical covariates**: one forest per level, coded symmetrically the way
  `multinomial()` codes its predictors. Mean-centered, the columns sum to zero
  within a row, which is the single function-valued redundancy the coding
  carries -- the same one symmetric multinomial coding has. `coef()` recentres
  to sum to zero across levels, exactly rather than approximately, which is what
  makes the reference a reporting choice: `center` names the level to report
  against and no refit is needed to change it. No other BART package supports a
  categorical varying covariate at all; flexBART silently coerces a factor to
  integer codes and fits a linear-in-codes model.
- **Forest names**: `(Intercept)` and the covariate's own name when the family
  has one additive parameter, `mean` and `mean:z` when it has several, so
  per-forest arguments are typable.
- **`coef()`**: the coefficient functions at each observation, on the identified
  scale, with `draws = TRUE` for the full posterior. stochtree makes you call
  `$mean_forests$predict_raw()` and rescale by hand; flexBART's return shape
  changes across four model types.
- **Prior scale**: a coefficient forest gets `sd(y)/sd(X_j)`, so a coefficient is
  shrunk the same amount whatever units its covariate is in. VCBART and flexBART
  both use a flat leaf scale and neither adapts. It turns out not to matter much
  (see above), but it is the right default for the wrong-units case nobody
  measured.
- **`bcf()`**: the causal wrapper. Fits the propensity score and puts it in the
  control function only, fewer trees for the effect, `sparsity = FALSE` on the
  outcome model per the measured guidance in `?bartisan_control`. The treatment's
  type decides the propensity model: `binomial()` for binary, `multinomial()` for
  $K$ categories with all $K$ probabilities going in, since the balancing score
  for a multi-valued treatment is the whole vector (Imbens 2000; Imai and van Dyk
  2004). A continuous treatment refuses the propensity with an explanation: its
  analogue is a conditional *density* (Hirano and Imbens 2004), which needs a
  density model. `dpm()` is one, so this is a natural extension and is not built.

## Three things found by auditing, not by testing

**The density was wrong for every varying-coefficient fit, and silently.**
`conditional_density()` handed the engine the stored basis and asked it to build
a varying-coefficient family, then passed the *combined* predictor -- so the
family expected one matrix per forest and got one. `predict(type = "density")`
died with "eta_draws must have one matrix per additive predictor", and so did
anything routed through it. Nothing in the suite covered a density on a `vc()`
fit, which is why it survived. The fix is one line of intent: the combination
already happened in R, so the density is the plain family's at that predictor.
That also makes it right under a drawn coding, where a single basis matrix could
not have carried it at all -- the stored column there is a zero placeholder, so
the version that "worked" would have priced every row as though the covariate had
no effect.

The lesson is the shape of the sweep that found it: every accessor against every
`vc()` shape, rather than the one accessor the feature was about. Five model
shapes times fifteen methods found one bug and confirmed the other seventy-four
combinations.

**An abbreviated `center` silently fitted a different model.** `arg::match_arg()`
completes abbreviations, so `center = "est"` came back as `"estimate"` inside the
`switch` -- and the guard below tested the *unmatched* string, so it fell through
to a fixed centring of zero and fitted the wrong model without a word. Matched
once and reused now.

**And the documentation demonstrated what its next paragraph forbade.** Both
`?vc` and `vignette("effects")` illustrated a coefficient varying across its own
covariate as `y ~ x1 + z + vc(z, ~ z + x1)`, with `z` in the fixed part as well
-- which is exactly the overlap the following paragraph tells the reader to avoid
and which `bartisan()` warns about. Verified by running it: the warning fires.

## Two things found by writing the vignettes, not by testing

**Augmented families discarded the wrapper.** `augmented_family()` builds a fresh
family from the name and the options, so the varying-coefficient wrapper
`make_family()` had put on was thrown away, and an augmented family reached the
sampler claiming one additive predictor while the rest of the fit expected
1 + J. Every `vc()` fit with `binomial()` failed with a confusing message about
the offset's shape. The gaussian-only tests never touched it. Fixed by wrapping
the rewritten family too, and there is now a test across `logit`, `probit` and
`poisson()`.

**A `bcf()` fit used a column the caller never named.** The propensity score is a
predictor of the control function, so `predict(fit, newdata = <their data>)`
failed on a missing `.propensity`. The fit now keeps the propensity model and
rebuilds the score. One consequence is documented rather than hidden: a predictor
goes through the quantile transform, which is a step function, so a score that
rebuilds to 1e-11 can still land on the other side of a step, and predictions for
the *training* data come back to within a few percent of the response's spread
rather than exactly. Supplying the stored score in `newdata` removes the
reconstruction and reproduces the fit to 2e-10, which is what pins the cause.

## Against VCBART

`_dev/vcbart-benchmark.R`, 8 replicates, $n = 500$, 50 trees, 1000 draws, five
coefficient functions of deliberately different kinds. Paired within replicate;
the standard error is of the paired difference.

| coefficient | RMSE, bartisan | RMSE, VCBART | difference (se) | width, bartisan | width, VCBART |
|---|---|---|---|---|---|
| $\beta_0$ control | 0.204 | 0.239 | -0.036 (0.006) | 0.840 | 1.104 |
| $\beta_1$ strong | 0.125 | 0.162 | -0.037 (0.004) | 0.655 | 0.922 |
| $\beta_2$ mild | 0.133 | 0.194 | -0.061 (0.008) | 0.568 | 0.972 |
| $\beta_3$ constant 1 | 0.107 | 0.192 | -0.085 (0.011) | 0.626 | 1.030 |
| $\beta_4$ null 0 | **0.052** | 0.200 | **-0.148 (0.015)** | **0.399** | 1.098 |

Better on every function, and the margin grows as the function gets simpler:
four times better on the null coefficient, where the intervals are also 64%
narrower while still covering (1.000 against 0.988, both above the nominal 0.95).
Coverage is at or above nominal everywhere for both, so the narrower intervals
are not bought by under-covering.

The cost is speed: 36.4 seconds per fit against 11.3, so VCBART is a bit over
three times faster. That is the soft-rule tax the package pays everywhere, and
`?bartisan_control` documents where to get it back.

## Ghosh et al. (2025), and whether sparsity already does it

Ghosh, Bhogale and Deshpande (2025) is sparseVCBART, which adds two things to
VCBART: a DART sparsity prior over the modifiers, which bartisan already has and
per forest rather than globally, and a regularized horseshoe on the leaf jumps,
where bartisan has only the local half -- a per-forest half-Cauchy `sigma_mu`,
with no global parameter tying the ensembles together and no slab. Their headline
claim is narrower and better-calibrated intervals "especially for null covariate
effects".

The first half of that is already delivered by a wide margin. On the null
coefficient above, bartisan's interval is 0.399 wide against VCBART's 1.098 --
64% narrower -- while covering at least as often. So whatever the horseshoe is
worth, it is not worth the gap Ghosh et al. measured against the VCBART they
compared to, because bartisan is already most of the way across it.

**But sparsity is not what carries that, and this is the useful finding.**
`_dev/sparsity-null.R`, 10 replicates on the same design, paired differences from
the default:

| arm | $\beta_1$ strong | $\beta_2$ mild | $\beta_3$ constant | $\beta_4$ null |
|---|---|---|---|---|
| RMSE, `sparsity = TRUE` | 0.141 | 0.131 | 0.104 | 0.047 |
| RMSE, `sparsity = FALSE` | +0.092 | +0.042 | +0.017 | +0.006 (se 0.006) |
| width, `sparsity = TRUE` | 0.669 | 0.577 | 0.615 | 0.384 |
| width, `sparsity = FALSE` | +0.473 | +0.225 | +0.090 | +0.043 (se 0.026) |

Sparsity earns its keep on the *complex* coefficient -- turning it off costs 65%
on $\beta_1$'s RMSE and 71% on its interval width -- and does almost nothing on
the null one, where the paired differences are inside two standard errors.
`"strong"` is a wash against the default everywhere and `"weak"` is worse on the
null. So the default is right and there is nothing to tune here.

**What is left on the table is conservatism, not error.** Coverage of the null
coefficient is **1.000 against a nominal 0.95 in every arm**, and no sparsity
setting moves it: the widths across all four arms span 0.377 to 0.464 while
coverage stays pinned at 1. That residual over-coverage is precisely what a
global-local prior addresses and what a per-forest local scale cannot, because
narrowing a null coefficient's interval towards nominal needs the *other*
ensembles to say that this one is null, which is what a shared global parameter
carries.

So the assessment is: worth knowing about, not worth building yet. The gain is
narrower intervals on coefficients that are already correct and already far
narrower than the comparator's, and the cost is coupling every coefficient
forest's leaf scale through a shared parameter in a design where they are
currently independent. Revisit if calibrated intervals on null coefficients
become the point rather than a by-product -- the number to beat is the 1.000.

## The forest space is two-dimensional

The plan proposed `list(mean = y ~ x1 + vc(z), log_sd = ~ x1)` and the first
implementation refused it: `expand_for_vc()` insisted on one additive predictor.
That is now built, and it is a smaller change than it looked, because the whole
of what the engine needs to know about the shape is two index vectors.

`param[h]` is the additive predictor forest `h` feeds and `column[h]` is the
basis column it is multiplied by, or none for a control function. Then

```
mu[param[h]] += slope(i, h) * eta[h]
```

for every forest, and every other method is the chain rule through that one
line: `score_info` asks the wrapped family about `param(h)` and scales by
`slope(i, h)`, and a target that is quadratic in that predictor is quadratic in
`eta_h`. The single-parameter case is the special case where `param` is all
zeros, which is why the old path came through unchanged.

Three things fell out of doing it per parameter rather than globally:

- **`center = "estimate"` is judged per parameter, and has to be.**
  `location_scale()` is quadratic in the mean and not in the log standard
  deviation, so the drawn coding's conjugate step is exact on one and would be an
  uncorrected Laplace step on the other. Asking `target_form(param(h))` accepts
  the first and refuses the second, and the refusal names the forest.
- **The duplicate-covariate check had to move.** `y ~ vc(z) + vc(z)` is a
  mistake, but the mean and the log standard deviation may each give `z` a
  coefficient -- and the union across the formulas, which is what settles the
  design, has `z` wrapped twice either way. So the check is per formula and the
  union call switches it off.
- **`num_pinned()` had to be forwarded.** A custom family's nuisance parameter is
  a trailing forest pinned at depth zero. The wrapper did not say so, so the
  engine treated it as an ordinary forest and the parameter was never drawn --
  invisible before, because a custom family with several predictors could not
  have a varying coefficient at all.

`multinomial()` and `mnp()` still refuse. Their forests are the levels of one
parameter, identified only up to a function they all share, and `report_shift()`
removes it at reporting; a coefficient forest per level would add one such
direction per coefficient and `report_shift()` does not carry them. Supporting
them means generalizing the recentring, not the forest layout.

## What is not built

- `lin()` for non-varying parametric coefficients, which was out of scope in the
  plan and remains so.
- The generalized propensity score for a continuous treatment via `dpm()`.
- Woody et al.'s heuristics for diagnosing the linearity assumption, which would
  be a good companion to the `vc(z, ~ z)` escape hatch.
- A global-local prior on the leaf jumps across coefficient forests, per Ghosh et
  al. See above for the measurement that says why it is not urgent.
- Nothing further on the forest layout; see below.
