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
original proof was the centring.

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

### 3. `center = "estimate"` is not built

Hahn et al. section 5.3 replaces the fixed coding with $b_{z_i}$ drawn from
$N(0, 1/2)$ per level, which makes the model invariant to how the treatment is
coded. The plan had it, and the analysis in the plan stands -- it extends to $K$
categories with every pairwise contrast calibrated as $N(0,1)$, and it is a
rank-one restriction for $K > 2$. It needs a conjugate Gibbs step the sampler
does not have, and it is not implemented. `center = "estimate"` is refused rather
than accepted and ignored.

Worth noting that the symmetric per-level coding for factors, which *is* built,
already delivers the invariance Hahn was after by a different route: relabelling
the levels permutes exchangeably-priored forests and leaves the prior unchanged.
What `"estimate"` would add for a factor is parsimony, not invariance. For a
numeric binary covariate it would add both.

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
  `multinomial()` codes its predictors. Mean-centred, the columns sum to zero
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

## What is not built

- `center = "estimate"`, above.
- `lin()` for non-varying parametric coefficients, which was out of scope in the
  plan and remains so.
- The generalized propensity score for a continuous treatment via `dpm()`.
- Woody et al.'s heuristics for diagnosing the linearity assumption, which would
  be a good companion to the `vc(z, ~ z)` escape hatch.
- A benchmark against VCBART's own simulation and against flexBART. The paired
  design with standard errors on the paired difference is the way to run it;
  means across a handful of replicates misled twice this session.
