# Causal estimands and plot methods: what to build and why

Status: agreed in outline, `bcf_fit` class added, nothing else written.

## The problem that decides the design

A BCF's varying coefficient is a contrast on the *link* scale. On a probit or
logit fit, `coef(fit)[, "z"]` is a per-unit log-odds difference, and the average
of those is the average conditional odds ratio, which is not the marginal odds
ratio and not an estimand anyone asked for. `vignette("causal")` already says
so. A `plot()` method that drew the coefficient forest would therefore be a trap,
and the worst kind: the easy thing to plot, wrong in a way that looks fine.

So the rule the whole design follows: **estimands are computed on the response
scale by per-draw g-computation.** Predict every unit under $z = 1$ and under
$z = 0$, difference within each draw, then aggregate over units. That is
interpretable for every family, and for a difference contrast the average of the
unit-level effects *is* the marginal effect. `_dev/hetero-dpm.R`'s `ate_draws()`
is the same computation, and it is what
`marginaleffects::avg_comparisons()` does internally.

The link scale is still offered, as `type = "link"`, because it is the right
thing for inspecting effect *heterogeneity* on the scale the forest actually
models. It is never the default and its documentation has to say what its
average is not.

## Why in the package rather than through *marginaleffects*

Two reasons, neither of them that *marginaleffects* is inadequate.

- It is in **Suggests**. A causal workflow is the reason a good part of this
  package exists, and it should not need a soft dependency to report an ATE.
- What a reader wants on first contact with a fit (the ATE, both potential
  outcomes, the spread of conditional effects, and whether the arms overlap at
  all) is four or five calls today, plus care about which scale each is on. One
  call that gets it right by default is where the scale warning belongs.

The risk is divergence between two routes to the same number, and the mitigation
is a test rather than discipline: `estimate_effect(fit, "ATE")` against
`avg_comparisons()`. These should be identical, with difference arising only due to how each package computes the posterior estimate. In *bartisan* we default to the posterior mean, while *marginaleffects* uses the posterior median by default. For any checks, make sure to set *marginaleffects* to use the posterior mean.

## Architecture: one engine, two entry points

`_dev/TASKS.md` settled this pattern for `variable_importance()` and
`error_density()`: the plot method dispatches on the *result* object, and the
computing function takes a `plot` argument that **calls the method** so two
entry points to one drawing cannot drift. Same here.

### `estimate_effect()`, exported

```r
estimate_effect(object,
              estimand = "ATE",           # "ATE", "ATT", "ATC", "CATE"
              comparison = "difference",  # "difference", "ratio", "lnratio", "or", "lnor"
              by = NULL,                  # variable or formula: subgroup effects
              newdata = NULL,
              level = 0.95,
              interval = "eti",           # "eti" or "hpdi"
              focal = NULL,               # For ATT with multi-categorical treatments
              type = "response"           # prediction type
              plot = FALSE)
```

Returns `<bartisan_effect>`: a data frame of `estimate`, `lower`, `upper`, and
either a unit or a subgroup identifier, with the estimand's posterior draws kept
as an attribute so a caller can re-contrast without refitting.

The simplification that makes this one code path: each estimand is an average over a subgroup of units: ATE is all, ATT treated only, ATC control only.

For bartisan fits that are not from bcf(), the treatment variables needs to be named. Otherwise it should be extracted from the bcf object.

Two things the documentation has to be explicit about:

- For a **ratio** contrast the potential outcomes are averaged first and then
  contrasted, which is the marginal ratio. Contrasting first and averaging after
  gives the average conditional ratio, a different quantity.
- For `estimand = "CATE"` there is no averaging, so a ratio there *is* a
  conditional ratio. The print method should label it as one.
  
For multi-category treatments, if the estimand is ATT or ATC, `focal` must be supplied to identify the focal (treated) group, just as it is in *WeightIt* and *cobalt*. In that case, only display contrasts between the focal group and the other groups, but compute all pairwise comparisons. For the ATE, compute and display all pairwise comparisons. A display option (possibly in the `print()` method) controls whether to display all pairwise comparisons, only those with a specific treatment group, or only those specifically named.

### `summary()` on a `<bcf_fit>`

Calls `estimate_effect()` with defaults and prints a compact block: the ATE, both
potential outcomes, the quartiles of the conditional effects, and a line naming
the interval type. Reporting the potential outcomes beside the difference is the
point `vignette("causal")` makes about lalonde, that a few hundred dollars means
something different against a baseline of six thousand than against six hundred.

### Dispatch

Done: `bcf()` now returns `c("bcf_fit", "bartisan_fit")`. Prepended rather than
replacing, so `predict()`, `pp_check()`, the *marginaleffects* methods and the
`arg::arg_is()` guards all still dispatch. The audit found no exact-class
comparison anywhere in `R/` or `tests/`, and `test-bcf.R:208`'s
`expect_s3_class(fit, "bartisan_fit")` passes unchanged because it uses
`inherits()`.

`print()` should gain a `bcf_fit` variant naming the treatment and
the estimand it would report by default, and instructions on how to extract the treatment effect estimate with `estimate_effect()`.

## What the plot methods draw

### `plot.bartisan_effect()`

Rendering follows what was asked for, since the object records it:

| estimand | drawing |
|---|---|
| `"CATE"` | units ordered by estimate, credible bars, a line at 0, and a band behind them for the marginal estimand. This is the forest plot `vignette("causal")` builds by hand. |
| `by = ~ x` | the same forest, one row per subgroup level |
| a scalar estimand | the posterior density of the estimand with its interval marked |

### `plot.bcf_fit()`

The CATE forest by default, that being the question a BCF is asked.

### `plot.bartisan_fit()`

For a fit with no named treatment the useful thing is **partial dependence**:
`plot(fit, ~ x1)` for the posterior mean and band of $E[Y \mid x_1]$
marginalizing over the rest, `~ x1 + x2` for a surface or a family of curves.
Following the settled pattern that means a `partial_dependence()` returning
`<bartisan_partial>` with its own `plot()` method, and `plot.bartisan_fit()` as
sugar over it. This is the row `vignette("implementation")`'s comparison table
marks as covered by a helper package rather than by this one.

### Considered and not planned

- A **random-effect caterpillar** for `(1 | group)` terms is worth having and is
  cheap. Not in the first version, but no reason against it.
- **Calibration**, binned observed against fitted, largely duplicates
  `pp_check()`.
- **Trace and rank plots** belong to *bayesplot* and *posterior*. `diagnose()`
  already routes there and should keep routing rather than reimplementing.

## Scope of a first version

Binary and multi-category treatment only, with an informative error for continuous. A continuous one wants
an average slope or a dose-response curve, which is a different plot with a
different argument. A multinomial or ordinal *response* needs a category
argument if the prediction type is probabilities before "the effect" means anything, so those error too unless a different prediction type is used (e.g., stdlv or mean). A narrow
version that is right beats a general one with quiet traps.

## Tests

1. `estimate_effect(fit, "ATE")` against `marginaleffects::avg_comparisons()`. The two routes must not diverge and the check is cheap. Remember to standardize the posterior summary between the two.
2. For an identity link, the ATT equals the mean of the conditional effects among
   the treated. For a logit link, the marginal odds ratio does **not** equal the
   mean conditional odds ratio. Asserting the inequality is what keeps the scale
   bug from coming back.
3. The `plot` argument and the `plot()` method produce the same `data` and the
   same `labels`. ggplot objects carry environments and do not compare equal,
   which is why the comparison is of the pieces.
