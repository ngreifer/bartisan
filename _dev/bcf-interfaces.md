# Varying-coefficient BART: what works, and four interfaces to choose among

## The model

$$g(\mu_i) = f_0(x_i) + \sum_{j=1}^{J} a_{ij}\, f_j(x_i)$$

Each $f_j$ is its own forest. The $a_j$ enter linearly with a coefficient that
varies with $x$. With $J = 1$ and $a_1$ a binary treatment this is the Bayesian
causal forest of Hahn, Murray and Carvalho (2020): a prognostic forest $f_0$ and
a treatment-effect forest $f_1 = \tau$, each with its own prior.

What that buys over putting $z$ among the predictors of a single forest is a
prior on the effect surface directly. A single forest shrinks the *fitted
surface* toward a constant, which shrinks $\tau$ only incidentally and by an
amount nobody chose.

## It works, and the proof cost nothing structural

`_dev/bcf-proof.R` writes the likelihood as a `custom_family()` with
`num_predictors = 2`, computes $\mu_i = \eta_{i1} + z_i \eta_{i2}$ inside it,
and supplies the chain rule as the analytic derivatives:

```r
derivatives = function(y, eta, h) {
  z <- y %/% 2
  prob <- plogis(eta[, 1] + z * eta[, 2])
  slope <- if (h == 1L) 1 else z          # d mu / d eta_h
  list(score = (y %% 2 - prob) * slope,
       info  = prob * (1 - prob) * slope^2)
}
```

n = 1000, confounded assignment, `tau(x) = 1 + 0.8 * x3`, two chains,
`num_trees = c(50, 25)`:

| | cor with true tau | RMSE | ATE (truth 0.981) |
|---|---|---|---|
| varying coefficient | 0.986 | 0.150 | 0.964 [0.827, 1.102] |
| one forest with `z` among its predictors | 0.969 | 0.233 | 0.967 [0.835, 1.097] |

The effect surface is recovered, and recovered better than by the thing a user
would do today. RMSE falls by a third. The binomial version behaves the same
way on the logit scale (0.859 vs 0.796, RMSE 0.498 vs 0.610).

**The device the proof needed, and the gap it exposes.** A custom likelihood is
called once per leaf with the rows reaching that leaf and is handed `y` for
those rows and nothing else, so a treatment vector held outside cannot be lined
up with them. The proof carries `z` inside the response and unpacks it on
arrival. That is a trick, and the reason it was necessary is a real hole:
**`custom_family()` has no way to see per-observation covariates.** Any censored
or truncated likelihood a user writes by hand hits the same wall. Worth closing
on its own merits, separately from this.

**Mixing: the varying coefficient is not the problem.** The proof reported rhat
around 1.39 on the prognostic forest and its leaf scale, with the
treatment-effect forest at 1.10. That looked like the weak separation of the two
surfaces that the BCF literature warns about, and it is not.
`_dev/bcf-mixing.R` puts the same custom-family route with and without the
varying coefficient against the compiled `gaussian()` on the same data, four
chains each:

| Fit | worst `sigma_mu` rhat | worst `eta` rhat |
|---|---|---|
| custom, one forest, no varying coefficient | 1.54 | 1.22 |
| custom, two forests, varying coefficient on `z` | 1.39 | 1.37 |
| compiled `gaussian()`, `z` an ordinary predictor | 1.36 | 1.35 |

The compiled family with no varying coefficient anywhere mixes exactly as badly.
Whatever this is, it belongs to the design or to the leaf-scale sampler, and
adding a varying coefficient does not make it worse. `sigma_mu` is the worst row
in all three at `ess_bulk` under 10 out of 3000 draws, which was followed up
separately in `_dev/sigma-mu-mixing.R` and `_dev/sigma-mu-cause.R`; it is not a
varying-coefficient question.

Adding an estimated propensity score to the covariates did not help either
(rhat 1.588, tau recovery unchanged at 0.986 / 0.152). That is not evidence
against Hahn et al.'s recommendation, because nothing here can give a covariate
to the prognostic forest and withhold it from the effect forest, which is what
they actually recommend. It is evidence that the recommendation cannot be tested
without design question (2) below being answered.

Two smaller findings, both filed:

- `derivatives` is called as `f(y, eta, h)` with no `aux`, so a custom family
  with a nuisance parameter cannot supply analytic derivatives that depend on
  it. The Gaussian half of the proof differences instead, at three R calls per
  leaf visit rather than one.
- `warn_runaway_scale()` compared `out$sigma_mu`, which has one column per
  reported forest, against a target with one entry per forest the engine builds,
  which includes the pinned forests standing in for nuisance parameters. Two
  predictors plus one nuisance parameter warned about the length; one predictor
  plus one nuisance parameter silently compared against the wrong target. Fixed.

## What the engine already gives us

Two things make this cheap to build properly, and both were checked rather than
assumed.

**The decorator shape already exists.** `LinkFamily` in `src/family.cpp` wraps
an inner family, transforms the predictor on the way in, and applies the chain
rule on the way out:

```cpp
inner->score_info_block(idx, n, th.memptr(), 0, d1, d2);
arma::vec tp = slope(block, n);
for (int k = 0; k < n; k++) { d1[k] *= tp(k); d2[k] *= tp(k) * tp(k); }
```

A varying-coefficient decorator is the same object with a different map:
present $1 + J$ predictors outward, hand $\mu$ inward, and scale $d_1$ by
$a_{ij}$ and $d_2$ by $a_{ij}^2$. Because it wraps the family rather than
replacing it, **every existing family gets varying coefficients with no
per-family work** -- `gaussian()`, `binomial()`, `dpm()`, `ph()`, the AFT
families, all of them.

**The fast leaf draws survive, for the families that have them.** The map is
linear in each $\eta_j$, so a target that is quadratic in $\mu$ is quadratic in
each $\eta_j$. `TARGET_QUADRATIC` is preserved, which means `gaussian()`,
`dpm()`, and every augmented family (binomial, ordinal, multinomial,
zero-inflated) keep their closed-form leaf draws at acceptance one. The
exponential form does not survive: `exp_rate()` returns one scalar per
predictor, and a varying coefficient makes the rate $a_i$, which varies by
observation. `poisson()` and `Gamma()` would fall back to `TARGET_GENERAL` and
pay the five-to-ten times that costs. Fixable by generalizing `exp_rate()` to a
vector, and not on the critical path.

## What any interface has to answer

1. Which variables get varying coefficients.
2. What each forest splits on. All the same covariates, or one set per forest?
   Hahn et al. specifically want an estimated propensity score in $f_0$ and not
   in $\tau$, so "one set per forest" is not hypothetical.
3. Per-forest settings. `num_trees` already takes a vector; `k`, `sparsity` and
   `sigma_mu` do not.
4. Composition with `family`, with `(1 | group)`, and with the families that
   already use several additive predictors.
5. What `predict()` and *marginaleffects* do with it.

On (5) there is a good answer available for free under every option below. If
the varying-coefficient combination happens inside the model rather than in the
family, then `predict()` on modified `newdata` is well defined, so
`avg_comparisons(fit, variables = "z")` gives the ATE with a credible interval,
`newdata = subset(z == 1)` gives the ATT, `by =` gives subgroups, and
`hypothesis = ~pairwise` tests moderation. Everything in `vignette("causal")`
keeps working, and $\tau(x)$ is also readable directly as the second forest.

On (4), the multi-predictor families are the one genuine conflict. A decorator
that presents $\mu$ to the inner family assumes the inner family wants one
predictor. `location_scale()`, `multinomial()` and the zero-inflated families
want several, and "a varying coefficient on which one" has no obvious default.
Refusing the combination at first is defensible.

---

## Option A: a bar in the formula

```r
bartisan(y ~ z | x1 + x2 + x3, data = d, family = binomial())
bartisan(y ~ z1 + z2 | x1 + x2 + x3, data = d)
```

Left of the bar, the varying coefficients; right of the bar, the covariates,
which are also what $f_0$ splits on. This is the `Formula` package's convention
and will read correctly to anyone who has used *plm* or *betareg*.

Against it:

- `|` is already spoken for. `(1 | group)` is random effects, and this package
  supports them. `y ~ z | x1 + x2 + (1 | g)` is ambiguous to read even where it
  is unambiguous to parse, and `y ~ z | x1 + (1 | g) + x2` is worse.
- No room for per-coefficient modifiers. Every $f_j$ splits on the same
  right-hand side, so question (2) is answered "no" and the propensity-score
  recommendation cannot be expressed.
- No room for per-coefficient settings without a separate argument keyed by
  position.

## Option B: a marker on the term

```r
bartisan(y ~ x1 + x2 + x3 + vc(z), data = d, family = binomial())

# per-forest modifiers and settings, in the place the term is named
bartisan(y ~ x1 + x2 + x3 + ps_hat + vc(z, ~ x1 + x2 + x3, num_trees = 25),
         data = d)

# more than one
bartisan(y ~ x1 + x2 + vc(z1) + vc(z2), data = d)
```

Unmarked terms are what $f_0$ splits on, as now. `vc(a)` says that `a`'s
coefficient is a forest, splitting by default on the same covariates and
optionally on a named subset. Extra arguments are that forest's settings.

This is how `s()` works in *mgcv*, `strata()` in *survival*, `offset()` in base
R. It answers all five questions in the place the term is written, composes with
`(1 | group)` without thinking about it, and leaves the rest of the formula
meaning exactly what it means today.

Against it: `vc(z)` inside a formula has to be recognized before
`model.frame()`, alongside the random-effect split that already happens there.
That machinery exists, and this is the second customer for it, but it is the
option with the most parsing.

## Option C: an argument

```r
bartisan(y ~ x1 + x2 + x3, varying = ~ z, data = d, family = binomial())

bartisan(y ~ x1 + x2 + x3 + ps_hat, data = d,
         varying = list(z = ~ x1 + x2 + x3))
```

No formula parsing at all, and the simple case is one short argument. The cost
is that the model specification now lives in two places, which is the thing a
formula interface exists to prevent, and the list form for per-forest modifiers
is noticeably clumsier than Option B's.

Cheapest of the four to build. A reasonable staging post even if the eventual
answer is B.

## Option D: a separate function

```r
bcf(y ~ x1 + x2 + x3, treatment = ~ z, data = d, family = binomial(),
    propensity = TRUE)
```

A named entry point for the causal case specifically, which could do the things
BCF wants done and `bartisan()` should not do by default: fit the propensity
score, put it in the prognostic forest only, use fewer trees for $\tau$, and
return an object whose `summary()` reports an ATE rather than a forest.

Against it as the *only* interface: it duplicates every argument of
`bartisan()`, and it answers the causal question while leaving the general
varying-coefficient model unreachable. A varying coefficient is not always a
treatment effect.

The better reading of D is as a **wrapper over B or C** rather than an
alternative to them. That also keeps the general machinery honest, because the
wrapper has to be expressible in it.

## A more radical option, for completeness

Make the additive structure explicit and drop the special case:

```r
bartisan(y ~ forest(x1 + x2 + x3) + z:forest(x1 + x2), data = d)
```

Fully general, and it says what the model is rather than naming a pattern. It
also abandons the `glm()`-shaped formula that is the package's whole premise,
and it makes the common case unreadable. Recorded and not recommended.

---

## Recommendation

**B, with D as a thin wrapper over it.**

B is the only option that answers question (2), and question (2) is not a
detail: the propensity-score-in-the-prognostic-forest recommendation is the main
practical advice in the BCF literature, and the proof could not test it because
the current interface cannot express it. B also costs nothing at the other four
questions and reads like the rest of the package.

C is worth building first if the parsing in B turns out to be more than it
looks, since the two are compatible: `varying = ~ z` and `vc(z)` can coexist,
and C's list form is the escape hatch B would not need.

Suggested order:

1. The `VaryingCoefficient` decorator in C++, on the `LinkFamily` pattern. This
   is the whole of the modeling work and it is family-agnostic.
2. Per-observation covariates for `custom_family()`. Independent of the above,
   and the proof showed it is a hole regardless.
3. Whichever of B or C you pick, plus `predict()` and the *marginaleffects*
   hooks, which mostly follow from the combination happening inside the model.
4. `exp_rate()` as a vector, so `poisson()` and `Gamma()` do not fall off the
   fast path.
5. D, and a section in `vignette("causal")`.

## Open questions

- Should `vc(z)` put `z` in $f_0$'s splitters as well? Hahn et al. say no for
  the treatment and yes for its propensity score. Default no, overridable.
- Should the varying coefficients be centered, so $f_0$ is interpretable as the
  surface at $a = 0$? Only matters for reporting.
- What happens for `location_scale()`, `multinomial()` and the zero-inflated
  families. Refuse, or apply to the first predictor and say so?
- Does `sparsity` operate per forest? It already does, and for BCF that is
  probably right: the prognostic and effect surfaces need not use the same
  covariates.
