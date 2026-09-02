# bartisan: Generalized Bayesian Additive Regression Trees

Fits Bayesian additive regression trees (BART) models for likelihoods
outside the conditionally conjugate Gaussian case, using the
Laplace-approximation reversible-jump sampler of Linero (2025).
Supported response families include Gaussian, binomial, Poisson,
negative binomial, gamma, ordinal (cumulative link), multinomial,
accelerated failure time models for right-censored survival data, and
location-scale Gaussian regression. Decision rules may be hard, as in
standard BART, or soft, as in the SoftBart model of Linero and Yang
(2018), which yields smoother fits. The interface follows that of
'glm()', so that a model is specified with a formula, a data frame, and
a family.

## Details

Fits Bayesian additive regression trees the way
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) fits a generalized
linear model: a formula, a data frame, and a family. The forest replaces
the linear predictor, so nothing has to be said about which terms enter,
which are curved, or which interact. Everything else about the workflow
stays where it was.

    fit <- bartisan(y ~ ., data = d)

## What to reach for

Most of what you will want to do with a fitted model lives in a package
that already does that job well, and *bartisan* registers the methods
those packages need rather than reimplementing them. This table is the
map.

|  |  |
|----|----|
| You want to | Use |
| fit a model | [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md) |
| choose a likelihood | [bartisan-families](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), [`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md) |
| change the sampler's settings | [`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md) |
| predict for new data | [`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md) |
| a prediction **with an interval** | [`marginaleffects::predictions()`](https://rdrr.io/pkg/marginaleffects/man/predictions.html) |
| an interval for a **new observation**, noise included | [posterior_predict()](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md) |
| how much a predictor moves the outcome | [`marginaleffects::avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html) |
| a partial dependence plot | [`marginaleffects::plot_predictions()`](https://rdrr.io/pkg/marginaleffects/man/plot_predictions.html) |
| which predictors the forest uses | [`variable_importance()`](https://ngreifer.github.io/bartisan/reference/variable_importance.md) |
| to tell the prior which predictors matter | `split_prior` in [`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md) |
| to give one forest its own predictors or settings | a list of formulas, and per-forest arguments; see [bartisan-families](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) |
| a coefficient that varies with the other predictors | [`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) in the formula, then [`coef()`](https://rdrr.io/r/stats/coef.html) |
| a treatment effect with its own prior | [`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md), [`vignette("causal")`](https://ngreifer.github.io/bartisan/articles/causal.md) |
| to check it converged and mixed | [`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md), then [as_draws()](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md) with bayesplot |
| to check it fits | [pp_check()](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md), [residuals()](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md) |
| to compare two models | [loo()](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md) |
| survival data | [`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), [`dpm_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), [`vignette("survival")`](https://ngreifer.github.io/bartisan/articles/survival.md) |
| a likelihood of your own | [`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) |

## If you are new to this

Three things are worth knowing before the first fit, and none of them
requires knowing anything about Bayesian statistics.

**There are no coefficients.** A forest has no slope to read off, so the
question "what is the effect of `x`" is answered by asking the fitted
model what it predicts under one value of `x` and under another, and
taking the difference.
`marginaleffects::avg_comparisons(fit, variables = "x")` does exactly
that, and reports an interval with it. This is a better habit than
reading coefficients even when coefficients exist, and here it is the
only habit available.

**The intervals mean what you would hope.** A 95% interval from any of
the above is the range the model considers most plausible, given the
data and the model. It already includes the uncertainty from not knowing
the shape of the relationship, which is the part a linear model leaves
out by assuming it away.

**The defaults are meant to be used.** The settings in
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
are there for people who need them; the priors and the number of trees
are chosen to work across a wide range of problems, and tuning them is
rarely where the gains are. Choosing the right
[family](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
matters much more.

## What the fit does not do for you

It is flexible about the shape of the relationship, not about anything
else. It will not tell you that a predictor is a cause, that the sample
represents the population, or that the outcome was measured well. A
forest fitted to confounded data returns a confounded answer with a
tight interval around it.

## See also

[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
to fit,
[bartisan-families](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
to choose a likelihood,
[bartisan-marginaleffects](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
and
[bartisan-interop](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
for the packages that read a fit, and
[`vignette("bartisan")`](https://ngreifer.github.io/bartisan/articles/bartisan.md)
for how the sampler works.

## Author

**Maintainer**: Noah Greifer <noah.greifer@gmail.com>
([ORCID](https://orcid.org/0000-0003-3067-7154))

Authors:

- Noah Greifer <noah.greifer@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-3067-7154))

Other contributors:

- Antonio R. Linero (Author of the FlexBart reference implementation
  from which the MCMC engine is adapted) \[contributor, copyright
  holder\]
