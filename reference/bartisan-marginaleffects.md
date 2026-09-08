# Counterfactual estimands with marginaleffects

Methods that let a `<bartisan_fit>` object be used by the
marginaleffects package, which computes predictions, comparisons and
slopes (and hypothesis tests on any of them) from the posterior draws.
There is nothing to set up: load marginaleffects and call its functions
on the fit.

## Usage

``` r
# S3 method for class 'bartisan_fit'
formula(x, ...)

# S3 method for class 'bartisan_fit'
terms(x, ...)

# S3 method for class 'bartisan_fit'
model.frame(formula, ...)

# S3 method for class 'bartisan_fit'
nobs(object, ...)

# S3 method for class 'bartisan_fit'
family(object, ...)

# S3 method for class 'bartisan_fit'
get_predict(model, newdata = NULL, type = NULL, ...)

# S3 method for class 'bartisan_fit'
get_group_names(model, ...)

# S3 method for class 'bartisan_fit'
get_coef(model, ...)

# S3 method for class 'bartisan_fit'
set_coef(model, coefs, ...)

# S3 method for class 'bartisan_fit'
get_vcov(model, ...)

# S3 method for class 'bartisan_fit'
get_data(x, ...)
```

## Arguments

- ...:

  further arguments. `values`, `iterations`, `offset`, `weights`, `log`
  and `times` are passed on to
  [`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md);
  anything else is ignored, since marginaleffects puts arguments of its
  own here too. Note that marginaleffects warns that it does not
  recognize `values`, which is expected (it is this package's argument,
  not one of its own), and the value is used regardless.

- model, x, object, formula:

  a `<bartisan_fit>` object; the output of a call to
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).
  The last is named `formula` only because
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html)
  names its first argument that way.

- newdata:

  optional; a data frame at which to evaluate the fit. Default is `NULL`
  to use the data the model was fitted to, which is retained in the fit
  for this purpose.

- type:

  string; the scale to work on. Allowable options include `"response"`
  (the default), `"link"`, `"prob"`, `"mean"`, `"stdlv"`, and
  `"survival"`. `"response"` is the fitted mean, `"link"` the additive
  predictor, `"prob"` the per-category probabilities of a categorical
  family, `"mean"` the mean of a categorical response with its labels
  read as numbers, `"stdlv"` the standardized latent variable of an
  ordinal fit, and `"survival"` the survival function at the times given
  in `times`; all but the first two are described under
  [`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md).
  An
  [`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  or
  [`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  response has no single mean, so `"response"` gives `"prob"` for those
  families, while a binomial fit reports the probability of its second
  level, as [`stats::glm()`](https://rdrr.io/r/stats/glm.html) does, and
  `"prob"` is what asks for both of its columns. `"probs"`, `"lp"`,
  `"lv"` and `"surv"` are accepted as aliases for `"prob"`, `"link"`,
  `"link"` and `"survival"`, since those are the names the same
  quantities go by for other ordinal fits in marginaleffects.

- coefs:

  ignored; a forest has no coefficient vector.

## Value

`get_predict()` returns a data frame with columns `rowid`, `group` and
`estimate`, carrying the draws in a `"posterior_draws"` attribute of
observations by draws, which is the interface marginaleffects documents.
The other methods exist to satisfy the generic and return what their
names suggest.

## Details

### How the Uncertainty Is Computed

A forest has no coefficient vector and no variance-covariance matrix, so
the delta method marginaleffects uses for a frequentist model has
nothing to work with. It has something better here: the posterior draws.
Every estimand is computed by pushing all of the draws through the same
transformation and summarizing at the end, so an interval is a posterior
quantile rather than a normal approximation, and a nonlinear estimand
needs no approximation at all. This is the same path marginaleffects
takes for brms and rstanarm fits.

One consequence worth knowing: marginaleffects centers a posterior at
its **median**, where
[`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
reports its **mean**. The two are summarizing the same draws, so a
difference between them is the skewness of the posterior and not a
disagreement.

### A Contrast of Exactly Zero Is Usually Real

`avg_comparisons()` reporting an estimate of exactly `0` is the most
common surprise here, and it is neither package computing anything
wrong. Two facts meet to produce it.

The posterior of a contrast has an **atom at exactly zero**. In any draw
where no tree in the forest splits on the variable being contrasted, the
fit does not depend on that variable at all, so the two counterfactual
predictions are identical to the last bit and their difference is
exactly zero. That is not a near-zero value that rounding flattered; it
is a point mass. The Dirichlet sparsity prior on the splitting
proportions, `sparsity` in
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md),
is what makes those draws common: it is a variable selection prior, and
dropping a weak predictor from every tree is what it is for.

And marginaleffects centers a posterior at its **median**. So once the
atom holds more than half the mass, the reported estimate is exactly
zero however large the rest of the posterior is. On `MatchIt::lalonde`
with the default settings, `treat` was absent from all 50 trees in 64%
of draws and the contrast came out exactly zero in 65%, which put the
median at 0 while the posterior mean was 197 and the upper limit was
above 2000.

Four things are worth doing about it, in the order given.

**Look at the inclusion probability**, which is what the zero reports.
[`summary()`](https://rdrr.io/r/base/summary.html) gives it as
`prop_used`: the posterior probability that each predictor group appears
anywhere in the forest. A contrast whose median is zero is a predictor
the model is not sure belongs.

**Ask for the mean instead**, with
`options(marginaleffects_posterior_center = mean)`. The mean is the
summary
[`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
reports, and it is the one that behaves sensibly against an atom.

**Reconsider the sparsity prior** when variable selection is not what
the fit is for. `sparsity = FALSE` in
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
removes the atom almost entirely: on the example above it fell from 20%
of draws to none by 50 trees. A larger `num_trees` does *not* remove it,
which is worth stating because it looks as though it should: with the
prior on, the contrast was exactly zero in 20% of draws at 50 trees and
18% at 200. Turning the prior off is a modeling choice rather than a
fix, so it wants a reason; it is the right one when a contrast on a
particular predictor is the estimand, and the wrong one when there are
many predictors and most are irrelevant.

**Run several chains and compare them.** The variable selection state
mixes slowly, because a predictor whose splitting proportion has gone
small is rarely proposed and so is hard to get back in. On the example
above, four chains disagreed by more than 100% of the estimate at every
tree count from 20 to 200 with the prior on; with `sparsity = FALSE` and
50 trees they agreed to within 9%. A single chain can look much more
settled than the posterior is.

### Slopes Need a Linear Predictor Transform

A slope is a numerical derivative, and taking one requires the fitted
function to be differentiable in the predictor *as the caller supplies
it*. The default `x_transform = "quantile"` maps each predictor through
its empirical distribution function before any rule sees it, and an
empirical distribution function is a step function; the fit is therefore
a step function of the original predictor whatever the decision rules
are, and its difference quotient grows without bound as the step
shrinks. Measured on a smooth surface where the average derivative is
zero:

|      |                            |                         |
|------|----------------------------|-------------------------|
| step | `x_transform = "quantile"` | `x_transform = "range"` |
| 1e-4 | -4.79                      | -0.28                   |
| 1e-2 | -0.40                      | -0.30                   |
| 5e-2 | -0.29                      | -0.25                   |

So **`x_transform = "range"` is what slopes want**, since it maps each
predictor linearly and leaves a soft-rule fit differentiable. Hard rules
give a piecewise-constant fit under either transform, and a derivative
of one is not a meaningful quantity however it is computed.

None of this affects
[`marginaleffects::predictions()`](https://rdrr.io/pkg/marginaleffects/man/predictions.html)
or
[`marginaleffects::comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
, which evaluate the fit at two points a substantive distance apart
rather than dividing by a vanishing one. Those are the estimands to
reach for with the default transform.

### The Usual Survival Estimand

For a survival family, the estimand is usually a contrast in survival at
a horizon rather than in the predictor. `type = "survival"` with `times`
gives it:

    # The difference in one-year survival between treated and untreated
    avg_comparisons(fit, variables = "trt", type = "survival", times = 1)

One time per call. marginaleffects checks the dots against a whitelist
of its own, hardcoded per model class, so it warns that it does not
recognize `times` while passing it through, which is what the warning
says. There is no hook for registering an argument with it, so the
warning is expected and the result is correct.

### What Is Not Covered

`type = "class"` and `type = "density"` are not available, because
neither is one number per observation that an average or a contrast
could be taken of: a class is a factor, and a density needs the outcome,
which a counterfactual grid does not have.
[`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
computes those.

`type = "link"` is refused for a family with more than one additive
predictor
([`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
the zero-inflated families,
[`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)),
because there is no single link there to be talking about. Those
families work on the response scale, which is one number per observation
whatever the family, and
[`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
and
[`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
work on the probability scale, which gives one group per category. To
reach a *particular* predictor of a multi-predictor family, call
[`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
directly.

Extrapolation is worth keeping in mind for `comparisons()`: a forest is
constant outside the range of the predictor it was fitted on, so a
contrast that steps a predictor beyond that range reports the boundary
value rather than an extrapolated one.

## See also

[`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
for the prediction scales these estimands are computed on;
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
for `sparsity` and `x_transform`;
[bartisan-interop](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
for the methods that let other packages assess the fit

## Examples

``` r
data("rhc")
set.seed(123)

# Whether a patient died, with every other variable a candidate predictor.
# The sparsity prior is turned off because a contrast on one predictor is
# the estimand rather than variable selection
fit <- bartisan(death ~ . - days, data = rhc, sparsity = FALSE,
                num_trees = 10, num_burn = 50, num_draws = 50)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` to choose another, which also silences this message.

# The effect of catheterization on the probability of death, as an average
# contrast of counterfactual predictions
marginaleffects::avg_comparisons(fit, variables = "rhc")
#> 
#>  Estimate  2.5 % 97.5 %
#>    0.0749 0.0211  0.112
#> 
#> Term: rhc
#> Type: response
#> Comparison: 1 - 0
#> 

# The same contrast within each sex, and a test that the two are equal
marginaleffects::avg_comparisons(fit, variables = "rhc", by = "sex",
                                 hypothesis = ~pairwise)
#> 
#>         Hypothesis Estimate    2.5 %  97.5 %
#>  (male) - (female) 0.000296 -0.00827 0.00302
#> 
#> Type: response
#> 

# Centering each posterior at its mean rather than its median, which is the
# summary `predict()` reports
op <- options(marginaleffects_posterior_center = mean)
marginaleffects::avg_comparisons(fit, variables = "rhc")
#> 
#>  Estimate  2.5 % 97.5 %
#>     0.072 0.0211  0.112
#> 
#> Term: rhc
#> Type: response
#> Comparison: 1 - 0
#> 
options(op)
```
