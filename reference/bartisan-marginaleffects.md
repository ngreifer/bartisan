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
  `"prob"` asks for both of its columns. `"probs"`, `"lp"`, `"lv"` and
  `"surv"` are accepted as aliases for `"prob"`, `"link"`, `"link"` and
  `"survival"`, since those are the names the same quantities go by for
  other ordinal fits in marginaleffects.

- coefs:

  ignored; a forest has no coefficient vector.

## Value

`get_predict()` returns a data frame with columns `rowid`, `group` and
`estimate`, carrying the draws in a `"posterior_draws"` attribute of
observations by draws, which is the interface marginaleffects documents.
The other methods exist to satisfy the generic and return what their
names suggest.

## Details

### Uncertainty

A forest has no coefficient vector and no variance-covariance matrix, so
the delta method marginaleffects uses for a frequentist model has
nothing to work with. It has the posterior draws instead: every estimand
is computed by pushing all of the draws through the same transformation
and summarizing at the end, so an interval is a posterior quantile
rather than a normal approximation, and a nonlinear estimand needs no
approximation at all. This is the same path marginaleffects takes for
brms and rstanarm fits.

One consequence is worth knowing: marginaleffects centers a posterior at
its median, where
[`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
reports its mean. The two summarize the same draws, so a difference
between them is the skewness of the posterior rather than a
disagreement.

### Contrasts of Exactly Zero

`avg_comparisons()` reporting an estimate of exactly `0` is the most
common surprise here, and it comes from the prior and the posterior
summary acting together.

The posterior of a contrast has an atom at exactly zero. In any draw
where no tree splits on the variable being contrasted, the fit does not
depend on that variable at all, so the two counterfactual predictions
are identical to the last bit and their difference is exactly zero. The
Dirichlet sparsity prior, set by the `sparsity` argument of
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md),
makes those draws common, since dropping a weak predictor from every
tree is the purpose of a variable selection prior. And marginaleffects
centers a posterior at its median, so once the atom holds more than half
the mass the reported estimate is exactly zero however large the rest of
the posterior is: the median can land on the atom while the posterior
mean and the upper limit of the interval are both far from zero.

There are a few things worth doing about it, in the order given. Look at
`prop_used` in [`summary()`](https://rdrr.io/r/base/summary.html), the
posterior probability that the predictor appears anywhere in the forest,
the quantity the zero reports. Ask for the mean instead, with
`options(marginaleffects_posterior_center = mean)`, which is the summary
[`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
reports and the one that behaves sensibly against an atom. Reconsider
the sparsity prior, since setting `sparsity = FALSE` removes the atom
almost entirely where a larger `num_trees` leaves it in place; that is a
modeling choice rather than a fix, and it is the right one when a
contrast on a particular predictor is the estimand. And run several
chains and compare them, because the variable selection state mixes
slowly and a single chain can look much more settled than the posterior
is.
[`vignette("effects")`](https://ngreifer.github.io/bartisan/articles/effects.md)
works the choice through.

### Slopes and the Predictor Transform

A slope is a numerical derivative, and taking one requires the fitted
function to be differentiable in the predictor as the caller supplies
it. Whether it is depends on `x_transform` in
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md),
because the fit is a smooth function of the transformed predictor rather
than of the original one.

Setting `x_transform = "range"` suits slopes, and refitting with it is
worthwhile when a derivative is the quantity being reported rather than
a prediction. It is the only affine map of the three, so it is the only
one through which a slope of the fit is a slope of the original
predictor. The default, `"smoothcdf"`, carries an estimated density
along with it and lands close; `"quantile"` maps each predictor through
a step function and so has no derivative to take, and the number
returned there is a property of the step size rather than of the fit.
[`vignette("effects")`](https://ngreifer.github.io/bartisan/articles/effects.md)
shows what each does as the step shrinks.

None of this affects
[`marginaleffects::predictions()`](https://rdrr.io/pkg/marginaleffects/man/predictions.html)
or
[`marginaleffects::comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
, which evaluate the fit at two points a substantive distance apart
rather than dividing by a vanishing one.

### The Usual Survival Estimand

For a survival family, the estimand is usually a contrast in survival at
a horizon rather than in the predictor. `type = "survival"` with `times`
gives it:

    # The difference in one-year survival between treated and untreated
    avg_comparisons(fit, variables = "trt", type = "survival", times = 1)

One time per call. marginaleffects checks the dots against a whitelist
of its own, hardcoded per model class, so it warns that it does not
recognize `times` while passing it through, as the warning says. There
is no hook for registering an argument with it, so the warning is
expected and the result is correct.

### Prediction Types

Setting `type = "link"` needs a family with a single additive predictor,
since there is no one link to be talking about in
[`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
the zero-inflated families or
[`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md).
Those families work on the response scale, which is one number per
observation whatever the family, and
[`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
and
[`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
work on the probability scale, which gives one group per category.
[`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
reaches a particular predictor of a multi-predictor family, and computes
`type = "class"` and `type = "density"`, which are not one number per
observation that an average could be taken of.

Extrapolation is worth keeping in mind for `comparisons()`: a forest is
constant outside the range of the predictor it was fitted on, so a
contrast that steps a predictor beyond that range reports the boundary
value rather than an extrapolated one.

## See also

[`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
for the prediction scales these estimands are computed on;
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
for `sparsity` and `x_transform`;
[`bartisan-interop`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
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
#> ℹ Set `family` explicitly to silence this message.

# The effect of catheterization on the probability of death, as an average
# contrast of counterfactual predictions
marginaleffects::avg_comparisons(fit, variables = "rhc")
#> 
#>  Estimate 2.5 % 97.5 %
#>    0.0453     0  0.101
#> 
#> Term: rhc
#> Type: response
#> Comparison: 1 - 0
#> 

# The same contrast within each sex, and a test that the two are equal
marginaleffects::avg_comparisons(fit, variables = "rhc", by = "sex",
                                 hypothesis = ~pairwise)
#> 
#>         Hypothesis Estimate   2.5 %  97.5 %
#>  (male) - (female)        0 -0.0388 0.00448
#> 
#> Type: response
#> 

# Centering each posterior at its mean rather than its median, which is the
# summary `predict()` reports
op <- options(marginaleffects_posterior_center = mean)
marginaleffects::avg_comparisons(fit, variables = "rhc")
#> 
#>  Estimate 2.5 % 97.5 %
#>    0.0463     0  0.101
#> 
#> Term: rhc
#> Type: response
#> Comparison: 1 - 0
#> 
options(op)
```
