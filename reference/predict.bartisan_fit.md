# Predictions from a generalized BART model

Evaluates the stored posterior draws of the forests, either at the data
used to fit the model or at new data. Because every draw of every tree
is kept, predictions carry full posterior uncertainty rather than being
a single point estimate.

## Usage

``` r
# S3 method for class 'bartisan_fit'
predict(
  object,
  newdata = NULL,
  type = "response",
  draws = FALSE,
  iterations = NULL,
  offset = NULL,
  weights = NULL,
  values = NULL,
  log = FALSE,
  times = NULL,
  ...
)
```

## Arguments

- object:

  a fitted model from
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).

- newdata:

  optional data frame at which to predict. Omit it to use the data the
  model was fit to. Missing predictor values are allowed in the columns
  that had them when the model was fit, since only those columns'
  splitting rules carry an answer for one; a missing value anywhere else
  is an error. See the Missing predictor values section of
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).

- type:

  the scale of the prediction:

  `"link"`

  :   the additive predictor, one column per predictor for families that
      have more than one.

  `"response"`

  :   the mean of the response; the median survival time for the
      accelerated failure time families; and, for a response with
      categories, the category probabilities, since there is no single
      mean to report.

  `"prob"`

  :   category probabilities, for the binomial, ordinal and multinomial
      families.

  `"class"`

  :   the most probable category, as a factor, for the same families.

  `"mean"`

  :   the mean of the response with the category labels read as numbers,
      for the same families. `"4"` counts as four. This is the summary
      an ordinal outcome with numeric labels usually wants, and it needs
      no assumption at the modeling stage: the model is still ordinal
      and only the reporting treats the categories as numbers. Use
      `values` to say what the categories are worth when the labels are
      not numbers, or are not the numbers you mean.

  `"stdlv"`

  :   the additive predictor divided by the standard deviation of the
      latent variable it indexes, for the ordinal and binomial families.
      Either response can be written as a threshold crossing of a
      continuous `y* = eta + e`, and the link fixes the distribution of
      `e` and so its variance; dividing by the standard deviation of
      `y*` puts fits with different links, or different amounts of
      signal, on one scale, which is what a standardized effect size on
      such an outcome needs. Available for the probit, logit and
      complementary log-log links, which are the ones with a latent
      distribution to name. See Details.

  `"density"`

  :   the conditional density of the outcome given the predictors,
      evaluated at the observed outcome. This requires the outcome, so
      `newdata` must contain it; omit `newdata` to use the data the
      model was fit to. The value is the likelihood contribution of the
      observation, so it is a density for a continuous response, a
      probability for a discrete one, and a survival probability for a
      censored survival time. Useful for held-out log scores and for
      posterior predictive checks. **The measure differs across the
      survival families**: the accelerated failure time families,
      [`dpm_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
      included, report the density of \\\log T\\, while
      [`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
      reports the density of \\T\\. The two differ by \\\sum \log t\\,
      so log scores are comparable within each group and not across
      them; `type = "survival"` is comparable throughout.

  `"survival"`

  :   the survival function \\S(t \mid x)\\ at the times given in
      `times`, for the accelerated failure time families and
      [`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md).
      This is the predictive distribution of a survival response, and
      the analogue of `"prob"` for a categorical one: `"response"`
      reports only the median. Returns one column per time, or a draws
      by rows by times array when `draws = TRUE`. It is also what makes
      the usual survival estimand reachable through marginaleffects – a
      contrast in \\t\\-year survival; see
      [bartisan-marginaleffects](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md).

- draws:

  return every posterior draw rather than the posterior mean. The result
  then gains a leading dimension indexing draws.

- iterations:

  optional integer vector selecting which stored draws to use. Defaults
  to all of them.

- offset:

  an offset for `newdata`, on the link scale. Required when the model
  was fit with an observation-level offset, because the offset is not a
  function of the predictors and so cannot be reconstructed.

- weights:

  prior weights for `newdata`, used only by `type = "density"`. For a
  binomial response given as proportions these are the numbers of
  trials. Ignored otherwise.

- values:

  for `type = "mean"`, a numeric vector named for every response level,
  giving what each category is worth. Defaults to the level labels read
  as numbers, which fails with an error rather than a guess when they
  are not numbers.

- log:

  for `type = "density"`, return the log of the value instead.

- times:

  for `type = "survival"`, the times at which to report the survival
  function. Required, because the horizon is a choice rather than a
  property of the fit. Summing across observations then gives a log
  score. With `draws = FALSE` the density is averaged over draws before
  the log is taken, so the result is the pointwise predictive density
  rather than the average log density.

- ...:

  ignored, present for compatibility with the generic.

## Value

With `draws = FALSE`, a vector for a single-predictor family on the
`"link"` or `"response"` scale, a matrix of observations by categories
for `"prob"`, a factor for `"class"`, and a matrix of observations by
predictors otherwise. With `draws = TRUE`, a matrix of draws by
observations, or a list of such matrices when the family has several
additive predictors, or an array of draws by observations by categories
for `"prob"`.

## Multinomial probit probabilities are simulated

The likelihood of a `multinomial("probit")` fit has no closed form: the
probability of a category is the chance that the largest of several
correlated Gaussian variables is the one belonging to it, which is a
multivariate orthant probability. So `"prob"`, `"class"`, `"response"`
and `"density"` are all simulated, using the number of replicates the
family was given. Fresh draws are taken on every call, so two calls
differ by Monte Carlo error; that error is per posterior draw and
averages down over them, which makes `draws = FALSE` much more accurate
than any single row of `draws = TRUE`.

## The standardized latent variable

`type = "stdlv"` reports `(eta - E[e]) / sd(y*)` for the latent
`y* = eta + e`, following
[`WeightIt::predict.ordinal_weightit()`](https://ngreifer.github.io/WeightIt/reference/predict.glm_weightit.html)
. Three parts of that need saying.

The **scale** is `sd(y*) = sqrt(var(eta) + var(e))`, where `var(eta)` is
taken over the sample the model was fitted to, per draw, so it is a
property of the model rather than of whatever is being predicted; the
same divisor is used when predicting new data. `var(e)` is whatever the
link implies: 1 for the probit link, `pi^2 / 3` for the logit,
`pi^2 / 6` for the complementary log-log.

The **location** subtracts the latent error's mean, which shifts `y*` so
that its error is centered. That is invisible for the logit and probit
links, whose errors are already centered, and is the whole of the
difference for the complementary log-log link, whose error is a smallest
extreme value variate.

The **sign of that shift differs between the two families**, and only
for the complementary log-log link. A normal or logistic error is
symmetric, so it does not matter whether `e` or `-e` is the thing added
to the index. A smallest extreme value error is not symmetric, and the
two families add it with opposite signs: an ordinal model has
`P(Y <= k) = G(c_k - eta)`, which is `P(eta + e <= c_k)`, so its error
has mean `-gamma`; a binomial model has `P(Y = 1) = G(eta)`, which is
`P(e <= eta)`, so its latent is `eta - e` and the error has mean
`+gamma`. The two are different models rather than the same one written
twice, which is also why a two-category ordinal complementary log-log
fit is not the same as a binomial one.

Such a model is identified only up to a common shift of its thresholds
and its predictor, so the location of this quantity is a convention
rather than a fact, and the one used here is the same one the cutpoints
use – a predictor centered over the fitted sample. Against WeightIt,
which identifies by dropping the intercept column instead, the two agree
on the scale and differ by a constant; measured on a linear truth, the
standard deviations agree to under 1% and the difference is constant to
three decimals. Differences on this scale, which is what a standardized
quantity is for, are unaffected.

## See also

[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)

## Examples

``` r
set.seed(1)

n <- 150
d <- data.frame(x1 = runif(n), x2 = runif(n))
d$y <- rpois(n, exp(0.5 + d$x1))

fit <- bartisan(y ~ x1 + x2, data = d, family = poisson(),
               control = bartisan_control(num_trees = 10, num_burn = 50,
                                         num_draws = 50, verbose = FALSE))

head(predict(fit))
#> [1] 2.676084 2.948299 3.282852 3.858021 2.468973 3.602617

nd <- data.frame(x1 = c(0.1, 0.9), x2 = c(0.5, 0.5))
predict(fit, newdata = nd)
#> [1] 2.306964 3.882649

# Posterior uncertainty for the two new points
apply(predict(fit, newdata = nd, draws = TRUE), 2,
      quantile, c(0.025, 0.975))
#>           [,1]     [,2]
#> 2.5%  1.510869 3.242627
#> 97.5% 3.262926 4.795421

# Held-out log score: the conditional density needs the outcome, so it must
# be present in `newdata`.
held_out <- data.frame(x1 = runif(20), x2 = runif(20))
held_out$y <- rpois(20, exp(0.5 + held_out$x1))
sum(log(predict(fit, newdata = held_out, type = "density")))
#> [1] -33.83554
```
