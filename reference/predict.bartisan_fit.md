# Predictions from a generalized BART model

Computes predictions from a
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
fit, at the data the model was fit to or at new data, on any of several
scales. Because every posterior draw of every tree is retained, a
prediction can be returned as its posterior mean or as the draws
themselves.

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

  a `<bartisan_fit>` object; the output of a call to
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).

- newdata:

  optional; a data frame at which to predict. Default is `NULL` to use
  the data the model was fit to. Missing predictor values are allowed in
  the columns that had them when the model was fit, since only those
  columns' splitting rules carry an answer for one; a missing value
  anywhere else is an error. See section *Missing Predictor Values* in
  the Details of
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).

- type:

  string; the scale of the prediction. Allowable options include
  `"link"`, `"response"` (the default), `"prob"`, `"class"`, `"mean"`,
  `"stdlv"`, `"density"`, `"survival"`, and `"hazard"`.

  `"link"`

  :   the additive predictor, one column per predictor for families that
      have more than one. For the accelerated failure time families it
      is on the scale of log time, so a difference between two units is
      a log time ratio. For
      [`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
      it is reported centered over the fitted sample, so a single value
      is a log hazard ratio against a unit whose predictor sits at that
      average and a difference between two units is the log hazard ratio
      between them.

  `"response"`

  :   the mean of the response; the median survival time for every
      survival family,
      [`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
      included; and, for a response with categories, the category
      probabilities, there being no single mean to report.

  `"prob"`

  :   category probabilities, for the binomial, ordinal and multinomial
      families.

  `"class"`

  :   the most probable category, as a factor, for the same families.

  `"mean"`

  :   the mean of the response with the category labels read as numbers,
      for the same families, so that `"4"` counts as four. The `values`
      argument says what the categories are worth where the labels are
      not the numbers intended.

  `"stdlv"`

  :   the additive predictor on the scale of the latent variable it
      indexes, for the ordinal and binomial families under the probit,
      logit and complementary log-log links. This puts fits with
      different links, or different amounts of signal, on one scale. See
      Details.

  `"density"`

  :   the conditional density of the outcome given the predictors,
      evaluated at the observed outcome, so the `newdata` argument must
      carry the outcome and leaving it empty uses the fitted data. This
      is the observation's likelihood contribution, so it is a
      probability for a discrete response and a survival probability for
      a censored time. See Details.

  `"survival"`

  :   the survival function \\S(t \mid x)\\ at the times given in
      `times`, for the accelerated failure time families and
      [`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md).
      Returns one column per time, or a draws by rows by times array
      when `draws = TRUE`. This makes the usual survival estimand
      reachable through marginaleffects; see
      [`bartisan-marginaleffects`](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md).

  `"hazard"`

  :   the hazard \\h(t \mid x) = f(t \mid x) / S(t \mid x)\\ at the
      times given in `times`, the rate at which the event occurs at
      \\t\\ among those who have survived to it, for the same families
      and in the same shape as `"survival"`. The cumulative hazard is
      \\-\log S(t \mid x)\\ and so comes from `"survival"`.

- draws:

  `logical`; whether to return every posterior draw rather than the
  posterior mean. Default is `FALSE` to return the mean. If `TRUE`, the
  result gains a leading dimension indexing the draws.

- iterations:

  `numeric`; optional indices of the stored draws to use, between 1 and
  the number of draws the fit retains. Default is `NULL` to use all of
  them.

- offset:

  `numeric`; an offset for the rows of `newdata`, on the link scale,
  given either as one value per row or as a matrix with one column per
  additive predictor. Used only when `newdata` is supplied, and required
  there when the model was fit with an observation-level offset, since
  an offset is not a function of the predictors and so cannot be
  reconstructed. Default is `NULL` to add nothing.

- weights:

  `numeric`; prior weights for the rows of `newdata`, used only by
  `type = "density"` and ignored otherwise. For a binomial response
  given as proportions these are the numbers of trials, as in
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html). Default is `NULL`,
  which uses the weights the model was fit with when `newdata` is
  omitted and a weight of 1 for every row when it is supplied.

- values:

  `numeric`; for `type = "mean"`, what each response category is worth,
  given as a vector named for every level of the response. Default is
  `NULL` to read the level labels as numbers, which is an error rather
  than a guess when they cannot be read that way. Ignored with a warning
  for every other `type`.

- log:

  `logical`; for `type = "density"`, whether to return the log of the
  value. Default is `FALSE`. Summing the log density across observations
  gives a log score. Note that with `draws = FALSE` the density is
  averaged over the draws before the log is taken, so the result is the
  pointwise predictive density rather than the average log density.

- times:

  `numeric`; for `type = "survival"` and `type = "hazard"`, the times at
  which to report the survival function or the hazard, which must be
  finite and strictly positive. It has no default, because the horizon
  is a choice rather than a property of the fit. Ignored with a warning
  for every other `type`.

- ...:

  ignored; present for compatibility with the generic.

## Value

With `draws = FALSE`, a vector with one element per observation for
`"response"`, `"mean"` and `"density"`, and for `"link"` and `"stdlv"`
when the family has a single additive predictor; a matrix of
observations by additive predictors for `"link"` and `"stdlv"` when it
has more than one; a matrix of observations by categories for `"prob"`,
and for `"response"` with an ordinal or multinomial family; a matrix of
observations by times for `"survival"` and `"hazard"`; and a factor for
`"class"`, ordered when the family is
[`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md).

With `draws = TRUE` each of these gains a leading dimension indexing the
posterior draws, so that a vector becomes a matrix of draws by
observations and a matrix becomes an array of draws by observations by
categories or by times. On the `"link"` and `"stdlv"` scales a family
with more than one additive predictor returns a list with one
draws-by-observations matrix per predictor. `"class"` is a factor either
way.

## Details

### Simulated Multinomial Probit Probabilities

The likelihood of a `multinomial("probit")` fit has no closed form: the
probability of a category is the chance that the largest of several
correlated Gaussian variables is the one belonging to it, which is a
multivariate orthant probability. So `"prob"`, `"class"`, `"response"`
and `"density"` are all simulated, using the number of replicates the
family was given. Fresh draws are taken on every call, so two calls
differ by Monte Carlo error; that error is per posterior draw and
averages down over them, which makes `draws = FALSE` much more accurate
than any single row of `draws = TRUE`.

### `type = "stdlv"`

This reports `(eta - E[e]) / sd(y*)` for the latent `y* = eta + e`,
following
[`WeightIt::predict.ordinal_weightit()`](https://ngreifer.github.io/WeightIt/reference/predict.glm_weightit.html)
, and makes an ordinal or binary fit comparable across links.

The scale is `sd(y*) = sqrt(var(eta) + var(e))`, where `var(eta)` is
taken over the sample the model was fitted to, per draw, so it is a
property of the model rather than of whatever is being predicted and the
same divisor is used for new data. `var(e)` is whatever the link
implies: 1 for the probit link, \\\pi^2 / 3\\ for the logit, \\\pi^2 /
6\\ for the complementary log-log. The location subtracts the latent
error's mean, which is invisible for the logit and probit links, whose
errors are already centered, and is the whole of the difference for the
complementary log-log.

Such a model is identified only up to a common shift of its thresholds
and its predictor, so the location here is a convention, and the one
used is the same the cutpoints use: a predictor centered over the fitted
sample. Differences on this scale, which a standardized quantity exists
to support, are unaffected by that choice.

### `type = "density"`

The measure differs across the survival families: the accelerated
failure time families,
[`dpm_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
included, report the density of \\\log T\\, where
[`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
reports the density of \\T\\. The two differ by \\\sum \log t\\, so log
scores are comparable within each group rather than across them, and the
`type = "survival"` prediction is comparable throughout.
[`loo()`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
takes a `scale` argument that puts them on one measure.

Setting `type = "density"` returns `NaN` for an observation whose
density is undefined at any of the draws, which happens when a link the
package does not compile has been composed onto the family's own scale
and its inverse does not cover the whole additive predictor (e.g.,
`poisson("identity")`, which gives a positive mean only where the
predictor is positive). A warning reports how many draw-by-observation
values were undefined and how many returned values that made `NaN`. The
draws are averaged before the log is taken, so one undefined draw is
enough to make an observation `NaN`.

## See also

[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
for fitting the model;
[`bartisan-marginaleffects`](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
for averages and contrasts of these predictions;
[`bartisan-interop`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
for the methods that let other packages assess the fit

## Examples

``` r
data("rhc")
set.seed(123)

# Whether a patient died, with the treatment and every
# covariate recorded on admission as candidate
# predictors
model <- death ~ rhc + age + sex + race + edu + aps +
  meanbp + resp + hema + pafi + paco2 + crea + surv2m +
  card

fit <- bartisan(model, data = rhc, num_trees = 10,
                num_burn = 50, num_draws = 50)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` explicitly to silence this message.

# Fitted probabilities of death, averaged over the draws
head(predict(fit, type = "response"))
#> [1] 0.7959859 0.7976197 0.2008229 0.5216073 0.4393714 0.5232934

# The whole posterior for the first five patients rather
# than its mean
post <- predict(fit, newdata = rhc[1:5, ], draws = TRUE)
apply(post, 2, quantile, c(.025, .5, .975))
#>            [,1]      [,2]      [,3]      [,4]      [,5]
#> 2.5%  0.7378238 0.6893544 0.1062475 0.4105328 0.3025466
#> 50%   0.7921700 0.8089335 0.1979082 0.5251180 0.4428998
#> 97.5% 0.8489815 0.8587385 0.3142898 0.6063779 0.5359435

# A held-out log score, which needs the outcome, so the
# `newdata` argument carries it
dens <- predict(fit, newdata = rhc[1:100, ],
                type = "density")
sum(log(dens))
#> [1] -58.29588
```
