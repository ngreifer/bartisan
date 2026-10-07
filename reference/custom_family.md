# Fit a model to a likelihood written in R

`custom_family()` takes the log density itself, as an R function, and
returns a family object
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
fits the model for. It is the route to a response distribution none of
the compiled families in
[`bartisan-families`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
cover.

## Usage

``` r
custom_family(
  logdens,
  num_predictors = 1L,
  start = 0,
  derivatives = NULL,
  aux_names = NULL,
  aux_start = NULL,
  name = "custom"
)
```

## Arguments

- logdens:

  the log density, given as a function of the response and the additive
  predictors, `function(y, eta)`, where `y` holds the response of `n`
  observations and `eta` is an `n` by `num_predictors` matrix of their
  additive predictors, returning a numeric vector of length `n`. `y` is
  the response as it appears in the model formula, with only those
  observations kept: a vector or factor of length `n`, or a matrix with
  `n` rows, such as one made with
  [`cbind()`](https://rdrr.io/r/base/cbind.html) or
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) . See
  Details. With nuisance parameters it must take a third argument,
  `function(y, eta, aux)`, where `aux` is a numeric vector of their
  current values. It is the log density of one unit of prior weight, so
  that `weights` behave as they do elsewhere, and terms free of `eta`
  may be dropped.

- num_predictors:

  `numeric`; how many additive predictors the density has (i.e., how
  many forests to fit). Default is 1.

- start:

  `numeric`; the value each additive predictor starts at, in place of
  the intercept-only fit the compiled families use. One value or one per
  predictor. Default is 0.

- derivatives:

  optional; a `function(y, eta, h)` returning a list with elements
  `score` and `info`, the first derivative of `logdens` with respect to
  the `h`th predictor and minus its second derivative, each a vector of
  length `n`. With nuisance parameters it takes a fourth argument,
  `function(y, eta, h, aux)`, which receives their current values as
  `logdens` does. Default is `NULL` to take central differences of
  `logdens`. It covers the additive predictors only: a nuisance
  parameter is always differenced, which costs three calls per sweep
  rather than three per leaf.

- aux_names:

  optional `character`; the names of the nuisance parameters to draw, if
  any. Naming them declares them, because the names label the columns of
  `fit$aux` and are the labels
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
  report them under. They must be distinct and non-empty. Default is
  `NULL` for none, unless `aux_start` is given, in which case the
  parameters are named by `names(aux_start)` when it carries names and
  positionally (`"aux1"`, `"aux2"`, and so on) when it does not.

- aux_start:

  optional `numeric`; the value each nuisance parameter starts at, given
  as one value or one per parameter. Default is `NULL`, which is 0 for
  each. The sampler will walk to the posterior from a poor start, so
  this need only be the right order of magnitude. Supplying it is a
  second way to declare the parameters, so setting
  `aux_start = c(shape = 1)` both names one and starts it at 1.

- name:

  string; a label used when printing the fit. Default is `"custom"`.

## Value

A `<bartisan_family>` object, which is a list containing at least the
elements `family` and `link` and which inherits from `family`, so that
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
recognizes it wherever it recognizes an ordinary
[`stats::family`](https://rdrr.io/r/stats/family.html) object.

## Details

Nothing else about the sampler changes: the leaf-level Laplace proposal
needs the first two derivatives of the log density with respect to each
additive predictor and nothing more, and central differences of the
supplied function produce both.

    # A Poisson model written out by hand. Terms free of eta may be dropped;
    # they cancel from every acceptance ratio.
    pois <- custom_family(function(y, eta) y * eta[, 1] - exp(eta[, 1]),
                          start = log(mean(d$y)))

    # Two predictors: a mean and a log standard deviation.
    ls <- custom_family(function(y, eta) dnorm(y, eta[, 1], exp(eta[, 2]),
                                               log = TRUE),
                        num_predictors = 2, start = c(0, 0))

The function is called once per leaf per Fisher-scoring step with the
observations reaching that leaf, so it must be vectorized over `y` and
the rows of `eta`; it must not be vectorized *within* an observation,
and it must return exactly one value per row. Supplying `derivatives`
cuts three calls to one and removes the differencing error.

Because it sees a subset rather than the whole sample, a per-observation
vector captured from the enclosing environment will not line up with the
rows it is handed, and nothing can detect that for you: the fit runs and
is wrong. Anything the density needs that varies across observations
should instead be part of the response, which is subset along with
`eta`. A varying number of trials, for example, can be supplied as a
second column, as in `cbind(hits, size) ~ x`, and read inside `logdens`
as `y[, 2]`. An offset belongs in the formula, as it does for the other
families.

### Nuisance Parameters

These are drawn alongside the trees when the `aux_names` argument names
them, and the function supplied as `logdens` then takes a third argument
holding their current values:

    # A Gaussian written out by hand, with its scale drawn rather than fixed.
    by_hand <- custom_family(
      logdens = function(y, eta, aux) dnorm(y, eta[, 1], exp(aux[1]), log = TRUE),
      aux_names = "log_sigma", aux_start = 0)

They are reported in `fit$aux` under those names, and covered by
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
like any other family's. There is no prior argument and no bounds
argument, because a nuisance parameter here is carried as an additive
predictor whose forest is pinned at depth zero (one tree that can never
split, so the forest is a single scalar), and it is drawn by the same
Laplace-plus-Metropolis step as any leaf, under that step's Gaussian
leaf prior. So a parameter with a restricted range is handled the way it
would be for a real predictor, by writing the transform into `logdens`:
the [`exp()`](https://rdrr.io/r/base/Log.html) above keeps the scale
positive.

### Responses Other Than a Numeric Vector

The response reaches `logdens` as it appears in the model formula, so it
can be a factor, a matrix with a column per outcome, or a
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) object,
which keeps its `time` and `status` columns. A factor arrives with the
levels it had when the model was fit, in
[`predict()`](https://rdrr.io/r/stats/predict.html) with `newdata` as
well, so a density can read `as.integer(y)` as a category code. A matrix
response makes it possible to fit a multivariate BART model, with one
additive predictor per outcome and the outcomes related through the
error distribution:

    # A bivariate normal, with the two log standard deviations and the Fisher
    # transform of the correlation drawn as nuisance parameters.
    bvn <- custom_family(
      logdens = function(y, eta, aux) {
        s <- exp(aux[1:2])
        rho <- tanh(aux[3])
        z1 <- (y[, 1] - eta[, 1]) / s[1]
        z2 <- (y[, 2] - eta[, 2]) / s[2]
        -sum(log(s)) - log1p(-rho^2) / 2 -
          (z1^2 - 2 * rho * z1 * z2 + z2^2) / (2 * (1 - rho^2))
      },
      num_predictors = 2, start = c(mean(d$y1), mean(d$y2)),
      aux_names = c("log_sd1", "log_sd2", "atanh_rho"))

    fit <- bartisan(cbind(y1, y2) ~ x1 + x2, data = d, family = bvn)

Each forest's prior is scaled for an additive predictor that varies by a
few units around its `start`, so outcomes measured on a large scale are
best standardized first. A row with any missing value in the response is
dropped, as it is for the other families.

### Limits of a Log Density

A density corresponds to how likely an observed value is, not how to
draw a new one, so a `custom_family()` fit has no posterior predictive
distribution, which
[`simulate()`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md),
[`pp_check()`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
and
[`r2()`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
rely on. Each of them errors on such a fit rather than returning
something it cannot support.
[loo()](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
and
[`waic()`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
are unaffected, since both read the pointwise log likelihood the family
already computes. Note that terms dropped from `logdens` because they
are free of `eta` are missing from that log likelihood too, so its
values can be compared across fits only when the same terms are dropped
from each (e.g., not against a compiled family, whose density is
complete).

There is a smaller limit as well. Since the package cannot know what the
mean of the density is, `predict(type = "response")` returns the
additive predictors rather than a fitted mean.

## See also

[`bartisan-families`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
for the compiled families, one of which is usually the better answer;
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
for fitting a model with the result;
[`vignette("families", package = "bartisan")`](https://ngreifer.github.io/bartisan/articles/families.md)
for the long form

## Examples

``` r
set.seed(123)

d <- data.frame(x1 = runif(300), x2 = runif(300))
d$y <- rpois(300, exp(1 + sin(pi * d$x1)))

# A Poisson likelihood written out by hand, with its
# derivatives. Terms free of eta may be dropped, since
# they cancel from every acceptance ratio, and the same
# terms are absent from the score.
pois <- custom_family(
  function(y, eta) y * eta[, 1] - exp(eta[, 1]),
  derivatives = function(y, eta, h) {
    list(score = y - exp(eta[, 1]),
         info = exp(eta[, 1]))
  },
  start = log(mean(d$y)))

fit <- bartisan(y ~ x1 + x2,
                data = d, family = pois,
                num_trees = 10, num_burn = 50,
                num_draws = 50)

fit
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = y ~ x1 + x2, data = d, family = pois, num_trees = 10, 
#>     num_burn = 50, num_draws = 50)
#> 
#> Family: "custom" (supplied from R)
#> Observations: 300
#> Structure: 1 forest of 10 trees, soft decision rules
#> Draws: 50 kept after 50 warmup

# \donttest{
# A beta-binomial, which no built-in family covers:
# counts out of a known number of trials, overdispersed
# relative to a binomial. The number of trials varies,
# so it is the response's second column, which reaches
# `logdens` row for row with `eta`. `phi` is drawn
# alongside the trees.
d$size <- sample(10:40, 300, replace = TRUE)
p_true <- plogis(d$x1)
d$hits <- rbinom(300, d$size,
                 rbeta(300, p_true * 6, (1 - p_true) * 6))

bb <- custom_family(
  logdens = function(y, eta, aux) {
    hits <- y[, 1]
    size <- y[, 2]
    p <- plogis(eta[, 1])
    phi <- exp(aux[1])
    lbeta(hits + p * phi, size - hits + (1 - p) * phi) -
      lbeta(p * phi, (1 - p) * phi)
  },
  aux_names = "log_phi", aux_start = log(5),
  name = "beta-binomial")

fit_bb <- bartisan(cbind(hits, size) ~ x1 + x2,
                   data = d, family = bb,
                   num_trees = 10, num_burn = 50,
                   num_draws = 50)

# The drawn precision, on the scale it was written on.
exp(mean(fit_bb$aux[, "log_phi"]))
#> [1] 5.522993
# }
```
