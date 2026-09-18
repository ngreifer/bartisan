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
  predictors, `function(y, eta)`, where `y` is a numeric vector of
  length `n` and `eta` an `n` by `num_predictors` matrix, returning a
  numeric vector of length `n`. With nuisance parameters it takes a
  third argument, `function(y, eta, aux)`, where `aux` is a numeric
  vector of their current values. It is the log density of one unit of
  prior weight, so that `weights` behave as they do elsewhere, and terms
  free of `eta` may be dropped.

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
  length `n`. Default is `NULL` to take central differences of
  `logdens`. It covers the additive predictors only: a nuisance
  parameter is always differenced, which costs three calls per sweep
  rather than three per leaf.

- aux_names:

  optional `character`; the names of the nuisance parameters to draw, if
  any. Naming them is what declares them, because the names label the
  columns of `fit$aux` and are what
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
  second way to declare the parameters, so `aux_start = c(shape = 1)`
  both names one and starts it at 1.

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

Because it sees a subset rather than the whole sample, anything else the
density needs has to be a scalar or reach it through `eta`. A
per-observation vector captured from the enclosing environment will not
line up with the rows it is handed, and nothing can detect that for you:
the fit runs and is wrong. Where such a quantity is genuinely needed, an
offset belongs in the formula and a varying trial count or exposure
belongs in a family written for it.

### Nuisance Parameters

These are drawn alongside the trees when `aux_names` names them, and
`logdens` then takes a third argument holding their current values:

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
the [`exp()`](https://rdrr.io/r/base/Log.html) above is what keeps the
scale positive.

### What a Log Density Cannot Supply

A density says how likely an observed value is, not how to draw a new
one, so a `custom_family()` fit has no posterior predictive
distribution, which is what
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
already computes.

Two smaller limits. The response must be numeric, so a factor has to be
coded first. And since the package cannot know what the mean of the
density is, `predict(type = "response")` returns the additive predictors
rather than a fitted mean.

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

# A Poisson likelihood written out by hand, with its derivatives. Terms
# free of eta may be dropped, since they cancel from every acceptance
# ratio, and the same terms are absent from the score.
pois <- custom_family(
  function(y, eta) y * eta[, 1] - exp(eta[, 1]),
  derivatives = function(y, eta, h) {
    list(score = y - exp(eta[, 1]), info = exp(eta[, 1]))
  },
  start = log(mean(d$y)))

fit <- bartisan(y ~ x1 + x2, data = d, family = pois,
                num_trees = 20, num_burn = 100, num_draws = 100,
                verbose = FALSE)

fit
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = y ~ x1 + x2, data = d, family = pois, num_trees = 20, 
#>     num_burn = 100, num_draws = 100, verbose = FALSE)
#> 
#> Family: "custom" (supplied from R)
#> Observations: 300
#> Structure: 1 forest of 20 trees, soft decision rules
#> Draws: 100 kept after 100 warmup

# A beta-binomial, which no built-in family covers: counts out of a known
# number of trials, overdispersed relative to a binomial. `phi` is drawn
# alongside the trees, and `size` is a scalar, so closing over it is safe.
size <- 30
d$hits <- rbinom(300, size, rbeta(300, plogis(d$x1) * 6,
                                  (1 - plogis(d$x1)) * 6))

bb <- custom_family(
  logdens = function(y, eta, aux) {
    p <- plogis(eta[, 1])
    phi <- exp(aux[1])
    lbeta(y + p * phi, size - y + (1 - p) * phi) -
      lbeta(p * phi, (1 - p) * phi)
  },
  aux_names = "log_phi", aux_start = log(5), name = "beta-binomial")

fit_bb <- bartisan(hits ~ x1 + x2, data = d, family = bb,
                   num_trees = 20, num_burn = 100, num_draws = 100,
                   verbose = FALSE)

# The drawn precision, on the scale it was written on.
exp(mean(fit_bb$aux[, "log_phi"]))
#> [1] 5.889742
```
