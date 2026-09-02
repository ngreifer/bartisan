# Fit a generalized Bayesian additive regression trees (BART) model

Fits a BART model in which the response distribution is arbitrary rather
than restricted to the conditionally conjugate cases, using the
Laplace-approximation reversible-jump sampler of Linero (2025). Decision
rules may be soft, as in Linero and Yang (2018), which gives smoother
fits than the step functions of standard BART.

The interface deliberately mirrors
[`stats::glm()`](https://rdrr.io/r/stats/glm.html): a formula, a data
frame and a family. Ordinary
[stats::family](https://rdrr.io/r/stats/family.html) objects work
unchanged, including their links, and the extra families that
[`glm()`](https://rdrr.io/r/stats/glm.html) has no counterpart for are
documented at
[bartisan-families](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
along with
[`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
for a likelihood of your own.

## Usage

``` r
bartisan(
  formula,
  data,
  family = NULL,
  weights = NULL,
  offset = NULL,
  subset = NULL,
  na.action = stats::na.pass,
  control = bartisan_control(),
  ...
)
```

## Arguments

- formula:

  a model formula. The right-hand side lists candidate predictors; the
  model finds interactions and nonlinearity on its own, so `y ~ .` is
  usually the right specification. Survival families take a
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)
  object on the left. A `(1 | group)` term adds a group-level random
  intercept, in the notation of lme4; see Details.

  For a family with more than one additive predictor this may be a
  *list* of formulas, one per forest, to give each one its own
  predictors. The first is the model for the main parameter and carries
  the response; the rest need no response, and follow the order in
  [bartisan-families](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
  under "Several additive predictors", which also gives the name of each
  forest so the list can be named instead of ordered:

      bartisan(list(y ~ x1 + x2, ~ x2 + x3), data = d, family = location_scale())
      bartisan(list(mean = y ~ x1 + x2, log_sd = ~ x2), data = d,
               family = location_scale())

  One formula applies to every forest, which is the ordinary case. A
  predictor left out of one forest's formula is still in the data and is
  never split on by that forest.

  [`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) terms
  are read out of each formula in turn, so a parameter has the varying
  coefficients its own formula asks for and no others. That makes the
  forests two-dimensional – one axis the parameter, the other the
  coefficient – and the names below are what per-forest settings are
  keyed by:

      # forests: mean, mean:z, log_sd
      bartisan(list(mean = y ~ x1 + x2 + vc(z), log_sd = ~ x1 + x2), data = d,
               family = location_scale())

      # one formula reaches every parameter, so both get a coefficient of `z`
      bartisan(y ~ x1 + x2 + vc(z), data = d, family = location_scale())

- data:

  a data frame containing the variables in `formula`.

- family:

  the response distribution, as a
  [stats::family](https://rdrr.io/r/stats/family.html) object, one of
  the families in
  [bartisan-families](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
  or a name. The default, `NULL`, reads one off the response and says
  which it chose; see Details for the rules and for what is supported.

- weights:

  optional prior weights. For a binomial response given as proportions,
  these are the numbers of trials, as in
  [`glm()`](https://rdrr.io/r/stats/glm.html).

- offset:

  optional known component of the additive predictor, on the link scale.

- subset:

  optional vector specifying a subset of rows to use.

- na.action:

  how to handle missing values. The default,
  [stats::na.pass](https://rdrr.io/r/stats/na.fail.html), keeps rows
  whose *predictors* are missing and lets the splitting rules decide
  where they go, which is what the trees are able to do and
  [`lm()`](https://rdrr.io/r/stats/lm.html) and
  [`glm()`](https://rdrr.io/r/stats/glm.html) are not; see Details. Pass
  [stats::na.omit](https://rdrr.io/r/stats/na.fail.html) to drop any row
  with a missing value anywhere instead. Rows with a missing response,
  weight or offset are dropped either way, with a warning, since there
  is nothing to fit them to.

- control:

  a list of sampler and prior settings from
  [`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md).

- ...:

  further arguments to
  [`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md).
  They are merged into `control`, overriding any value given there, so
  that `bartisan(..., num_trees = 20)` and
  `bartisan(..., control = bartisan_control(num_trees = 20))` are the
  same call. Names that are not arguments of
  [`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
  are an error rather than being silently ignored.

## Value

An object of class `bartisan`, a list with elements including:

- `eta`:

  a list with one matrix per additive predictor, each of posterior draws
  by observation, on the link scale.

- `fitted`:

  fitted values on the response scale, averaged over draws.

- `counts`:

  a list with one matrix per additive predictor, of the number of
  splitting rules using each predictor group in each draw. Useful for
  variable selection.

- `aux`:

  draws of the nuisance parameters, such as the residual standard
  deviation or the ordinal cutpoints, when the family has any.

- `has_na`:

  which predictor columns contained a missing value, which is what
  determines where [`predict()`](https://rdrr.io/r/stats/predict.html)
  will accept one.

- `rhat`:

  a data frame of convergence diagnostics, when more than one chain was
  run: rank-normalized folded split R-hat and the bulk and tail
  effective sample sizes (Vehtari et al. 2021) for the log likelihood,
  the leaf scales, the nuisance parameters and the additive predictor.
  R-hat above about 1.01 says the chains have not agreed; an effective
  sample size below about 400 says the run is too short for the quantity
  it belongs to, and the tail column is the one that governs interval
  endpoints.

- `sigma_mu`, `bandwidth`:

  draws of the leaf standard deviation and, for soft rules, the per-tree
  gate bandwidths.

- `loglik`:

  the log likelihood at each draw.

## What the sampler does

Standard BART relies on the leaf parameters being integrable in closed
form, which restricts it to a Gaussian response, or to models that can
be reduced to one by data augmentation. Linero's algorithm removes that
restriction. At each candidate move it builds a Gaussian approximation
to the conditional posterior of the affected leaf parameters, by Fisher
scoring, and uses that approximation as the proposal in a
reversible-jump Metropolis step. The approximation only has to be good
enough to be accepted often; the stationary distribution is the exact
posterior either way.

What a new family therefore has to supply is only the log density of one
observation and its first two derivatives with respect to the additive
predictor. Families whose response has more than one unconstrained
parameter, such as
[`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
and
[`location_scale()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
carry one forest per parameter. Because that is the whole interface, it
can be reached from R:
[`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
takes the log density as an R function and differences it for the
derivatives, and a link the package does not compile is composed onto
the scale its family works on the same way.

## The family is inferred when you do not name one

`family` may be left alone, in which case it is read off the response:

|  |  |
|----|----|
| Response | Family |
| [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) object, or a two-column matrix of times and events | [`dpm_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) |
| ordered factor | [`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) |
| logical, or two levels, or numeric zeros and ones | [`binomial()`](https://rdrr.io/r/stats/family.html) |
| factor or character with more than two levels | [`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) |
| two-column matrix of successes and failures | [`binomial()`](https://rdrr.io/r/stats/family.html) |
| anything else | [`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) |

A message reports the choice. Naming `family` yourself is what silences
it, which is the same thing you would do to change the choice.

Two of these are worth saying out loud. A **count** is not inferred as
[`poisson()`](https://rdrr.io/r/stats/family.html): a non-negative
integer response is often Poisson and often not, and the Poisson
variance assumption is strong enough that making it silently would be a
modeling decision taken on the caller's behalf. Gaussian is the weaker
guess and the one whose failure is easy to see. And a numeric response
with exactly two values that are *not* zero and one – `c(1, 2)`, say –
is Gaussian rather than binomial, because which of the two counts as the
success is not something to guess at.

## Soft decision rules

By default a decision rule is a smooth gate rather than a step, so an
observation reaches every leaf with some weight and the fitted function
is smooth. `gate` in
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
chooses both whether the rules are soft and, if they are, the gate's
shape; the default is the bounded `"smoothstep"`, and `"logistic"` is
Linero and Yang's (2018) original. This costs more per iteration, since
a leaf now touches every observation rather than the ones inside its
cell, and it makes the leaf parameters of a tree dependent on one
another. Combining soft rules with a non-conjugate likelihood is an
extension of Linero (2025), which leaves it as an open problem; it is
handled here by giving the reversible-jump move a bivariate Laplace
proposal for the pair of child leaves, which reduces to Linero's
independent pair exactly when the rules are hard.

Set `gate = "hard"` in
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
for the faster hard-rule sampler.

## Random intercepts

A `(1 | group)` term in the formula adds an intercept per level of
`group`, drawn from a common mean-zero normal whose standard deviation
is itself drawn under the same half-Cauchy prior the leaf scale uses.
Several grouping factors are allowed, and `(1 | a/b)` expands to nesting
as it does in lme4:

    bartisan(y ~ x1 + x2 + (1 | school), data = d)
    bartisan(y ~ x1 + (1 | school) + (1 | year), data = d)

The intercepts are in `fit$ranef` and their standard deviations in
`fit$tau`, one matrix per additive predictor. A family with several
predictors gets a separate set for each – a zero-inflated count model
has a group effect on the count part and another on the inflation part –
and they are independent of one another.

Only random *intercepts* are supported, and a random slope is refused
rather than ignored. The reason is that a random intercept is a scalar
entering the predictor with weight one for the observations in its
level, which is what a leaf is once its gate is removed, so the
sampler's leaf machinery handles it exactly; a slope is a different
shape of parameter. A variable whose effect varies by group belongs in
the fixed part of the formula, where a tree can split on the group and
on the variable together and get an interaction of any shape.

**When to reach for this rather than putting the group in as a
predictor.** A grouping factor can also go in the fixed part, where a
tree splits on it like anything else, and with few large groups that is
the better choice – measured, it beats a random intercept, because the
group means are well determined without pooling and a split can interact
the group with the covariates. The random intercept wins when there are
many small groups, which is where partial pooling earns its keep: at 250
groups of four observations it cut held-out error by 30% against the
factor route, and at five groups of a hundred it lost to it.

A level of `group` that was not present at fitting time is given the
prior mean of zero when predicting, with a warning.

## Missing predictor values

A missing predictor is not imputed and its row is not dropped, which is
the default here because a tree can do something better with a missing
value than either. Instead each splitting rule carries the answer for
itself. A rule on a variable that has missing values is drawn as one of
three, with equal probability:

- `x < c`, or missing, goes left;

- `x < c` goes left, missing goes right;

- missing goes left, present goes right.

This is *missingness incorporated in attributes* (Twala, Jones and Hand
2008; for BART, Kapelner and Bleich 2015). The third rule is what lets
the model split on missingness itself, so a variable whose *absence*
carries the signal is usable even if its observed values say nothing.
Since the choice is drawn from its prior along with the variable and the
cutpoint, it cancels from every acceptance ratio, and a variable with no
missing values is not given the extra draw at all: complete data
reproduces the sampler exactly as it was.

A missing value takes a hard path through the tree even when the rules
are soft, which is the right thing – there is nothing about being absent
to smooth over – and it keeps the leaf weights summing to one.

Two consequences to be clear about.
[`predict()`](https://rdrr.io/r/stats/predict.html) accepts missing
values only in columns that had them at fitting time, because only those
columns' rules carry an answer; elsewhere every rule would send the
value the same arbitrary way, so it is an error instead. And what this
estimates is the mean of the response given the predictors *and the
pattern of missingness*. That is what you want for prediction. If the
estimand is a regression or causal effect defined on complete data,
multiple imputation is the right tool and this is not.

## Preprocessing

Predictors are mapped to the unit interval, because the cutpoint prior
is uniform on a node's live range and the soft-rule bandwidth is
measured on the predictor scale. Factors are expanded to an indicator
per level and share a single weight in the sparsity prior, so that a
factor is selected or not as a whole rather than one level at a time.
The additive predictor starts from an intercept-only fit, so the leaf
prior describes departures from that fit rather than the absolute level
of the response. That starting value is the exact null-model estimate
for most families; for the accelerated failure time families, where
censoring makes the sample mean of the log times biased, and for the
zero-inflated and ordered beta families, it is a moment approximation,
which the sampler then moves away from.

## References

Linero, A. R. (2025). Generalized Bayesian additive regression trees
models: beyond conditional conjugacy. *Journal of the American
Statistical Association*, 120(549), 356–369.
[doi:10.1080/01621459.2024.2337156](https://doi.org/10.1080/01621459.2024.2337156)

Linero, A. R., & Yang, Y. (2018). Bayesian regression tree ensembles
that adapt to smoothness and sparsity. *Journal of the Royal Statistical
Society Series B*, 80(5), 1087–1110.
[doi:10.1111/rssb.12293](https://doi.org/10.1111/rssb.12293)

Albert, J. H., & Chib, S. (1993). Bayesian analysis of binary and
polychotomous response data. *Journal of the American Statistical
Association*, 88(422), 669–679.
[doi:10.1080/01621459.1993.10476321](https://doi.org/10.1080/01621459.1993.10476321)

Polson, N. G., Scott, J. G., & Windle, J. (2013). Bayesian inference for
logistic models using Polya-Gamma latent variables. *Journal of the
American Statistical Association*, 108(504), 1339–1349.
[doi:10.1080/01621459.2013.829001](https://doi.org/10.1080/01621459.2013.829001)

Kapelner, A., & Bleich, J. (2015). Prediction with missing data via
Bayesian additive regression trees. *Canadian Journal of Statistics*,
43(2), 224–239.
[doi:10.1002/cjs.11248](https://doi.org/10.1002/cjs.11248)

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Buerkner, P.-C.
(2021). Rank-normalization, folding, and localization: an improved
\\\widehat{R}\\ for assessing convergence of MCMC. *Bayesian Analysis*,
16(2), 667–718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

Twala, B. E. T. H., Jones, M. C., & Hand, D. J. (2008). Good methods for
coping with missing data in decision trees. *Pattern Recognition
Letters*, 29(7), 950–956.
[doi:10.1016/j.patrec.2008.01.010](https://doi.org/10.1016/j.patrec.2008.01.010)

## See also

[`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md),
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md),
[bartisan-families](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
and
[`vignette("families", package = "bartisan")`](https://ngreifer.github.io/bartisan/articles/families.md)
for a family-by-family guide.

## Examples

``` r
set.seed(1)

n <- 200
d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
d$y <- rbinom(n, 1, plogis(3 * sin(pi * d$x1 * d$x2) - 1))

fit <- bartisan(y ~ x1 + x2 + x3, data = d, family = binomial(),
               control = bartisan_control(num_trees = 10, num_burn = 50,
                                         num_draws = 50, verbose = FALSE))
fit
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = y ~ x1 + x2 + x3, data = d, family = binomial(), 
#>     control = bartisan_control(num_trees = 10, num_burn = 50, 
#>         num_draws = 50, verbose = FALSE))
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 200
#> Structure: 1 forest of 10 trees, soft decision rules
#> Draws: 50 kept after 50 warmup

head(predict(fit, type = "response"))
#> [1] 0.3255948 0.4875942 0.8162769 0.6190264 0.1925652 0.7843789
```
