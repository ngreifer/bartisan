# Fit a generalized Bayesian additive regression trees (BART) model

Fits a BART model in which the response distribution is arbitrary rather
than restricted to the conditionally conjugate cases, using the
Laplace-approximation reversible-jump sampler of Linero (2025). Decision
rules may be soft, as in Linero and Yang (2018), which gives smoother
fits than the step functions of standard BART. The interface mirrors
that of [`stats::glm()`](https://rdrr.io/r/stats/glm.html): a formula, a
data frame, and a family, with the families that
[`glm()`](https://rdrr.io/r/stats/glm.html) has no counterpart for
documented at
[`bartisan-families`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md).

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
  prior_only = FALSE,
  ...
)
```

## Arguments

- formula:

  a model formula. The right-hand side lists candidate predictors; the
  model finds interactions and nonlinearity on its own, so
  `y ~ x1 + x2 + x3` is usually the right specification. Survival
  families take a
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)
  object on the left. A `(1 | group)` term adds a group-level random
  intercept, in the notation of lme4; see Details.

  For a family with more than one additive predictor this may be a
  *list* of formulas, one per forest, to give each one its own
  predictors. The first is the model for the main parameter and carries
  the response; the rest need no response, and follow the order in
  [`bartisan-families`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
  under "Several additive predictors", which also gives the name of each
  forest so the list can be named instead of ordered:

      bartisan(list(y ~ x1 + x2, ~ x2 + x3), data = d,
               family = gaussian_ls())
      bartisan(list(mean = y ~ x1 + x2, log_sd = ~ x2), data = d,
               family = gaussian_ls())

  One formula applies to every forest, which is the ordinary case. A
  predictor left out of one forest's formula is still in the data and is
  never split on by that forest.

  A formula naming no predictor at all makes that parameter a constant.
  `~ 1` leaves its forest nothing to split on, so every tree in it is a
  stump and the forest is a single drawn scalar rather than a function
  of the predictors.

  [`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) terms
  are read out of each formula in turn, so a parameter has the varying
  coefficients its own formula asks for and no others, which makes the
  forests two-dimensional (one axis the parameter, the other the
  coefficient).
  [`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md)
  documents how they are then named and keyed.

- data:

  a data frame containing the variables named in `formula`.

- family:

  the response distribution, given as a
  [`stats::family`](https://rdrr.io/r/stats/family.html) object, as one
  of the families in
  [`bartisan-families`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
  or as the name of either. A `family` object is accepted when the
  distribution it names is one this package implements, since the
  likelihood is the package's rather than the object's: a `family`
  object carries a link and a variance function and not a density, so
  one naming anything else (e.g.,
  [`stats::inverse.gaussian`](https://rdrr.io/r/stats/family.html), or a
  Tweedie from another package) is an error rather than something a
  likelihood can be built from, and
  [`custom_family()`](https://ngreifer.github.io/bartisan/reference/custom_family.md)
  is the route for those. Links are used as supplied, and a link the
  package does not compile is composed onto the scale its family works
  on. Default is `NULL`, in which case the family is read off the
  response and a message reports the choice; see Details for the rules.

- weights:

  optional; prior weights, one per observation. For a binomial response
  given as proportions, these are the numbers of trials, as in
  [`glm()`](https://rdrr.io/r/stats/glm.html).

- offset:

  optional; a known component of the additive predictor, on the link
  scale. One value per observation, or a matrix with one column per
  additive predictor to give each its own; see Details.

- subset:

  optional; a vector specifying the subset of rows to use.

- na.action:

  how missing values are handled. Default is
  [`stats::na.pass`](https://rdrr.io/r/stats/na.fail.html), which keeps
  rows whose *predictors* are missing and lets the splitting rules
  decide where they go, which is something the trees can do and
  [`lm()`](https://rdrr.io/r/stats/lm.html) and
  [`glm()`](https://rdrr.io/r/stats/glm.html) cannot; see Details. Pass
  [`stats::na.omit`](https://rdrr.io/r/stats/na.fail.html) to drop any
  row with a missing value anywhere instead. Note that rows with a
  missing response, weight, or offset are dropped either way, with a
  warning, since there is nothing to fit them to.

- control:

  a `<bartisan_control>` object; the output of a call to
  [`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md),
  containing the sampler and prior settings.

- prior_only:

  `logical`; whether to draw from the prior rather than the posterior,
  which is what a prior predictive check reads. Default is `FALSE`. Not
  available for every family; see Details.

- ...:

  further arguments to
  [`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md),
  which are merged into `control` and override any value given there, so
  that `bartisan(..., num_trees = 20)` and
  `bartisan(..., control = bartisan_control(num_trees = 20))` are the
  same call. A name that is not an argument of
  [`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
  is an error rather than being silently ignored.

## Value

A `<bartisan_fit>` object, a list with the following components among
others.

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

- `sigma_mu`, `bandwidth`:

  draws of the leaf standard deviation and, for soft rules, the per-tree
  gate bandwidths.

- `loglik`:

  the log likelihood at each draw.

- `control`:

  the `<bartisan_control>` object the fit used, with any settings given
  in `...` merged in.

## Details

### The Sampler

Standard BART relies on the leaf parameters being integrable in closed
form, which restricts it to a Gaussian response, or to models that can
be reduced to one by data augmentation. Linero's (2025) algorithm
removes that restriction. At each candidate move it builds a Gaussian
approximation to the conditional posterior of the affected leaf
parameters, by Fisher scoring, and uses that approximation as the
proposal in a reversible-jump Metropolis step. The approximation only
has to be good enough to be accepted often; the stationary distribution
is the exact posterior either way.

What a new family therefore has to supply is only the log density of one
observation and its first two derivatives with respect to the additive
predictor. Families whose response has more than one unconstrained
parameter, such as
[`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
and
[`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
carry one forest per parameter. Because that is the whole interface, it
can be reached from R:
[`custom_family()`](https://ngreifer.github.io/bartisan/reference/custom_family.md)
takes the log density as an R function and differences it for the
derivatives.

### Inferring the Family

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

A message reports the choice, and naming `family` is what silences it,
which is also what changes it.

Two scenarios are worth noting. A count is read as a numeric variable
and therefore has
[`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
as its default. And a numeric response with exactly two values other
than zero and one (e.g., `c(1, 2)`) is also given
[`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
as its default rather than
[`binomial()`](https://rdrr.io/r/stats/family.html), which of the two
counts as the success not being something to guess at.

### Soft Decision Rules

By default a decision rule is a smooth gate rather than a step, so an
observation reaches every leaf with some weight and the fitted function
is smooth. `gate` in
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
chooses both whether the rules are soft and, if they are, the gate's
shape; the default is the bounded `"smoothstep"`, and `"logistic"` is
Linero and Yang's (2018) original. Soft rules cost more per iteration,
since a leaf now touches every observation rather than only the ones
inside its cell, and they make the leaf parameters of a tree dependent
on one another. Combining them with a non-conjugate likelihood is an
extension of Linero (2025), which leaves it as an open problem; it is
handled here by giving the reversible-jump move a bivariate Laplace
proposal for the pair of child leaves, which reduces to Linero's
independent pair exactly when the rules are hard. Set `gate = "hard"` in
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
for the faster hard-rule sampler.

### Random Intercepts

A `(1 | group)` term in the formula adds an intercept per level of
`group`, drawn from a common mean-zero normal whose standard deviation
is itself drawn under the same half-Cauchy prior the leaf scale uses.
Several grouping factors are allowed, and `(1 | a/b)` expands to nesting
as it does in lme4:

    bartisan(y ~ x1 + x2 + (1 | school), data = d)
    bartisan(y ~ x1 + (1 | school) + (1 | year), data = d)

The intercepts are in `fit$ranef` and their standard deviations in
`fit$tau`, one matrix per additive predictor. A family with several
predictors gets a separate set for each (i.e., a zero-inflated count
model has a group effect on the count part and another on the inflation
part), and they are independent of one another.

Only random *intercepts* are supported, and a random slope is refused
rather than ignored. The reason is that a random intercept is a scalar
entering the predictor with weight one for the observations in its
level, which is what a leaf is once its gate is removed, so the
sampler's leaf machinery handles it exactly; a slope is a different
shape of parameter. A variable whose effect varies by group belongs in
the fixed part of the formula, where a tree can split on the group and
on the variable together and get an interaction of any shape.

A grouping factor can also go in the fixed part, where a tree splits on
it like anything else, and with few large groups that is the better
choice: the group means are well determined without pooling and a split
can interact the group with the covariates. The random intercept wins
where there are many small groups, which is where partial pooling earns
its keep.

A level of `group` that was not present at fitting time is given the
prior mean of zero when predicting, with a warning.

### Offsets

An offset is a known part of the additive predictor, supplied on the
link scale and not estimated: the log of an exposure for a count, say.
It is given as one value per observation, through `offset` or as an
[`offset()`](https://rdrr.io/r/stats/offset.html) term in the formula.

A family with several additive predictors takes the same offset on each
of them, unless it is given a matrix with one column per predictor,
which offsets each separately: a category-specific exposure for
[`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
say, or an offset on the count part of
[`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
and not on its zero-inflation part. `?bartisan-families` lists each
family's forests in order, which is the column order.
[`predict()`](https://rdrr.io/r/stats/predict.html) takes the same two
forms.

A vector is worth thinking about before reaching for under
[`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
because the default symmetric coding gives every category a forest and a
shift common to all of them cancels out of the softmax. A vector offset
there leaves the fitted probabilities unchanged; with `reference` set it
does not cancel, and moves every non-reference category against the
reference one. A matrix is what expresses a per-category offset either
way.

An offset is not a function of the predictors, so it cannot be rebuilt
for rows the fit has not seen: a model fitted with one requires `offset`
at [`predict()`](https://rdrr.io/r/stats/predict.html) time.

### Missing Predictor Values

A missing predictor is not imputed and its row is not dropped, which is
the default here because a tree can do something better with a missing
value than either. Instead each splitting rule carries the answer for
itself. A rule on a variable that has missing values is drawn as one of
three, with equal probability:

- `x < c`, or missing, goes left;

- `x < c` goes left, missing goes right;

- missing goes left, present goes right.

This is missingness incorporated in attributes, and the third rule is
what lets the model split on missingness itself, so a variable whose
absence carries the signal is usable even where its observed values say
nothing.

Two consequences are worth being clear about.
[`predict()`](https://rdrr.io/r/stats/predict.html) accepts missing
values in a column that had them at fitting time, those being the
columns whose rules carry an answer. And what the model estimates is the
mean of the response given the predictors and the pattern of
missingness, which is the quantity prediction calls for; where the
estimand is a regression or causal effect defined on complete data,
multiple imputation is the right tool.

### Preprocessing

Predictors are mapped to the unit interval, because the cutpoint prior
is uniform on a node's live range and the soft-rule bandwidth is
measured on the predictor scale. Factors share a single weight in the
sparsity prior, so that a factor is selected or not as a whole rather
than one level at a time. The additive predictor starts from an
intercept-only fit, so the leaf prior describes departures from that fit
rather than the absolute level of the response. That starting value is
the exact null-model estimate for most families; for the accelerated
failure time families, where censoring makes the sample mean of the log
times biased, and for the zero-inflated and ordered beta families, it is
a moment approximation, which the sampler then moves away from.

### Drawing From the Prior (`prior_only`)

`prior_only = TRUE` fits the same model to no data. Every observation is
given a weight of zero, and since the weight multiplies that
observation's log density, its gradient, and its curvature, the
likelihood is flat: each tree move is accepted or rejected on the prior
alone and each leaf is drawn from its prior. A family that draws an
auxiliary parameter from the response directly rather than through the
weighted density (the mixture atoms under
[`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
the latent utilities under `multinomial(link = "probit")`) is told
separately that a weightless observation carries no information, and
draws that parameter from its own prior instead. The sampler is
otherwise untouched, so what comes back is an ordinary fit whose draws
are prior draws, and
[`rstantools::posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
on it gives the prior predictive distribution.

It answers a question the priors themselves cannot. `k`, `gamma` and
`beta` are statements about trees and leaves, and what they imply about
an outcome is opaque. The replicates put that on the response's own
scale, where it can be judged: a prior predictive that puts its mass
where the outcome cannot go, or spread over an implausible range, is a
prior worth changing before the data are seen and before any of their
information is spent.

The additive predictor is anchored at an intercept-only fit on the link
scale and the leaf scale is calibrated from the response, which is how a
BART prior is specified. The replicates therefore take their location
and scale from the response and everything else from the prior: which
predictors are split on, how deep, how far the fitted function departs
from that anchor. Read them for shape and spread rather than for level,
and note that the wider a prior is the wider its replicates, so
[`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
and
[`Gamma_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
which put a log scale in a second forest, run far wider than the
response ever does.

`loo()`, `waic()` and `kfold()` score a fit against data, so they refuse
a prior-only fit, and
[`performance::model_performance()`](https://easystats.github.io/performance/reference/model_performance.html)
leaves out the three columns built on them.

## References

Linero, A. R. (2025). Generalized Bayesian additive regression trees
models: beyond conditional conjugacy. *Journal of the American
Statistical Association*, 120(549), 356–369.
[doi:10.1080/01621459.2024.2337156](https://doi.org/10.1080/01621459.2024.2337156)

Linero, A. R., & Yang, Y. (2018). Bayesian regression tree ensembles
that adapt to smoothness and sparsity. *Journal of the Royal Statistical
Society Series B*, 80(5), 1087–1110.
[doi:10.1111/rssb.12293](https://doi.org/10.1111/rssb.12293)

## See also

[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
for the sampler and prior settings;
[`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
for prediction;
[`bartisan-families`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
for the likelihoods, and
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
for a family-by-family guide;
[`bartisan-marginaleffects`](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
for reading effects off a fit

## Examples

``` r
data("rhc")
set.seed(123)

# Whether a patient died, with every other variable a candidate predictor
# and the family read off the response. `days` is the timing of the same
# event, so it is excluded rather than conditioned on
fit <- bartisan(death ~ . - days, data = rhc,
                num_trees = 10, num_burn = 50, num_draws = 50,
                verbose = FALSE)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` to choose another, which also silences this message.
fit
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = death ~ . - days, data = rhc, num_trees = 10, 
#>     num_burn = 50, num_draws = 50, verbose = FALSE)
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 1500
#> Structure: 1 forest of 10 trees, soft decision rules
#> Draws: 50 kept after 50 warmup

# Fitted probabilities
head(predict(fit, type = "response"))
#> [1] 0.7685691 0.8202057 0.2465819 0.3451376 0.3927522 0.5274813

# The forest has no coefficients, so an effect is a contrast of
# predictions, here of catheterization on the probability of death
if (rlang::is_installed("marginaleffects")) {
  marginaleffects::avg_comparisons(fit, variables = "rhc")
}
#> 
#>  Estimate 2.5 % 97.5 %
#>    0.0585     0  0.102
#> 
#> Term: rhc
#> Type: response
#> Comparison: 1 - 0
#> 
```
