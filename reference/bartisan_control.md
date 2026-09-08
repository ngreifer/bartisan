# Sampler and prior settings for `bartisan()`

Collects the tuning constants of the sampler and the hyperparameters of
the tree prior. The defaults follow Linero (2025) and, for the soft
decision rules, Linero and Yang (2018), and are intended to be usable
without adjustment.

The arguments fall into three groups, and the first group is the one
worth reading: `num_trees`, `gate`, `sparsity`, `k`, `bandwidth`, the
three chain lengths, and `x_transform` are **modeling decisions**, in
that changing one changes what is being fitted or how long it is fitted
for. Everything from `augment` to `num_print` is an **advanced
setting**: a hyperparameter of a prior the first group summarizes, or a
switch whose default is almost always right. The last three
(`block_eval`, `exact_quadratic`, and `generic_accumulate`) exist **for
internal validation** and are documented so that the checks that use
them can be read; they compute the same posterior more slowly.

## Usage

``` r
bartisan_control(
  num_trees = NULL,
  gate = "smoothstep",
  sparsity = TRUE,
  share_sparsity = FALSE,
  split_prior = NULL,
  categorical = "subset",
  k = 2,
  bandwidth = 0.1,
  chains = 1L,
  num_burn = 200L,
  num_draws = 800L,
  num_thin = 1L,
  augment = TRUE,
  x_transform = "quantile",
  gamma = 0.95,
  beta = 2,
  sigma_mu = NULL,
  update_sigma_mu = TRUE,
  sigma_mu_ramp = 0.25,
  update_tau = TRUE,
  update_bandwidth = TRUE,
  bandwidth_every = 1L,
  alpha = NULL,
  alpha_scale = NULL,
  alpha_shape_1 = NULL,
  alpha_shape_2 = NULL,
  update_s = NULL,
  update_alpha = NULL,
  verbose = FALSE,
  num_print = 100L,
  block_eval = FALSE,
  exact_quadratic = TRUE,
  generic_accumulate = FALSE
)
```

## Arguments

- num_trees:

  `numeric`; the number of trees, given either as one number for every
  forest or as one number per additive predictor. Default is `NULL`,
  which is 50 for both kinds of rule. The two do not want the same count
  (a soft rule makes each tree more expressive, so soft rules reach
  their best held-out error at about 20 trees and get worse past that,
  while hard rules keep improving to 200), but a smaller forest mixes
  worse, which is why 50 is the compromise for both. A family with more
  than one additive predictor takes a vector, which is worth using:
  [`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  spends about 90% of its time on the scale forest, and a variance
  surface needs less capacity than a mean surface, so
  `num_trees = c(50, 10)` runs about 2.5 times faster at the same
  accuracy. All of this is measured in Details.

- gate:

  string; the shape of a decision rule, which is also how hard and soft
  rules are chosen between. Allowable options include `"smoothstep"`
  (the default), `"smootherstep"`, `"logistic"`, and `"hard"` (or
  equivalently `"step"`). `"hard"` gives the step functions of standard
  BART. The other three give soft rules, as in the SoftBart model, in
  which every observation reaches every leaf with some weight and the
  fitted function is smooth: `"smoothstep"` is the Beta(2, 2) cumulative
  distribution function, once differentiable and supported on a bounded
  interval; `"smootherstep"` is Beta(3, 3), twice differentiable and
  bounded; and `"logistic"` is the logistic function of Linero and Yang
  (2018), which is infinitely differentiable. Soft rules cost three to
  five times as much per iteration and cut held-out error by 35 to 40%,
  so they are the accuracy argument rather than a tax. The two bounded
  gates are about 1.4 times faster than the logistic and equally
  accurate, and are within noise of each other; see Details.

- sparsity:

  the prior on which predictors are split on, given as either a logical
  value or a string. `TRUE` (the default) is the Dirichlet sparsity
  prior of Linero (2018) (i.e., DART), which concentrates splits on the
  predictors that earn them and can drop the rest from the forest
  entirely, and `FALSE` gives every predictor the same splitting
  probability, which is classic BART. The strings `"none"`, `"weak"`,
  `"moderate"`, and `"strong"` name four strengths, with `"none"` equal
  to `FALSE`, `"moderate"` equal to `TRUE`, and the other two moving the
  prior on the concentration. Note that this argument sets `update_s`,
  `update_alpha`, `alpha_shape_1`, and `alpha_shape_2` together, and
  that supplying any of those directly overrides it. Read the trade-off
  in Details before turning it off or up.

- share_sparsity:

  `logical`; for a family with more than one additive predictor, whether
  the forests draw their splitting proportions from one pooled Dirichlet
  instead of one each. Default is `FALSE` for a separate prior per
  forest. `TRUE` says the same predictors are relevant to every
  component, which is an assumption about the data rather than a free
  improvement; see Details. It requires that the proportions be drawn at
  all, so `sparsity` must not be `FALSE` for the forests that are to
  share, and that those forests be able to split on the same predictors.
  Ignored rather than an error for a family with a single forest, where
  there is nothing to share.

- split_prior:

  `numeric`; the relative prior weight on each predictor, for use when
  some are expected to matter more than others. This is a named vector,
  keyed by the names the predictors have in the formula, and every
  predictor not named gets a weight of 1. The prior probability of
  splitting on a predictor is its weight divided by the total, so in a
  three-predictor model `split_prior = c(x1 = 3, x3 = 0.5)` gives `x1` a
  probability of `3 / 4.5`, `x2` one of `1 / 4.5`, and `x3` one of
  `0.5 / 4.5`. Weights must be finite and nonnegative, and naming a
  predictor the model does not have is an error rather than being
  silently ignored. A weight of zero is allowed and means the predictor
  is never split on: it stays in the model frame and out of every tree.
  Default is `NULL` to weight every predictor equally. Note that setting
  this overrides `sparsity`; see Details.

- categorical:

  string; how a splitting rule divides the levels of a factor. Allowable
  options include `"subset"` (the default), which draws a subset of the
  levels still available at the node and sends those left, as in
  Deshpande (2024), and `"onehot"`, which is what most BART
  implementations do: it splits on one indicator column, peeling a
  single level off the rest. The choice matters because `"onehot"`
  reaches only `2^K - K` of the `B_K` partitions of `K` levels (27 of 52
  at `K = 5`, and 1,014 of 115,975 at `K = 10`), and the partitions it
  can form all have at most one cell with more than one level in it, so
  the bulk of the levels is never divided. See Details.

- k:

  `numeric`; controls the leaf prior. The prior standard deviation of a
  forest is `3 / k` times the natural scale of its additive predictor,
  so larger values shrink the fit harder toward the intercept-only
  model. Default is 2.

- bandwidth:

  `numeric`; the prior mean of the gate bandwidth of a soft rule, on the
  scale of the transformed predictors, which lie in `[0, 1]`. Smaller
  values approach hard rules. Default is .1. Ignored when
  `gate = "hard"`.

- chains:

  `numeric`; how many independent chains to run. Default is 1. The draws
  are pooled and
  [split-R-hat](https://ngreifer.github.io/bartisan/reference/bartisan.md)
  is reported in the `rhat` element. With the future.apply package
  installed the chains run in parallel under whatever backend the caller
  has planned with
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  (`multisession`, `multicore`, a cluster, or mirai's
  `mirai_multisession`); without it they run one after another, which is
  slower and otherwise identical. One
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  reproduces the whole run either way, because each chain is given its
  own L'Ecuyer stream.

- num_burn:

  `numeric`; the number of warmup iterations to discard. Default is 200.
  Warmup is where the trees grow into the data and the hyperparameters
  find their scale, so raising it buys convergence rather than
  precision: increase it when `rhat` says the chains have not agreed. It
  was 500 through earlier development, which measurement says is several
  times longer than anything here needs; see Details for what was
  measured, and raise it for a
  [`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  fit, which is the one family that showed a cost.

- num_draws:

  `numeric`; the number of draws to keep. Default is 800. These are what
  every estimate and interval is computed from, so raising it narrows
  Monte Carlo error and does nothing about convergence: increase it when
  `ess_bulk` or `ess_tail` is small relative to what the reported
  quantity needs.

- num_thin:

  `numeric`; keep one draw in every `num_thin` after warmup. Default is
  1 to keep every draw. Thinning discards draws to make the kept ones
  less correlated, which costs information and is worth it only to hold
  down the memory a long chain would otherwise take: for a given amount
  of computing, more draws beat fewer less-correlated ones.

- augment:

  *Advanced.* Whether to rewrite the likelihood as the margin of a
  Gaussian one, or of a Poisson one, which makes the target a shape the
  sampler can exploit and the Laplace approximation exact or nearly so.
  **The posterior is the same either way**, so this is a sampling
  setting rather than a modeling one. Default is `TRUE`, which does it
  wherever it has been measured to pay: the binomial, ordinal,
  multinomial, zero-inflated, and survival families always, and the
  negative binomial when the rules are hard. `FALSE` never does it, and
  a character vector of engine family names (`"binomial"`, `"ordinal"`,
  `"multinomial"`, `"negbin"`, `"zip"`, `"zinb"`, or `"aft"`) asks for
  exactly those. Every rewriting trades speed for mixing, so the
  measured effect on effective sample size per second is what matters,
  and it differs by family; see Details. In a fitted model,
  `control$augment` is instead a `logical` recording whether a rewriting
  was applied, since whether the data admit one (a Bernoulli response,
  single trials) is settled only when the model is fitted.

- x_transform:

  string; how numeric predictors are mapped to `[0, 1]`. Allowable
  options include `"quantile"` (the default), which uses each
  predictor's empirical distribution function and so makes the cutpoint
  prior invariant to monotone reparameterization, and `"range"`, which
  rescales linearly and preserves the original spacing.

- gamma, beta:

  *Advanced.* `numeric`; the branching probability at depth `d` is
  `gamma * (1 + d)^(-beta)`. Defaults are .95 and 2.

- sigma_mu:

  *Advanced.* `numeric`; the prior median of the leaf standard
  deviation, one value per additive predictor. Default is `NULL` to
  derive it from `k` and that forest's own tree count.

- update_sigma_mu:

  *Advanced.* `logical`; whether to draw the leaf standard deviation
  under a half-Cauchy prior rather than fixing it. Default is `TRUE`.
  `FALSE` is worth reaching for when a binary fit mixes badly: where the
  predictors separate the response well the leaf scale is barely
  identified, and it wanders rather than settling, which drags the
  effective sample size of everything built on it down with it. On the
  propensity model of
  [`vignette("causal")`](https://ngreifer.github.io/bartisan/articles/causal.md)
  it takes the additive predictor's effective sample size from 91 to 301
  and the log likelihood's from 18 to 249, for about 1% of held-out AUC.

- sigma_mu_ramp:

  *Advanced.* `numeric`; the fraction of warmup over which the leaf
  standard deviation is raised from near zero to its target. Default is
  .25; set to 0 to disable. Linero (2025) describes this as essential:
  started at its full value, the sampler can settle early into a poor
  configuration and fail to move.

- update_tau:

  *Advanced.* `logical`; whether to draw the standard deviation of each
  random-effect term under the same half-Cauchy prior the leaf scale
  uses, rather than fixing it at that prior's median. Default is `TRUE`.
  Only relevant when the formula has a `(1 | group)` term.

- update_bandwidth:

  *Advanced.* `logical`; whether to draw the bandwidth of each tree
  rather than holding it at `bandwidth`. Default is `TRUE`.

- bandwidth_every:

  *Advanced.* `numeric`; how many sweeps between bandwidth draws for a
  given tree. Default is 1. The bandwidth is one scalar per tree, drawn
  by an adaptive random walk, and every attempt costs a full rebuild of
  the tree's memberships, which is the single largest item in a
  soft-rule fit. Drawing it less often than the trees themselves trades
  mixing in that one parameter for time; see Details.

- alpha:

  *Advanced.* `numeric`; the concentration of the Dirichlet prior on the
  splitting proportions, where smaller values concentrate splits on
  fewer predictors. Default is `NULL`, which uses 1 as the starting
  value for a parameter that is then drawn.

- alpha_scale, alpha_shape_1, alpha_shape_2:

  *Advanced.* `numeric`; the prior on `alpha`, in which
  `alpha / (alpha + alpha_scale)` is Beta(`alpha_shape_1`,
  `alpha_shape_2`). Default for `alpha_scale` is `NULL` to use the
  number of predictor groups; the two shapes default to whatever
  `sparsity` implies.

- update_s, update_alpha:

  *Advanced.* `logical`; whether to draw the splitting proportions and
  their concentration. Defaults are `NULL` to follow `sparsity`. Turning
  both off recovers a uniform prior over predictors, which is what
  `sparsity = FALSE` does.

- verbose:

  *Advanced.* `logical`; whether to print progress to the console while
  sampling. Default is `FALSE`. For a progress bar instead, see the
  Progress section below, which needs no argument here.

- num_print:

  *Advanced.* `numeric`; how many iterations between the reports
  `verbose` prints. Default is 100.

- block_eval:

  *Validation.* `logical`; whether to evaluate the likelihood one leaf
  at a time rather than one observation at a time. Default is `FALSE`. A
  family built by
  [`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  does this regardless, since it must call back into R; setting it for a
  compiled family produces the same draws from the same seed at somewhat
  greater cost.

- exact_quadratic:

  *Validation.* `logical`; whether to use the closed forms that a target
  quadratic in the additive predictor allows, in which one pass over a
  node determines the log target everywhere, so that the Laplace
  approximation is the conditional posterior rather than an
  approximation to it. Default is `TRUE`. This is what a conjugate
  sampler does, and it is what makes a Gaussian response, or any of the
  rewritings in `augment`, cheap. Setting it to `FALSE` falls back on
  the general path; the two agree, at greater cost.

- generic_accumulate:

  *Validation.* `logical`; whether to accumulate a leaf's sums through
  the family's virtual interface rather than through its own statically
  dispatched loop. Default is `FALSE`. The two compute the same thing;
  the second lets the compiler inline the family's arithmetic, which is
  most of the remaining per-observation cost.

## Value

A `<bartisan_control>` object, a list containing the supplied settings
and the defaults for those not supplied, for passing to the `control`
argument of
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).

## Details

### How Long Warmup Needs to Be

`num_burn` was 500 through earlier development and is now 200, because
500 is several times longer than the sampler takes to settle and warmup
is time that produces no draws.

**The transient is short.** Run with no warmup at all, so that every
sweep is retained, the log likelihood reaches within two standard
deviations of its eventual level by sweep 36 under soft rules, 55 under
hard ones, and 61 for a hard-rule
[`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
fit. On the harder case for warmup, a sparse one with 30 predictors of
which 25 are irrelevant, the share of splits landing on the five that
matter climbs from 0.20 to its plateau by sweep 34 under soft rules and
48 under hard, which is the variable-selection state settling and is the
slowest thing warmup has to do.

**Shortening it costs nothing measurable, and buys effective sample size
per second.** Out-of-sample root mean squared error against the true
regression function, four chains of 500 draws, against the same fit at
`num_burn = 500`. A positive difference is worse; the standard error is
of the paired difference.

|  |  |  |  |
|----|----|----|----|
| design | 100 | 200 | paired SE |
| Gaussian, soft, n = 1500, p = 10 | -0.001 | +0.001 | 0.002 |
| Gaussian, soft, n = 500, p = 30 | +0.003 | -0.002 | 0.003 |
| Gaussian, hard, n = 1500, p = 10 | -0.001 | -0.002 | 0.003 |
| [`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), hard | -0.021 | -0.007 | 0.011 |
| [`negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | -0.002 | +0.001 | 0.003 |
| [`binomial()`](https://rdrr.io/r/stats/family.html) | -0.019 | -0.016 | 0.009 |
| [`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | **+0.011** | **+0.006** | 0.002 |

Every design but the last is within a standard error or two of the
longer warmup, and several are better with the shorter one. Effective
sample size per second improves everywhere, by 1.4 to 2.0 times, because
the sweeps saved were producing nothing. `num_draws` went from 500 to
800 at the same time, which spends some of what warmup gave back on
draws that do count towards an effective sample size; the two together
still run in less time than the old pair did.

**[`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is the exception.** Its mixture carries a component-count state that
settles more slowly than a forest does, and the cost is real rather than
noise: 4 standard errors at 200 and 4 at 100. It is small in absolute
terms, about 2.5% of the error at 200, and `num_burn = 500` is the
setting to reach for when the error distribution itself is the object of
interest.

Raising `num_burn` is still what a large `rhat` calls for. What the
measurement says is that 500 is not a good place to start from.

### How Many Trees, and How Many per Forest

The tree count is the setting most worth thinking about after the
family, and the default of 50 is a compromise rather than an optimum.
Measured on the Friedman function with 1000 training and 1000 test
observations, four chains of 500 draws after 500 warmup iterations, as
held-out root mean squared error against the true regression function:

|            |         |       |           |           |       |       |
|------------|---------|-------|-----------|-----------|-------|-------|
| Rules      | 5 trees | 10    | 20        | 50        | 100   | 200   |
| Soft rules | 0.286   | 0.281 | **0.270** | 0.284     | 0.289 | 0.319 |
| Hard rules | 1.149   | 0.682 | 0.558     | **0.521** | 0.531 | 0.510 |

Two things to read off it. **Soft rules need far fewer trees than hard
ones**, which is what makes 200 (the default in most BART packages)
actively worse here than 20. And **the two kinds of rule want different
counts**, since hard rules are still improving at 200 where soft rules
peaked at 20.

The default is 50 for both anyway, for a reason the table cannot show: a
smaller forest mixes worse. On `MatchIt::lalonde`, four chains at 20
soft trees disagreed by 35% on an average contrast where 50 trees
disagreed by 9%, and the Friedman gain at 20 trees was 5%. So 50 buys
reliable inference at a small cost in point accuracy. Drop to 20 if
prediction is the only goal and the fit is soft; raise toward 200 with
hard rules if it is not.

**A vector is worth using when the family has more than one forest.**
They are not equally expensive and they do not need equal capacity.
[`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
spends about 90% of its time on the scale forest, because that forest's
target is not quadratic, and a variance surface carries much less
information than a mean surface. Measured on 1000 observations with a
smooth mean and a log-linear standard deviation:

|             |         |           |             |           |
|-------------|---------|-----------|-------------|-----------|
| `num_trees` | Seconds | Mean RMSE | Log-SD RMSE | Log score |
| `c(50, 50)` | 14.5    | 0.092     | 0.050       | -1188     |
| `c(50, 20)` | 8.3     | 0.094     | 0.047       | -1188     |
| `c(50, 10)` | 5.9     | 0.093     | 0.051       | -1188     |
| `c(50, 5)`  | **4.9** | 0.094     | 0.046       | -1187     |
| `c(20, 5)`  | **2.9** | 0.084     | 0.041       | -1184     |

A Gaussian fit on the same data takes 1.4 seconds, so `c(50, 50)` costs
10 times a Gaussian fit and `c(50, 5)` costs 3.5 times, at the same
accuracy to three decimal places. That is not the default, because how
many trees a variance surface needs depends on how complicated it is,
and silently under-parameterizing it would show up as intervals that are
wrong, which is the thing
[`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
exists to get right. It is worth setting by hand.

The leaf prior scale divides by the square root of each forest's own
tree count, so shrinking one forest does not change the prior on the sum
it forms.

### The Sparsity Prior, and What It Costs

`sparsity = TRUE` is the Dirichlet prior of Linero (2018) on the
splitting proportions. It is a genuine variable-selection prior: it can
and does drop a predictor from every tree in the forest at once, and
that is the point of it in the high-dimensional problems it was built
for.

It has a consequence that is easy to misread. **A contrast on a
predictor the prior has dropped is exactly zero, not nearly zero**,
because in that draw the fit does not depend on that predictor at all.
So the posterior of a contrast has an atom at zero whose mass is one
minus the predictor's inclusion probability, which
[`summary()`](https://rdrr.io/r/base/summary.html) reports as
`prop_used`. Any summary that reports a median lands on that atom
whenever it holds half the mass; see
[bartisan-marginaleffects](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md).

**More trees does not fix it**, which is worth stating because it looks
as though it should. On `MatchIt::lalonde`, four chains with the
sparsity prior on left `treat` out of every tree in 20% of draws at 50
trees and 18% at 200, and disagreed with each other by more than 100% of
the estimate at every count from 20 trees to 200. With
`sparsity = FALSE` the atom disappears by 50 trees and the four chains
agreed to within 9%. The variable-selection state mixes slowly: a
predictor whose splitting proportion has gone small is rarely proposed,
so it is hard to get back in.

**What each setting is worth, measured.** The trade is real in both
directions, and the strength barely matters in either: what matters is
whether the prior is on at all.

For prediction, any sparsity beats none. Scored against the true
regression function on held-out data, Friedman with five relevant
predictors:

|              |               |               |
|--------------|---------------|---------------|
| Sparsity     | 10 predictors | 50 predictors |
| `"none"`     | 0.446         | 0.465         |
| `"weak"`     | 0.385         | 0.346         |
| `"moderate"` | 0.400         | 0.372         |
| `"strong"`   | 0.374         | 0.362         |

For a contrast on a predictor whose signal is weak, any sparsity is
actively harmful, and not only in the atom-at-zero sense above. A binary
treatment among 20 predictors, continuous outcome, residual standard
deviation 1, n = 800, five replicates, with `covers` the share of
replicates whose 95% interval contains the truth:

|             |         |          |      |        |
|-------------|---------|----------|------|--------|
| true effect | setting | estimate | atom | covers |
| 0.05        | `FALSE` | 0.031    | 0.08 | 0.80   |
| 0.05        | `TRUE`  | 0.000    | 0.89 | 0.40   |
| 0.10        | `FALSE` | 0.131    | 0.03 | 1.00   |
| 0.10        | `TRUE`  | 0.029    | 0.69 | 1.00   |
| 0.20        | `FALSE` | 0.161    | 0.05 | 1.00   |
| 0.20        | `TRUE`  | 0.094    | 0.55 | 0.60   |
| 0.50        | `FALSE` | 0.475    | 0.00 | 1.00   |
| 0.50        | `TRUE`  | 0.474    | 0.00 | 1.00   |

The prior attenuates a weak effect by half or more and its interval
covers well below its nominal rate. A strong effect is untouched,
because the prior never has reason to drop a predictor that is earning
its splits, so this is a weak-signal failure rather than a general one.

So the setting is close to a switch, and which way to throw it follows
from the estimand rather than from the data.

**For prediction, or for variable selection**, keep the default.
`"weak"` is already worth most of what `"strong"` is worth, so reaching
past `TRUE` is warranted only when the predictors are many and nearly
all of them are expected to be irrelevant.

**For a contrast, a partial effect, or a treatment effect**, set
`sparsity = FALSE`, or use `split_prior`, which fixes the weights and so
cannot drop anything. `split_prior` is preferable when the other
predictors are numerous enough that weighting them all alike is
wasteful.

**For a varying-coefficient model, the two answers can differ by
forest.** What the prior can drop is a predictor a forest splits on, and
in a [`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) model
the treatment is not one of those: it is the coefficient, carried by a
forest of its own, so no splitting proportion can drop it. The prior on
that forest selects among the moderators instead, and dropping all of
them leaves a coefficient that does not vary rather than one that is
zero. So `sparsity = c(FALSE, TRUE)` is a coherent thing to ask for, and
in [`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) it is
the asymmetry worth considering: the control function keeps every
predictor, the propensity score included, and the effect forest is left
to work out which covariates moderate.

Run several chains either way, because one chain can look far more
settled than the posterior is.

### Sharing the Sparsity Prior Across Forests

A family with more than one additive predictor fits a forest for each,
and by default each forest draws its own splitting proportions. That is
the conservative choice and it wastes information whenever the
components are functions of the same predictors, which they usually are:
the zero-inflation probability and the count mean of
[`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
are two aspects of one process, and so are the mean and the spread of
[`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md).
One component is often far more informative about *which* predictors
matter than another, and with separate priors the less informative one
has to find them again from its own residual signal alone.

`share_sparsity = TRUE` pools the splitting counts of the forests and
draws one Dirichlet from the total, so every forest reaches for
predictors in the same proportions. The forests themselves stay
separate: they have their own trees, their own cut points, and their own
leaf scales, and only the prior over which predictors they may use is
held in common.

Whether it helps has a clear shape, and the trade is asymmetric.
Measured on
[`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
at `n = 400` over 20 replicates, scoring each component against the
truth on a held-out thousand, with the mean and the log standard
deviation driven either by the same five predictors or by disjoint sets
of five:

|  |  |  |  |
|----|----|----|----|
|  | predictors | weaker component | stronger component |
| [`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), same 5 | 5 | 1.02x | 0.99x |
| [`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), same 5 | 25 | 1.13x | 1.03x |
| [`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), same 5 | 100 | **1.23x** | 1.12x |
| [`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), disjoint | 100 | 0.98x | 0.96x |
| [`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), same 5 | 5 | 1.04x | 0.97x |
| [`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), same 5 | 25 | 1.30x | 1.03x |
| [`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), same 5 | 100 | **1.40x** | 1.02x |
| [`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), disjoint | 100 | 1.00x | 0.92x |

Ratios above 1 are reductions in root mean squared error against the
default. The weaker component is the log standard deviation and the
zero-inflation probability respectively, which are the ones with less
signal to find the relevant predictors from on their own. Sharing buys
nothing at five predictors, because there is no selection problem to
transfer, and takes a quarter to two fifths off the weaker component's
error at a hundred. When the assumption is false it costs 2% to 8%,
since the pooled prior pulls each forest towards the other's variables.
Turning it on is therefore a statement about the data, and one whose
downside is a good deal smaller than its upside.
[`variable_importance()`](https://ngreifer.github.io/bartisan/reference/variable_importance.md)
is where to check it: with a shared prior the forests report similar
`prop_splits`, and a fit that wants them different will show that under
the default.

This is the variable-selection content of the shared forests of Linero
et al. (2020) and not their model. Theirs shares the tree *topology*
across components, so the forests have the same partitions with
different leaf values, which fixes the cut points as well as the choice
of predictor. The pooled prior here is the weaker and cheaper
assumption;
[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) sits at
the other extreme, sharing an entire function between components rather
than a prior.

### Arguments That Vary by Forest

A family with several additive predictors has one forest per predictor,
each with its own prior. Every argument that could mean something
different for one of them may be given once, to apply to all, or once
per forest, either positionally or keyed by the forest names listed in
[bartisan-families](https://ngreifer.github.io/bartisan/reference/bartisan-families.md).
Note that a forest a named argument does not mention keeps that
argument's default rather than borrowing the value chosen for another
forest.

That covers `num_trees`, `k`, `sigma_mu`, `sparsity`, `split_prior`,
`bandwidth`, `gamma`, `beta`, the four `alpha` arguments, and the three
`update_` flags for the leaf scale, the splitting proportions and the
bandwidth. `formula` works the same way; see
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).

    # A scale forest with less capacity than the mean forest, and no sparsity
    # prior on it.
    bartisan_control(num_trees = c(mean = 50, log_sd = 10),
                     sparsity = c(mean = TRUE, log_sd = FALSE))

The multinomial families are the exception, for the reason given in
[bartisan-families](https://ngreifer.github.io/bartisan/reference/bartisan-families.md):
their forests act as one, so these arguments take a single value.

### Splitting a Factor

A rule on one indicator column of a factor can only separate one level
from the rest. Applied repeatedly down a path that produces a partition
with some number of singleton levels and one cell holding everything
else, and those are the only partitions available: `2^K - K` of the
`B_K` partitions of `K` levels, which is 27 of 52 at `K = 5` and under
1% at `K = 10`.

What that costs is partial pooling. Simulating the tree prior directly
at `K = 10`, the probability that two levels land in the same leaf is
0.77 under `"onehot"` against 0.46 under `"subset"`, and a typical
`"onehot"` tree has 1.18 singleton levels out of 2.18 leaves: one level
alone, the rest together, whether the data want that or not. A rule that
takes a subset of the levels reaches every partition, and the
co-clustering probabilities come out far more even.

**What it is worth depends on the decision rules.** Twenty levels in
four groups of five sharing a mean, twelve replicates, paired within
replicate, RMSE against the true mean function:

|              |              |              |               |
|--------------|--------------|--------------|---------------|
| Rule         | 10 per level | 25 per level | 100 per level |
| subset, hard | 0.3705       | 0.2617       | 0.1421        |
| onehot, hard | 0.3998       | 0.2725       | 0.1488        |
| subset, soft | 0.3435       | 0.2265       | 0.1098        |
| onehot, soft | 0.3445       | 0.2201       | 0.1059        |

Under hard rules `"subset"` is better at every size, clearly so at ten
observations per level, where the gap is 7% against a standard error of
2%, and by 4% at the two larger sizes, where the standard error is
around half the gap. Under soft rules, which is the default, **the two
are indistinguishable**: the largest gap is 0.006 against a standard
error of 0.005. Soft rules are worth far more than either choice, which
is the biggest number in the table and the one to act on.

So `"subset"` is the default because it is right about the prior and
never loses beyond noise, not because it will visibly improve a fit.

Two consequences of `"subset"` worth knowing. A rule on a factor is
always hard, even in a soft tree: a gate is a smooth function of the
distance from a cutpoint and there is no distance between two levels.
Under `"onehot"` the distance-based gate did apply to the 0 and 1 of an
indicator, and at the default bandwidth it left one level with a
fractional membership weight for 81% of cutpoints, which is not a
smoothing anyone asked for. And a rule can no longer land on an
indicator that a path has already used up, which `"onehot"` did for
between 1% and 7% of its draws on a factor depending on `K`.

A two-level factor is unaffected either way, since `2^2 - 2` is `B_2`.

### Telling the Prior What Is Already Known

`sparsity` and `split_prior` answer different questions and cannot both
be in force, so giving `split_prior` turns `sparsity` off. This is about
weights a caller supplies, and not about a forest being held to the
predictors its own formula names: a
[`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) term whose
moderators are only some of the covariates still draws its splitting
proportions, over those moderators. `sparsity` is for the case in which
which predictors matter is unknown and the prior is to work it out from
the data; the splitting proportions are drawn, and a predictor can be
dropped from the forest entirely. `split_prior` is for the case in which
something is known and is to be honored; the proportions are held at the
supplied weights and nothing draws over them.

A weight is a statement about relative attention, not about effect size.
It changes how often the sampler proposes a split on a predictor, which
is a prior, so the data can still overrule it: a predictor given a large
weight whose splits do not improve the fit will collect rules that go
nowhere, and one given a small weight that genuinely matters will still
be found, more slowly. On pure noise, where nothing in the data prefers
any predictor, the realized share of splitting rules matches the weights
closely, which is the cleanest way to see what the argument does and the
only case where the weights fully determine the outcome.

Because the weights are fixed, `split_prior` does not accumulate the
atom-at-zero mass described above. A predictor can still miss out on a
rule in some draw of a small forest, but that is sampling variation at a
fixed probability rather than a prior state that persists, and unlike
the sparsity prior it does go away as trees are added: on four noise
predictors with `split_prior = c(x1 = 8)`, `prop_used` runs from 0.89 to
1 at 20 trees and is 1 for every predictor at 50 and at 200. The
sparsity prior on the same data at 50 trees left the most heavily
weighted predictor out of 41% of draws. So `split_prior` is a reasonable
middle course when a particular contrast is the estimand but the
predictors are still too many to treat alike.

One weight per term in the formula, not per column of the design matrix,
so a factor is named once and its levels share the weight, the way they
already share one entry of the sparsity prior.

### Gaussian Rewritings

The expensive part of this sampler is not the non-conjugacy of the
likelihood but that the leaf-level target is not quadratic in the
additive predictor. Where a data augmentation makes it quadratic, Fisher
scoring reaches the mode in one step, the fitted normal *is* the
conditional posterior, and the leaf refresh becomes a Gibbs step with
acceptance one.

Every such rewriting costs mixing, because the chain now has to move the
latent variables and the predictor in alternation. Measured on 500
observations, 50 trees and 800 draws, with the ratio that matters being
the last column:

|  |  |  |  |  |
|----|----|----|----|----|
| Family | Rules | Speed | Effective sample size | ESS per second |
| `binomial("probit")` | either | 5.7 to 7.2x | 0.66 to 0.75x | **3.8 to 5.3x** |
| `binomial("logit")` | either | 2.6 to 3.7x | 0.81 to 1.8x | **2.1 to 6.7x** |
| `ordinal("probit")` | soft | 14x | 0.79 to 0.95x | **11 to 13x** |
| `ordinal("probit")` | hard | 26 to 30x | 0.73 to 0.90x | **22 to 24x** |
| `ordinal("logit")` | soft | 7.0 to 7.2x | 0.70 to 0.82x | **5.0 to 5.9x** |
| `ordinal("logit")` | hard | 14.5 to 15.3x | 0.87 to 1.02x | **13 to 15x** |
| `ordinal("cloglog")` | soft | 2.8x | 0.57x | 1.6x |
| `ordinal("cloglog")` | hard | 5.1x | 1.05x | **5.2x** |
| [`negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | hard | 1.7 to 1.9x | 0.61 to 1.14x | 1.2 to 2.0x |
| [`negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | soft | 1.1x | 0.71x | 0.8x |
| [`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | soft | 9.3x | 1.09x | **10.1x** |
| [`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | hard | 14.5x | 0.66x | **9.6x** |
| [`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | soft | 4.6x | 0.85x | **3.9x** |
| [`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | hard | 7.1x | 1.42x | **10.1x** |
| [`zi_negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | soft | 5.9x | 0.98x | **5.8x** |
| [`zi_negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | hard | 9.4x | 0.84x | **7.9x** |
| [`lognormal_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | soft | 11.2 to 19.7x | 0.83 to 1.27x | **14 to 16x** |
| [`lognormal_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | hard | 19.7 to 29.3x | 0.74 to 1.01x | **20 to 22x** |
| [`loglogistic_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | soft | 8.4 to 9.3x | 0.58 to 0.84x | **4.9 to 7.7x** |
| [`loglogistic_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | hard | 11.3 to 12.3x | 0.87 to 0.89x | **9.8 to 11x** |

The ranges are two problems of different size and shape, which is a fair
picture of how much this varies: what an augmentation costs in mixing
depends on the data, not only on the family. The negative binomial is
the marginal case (a clear gain on one problem and a slight one on the
other) and is worth turning off if its diagnostics look poor.

The two survival rows are the other end of the range from the negative
binomial, and for a reason worth stating: right-censoring is what makes
the direct likelihood expensive. An observed failure contributes a
density in the predictor and a censored observation contributes a
survival function, and the two have different shapes, so the target has
no exploitable form and every trial value of a leaf parameter costs its
own pass over the data. Imputing each censored failure time above its
censoring time replaces the survival term with a density, and then every
observation contributes the same quadratic shape.
[`weibull_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is absent from the table because it needs none of this: its likelihood
is already of the exponential form, censoring included, so it gets the
single-pass treatment directly under hard rules.

The binomial family uses Polya-Gamma augmentation for a logit link and
Albert and Chib's latent normal for a probit one. The negative binomial
uses neither: it is written as a Poisson whose rate is drawn from a
gamma, which costs one gamma draw per observation and leaves the target
in the exponential form the sampler can collapse to a single pass, but
only under hard rules, which is why the gain is there and not under soft
ones.

`ordinal("probit")` is the largest gain of any of them, because the
target it replaces is the most expensive: two cumulative-normal
evaluations per observation per pass, and a target with no exploitable
shape, so every value of a leaf parameter costs its own pass.
Conditional on the latent normal all of that collapses to one pass of
arithmetic. The cutpoints are *not* drawn from their conditional given
the latent variables, which would be uniform between the two order
statistics bracketing each one and would mix worse and worse as the
sample grows; they are drawn from the ordinal likelihood with the latent
variables integrated out, and the latent variables are redrawn
immediately afterwards. That is a partially collapsed Gibbs sampler (Van
Dyk and Park 2008) and is the standard remedy (Cowles 1996).

`ordinal("logit")` gets the same treatment by a different route, because
its latent variable is logistic rather than normal. A logistic variate
is a normal whose precision is itself random, and the mixing
distribution is usually reached through the Kolmogorov-Smirnov density,
which needs a sampler of its own. It does not have to be: Polson, Scott
and Windle's (2013) Theorem 1 at `a = 1`, `b = 2` says the standard
logistic density is `(1/4) E[exp(-w x^2 / 2)]` with `w` Polya-Gamma(2,
0), so the conditional of the precision given a residual `r` is exactly
Polya-Gamma(2, \|r\|), an integer-parameter draw, which the exact
Devroye sampler already used for the logistic family covers. Nothing
approximate enters. It costs one extra draw per observation per sweep
relative to the probit, which is why its gain is about half as large.

`ordinal("cloglog")` has a latent variable of a different kind, and it
is the one the model is usually named for: the cumulative complementary
log-log model is the discrete proportional hazards model, so its
survivor function `exp(-exp(c_k - eta))` says exactly that a waiting
time with exponential rate `exp(-eta)` has passed `exp(c_k)`.
Conditional on that time the log density is `-eta - T * exp(-eta)`,
which is the *exponential* form rather than the quadratic one (i.e., the
same shape as the gamma family). That is worth 5.1x under hard rules,
where the exponential form applies; under soft rules, where it does not,
what is left is one [`exp()`](https://rdrr.io/r/base/Log.html) per
observation instead of a difference of two extreme-value distribution
functions, which is a smaller gain and costs more mixing.

### Rare Events Mix Slowly

With very few events the whole fit mixes slowly, augmented or not.
Measured at 2000 observations, 50 trees and 1000 draws, the effective
sample size of the level of the predictor falls from about 860 at half
positives to 75 at 1.4% and to 8 at 0.2%. That is worth knowing about
but it is not a fault of the rewritings: the *shape* of the fit mixes at
the same rate as its level, moving the anchor by three units on the
probit scale does not change the figure, and dbarts (an independent
implementation of the same latent normal, with the same anchoring)
reproduces it to three digits. It is the information in a handful of
events. Lengthen the chain, and read the `rhat` element of the fit
rather than assuming the default draw count is enough.

[`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
uses the same Polya-Gamma identity as the binomial logit, and can,
because conditional on the other categories the likelihood of category
*j* is *exactly* binomial-logistic in `eta_j - log C_j`, where `C_j`
sums `exp(eta)` over the others. So no stick-breaking decomposition is
needed: the structure the sampler already has, one forest at a time, is
the structure the augmentation wants. The latent variables are redrawn
per forest rather than per sweep, since `C_j` moves during a sweep.

The **zero-inflated** families need two latent variables, because what
blocks them is a mixture rather than a link. The zero contributes
`log[pi + (1 - pi) P_0]`, a log-sum-exp of the two components, so
neither predictor has a shape. Introducing the indicator the mixture is
taken over (i.e., whether the observation is a structural zero)
separates them: conditional on it the count forest sees a plain Poisson
or negative binomial, and the inflation forest sees a Bernoulli logistic
likelihood, which Polya-Gamma handles. The indicator is drawn from its
exact conditional, which is zero whenever the count is positive. For
[`zi_negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
the Poisson-gamma rate goes on top of that, so the count forest gets the
exponential form as well; the indicator is drawn with the rate
integrated out and the rate redrawn afterwards, which is a valid
partially collapsed Gibbs step in that order.

### Progress

`verbose = TRUE` prints a line to the console every `num_print`
iterations, which is the whole of what this package decides about
progress. A progress bar is
[progressr](https://CRAN.R-project.org/package=progressr)'s business,
and the sampler reports to it unconditionally: nothing is shown unless a
handler is active, so there is no argument to switch on.

    progressr::with_progress(
      bartisan(y ~ ., data = d, family = gaussian())
    )

    # or once, for the session
    progressr::handlers(global = TRUE)

The bar is sized for the whole fit, so `chains = 4` fills one bar once
rather than four in sequence, and progress from chains running in
parallel under [future](https://CRAN.R-project.org/package=future) is
relayed back as it arrives. Warmup and sampling are one run for this
purpose, because they cost the same per iteration and a bar that
restarts halfway is a worse report than one that does not.

The bar covers the sampling, which is what a fit spends its time on.
Convergence diagnostics are not part of it because they are not part of
a fit: computing R-hat and the two effective sample sizes for every
observation cost about as much as all four chains' sampling at 2000
observations, so they run in
[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
instead, which reports its own progress the same way.

Reporting is capped at 50 steps per chain rather than one per iteration.
The callback itself is nothing next to a sweep, but a handler that
redraws a bar is not, and it would otherwise be possible for the
reporting to cost more than the sampling. Progress does not touch the
draws: the same seed gives the same fit whether or not anything was
listening.

### Soft Rules and the Cost of a Gate

A soft rule is what makes the fit smooth, and it is charged for in two
places. Every observation reaches more than one leaf, so a pass over a
node covers more than the sample (measured at 2.5 times, at the default
bandwidth, for the logistic gate). And the bandwidth is itself a
parameter with a Metropolis step per tree per sweep, each of which
rebuilds every membership weight in the tree, twice when it is rejected.
On a Gaussian response with a thousand observations that one move is
about half the total time.

A bounded gate addresses the first: past its half-width from the
cutpoint the gate is exactly zero or one, so the observation takes one
side outright, the other subtree is never visited, and the gate is a
polynomial rather than an [`exp()`](https://rdrr.io/r/base/Log.html).
Measured at 1.4x with no loss of accuracy, on three test functions and
six replicates. The second is cheaper than it was, since a rejected
proposal is rolled back from a snapshot rather than by evaluating every
gate again, and a tree with no splits has no gate at all, so its
bandwidth is drawn straight from its prior, which is exactly its full
conditional; it is still about half of a soft-rule fit.

The half-width is `pi * sqrt((2a + 1) / 3)` times `bandwidth` for the
Beta(a, a) gate (4.06 for `"smoothstep"` and 4.80 for `"smootherstep"`),
which equates the gates' standard deviations, so `bandwidth` means the
same amount of smoothing whichever is chosen.

*Which* bounded gate makes almost no difference, and not for the reason
it looks like it should. At a bandwidth wide enough that a bounded gate
truncates nothing at all it is still 1.45 times faster than the
logistic: what the bounded gates save is the
[`exp()`](https://rdrr.io/r/base/Log.html), not the work on the far side
of the cutpoint. So the two come out within noise of each other on time
and on accuracy, and the reason to prefer `"smootherstep"` is that it
gives a twice-differentiable fit.

`bandwidth_every` addresses the second, and is a real trade rather than
a free one. Fixing the bandwidth entirely (`update_bandwidth = FALSE`)
is 2.4x faster again and *more* accurate on smooth functions, but much
worse on functions with jumps, where what the update is for is letting
the rules sharpen toward hard ones. On a three-step function, measured
RMSE was .42 fixed against .19 drawn. So the default draws it every
sweep, and a larger `bandwidth_every` buys time at the cost of how fast
that adaptation happens.

There is nothing to rewrite for the Poisson and gamma families, whose
targets are already in the exponential form, and no known rewriting for
the accelerated failure time, ordered beta, or location-scale families.

## References

Albert, J. H., & Chib, S. (1993). Bayesian analysis of binary and
polychotomous response data. *Journal of the American Statistical
Association*, 88(422), 669–679.

Cowles, M. K. (1996). Accelerating Monte Carlo Markov chain convergence
for cumulative-link generalized linear models. *Statistics and
Computing*, 6(2), 101–111.

Linero, A. R. (2018). Bayesian regression trees for high-dimensional
prediction and variable selection. *Journal of the American Statistical
Association*, 113(522), 626–636.
[doi:10.1080/01621459.2016.1264957](https://doi.org/10.1080/01621459.2016.1264957)

Linero, A. R. (2025). Generalized Bayesian additive regression trees
models: beyond conditional conjugacy. *Journal of the American
Statistical Association*, 120(549), 356–369.

Linero, A. R., & Yang, Y. (2018). Bayesian regression tree ensembles
that adapt to smoothness and sparsity. *Journal of the Royal Statistical
Society Series B*, 80(5), 1087–1110.

Linero, A. R., Sinha, D., & Lipsitz, S. R. (2020). Semiparametric
mixed-scale models using shared Bayesian forests. *Biometrics*, 76(1),
131–144. [doi:10.1111/biom.13107](https://doi.org/10.1111/biom.13107)

Van Dyk, D. A., & Park, T. (2008). Partially collapsed Gibbs samplers:
theory and methods. *Journal of the American Statistical Association*,
103(482), 790–796.

## See also

[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md),
which takes the result as its `control` argument;
[bartisan-families](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
for the forest names the per-forest arguments are keyed by

## Examples

``` r
data("rhc")

# Settings can be built up once and reused across fits
ctrl <- bartisan_control(num_trees = 20, gate = "hard", num_burn = 50,
                         num_draws = 50)

fit <- bartisan(death ~ . - days, data = rhc, control = ctrl)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` to choose another, which also silences this message.

# The same call, with the settings passed through `...` instead
fit2 <- bartisan(death ~ . - days, data = rhc, num_trees = 20,
                 gate = "hard", num_burn = 50, num_draws = 50)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` to choose another, which also silences this message.

# A setting given once applies to every forest, and a vector gives each
# forest its own value. A variance surface needs less capacity than a
# mean surface
bartisan_control(num_trees = c(mean = 50, log_sd = 10))
#> $num_trees
#>   mean log_sd 
#>     50     10 
#> 
#> $gate
#> [1] "smoothstep"
#> 
#> $soft
#> [1] TRUE
#> 
#> $sparsity
#> [1] TRUE
#> 
#> $share_sparsity
#> [1] FALSE
#> 
#> $share_forests
#> [1] FALSE
#> 
#> $split_prior
#> NULL
#> 
#> $categorical
#> [1] "subset"
#> 
#> $k
#> [1] 2
#> 
#> $bandwidth
#> [1] 0.1
#> 
#> $chains
#> [1] 1
#> 
#> $num_burn
#> [1] 200
#> 
#> $num_draws
#> [1] 800
#> 
#> $num_thin
#> [1] 1
#> 
#> $augment
#> [1] "binomial"    "ordinal"     "multinomial" "zip"         "zinb"       
#> [6] "aft"        
#> 
#> $x_transform
#> [1] "quantile"
#> 
#> $gamma
#> [1] 0.95
#> 
#> $beta
#> [1] 2
#> 
#> $sigma_mu
#> NULL
#> 
#> $update_sigma_mu
#> [1] TRUE
#> 
#> $sigma_mu_ramp
#> [1] 0.25
#> 
#> $update_tau
#> [1] TRUE
#> 
#> $update_bandwidth
#> [1] TRUE
#> 
#> $bandwidth_every
#> [1] 1
#> 
#> $alpha
#> [1] 1
#> 
#> $alpha_scale
#> [1] 0
#> 
#> $alpha_shape_1
#> [1] 0.5
#> 
#> $alpha_shape_2
#> [1] 1
#> 
#> $update_s
#> [1] TRUE
#> 
#> $update_alpha
#> [1] TRUE
#> 
#> $verbose
#> [1] FALSE
#> 
#> $num_print
#> [1] 100
#> 
#> $block_eval
#> [1] FALSE
#> 
#> $exact_quadratic
#> [1] TRUE
#> 
#> $generic_accumulate
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "bartisan_control"
#> attr(,"supplied")
#> attr(,"supplied")$num_trees
#>   mean log_sd 
#>     50     10 
#> 

# Weighting the splitting prior toward the treatment, which the sparsity
# prior would otherwise be free to drop
bartisan_control(split_prior = c(rhc = 10))
#> $num_trees
#> NULL
#> 
#> $gate
#> [1] "smoothstep"
#> 
#> $soft
#> [1] TRUE
#> 
#> $sparsity
#> [1] FALSE
#> 
#> $share_sparsity
#> [1] FALSE
#> 
#> $share_forests
#> [1] FALSE
#> 
#> $split_prior
#> rhc 
#>  10 
#> 
#> $categorical
#> [1] "subset"
#> 
#> $k
#> [1] 2
#> 
#> $bandwidth
#> [1] 0.1
#> 
#> $chains
#> [1] 1
#> 
#> $num_burn
#> [1] 200
#> 
#> $num_draws
#> [1] 800
#> 
#> $num_thin
#> [1] 1
#> 
#> $augment
#> [1] "binomial"    "ordinal"     "multinomial" "zip"         "zinb"       
#> [6] "aft"        
#> 
#> $x_transform
#> [1] "quantile"
#> 
#> $gamma
#> [1] 0.95
#> 
#> $beta
#> [1] 2
#> 
#> $sigma_mu
#> NULL
#> 
#> $update_sigma_mu
#> [1] TRUE
#> 
#> $sigma_mu_ramp
#> [1] 0.25
#> 
#> $update_tau
#> [1] TRUE
#> 
#> $update_bandwidth
#> [1] TRUE
#> 
#> $bandwidth_every
#> [1] 1
#> 
#> $alpha
#> [1] 1
#> 
#> $alpha_scale
#> [1] 0
#> 
#> $alpha_shape_1
#> [1] 0.5
#> 
#> $alpha_shape_2
#> [1] 1
#> 
#> $update_s
#> [1] FALSE
#> 
#> $update_alpha
#> [1] FALSE
#> 
#> $verbose
#> [1] FALSE
#> 
#> $num_print
#> [1] 100
#> 
#> $block_eval
#> [1] FALSE
#> 
#> $exact_quadratic
#> [1] TRUE
#> 
#> $generic_accumulate
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "bartisan_control"
#> attr(,"supplied")
#> attr(,"supplied")$split_prior
#> rhc 
#>  10 
#> 
```
