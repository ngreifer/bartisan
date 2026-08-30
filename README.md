
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bartisan

Bayesian additive regression trees for response distributions that
standard BART implementations cannot reach, with an interface that looks
like `glm()`.

``` r
bartisan(y ~ x1 + x2, data = d, family = binomial("logit"))
```

## Why this exists

BART is normally tied to a Gaussian likelihood. The Gibbs sampler
integrates the leaf parameters out of the tree likelihood in closed
form, and that step needs conjugacy. Getting anything else — a logistic
model, a survival model — usually means finding a data-augmentation
scheme that turns the problem back into a Gaussian one, and for many
likelihoods no such scheme exists.

*bartisan* implements the approach of Linero (2025), which removes the
restriction. Rather than integrating the leaf parameters out, the
sampler builds a Gaussian approximation to their conditional posterior
by Fisher scoring, and uses that as the proposal in a reversible-jump
Metropolis step. The approximation only has to be good enough to be
accepted often; the chain targets the exact posterior regardless of how
good it is.

The practical consequence is that adding a family requires only the log
density of one observation and its first two derivatives with respect to
the additive predictor. That is why the list below includes models that
are otherwise hard to find in a BART package.

| Family | Compiled links | Nuisance parameters drawn |
|----|----|----|
| `gaussian()` | identity | residual standard deviation |
| `binomial()` | logit, probit, cloglog | — |
| `poisson()` | log | — |
| `negbin()` | log | dispersion |
| `Gamma("log")` | log | shape |
| `ordinal()` | logit, probit, cloglog | cutpoints |
| `multinomial()` | logit, probit | latent covariance, for probit |
| `weibull_aft()`, `loglogistic_aft()`, `lognormal_aft()` | — | scale |
| `ph()` | — | baseline hazard per bin |
| `location_scale()` | identity | — |
| `zi_poisson()` | log | — |
| `zi_negbin()` | log | dispersion |
| `Beta()` | logit, probit, cloglog | precision |
| `ordbeta()` | logit | 2 cutpoints, precision |
| `dpm()` | identity | error mixture, concentration |

`family` can also be left out, in which case it is read off the response
— a `Surv` object gives `dpm_aft()`, an ordered factor `ordinal()`, a
logical or a two-level factor `binomial()`, a factor with more levels
`multinomial()`, and any other numeric response `dpm()`, whatever its
number of distinct values. A message says which was chosen; naming
`family` silences it. A count is deliberately *not* read as `poisson()`,
and neither is a numeric response whose two values are not zero and one:
in both cases the guess would be a modeling decision rather than a
reading of the response’s type.

Ordinary `stats::family` objects work unchanged, so moving a model from
`glm()` to `bartisan()` is a one-word change. `multinomial()`,
`location_scale()` and the two zero-inflated families fit one forest per
unconstrained parameter, so the excess-zero mechanism of a zero-inflated
count model is itself nonparametric rather than a single constant.
`ordbeta()` is the ordered beta regression of Kubinec (2023), for a
proportion or slider response with mass piled up at zero and one.
`multinomial()` follows Murray (2021) in fitting one forest per category
and leaving the model unidentified, so that the prior does not depend on
which category was singled out as the reference.

`dpm()` is the default worth reaching for on a numeric response, and it
is not a distribution but a way of not choosing one. It is DPMBART
(George et al. 2019): the sum of trees for the mean, as `gaussian()`
has, and a Dirichlet process mixture of normals for the errors instead
of one normal, so the error distribution comes out as whatever shape the
data ask for. `error_density()` reports that shape with a pointwise
interval. Against `gaussian()` on held-out predictive log score it is
level when the errors are normal (−2097 against −2097), and ahead by 175
and 224 log points when they are heavy tailed or skewed — on the skewed
case it also cuts the interval width for the regression function from
1.59 to 0.99 while holding coverage. It is *not* the family for
heteroskedasticity: a mixture makes the error distribution flexible but
keeps it the same at every `x`, which is what `location_scale()` is for,
and which wins that case by 119 log points. Conditional on the mixture
the target is still quadratic, so a fit costs 1.3× a Gaussian one
against `location_scale()`’s 9×. It also does not pay for its
flexibility: on normal errors, where `gaussian()` is exactly right, it
came out slightly *ahead* on held-out error and log score, which is why
it is what a numeric response gets by default. The reason to prefer
`gaussian()` is that it takes prior weights, which `dpm()` refuses. The
additive predictor is the conditional mean for both: nothing forces the
mixture to be centered, so the sampler works in a chart where only their
sum is identified and reports in the one where the mixture has mean
zero.

`multinomial("probit")` is the same categorical model with *correlated*
latent utilities — the thing a multinomial logit cannot express. Its
covariance matrix is drawn and reported, normalized by the trace
constraint so that a two-category fit is exactly binary probit.
Conditional on the latent utilities the target is exactly quadratic, so
the leaf draw is a closed form and a fit costs about twice a multinomial
logit’s. Its likelihood is a Gaussian orthant probability with no closed
form, so probabilities are simulated. Its latent correlations are weakly
identified, so read the fitted probabilities rather than the covariance.

`vignette("families")` goes through every family: what each is and is
not for, how to choose among the links a family offers, how to choose
among the zero-inflated models, and what a family inferred from the
response will be if you name none. `vignette("survival")` takes the five
censored-response families on their own: what each one is estimating,
which shapes of hazard each can and cannot represent, how they compare
across six data-generating truths and across censoring from none to 70%,
and what to do when hazards are not proportional.

## Links and likelihoods of your own

The links in the table are the ones the sampler evaluates in compiled
code. Any other link works for the five families with a single mean and
a conventional link, applied from R by composing it onto the scale the
compiled family works on, so `binomial("cauchit")` needs no new syntax
and neither does a link you write yourself and pass the way `glm()`
takes one.

Because the whole interface a family has to satisfy is a log density and
two derivatives, it can also be supplied from R outright:

``` r
# A Poisson model written by hand. Terms free of eta may be dropped.
bartisan(y ~ ., data = d,
         family = custom_family(function(y, eta) y * eta[, 1] - exp(eta[, 1])))

# Two additive predictors: a mean and a log standard deviation.
bartisan(y ~ ., data = d,
         family = custom_family(
          function(y, eta) dnorm(y, eta[, 1], exp(eta[, 2]), log = TRUE),
          num_predictors = 2))
```

The function is called once per leaf, not once per observation, so it
must be vectorized over its arguments. Central differences supply the
derivatives unless you pass them; a custom family cannot draw a nuisance
parameter.

## Missing predictor values

A missing predictor is **kept by default**, not dropped and not imputed.
Each splitting rule carries the answer for itself: a rule on a variable
that has any missing values is drawn as one of three, with equal
probability, sending them left with the low values, right with the high
ones, or splitting on whether the value is there at all. Pass
`na.action = na.omit` for the usual behavior of dropping the row.

``` r
bartisan(y ~ ., data = d)                          # missing predictors kept
bartisan(y ~ ., data = d, na.action = na.omit)     # rows dropped instead
```

A missing response, weight or offset is dropped either way, with a
warning, since there is nothing to fit those rows to.

This is *missingness incorporated in attributes* (Twala et al. 2008; for
BART, Kapelner and Bleich 2015). The third rule is the one that matters:
it lets the model use a variable whose **absence** is the signal, which
no imputation scheme can recover. On a problem where the response
depends on whether a predictor was recorded and on nothing else about
it, the fitted means came out at 1.98 and -0.01 against a truth of 2 and
0, with 29.6 of the 30 splitting rules landing on that variable.

Because the choice is drawn from its prior alongside the variable and
the cutpoint, it cancels from every acceptance ratio, and a variable
with no missing values is not given the extra draw at all – complete
data reproduces the sampler exactly as it was. What this estimates is
the mean of the response given the predictors *and* the pattern of
missingness, which is what you want for prediction; if the estimand is
defined on complete data, multiple imputation is the right tool instead.

## Smoother fits

By default the decision rules are *soft*, following Linero and Yang
(2018): a rule is a smooth gate rather than a step, so an observation
reaches every leaf with some weight and the fitted function is smooth
rather than piecewise constant. The default gate is the bounded
`"smoothstep"`, which is about 1.4 times faster than Linero and Yang’s
logistic one and equally accurate; `gate` in `bartisan_control()`
chooses. Combining soft rules with a non-conjugate likelihood is an
extension of Linero (2025), which leaves it as an open problem — soft
rules make the leaf parameters of a tree dependent on one another, so
the reversible-jump move needs a bivariate Laplace proposal for the pair
of child leaves. That proposal reduces exactly to Linero’s independent
pair when the rules are hard, so the two cases share one implementation.

Set `gate = "hard"` in `bartisan_control()` for the faster hard-rule
sampler.

## Settings worth knowing about

`bartisan_control()` groups its arguments so that the ones that matter
come first. Three are worth a look before any fit.

`gate` is the one argument that chooses between hard and soft rules and,
when soft, the shape of the gate. `sparsity` turns the Dirichlet
variable-selection prior of Linero (2018) on and off, or names one of
four strengths, in place of the four hyperparameters that actually
parameterize it. And `num_trees` takes one value per additive predictor,
not just one for the fit:

``` r
bartisan(y ~ ., data = d, family = location_scale(), num_trees = c(50, 10))
```

That last one is worth using. `location_scale()` spends about 90% of its
time on the scale forest, and a variance surface needs less capacity
than a mean surface, so a smaller second forest cuts the fit from ten
times a Gaussian one to three and a half at the same accuracy.
`?bartisan_control` has the measurements, along with the tree-count
curves for the two kinds of rule — soft rules peak at around 20 trees
and get worse past that, hard rules keep improving to 200 — and what the
sparsity prior costs in exchange for dropping predictors.

## Example

``` r
library(bartisan)

set.seed(1)

n <- 400
d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n), x4 = runif(n))
d$y <- rbinom(n, 1, plogis(3 * sin(pi * d$x1 * d$x2) - 1))

fit <- bartisan(y ~ ., data = d, family = binomial())

fit
```

    ## Generalized BART
    ## 
    ## Call: `bartisan(formula = y ~ ., data = d, family = binomial())`
    ## 
    ## Family: "binomial" with the "logit" link
    ## Observations: 400
    ## Structure: 1 forest of 50 trees, soft decision rules
    ## Draws: 500 kept after 500 warmup

``` r
summary(fit)          # includes how often each predictor was split on
```

    ## Generalized BART
    ## 
    ## Call: `bartisan(formula = y ~ ., data = d, family = binomial())`
    ## 
    ## Family: "binomial" with the "logit" link
    ## Observations: 400
    ## Structure: 1 forest of 50 trees, soft decision rules
    ## Draws: 500
    ## 
    ## Predictor usage
    ## Splitting rules per draw, and how often used at all.
    ##      mean     sd lower upper prop_used
    ## x1 39.410 12.496    19    65     1.000
    ## x2 32.164 12.663     9    52     1.000
    ## x4  2.488  3.047     0    10     0.610
    ## x3  2.156  2.823     0    10     0.602

``` r
head(predict(fit, type = "response"))
```

    ## [1] 0.6334514 0.4320595 0.7682428 0.8133213 0.5981106 0.8445637

``` r
predict(fit, newdata = d[1:5, ], type = "prob")
```

    ##              0         1
    ## [1,] 0.3665486 0.6334514
    ## [2,] 0.5679405 0.4320595
    ## [3,] 0.2317572 0.7682428
    ## [4,] 0.1866787 0.8133213
    ## [5,] 0.4018894 0.5981106

Every draw of every tree is retained, so `predict()` can return the full
posterior rather than a point estimate:

``` r
draws <- predict(fit, newdata = d[1:5, ], type = "response", draws = TRUE)
apply(draws, 2, quantile, c(0.025, 0.975))
```

    ##            [,1]      [,2]      [,3]      [,4]      [,5]
    ## 2.5%  0.5063507 0.2833474 0.6221387 0.6946875 0.4310108
    ## 97.5% 0.7655796 0.5800427 0.8922670 0.9007625 0.7366930

`type = "density"` evaluates the conditional density of the outcome
given the predictors, at the outcome’s own value, so the outcome has to
be present in `newdata`. Summed on the log scale that is a held-out log
score, which is a fair basis for comparing models:

``` r
held_out <- data.frame(x1 = runif(50), x2 = runif(50), x3 = runif(50),
                       x4 = runif(50))
held_out$y <- rbinom(50, 1, plogis(3 * sin(pi * held_out$x1 * held_out$x2) - 1))

sum(predict(fit, newdata = held_out, type = "density", log = TRUE))
```

    ## [1] -31.79423

A right-censored survival model, with the response given as a `Surv`
object:

``` r
bartisan(survival::Surv(time, status) ~ ., data = d, family = weibull_aft())
```

A zero-inflated count model, where both the count mean and the
structural-zero probability get their own forest:

``` r
bartisan(count ~ ., data = d, family = zi_negbin())
```

## Random intercepts

A `(1 | group)` term adds an intercept per level of `group`, drawn from
a common mean-zero normal whose standard deviation is itself drawn:

``` r
bartisan(y ~ x1 + x2 + (1 | school), data = d)
bartisan(y ~ x1 + (1 | school) + (1 | year), data = d)   # several factors
bartisan(y ~ x1 + (1 | district / school), data = d)     # nesting
```

Each additive predictor gets its own set, so a zero-inflated count model
has a group effect on the count part and another on the inflation part.
Only random *intercepts* are supported; a slope is refused rather than
ignored, because a variable whose effect varies by group belongs in the
fixed part of the formula, where a tree can split on the group and the
variable together and get an interaction of any shape.

Whether to reach for this at all depends on the shape of the data, and
the answer is not “always”. A grouping factor can go in the fixed part
instead, where a tree splits on it like anything else. Measured against
`dbarts::rbart_vi()` on data with a true group effect, held-out RMSE:

| groups | per group | group as a predictor | random intercept |
|--------|-----------|----------------------|------------------|
| 5      | 100       | **0.190**            | 0.279            |
| 25     | 20        | **0.306**            | 0.341            |
| 100    | 5         | 0.541                | **0.468**        |
| 250    | 4         | 0.708                | **0.493**        |

With few large groups the predictor route is better — the group means
are well determined without pooling, and a split can interact the group
with the covariates. The random intercept wins once there are many small
groups, which is where partial pooling earns its keep.

## Counterfactual estimands

A fit works with **marginaleffects**, so predictions, comparisons and
slopes — and hypothesis tests on any of them — come out without
extracting draws by hand:

``` r
library(marginaleffects)

avg_predictions(fit)
avg_comparisons(fit)                        # average effect of each predictor
avg_predictions(fit, by = "group")          # by a factor
comparisons(fit, variables = "treatment",   # at a grid of covariate values
            newdata = datagrid(age = c(30, 50, 70)))
```

A forest has no coefficient vector and no variance-covariance matrix, so
there is no delta method to apply. There is something better: every
estimand is computed by pushing all of the posterior draws through the
same transformation, so the intervals are posterior quantiles rather
than a normal approximation, and a nonlinear estimand needs no
approximation at all.

One caveat, documented in `?bartisan-marginaleffects`. `avg_slopes()` is
a numerical derivative, and the default `x_transform = "quantile"` maps
each predictor through its empirical distribution function — a step
function — so the fit is a step function of the original predictor
whatever the decision rules are. Use `x_transform = "range"` if slopes
are the estimand. Predictions and comparisons evaluate the fit at two
points a substantive distance apart and are unaffected.

## Posterior predictive checks and fit statistics

Replicate outcomes can be drawn from the fitted model, which is what the
packages that assess a Bayesian fit want:

``` r
yrep <- posterior_predict(fit)          # draws by observations
ll   <- log_lik(fit)                    # pointwise log likelihood

loo::loo(ll)                            # leave-one-out, with its diagnostics
bayesplot::pp_check(fit, "dens_overlay")
performance::model_performance(fit)     # ELPD, LOOIC, WAIC, R2, RMSE, sigma
performance::r2(fit)                    # Bayesian R2, with a posterior interval
posterior::summarise_draws(posterior::as_draws(fit))
```

`fitted()`, `residuals()`, `weights()`, `sigma()` and `simulate()`
behave as they do for a `glm`, and `insight` reads the fit through them,
so the **easystats** packages see it too. `?bartisan-interop` documents
the scale each replicate comes back on — a binomial replicate is a
proportion, a categorical one a category index, an
accelerated-failure-time one an uncensored event time — and why a
leave-one-out estimate is strained for a model as flexible as a forest.

## Variable selection

Splitting variables are drawn from a Dirichlet prior that concentrates
on few predictors, following Linero (2018), so the model filters out
irrelevant ones. `summary()` reports the posterior distribution of the
number of splitting rules using each predictor and the proportion of
draws using it at all. The indicator columns of a factor share one
weight, so a factor is selected as a whole rather than one level at a
time.

## Correctness

Three checks are in the test suite and worth knowing about.

Shrinking the prior weights to nothing makes the likelihood constant, so
the posterior collapses to the tree prior. The sampled trees then
reproduce the prior’s branching process to within a fraction of a
percent, which is a direct check that the birth, death and change moves
satisfy detailed balance.

Every family’s analytic score is checked against a central difference of
its own log density, which is an independent route to the same quantity.
This is the check that matters for the sampler: a wrong score sends the
Fisher-scoring proposal to the wrong place.

Across replicate datasets at the default settings, 95% credible
intervals for the additive predictor cover the truth 0.95 (Gaussian),
0.91 (binomial), 0.96 (Poisson) and 0.96 (gamma) of the time. The
binomial case is the weakest because a binary response carries the least
information per observation. The shortfall tracks the ratio of absolute
bias to posterior standard deviation, which sits near 0.8: an interval
centered on a shrunken estimate loses coverage in proportion to how far
that center sits from the truth. Running more or longer chains does not
help, which rules out mixing as the cause at this scale. This is in line
with what is reported for BART generally; treat pointwise intervals for
the regression function as approximate.

## Speed

A non-conjugate sampler costs more than a conjugate one, and it is worth
being precise about how much – and about what the extra buys. From
`_dev/benchmark.Rmd`, which you can run yourself: the Friedman function,
n = 1000, 10 predictors, 50 trees, 1000 warmup and 1000 saved draws,
scored against the true regression function on a held-out thousand.

| Task | Package | Time | Effective sample size | Held-out RMSE |
|----|----|----|----|----|
| Gaussian | `dbarts` | 0.24 s | 20 | 0.221 |
|  | **bartisan, hard rules** | **0.44 s** | 18 | 0.219 |
|  | `stochtree` | 0.70 s | 18 | 0.218 |
|  | `bartMachine` | 0.87 s | — | 0.230 |
|  | `BART::wbart` | 1.06 s | 21 | 0.228 |
|  | bartisan, soft, smoothstep gate (default) | 1.40 s | 42 | **0.135** |
|  | bartisan, soft, smootherstep gate | 1.42 s | 43 | 0.140 |
|  | bartisan, soft, logistic gate | 2.00 s | 44 | 0.145 |
| probit | `dbarts` | 0.27 s | 36 | 0.134 |
|  | **bartisan, hard rules** | **0.49 s** | 35 | 0.122 |
|  | `stochtree` | 0.99 s | 29 | 0.134 |
|  | `BART::pbart` | 1.13 s | 36 | 0.125 |
|  | bartisan, soft, logistic gate | 2.05 s | 33 | **0.111** |
|  | bartisan, `augment = FALSE` | 23.6 s | 65 | 0.106 |
| logit | bartisan, soft, logistic gate | 2.11 s | 129 | **0.088** |
|  | bartisan, `augment = FALSE` | 13.9 s | 104 | 0.085 |
|  | `BART::lbart` | 20.8 s | 38 | 0.109 |
| ordinal | bartisan, hard, probit | 0.87 s | 35 | — |
|  | bartisan, hard, logit | 0.93 s | 31 | — |
|  | bartisan, soft, probit | 2.38 s | 51 | — |
|  | bartisan, soft, logit | 2.52 s | 48 | — |
|  | **bartisan, hard, cloglog** | **3.03 s** | 24 | — |
|  | `stochtree` (cloglog) | 4.01 s | 20 | — |
|  | bartisan, hard, logit, `augment = FALSE` | 14.9 s | 42 | — |
|  | bartisan, hard, cloglog, `augment = FALSE` | 17.0 s | 37 | — |
|  | bartisan, hard, probit, `augment = FALSE` | 25.2 s | 36 | — |

With hard rules bartisan is **within a factor of 1.8 of `dbarts`** on
the two tasks `dbarts` supports — 0.426 s against 0.241 s on the
Gaussian task and 0.482 s against 0.272 s on probit, best of six runs
each — at the same mixing and better accuracy, and it is faster *and*
more accurate than every other package here on both. On the
complementary log-log ordinal model, the one task `stochtree` supports
and `dbarts` does not, bartisan is now the faster of the two.

The other rows come from `_dev/benchmark.Rmd` at two replicates, which
is noisy at the ten to thirty percent level; the two ratios above are
quoted from a longer run because a claim about a factor of two should
not rest on a best-of-two. On a logit link it is 9.9 times faster than
`BART::lbart` and mixes three times better. Soft rules – the default –
cost a further four to five times and cut the held-out error by 30 to
40%, which makes the default configuration the most accurate fit in the
table; with a bounded gate that falls to about three times.

The ordinal rows have no RMSE column because the links put the additive
predictor on different scales, so the numbers would compare scales
rather than fits; the timings and the mixing are comparable.

The remaining families are the expensive ones, at 2 to 8 seconds, which
is the price of the general machinery where no rewriting is available.

Both bartisan figures use `augment = TRUE`. The gap decomposes into
factors that were measured rather than guessed, and profiling the
sampler by phase says where each one lives.

With **hard rules**, what is left against `dbarts` is tree bookkeeping,
not the family. A birth, death or change move makes five or six passes
over the node it touches – materializing the predictor with the node’s
own contribution removed, fitting the proposal, splitting the node’s
membership, recomputing the child weights, committing – against about
two for a sampler that only ever needs a sum of residuals and a count.
Nothing in that gap is the likelihood: the acceptance ratio bartisan
forms for a Gaussian response is algebraically the same
marginal-likelihood ratio `dbarts` forms, because a quadratic target
makes the Laplace approximation exact.

What has actually paid, in order, is worth knowing before reaching for
the obvious ideas. **Sizing the children’s index vectors and filling
them by index, instead of `push_back`,** is worth 20% of a soft fit on
its own – four vectors’ worth of per-element capacity tests in the
innermost loop of the sampler. **Recycling tree nodes** rather than
allocating two per birth proposal, about two thirds of which are
rejected, is worth 5%. **Folding the child weights into the split that
produces them** removes a whole pass over the same gates.

Two things that looked like they should help did not. A hard tree used
to carry a membership weight per observation per node, every one exactly
1.0; and the leaf-level sums used to write a node’s predictors into a
buffer and read them back. Removing both is bit-identical – multiplying
by 1.0 is exact – and together they are worth about 10%, against a
predicted factor of several. Memory bandwidth was never the constraint:
a node holds a few thousand observations at most, so the buffer stays in
L1. Nor is dropping the remaining materialized base array worth much:
measured directly, at about 2%, and it would trade a sequential read for
a gather.

With **soft rules**, half the cost is a single move. The bandwidth of
each tree is a parameter with a Metropolis step per tree per sweep, and
every attempt rebuilds every membership weight in the tree. A rejected
attempt – 58% of them – used to rebuild a second time to get back where
it started; it is now rolled back from a snapshot, and a tree with no
splits skips the move entirely and draws its bandwidth from the prior,
which is its exact full conditional. The other half is that a logistic
gate never saturates: dropping a weight needs it below 1e-10, which
needs the observation 23 bandwidths from the cutpoint, further than the
whole unit interval. So every observation reaches every leaf and a pass
over a node covers 2.5 times the sample. `gate = "smoothstep"`, the
default, and `gate = "smootherstep"` are the fix for the second half and
are in the table above – though not for the reason they were built. At a
bandwidth wide enough that a bounded gate truncates *nothing*, it is
still 1.45 times faster than the logistic: what a bounded gate saves is
the `exp()`, not the work on the far side of the cutpoint. So the two
bounded gates come out within noise of each other and the choice between
them is about smoothness – one derivative against two – rather than
speed. Fixing the bandwidth outright is the fix for the first, and is
*not* the default, because on a function with jumps it more than doubles
the error – letting the rules sharpen towards hard ones is what the move
is for.

The passes over the data are no longer part of it. A leaf value enters
the additive predictor linearly, so the log target over a leaf inherits
the shape of the log density – and where that shape is known, one pass
over a node determines the whole function and everything after it is
arithmetic. Two shapes are known. A *quadratic* log density (a Gaussian
response, and anything `augment` rewrites into one) makes the Laplace
approximation the conditional posterior exactly. An *exponential* one,
`a*eta + b*exp(r*eta)`, covers the Poisson, the gamma, the Weibull
survival model, the augmented negative binomial and `location_scale()`’s
log-scale forest under hard rules; there the fit still has to be
iterated, but on three numbers rather than on the data. The rate `r`
matters: it is what brings in the last two, at `-1/sigma` and `-2`.
Either way a birth move goes from six passes over the node to two.

What closes most of the first factor is a data augmentation that makes
the conditional Gaussian, since then the Laplace approximation is exact
rather than approximate. `augment` does this: Albert and Chib (1993) for
a probit link and for an ordinal probit, Pólya-Gamma augmentation
(Polson, Scott and Windle 2013) for a logit one, an ordinal logit, or a
multinomial, and imputation of the censored failure times for the
log-normal and log-logistic survival models. Every one of them trades
speed for mixing, so the ratio to judge is effective samples per second,
and it differs a lot by family:

| Family                          | Speed | Effective sample size | ESS per second |
|---------------------------------|-------|-----------------------|----------------|
| `lognormal_aft()`, hard rules   | 29x   | 0.87x                 | **21x**        |
| `ordinal("probit")`, hard rules | 30x   | 0.76x                 | **23x**        |
| `lognormal_aft()`, soft rules   | 20x   | 1.05x                 | **15x**        |
| `loglogistic_aft()`, hard rules | 12x   | 0.88x                 | **10x**        |
| `loglogistic_aft()`, soft rules | 8.9x  | 0.71x                 | **6.3x**       |
| `ordinal("logit")`, hard rules  | 15x   | 0.94x                 | **14x**        |
| `ordinal("probit")`, soft rules | 14x   | 0.87x                 | **12x**        |
| `ordinal("logit")`, soft rules  | 7.1x  | 0.76x                 | **5.4x**       |
| `binomial("probit")`            | 5.7x  | 0.66x                 | **3.8x**       |
| `binomial("logit")`             | 2.6x  | 0.81x                 | **2.1x**       |
| `multinomial()`, soft rules     | 9.3x  | 1.09x                 | **10.1x**      |
| `multinomial()`, hard rules     | 14.5x | 0.66x                 | **9.6x**       |
| `zi_poisson()`, hard rules      | 7.1x  | 1.42x                 | **10.1x**      |
| `zi_negbin()`, hard rules       | 9.4x  | 0.84x                 | **7.9x**       |
| `zi_poisson()`, soft rules      | 4.6x  | 0.85x                 | **3.9x**       |
| `negbin()`                      | 0.8x  | 0.56x                 | 0.5x           |

`augment` is on by default and covers the binomial, ordinal,
multinomial, zero-inflated and survival families, plus the negative
binomial when the rules are hard – there it is written as a Poisson
whose rate comes from a gamma, which needs one gamma draw per
observation rather than a Pólya-Gamma one.

The survival families are the clearest case of what the rewriting is
for. Right-censoring is what makes their likelihood expensive: an
observed failure contributes a density in the additive predictor and a
censored one contributes a survival function, and the two have different
shapes, so nothing about the target is exploitable. Imputing each
censored failure time above its censoring time replaces the survival
term with a density, and then every observation contributes the same
quadratic shape. `weibull_aft()` is absent from the table because it
needs none of this – its likelihood is already exponential in the sense
above, censoring included – which also makes it the slowest of the three
at the default soft gate, where that form cannot be used.

The zero-inflated families are the odd ones out, because what blocks
them is a mixture rather than a link: the zero contributes
`log[π + (1 − π)P₀]`, so neither the count predictor nor the inflation
predictor has a shape to exploit. Introducing the indicator of *which*
component produced each observation separates them — the count forest
then sees a plain Poisson or negative binomial, and the inflation forest
a Bernoulli logistic likelihood that Pólya-Gamma handles.

The ordinal probit is the largest of these because the target it
replaces is the most expensive one in the package: two cumulative-normal
evaluations per observation per pass, and no exploitable shape, so every
trial value of a leaf parameter costs its own pass. It comes with one
subtlety worth naming. The cutpoints are *not* drawn from their
conditional given the latent normals – uniform between the two order
statistics that bracket each one, which pins them to an interval of
width O(1/n) and mixes worse the more data there is. They are drawn from
the ordinal likelihood with the latent normals integrated out, and the
latent normals are redrawn immediately afterwards: a partially collapsed
Gibbs sampler (Van Dyk and Park 2008), which is the standard remedy
(Cowles 1996).

The cumulative **logit** has no latent normal – its latent variable is
logistic – and the usual route to one goes through the
Kolmogorov-Smirnov density, which needs a sampler of its own. It does
not have to. Polson, Scott and Windle’s Theorem 1 at *a* = 1, *b* = 2
says that the standard logistic density is `(1/4) E[exp(-w x²/2)]` with
*w* ~ PG(2, 0), so a logistic residual is a normal whose precision is
Pólya-Gamma – and because PG(*b*, *c*) is PG(*b*, 0) tilted by
`exp(-c²w/2)`, the precision’s conditional given a residual *r* is
exactly PG(2, \|*r*\|). That is an integer-parameter draw, which the
exact Devroye sampler already in the package covers, so the logit link
gets the same quadratic target with nothing approximate anywhere. The
Poisson and gamma families need no rewriting: their targets are already
in a shape the sampler collapses to a single pass.

One warning worth heeding: an unoptimized build of the compiled code is
5 to 20 times slower and looks no different. `devtools::load_all()` and
`devtools::install()` compile without optimization and leave the object
files for a later `R CMD INSTALL` to reuse. `bartisan()` warns once per
session when it detects one; `bartisan:::.bartisan_optimized()` reports
the state directly.

Independent chains run in parallel:

``` r
future::plan(future::multisession, workers = 4)

fit <- bartisan(y ~ ., data = d, chains = 4)
fit$rhat
```

The chain is the only parallel axis this sampler has – a sweep
conditions on the last one – and it is also what makes a convergence
diagnostic possible. Any `future` backend works, including mirai’s. One
`set.seed()` reproduces the whole run whatever the backend.

## Notes and limitations

Cox’s *partial* likelihood is not available here: it couples
observations through risk sets, so it does not decompose into a sum over
the observations reaching a leaf, which is what the leafwise Laplace
approximation requires. What `ph()` fits instead is the full likelihood
of the piecewise-exponential proportional hazards model, which does
decompose and which approaches the partial likelihood as the bins
shrink. Its `num_bins` is not a tuning knob: over a sixty-fold range the
estimates are flat, and what the bin count changes is the effective
number of parameters, which is what keeps `loo()` and `waic()` usable.

Soft rules cost more per iteration than hard ones, because a leaf
touches every observation rather than only those in its cell. Negligible
weights are pruned, which holds the gap to roughly a factor of three
rather than the factor of *L* it would otherwise be.

Absolute timings are machine- and load-dependent, so treat these as
indicative rather than exact; repeated runs on the same laptop varied by
about 20%. As an anchor, a Gaussian response with 10 predictors, 50 soft
trees and the default 1000 warmup plus 1000 saved draws takes roughly 6
seconds at n = 500 and roughly a minute at n = 5000, measured at steady
state on one core of an M-series Mac. Hard rules are about 3 times
faster, and doubling the tree count roughly doubles the cost.

The ratios are stabler than the absolute times. Cost relative to a
Gaussian fit at the same size is set by how expensive the family’s log
density and its derivatives are:

| Family                      | Relative cost                            |
|-----------------------------|------------------------------------------|
| `gaussian()`                | 1x                                       |
| `Gamma("log")`, `poisson()` | ~2x                                      |
| `negbin()`, `binomial()`    | 2–3x                                     |
| `ordinal()`                 | ~2x with a latent variable, ~25x without |
| `zi_poisson()`              | ~10x                                     |
| `ordbeta()`                 | ~25x                                     |

`ordbeta()` is the expensive one because every evaluation of its log
density needs two log-gamma calls that no amount of restructuring
removes. `ordinal()` spans a wide range because all three of its links
have a latent-variable representation that the sampler uses by default —
a normal for the probit, a normal with a Pólya-Gamma precision for the
logit, an exponential waiting time for the complementary log-log — and
`augment = FALSE` turns them off.

An ordinal model is identified only up to a common shift of its
cutpoints and its predictor. With three or more categories the draws are
reported in the chart where **the predictor has mean zero over the
fitted sample and every cutpoint is free**, which makes the cutpoints
readable as category boundaries. With two categories the single boundary
is folded into the intercept instead, so the fit is on the same scale as
`binomial()`.

That is `MASS::polr()`’s chart with its predictors centered. `polr()`
identifies the location by leaving the intercept out of the design
matrix rather than by centering, so its `zeta` is shifted by the mean of
its own linear predictor. Either of these lines puts the two side by
side, and on a linear truth they agree to Monte Carlo error:

``` r
p <- MASS::polr(y ~ x1 + x2, data = d, method = "probit")
p$zeta - mean(p$lp)                      # same chart as colMeans(fit$aux)
```

Two prediction types exist for reporting an ordinal fit on a single
scale, both following **WeightIt**. `predict(type = "mean")` weights the
category probabilities by the labels read as numbers, so levels `"1"`,
`"2"`, `"4"` give a mean between one and four; `values` says what the
categories are worth when the labels are not numbers.
`predict(type = "stdlv")` divides the predictor by the standard
deviation of the latent variable it indexes, which puts fits with
different links, or different amounts of signal, on one scale; it works
for `binomial()` too, since a binary response is the same construction
with one threshold. The complementary log-log error enters the two
families with opposite signs, and `?predict.bartisan` says why.

The leaf scale is raised from near zero over the first quarter of
warmup. Linero (2025) describes this as essential: started at its full
value the sampler can settle early into a poor configuration and fail to
move. It is controlled by `sigma_mu_ramp` in `bartisan_control()`.

Nuisance parameters that are weakly identified — the negative binomial
dispersion, the shape of a survival model — mix more slowly than the
regression function, and the data often do not pin them down. Check
their draws before interpreting them.

The leaf scale is drawn under a half-Cauchy prior, which has no upper
bound. When the predictors nearly separate a binary response the
likelihood rewards an unbounded predictor and that prior does not hold
the scale down: in a fully separated example the drawn scale wandered
between 3 and 9 times its prior median over 1600 draws without settling.
`bartisan()` warns when the scale settles more than five times above its
prior median; `update_sigma_mu = FALSE` pins it.

A link supplied from R costs a call into the interpreter for every leaf
the sampler visits, and `custom_family()` costs three unless you supply
the derivatives. Both are usable at the sizes a BART model is usually
fit at, and both are slower than the compiled equivalent.

## Installation

``` r
# install.packages("pak")
pak("ngreifer/bartisan")
```

Compilation requires a C++ toolchain, plus *Rcpp* and *RcppArmadillo*.

## References

Albert, J. H., & Chib, S. (1993). Bayesian analysis of binary and
polychotomous response data. *Journal of the American Statistical
Association*, 88(422), 669–679.

Kapelner, A., & Bleich, J. (2015). Prediction with missing data via
Bayesian additive regression trees. *Canadian Journal of Statistics*,
43(2), 224–239.

Polson, N. G., Scott, J. G., & Windle, J. (2013). Bayesian inference for
logistic models using Polya-Gamma latent variables. *Journal of the
American Statistical Association*, 108(504), 1339–1349.

Murray, J. S. (2021). Log-linear Bayesian additive regression trees for
multinomial logistic and count regression models. *Journal of the
American Statistical Association*, 116(534), 756–769.

Twala, B. E. T. H., Jones, M. C., & Hand, D. J. (2008). Good methods for
coping with missing data in decision trees. *Pattern Recognition
Letters*, 29(7), 950–956.

Linero, A. R. (2025). Generalized Bayesian additive regression trees
models: beyond conditional conjugacy. *Journal of the American
Statistical Association*, 120(549), 356–369.
<https://doi.org/10.1080/01621459.2024.2337156>

Linero, A. R., & Yang, Y. (2018). Bayesian regression tree ensembles
that adapt to smoothness and sparsity. *Journal of the Royal Statistical
Society Series B*, 80(5), 1087–1110.
<https://doi.org/10.1111/rssb.12293>

Linero, A. R. (2018). Bayesian regression trees for high-dimensional
prediction and variable selection. *Journal of the American Statistical
Association*, 113(522), 626–636.
<https://doi.org/10.1080/01621459.2016.1264957>

Kubinec, R. (2023). Ordered beta regression: a parsimonious,
well-fitting model for continuous data with lower and upper bounds.
*Political Analysis*, 31(4), 519–536.
<https://doi.org/10.1017/pan.2022.20>

The MCMC engine is adapted from Antonio Linero’s `FlexBart` reference
implementation, distributed with Linero (2025) and licensed GPL-2.
