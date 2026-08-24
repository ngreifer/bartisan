# genbart

Bayesian additive regression trees for response distributions that standard BART
implementations cannot reach, with an interface that looks like `glm()`.

``` r
genbart(y ~ x1 + x2, data = d, family = binomial("logit"))
```

## Why this exists

BART is normally tied to a Gaussian likelihood. The Gibbs sampler integrates the
leaf parameters out of the tree likelihood in closed form, and that step needs
conjugacy. Getting anything else — a logistic model, a survival model — usually
means finding a data-augmentation scheme that turns the problem back into a
Gaussian one, and for many likelihoods no such scheme exists.

*genbart* implements the approach of Linero (2025), which removes the
restriction. Rather than integrating the leaf parameters out, the sampler builds
a Gaussian approximation to their conditional posterior by Fisher scoring, and
uses that as the proposal in a reversible-jump Metropolis step. The
approximation only has to be good enough to be accepted often; the chain
targets the exact posterior regardless of how good it is.

The practical consequence is that adding a family requires only the log density
of one observation and its first two derivatives with respect to the additive
predictor. That is why the list below includes models that are otherwise hard to
find in a BART package.

| Family | Compiled links | Nuisance parameters drawn |
|---|---|---|
| `gaussian()` | identity | residual standard deviation |
| `binomial()` | logit, probit, cloglog | — |
| `poisson()` | log | — |
| `negbin()` | log | dispersion |
| `Gamma()` | log | shape |
| `ordinal()` | logit, probit, cloglog | cutpoints |
| `multinomial()` | logit | — |
| `weibull_aft()`, `loglogistic_aft()`, `lognormal_aft()` | — | scale |
| `location_scale()` | identity | — |
| `zi_poisson()` | log | — |
| `zi_negbin()` | log | dispersion |
| `ordbeta()` | logit | 2 cutpoints, precision |

Ordinary `stats::family` objects work unchanged, so moving a model from `glm()`
to `genbart()` is a one-word change. `multinomial()`, `location_scale()` and the
two zero-inflated families fit one forest per unconstrained parameter, so the
excess-zero mechanism of a zero-inflated count model is itself nonparametric
rather than a single constant. `ordbeta()` is the ordered beta regression of
Kubinec (2023), for a proportion or slider response with mass piled up at zero
and one. `multinomial()` follows Murray (2021) in fitting one forest per
category and leaving the model unidentified, so that the prior does not depend
on which category was singled out as the reference.

## Links and likelihoods of your own

The links in the table are the ones the sampler evaluates in compiled code. Any
other link works for the five families with a single mean and a conventional
link, applied from R by composing it onto the scale the compiled family works
on, so `binomial("cauchit")` needs no new syntax and neither does a link you
write yourself and pass the way `glm()` takes one.

Because the whole interface a family has to satisfy is a log density and two
derivatives, it can also be supplied from R outright:

```r
# A Poisson model written by hand. Terms free of eta may be dropped.
genbart(y ~ ., data = d,
        family = custom_family(function(y, eta) y * eta[, 1] - exp(eta[, 1])))

# Two additive predictors: a mean and a log standard deviation.
genbart(y ~ ., data = d,
        family = custom_family(
          function(y, eta) dnorm(y, eta[, 1], exp(eta[, 2]), log = TRUE),
          num_predictors = 2))
```

The function is called once per leaf, not once per observation, so it must be
vectorized over its arguments. Central differences supply the derivatives unless
you pass them; a custom family cannot draw a nuisance parameter.

## Missing predictor values

A missing predictor is **kept by default**, not dropped and not imputed. Each
splitting rule carries the answer for itself: a rule on a variable that has any
missing values is drawn as one of three, with equal probability, sending them left
with the low values, right with the high ones, or splitting on whether the value
is there at all. Pass `na.action = na.omit` for the usual behavior of dropping
the row.

```r
genbart(y ~ ., data = d)                          # missing predictors kept
genbart(y ~ ., data = d, na.action = na.omit)     # rows dropped instead
```

A missing response, weight or offset is dropped either way, with a warning, since
there is nothing to fit those rows to.

This is *missingness incorporated in attributes* (Twala et al. 2008; for BART,
Kapelner and Bleich 2015). The third rule is the one that matters: it lets the
model use a variable whose **absence** is the signal, which no imputation scheme
can recover. On a problem where the response depends on whether a predictor was
recorded and on nothing else about it, the fitted means came out at 1.98 and
-0.01 against a truth of 2 and 0, with 29.6 of the 30 splitting rules landing on
that variable.

Because the choice is drawn from its prior alongside the variable and the
cutpoint, it cancels from every acceptance ratio, and a variable with no missing
values is not given the extra draw at all -- complete data reproduces the
sampler exactly as it was. What this estimates is the mean of the response given
the predictors *and* the pattern of missingness, which is what you want for
prediction; if the estimand is defined on complete data, multiple imputation is
the right tool instead.

## Smoother fits

By default the decision rules are *soft*, following Linero and Yang (2018): a
rule is a logistic gate rather than a step, so an observation reaches every leaf
with some weight and the fitted function is smooth rather than piecewise
constant. Combining soft rules with a non-conjugate likelihood is an extension
of Linero (2025), which leaves it as an open problem — soft rules make the leaf
parameters of a tree dependent on one another, so the reversible-jump move needs
a bivariate Laplace proposal for the pair of child leaves. That proposal reduces
exactly to Linero's independent pair when the rules are hard, so the two cases
share one implementation.

Set `soft = FALSE` in `genbart_control()` for the faster hard-rule sampler.

## Example


``` r
library(genbart)

set.seed(1)

n <- 400
d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n), x4 = runif(n))
d$y <- rbinom(n, 1, plogis(3 * sin(pi * d$x1 * d$x2) - 1))

fit <- genbart(y ~ ., data = d, family = binomial())

fit
```

```
## Generalized BART
## 
## Call: `genbart(formula = y ~ ., data = d, family = binomial())`
## 
## Family: "binomial" with the "logit" link
## Observations: 400
## Forests: 1 of 50 trees, soft decision rules
## Draws: 500 kept after 500 warmup
```

``` r
summary(fit)          # includes how often each predictor was split on
```

```
## Generalized BART
## 
## Call: `genbart(formula = y ~ ., data = d, family = binomial())`
## 
## Family: "binomial" with the "logit" link
## Observations: 400
## Forests: 1 of 50 trees, soft decision rules
## Draws: 500
## 
## Predictor usage
## Splitting rules per draw, and how often used at all.
##      mean     sd lower  upper prop_used
## x1 25.866  8.745     9 43.525     1.000
## x2 20.002  8.914     6 39.000     1.000
## x3 17.578 12.292     0 41.000     0.964
## x4 13.182  9.156     0 35.525     0.956
```

``` r
head(predict(fit, type = "response"))
```

```
## [1] 0.6353835 0.4244375 0.7686371 0.8328150 0.6070707 0.8538756
```

``` r
predict(fit, newdata = d[1:5, ], type = "prob")
```

```
##              0         1
## [1,] 0.3646165 0.6353835
## [2,] 0.5755625 0.4244375
## [3,] 0.2313629 0.7686371
## [4,] 0.1671850 0.8328150
## [5,] 0.3929293 0.6070707
```

Every draw of every tree is retained, so `predict()` can return the full
posterior rather than a point estimate:


``` r
draws <- predict(fit, newdata = d[1:5, ], type = "response", draws = TRUE)
apply(draws, 2, quantile, c(0.025, 0.975))
```

```
##            [,1]      [,2]      [,3]      [,4]      [,5]
## 2.5%  0.4741241 0.2493829 0.5849570 0.7048386 0.4309315
## 97.5% 0.7808235 0.6162683 0.9010255 0.9229802 0.7688138
```

`type = "density"` evaluates the conditional density of the outcome given the
predictors, at the outcome's own value, so the outcome has to be present in
`newdata`. Summed on the log scale that is a held-out log score, which is a fair
basis for comparing models:


``` r
held_out <- data.frame(x1 = runif(50), x2 = runif(50), x3 = runif(50),
                       x4 = runif(50))
held_out$y <- rbinom(50, 1, plogis(3 * sin(pi * held_out$x1 * held_out$x2) - 1))

sum(predict(fit, newdata = held_out, type = "density", log = TRUE))
```

```
## [1] -27.60962
```

A right-censored survival model, with the response given as a `Surv` object:

``` r
genbart(survival::Surv(time, status) ~ ., data = d, family = weibull_aft())
```

A zero-inflated count model, where both the count mean and the structural-zero
probability get their own forest:

``` r
genbart(count ~ ., data = d, family = zi_negbin())
```

## Random intercepts

A `(1 | group)` term adds an intercept per level of `group`, drawn from a common
mean-zero normal whose standard deviation is itself drawn:

```r
genbart(y ~ x1 + x2 + (1 | school), data = d)
genbart(y ~ x1 + (1 | school) + (1 | year), data = d)   # several factors
genbart(y ~ x1 + (1 | district / school), data = d)     # nesting
```

Each additive predictor gets its own set, so a zero-inflated count model has a
group effect on the count part and another on the inflation part. Only random
*intercepts* are supported; a slope is refused rather than ignored, because a
variable whose effect varies by group belongs in the fixed part of the formula,
where a tree can split on the group and the variable together and get an
interaction of any shape.

Whether to reach for this at all depends on the shape of the data, and the answer
is not "always". A grouping factor can go in the fixed part instead, where a tree
splits on it like anything else. Measured against `dbarts::rbart_vi()` on data
with a true group effect, held-out RMSE:

| groups | per group | group as a predictor | random intercept |
|---|---|---|---|
| 5 | 100 | **0.190** | 0.279 |
| 25 | 20 | **0.306** | 0.341 |
| 100 | 5 | 0.541 | **0.468** |
| 250 | 4 | 0.708 | **0.493** |

With few large groups the predictor route is better — the group means are well
determined without pooling, and a split can interact the group with the
covariates. The random intercept wins once there are many small groups, which is
where partial pooling earns its keep.

## Counterfactual estimands

A fit works with **marginaleffects**, so predictions, comparisons and slopes — and
hypothesis tests on any of them — come out without extracting draws by hand:

```r
library(marginaleffects)

avg_predictions(fit)
avg_comparisons(fit)                        # average effect of each predictor
avg_predictions(fit, by = "group")          # by a factor
comparisons(fit, variables = "treatment",   # at a grid of covariate values
            newdata = datagrid(age = c(30, 50, 70)))
```

A forest has no coefficient vector and no variance-covariance matrix, so there is
no delta method to apply. There is something better: every estimand is computed by
pushing all of the posterior draws through the same transformation, so the
intervals are posterior quantiles rather than a normal approximation, and a
nonlinear estimand needs no approximation at all.

One caveat, documented in `?genbart-marginaleffects`. `avg_slopes()` is a
numerical derivative, and the default `x_transform = "quantile"` maps each
predictor through its empirical distribution function — a step function — so the
fit is a step function of the original predictor whatever the decision rules are.
Use `x_transform = "range"` if slopes are the estimand. Predictions and
comparisons evaluate the fit at two points a substantive distance apart and are
unaffected.

## Variable selection

Splitting variables are drawn from a Dirichlet prior that concentrates on few
predictors, following Linero (2018), so the model filters out irrelevant ones.
`summary()` reports the posterior distribution of the number of splitting rules
using each predictor and the proportion of draws using it at all. The indicator
columns of a factor share one weight, so a factor is selected as a whole rather
than one level at a time.

## Correctness

Three checks are in the test suite and worth knowing about.

Shrinking the prior weights to nothing makes the likelihood constant, so the
posterior collapses to the tree prior. The sampled trees then reproduce the
prior's branching process to within a fraction of a percent, which is a direct
check that the birth, death and change moves satisfy detailed balance.

Every family's analytic score is checked against a central difference of its own
log density, which is an independent route to the same quantity. This is the
check that matters for the sampler: a wrong score sends the Fisher-scoring
proposal to the wrong place.

Across replicate datasets at the default settings, 95% credible intervals for
the additive predictor cover the truth 0.95 (Gaussian), 0.91 (binomial), 0.96
(Poisson) and 0.96 (gamma) of the time. The binomial case is the weakest because
a binary response carries the least information per observation. The shortfall
tracks the ratio of absolute bias to posterior standard deviation, which sits
near 0.8: an interval centered on a shrunken estimate loses coverage in
proportion to how far that centre sits from the truth. Running more or longer
chains does not help, which rules out mixing as the cause at this scale. This is
in line with what is reported for BART generally; treat pointwise intervals for the
regression function as approximate.

## Speed

A non-conjugate sampler costs more than a conjugate one, and it is worth being
precise about how much -- and about what the extra buys. From
`_dev/benchmark.Rmd`, which you can run yourself: the Friedman function,
n = 1000, 10 predictors, 50 trees, 1000 warmup and 1000 saved draws, scored
against the true regression function on a held-out thousand.

| Task | Package | Time | Effective sample size | Held-out RMSE |
|---|---|---|---|---|
| Gaussian | `dbarts` | 0.24 s | 20 | 0.221 |
| | **genbart, hard rules** | **0.44 s** | 18 | 0.219 |
| | `stochtree` | 0.70 s | 18 | 0.218 |
| | `bartMachine` | 0.87 s | — | 0.230 |
| | `BART::wbart` | 1.06 s | 21 | 0.228 |
| | genbart, soft, smoothstep gate | 1.40 s | 42 | **0.135** |
| | genbart, soft, smootherstep gate | 1.42 s | 43 | 0.140 |
| | genbart, soft rules (default) | 2.00 s | 44 | 0.145 |
| probit | `dbarts` | 0.27 s | 36 | 0.134 |
| | **genbart, hard rules** | **0.49 s** | 35 | 0.122 |
| | `stochtree` | 0.99 s | 29 | 0.134 |
| | `BART::pbart` | 1.13 s | 36 | 0.125 |
| | genbart, soft rules (default) | 2.05 s | 33 | **0.111** |
| | genbart, `augment = FALSE` | 23.6 s | 65 | 0.106 |
| logit | genbart, soft rules | 2.11 s | 129 | **0.088** |
| | genbart, `augment = FALSE` | 13.9 s | 104 | 0.085 |
| | `BART::lbart` | 20.8 s | 38 | 0.109 |
| ordinal | genbart, hard, probit | 0.87 s | 35 | — |
| | genbart, hard, logit | 0.93 s | 31 | — |
| | genbart, soft, probit | 2.38 s | 51 | — |
| | genbart, soft, logit | 2.52 s | 48 | — |
| | **genbart, hard, cloglog** | **3.03 s** | 24 | — |
| | `stochtree` (cloglog) | 4.01 s | 20 | — |
| | genbart, hard, logit, `augment = FALSE` | 14.9 s | 42 | — |
| | genbart, hard, cloglog, `augment = FALSE` | 17.0 s | 37 | — |
| | genbart, hard, probit, `augment = FALSE` | 25.2 s | 36 | — |

With hard rules genbart is **within a factor of 1.8 of `dbarts`** on the two
tasks `dbarts` supports — 0.426 s against 0.241 s on the Gaussian task and 0.482 s
against 0.272 s on probit, best of six runs each — at the same mixing and better
accuracy, and it is faster *and* more accurate than every other package here on
both. On the complementary log-log ordinal model, the one task `stochtree`
supports and `dbarts` does not, genbart is now the faster of the two.

The other rows come from `_dev/benchmark.Rmd` at two replicates, which is noisy at
the ten to thirty percent level; the two ratios above are quoted from a longer run
because a claim about a factor of two should not rest on a best-of-two. On a logit link it is
9.9 times faster than `BART::lbart` and mixes three times better. Soft rules --
the default -- cost a further four to five times and cut the held-out error by 30
to 40%, which makes the default configuration the most accurate fit in the table;
with a bounded gate that falls to about three times.

The ordinal rows have no RMSE column because the links put the additive predictor
on different scales, so the numbers would compare scales rather than fits; the
timings and the mixing are comparable.

The remaining families are the expensive ones, at 2 to 8 seconds, which is the
price of the general machinery where no rewriting is available.

Both genbart figures use `augment = TRUE`. The gap decomposes into factors that
were measured rather than guessed, and profiling the sampler by phase says where
each one lives.

With **hard rules**, what is left against `dbarts` is tree bookkeeping, not the
family. A birth, death or change move makes five or six passes over the node it
touches -- materializing the predictor with the node's own contribution removed,
fitting the proposal, splitting the node's membership, recomputing the child
weights, committing -- against about two for a sampler that only ever needs a sum
of residuals and a count. Nothing in that gap is the likelihood: the acceptance
ratio genbart forms for a Gaussian response is algebraically the same
marginal-likelihood ratio `dbarts` forms, because a quadratic target makes the
Laplace approximation exact.

What has actually paid, in order, is worth knowing before reaching for the
obvious ideas. **Sizing the children's index vectors and filling them by index,
instead of `push_back`,** is worth 20% of a soft fit on its own -- four vectors'
worth of per-element capacity tests in the innermost loop of the sampler.
**Recycling tree nodes** rather than allocating two per birth proposal, about two
thirds of which are rejected, is worth 5%. **Folding the child weights into the
split that produces them** removes a whole pass over the same gates.

Two things that looked like they should help did not. A hard tree used to carry a
membership weight per observation per node, every one exactly 1.0; and the
leaf-level sums used to write a node's predictors into a buffer and read them
back. Removing both is bit-identical -- multiplying by 1.0 is exact -- and
together they are worth about 10%, against a predicted factor of several. Memory
bandwidth was never the constraint: a node holds a few thousand observations at
most, so the buffer stays in L1. Nor is dropping the remaining materialized base
array worth much: measured directly, at about 2%, and it would trade a sequential
read for a gather.

With **soft rules**, half the cost is a single move. The bandwidth of each tree
is a parameter with a Metropolis step per tree per sweep, and every attempt
rebuilds every membership weight in the tree. A rejected attempt -- 58% of them --
used to rebuild a second time to get back where it started; it is now rolled back
from a snapshot, and a tree with no splits skips the move entirely and draws its
bandwidth from the prior, which is its exact full conditional. The other half is that a logistic gate never saturates:
dropping a weight needs it below 1e-10, which needs the observation 23 bandwidths
from the cutpoint, further than the whole unit interval. So every observation
reaches every leaf and a pass over a node covers 2.5 times the sample.
`gate = "smoothstep"` and `gate = "smootherstep"` are the fix for the second half
and are in the table above -- though not for the reason they were built. At a
bandwidth wide enough that a bounded gate truncates *nothing*, it is still 1.45
times faster than the logistic: what a bounded gate saves is the `exp()`, not the
work on the far side of the cutpoint. So the two bounded gates come out within
noise of each other and the choice between them is about smoothness -- one
derivative against two -- rather than speed.
Fixing the bandwidth outright is the fix for the first, and is *not* the default,
because on a function with jumps it more than doubles the error -- letting the
rules sharpen towards hard ones is what the move is for.

The passes over the data are no longer part of it. A leaf value enters the
additive predictor linearly, so the log target over a leaf inherits the shape of
the log density -- and where that shape is known, one pass over a node
determines the whole function and everything after it is arithmetic. Two shapes
are known. A *quadratic* log density (a Gaussian response, and anything `augment`
rewrites into one) makes the Laplace approximation the conditional posterior
exactly. An *exponential* one, `a*eta + b*exp(s*eta)`, covers the Poisson, the
gamma and the augmented negative binomial under hard rules; there the fit still
has to be iterated, but on three numbers rather than on the data. Either way a
birth move goes from six passes over the node to two.

What closes most of the first factor is a data augmentation that makes the
conditional Gaussian, since then the Laplace approximation is exact rather than
approximate. `augment` does this: Albert and Chib (1993) for a probit link and for an ordinal
probit, Pólya-Gamma augmentation (Polson, Scott and Windle 2013) for a logit
one, an ordinal logit, or a multinomial. Every one of them trades speed for mixing,
so the ratio to judge is effective samples per second, and it differs a lot by
family:

| Family | Speed | Effective sample size | ESS per second |
|---|---|---|---|
| `ordinal("probit")`, hard rules | 30x | 0.76x | **23x** |
| `ordinal("logit")`, hard rules | 15x | 0.94x | **14x** |
| `ordinal("probit")`, soft rules | 14x | 0.87x | **12x** |
| `ordinal("logit")`, soft rules | 7.1x | 0.76x | **5.4x** |
| `binomial("probit")` | 5.7x | 0.66x | **3.8x** |
| `binomial("logit")` | 2.6x | 0.81x | **2.1x** |
| `multinomial()` | 4.2x | 0.37x | 1.6x |
| `negbin()` | 0.8x | 0.56x | 0.5x |

`augment` is on by default and covers the binomial and ordinal families, plus the
negative binomial when the rules are hard -- there it is written as a Poisson whose rate
comes from a gamma, which needs one gamma draw per observation rather than a
Pólya-Gamma one. `augment = "multinomial"` asks for a rewriting that pays only if
the chain is lengthened to match.

The ordinal probit is the largest of these because the target it replaces is the
most expensive one in the package: two cumulative-normal evaluations per
observation per pass, and no exploitable shape, so every trial value of a leaf
parameter costs its own pass. It comes with one subtlety worth naming. The
cutpoints are *not* drawn from their conditional given the latent normals --
uniform between the two order statistics that bracket each one, which pins them
to an interval of width O(1/n) and mixes worse the more data there is. They are
drawn from the ordinal likelihood with the latent normals integrated out, and the
latent normals are redrawn immediately afterwards: a partially collapsed Gibbs
sampler (Van Dyk and Park 2008), which is the standard remedy (Cowles 1996).

The cumulative **logit** has no latent normal -- its latent variable is logistic --
and the usual route to one goes through the Kolmogorov-Smirnov density, which
needs a sampler of its own. It does not have to. Polson, Scott and Windle's
Theorem 1 at *a* = 1, *b* = 2 says that the standard logistic density is
`(1/4) E[exp(-w x²/2)]` with *w* ~ PG(2, 0), so a logistic residual is a normal
whose precision is Pólya-Gamma -- and because PG(*b*, *c*) is PG(*b*, 0) tilted by
`exp(-c²w/2)`, the precision's conditional given a residual *r* is exactly
PG(2, |*r*|). That is an integer-parameter draw, which the exact Devroye sampler
already in the package covers, so the logit link gets the same quadratic target
with nothing approximate anywhere. The Poisson and gamma families need no
rewriting: their targets are already in a shape the sampler collapses to a single
pass.

One warning worth heeding: an unoptimized build of the compiled code is 5 to 20
times slower and looks no different. `devtools::load_all()` and
`devtools::install()` compile without optimization and leave the object files
for a later `R CMD INSTALL` to reuse. `genbart()` warns once per session when it
detects one; `genbart:::.genbart_optimized()` reports the state directly.

Independent chains run in parallel:

```r
future::plan(future::multisession, workers = 4)

fit <- genbart(y ~ ., data = d, chains = 4)
fit$rhat
```

The chain is the only parallel axis this sampler has -- a sweep conditions on
the last one -- and it is also what makes a convergence diagnostic possible.
Any `future` backend works, including mirai's. One `set.seed()` reproduces the
whole run whatever the backend.

## Notes and limitations

Proportional hazards regression is deliberately absent. Its partial likelihood
couples observations through risk sets, so it does not decompose into a sum over
the observations reaching a leaf, which is what the leafwise Laplace
approximation requires. The accelerated failure time families cover the same
ground within this framework.

Soft rules cost more per iteration than hard ones, because a leaf touches every
observation rather than only those in its cell. Negligible weights are pruned,
which holds the gap to roughly a factor of three rather than the factor of *L*
it would otherwise be.

Absolute timings are machine- and load-dependent, so treat these as indicative
rather than exact; repeated runs on the same laptop varied by about 20%. As an
anchor, a Gaussian response with 10 predictors, 50 soft trees and the default
1000 warmup plus 1000 saved draws takes roughly 6 seconds at n = 500 and roughly
a minute at n = 5000, measured at steady state on one core of an M-series Mac.
Hard rules are about 3 times faster, and doubling the tree count roughly doubles
the cost.

The ratios are stabler than the absolute times. Cost relative to a Gaussian fit
at the same size is set by how expensive the family's log density and its
derivatives are:

| Family | Relative cost |
|---|---|
| `gaussian()` | 1x |
| `Gamma()`, `poisson()` | ~2x |
| `negbin()`, `binomial()` | 2–3x |
| `ordinal()` | ~2x with a latent variable, ~25x without |
| `zi_poisson()` | ~10x |
| `ordbeta()` | ~25x |

`ordbeta()` is the expensive one because every evaluation of its log density
needs two log-gamma calls that no amount of restructuring removes. `ordinal()`
spans a wide range because all three of its links have a latent-variable
representation that the sampler uses by default — a normal for the probit, a
normal with a Pólya-Gamma precision for the logit, an exponential waiting time
for the complementary log-log — and `augment = FALSE` turns them off.

An ordinal model is identified only up to a common shift of its cutpoints and its
predictor. With three or more categories the draws are reported in the chart where
**the predictor has mean zero over the fitted sample and every cutpoint is free**,
which makes the cutpoints readable as category boundaries. With two categories the
single boundary is folded into the intercept instead, so the fit is on the same
scale as `binomial()`.

That is `MASS::polr()`'s chart with its predictors centered. `polr()` identifies
the location by leaving the intercept out of the design matrix rather than by
centering, so its `zeta` is shifted by the mean of its own linear predictor.
Either of these lines puts the two side by side, and on a linear truth they agree
to Monte Carlo error:

```r
p <- MASS::polr(y ~ x1 + x2, data = d, method = "probit")
p$zeta - mean(p$lp)                      # same chart as colMeans(fit$aux)
```

Two prediction types exist for reporting an ordinal fit on a single scale, both
following **WeightIt**. `predict(type = "mean")` weights the category
probabilities by the labels read as numbers, so levels `"1"`, `"2"`, `"4"` give a
mean between one and four; `values` says what the categories are worth when the
labels are not numbers. `predict(type = "stdlv")` divides the predictor by the
standard deviation of the latent variable it indexes, which puts fits with
different links, or different amounts of signal, on one scale.

The leaf scale is raised from near zero over the first quarter of warmup.
Linero (2025) describes this as essential: started at its full value the sampler
can settle early into a poor configuration and fail to move. It is controlled by
`sigma_mu_ramp` in `genbart_control()`.

Nuisance parameters that are weakly identified — the negative binomial
dispersion, the shape of a survival model — mix more slowly than the regression
function, and the data often do not pin them down. Check their draws before
interpreting them.

The leaf scale is drawn under a half-Cauchy prior, which has no upper bound.
When the predictors nearly separate a binary response the likelihood rewards an
unbounded predictor and that prior does not hold the scale down: in a fully
separated example the drawn scale wandered between 3 and 9 times its prior
median over 1600 draws without settling. `genbart()` warns when the scale
settles more than five times above its prior median; `update_sigma_mu = FALSE`
pins it.

A link supplied from R costs a call into the interpreter for every leaf the
sampler visits, and `custom_family()` costs three unless you supply the
derivatives. Both are usable at the sizes a BART model is usually fit at, and
both are slower than the compiled equivalent.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("ngreifer/genbart")
```

Compilation requires a C++ toolchain, plus *Rcpp* and *RcppArmadillo*.

## References

Albert, J. H., & Chib, S. (1993). Bayesian analysis of binary and polychotomous
response data. *Journal of the American Statistical Association*, 88(422),
669--679.

Kapelner, A., & Bleich, J. (2015). Prediction with missing data via Bayesian
additive regression trees. *Canadian Journal of Statistics*, 43(2), 224--239.

Polson, N. G., Scott, J. G., & Windle, J. (2013). Bayesian inference for logistic
models using Polya-Gamma latent variables. *Journal of the American Statistical
Association*, 108(504), 1339--1349.

Murray, J. S. (2021). Log-linear Bayesian additive regression trees for
multinomial logistic and count regression models. *Journal of the American
Statistical Association*, 116(534), 756--769.

Twala, B. E. T. H., Jones, M. C., & Hand, D. J. (2008). Good methods for coping
with missing data in decision trees. *Pattern Recognition Letters*, 29(7),
950--956.

Linero, A. R. (2025). Generalized Bayesian additive regression trees models:
beyond conditional conjugacy. *Journal of the American Statistical Association*,
120(549), 356–369. <https://doi.org/10.1080/01621459.2024.2337156>

Linero, A. R., & Yang, Y. (2018). Bayesian regression tree ensembles that adapt
to smoothness and sparsity. *Journal of the Royal Statistical Society Series B*,
80(5), 1087–1110. <https://doi.org/10.1111/rssb.12293>

Linero, A. R. (2018). Bayesian regression trees for high-dimensional prediction
and variable selection. *Journal of the American Statistical Association*,
113(522), 626–636. <https://doi.org/10.1080/01621459.2016.1264957>

Kubinec, R. (2023). Ordered beta regression: a parsimonious, well-fitting model
for continuous data with lower and upper bounds. *Political Analysis*, 31(4),
519–536. <https://doi.org/10.1017/pan.2022.20>

The MCMC engine is adapted from Antonio Linero's `FlexBart` reference
implementation, distributed with Linero (2025) and licensed GPL-2.
