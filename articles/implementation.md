# Generalized BART with bartisan

``` r

library(bartisan)

set.seed(2026)
```

## Introduction

This vignette describes what the model is and how it is fitted. It is
the reference document rather than the entry point:
[`vignette("bartisan")`](https://ngreifer.github.io/bartisan/articles/bartisan.md)
shows how to use the package on a real analysis, and
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
covers choosing a likelihood. This is the one to read to find out what
the fit is doing, why the defaults are what they are, or what an
advanced setting changes.

In this guide we will first set out the model, which is a sum of trees
supplying an additive predictor to an arbitrary likelihood, and then the
prior that regularizes it. Next we will describe how it is fitted, which
is where the generalization actually lives, and what data augmentation
buys where it is available. We then turn to the practical matters that
come up in a real fit, to what the sampler costs against the other BART
packages, and to the checks that say the implementation is right.

Fits below use a deliberately small chain so that the vignette builds
quickly: 20 trees and 300 draws after 300 warmup iterations, on 400
observations. The defaults are 50 trees and 500 draws after 500 warmup
iterations.

``` r

ctrl <- bartisan_control(num_trees = 20, num_burn = 300, num_draws = 300)

n <- 400
```

## The Model

### A Sum of Trees

A binary regression tree is a pair \\(T, M)\\, where \\T\\ is the tree
structure, meaning the splitting rules at the interior nodes, and \\M =
\\\mu_1, \dots, \mu_B\\\\ collects one parameter for each of its \\B\\
leaves. The tree defines a function \\g(x; T, M)\\ that sends \\x\\ down
the tree according to the rules and returns the parameter in the leaf it
reaches. A single tree is a step function.

BART sums many of them ([Chipman et al. 2010](#ref-chipman2010)):

\\f(x) = \sum\_{m=1}^{M} g(x; T_m, M_m).\\

Each tree is kept small by its prior, so no single one explains much.
The sum is flexible while each term is a weak learner, which is what
makes the fit stable. Interactions come free: a path through a tree that
splits on \\x_1\\ and then on \\x_2\\ is an interaction between them,
and nothing had to be specified for it to appear.

In ordinary BART the sum is the conditional mean of a Gaussian outcome,

\\Y = f(x) + \varepsilon, \qquad \varepsilon \sim \mathrm{N}(0,
\sigma^2).\\

### The Generalization

*bartisan* keeps the sum of trees and replaces what sits on top of it.
The forest supplies an additive predictor \\\eta(x)\\, and the outcome
follows whatever density the family names,

\\\eta(x) = \sum\_{m=1}^{M} g(x; T_m, M_m), \qquad Y \mid x \sim p\\y
\mid \eta(x), \theta\\,\\

with \\\theta\\ collecting any parameters of the density that are not
functions of \\x\\, such as a residual standard deviation or a shape.
These are the nuisance parameters, and they are drawn alongside the
trees.

Some families need more than one function of \\x\\. A location-scale
model has a mean and a log standard deviation; a zero-inflated count has
a rate and a probability of a structural zero. Those get one forest
each:

\\\eta_h(x) = \sum\_{m=1}^{M_h} g(x; T\_{hm}, M\_{hm}), \qquad h = 1,
\dots, H,\\

and the sampler cycles over all the forests.
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
says how many each family uses.

### The Prior

Three pieces, all with defaults chosen so that the fit is regularized
without being tuned ([Chipman et al. 2010](#ref-chipman2010)).

**Tree structure.** A node at depth \\d\\ is split rather than left as a
leaf with probability

\\\alpha (1 + d)^{-\beta}, \qquad \alpha = 0.95, \\ \beta = 2.\\

This falls off quickly, so trees are shallow: most have two or three
leaves. It is the main reason a forest of 50 trees does not overfit.

**Leaf values.** Each \\\mu_b\\ is normal with mean zero and a standard
deviation that shrinks as trees are added, so that the sum has a
sensible scale whatever \\M\\ is. The `k` argument in
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
sets how hard this shrinks: larger values of `k` pull the fit toward the
intercept-only model. The default is a compromise that works across a
wide range of problems, and it is rarely worth changing.

**Nuisance parameters.** Each family sets its own, for example a
half-Cauchy on a scale parameter.

The prior does the work that cross-validation does in other tree
ensembles. There is no learning rate, no depth limit, and no early
stopping.

### Soft Rules (`gate`)

A rule in classic BART is a step: an observation with \\x_j \le c\\ goes
left and everything else goes right. The fitted function is then a step
function, which is a poor description of a smooth truth.

Soft rules replace the step with a gate (i.e., a smooth function of the
distance from the cutpoint) ([Linero and Yang 2018](#ref-linero2018)).
An observation goes left with weight

\\\psi\left(\frac{x_j - c}{\tau}\right),\\

where \\\psi\\ is a cumulative distribution function and \\\tau\\ is a
bandwidth, and right with the remaining weight. Every observation
reaches every leaf with some weight, and the fitted function is smooth.
As \\\tau \to 0\\ the gate becomes a step and the hard rule is
recovered.

The `gate` argument chooses \\\psi\\: `"smoothstep"` (the default),
`"smootherstep"`, `"logistic"`, or `"hard"`. The first two are
polynomial and reach exactly zero and one outside a finite window, which
is faster than the logistic and usually just as good. `bandwidth` sets
the prior mean of \\\tau\\, which is drawn rather than fixed.

Soft rules cost about three times as much per iteration and are usually
worth it:

``` r

friedman <- function(n) {
  x <- as.data.frame(matrix(runif(n * 10), n, 10))
  names(x) <- paste0("x", 1:10)
  x$eta <- 10 * sin(pi * x$x1 * x$x2) + 20 * (x$x3 - 0.5)^2 +
    10 * x$x4 + 5 * x$x5
  x
}

train <- friedman(n)
test  <- friedman(1000)
train$y <- train$eta + rnorm(n)

timed <- function(gate) {
  started <- proc.time()[["elapsed"]]
  fit <- bartisan(y ~ . - eta, data = train, family = gaussian(),
                  control = bartisan_control(num_trees = 20, num_burn = 300,
                                             num_draws = 300, gate = gate))
  data.frame(rules = gate,
             test_rmse = sqrt(mean((predict(fit, newdata = test) - test$eta)^2)),
             seconds = round(proc.time()[["elapsed"]] - started, 1))
}

rbind(timed("smoothstep"), timed("hard"))
#>        rules test_rmse seconds
#> 1 smoothstep     0.434     1.3
#> 2       hard     1.025     0.4
```

The true function has a standard deviation of about 4.9, so both are
fitting real structure and the soft fit is the more accurate.

### Sparsity (`sparsity`)

By default the variable a rule splits on is drawn from a categorical
distribution whose probabilities have a Dirichlet prior ([Linero
2018](#ref-linero2018sparse)):

\\s \sim \mathrm{Dirichlet}(\alpha/p, \dots, \alpha/p),\\

with \\\alpha\\ itself drawn. Small \\\alpha\\ concentrates the
probability on a few predictors, so the forest can stop splitting on the
rest entirely. This is DART, and it is what `sparsity = TRUE` means. It
is the default.

The alternative, `sparsity = FALSE`, gives every predictor the same
splitting probability, which is classic BART.

The difference matters for variable selection. Chipman et al.
([2010](#ref-chipman2010)) note that counting splits works poorly when
there are many trees, “because the redundancy offered by so many trees
tends to mix many irrelevant predictors in with the relevant ones”, and
recommend reducing the number of trees so that predictors compete. The
Dirichlet prior addresses the same problem directly, by letting
irrelevant predictors be dropped rather than crowded out.

Both effects are visible on the Friedman function, where `x1` to `x5`
matter and `x6` to `x10` are noise. The table reports the average
posterior probability of being used, over three replicates at \\n =
500\\:

| Sparsity           | 10 trees | 20 trees | 50 trees | 100 trees |
|:-------------------|---------:|---------:|---------:|----------:|
| `sparsity = FALSE` |     0.28 |     0.50 |     0.95 |      1.00 |
| `sparsity = TRUE`  |     0.09 |     0.08 |     0.09 |      0.14 |

Mean prop_used for the five noise predictors. The five real predictors
sit at 1.00 in every cell. {.table}

Without the sparsity prior the advice to use fewer trees is essential:
at 50 trees the noise predictors are used in 95% of draws and are
indistinguishable from the real ones. With it, the noise predictors stay
near zero at every tree count, and reducing the trees buys little. The
recommendation is a workaround for the absence of the prior rather than
a property of variable selection in general.

Two things this does not fix. Sparsity needs signal, so on a small
sample with a weak fit nothing separates, and the honest reading of a
flat table is that the data cannot rule anything out; and it says
nothing about effect size, since a predictor can be split on constantly
and move the prediction very little.
[`vignette("importance")`](https://ngreifer.github.io/bartisan/articles/importance.md)
covers both.

## How It Is Fitted

### Bayesian Backfitting

The trees are updated one at a time. To update tree \\m\\, subtract the
contribution of every other tree from the predictor, leaving a partial
residual that only tree \\m\\ has to explain, then draw \\(T_m, M_m)\\
given that residual. Cycling over all trees, and then over the nuisance
parameters, is one sweep. The draws saved after warmup are the posterior
sample.

### Moving the Trees

A tree is proposed by a small local change: growing a leaf into two,
pruning a pair of leaves back into one, or changing a splitting rule.
The proposal is accepted or rejected by a Metropolis-Hastings step.
Because growing and pruning change the number of parameters, this is a
reversible-jump sampler.

The acceptance ratio needs the likelihood of the tree with the leaf
parameters integrated out,

\\p(\text{residual} \mid T) = \int p(\text{residual} \mid T, M) \\ p(M)
\\ dM.\\

### The General Case, and the Laplace Approximation

For a Gaussian outcome with a normal leaf prior that integral is a
normal integral and has a closed form. This is what confines classic
BART to Gaussian likelihoods, or to likelihoods that can be made
Gaussian by augmentation.

For anything else the integral has no closed form. Linero
([2025](#ref-linero2025)) replaces it with a Laplace approximation:
expand the log integrand around its mode, integrate the resulting
quadratic, and use that both to score the tree and to propose the leaf
values. A Metropolis correction then makes the chain target the exact
posterior. The approximation only has to be good enough to be accepted
often; it does not have to be right.

The practical consequence is the reason the package exists. A family
needs only the log density of one observation and its first two
derivatives with respect to the predictor. Anything that can supply
those can be fitted, which is why
[`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
takes an R function and works.

The sampler recognizes three shapes of the leaf-level target and uses
the cheapest machinery each admits:

| Shape | Meaning | Cost |
|----|----|----|
| Quadratic | the target is exactly Gaussian, so the Laplace step is exact | closed form, no rejection |
| Exponential | of the form \\a\eta + b e^{r\eta}\\, as for a count or a Weibull | one mode-solve |
| General | anything else | Newton iteration |

A family that reaches the quadratic shape is several times faster than
one that does not. This is what data augmentation is for.

### Data Augmentation (`augment`)

Some likelihoods are not Gaussian but become Gaussian once an extra
latent variable is imagined and drawn alongside everything else. Three
such schemes are in use here.

**A probit model** is a Gaussian one on a latent scale that is truncated
by the observed category ([Albert and Chib 1993](#ref-albert1993)), so
drawing the latent value makes the target exactly quadratic. **A
logistic model** becomes Gaussian conditional on a Pólya-Gamma variable
([Polson et al. 2013](#ref-polson2013)), which acts as a per-observation
precision. And **a censored outcome** becomes an ordinary one once the
unobserved value is imputed from its truncated distribution.

`augment = TRUE`, which is the default wherever a scheme exists, turns
these on. Note that they change the sampler and not the model: the
posterior is the same, reached faster.
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
says which families have one.

### Nuisance Parameters

Anything in the density that is not a function of \\x\\ is drawn once
per sweep, after the trees. A residual standard deviation, a negative
binomial dispersion, the cutpoints of an ordinal model, and the baseline
hazard of a proportional hazards model are all handled this way. They
come back in `fit$aux` and are summarized by
[`summary()`](https://rdrr.io/r/base/summary.html).

A likelihood written with
[`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
can have them too, by naming them in `aux_names`. There is no prior
argument: a parameter with a restricted range is handled by writing the
transform into the density, exactly as it would be for a real predictor.

### The Dirichlet Process Mixture

Two families,
[`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
and
[`dpm_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
do not assume a shape for the error at all. They model it as a mixture
of normals with a Dirichlet process prior on the mixing distribution
([Escobar and West 1995](#ref-escobar1995); [George et al.
2019](#ref-george2019)), which lets the number of components grow with
the data. The mixture is constrained to have mean zero so that the
forest still carries the conditional mean.

Conditional on which component each observation belongs to, the model is
Gaussian and the target is quadratic, so the extra flexibility is cheap.
[`error_density()`](https://ngreifer.github.io/bartisan/reference/error_density.md)
reports the fitted error distribution, which is the diagnostic the
family exists to provide. Because the mixture can collapse to a single
component, it costs little when a single normal was right.

## Practical Matters

### Missing Predictor Values

Missing predictor values are handled natively, by treating missingness
as something to split on ([Twala et al. 2008](#ref-twala2008)). A rule
can send missing values left, send them right, or split on missingness
itself, and which of these is used is part of the posterior. No
imputation is required and rows are not dropped.

``` r

d_miss <- train
d_miss$x1[sample(n, 60)] <- NA

fit_miss <- bartisan(y ~ . - eta, data = d_miss, family = gaussian(),
                     control = ctrl)

sqrt(mean((predict(fit_miss, newdata = test) - test$eta)^2))
#> [1] 0.662
```

The fit uses the incomplete variable rather than discarding it. Note
that this is a prediction tool and not a substitute for thinking about
why values are missing: what it estimates is the mean of the response
given the predictors and the pattern of missingness, which is the
quantity prediction calls for, and which is not the regression or causal
effect defined on complete data that multiple imputation targets.

### Several Chains (`chains`)

`chains` runs the sampler more than once from different starting points.
The default is one, which is enough for estimates. Running several is
what makes a between-chain comparison possible, which is the only way to
see whether the sampler has converged, and
[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
is where that comparison is computed. It is not part of the fit: the
per-observation statistics cost more than the sampling does, so they run
when they are asked for.

``` r

fit_chains <- bartisan(y ~ . - eta, data = train, family = gaussian(),
                       chains = 4, control = ctrl)

diagnose(fit_chains)$table
#>                              quantity rhat rhat_late ess_bulk ess_tail ess_frac
#> 1                              loglik 1.39      1.27     8.87     52.0  0.00739
#> 2                           aux.sigma 1.06      1.04    71.06    488.2  0.05922
#> 3                          splits.eta 1.63      1.89     6.79     23.6  0.00566
#> 4 eta.eta (average over observations) 1.00      1.01  1031.88   1017.2  0.85990
#> 5  eta.eta (worst 5% of observations) 1.35      1.46     9.47     31.3  0.00789
#>   rhat_bad late_bad
#> 1        1        1
#> 2        1        1
#> 3        1        1
#> 4        0        0
#> 5        1        1
```

This table is read selectively. `aux.sigma` and `loglik` are close to 1,
which is what convergence looks like, and they respond to longer chains
in the usual way.

The `eta` row does not, and it is worth knowing why before it causes
alarm. It is the worst value over all 400 observations, so it is a
maximum by construction, and forests mix slowly on their fitted values:
two chains can visit quite different sets of trees that imply very
similar predictions, which inflates a between-chain statistic without
meaning that the answers disagree. Adding draws does not reliably bring
it down. Values above 1.05 on this row are ordinary for a BART fit, and
are characteristic of the method rather than of this implementation.

What to check instead is that the quantities to be reported are stable.
If `sigma` and the log likelihood have converged and an
[`avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
estimate is the same across chains, the fit is usable whatever the `eta`
row says.

Chains run sequentially unless a `future` plan is set, in which case
they run in parallel:

``` r

future::plan(future::multisession, workers = 4)

fit <- bartisan(y ~ ., data = d, chains = 4)
diagnose(fit)
```

The chain is the only parallel axis this sampler has, since a sweep
conditions on the one before it, and it is also what makes a convergence
diagnostic possible. Any *future* backend works, including *mirai*’s,
and one [`set.seed()`](https://rdrr.io/r/base/Random.html) reproduces
the whole run whatever the backend.

[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
is the one call for whether the fit can be reported. It computes
split-R-hat and the bulk and tail effective sample sizes for every
scalar the sampler draws, for the fitted function over its worst 5% of
observations, and for the size of the forest itself, then says which of
them fall short and what to change. Two of its choices are worth knowing
about. It repeats R-hat on the second half of the draws alone, which is
what distinguishes a warmup that ended too early from chains that have
each settled somewhere different: discarding the early draws is what
more `num_burn` would have done, so if that fixes R-hat then warmup was
the problem, and the advice says so rather than listing everything a
reader might try. And it leaves the leaf scale out, because that
parameter mixes badly in every implementation (1.12 in *dbarts* and 1.16
in *stochtree* against 1.19 here, on the same data) and its disagreement
never reaches the fitted function.

A within-chain drift test would answer the warmup question more
directly, and is not offered, because it cannot be made to work at a
forest’s autocorrelation. Three versions were calibrated against
stationary series, where there is nothing to find by construction:
taking each half’s Monte Carlo error from that half alone fires 31% of
the time at an autocorrelation of .995 against a nominal 5%, taking it
from the whole chain holds specificity under 3% but then misses a
six-standard-deviation trend three times in four, and batch means catch
everything and fire 92% of the time on a chain that has converged.
Re-running a statistic that is already calibrated sidesteps the choice.

[`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md)
covers all of this from the reader’s side, including which statistic to
look at when the fitted values are the estimand.

### Choosing the Settings

Most of
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
does not need to be touched. In rough order of how often changing one is
worth it:

| Setting | When to change it |
|----|----|
| `num_burn`, `num_draws` | when `rhat` says the chain has not converged |
| `num_trees` | more for a complex function and a large sample; fewer to speed up |
| `gate` | `"hard"` when the truth really is a step function, or for speed |
| `sparsity` | `FALSE` to recover classic BART, which is rarely the goal |
| `k` | to shrink harder toward the mean, in a very small sample |

Everything else exists so that the checks in the test suite can be
written, or because a specific family needs it.

The single decision that matters more than all of these is the family.
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
covers it.

## Speed

A non-conjugate sampler costs more than a conjugate one, and it is worth
being precise about how much, and about what the extra buys. Cost per
sweep is roughly linear in the number of observations and in the number
of trees, but the family matters more than either: one that reaches the
quadratic target is several times faster than one that does not, which
is why
[`lognormal_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
fits in about a second where
[`weibull_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
takes eight on the same data.

### The Benchmark

From `_dev/benchmark.Rmd`, which is reproducible: the Friedman function,
\\n = 1000\\, 10 predictors, 50 trees, 1000 warmup and 1000 saved draws,
scored against the true regression function on a held-out thousand, best
of three runs. Every package fits training data only and predicts
afterward, so what is timed is the sampler rather than the sampler plus
a thousand test evaluations, and *bartisan* names `family = gaussian()`
rather than taking the default.

| Task | Package | Time | Effective sample size | Held-out RMSE |
|----|----|----|----|----|
| Gaussian | `dbarts` | **0.28 s** | 24 | 0.225 |
|  | bartisan, hard rules | 0.91 s | 22 | 0.193 |
|  | `stochtree` | 1.19 s | 20 | 0.245 |
|  | `BART::wbart` | 1.93 s | 23 | 0.231 |
|  | `bartMachine` | 2.08 s | (not reported) | 0.234 |
|  | bartisan, soft, smoothstep gate (default) | 2.89 s | 50 | 0.130 |
|  | bartisan, soft, smootherstep gate | 2.98 s | 56 | **0.126** |
|  | bartisan, soft, logistic gate | 4.21 s | 54 | 0.130 |
| probit | `dbarts` | **0.36 s** | 51 | 0.113 |
|  | bartisan, hard rules | 1.09 s | 36 | 0.107 |
|  | `stochtree` | 1.88 s | 34 | 0.122 |
|  | `BART::pbart` | 2.06 s | 39 | 0.117 |
|  | bartisan, soft, logistic gate | 3.10 s | 64 | 0.086 |
|  | bartisan, `augment = FALSE` | 42.1 s | 73 | **0.085** |
| logit | bartisan, soft, logistic gate | **3.97 s** | 179 | **0.082** |
|  | bartisan, `augment = FALSE` | 25.0 s | 253 | 0.088 |
|  | `BART::lbart` | 48.9 s | 69 | 0.108 |
| ordinal | bartisan, hard, logit | **1.87 s** | 44 | (not comparable) |
|  | bartisan, hard, probit | 2.41 s | 44 | (not comparable) |
|  | bartisan, soft, logit | 3.92 s | 52 | (not comparable) |
|  | bartisan, soft, probit | 5.15 s | 100 | (not comparable) |
|  | bartisan, hard, cloglog | 6.42 s | 36 | (not comparable) |
|  | `stochtree` (cloglog) | 8.10 s | 27 | (not comparable) |
|  | bartisan, hard, logit, `augment = FALSE` | 32.1 s | 47 | (not comparable) |
|  | bartisan, hard, cloglog, `augment = FALSE` | 38.2 s | 33 | (not comparable) |
|  | bartisan, hard, probit, `augment = FALSE` | 51.7 s | 55 | (not comparable) |

With hard rules *bartisan* is **about three times slower than `dbarts`**
on the two tasks `dbarts` supports (0.91 s against 0.28 s on the
Gaussian task and 1.09 s against 0.36 s on probit), at comparable mixing
and better accuracy, and it is faster *and* more accurate than every
other package here on both. On the complementary log-log ordinal model,
the one task `stochtree` supports and `dbarts` does not, *bartisan* is
the faster of the two, and on a logit link it is 12 times faster than
`BART::lbart` and mixes two and a half times better. Soft rules, which
are the default, cost a further three times and cut the held-out error
by about a third, which makes the default configuration the most
accurate fit in the table.

The ordinal rows have no RMSE because the three links put the additive
predictor on different scales, so the numbers would compare scales
rather than fits; the timings and the mixing are comparable. The
remaining families are the expensive ones (2.1 s for a log-logistic
accelerated failure time model, 5.1 s for a Poisson, 10.0 s for a
negative binomial, and 10.5 s for a Gamma), which is the price of the
general machinery where no rewriting is available.

**Two earlier versions of this table were not measuring what they
claimed**, and both errors are worth stating because they cut in
opposite directions and the net was a gap reported as 1.8x that is
really about 3x. The first is that *bartisan*’s Gaussian rows called
`bartisan(y ~ ., data = dtr)` with no `family`; a numeric response with
no family reaches
[`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
a Dirichlet process mixture for the error distribution rather than a
Gaussian, so those rows timed a strictly more expensive model than the
three packages they were compared against, and their accuracy column
described a different fit. Naming the family is worth about 35%. The
second, and the larger, is that the competing packages were called as
`bart(xtr, ytr, xte, ...)`, which evaluates a thousand test points at
every draw *inside* the timed call, while *bartisan*’s
[`predict()`](https://rdrr.io/r/stats/predict.html) ran after the timer
stopped; on this task that is 0.28 s against 0.47 s for `dbarts`, so the
comparison charged `dbarts` for two thirds again as much work as
*bartisan* was doing. `_dev/parity-benchmark.R` isolates each of these,
along with two settings that turned out not to matter: *bartisan*’s
Dirichlet sparsity prior and its drawn leaf scale, neither of which
`dbarts` has, are both free to within measurement noise.

### Where the Time Goes

With **hard rules**, what is left against `dbarts` is tree bookkeeping
rather than the family. A birth, death, or change move makes five or six
passes over the node it touches, materializing the predictor with the
node’s own contribution removed, fitting the proposal, splitting the
node’s membership, recomputing the child weights, and committing,
against about two passes for a sampler that only ever needs a sum of
residuals and a count. Nothing in that gap is the likelihood: the
acceptance ratio *bartisan* forms for a Gaussian response is
algebraically the same marginal-likelihood ratio `dbarts` forms, because
a quadratic target makes the Laplace approximation exact.

Reading `dbarts`’s source says where the rest of it is, and all three of
its techniques are things soft rules give up. Its nodes hold a pointer
into **one shared index array** rather than a vector of their own, so a
split is an in-place partition and costs no allocation and no copy; each
node caches its **sufficient statistics** (a mean and an effective
count) computed in the same recursion as that partition, so the marginal
likelihood reads two doubles and never touches the data; and the
predictors are pre-discretized to **small integer cutpoint codes**, so
the split test is an integer compare over one or two bytes instead of a
gather over doubles, and the partition itself is hand-vectorized with
runtime dispatch to AVX2, SSE4.1, SSE2, or NEON.

None of the three survives a soft rule, where an observation reaches
both children with a weight rather than going to one of them: there is
no partition to do in place, the leaf target is not a function of two
summaries, and the gate needs the covariate’s value rather than which
side of a cutpoint it fell. They are all available under
`gate = "hard"`, but only by specializing a code path that soft and hard
rules currently share, which is the trade this package has made
deliberately, since the soft default is the configuration that wins on
accuracy. `BART` has none of the three and reaches every observation on
every move, which is most of why it is six times slower than `dbarts`
here.

What has actually paid, in order, is worth knowing before reaching for
the obvious ideas. **Sizing the children’s index vectors and filling
them by index, instead of `push_back`,** is worth 20% of a soft fit on
its own, being four vectors’ worth of per-element capacity tests in the
innermost loop of the sampler. **Recycling tree nodes** rather than
allocating two per birth proposal, about two thirds of which are
rejected, is worth 5%. And **folding the child weights into the split
that produces them** removes a whole pass over the same gates.

Two things that looked as though they should help did not. A hard tree
used to carry a membership weight per observation per node, every one
exactly 1.0, and the leaf-level sums used to write a node’s predictors
into a buffer and read them back. Removing both is bit-identical
(multiplying by 1.0 is exact) and together they are worth about 10%,
against a predicted factor of several. Memory bandwidth was never the
constraint: a node holds a few thousand observations at most, so the
buffer stays in L1. Nor is dropping the remaining materialized base
array worth much, measured directly at about 2%, and it would trade a
sequential read for a gather.

With **soft rules**, half the cost is a single move. The bandwidth of
each tree is a parameter with a Metropolis step per tree per sweep, and
every attempt rebuilds every membership weight in the tree. A rejected
attempt, which is 58% of them, used to rebuild a second time to get back
where it started; it is now rolled back from a snapshot, and a tree with
no splits skips the move entirely and draws its bandwidth from the
prior, which is its exact full conditional. The other half is that a
logistic gate never saturates: dropping a weight needs it below 1e-10,
which needs the observation 23 bandwidths from the cutpoint, further
than the whole unit interval, so every observation reaches every leaf
and a pass over a node covers 2.5 times the sample.

`gate = "smoothstep"`, the default, and `gate = "smootherstep"` are the
fix for the second half, though not for the reason they were built. At a
bandwidth wide enough that a bounded gate truncates *nothing*, it is
still 1.45 times faster than the logistic: what a bounded gate saves is
the [`exp()`](https://rdrr.io/r/base/Log.html), not the work on the far
side of the cutpoint. So the two bounded gates come out within noise of
each other, and the choice between them is about smoothness (one
derivative against two) rather than speed. Fixing the bandwidth outright
is the fix for the first half, and is *not* the default, because on a
function with jumps it more than doubles the error; letting the rules
sharpen toward hard ones is what the move is for.

The passes over the data are no longer part of it. A leaf value enters
the additive predictor linearly, so the log target over a leaf inherits
the shape of the log density, and where that shape is known, one pass
over a node determines the whole function and everything after it is
arithmetic. Two shapes are known. A *quadratic* log density (a Gaussian
response, and anything `augment` rewrites into one) makes the Laplace
approximation the conditional posterior exactly. An *exponential* one,
\\a\eta + b e^{r\eta}\\, covers the Poisson, the gamma, the Weibull
survival model, the augmented negative binomial, and
[`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)’s
log-scale forest under hard rules; there the fit still has to be
iterated, but on three numbers rather than on the data. The rate \\r\\
matters, since it is what brings in the last two, at \\-1/\sigma\\ and
\\-2\\. Either way a birth move goes from six passes over the node to
two.

### What Data Augmentation Buys (`augment`)

What closes most of the first factor is a data augmentation that makes
the conditional Gaussian, since then the Laplace approximation is exact
rather than approximate. `augment` does this: Albert and Chib
([1993](#ref-albert1993)) for a probit link and for an ordinal probit,
Pólya-Gamma augmentation ([Polson et al. 2013](#ref-polson2013)) for a
logit one, an ordinal logit, or a multinomial, and imputation of the
censored failure times for the log-normal and log-logistic survival
models. Every one of them trades speed for mixing, so the ratio to judge
is effective samples per second, and it differs a great deal by family:

| Family | Speed | Effective sample size | ESS per second |
|----|----|----|----|
| [`lognormal_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), hard rules | 29x | 0.87x | **21x** |
| `ordinal("probit")`, hard rules | 30x | 0.76x | **23x** |
| [`lognormal_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), soft rules | 20x | 1.05x | **15x** |
| [`loglogistic_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), hard rules | 12x | 0.88x | **10x** |
| [`loglogistic_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), soft rules | 8.9x | 0.71x | **6.3x** |
| `ordinal("logit")`, hard rules | 15x | 0.94x | **14x** |
| `ordinal("probit")`, soft rules | 14x | 0.87x | **12x** |
| `ordinal("logit")`, soft rules | 7.1x | 0.76x | **5.4x** |
| `binomial("probit")` | 5.7x | 0.66x | **3.8x** |
| `binomial("logit")` | 2.6x | 0.81x | **2.1x** |
| [`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), soft rules | 9.3x | 1.09x | **10.1x** |
| [`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), hard rules | 14.5x | 0.66x | **9.6x** |
| [`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), hard rules | 7.1x | 1.42x | **10.1x** |
| [`zi_negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), hard rules | 9.4x | 0.84x | **7.9x** |
| [`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), soft rules | 4.6x | 0.85x | **3.9x** |
| [`negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | 0.8x | 0.56x | 0.5x |

`augment` is on by default and covers the binomial, ordinal,
multinomial, zero-inflated, and survival families, plus the negative
binomial when the rules are hard; there it is written as a Poisson whose
rate comes from a gamma, which needs one gamma draw per observation
rather than a Pólya-Gamma one.

The survival families are the clearest case of what the rewriting is
for. Right-censoring is what makes their likelihood expensive: an
observed failure contributes a density in the additive predictor and a
censored one contributes a survival function, and the two have different
shapes, so nothing about the target is exploitable. Imputing each
censored failure time above its censoring time replaces the survival
term with a density, and then every observation contributes the same
quadratic shape.
[`weibull_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is absent from the table because it needs none of this, its likelihood
being already exponential in the sense above, censoring included, which
also makes it the slowest of the three at the default soft gate, where
that form cannot be used.

The zero-inflated families are the odd ones out, because what blocks
them is a mixture rather than a link: the zero contributes \\\log\[\pi +
(1 - \pi) P_0\]\\, so neither the count predictor nor the inflation
predictor has a shape to exploit. Introducing the indicator of *which*
component produced each observation separates them, and the count forest
then sees a plain Poisson or negative binomial, and the inflation forest
a Bernoulli logistic likelihood that Pólya-Gamma handles.

The ordinal probit is the largest of these because the target it
replaces is the most expensive one in the package: two cumulative-normal
evaluations per observation per pass, and no exploitable shape, so every
trial value of a leaf parameter costs its own pass. It comes with one
subtlety worth naming. The cutpoints are *not* drawn from their
conditional given the latent normals, which would be uniform between the
two order statistics that bracket each one, pinning them to an interval
of width \\O(1/n)\\ and mixing worse the more data there is. They are
drawn from the ordinal likelihood with the latent normals integrated
out, and the latent normals are redrawn immediately afterward: a
partially collapsed Gibbs sampler ([van Dyk and Park
2008](#ref-vandyk2008)), which is the standard remedy ([Cowles
1996](#ref-cowles1996)).

The cumulative **logit** has no latent normal, its latent variable being
logistic, and the usual route to one goes through the Kolmogorov-Smirnov
density, which needs a sampler of its own. It does not have to. Polson
et al. ([2013](#ref-polson2013))’s Theorem 1 at \\a = 1\\, \\b = 2\\
says that the standard logistic density is \\\frac{1}{4} E\[\exp(-w x^2
/ 2)\]\\ with \\w \sim \mathrm{PG}(2, 0)\\, so a logistic residual is a
normal whose precision is Pólya-Gamma; and because \\\mathrm{PG}(b, c)\\
is \\\mathrm{PG}(b, 0)\\ tilted by \\\exp(-c^2 w / 2)\\, the precision’s
conditional given a residual \\r\\ is exactly \\\mathrm{PG}(2, \lvert r
\rvert)\\. That is an integer-parameter draw, which the exact Devroye
sampler already in the package covers, so the logit link gets the same
quadratic target with nothing approximate anywhere. The Poisson and
gamma families need no rewriting, their targets being already in a shape
the sampler collapses to a single pass.

### Unoptimized Builds

One warning is worth heeding: an unoptimized build of the compiled code
is 5 to 20 times slower and looks no different. `devtools::load_all()`
and `devtools::install()` compile without optimization and leave the
object files for a later `R CMD INSTALL` to reuse.
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
warns once per session when it detects one, and
`bartisan:::.bartisan_optimized()` reports the state directly.

## Correctness

Three checks are in the test suite and worth knowing about.

Shrinking the prior weights to nothing makes the likelihood constant, so
the posterior collapses to the tree prior. The sampled trees then
reproduce the prior’s branching process to within a fraction of a
percent, which is a direct check that the birth, death, and change moves
satisfy detailed balance.

Every family’s analytic score is checked against a central difference of
its own log density, which is an independent route to the same quantity.
This is the check that matters for the sampler, since a wrong score
sends the Fisher-scoring proposal to the wrong place.

Across replicate datasets at the default settings, 95% credible
intervals for the additive predictor cover the truth .95 (Gaussian), .91
(binomial), .96 (Poisson), and .96 (gamma) of the time. The binomial
case is the weakest because a binary response carries the least
information per observation. The shortfall tracks the ratio of absolute
bias to posterior standard deviation, which sits near .8: an interval
centered on a shrunken estimate loses coverage in proportion to how far
that center sits from the truth. Running more or longer chains does not
help, which rules out mixing as the cause at this scale. This is in line
with what is reported for BART generally, so pointwise intervals for the
regression function are best treated as approximate.

## Comparison With Other BART Packages

Checked against the installed versions of each package rather than from
memory: *dbarts* 0.9.34, *BART* 2.9.10, *flexBART* 2.0.3, *SoftBart*
1.0.3, *bartMachine* 1.4.2, and *stochtree* 0.4.5. A dash means the
package does not offer the feature, not that it fits it badly.

|  | bartisan | dbarts | BART | flexBART | SoftBart | bartMachine | stochtree |
|----|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| **Responses** |  |  |  |  |  |  |  |
| Gaussian | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Binary, probit | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Binary, logit | ✓ | — | ✓ | ✓ | — | — | — |
| Count: Poisson, negative binomial | ✓ | — | — | — | — | — | — |
| Gamma, beta, ordered beta | ✓ | — | — | — | — | — | — |
| Ordinal | ✓ 3 links | — | — | — | — | — | ✓ cloglog |
| Multinomial | ✓ logit, probit | — | ✓ **separate binary fits** | — | — | — | — |
| Zero-inflated counts | ✓ | — | — | — | — | — | — |
| Survival | ✓ 3 AFT, PH | — | ✓ **+ recurrent, competing risks** | — | — | — | — |
| Heteroskedastic | ✓ | — | — | ✓ | — | — | ✓ |
| Error distribution itself modeled | ✓ | — | — | — | — | — | — |
| A likelihood written by the user | ✓ | — | — | — | — | — | — |
| **Rules and priors** |  |  |  |  |  |  |  |
| Hard decision rules | ✓ | ✓ | ✓ | ✓ | — | ✓ | ✓ |
| Soft decision rules | ✓ 4 gates | — | — | — | ✓ | — | — |
| Dirichlet sparsity prior (DART) | ✓ | — | ✓ | ✓ | ✓ | — | — |
| Splitting weights fixed by the user | ✓ | — | — | — | ✓ groups | ✓ | — |
| Categorical splits on level subsets | ✓ | — | — | ✓ **+ nested, network** | — | — | — |
| Missing predictors, no imputation | ✓ default | — | — | — | — | ✓ | — |
| **Structure** |  |  |  |  |  |  |  |
| Random intercepts | ✓ | ✓ | — | — | — | — | ✓ |
| Varying coefficients | ✓ | — | — | ✓ | ✓ | — | ✓ |
| Treatment-effect (BCF) structure | ✓ | — | — | — | — | — | ✓ |
| Formula interface | ✓ | ✓ | — | ✓ | — | — | — |
| **Running it** |  |  |  |  |  |  |  |
| Several chains | ✓ | ✓ | ✓ | ✓ | — | ✓ | ✓ |
| Threads inside one chain | — | ✓ | ✓ | — | — | ✓ | ✓ |
| Grow-from-root warm start | — | — | — | — | — | — | ✓ |
| Cross-validation over hyperparameters | — | ✓ | — | — | — | ✓ | — |
| Save and reload a fitted model | RDS | — | — | — | — | ✓ | ✓ JSON |
| **Reading the fit** |  |  |  |  |  |  |  |
| [`predict()`](https://rdrr.io/r/stats/predict.html) on new data | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Convergence diagnostics built in | ✓ | — | ✓ | — | — | ✓ | — |
| Variable importance | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — |
| Formal variable-selection test | ✎ | — | ✓ | — | ✓ | ✓ | — |
| Partial dependence | ✎ | ✓ | — | — | ✓ | ✓ | — |
| Interaction detection | ✎ | — | — | — | — | ✓ | — |
| Counterfactual estimands with intervals | ✓ | — | — | — | — | — | — |
| Cross-validated model comparison | ✓ | — | — | — | — | — | — |

**On the multinomial row, *BART*’s support is a set of binary fits
rather than a multinomial model**, which is worth knowing before the
checkmark is read as a like-for-like. `BART::mbart2()` fits \\K\\
independent one-vs-rest binary probit or logit BARTs, each to the
indicator `(y == h)`, and then normalizes their latent predictors
afterward, so no forest is ever fitted to a multinomial likelihood.
`BART::mbart()` instead fits \\K - 1\\ binary models to nested subsets,
the \\h\\th to `y == h` among the observations with `y >= h`, and
multiplies them as a continuation-ratio. That second one is an exact
factorization of the multinomial mass function, so it is a coherent
model, but each conditional carries its own independent forest and
prior, and the factorization runs over the categories in sorted order,
so the prior it puts on a probability vector is not exchangeable in the
categories. *bartisan*’s
[`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
couples the forests through the likelihood itself instead: conditional
on the others, the likelihood for category \\j\\ is exactly
binomial-logistic in \\\eta_j - \log C_j\\ with \\C_j\\ summing
\\e^{\eta}\\ over the rest, which is what the Pólya-Gamma augmentation
exploits and what keeps the prior symmetric in the categories ([Murray
2021](#ref-murray2021)).

**✎ means a helper package covers it, not this one.** All three of those
rows are things *bartMachine* does natively and *bartisan* does through
*marginaleffects*, because a fit works with it and every estimand there
is computed by pushing the draws through:
[`plot_predictions()`](https://rdrr.io/pkg/marginaleffects/man/plot_predictions.html)
is a partial dependence plot, and
`avg_comparisons(variables = "x", by = "z")` and `hypotheses(~pairwise)`
are interaction detection, neither of which needs code here. What is
genuinely missing is a *formal* variable-selection test: this package
reports split counts and the share of draws that used a predictor, which
is not a test, where *bartMachine* permutes the response and *SoftBart*
reports posterior inclusion probabilities.

The other columns worth reading as gaps rather than as differences are
these. There are **no threads inside a chain**, so a single chain is
single-core here where *dbarts*, *BART*, *bartMachine*, and *stochtree*
all use several; chains do run in parallel, which is the cheaper win,
but it does not help one chain. There is **no cross-validation** over
`k`, `num_trees`, and the tree prior, which `dbarts::xbart()` and
`bartMachine::bartMachineCV()` both automate. There is **no
grow-from-root warm start**, which is *stochtree*’s way of shortening
burn-in. And there is **no JSON serialization**, so a fit round-trips
through [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and not into
another language. On survival, *BART* is ahead: recurrent events and
competing risks are there and here they are not.

Three more packages are single-purpose rather than general, so they are
not columns above: [*bcf*](https://CRAN.R-project.org/package=bcf) fits
Bayesian causal forests only,
[*VCBART*](https://github.com/skdeshpande91/VCBART) varying-coefficient
models only, and *bartCause* wraps *dbarts* for causal estimands.

## Notes and Limitations

Soft rules cost more per iteration than hard ones, because a leaf
touches every observation rather than only those in its cell. Negligible
weights are pruned, which holds the gap to roughly a factor of three
rather than the factor of \\B\\ it would otherwise be.

Absolute timings are machine- and load-dependent, so these are
indicative rather than exact; repeated runs on the same laptop varied by
about 20%. As an anchor, a Gaussian response with 10 predictors, 50 soft
trees, and the default 1000 warmup plus 1000 saved draws takes roughly 6
seconds at \\n = 500\\ and roughly a minute at \\n = 5000\\, measured at
steady state on one core of an M-series Mac. Hard rules are about 3
times faster, and doubling the tree count roughly doubles the cost.

The ratios are stabler than the absolute times. Cost relative to a
Gaussian fit at the same size is set by how expensive the family’s log
density and its derivatives are:

| Family | Relative cost |
|----|----|
| [`gaussian()`](https://rdrr.io/r/stats/family.html) | 1x |
| `Gamma("log")`, [`poisson()`](https://rdrr.io/r/stats/family.html) | ~2x |
| [`negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), [`binomial()`](https://rdrr.io/r/stats/family.html) | 2–3x |
| [`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | ~2x with a latent variable, ~25x without |
| [`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | ~10x |
| [`ordbeta()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | ~25x |
| [`tweedie()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | ~2x, half again as much with the power drawn |

[`tweedie()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
earns a note because it looks as though it should be far worse: its
density has no closed form at a positive response, and normalizing it
takes an infinite series. Measured against a Gaussian fit of the same
size it costs 5.5 times as much where `Gamma("log")` costs 5.2, so the
series is not visible in the total. The reason is where the series sits.
Writing the log density in exponential-dispersion form separates it into
a part that moves with the predictor, which is closed form, and a
normalizing term that does not depend on the predictor at all; the
second goes in the eta-free part, which is evaluated once per sweep
rather than at every leaf, and cancels from every acceptance ratio in
between. Its length depends on the response and the dispersion but never
on the mean, so it does not grow as the forest moves. Drawing `power`
costs about a third again, because the slice sampler has to re-sum the
series at each candidate value and cannot use the table of log-gammas
that a fixed power allows.

[`ordbeta()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is the expensive one because every evaluation of its log density needs
two log-gamma calls that no amount of restructuring removes.
[`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
spans a wide range because all three of its links have a latent-variable
representation that the sampler uses by default (a normal for the
probit, a normal with a Pólya-Gamma precision for the logit, and an
exponential waiting time for the complementary log-log), and
`augment = FALSE` turns them off.

An ordinal model is identified only up to a common shift of its
cutpoints and its predictor. With three or more categories the draws are
reported in the chart where **the predictor has mean zero over the
fitted sample and every cutpoint is free**, which makes the cutpoints
readable as category boundaries. With two categories the single boundary
is folded into the intercept instead, so the fit is on the same scale as
[`binomial()`](https://rdrr.io/r/stats/family.html).

That is [`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html)’s chart
with its predictors centered. `polr()` identifies the location by
leaving the intercept out of the design matrix rather than by centering,
so its `zeta` is shifted by the mean of its own linear predictor. Either
of these lines puts the two side by side, and on a linear truth they
agree to Monte Carlo error:

``` r

p <- MASS::polr(y ~ x1 + x2, data = d, method = "probit")
p$zeta - mean(p$lp)                      # same chart as colMeans(fit$aux)
```

Two prediction types exist for reporting an ordinal fit on a single
scale, both following *WeightIt*. `predict(type = "mean")` weights the
category probabilities by the labels read as numbers, so levels `"1"`,
`"2"`, and `"4"` give a mean between one and four, and `values` says
what the categories are worth when the labels are not numbers.
`predict(type = "stdlv")` divides the predictor by the standard
deviation of the latent variable it indexes, which puts fits with
different links, or different amounts of signal, on one scale; it works
for [`binomial()`](https://rdrr.io/r/stats/family.html) too, since a
binary response is the same construction with one threshold. Note that
the complementary log-log error enters the two families with opposite
signs, and
[`?predict.bartisan_fit`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
says why.

The leaf scale is raised from near zero over the first quarter of
warmup. Linero ([2025](#ref-linero2025)) describes this as essential:
started at its full value, the sampler can settle early into a poor
configuration and fail to move. It is controlled by `sigma_mu_ramp` in
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md).

Nuisance parameters that are weakly identified, such as the negative
binomial dispersion or the shape of a survival model, mix more slowly
than the regression function, and the data often do not pin them down.
Their draws are worth checking before they are interpreted.

The leaf scale is drawn under a half-Cauchy prior, which has no upper
bound. When the predictors nearly separate a binary response the
likelihood rewards an unbounded predictor and that prior does not hold
the scale down: in a fully separated example the drawn scale wandered
between 3 and 9 times its prior median over 1600 draws without settling.
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
warns when the scale settles more than five times above its prior
median, and `update_sigma_mu = FALSE` pins it.

A link supplied from R costs a call into the interpreter for every leaf
the sampler visits, and
[`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
costs three unless the derivatives are supplied. Both are usable at the
sizes a BART model is usually fit at, and both are slower than the
compiled equivalent.

## Further Reading

Hill et al. ([2020](#ref-hill2020)) is a review of BART covering the
model, its extensions and its applications. Chipman et al.
([2010](#ref-chipman2010)) is the original. Linero and Yang
([2018](#ref-linero2018)) introduces soft rules and Linero
([2018](#ref-linero2018sparse)) the sparsity prior. Linero
([2025](#ref-linero2025)) is the sampler this package implements. For
causal inference with BART, Hill ([2011](#ref-hill2011)) is the starting
point and Hahn et al. ([2020](#ref-hahn2020)) the standard treatment of
what goes wrong when treatment assignment is strongly predicted by the
covariates. Murray ([2021](#ref-murray2021)) gives an alternative route
to count and multinomial outcomes.

## References

Albert, James H., and Siddhartha Chib. 1993. “Bayesian Analysis of
Binary and Polychotomous Response Data.” *Journal of the American
Statistical Association* 88 (422): 669–79.
<https://doi.org/10.1080/01621459.1993.10476321>.

Chipman, Hugh A., Edward I. George, and Robert E. McCulloch. 2010.
“BART: Bayesian Additive Regression Trees.” *The Annals of Applied
Statistics* 4 (1): 266–98. <https://doi.org/10.1214/09-AOAS285>.

Cowles, Mary Kathryn. 1996. “Accelerating Monte Carlo Markov Chain
Convergence for Cumulative-Link Generalized Linear Models.” *Statistics
and Computing* 6 (2): 101–11. <https://doi.org/10.1007/BF00162520>.

Escobar, Michael D., and Mike West. 1995. “Bayesian Density Estimation
and Inference Using Mixtures.” *Journal of the American Statistical
Association* 90 (430): 577–88.
<https://doi.org/10.1080/01621459.1995.10476550>.

George, Edward, Purushottam Laud, Brent Logan, Robert McCulloch, and
Rodney Sparapani. 2019. “Fully Nonparametric Bayesian Additive
Regression Trees.” In *Topics in Identification, Limited Dependent
Variables, Partial Observability, Experimentation, and Flexible
Modeling: Part b*, vol. 40B. Advances in Econometrics. Emerald
Publishing Limited. <https://doi.org/10.1108/S0731-90532019000040B006>.

Hahn, P. Richard, Jared S. Murray, and Carlos M. Carvalho. 2020.
“Bayesian Regression Tree Models for Causal Inference: Regularization,
Confounding, and Heterogeneous Effects (with Discussion).” *Bayesian
Analysis* 15 (3): 965–1056. <https://doi.org/10.1214/19-BA1195>.

Hill, Jennifer L. 2011. “Bayesian Nonparametric Modeling for Causal
Inference.” *Journal of Computational and Graphical Statistics* 20 (1):
217–40. <https://doi.org/10.1198/jcgs.2010.08162>.

Hill, Jennifer, Antonio Linero, and Jared Murray. 2020. “Bayesian
Additive Regression Trees: A Review and Look Forward.” *Annual Review of
Statistics and Its Application* 7 (1): 251–78.
<https://doi.org/10.1146/annurev-statistics-031219-041110>.

Linero, Antonio R. 2018. “Bayesian Regression Trees for High-Dimensional
Prediction and Variable Selection.” *Journal of the American Statistical
Association* 113 (522): 626–36.
<https://doi.org/10.1080/01621459.2016.1264957>.

Linero, Antonio R. 2025. “Generalized Bayesian Additive Regression Trees
Models: Beyond Conditional Conjugacy.” *Journal of the American
Statistical Association* 120 (549): 356–69.
<https://doi.org/10.1080/01621459.2024.2337156>.

Linero, Antonio R., and Yun Yang. 2018. “Bayesian Regression Tree
Ensembles That Adapt to Smoothness and Sparsity.” *Journal of the Royal
Statistical Society Series B: Statistical Methodology* 80 (5): 1087–110.
<https://doi.org/10.1111/rssb.12293>.

Murray, Jared S. 2021. “Log-Linear Bayesian Additive Regression Trees
for Multinomial Logistic and Count Regression Models.” *Journal of the
American Statistical Association* 116 (534): 756–69.
<https://doi.org/10.1080/01621459.2020.1813587>.

Polson, Nicholas G., James G. Scott, and Jesse Windle. 2013. “Bayesian
Inference for Logistic Models Using Pólya–Gamma Latent Variables.”
*Journal of the American Statistical Association* 108 (504): 1339–49.
<https://doi.org/10.1080/01621459.2013.829001>.

Twala, B. E. T. H., M. C. Jones, and D. J. Hand. 2008. “Good Methods for
Coping with Missing Data in Decision Trees.” *Pattern Recognition
Letters* 29 (7): 950–56. <https://doi.org/10.1016/j.patrec.2008.01.010>.

van Dyk, David A., and Taeyoung Park. 2008. “Partially Collapsed Gibbs
Samplers: Theory and Methods.” *Journal of the American Statistical
Association* 103 (482): 790–96.
<https://doi.org/10.1198/016214508000000409>.
