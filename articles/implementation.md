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
observations. The defaults are 50 trees and 800 draws after 200 warmup
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

In ordinary BART, the sum is the conditional mean of a Gaussian outcome,

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
says how many forests each family uses.

### The Prior

The model prior has three pieces, all with defaults chosen so that the
fit is regularized without being tuned ([Chipman et al.
2010](#ref-chipman2010)).

**Tree structure.** A node at depth \\d\\ is split rather than left as a
leaf with probability

\\\gamma (1 + d)^{-\beta}, \qquad \gamma = 0.95, \\ \beta = 2,\\

which are the `gamma` and `beta` arguments of
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md).

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
bandwidth, and right with the remaining weight. An observation near a
cutpoint reaches both of that node’s children, so the fitted function is
smooth rather than a step. As \\\tau \to 0\\ the gate becomes a step and
the hard rule is recovered.

The `gate` argument chooses \\\psi\\: `"smoothstep"` (the default),
`"smootherstep"`, `"logistic"`, or `"hard"`. How far the smoothing
reaches is what separates them. The logistic has unbounded support, so
under it every observation really does reach every leaf with some
weight. The two polynomial gates are the cumulative distribution
functions of symmetric Beta kernels and reach exactly zero and one
outside a finite window, so an observation far enough from a cutpoint
takes one side only and the subtree on the other side is never visited.
That skipped work turns out to be worth less than it sounds; what makes
them faster is that a polynomial needs no
[`exp()`](https://rdrr.io/r/base/Log.html), and they keep that advantage
even at a bandwidth wide enough to truncate nothing at all. `bandwidth`
sets the prior mean of \\\tau\\, which is drawn rather than fixed, and
it means the same amount of smoothing under every gate, the kernels
being matched on their standard deviation rather than on their width.

Soft rules cost about three times as much per iteration but are usually
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
#> 1 smoothstep     0.416     1.3
#> 2       hard     1.079     0.4
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

Without the sparsity prior, the advice to use fewer trees is essential:
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

## How the Model Is Fit

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
\\ dM\\

### The General Case, and the Laplace Approximation

For a Gaussian outcome with a normal leaf prior, that integral is a
normal integral and has a closed form. This is what confines classic
BART to Gaussian likelihoods, or to likelihoods that can be made
Gaussian by augmentation.

For anything else, the integral has no closed form. Linero
([2025](#ref-linero2025)) replaces it with a Laplace approximation:
expand the log integrand around its mode, integrate the resulting
quadratic, and use that both to score the tree and to propose the leaf
values. A Metropolis correction then makes the chain target the exact
posterior. The approximation only has to be good enough to be accepted
often; it does not have to be right.

The practical consequence is the reason the package exists. A family
needs only the log density of one observation and its first two
derivatives with respect to the predictor. Anything that can supply
those can be fit, which is why
[`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
takes an R function and works.

The sampler recognizes three shapes of the leaf-level target and uses
the cheapest machinery each admits:

| Shape | Meaning | Cost |
|----|----|----|
| Quadratic | the target is exactly Gaussian, so the Laplace step is exact | closed form, no rejection |
| Exponential | of the form \\a\eta + b e^{r\eta}\\, as for a count or a Weibull | one mode-solve |
| General | anything else | Newton iteration |

A family that lands higher in this table is several times faster than
one that does not. This is what data augmentation is for.

### Data Augmentation (`augment`)

Some likelihoods move up the table once an extra latent variable is
imagined and drawn alongside everything else. The variable is not part
of the model; it is machinery, and the posterior over everything else is
unchanged.

Several such rewritings are in use here, and they do not all arrive at
the same shape. **A probit model** is a Gaussian one on a latent scale
truncated by the observed category ([Albert and Chib
1993](#ref-albert1993)), so drawing the latent value makes the target
exactly quadratic. **A logistic model** is Gaussian conditional on a
Pólya-Gamma variable ([Polson et al. 2013](#ref-polson2013)), which acts
as a per-observation precision, and an ordinal complementary log-log
model is Gaussian conditional on an exponential waiting time. **A
censored outcome** becomes an uncensored one once the unobserved value
is imputed from its truncated distribution, which for the log-normal and
log-logistic survival families leaves a quadratic target. **A negative
binomial** is a Poisson whose rate carries a gamma prior, so drawing the
rate leaves the *exponential* shape rather than the quadratic one, which
is a smaller gain. And **a zero-inflated count** is blocked by a mixture
rather than by a link, so what is drawn there is the indicator of which
component produced each observation; that separates the two forests, and
for a zero-inflated negative binomial a second rewriting goes on top of
the first.

`augment = TRUE`, which is the default wherever a rewriting exists,
turns these on. They change the sampler and not the model: the posterior
is the same, reached faster.
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
says which families have one, and the section on speed below measures
what each is worth.

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
#> [1] 0.504
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
#> 1                              loglik 1.07     1.064    58.56    325.1  0.04880
#> 2                           aux.sigma 1.01     1.018   408.53    496.5  0.34044
#> 3                          splits.eta 1.27     1.346    12.56     18.7  0.01046
#> 4 eta.eta (average over observations) 1.00     0.999  1060.22   1172.6  0.88352
#> 5  eta.eta (worst 5% of observations) 1.44     1.607     7.98     27.4  0.00665
#>   rhat_bad late_bad
#> 1        1        1
#> 2        0        1
#> 3        1        1
#> 4        0        0
#> 5        1        1
```

This table is read selectively. `aux.sigma` and `loglik` are close to 1,
which is what convergence looks like, and they respond to longer chains
in the usual way.

The `eta` row does not, and it is worth knowing why before it causes
alarm. It summarizes the worst 5% of the observations rather than the
single worst, so it is a high percentile by construction, and forests
mix slowly on their fitted values: two chains can visit quite different
sets of trees that imply very similar predictions, which inflates a
between-chain statistic without meaning that the answers disagree.
Adding draws does not reliably bring it down. Values above 1.05 on this
row are ordinary for a BART fit, and are characteristic of the method
rather than of this implementation.

What to check instead is the quantity that will be reported, on its own
draws.
[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
takes the output of
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
for exactly this, and
[`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html)
takes anything else that can be arranged as draws by chains. A high
`eta` row is not on its own a reason to discard a fit, but neither is a
low one a reason to trust an effect read off it: a contrast can mix
badly where the function it contrasts mixes well, so it has to be
diagnosed rather than inferred.

Chains run sequentially unless a `future` plan is set, in which case
they run in parallel:

``` r

future::plan(future::multisession, workers = 4)

fit <- bartisan(y ~ ., data = d, chains = 4)
diagnose(fit)
```

The chain is the only parallel axis the *sampler* has, since a sweep
conditions on the one before it, and it is also what makes a convergence
diagnostic possible. Two things outside the sampler use a plan when one
is set: the convergence pass splits its per-observation columns across
workers, and
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
predicts the treatment levels in parallel. Any *future* backend works,
including *mirai*’s, and one
[`set.seed()`](https://rdrr.io/r/base/Random.html) reproduces the whole
run whatever the backend.

[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
is the one call for whether the fit can be reported. It computes
split-R-hat and the bulk and tail effective sample sizes for every
scalar the sampler draws, for the fitted function twice over, once
averaged across observations and once over its worst 5% of them, and for
the size of the forest itself, then says which of them fall short and
what to change. Handed the output of
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
instead of a fit, it reports on the contrast and on the two potential
outcomes it is a contrast of. Two of its choices are worth knowing
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

### What Data Augmentation Buys (`augment`)

A major contributor to *bartisan*’s speed is data augmentation. Where
the rewriting reaches a Gaussian conditional the Laplace approximation
stops being an approximation at all, which is the largest of the gains:
`augment` uses Albert and Chib ([1993](#ref-albert1993)) for a probit
link and for an ordinal probit, Pólya-Gamma augmentation ([Polson et al.
2013](#ref-polson2013)) for a logit link, an ordinal logit, or a
multinomial, and imputation of the censored failure times for the
log-normal and log-logistic survival models. The negative binomial and
the zero-inflated families are rewritten too, but to a Poisson rather
than to a Gaussian, and the two paragraphs after the table say what that
costs them.

Each of them is faster, from 1.2 times for the negative binomial under
soft rules to 31 for an ordinal probit with hard rules, with a median of
8.8. What each costs in mixing is the question, and the answer is not
that they all cost something: over 15 replicate datasets the
worst-mixing quantity’s effective sample size came out a median of 1.05
times what the direct likelihood gives, and of the seventeen cases
below, seven have an 80% interval lying entirely below 1, two lie
entirely above it, and eight include it. Effective draws per second,
which is the ratio to judge, favors augmentation in all seventeen:

| Family | Speed | ESS, worst quantity | ESS, median quantity | ESS per second |
|----|----|----|----|----|
| `ordinal("probit")`, hard rules | 31x | 0.72x \[0.54, 0.95\] | 0.73x | **22x** |
| `ordinal("probit")`, soft rules | 24x | 0.65x \[0.52, 0.84\] | 0.61x | **16x** |
| [`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), hard rules | 12x | 1.06x \[0.73, 1.48\] | 0.97x | **13x** |
| [`lognormal_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), soft rules | 11x | 1.09x \[0.62, 1.95\] | 1.05x | **12x** |
| [`loglogistic_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), hard rules | 8.8x | 1.34x \[0.72, 2.41\] | 0.96x | **12x** |
| [`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), soft rules | 8.2x | 1.43x \[1.00, 1.96\] | 0.91x | **12x** |
| `binomial("probit")` | 9.7x | 1.19x \[0.77, 2.01\] | 1.12x | **12x** |
| [`lognormal_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), hard rules | 14x | 0.74x \[0.57, 0.97\] | 0.96x | **10x** |
| `binomial("logit")` | 6.7x | 1.44x \[1.02, 1.98\] | 1.11x | **9.6x** |
| [`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), hard rules | 5.8x | 1.58x \[1.12, 2.09\] | 1.05x | **9.1x** |
| `ordinal("logit")`, soft rules | 12x | 0.57x \[0.45, 0.76\] | 0.65x | **7.1x** |
| [`zi_negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), hard rules | 6.7x | 1.04x \[0.81, 1.37\] | 0.99x | **6.9x** |
| `ordinal("logit")`, hard rules | 15x | 0.37x \[0.30, 0.46\] | 0.52x | **5.5x** |
| [`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), soft rules | 3.7x | 1.24x \[0.86, 1.79\] | 1.15x | **4.6x** |
| [`loglogistic_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), soft rules | 6.9x | 0.58x \[0.46, 0.80\] | 0.69x | **4.0x** |
| [`negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), soft rules | 1.2x | 1.11x \[0.68, 1.64\] | 0.86x | 1.3x |
| [`negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), hard rules | 1.8x | 0.61x \[0.47, 0.79\] | 0.84x | 1.1x |

The bracketed range is an 80% bootstrap interval over the replicates,
and it is in the table because without it the column invites conclusions
it cannot support: at three replicates the same quantity varied by a
median factor of 4.5 within a single cell, since the worst-quantity
effective sample size is a minimum over many quantities and so is a
high-variance thing to estimate. The median quantity is the stabler
reading, and comparing the two columns is informative on its own: those
ratios run .52 to 1.15 where the worst-quantity ones run .37 to 1.58, so
a family whose worst corner mixes noticeably better or worse under
augmentation has usually moved less than that as a whole.

`augment` is on by default for every family that has a rewriting: the
binomial, ordinal, multinomial, negative binomial, zero-inflated, and
survival families. The negative binomial’s is a Poisson whose rate comes
from a gamma, which needs one gamma draw per observation rather than a
Pólya-Gamma one, and it is the one family whose two rows differ in kind
rather than degree. Under hard rules the target reaches the exponential
form and the rewriting buys time at a real cost in mixing; under soft
rules it cannot, and the rewriting instead buys a little time at no
measurable cost in mixing. Both come out ahead on effective draws per
second, but by the least of any family here, so `augment = FALSE` is
worth trying if a negative binomial fit’s diagnostics look poor.

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

Over 40 replicate datasets at \\n = 500\\ with 10 predictors and the
default settings, 95% credible intervals for the additive predictor
cover the truth .964 of the time for a Gaussian response, .961 binomial,
.973 Poisson, and .970 gamma, each with a Monte Carlo standard error
near .006. Nominal coverage is .95, so these intervals are mildly
conservative rather than deficient, and the binomial sits .008 below the
other three, which at this replicate count is not a difference (p =
.35). The ratio of absolute bias to posterior standard deviation runs
.69 to .74, so the posterior spread is comfortably larger than the
distance from the truth, which is the same fact read a second way. Four
chains and five times the draws move neither number (p = .20 and p =
.18), so mixing is not what sets them. All of this is one
data-generating process: a mean function of moderate spread mapped to
each family’s link, which for the binomial keeps the success
probabilities informative. A response whose probabilities sit near 0 or
1 carries much less information per observation, and pointwise intervals
for the regression function are best treated as approximate on that
account rather than because of anything measured here.

## Comparison With Other BART Packages

Below we compare *bartisan*’s capabilities to those of other popular
BART R packages: *dbarts* 0.9.34, *BART* 2.9.10, *flexBART* 2.0.3,
*SoftBart* 1.0.3, *bartMachine* 1.4.2, and *stochtree* 0.4.5. A dash
means the package does not offer the feature.

|  | bartisan | dbarts | BART | flexBART | SoftBart | bartMachine | stochtree |
|----|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| **Responses** |  |  |  |  |  |  |  |
| Gaussian | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Binary, probit | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Binary, logit | ✓ | — | ✓ | ✓ | — | — | — |
| Count: Poisson, negative binomial | ✓ | — | — | — | — | — | — |
| Gamma, beta, ordered beta, Tweedie | ✓ | — | — | — | — | — | — |
| Ordinal | ✓ 3 links | — | — | — | — | — | ✓ cloglog |
| Multinomial | ✓ logit, probit | — | ✓[^1] | — | — | — | — |
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
| Formal variable-selection test | ✓ inclusion | — | ✓ permutation | — | ✓ inclusion | ✓ permutation | — |
| Partial dependence | ✓ | ✓ | — | — | ✓ | ✓ | — |
| Interaction detection | ✎ | — | — | — | — | ✓ | — |
| Counterfactual estimands with intervals | ✓ | — | — | — | — | — | — |
| ATE, ATT or per-unit treatment effects | ✓ | — | — | — | — | — | — |
| Cross-validated model comparison | ✓ | — | — | — | — | — | — |

**✎ means a helper package covers it, not this one.** Interaction
detection is the remaining one, and *bartisan* does it through
*marginaleffects*, because a fit works with it and every estimand there
is computed by pushing the draws through:
`avg_comparisons(variables = "x", by = "z")` and `hypotheses(~pairwise)`
are interaction detection, and neither needs code here. Partial
dependence used to sit in the same row and now does not:
[`partial_dependence()`](https://ngreifer.github.io/bartisan/reference/partial_dependence.md)
averages the fitted surface over the sample at each value of one or two
predictors, and [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
on a fit is the same thing drawn.
[`marginaleffects::plot_predictions()`](https://rdrr.io/pkg/marginaleffects/man/plot_predictions.html)
remains the one to reach for when the grid needs more control than that.

The variable-selection row divides differently, into two things that
both get called a test. One is the posterior inclusion probability under
the Dirichlet sparsity prior, thresholded at .5 to give the median
probability model of Barbieri and Berger ([2004](#ref-barbieri2004)).
*SoftBart* reports it as `posterior_probs()`, and it is
`colMeans(var_counts > 0)` on the splitting counts: the `prop_used`
column of
[`variable_importance()`](https://ngreifer.github.io/bartisan/reference/variable_importance.md)
is that same quantity, under a prior that is on by default here as it is
there. The other is a permutation test ([Bleich et al.
2014](#ref-bleich2014)), which refits the model on responses whose link
to the predictors has been broken and thresholds the splitting shares
against the resulting null rather than against .5. *bartMachine* offers
all three of its thresholds and *BART* offers the third one,
`BART::mc.wbart.gse()`; *bartisan* offers none of them, and that is the
actual gap.

It is a gap worth closing, because the two answer differently. Measured
over 20 replicates at \\n = 300\\ with 25 predictors, 3 of them real,
and 20 trees: on data where nothing matters, the .5 cut selected a null
predictor in every replicate, 7 of the 25 on average with the sparsity
prior and all 25 without, while the two simultaneous permutation
thresholds selected one in 1 replicate of 10 against a nominal .05. The
.5 cut is a rule for choosing a predictive submodel and behaves like
one; only the permutation null controls a false-positive rate. The cost
is not the obstacle either, since a hundred fits of a 20-tree model runs
in seconds for a cheap family. The obstacle is that the two corrections
do not compose: run the permutation test with `sparsity = TRUE` and
power collapses, to .12 for the max threshold and .37 for the SE one
against .75 and .82 with `sparsity = FALSE`, because the prior shrinks
the null distribution by the same mechanism it shrinks the observed
shares. Anything implemented here would have to fit its null with the
sparsity prior off, and say so.

The other columns worth reading as gaps rather than as differences are
these. There are no threads inside a chain, so a single chain is
single-core here where *dbarts*, *BART*, *bartMachine*, and *stochtree*
all use several; chains do run in parallel, which is the cheaper win,
but it does not help one chain. There is no cross-validation over `k`,
`num_trees`, or the tree prior, which `dbarts::xbart()` and
`bartMachine::bartMachineCV()` both automate. There is no grow-from-root
warm start, which is *stochtree*’s way of shortening burn-in. And there
is no JSON serialization, so a fit round-trips through
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and not into another
language. On survival, *BART* is ahead: recurrent events and competing
risks are there and here they are not.

Three more packages are single-purpose rather than general, so they are
not columns above: [*bcf*](https://CRAN.R-project.org/package=bcf) fits
Bayesian causal forests only,
[*VCBART*](https://github.com/skdeshpande91/VCBART) varying-coefficient
models only, and *bartCause* wraps *dbarts* for causal estimands.

## Notes and Limitations

Soft rules cost more per iteration than hard ones, because an
observation near a cutpoint has to be followed down both sides of it
rather than into one cell, so a leaf that a hard rule would skip still
has to be evaluated. The gap is about a factor of \\B\\, the number of
leaves, and the reason it comes out near three rather than near a
hundred is the tree prior: it keeps trees to about two and a half leaves
on average, so there is not much for the softening to multiply.

Absolute timings are machine- and load-dependent, so these are
indicative rather than exact; repeated runs on the same laptop varied by
about 20%. As an anchor, a Gaussian response with 10 predictors, 50 soft
trees, and the default 200 warmup plus 800 saved draws takes about a
second at \\n = 500\\ and about ten seconds at \\n = 5000\\, measured at
steady state on one core of an M-series Mac. Hard rules are about 3
times faster (2.9 times at \\n = 500\\ and 3.2 at \\n = 5000\\), and
doubling the tree count roughly doubles the cost.

The ratios are stabler than the absolute times. Cost relative to a
Gaussian fit at the same size is set by how expensive the family’s log
density and its derivatives are:

| Family | Relative cost |
|----|----|
| [`gaussian()`](https://rdrr.io/r/stats/family.html) | 1x |
| [`binomial()`](https://rdrr.io/r/stats/family.html), either link | 1.1x |
| [`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | 1.2x |
| [`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), logit or probit | 1.2x |
| [`poisson()`](https://rdrr.io/r/stats/family.html) | 3.6x |
| `Gamma("log")` | 5.6x |
| [`tweedie()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | 6.5x |
| [`negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | 6.8x |
| [`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | 7.8x |
| [`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | 8.5x |
| [`tweedie()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), with `power` drawn | 9.1x |
| [`Beta()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | 18x |
| [`ordbeta()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | 23x |

[`tweedie()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
earns a note because it looks as though it should be far worse: its
density has no closed form at a positive response, and normalizing it
takes an infinite series. Measured against a Gaussian fit of the same
size it costs 6.5 times as much where `Gamma("log")` costs 5.6, so the
series is not visible in the total. The reason is where the series sits.
Writing the log density in exponential-dispersion form separates it into
a part that moves with the predictor, which is closed form, and a
normalizing term that does not depend on the predictor at all; the
second goes in the eta-free part, which is evaluated once per sweep
rather than at every leaf, and cancels from every acceptance ratio in
between. Its length depends on the response and the dispersion but never
on the mean, so it does not grow as the forest moves. Drawing `power`
costs about 40% again, because the slice sampler has to re-sum the
series at each candidate value and cannot use the table of log-gammas
that a fixed power allows.

[`ordbeta()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
and
[`Beta()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
are the expensive ones because every evaluation of their log density
needs two log-gamma calls that no amount of restructuring removes.
[`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is at the other end, and only because of augmentation: all three of its
links have a latent-variable representation that the sampler uses by
default (a normal for the probit, a normal with a Pólya-Gamma precision
for the logit, and an exponential waiting time for the complementary
log-log), which is what puts it beside
[`binomial()`](https://rdrr.io/r/stats/family.html) here.
`augment = FALSE` turns them off, and the table above is then no guide
at all; the augmentation table earlier in this section is.

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

Barbieri, Maria Maddalena, and James O. Berger. 2004. “Optimal
Predictive Model Selection.” *The Annals of Statistics* 32 (3).
<https://doi.org/10.1214/009053604000000238>.

Bleich, Justin, Adam Kapelner, Edward I. George, and Shane T. Jensen.
2014. “Variable Selection for BART: An Application to Gene Regulation.”
*The Annals of Applied Statistics* 8 (3).
<https://doi.org/10.1214/14-AOAS755>.

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

[^1]: *BART*’s support is a set of binary fits rather than a multinomial
    model, which is worth knowing before the checkmark is read as a
    like-for-like. `BART::mbart2()` fits \\K\\ independent one-vs-rest
    binary probit or logit BARTs, each to the indicator `(y == h)`, and
    then normalizes their latent predictors afterward, so no forest is
    ever fitted to a multinomial likelihood. `BART::mbart()` instead
    fits \\K - 1\\ binary models to nested subsets, the \\h\\th to
    `y == h` among the observations with `y >= h`, and multiplies them
    as a continuation-ratio. That second one is an exact factorization
    of the multinomial mass function, so it is a coherent model, but
    each conditional carries its own independent forest and prior, and
    the factorization runs over the categories in sorted order, so the
    prior it puts on a probability vector is not exchangeable in the
    categories. *bartisan*’s
    [`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
    couples the forests through the likelihood itself instead:
    conditional on the others, the likelihood for category \\j\\ is
    exactly binomial-logistic in \\\eta_j - \log C_j\\ with \\C_j\\
    summing \\e^{\eta}\\ over the rest, which is what the Pólya-Gamma
    augmentation exploits and what keeps the prior symmetric in the
    categories ([Murray 2021](#ref-murray2021)).
