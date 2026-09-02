# Generalized BART with bartisan

``` r

library(bartisan)

set.seed(2026)
```

This vignette describes what the model is and how it is fitted. It is
the reference document rather than the entry point:
[`vignette("bartisan")`](https://ngreifer.github.io/bartisan/articles/bartisan.md)
shows how to use the package on a real analysis, and
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
covers choosing a likelihood. Read this one when you want to know what
the fit is doing, why the defaults are what they are, or what an
advanced setting changes.

Fits below use a deliberately small chain so the vignette builds
quickly: 20 trees and 300 draws after 300 warmup, on 400 observations.
The defaults are 50 trees and 500 draws after 500 warmup.

``` r

ctrl <- bartisan_control(num_trees = 20, num_burn = 300, num_draws = 300)

n <- 400
```

## The model

### A sum of trees

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

### The generalization

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

### The prior

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
sets how hard this shrinks: larger `k` pulls the fit toward the
intercept-only model. The default is a compromise that works across a
wide range of problems, and it is rarely worth changing.

**Nuisance parameters.** Each family sets its own, for example a
half-Cauchy on a scale parameter.

The prior does the work that cross-validation does in other tree
ensembles. There is no learning rate, no depth limit, and no early
stopping.

### Soft rules

A rule in classic BART is a step: an observation with \\x_j \le c\\ goes
left and everything else goes right. The fitted function is then a step
function, which is a poor description of a smooth truth.

Soft rules replace the step with a gate ([Linero and Yang
2018](#ref-linero2018)). An observation goes left with weight

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
#> 1 smoothstep     0.434     1.4
#> 2       hard     1.025     0.4
```

The true function has a standard deviation of about 4.9, so both are
fitting real structure and the soft fit is the more accurate.

### Sparsity

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

| Trees | noise, sparsity = FALSE | noise, sparsity = TRUE |
|------:|------------------------:|-----------------------:|
|    10 |                    0.28 |                   0.09 |
|    20 |                    0.50 |                   0.08 |
|    50 |                    0.95 |                   0.09 |
|   100 |                    1.00 |                   0.14 |

Mean prop_used for the five noise predictors. The five real predictors
sit at 1.00 in every cell. {.table}

Without the sparsity prior the advice to use fewer trees is essential:
at 50 trees the noise predictors are used in 95% of draws and are
indistinguishable from the real ones. With it, the noise predictors stay
near zero at every tree count, and reducing the trees buys little. The
recommendation is a workaround for the absence of the prior rather than
a property of variable selection in general.

Two things this does not fix. Sparsity needs signal: on a small sample
with a weak fit nothing separates, and the honest reading of a flat
table is that the data cannot rule anything out. And it says nothing
about effect size.
[`vignette("importance")`](https://ngreifer.github.io/bartisan/articles/importance.md)
covers both.

## How it is fitted

### Bayesian backfitting

The trees are updated one at a time. To update tree \\m\\, subtract the
contribution of every other tree from the predictor, leaving a partial
residual that only tree \\m\\ has to explain, then draw \\(T_m, M_m)\\
given that residual. Cycling over all trees, and then over the nuisance
parameters, is one sweep. The draws saved after warmup are the posterior
sample.

### Moving the trees

A tree is proposed by a small local change: growing a leaf into two,
pruning a pair of leaves back into one, or changing a splitting rule.
The proposal is accepted or rejected by a Metropolis-Hastings step.
Because growing and pruning change the number of parameters, this is a
reversible-jump sampler.

The acceptance ratio needs the likelihood of the tree with the leaf
parameters integrated out,

\\p(\text{residual} \mid T) = \int p(\text{residual} \mid T, M) \\ p(M)
\\ dM.\\

### Why the general case is hard, and what this package does

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

### Data augmentation

Some likelihoods are not Gaussian but become Gaussian once an extra
latent variable is imagined and drawn alongside everything else.

- A **probit** model is a Gaussian one on a latent scale that is
  truncated by the observed category ([Albert and Chib
  1993](#ref-albert1993)). Draw the latent value and the target is
  exactly quadratic.
- A **logistic** model becomes Gaussian conditional on a Pólya-Gamma
  variable ([Polson et al. 2013](#ref-polson2013)), which acts as a
  per-observation precision.
- A **censored** outcome becomes an ordinary one once the unobserved
  value is imputed from its truncated distribution.

`augment = TRUE`, the default where a scheme exists, turns these on.
They change the sampler and not the model: the posterior is the same,
reached faster.
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
says which families have one.

### Nuisance parameters

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

### The Dirichlet process mixture

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

## Practical matters

### Missing predictor values

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

The fit uses the incomplete variable rather than discarding it. This is
a prediction tool and not a substitute for thinking about why values are
missing; it assumes the pattern of missingness is itself informative,
which is sometimes true and sometimes not.

### Several chains

`chains` runs the sampler more than once from different starting points.
The default is one, which is enough for estimates. Running several is
what makes `fit$rhat` available, which is the only way to see whether
the sampler has converged.

``` r

fit_chains <- bartisan(y ~ . - eta, data = train, family = gaussian(),
                       chains = 4, control = ctrl)

fit_chains$rhat
#>                            quantity rhat ess_bulk ess_tail
#> 1                            loglik 1.39     8.87     52.0
#> 2                         aux.sigma 1.06    71.06    488.2
#> 3 eta.eta (worst over observations) 1.48     7.95     20.6
```

Read this table selectively. `aux.sigma` and `loglik` are close to 1,
which is what you want, and they respond to longer chains in the usual
way.

The `eta` row does not, and it is worth knowing why before it alarms
you. It is the worst value over all 400 observations, so it is a maximum
by construction, and forests mix slowly on their fitted values: two
chains can visit quite different sets of trees that imply very similar
predictions, which inflates a between-chain statistic without meaning
the answers disagree. Adding draws does not reliably bring it down.
Values above 1.05 on this row are ordinary for a BART fit, and it is
characteristic of the method rather than of this implementation.

What to do instead is check that the quantities you will report are
stable. If `sigma` and the log likelihood have converged and an
[`avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
estimate is the same across chains, the fit is usable whatever the `eta`
row says.

Chains run sequentially unless a `future` plan is set, in which case
they run in parallel.
[`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md)
covers this properly, including which statistic to look at when the
fitted values are the estimand.

### Choosing the settings

Most of
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
does not need to be touched. In rough order of how often changing one is
worth it:

| Setting | When to change it |
|----|----|
| `num_burn`, `num_draws` | when `rhat` says the chain has not converged |
| `num_trees` | more for a complex function and a large sample; fewer to speed up |
| `gate` | `"hard"` when the truth really is a step function, or for speed |
| `sparsity` | `FALSE` to recover classic BART, which is rarely what you want |
| `k` | to shrink harder toward the mean, in a very small sample |

Everything else exists so that the checks in the test suite can be
written, or because a specific family needs it.

The single decision that matters more than all of these is the family.
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
covers it.

### Where the time goes

Cost per sweep is roughly linear in the number of observations and the
number of trees. The family matters more than either: a family that
reaches the quadratic target is several times faster than one that does
not, which is why
[`lognormal_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
fits in about a second where
[`weibull_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
takes eight on the same data.

## Further reading

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
