# Bayesian causal forests

A varying-coefficient model set up for estimating a treatment effect: a
control function for the outcome under no treatment and a separate
forest for the effect, with the prior on the effect regularized more
heavily than the prior on the control function. This is Hahn, Murray and
Carvalho (2020) for a binary treatment and Woody, Carvalho, Hahn and
Murray (2020) for a continuous one.

## Usage

``` r
bcf(
  formula,
  treatment,
  data,
  family = NULL,
  moderators = NULL,
  propensity = TRUE,
  propensity_args = list(),
  ...
)
```

## Arguments

- formula:

  a model formula. The right-hand side lists the covariates; the
  treatment is named in `treatment` rather than here, and is removed
  from the covariates if it appears among them.

- treatment:

  a one-sided formula naming the treatment, as in `~ z`.

- data:

  a data frame.

- family:

  the outcome distribution, as in
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).

- moderators:

  a one-sided formula naming the covariates the treatment effect may
  vary with. The default, `NULL`, is all of them.

- propensity:

  what to do about the probability of treatment. `TRUE`, the default,
  fits a model for it and adds the fitted values to the **control
  function only**; `FALSE` fits nothing; a numeric vector or matrix is
  used as given; a one-sided formula fits it with that model rather than
  the outcome's covariates.

- propensity_args:

  a list of
  [`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
  settings for the propensity model, which is a prediction problem and
  so keeps the sparsity prior the outcome model turns off.

- ...:

  passed to
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md),
  including
  [`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
  settings.

## Value

A `bartisan_fit`, with the treatment's coefficient forest named for the
treatment. [`coef()`](https://rdrr.io/r/stats/coef.html) gives the
conditional effect for each observation and
[`marginaleffects::avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
the average.

## What the wrapper decides

Four things, all of which can be written out in
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
directly:

- The treatment gets a
  [`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) term, so
  the effect is a forest of its own with its own prior rather than
  whatever difference a single forest with the treatment among its
  predictors happens to produce.

- The propensity score goes in the control function and **not** in the
  effect forest. That is Hahn et al.'s recommendation and the flag
  `bcf`, `stochtree` and this package all provide; the point is to let
  the control function absorb the selection without letting the effect
  vary with it.

- The effect forest gets fewer trees than the control function, since
  patterns of effect heterogeneity are usually simpler than prognostic
  surfaces.

- A **binary** treatment has its coding drawn rather than fixed, which
  is `center = "estimate"` in
  [`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) and the
  parameter expansion of Hahn et al.'s section 5.3. At two levels it
  restricts nothing and costs nothing in recovery, and it removes the
  dependence on which level was written as 1. A treatment with more
  levels keeps the symmetric per-level coding, because there the drawn
  coding gives every contrast one shared shape;
  [`?vc`](https://ngreifer.github.io/bartisan/reference/vc.md) has the
  numbers on both.

- `sparsity = FALSE` on the outcome model. A variable-selection prior on
  the variable whose contrast is the estimand puts a point mass at
  exactly zero in the posterior of the effect; see the measurements in
  [`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md).
  The propensity model keeps the default, because predicting who was
  treated is a prediction problem.

## The treatment's type

The treatment decides the model for the propensity score and what that
score even is.

|  |  |  |
|----|----|----|
| treatment | propensity score | model |
| binary | one column, the probability of treatment | [`binomial()`](https://rdrr.io/r/stats/family.html) |
| `K` categories | the whole vector of assignment probabilities | [`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) |
| continuous | a conditional density, not a regression | not fitted |

For a treatment with more than two categories the balancing score is the
*vector* of assignment probabilities (Imbens, 2000; Imai and van Dyk,
2004), not any one of them, so all of them go into the control function.
They sum to one and are therefore collinear, which costs a tree ensemble
nothing.

For a continuous treatment the analogue is the conditional density of
the treatment given the covariates evaluated at the observed dose
(Hirano and Imbens, 2004), which needs a density model rather than a
regression, so `propensity = TRUE` is refused and a score supplied as a
number is used as given.

## The assumption a continuous treatment carries

With a continuous treatment this fits `f0(x) + z * f1(x)`: a dose
response that is **linear in the dose**, with a slope that varies. For a
binary treatment that is no assumption at all. For a continuous one it
is a real one, and it is the assumption Woody et al. (2020) make and
diagnose. If the dose response itself might be curved, either put the
treatment in as an ordinary predictor or give `f1` the treatment among
its moderators, which makes the effect vary across the dose; see
[`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md).

## Predicting for new data

The propensity score is a predictor of the control function, and it is
one the caller never named, so `newdata` taken from their own frame does
not carry it. [`predict()`](https://rdrr.io/r/stats/predict.html)
rebuilds it from the model this kept, which is what makes that work at
all.

One consequence is worth knowing. A predictor goes through the quantile
transform, which is a step function, so a score that rebuilds to within
1e-11 can still land on the other side of a step: predictions for the
*training* data come back to within a few percent of the response's
spread rather than exactly. Passing the stored score in `newdata` – it
is in `fit$bcf$propensity` – removes the reconstruction and reproduces
the fit to machine precision. Supplying `propensity` as a number rather
than fitting it has the same effect, and then `newdata` must carry the
column.

## References

Hahn, P. R., Murray, J. S., & Carvalho, C. M. (2020). Bayesian
regression tree models for causal inference: regularization,
confounding, and heterogeneous effects. *Bayesian Analysis*, 15(3),
965–1056. [doi:10.1214/19-BA1195](https://doi.org/10.1214/19-BA1195)

Woody, S., Carvalho, C. M., Hahn, P. R., & Murray, J. S. (2020).
Estimating heterogeneous effects of continuous exposures using Bayesian
tree ensembles.
[doi:10.48550/arXiv.2007.09845](https://doi.org/10.48550/arXiv.2007.09845)

## See also

[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
and [`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) for
the general interface this is written in terms of, and
[`vignette("causal")`](https://ngreifer.github.io/bartisan/articles/causal.md).

## Examples

``` r
set.seed(1)
n <- 300
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
d$z <- rbinom(n, 1, plogis(d$x1))
d$y <- d$x1 + d$z * (1 + d$x2) + rnorm(n)

fit <- bcf(y ~ x1 + x2, treatment = ~ z, data = d, family = gaussian(),
           num_trees = c(10, 5), num_burn = 50, num_draws = 50,
           verbose = FALSE)

head(coef(fit))
#>                z
#> [1,]  2.02330640
#> [2,] -0.03702258
#> [3,]  2.28847658
#> [4,]  0.97773484
#> [5,]  2.19884085
#> [6,]  2.19558258
```
