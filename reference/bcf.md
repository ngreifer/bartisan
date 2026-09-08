# Bayesian causal forests

Fits a varying-coefficient model set up for estimating a treatment
effect: a control function for the outcome under no treatment and a
separate forest for the effect, with the prior on the effect regularized
more heavily than the prior on the control function. This is Hahn,
Murray and Carvalho (2020) for a binary treatment and Woody, Carvalho,
Hahn and Murray (2020) for a continuous one.

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
  from the covariates if it appears among them, so `y ~ .` is usually
  the right specification.

- treatment:

  a one-sided formula naming the treatment, as in `~ z`. The treatment
  may be binary, categorical, or continuous, and which it is decides
  what the propensity score is and how it is modeled; see Details.

- data:

  a data frame containing the variables named in `formula`.

- family:

  the response distribution, given as a
  [stats::family](https://rdrr.io/r/stats/family.html) object, as one of
  the families in
  [bartisan-families](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
  or as the name of either. A `family` object is accepted when the
  distribution it names is one this package implements, since the
  likelihood is the package's rather than the object's: a `family`
  object carries a link and a variance function and not a density, so
  one naming anything else (e.g.,
  [stats::inverse.gaussian](https://rdrr.io/r/stats/family.html), or a
  Tweedie from another package) is an error rather than something a
  likelihood can be built from, and
  [`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  is the route for those. Links are used as supplied, and a link the
  package does not compile is composed onto the scale its family works
  on. Default is `NULL`, in which case the family is read off the
  response and a message reports the choice; see Details for the rules.

- moderators:

  a one-sided formula naming the covariates the treatment effect may
  vary with. Default is `NULL` to let the effect vary with every
  covariate.

- propensity:

  what to do about the probability of treatment, given as either a
  logical value, a numeric vector or matrix, or a one-sided formula.
  Default is `TRUE` to fit a model for it and add the fitted values to
  the **control function only**. `FALSE` fits nothing; a numeric vector
  or matrix is used as given; and a one-sided formula fits it with the
  predictors that formula names rather than with the outcome's
  covariates. Note that a continuous treatment has no propensity score
  that is a probability, so `TRUE` is refused for one; see Details.

- propensity_args:

  a list of
  [`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
  settings for the propensity model. Default is
  [`list()`](https://rdrr.io/r/base/list.html) to leave every setting at
  its own default.

- ...:

  passed to
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md),
  including
  [`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
  settings.

## Value

A `<bartisan_fit>` object, as
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
returns, with the treatment's coefficient forest named for the
treatment. [`coef()`](https://rdrr.io/r/stats/coef.html) gives the
conditional effect for each observation and
[`marginaleffects::avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
the average; see
[bartisan-marginaleffects](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md).

## Details

### What the Wrapper Decides

Five things, all of which can be written out in
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
directly.

**The treatment gets a
[`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) term**, so
the effect is a forest of its own with its own prior rather than
whatever difference a single forest with the treatment among its
predictors happens to produce. **The propensity score goes in the
control function and not in the effect forest**, which is Hahn et al.'s
recommendation and the flag bcf, stochtree and this package all provide;
the point is to let the control function absorb the selection without
letting the effect vary with it. **The effect forest gets fewer trees
than the control function**, since patterns of effect heterogeneity are
usually simpler than prognostic surfaces.

**A binary treatment has its coding drawn rather than fixed**, which is
`center = "estimate"` in
[`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) and the
parameter expansion of Hahn et al.'s section 5.3. At two levels it
restricts nothing and costs nothing in recovery, and it removes the
dependence on which level was written as 1. A treatment with more levels
keeps the symmetric per-level coding, because there the drawn coding
gives every contrast one shared shape;
[`?vc`](https://ngreifer.github.io/bartisan/reference/vc.md) has the
numbers on both.

**The sparsity prior is left at its default**, which is on. Turning it
off is what a contrast on a predictor calls for, because a
variable-selection prior can drop that predictor from the forest
entirely and put a point mass at exactly zero in the posterior of the
effect; see the measurements in
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md).
That is a reason to reach for `sparsity = FALSE` when the treatment is
one predictor among many in a single forest, which is how a contrast is
written without
[`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md).

It is not the situation here. The treatment is not split on at all: it
is the coefficient, carried by a forest of its own, so no splitting
proportion can drop it and there is no mass to pile at zero. What the
prior selects among on that forest is the moderators, and dropping all
of them leaves an effect that does not vary rather than one that is
zero, which is the ordinary shrinkage a heterogeneity model wants.
Measured on a truth with one moderator among twenty covariates, the
point mass is absent at every setting, the average effect is recovered
either way, and the conditional effect is recovered better with the
prior on than off.

Either forest can still be set on its own, `sparsity = c(FALSE, TRUE)`
leaving the control function every predictor and asking only the effect
forest to select.

### The Treatment's Type

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

### The Assumption a Continuous Treatment Carries

With a continuous treatment this fits `f0(x) + z * f1(x)`: a dose
response that is **linear in the dose**, with a slope that varies. For a
binary treatment that is no assumption at all. For a continuous one it
is a real one, and it is the assumption Woody et al. (2020) make and
diagnose. If the dose response itself might be curved, either put the
treatment in as an ordinary predictor or give `f1` the treatment among
its moderators (i.e., write the term as `vc(z, ~ z + ...)` in
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
directly), which makes the effect vary across the dose; see
[`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md).

### Predicting for New Data

The propensity score is a predictor of the control function, and it is
one the caller never named, so `newdata` taken from their own frame does
not carry it. [`predict()`](https://rdrr.io/r/stats/predict.html)
rebuilds it from the model this kept, which is what makes that work at
all.

One consequence is worth knowing. A predictor goes through the quantile
transform, which is a step function, so a score that rebuilds to within
1e-11 can still land on the other side of a step: predictions for the
*training* data come back to within a few percent of the response's
spread rather than exactly. Passing the stored score in `newdata` (it is
in `fit$bcf$propensity`) removes the reconstruction and reproduces the
fit to machine precision. Supplying `propensity` as a number rather than
fitting it has the same effect, and then `newdata` must carry the
column.

## References

Hahn, P. R., Murray, J. S., & Carvalho, C. M. (2020). Bayesian
regression tree models for causal inference: regularization,
confounding, and heterogeneous effects. *Bayesian Analysis*, 15(3),
965–1056. [doi:10.1214/19-BA1195](https://doi.org/10.1214/19-BA1195)

Imai, K., & van Dyk, D. A. (2004). Causal inference with general
treatment regimes: generalizing the propensity score. *Journal of the
American Statistical Association*, 99(467), 854–866.

Imbens, G. W. (2000). The role of the propensity score in estimating
dose-response functions. *Biometrika*, 87(3), 706–710.

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
data("rhc")

set.seed(123)

# The effect of right heart catheterization on death, free to vary with
# every covariate, with the propensity score entering the control function
# alone
fit <- bcf(death ~ . - days, treatment = ~ rhc, data = rhc,
           family = binomial(), num_trees = c(10, 5), num_burn = 50,
           num_draws = 50,
           propensity_args = list(num_trees = 10, num_burn = 50,
                                  num_draws = 50))

# One conditional effect per patient, which is what the effect forest comes
# to at each observation
head(coef(fit))
#>            rhc
#> [1,] 0.2976305
#> [2,] 0.1495210
#> [3,] 0.2202776
#> [4,] 0.1666915
#> [5,] 0.5088872
#> [6,] 0.4356817

# The average effect over the sample
if (rlang::is_installed("marginaleffects")) {
  marginaleffects::avg_comparisons(fit, variables = "rhc")
}
#> 
#>  Estimate  2.5 % 97.5 %
#>    0.0532 0.0048  0.119
#> 
#> Term: rhc
#> Type: response
#> Comparison: 1 - 0
#> 
```
