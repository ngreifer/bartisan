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
  treat,
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
  treatment is named in `treat` rather than here, and is removed from
  the covariates if it appears among them.

- treat:

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

A `<bcf_fit>` object, which is the `<bartisan_fit>`
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
returns with a class in front of it and the treatment's coefficient
forest named for the treatment. Everything that works on a
`<bartisan_fit>` works here unchanged; the class exists so that methods
needing a named treatment have something to dispatch on.
[`coef()`](https://rdrr.io/r/stats/coef.html) gives the conditional
effect for each observation and
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

[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
for the average or conditional effect;
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
fit <- bcf(death ~ . - days, treat = ~ rhc, data = rhc,
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

# The average effect over the sample, on the response scale, which for a
# binomial fit makes it a risk difference rather than a log odds ratio
estimate_effect(fit)
#> Average treatment effect (difference)
#> 
#> Treatment: "rhc"
#> Averaged over 1500 units
#> 
#>     contrast estimate  lower upper    n
#>  Y[1] - Y[0]   0.0539 0.0048 0.119 1500
#> 
#> Average potential outcomes
#> 
#>  quantity estimate lower upper
#>      Y[0]    0.632 0.607 0.659
#>      Y[1]    0.686 0.652 0.720
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ Y[a] is the average response with "rhc" set to a.

# Or the whole picture at once
summary(fit)
#> Generalized BART
#> 
#> Call:
#> bcf(formula = death ~ . - days, treat = ~rhc, data = rhc, family = binomial(), 
#>     propensity_args = list(num_trees = 10, num_burn = 50, num_draws = 50), 
#>     num_trees = c(10, 5), num_burn = 50, num_draws = 50)
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 1500
#> Structure: 2 forests of 10 and 5 trees, soft decision rules
#> Draws: 50
#> 
#> Nuisance parameters
#>          mean    sd  lower upper
#> b.rhc.0 0.431 0.424 -0.595 0.965
#> b.rhc.1 1.148 0.327  0.732 1.690
#> 
#> Predictor usage
#> Splitting rules per draw, and how often used at all.
#> 
#> Predictor "(Intercept)":
#>             mean    sd lower  upper prop_used
#> age         5.46 2.305     2 10.550      1.00
#> surv2m      5.78 3.079     2 11.775      1.00
#> paco2       1.64 1.367     0  4.775      0.78
#> card        1.44 1.280     0  3.775      0.64
#> race        0.64 0.776     0  2.775      0.50
#> hema        0.34 0.593     0  1.000      0.30
#> resp        0.26 0.443     0  1.000      0.26
#> edu         0.30 0.735     0  2.775      0.18
#> aps         0.14 0.405     0  1.000      0.12
#> sex         0.08 0.274     0  1.000      0.08
#> pafi        0.12 0.435     0  1.775      0.08
#> meanbp      0.04 0.198     0  0.775      0.04
#> .propensity 0.04 0.198     0  0.775      0.04
#> crea        0.02 0.141     0  0.000      0.02
#> 
#> Predictor "rhc":
#>             mean    sd lower upper prop_used
#> pafi        1.52 0.789     0     3      0.90
#> crea        1.26 0.828     0     3      0.84
#> hema        1.42 1.357     0     4      0.66
#> age         0.70 0.707     0     2      0.56
#> aps         0.58 0.642     0     2      0.50
#> edu         0.70 0.814     0     2      0.48
#> sex         0.34 0.688     0     2      0.24
#> meanbp      0.26 0.487     0     1      0.24
#> race        0.24 0.476     0     1      0.22
#> paco2       0.22 0.418     0     1      0.22
#> card        0.18 0.388     0     1      0.18
#> surv2m      0.16 0.370     0     1      0.16
#> resp        0.14 0.351     0     1      0.14
#> .propensity 0.00 0.000     0     0      0.00
#> 
#> ℹ This fit has a treatment, "rhc". `estimate_effect()` reports its effect, with
#>   the average potential outcomes beside it.
```
