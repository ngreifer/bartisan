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

- moderators:

  a one-sided formula naming the covariates the treatment effect may
  vary with. Default is `NULL` to let the effect vary with every
  covariate.

- propensity:

  what to do about the propensity score, given as either a logical
  value, a numeric vector or matrix, or a one-sided formula. Default is
  `TRUE` to fit a model for it and add the fitted values to the control
  function. `FALSE` fits nothing, a numeric vector or matrix is used as
  given, and a one-sided formula fits it with the predictors that
  formula names. For a continuous treatment the fitted values are the
  treatment's conditional mean given the covariates; see Details.

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
[`bartisan-marginaleffects`](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md).

## Details

### The Model

`bcf()` writes several settings into a
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
call, all of which can be written out there directly. The treatment gets
a [`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) term, so
the effect is a forest of its own with its own prior rather than
whatever difference a single forest with the treatment among its
predictors happens to produce. The propensity score goes in the control
function and not in the effect forest, so that the control function
absorbs the selection while the effect stays free of it. The effect
forest gets fewer trees than the control function, since patterns of
effect heterogeneity are usually simpler than prognostic surfaces. A
binary treatment has its coding drawn rather than fixed, which is
`center = "estimate"` in
[`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) and makes
the answer the same whichever level was written as 1; a treatment with
more levels keeps the symmetric per-level coding. And `sparsity` is left
at its default, which is on.

The propensity score that enters the control function is the posterior
mean of a separate model for the treatment, fit before the outcome model
and then held fixed. Its uncertainty is therefore not carried into the
interval for the effect, and the outcome has no say in the score, as it
would in a joint model of the two. See
[`vignette("causal")`](https://ngreifer.github.io/bartisan/articles/causal.md)
for more on this choice. Including the score at all matters because,
with flexible priors on the outcome that are independent of the model
for the treatment, the implied prior on the amount of confounding bias
concentrates near zero (Linero, 2024).

### Setting `sparsity`

A variable-selection prior can drop a predictor from the forest entirely
and put a point mass at exactly zero in the posterior of a contrast on
it, which is why `sparsity = FALSE` is the setting to reach for when a
treatment is one predictor among many in a single forest.

Here the treatment is the coefficient rather than a predictor the forest
splits on, so no splitting proportion can drop it and there is no mass
to pile at zero. What the prior selects among on the effect forest is
the moderators, and dropping all of them leaves an effect that does not
vary rather than one that is zero, which is the shrinkage a
heterogeneity model wants. Either forest can still be set on its own,
`sparsity = c(FALSE, TRUE)` leaving the control function every predictor
and asking only the effect forest to select.

### Setting `treat`

The treatment decides the model for the propensity score and what that
score is.

|  |  |  |
|----|----|----|
| treatment | propensity score | model |
| binary | one column, the probability of treatment | [`binomial()`](https://rdrr.io/r/stats/family.html) |
| `K` categories | the whole vector of assignment probabilities | [`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) |
| continuous | the conditional mean of the treatment | [`gaussian()`](https://rdrr.io/r/stats/family.html) |

For a treatment with more than two categories the balancing score is the
whole vector of assignment probabilities rather than any one of them, so
all of them go into the control function. They sum to one and are
therefore collinear, which costs a tree ensemble nothing.

### Continuous Treatments

For a continuous treatment, the score added to the control function is
the treatment's conditional mean given the covariates, fit with a
Gaussian model. The balancing score proper would be the conditional
density of the treatment at the observed dose, but the purpose of the
score here is to let the control function absorb the confounding, and
Linero (2024) shows for linear models that it is the conditional mean of
the treatment whose absence leaves the prior on the confounding bias
concentrated near zero. A score supplied as a number is used as given.

A continuous treatment also carries an assumption. This fits
`f0(x) + z * f1(x)`, a dose response that is linear in the dose with a
slope that varies with the covariates, which for a binary treatment is
no assumption at all and for a continuous one is a real one. Where the
dose response itself might be curved, either put the treatment in as an
ordinary predictor or give `f1` the treatment among its moderators by
writing the term as `vc(z, ~ z + ...)` in
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
directly, which lets the effect vary across the dose; see
[`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md).

### Predicting for New Data

The propensity score is a predictor of the control function and one the
caller never named, so a `newdata` data frame taken from their own frame
does not carry it. [`predict()`](https://rdrr.io/r/stats/predict.html)
rebuilds it from the model the fit kept, which makes that work.
Supplying `propensity` as a number rather than fitting it removes the
reconstruction, and the data given as `newdata` then has to carry the
column itself.

## References

Hahn, P. R., Murray, J. S., & Carvalho, C. M. (2020). Bayesian
regression tree models for causal inference: regularization,
confounding, and heterogeneous effects. *Bayesian Analysis*, 15(3),
965–1056. [doi:10.1214/19-BA1195](https://doi.org/10.1214/19-BA1195)

Linero, A. R. (2024). In nonparametric and high-dimensional models,
Bayesian ignorability is an informative prior. *Journal of the American
Statistical Association*, 119(548), 2785–2798.
[doi:10.1080/01621459.2023.2278202](https://doi.org/10.1080/01621459.2023.2278202)

Woody, S., Carvalho, C. M., Hahn, P. R., & Murray, J. S. (2020).
Estimating heterogeneous effects of continuous exposures using Bayesian
tree ensembles.
[doi:10.48550/arXiv.2007.09845](https://doi.org/10.48550/arXiv.2007.09845)

## See also

[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
for the average or conditional effect;
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
and [`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) for
the general interface this is written in terms of;
[`vignette("varying")`](https://ngreifer.github.io/bartisan/articles/varying.md)
for the varying-coefficient model this is a case of and what it writes
into the call; and
[`vignette("causal")`](https://ngreifer.github.io/bartisan/articles/causal.md)
for the assumptions under which the effect is causal.

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

# One conditional effect per patient: the effect forest evaluated at each
# observation
head(coef(fit))
#>               rhc
#> [1,]  0.276254026
#> [2,]  0.220726838
#> [3,] -0.057821159
#> [4,]  0.025513898
#> [5,] -0.005822382
#> [6,]  0.071148221

# The average effect over the sample, on the response scale, which for a
# binomial fit makes it a risk difference rather than a log odds ratio
estimate_effect(fit)
#> Average treatment effect (difference)
#> 
#> Treatment: `rhc`
#> Averaged over 1500 units
#> 
#>     contrast estimate    lower  upper    n
#>  Y[1] - Y[0]   0.0289 -0.00395 0.0687 1500
#> 
#> Average potential outcomes
#> 
#>  quantity estimate lower upper
#>      Y[0]    0.643 0.624 0.667
#>      Y[1]    0.672 0.645 0.699
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ Y[a] is the average response with `rhc` set to "a".

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
#>          mean    sd lower upper
#> b.rhc.0 0.834 0.172 0.503 1.084
#> b.rhc.1 1.165 0.288 0.745 1.610
#> 
#> Predictor usage
#> Splitting rules per draw, and how often used at all.
#> 
#> Predictor "(Intercept)":
#>             mean    sd lower  upper prop_used
#> aps         2.02 1.672     1  6.000      1.00
#> paco2       2.42 0.702     2  4.000      1.00
#> surv2m      7.14 1.830     4 11.000      1.00
#> pafi        1.34 0.823     0  2.775      0.82
#> card        0.70 0.614     0  2.000      0.62
#> crea        0.80 0.990     0  3.000      0.46
#> meanbp      0.48 0.646     0  2.000      0.40
#> hema        0.50 0.814     0  2.775      0.34
#> age         0.34 0.519     0  1.000      0.32
#> edu         0.10 0.303     0  1.000      0.10
#> race        0.04 0.198     0  0.775      0.04
#> resp        0.04 0.198     0  0.775      0.04
#> sex         0.00 0.000     0  0.000      0.00
#> .propensity 0.00 0.000     0  0.000      0.00
#> 
#> Predictor "rhc":
#>             mean    sd lower upper prop_used
#> age         6.02 2.005 2.225 9.775      1.00
#> surv2m      1.46 1.705 0.000 5.775      0.54
#> resp        0.22 0.465 0.000 1.000      0.20
#> crea        0.10 0.303 0.000 1.000      0.10
#> sex         0.08 0.274 0.000 1.000      0.08
#> hema        0.10 0.463 0.000 1.000      0.06
#> edu         0.06 0.314 0.000 0.775      0.04
#> aps         0.04 0.198 0.000 0.775      0.04
#> pafi        0.04 0.198 0.000 0.775      0.04
#> paco2       0.04 0.198 0.000 0.775      0.04
#> race        0.02 0.141 0.000 0.000      0.02
#> meanbp      0.02 0.141 0.000 0.000      0.02
#> card        0.02 0.141 0.000 0.000      0.02
#> .propensity 0.00 0.000 0.000 0.000      0.00
#> 
#> ℹ This fit has a treatment, "rhc". `estimate_effect()` reports its effect, with
#>   the average potential outcomes beside it.
```
