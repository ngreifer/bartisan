# Effects, Curves, and Interactions

## Introduction

A forest has no coefficients: there is no table of slopes to read, and
no standard error to put beside one. This is the part of the workflow
that changes most when moving from
[`glm()`](https://rdrr.io/r/stats/glm.html) to BART, and it is the part
where the change is an improvement rather than a cost.

The replacement is to ask the fitted model questions about predictions.
What does it predict for these people? What would it predict if this
variable were different? How much does the answer differ between groups?
The *marginaleffects* package ([Arel-Bundock et al.
2024](#ref-arelbundock2024)) asks all of them, and returns posterior
intervals with the answers.

Two of those questions the package answers itself, and it is worth
knowing which before reading further.
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
contrasts one binary or factor treatment and averages it as the ATE, the
ATT, the ATC, or one effect per unit, optionally within subgroups;
[`partial_dependence()`](https://ngreifer.github.io/bartisan/reference/partial_dependence.md)
averages the fitted surface over the sample at each value of one or two
predictors. Neither needs a suggested package, and
[`vignette("causal")`](https://ngreifer.github.io/bartisan/articles/causal.md)
is the worked version of the first. Everything else in this vignette is
*marginaleffects*: slopes, contrasts between arbitrary covariate values,
a continuous treatment, hypotheses comparing one estimate to another,
and any grid that fixes the other covariates at chosen values rather
than averaging over them. That last class is the one to watch for,
because it is where the native route stops rather than where it is
merely less convenient.

This vignette covers the questions worth asking and how to phrase them.
[`vignette("bartisan")`](https://ngreifer.github.io/bartisan/articles/bartisan.md)
is the shorter tour, and this expands its section on interpreting the
fit.

In this guide, we will start from the three kinds of question the
package answers and then work through them in turn. First we’ll take
average effects and the choice of step for a numeric predictor, along
with the one prior setting (`sparsity`) that can quietly attenuate a
contrast. Next we’ll split those effects by subgroup and test a
moderation with the difference of the two, then plot the shape of a
fitted relationship and choose the scale on which an effect is reported.
Finally we’ll write the effect as a parameter rather than a contrast
with [`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md), and
close with what none of these estimates can be taken to mean.

``` r

library(bartisan)
library(marginaleffects)

data("rhc")

model <- death ~ rhc + age + sex + race + edu + aps + meanbp + resp +
  hema + pafi + paco2 + crea + surv2m + card

set.seed(2026)

fit <- bartisan(model, data = rhc, family = binomial(), chains = 4)
```

## Predictions, Comparisons, and Slopes

Everything below is one of three things.

A **prediction** is what the model expects for a set of covariate
values.
[`predictions()`](https://rdrr.io/pkg/marginaleffects/man/predictions.html)
gives one per row,
[`avg_predictions()`](https://rdrr.io/pkg/marginaleffects/man/predictions.html)
averages them.

A **comparison** is the difference between two predictions that differ
in one variable. This is the closest thing to a regression coefficient,
and it is usually what we want.

A **slope** is the derivative of the prediction with respect to a
numeric variable. It is the least useful of the three here, for a reason
given below.

## Average Effects (`avg_comparisons()`)

[`avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
with no `variables` argument gives every predictor at once, which for
this model is a long table. A few at a time is easier to read:

``` r

avg_comparisons(fit, variables = c("rhc", "age", "card"))
#> 
#>  Term Contrast Estimate    2.5 % 97.5 %
#>  age  +1        0.00308  0.00131 0.0048
#>  card yes - no  0.02794 -0.00266 0.0852
#>  rhc  1 - 0     0.05588  0.00000 0.1052
#> 
#> Type: response
```

This reads as a coefficient table. Each estimate is an average
difference in predicted probability, holding everything else at each
patient’s own values: for a factor, between the levels named in the
`Contrast` column, and for a numeric predictor, for an increase of one
unit. `rhc` is coded `0` and `1`, so its one-unit contrast is the
treatment effect.

### The Splitting Prior and a Contrast (`sparsity`)

An estimate here can come back as exactly zero, and an interval bound
with it; that is not a rounding artifact. The default splitting prior is
a variable-selection prior (i.e., one that can leave a predictor out of
the forest altogether), so in a draw where it uses the predictor in no
tree the prediction does not depend on it and the contrast is exactly
zero; the posterior of the contrast is a mixture with a point mass
there, holding whatever share of draws dropped the predictor.
*marginaleffects* centers a posterior at its median, so once that point
mass holds half of it the reported estimate is exactly zero however far
the rest of the posterior sits from zero. Setting
`options(marginaleffects_posterior_center = mean)` asks for the mean
instead, which is the summary
[`predict()`](https://rdrr.io/r/stats/predict.html) reports.

It matters more than it sounds: on a weak signal the prior attenuates
the estimate substantially and its interval covers the truth well below
its nominal rate. A strong effect is untouched, because the prior never
has reason to drop a predictor that is earning its splits, so this is a
weak-signal problem rather than a general one.

If a contrast is what we are reporting, we fit with `sparsity = FALSE`,
or with `split_prior`, which fixes the weights (i.e., the probability
that each predictor is chosen for a split) and so cannot drop anything.
[`?bartisan_control`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
explains both halves of this.

### Choosing the Step (`variables`)

One unit is the default and is often the wrong scale. One point of an
illness score is a small change; ten points is a difference someone
would notice.

``` r

avg_comparisons(fit, variables = list(aps = 10))
#> 
#>  Estimate 2.5 % 97.5 %
#>    0.0122     0 0.0291
#> 
#> Term: aps
#> Type: response
#> Comparison: +10
```

A numeric effect should always be reported together with the step it was
computed at. Unlike a linear model, the answer here is not ten times the
one-unit effect, because the relationship is not assumed to be a
straight line.

## Effects for Subgroups (`by`)

`by` splits the average by a grouping variable.

``` r

avg_comparisons(fit, variables = "rhc", by = "card")
#> 
#>  card Estimate 2.5 % 97.5 %
#>   no    0.0573     0  0.110
#>   yes   0.0525     0  0.104
#> 
#> Term: rhc
#> Type: response
#> Comparison: 1 - 0
```

The two subgroup estimates are close, and both intervals reach zero.
`estimate_effect(fit, treat = "rhc", by = ~ card)` asks the same
question natively and returns the same intervals; its point estimates
differ a little, because it centers each posterior at its mean where
*marginaleffects* centers at its median.

A common mistake is to stop here and conclude that the effect differs
between groups; that comparison is not a test. The question is whether
the two effects differ from each other, which needs the difference of
the two (i.e., a difference of differences) with an interval of its own.

``` r

avg_comparisons(fit, variables = "rhc", by = "card",
                hypothesis = ~pairwise)
#> 
#>    Hypothesis Estimate   2.5 % 97.5 %
#>  (yes) - (no) -0.00132 -0.0494 0.0122
#> 
#> Type: response
```

The difference is small with an interval covering zero, so there is no
evidence here that the effect of catheterization depends on
cardiovascular disease. This is how to test an interaction in a model
that never had an interaction term to test, and the interval on that
difference is the only thing that separates a real interaction from two
subgroup estimates that merely look different.

It is also the clearest place where the native route stops:
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
reports each subgroup effect with its own interval and has no way to ask
for the difference between them, so this call is the reason to reach for
*marginaleffects* even when the subgroup estimates came from the other
route.

[`avg_predictions()`](https://rdrr.io/pkg/marginaleffects/man/predictions.html)
does the same thing for predictions rather than differences, which is
useful for describing groups:

``` r

avg_predictions(fit, by = "card")
#> 
#>  card Estimate 2.5 % 97.5 %
#>   no     0.634 0.607  0.660
#>   yes    0.691 0.656  0.731
#> 
#> Type: response
```

## The Shape of a Relationship (`plot_predictions()`)

Averages hide shape. The fitted function is seen by plotting predictions
against one predictor while the others are either averaged over or held
fixed.
[`partial_dependence()`](https://ngreifer.github.io/bartisan/reference/partial_dependence.md)
averages them over the sample and draws the result natively, and
`plot(fit, ~ x)` is the short way to it; the section uses
`plot_predictions(draw = FALSE)` instead because it gives more control
over the grid and over what is held fixed, which is what the hand-built
*ggplot2* calls below are for:

``` r

library(ggplot2)

curve <- plot_predictions(fit, condition = "aps", draw = FALSE)

ggplot(curve, aes(aps, estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line() +
  labs(x = "APACHE III score on day 1", y = "fitted probability of death") +
  theme_bw(base_size = 9)
```

![](effects_files/figure-html/pdp-1.png)

The probability of death rises with the illness score, and the rise is
not a straight line on any scale the model was told about; nothing was
specified to find the shape.

The band is a credible interval and is wide at the top, where few
patients were that sick, so the ends of a curve deserve more caution
than the middle: with little data out there the forest shrinks its
predictions toward the overall mean, which flattens the curve at both
edges of the predictor’s range.

Adding a second variable shows how the shape differs across groups.

``` r

curve2 <- plot_predictions(fit, condition = c("aps", "rhc"), draw = FALSE)

ggplot(curve2, aes(aps, estimate, color = factor(rhc))) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = factor(rhc)),
              alpha = 0.15, color = NA) +
  geom_line() +
  labs(x = "APACHE III score on day 1", y = "fitted probability of death",
       color = "catheterized", fill = "catheterized") +
  theme_bw(base_size = 9)
```

![](effects_files/figure-html/pdp2-1.png)

The two curves run close together; if they diverged, that would be a
moderation worth reporting, and the difference of differences above is
how to put a number on it.

## Scales (`type`)

`type` chooses the scale on which predictions are made and therefore the
scale on which effects are reported.

``` r

avg_comparisons(fit, variables = "rhc", type = "link")
#> 
#>  Estimate 2.5 % 97.5 %
#>     0.296     0  0.571
#> 
#> Term: rhc
#> Type: link
#> Comparison: 1 - 0
```

On the link scale this is a difference in log-odds, which is what a
logistic regression coefficient is. It is the less useful of the two
here: a difference in probability is interpretable without reference to
the model and is the number a reader can act on, whereas a difference in
log-odds needs a baseline before it means anything.

For survival families, `type = "survival"` with a `times` argument gives
a difference in survival probability at a horizon. See
[`vignette("survival")`](https://ngreifer.github.io/bartisan/articles/survival.md).

## Slopes and the Predictor Transform (`avg_slopes()`)

[`avg_slopes()`](https://rdrr.io/pkg/marginaleffects/man/slopes.html)
reports a derivative, and whether that derivative means anything depends
on `x_transform` in
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md).
The fit is a smooth function of the *transformed* predictor, so the
question is what the transform does to the original one.

The default, `"smoothcdf"`, is differentiable, so
[`avg_slopes()`](https://rdrr.io/pkg/marginaleffects/man/slopes.html)
works. The way to check such a thing is to shrink the numerical step and
see whether the answer settles. On a surface whose average slope over
the sample is .5485, the default gives .553, .559, .560, .561, .561 as
the step runs from 1e-1 down to 1e-5. It settles.

Under `"quantile"` the same sequence is .577, .882, 3.73, 31.8, 311.
That transform maps each predictor through its empirical distribution
function, which is a step, so the fitted function is a step function of
the original predictor and the difference quotient grows without bound
as the step shrinks. There is no derivative there to estimate, and the
number
[`avg_slopes()`](https://rdrr.io/pkg/marginaleffects/man/slopes.html)
returns is a property of the step size rather than of the fit.

Under `"range"` the sequence is .544, .543, .543, .543, .543, which
settles closest to the truth. Writing the fit as \\f(T(x))\\, a slope is
\\f'(T(x))\\T'(x)\\: an affine \\T\\ has a known constant derivative, so
only \\f'\\ is estimated, where a smoothed distribution function
contributes an estimated density and the two errors multiply. So **we
refit with `x_transform = "range"` when a slope is the quantity being
reported**, and leave the default alone when it is not.

For everything else we use
[`comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
with a step we can interpret, as with `aps = 10` above, which evaluates
the fit at two points a substantive distance apart rather than dividing
by a vanishing one. With `gate = "hard"` the fit is piecewise constant
under any transform and has no derivative worth taking. This is
documented at `?bartisan-marginaleffects`.

## Varying Coefficients (`vc()`)

Everything above reads an effect out of a fitted surface by asking the
model what it predicts under two versions of the data. There is another
way to write the model, in which the effect is a parameter rather than a
contrast:

\\f_0(x) + z\\f_1(x)\\

Here \\f_1\\ is a forest of its own, and it *is* the effect of `z`: how
much the prediction moves per unit of `z`, as a function of the other
predictors.
[`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) asks for
it. This is the varying-coefficient model of Deshpande et al.
([2026](#ref-deshpande2026)), of which Hahn et al.
([2020](#ref-hahn2020)) is the case of one binary covariate and Woody et
al. ([2020](#ref-woody2020)) the case of one continuous one.

``` r

fit_vc <- bartisan(death ~ age + sex + race + edu + aps + meanbp + resp +
                     hema + pafi + paco2 + crea + surv2m + card + vc(rhc),
                   data = rhc, family = binomial(), chains = 4,
                   sparsity = FALSE)

head(coef(fit_vc))
#>         rhc
#> [1,] 0.4892
#> [2,] 0.1855
#> [3,] 0.2007
#> [4,] 0.3322
#> [5,] 0.2895
#> [6,] 0.2646
```

[`coef()`](https://rdrr.io/r/stats/coef.html) returns one value per
patient, which is what a coefficient becomes when it is allowed to vary.
It is on the link scale, so for this binary outcome it is a difference
in log odds rather than in probability;
[`avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
is still what reports an effect on the scale a reader can act on.

What the reparameterization buys is a prior on the effect itself. The
forest for \\f_1\\ is regularized separately from the forest for the
rest of the outcome, so shrinking the prognostic part does not shrink
the effect.
[`vignette("causal")`](https://ngreifer.github.io/bartisan/articles/causal.md)
covers why that matters and
[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) sets it
up for the causal case.

Two things worth knowing before reaching for it.

The effect is **linear in the covariate** unless we say otherwise. For a
binary treatment that is no assumption at all, since there are only two
values. For a continuous predictor it says the effect is proportional to
it, which is a real restriction. Letting the coefficient’s forest split
on the covariate itself removes it, and then the effect varies across
the covariate’s own range:

``` r

# The effect of `aps` may itself change across `aps`.
# A modifier has to be a predictor of the model, so `aps` reaches the formula
# through `.`, which also keeps it out of the control function.
y ~ . + vc(aps, ~ aps + age)
```

And a covariate whose coefficient varies should not also be a predictor
of the control function. With it in both, the two are not separately
identified: any function of it can move between them. Writing the
covariate only inside
[`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) is what
keeps them apart. Named outright in the fixed part as well, the model is
fitted as asked and
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
warns; reached through `.`, the covariate is dropped from the control
function without comment, since `.` did not name it.

For a family with several additive predictors, each parameter’s formula
carries its own
[`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) terms, so
a covariate can have a coefficient on more than one of them:

``` r

# The effect of `z` on the mean, and separately on the spread.
bartisan(list(mean = y ~ x1 + x2 + vc(z),
              log_sd = ~ x1 + x2 + vc(z)),
         data = d, family = gaussian_ls())
```

[`coef()`](https://rdrr.io/r/stats/coef.html) then returns one column
per coefficient, named `mean:z` and `log_sd:z` for the forests they come
from, which is also how per-forest settings like `num_trees` are keyed.
[`?vc`](https://ngreifer.github.io/bartisan/reference/vc.md) covers the
rest, including the two multinomial families that refuse this.

## Descriptions of the Fitted Model

Everything here is a description of the fitted model.
[`avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
reports what the model predicts would differ between two versions of the
data, which is a causal quantity only if the model contains enough
covariates to account for confounding. Patients were not randomized to
catheterization, so for this fit that is a strong assumption.
[`vignette("causal")`](https://ngreifer.github.io/bartisan/articles/causal.md)
covers what is needed.

The intervals are posterior credible intervals under the model: they
cover the uncertainty in the fitted function, and they do not cover the
possibility that the model is missing a confounder, that the outcome is
measured with bias, or that the sample is not the population of
interest.

## Where to Go Next

[`vignette("importance")`](https://ngreifer.github.io/bartisan/articles/importance.md)
covers which predictors the forest uses, which is a different question
from how much they move the outcome.
[`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md)
covers whether the fit can be trusted before any of this is read.
`?bartisan-marginaleffects` documents which *marginaleffects* functions
are supported and the arguments that are specific to this package, and
[`?estimate_effect`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
and
[`?partial_dependence`](https://ngreifer.github.io/bartisan/reference/partial_dependence.md)
document the two questions answered without it.

## References

Arel-Bundock, Vincent, Noah Greifer, and Andrew Heiss. 2024. “How to
Interpret Statistical Models Using marginaleffects for R and Python.”
*Journal of Statistical Software* 111 (9): 1–32.
<https://doi.org/10.18637/jss.v111.i09>.

Deshpande, Sameer K., Ray Bai, Cecilia Balocchi, Jennifer E. Starling,
and Jordan Weiss. 2026. “VCBART: Bayesian Trees for Varying
Coefficients.” *Bayesian Analysis* 21 (1): 281–308.
<https://doi.org/10.1214/24-BA1470>.

Hahn, P. Richard, Jared S. Murray, and Carlos M. Carvalho. 2020.
“Bayesian Regression Tree Models for Causal Inference: Regularization,
Confounding, and Heterogeneous Effects (with Discussion).” *Bayesian
Analysis* 15 (3): 965–1056. <https://doi.org/10.1214/19-BA1195>.

Woody, Spencer, Carlos M. Carvalho, P. Richard Hahn, and Jared S.
Murray. 2020. *Estimating Heterogeneous Effects of Continuous Exposures
Using Bayesian Tree Ensembles: Revisiting the Impact of Abortion Rates
on Crime*. <https://arxiv.org/abs/2007.09845>.
