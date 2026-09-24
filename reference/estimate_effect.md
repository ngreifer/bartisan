# Causal effects from a fitted model

Estimates the average or conditional effect of a treatment by
g-computation over the posterior draws. `estimate_effect()` predicts
every unit under every level of the treatment, contrasts those potential
outcomes within each draw, and averages over the units the estimand asks
for.

## Usage

``` r
estimate_effect(
  object,
  treat = NULL,
  estimand = "ATE",
  comparison = "difference",
  by = NULL,
  newdata = NULL,
  level = 0.95,
  interval = "eti",
  focal = NULL,
  type = "response"
)

# S3 method for class 'bartisan_effect'
print(x, digits = 3L, contrasts = NULL, potential_outcomes = TRUE, ...)

# S3 method for class 'bartisan_effect'
plot(x, marginal = TRUE, ...)
```

## Arguments

- object:

  a `<bartisan_fit>` object; the output of a call to
  [`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) or
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).

- treat:

  `string`; the name of the treatment variable. A fit from
  [`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md)
  carries its own and needs none, so this is for a fit from
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md),
  where nothing marks one predictor as the treatment.

- estimand:

  `string`; which units the effect is averaged over. `"ATE"` (the
  default) uses all of them, and `"ATT"` and `"ATC"` those in the focal
  group, which by default is the treatment's second level for the former
  and its first for the latter. With a binary treatment that makes them
  the treated and the untreated without `focal` being named. `"CATE"`
  does not average at all and returns one effect per unit. Abbreviations
  and lowercase spellings are allowed.

- comparison:

  `string`; how the two potential outcomes are contrasted. Allowable
  options include `"difference"` (the default), `"ratio"`, `"lnratio"`,
  `"or"`, and `"lnor"`. The last two are available only when the
  response is a probability.

- by:

  optional; a one-sided formula or a variable name naming a grouping
  variable, in which case the effect is averaged within each of its
  levels rather than over the whole sample. A formula is evaluated with
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html)
  rather than read for the names it mentions, so it can define a
  grouping the data has no column for, as in `by = ~ age > 50` or
  `by = ~ interaction(sex, region)`. It must give exactly one grouping
  variable, and the term as written names the column it occupies in the
  output.

- newdata:

  optional; a data frame of units to average over. Default is the data
  the model was fit to, which is what makes the default estimand the
  sample average effect.

- level:

  `numeric`; the level of the credible interval. Default is `.95`.

- interval:

  `string`; `"eti"` (the default) for an equal-tailed interval from the
  quantiles of the draws, or `"hpdi"` for the highest posterior density
  interval, which is the shortest interval containing `level` of the
  posterior mass.

- focal:

  optional; the treatment level whose units the effect is averaged over
  for `"ATT"` and `"ATC"`. With a binary treatment it defaults to the
  second level for `"ATT"` and the first for `"ATC"`, which makes them
  the effect among the treated and among the untreated. With a treatment
  of more than two levels it has no default and must be supplied, and
  then `"ATT"` and `"ATC"` differ only in the level named.

- type:

  `string`; the prediction scale the effect is computed on, passed to
  [`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md).
  Default is `"response"`, which is the scale on which an average of
  unit-level differences is the marginal effect. See Details before
  changing it.

- x:

  a `<bartisan_effect>` object; the output of a call to
  `estimate_effect()`.

- digits:

  `integer`; the number of significant digits to print.

- contrasts:

  `string`; which contrasts to display when the treatment has more than
  two levels. `"focal"` (the default when a focal group is known) shows
  only the contrasts involving it, `"all"` shows every pairwise
  contrast, and a character vector of contrast labels shows those. All
  of them are computed either way; this only decides what is printed.

- potential_outcomes:

  `logical`; whether to print the average response under each treatment
  level below the contrasts, those being what the contrasts were
  computed from. Default is `TRUE`. They are in the result's
  `"potential_outcomes"` attribute either way.

- ...:

  ignored.

- marginal:

  `logical`; for
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on
  conditional effects (i.e., `estimand = "CATE"`), whether to draw the
  marginal effect beside the units. Default is `TRUE`. Ignored by the
  other plots, which draw no marginal effect beside their estimates.

## Value

A `<bartisan_effect>` object, a data frame with one row per reported
effect and columns `contrast`, `estimate`, `lower`, and `upper`, plus
`unit` when `estimand = "CATE"` and the `by` variable's name when `by`
is used. The posterior draws of every reported quantity are kept in the
`"draws"` attribute, and the marginal mean of each potential outcome in
`"potential_outcomes"`, so a caller can re-contrast or re-summarize
without refitting.

[`print()`](https://rdrr.io/r/base/print.html) returns its input
invisibly, and prints the potential outcomes below the contrasts unless
`potential_outcomes = FALSE`.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns a
ggplot2 object.

## Details

### Setting `type`

A varying coefficient is a contrast on the link scale: on a
`binomial("logit")` fit, `coef(object)` is a per-unit difference in log
odds, and the average of those is the average conditional log odds ratio
rather than the marginal one a treatment question usually asks for. The
default is therefore `type = "response"`, where every unit's contrast is
on the scale the response is measured on and averaging them gives the
marginal effect.

`type = "link"` is the right choice for looking at how the effect
varies, since that is the scale the forest models it on, and the wrong
one for reporting an average.

### Setting `comparison`

For `"ratio"`, `"lnratio"`, `"or"` and `"lnor"` the potential outcomes
are averaged over units first and contrasted afterward, which gives the
marginal ratio. With `estimand = "CATE"` there is no averaging to do, so
a ratio reported there is a conditional ratio and the
[`print()`](https://rdrr.io/r/base/print.html) method says so.

For the same reason `"or"` and `"lnor"` are not each other's
[`exp()`](https://rdrr.io/r/base/Log.html) and
[`log()`](https://rdrr.io/r/base/Log.html). Each summarizes the
posterior of the quantity it names, and a posterior mean does not
survive a nonlinear transformation: `"or"` reports the mean of the odds
ratio and `"lnor"` the mean of its logarithm, which exponentiates to
something smaller. Report whichever scale the interval should be
symmetric on, which for a ratio is usually the log. The same holds of
`"ratio"` against `"lnratio"`.

### Setting `focal`

With a binary treatment and no `focal`, which level is the treated one
is worked out from the treatment's values by the rules
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html)
uses, so that the same variable is read the same way by both packages.
In order: a `logical` treatment is treated at `TRUE`; a `numeric` one is
untreated at 0; a `character` or `factor` one whose levels parse as
numbers is untreated at 0, or else treated at the largest; and otherwise
the conventional names are matched, `"t"`, `"tr"`, `"treat"`,
`"treated"` and `"exposed"` against `"c"`, `"co"`, `"ctrl"`, `"control"`
and `"unexposed"`. When none of those applies the second level is taken
as the treated one and a message says so, that being the case where the
guess can be wrong. `"ATC"` takes the other level as its focal group, so
both estimands average over the group `focal` names and only the default
differs.

### Multi-category Treatments

There is nothing to work out from the values, so `focal` is required for
`"ATT"` and `"ATC"`, and those two then name the same estimand: the
effect among the units in the level named.

Every pairwise contrast is computed. Which ones are shown is a display
choice, made by `contrasts` in the
[`print()`](https://rdrr.io/r/base/print.html) method: with a focal
group the default is to show only the contrasts involving it, matching
how
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html)
and
[`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html)
use `focal`.

The effect of a continuous treatment is a slope or a dose-response curve
rather than a contrast of levels, and the
[adrftools](https://CRAN.R-project.org/package=adrftools) package has
tools for summarizing and visualizing one.

## See also

[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md), which
fits the model this is usually called on;
[`print.bcf_fit()`](https://ngreifer.github.io/bartisan/reference/print.bcf_fit.md)
for its other methods;
[`bartisan-marginaleffects`](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
for the same estimands through marginaleffects, which also covers the
ones not offered here

## Examples

``` r
data("rhc")
set.seed(123)

fit <- bcf(death ~ age + sex + meanbp + aps, treat = ~ rhc,
           data = rhc, num_trees = 10, num_burn = 50, num_draws = 50,
           verbose = FALSE)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` explicitly to silence this message.

# The risk difference, averaged over everyone
estimate_effect(fit)
#> Average treatment effect (difference)
#> 
#> Treatment: `rhc`
#> Averaged over 1500 units
#> 
#>     contrast estimate    lower  upper    n
#>  Y[1] - Y[0]   0.0334 -0.00348 0.0878 1500
#> 
#> Average potential outcomes
#> 
#>  quantity estimate lower upper
#>      Y[0]    0.642 0.610 0.665
#>      Y[1]    0.675 0.638 0.704
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ Y[a] is the average response with `rhc` set to "a".

# Among the treated, and as a risk ratio rather than a difference
estimate_effect(fit, estimand = "ATT", comparison = "ratio")
#> Average treatment effect on the treated (ratio)
#> 
#> Treatment: `rhc`
#> Averaged over the 565 units in group "1"
#> 
#>     contrast estimate lower upper   n
#>  Y[1] / Y[0]     1.04 0.994  1.12 565
#> 
#> Average potential outcomes
#> 
#>  quantity estimate lower upper
#>      Y[0]    0.670 0.634 0.693
#>      Y[1]    0.699 0.668 0.728
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ Y[a] is the average response with `rhc` set to "a".

# One effect per unit, ordered, with the marginal effect beside them
cate <- estimate_effect(fit, estimand = "CATE")
plot(cate)


# By a subgroup, which is where a forest plot earns its keep
estimate_effect(fit, by = ~ sex)
#> Average treatment effect (difference)
#> 
#> Treatment: `rhc`
#> Averaged over 1500 units
#> Within levels of `sex`
#> 
#>     sex    contrast estimate    lower  upper   n
#>  female Y[1] - Y[0]   0.0490 -0.00818 0.1270 676
#>    male Y[1] - Y[0]   0.0206 -0.01960 0.0634 824
#> 
#> Average potential outcomes
#> 
#>  quantity estimate lower upper
#>      Y[0]    0.642 0.610 0.665
#>      Y[1]    0.675 0.638 0.704
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ Y[a] is the average response with `rhc` set to "a".
```
