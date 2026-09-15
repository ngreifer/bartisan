# Methods for Bayesian causal forest fits

A fit from
[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) is a
`<bartisan_fit>` with a treatment named, so it takes every method that
fit takes. [`print()`](https://rdrr.io/r/base/print.html) adds the
treatment and says where the effect comes from, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
conditional effects.

## Usage

``` r
# S3 method for class 'bcf_fit'
print(x, digits = 3L, ...)

# S3 method for class 'bcf_fit'
plot(x, level = 0.95, comparison = "difference", ...)
```

## Arguments

- x:

  a `<bcf_fit>` object; the output of a call to
  [`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md).

- digits:

  `integer`; the number of significant digits to print.

- ...:

  for [`plot()`](https://rdrr.io/r/graphics/plot.default.html), further
  arguments passed to
  [`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md);
  otherwise ignored.

- level:

  `numeric`; the level of the credible interval. Default is `.95`.

- comparison:

  passed to
  [`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md),
  which is what computes the numbers.

## Value

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns a
ggplot2 object. [`print()`](https://rdrr.io/r/base/print.html) returns
its input invisibly.

## Details

The effect itself is
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)'s
to report, and [`summary()`](https://rdrr.io/r/base/summary.html) on a
`<bcf_fit>` is the same summary of the forests it is on any other fit,
with a line at the end saying so. That way the same call prints the same
thing whether the model came from
[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) or from
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
with a [`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md)
term.

## See also

[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
for the effect and its estimands;
[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) for the
model

## Examples

``` r
data("rhc")
set.seed(123)

fit <- bcf(death ~ age + sex + meanbp + aps, treat = ~ rhc,
           data = rhc, num_trees = 10, num_burn = 50, num_draws = 50,
           verbose = FALSE)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` to choose another, which also silences this message.

fit
#> Generalized BART
#> 
#> Call:
#> bcf(formula = death ~ age + sex + meanbp + aps, treat = ~rhc, 
#>     data = rhc, num_trees = 10, num_burn = 50, num_draws = 50, 
#>     verbose = FALSE)
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 1500
#> Structure: 2 forests of 10 trees, soft decision rules
#> Draws: 50 kept after 50 warmup
#> 
#> Posterior means: b.rhc.0 = 0.239, b.rhc.1 = 0.455
#> 
#> Treatment: "rhc"
#> Effect moderators: "age", "sex", "meanbp", and "aps"
#> ℹ `estimate_effect()` reports the treatment effect, with the average potential
#>   outcomes beside it; `plot()` draws the conditional ones.

# The effect, with the potential outcomes it is a difference of
estimate_effect(fit)
#> Average treatment effect (difference)
#> 
#> Treatment: "rhc"
#> Averaged over 1500 units
#> 
#>     contrast estimate    lower  upper    n
#>  Y[1] - Y[0]   0.0331 -0.00319 0.0902 1500
#> 
#> Average potential outcomes
#> 
#>  quantity estimate lower upper
#>      Y[0]    0.642 0.615  0.67
#>      Y[1]    0.675 0.647  0.71
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ Y[a] is the average response with "rhc" set to a.

# The conditional effects, ordered, with the average beside them
plot(fit)

```
