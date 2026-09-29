# Group intercepts from a random-effect term

Extracts the intercept each level of a grouping factor was given by a
`(1 | group)` term in the formula, as a posterior mean or as every draw.

## Usage

``` r
# S3 method for class 'bartisan_fit'
ranef(object, draws = FALSE, ...)
```

## Arguments

- object:

  a `<bartisan_fit>` object; the output of a call to
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md),
  fitted with at least one `(1 | group)` term.

- draws:

  `logical`; whether to return every posterior draw of each intercept
  rather than its posterior mean. Default is `FALSE`.

- ...:

  not used.

## Value

With `draws = FALSE`, a named list with one entry per grouping factor,
each a data frame of one row per level of that factor and one column per
additive predictor that has a group intercept, with the levels as row
names. This is the shape
[`lme4::ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html)
returns, so anything that reads that reads this.

With `draws = TRUE`, a named list with one entry per grouping factor,
each itself a named list of draws-by-levels matrices, one per additive
predictor.

## Details

The column is named `(Intercept)`, as in lme4, when the family has a
single additive predictor. A family with more than one gets a group
intercept on each, independent of the others, and the columns are named
for the predictors instead;
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
lists them per family.

Only the intercepts are returned. The standard deviation each grouping
factor was drawn under is in `object$tau`, one column per factor and one
matrix per additive predictor, and
[`prior_summary()`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
reports the prior it was drawn from.

A posterior mean is the wrong summary for a level with few observations,
which is the case a group intercept exists for. Setting `draws = TRUE`
gives an interval, and a level whose interval covers zero is one the
data had little to say about.

The intercepts are shrunk towards zero by their prior and are deviations
from the additive predictor, so they come out approximately centered
without being constrained to sum to zero exactly. Nothing is lost by
that: the level of the fitted function is the additive predictor's, and
a shift common to every intercept is one the forest did not take.

## See also

[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
for the `(1 | group)` syntax and what it fits;
[`prior_summary()`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
for the prior on these

## Examples

``` r
set.seed(123)
sites <- factor(sample(letters[1:5], 200, TRUE))
d <- data.frame(x = runif(200), site = sites)
d$y <- rnorm(200, d$x + as.numeric(d$site) / 3)

fit <- bartisan(y ~ x + (1 | site), data = d,
                num_trees = 10, num_burn = 50,
                num_draws = 50, verbose = FALSE)
#> ℹ Using `family = dpm()`.
#> ℹ Set `family` explicitly to silence this message.

# One intercept per site, as posterior means. The
# generic is \pkg{nlme}'s, which \pkg{lme4} re-exports,
# so either qualification reaches this.
nlme::ranef(fit)
#> $site
#>   (Intercept)
#> a -0.41478565
#> b -0.10486625
#> c  0.38239000
#> d  0.08158241
#> e  0.62043215
#> 

# With the draws, so the intercepts come with intervals
re <- nlme::ranef(fit, draws = TRUE)
apply(re$site[["(Intercept)"]], 2L, quantile,
      c(.025, .975))
#>                 a          b            c          d         e
#> 2.5%  -0.72736163 -0.4163580 -0.006805168 -0.3149343 0.2263356
#> 97.5% -0.06789739  0.2126821  0.646165066  0.3988610 0.9472335
```
