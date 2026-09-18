# Give a predictor a varying coefficient

Marks a predictor inside a
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
formula as one whose coefficient is a function of the other predictors
rather than a constant, so that the coefficient gets a forest and a
prior of its own. `vc()` is a marker rather than a function: it is read
out of the formula it appears in and throws an error if it is called on
its own.

## Usage

``` r
vc(x, modifiers = NULL, center = "auto")
```

## Arguments

- x:

  the predictor whose coefficient varies, named as a bare variable
  rather than as an expression. A numeric variable gets one forest; a
  factor gets one per level, coded symmetrically (i.e., with no level
  held out as a reference) as
  [`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  codes its predictors.

- modifiers:

  a one-sided formula naming the predictors this coefficient's forest
  may split on. Default is `NULL` to allow every predictor in the model
  except `x` itself. Note that naming something that is not a predictor
  is an error rather than a silent restriction.

- center:

  the value of `x` at which the control function is read, given as
  either a string or a number. Allowable options include `"auto"` (the
  default), `"mean"`, `"zero"`, `"mid"`, and `"estimate"`. `"mean"`
  centers `x`, `"zero"` leaves it alone, `"mid"` uses the midpoint of
  its range, a number uses that number, and `"auto"` picks between
  `"zero"` and `"mean"` by the covariate. For a factor, `center` is
  `"mean"` or the name of a level to report against. `"estimate"` draws
  the coding rather than fixing it, and needs a covariate with between
  two and twenty distinct values. See Details.

## Value

Nothing; `vc()` is never evaluated, and calling it directly is an error.

## Details

The model is

\$\$g(\mu_i) = f_0(Z_i) + \sum_j (X\_{ij} - c_j) f_j(Z_i)\$\$

with a forest for the control function \\f_0\\ and one for each varying
coefficient \\f_j\\. Every forest is fitted at once, so the coefficient
has a prior of its own rather than being whatever difference a single
forest with `x` among its predictors happens to produce.

### Setting `modifiers`

By default a coefficient may vary with every predictor in the model
except `x` itself, and `modifiers` narrows that.

A covariate with a varying coefficient is kept out of the control
function, since \\f_0(Z) + z f_1(Z)\\ is not identified when `z` is
among \\f_0\\'s predictors: any function of `z` moves between the two.
Reached through `.` it is dropped quietly, because `.` did not name it;
named outright, the model is fitted as asked and a warning says why that
is a choice.

A numeric covariate may modify its own coefficient, and that is how the
effect stops being linear in it. `vc(z, ~ z + x1)` fits \\z f_1(z,
x_1)\\, so the slope moves across `z`'s range and the dose response is a
curve rather than a line, which is worth reaching for whenever the
effect of a continuous predictor might not be proportional to it. A
categorical covariate is removed from its own forests instead: a level's
indicator is nonzero only on the rows where that level holds, and the
variable is constant on exactly those rows, so such a split separates
rows that contribute from rows that contribute nothing.

At the other extreme, `~ 1` names nothing at all, which leaves the
coefficient's forest no predictor to split on: every tree in it is a
stump, so the coefficient is one drawn number and `x` enters as a linear
term while the rest of the model stays nonparametric. Comparing
`vc(z, ~ 1)` against `vc(z)` with
[loo()](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
is then a test of whether the effect of `z` varies at all, and the
constant fit reports it as a single coefficient, which under a logit
link is one conditional log odds ratio. It is drawn under the leaf prior
rather than a prior written for a regression coefficient, so it is
shrunk toward zero. See
[`vignette("comparison", package = "bartisan")`](https://ngreifer.github.io/bartisan/articles/comparison.md).

### Setting `center`

Centering is a reparameterization of \\f_0\\ alone: every coefficient
and every estimand is identical under any choice, and what changes is
what the control function means. `"auto"` picks by the covariate. For a
`0`/`1` covariate it uses zero, so \\f_0\\ is the surface among the
untreated, which is a quantity with its own meaning. For any other
numeric covariate it uses the mean, because zero may be nowhere near the
data and the control function there would be an extrapolation.

A factor is always fitted mean-centered and gets one forest per level,
coded symmetrically rather than as contrasts against a level that
happened to sort first. That coding leaves the reference a reporting
choice: `center` names the level
[`coef()`](https://rdrr.io/r/stats/coef.html) reports against, and no
refit is needed to change it.

### Drawing the Coding (`center = "estimate"`)

This is different in kind from the choices above. Rather than subtract a
number from `x`, it gives each of `x`'s values a coefficient of its own
and draws it, so that every contrast carries the same prior whatever the
number of values and no value is a reference.
[`coef()`](https://rdrr.io/r/stats/coef.html) returns the identified
contrasts.

At two values it restricts nothing, and it is what
[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) uses:
what it buys is that the answer stops depending on which level was
written as 1. Above two values every contrast becomes one shared shape
times a scalar, where the symmetric coding gives each level its own, so
it is the parsimonious model against a general one. Reach for it when
the levels plausibly differ in degree rather than in kind.

### Families with Several Additive Predictors

[`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
[`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
and the rest fit one forest per distributional parameter, and each
parameter's formula carries its own `vc()` terms. The forests are then
two-dimensional (a control function and its coefficients, for each
parameter) and named accordingly, which is what per-forest settings are
keyed by:

    # forests: mean, mean:z, log_sd
    bartisan(list(mean = y ~ x1 + x2 + vc(z), log_sd = ~ x1 + x2), data = d,
             family = gaussian_ls())

    # forests: mean, mean:z, log_sd, log_sd:z; one formula reaches every
    # parameter, which is the rule every per-forest argument follows
    bartisan(y ~ x1 + x2 + vc(z), data = d, family = gaussian_ls())

    # each coefficient with its own modifiers
    bartisan(list(mean = y ~ x1 + x2 + vc(z, ~ x2),
                  log_sd = ~ x1 + x2 + vc(z, ~ x1)), data = d,
             family = gaussian_ls())

So the same covariate may have a coefficient on more than one parameter:
\\z\\ shifting the mean and widening the spread are different questions,
and both are answered at once.
[`coef()`](https://rdrr.io/r/stats/coef.html) returns one column per
coefficient, named for its forest.

Two things follow from the parameters being different. A group intercept
from `(1 | g)` reaches every control function and no coefficient, since
a group-varying coefficient is a random slope. And `center = "estimate"`
is judged per parameter: the drawn coding needs a leaf target that is
quadratic in the predictor it feeds, which
[`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is in the mean and is not in the log standard deviation, so the same
request is accepted on one and refused on the other.

The two multinomial families are the exception and refuse `vc()`. Their
forests are the levels of one parameter rather than separate parameters,
identified only up to a function they all share, and reporting removes
it; a coefficient forest per level would add one such direction per
coefficient and the reporting does not carry them.

## See also

- [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
  for the formula interface

- [`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) for
  the causal case, which is this term with the priors and the propensity
  score set up for it

- [`coef.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/coef.bartisan_fit.md)
  for reading the coefficients out

- [`bartisan-families`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  for the order the forests come in

## Examples

``` r
data("rhc")

set.seed(123)

# The effect of right heart catheterization on death, free to vary with
# every other covariate. `rhc` reaches the fixed part through `.`, so it is
# dropped from the control function, which is what keeps the two identified
fit <- bartisan(death ~ . - days + vc(rhc), data = rhc,
                family = binomial(), num_trees = 10, num_burn = 50,
                num_draws = 50)

# One coefficient per patient, which is what a coefficient function comes to
head(coef(fit))
#>            rhc
#> [1,] 0.3852362
#> [2,] 0.2000175
#> [3,] 0.2619060
#> [4,] 0.2908669
#> [5,] 0.4621013
#> [6,] 0.3784168

# The same effect, free to vary with severity of illness alone
fit2 <- bartisan(death ~ . - days + vc(rhc, ~ aps), data = rhc,
                 family = binomial(), num_trees = 10, num_burn = 50,
                 num_draws = 50)

head(coef(fit2))
#>            rhc
#> [1,] 0.3902913
#> [2,] 0.2556783
#> [3,] 0.2023919
#> [4,] 0.3963084
#> [5,] 0.3731666
#> [6,] 0.3302474
```
