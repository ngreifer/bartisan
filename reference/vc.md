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

  the value of `x` at which the control function is read, which is both
  how the model is fitted and how it is reported, given as either a
  string or a number. Allowable options include `"auto"` (the default),
  `"mean"`, `"zero"`, `"mid"`, and `"estimate"`. `"auto"` uses `"zero"`
  for a `0`/`1` covariate and `"mean"` for any other numeric one, for
  the reason in Details. `"mean"` centers `x`, so the control function
  is the surface at its average; `"zero"` leaves `x` alone, so the
  control function is the surface at `x = 0`; `"mid"` uses the midpoint
  of `x`'s range (i.e., the average of its smallest and largest values);
  and a number uses that number. For a factor, `center` is `"mean"` or
  the name of a level to report against. `"estimate"` draws the coding
  rather than fixing it, which is the parameter expansion of Hahn,
  Murray and Carvalho (2020); see Details, and note that it needs a
  covariate with between two and twenty distinct values rather than a
  continuous one.

## Value

Nothing; `vc()` is never evaluated, and calling it directly is an error.

## Details

The model is

\$\$g(\mu_i) = f_0(Z_i) + \sum_j (X\_{ij} - c_j) f_j(Z_i)\$\$

with a forest for the control function \\f_0\\ and one for each varying
coefficient \\f_j\\. Every forest is fitted at once, so the coefficient
has a prior of its own rather than being whatever difference a single
forest with `x` among its predictors happens to produce.

### Which Predictors a Coefficient May Vary With

By default every predictor in the model except `x` itself. The
`modifiers` argument narrows that, and naming something that is not a
predictor is an error rather than a silent restriction.

**A varying covariate may not modify the control function.** With `z`
among \\f_0\\'s predictors, \\f_0(Z) + z f_1(Z)\\ is not identified: any
function of `z` moves between the two. Reached through `.`, the
covariate is dropped from the control function without comment, since
`.` did not name it. Named outright, the model is fitted as asked and a
warning says why that is a choice.

**A numeric covariate may modify its own coefficient, and doing so is
how the effect stops being linear.** `vc(z, ~ z + x1)` fits \\z f_1(z,
x_1)\\, so the slope itself moves across `z`'s range and the dose
response is a curve rather than a line. This is worth reaching for
whenever the effect of a continuous predictor might not be proportional
to it. On a simulation where the truth is \\y = 2x_1 + z^2\\, letting
\\f_1\\ split on `z` recovers the fitted surface to a root mean squared
error of 0.099 where forbidding it gives 1.187, and the fitted
coefficient traces the truth: -1.15, 0.15 and 1.09 at `z` of -1, 0 and
1, where the slope of \\z^2\\ is -1, 0 and 1.

**A categorical covariate may not**, and is removed from its own forests
quietly. A level's indicator is nonzero only on the rows where that
level holds, and the variable is constant on exactly those rows, so such
a split separates rows that contribute from rows that contribute
nothing. It is wasted rather than unidentified.

### Where the Control Function Sits

Centering is a reparameterization of \\f_0\\ alone: every coefficient
and every estimand is identical under any choice, and what changes is
what the control function means. `"auto"` picks by the covariate,
because neither answer wins everywhere. For a `0`/`1` covariate it uses
zero, so \\f_0\\ is the surface among the untreated (a quantity with its
own meaning, and the one that recovers the coefficient best). For any
other numeric covariate it uses the mean, because zero may be nowhere
near the data: with a covariate around 50 the control function at zero
is an extrapolation and recovery collapses to a correlation of 0.42.

A factor is always fitted mean-centered and gets one forest per level,
coded symmetrically the way
[`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
codes its predictors rather than as contrasts against a level that
happened to sort first. That coding carries one spare function-valued
dimension, which is what makes the reference a reporting choice:
`center` names the level [`coef()`](https://rdrr.io/r/stats/coef.html)
reports against, and no refit is needed to change it.

### Drawing the Coding Instead of Fixing It

`center = "estimate"` is different in kind from the choices above.
Rather than subtract a number from `x`, it gives each of `x`'s values a
coefficient of its own and draws it: the model is

\$\$g(\mu_i) = f_0(Z_i) + b\_{x_i}\\ \tilde f(Z_i), \qquad b_k \sim
\mathrm{N}(0, 1/2)\$\$

so the effect of moving from value \\j\\ to value \\k\\ is \\(b_k -
b_j)\tilde f(Z)\\, and because \\b_k - b_j\\ is marginally standard
normal every contrast carries the same prior whatever the number of
values and no value is a reference. This is the parameter expansion of
Hahn, Murray and Carvalho (2020, section 5.3), whose point is that a
fixed coding is never neutral: coding a treatment `0`/`1`, `1`/`0` or
\\\pm 1/2\\ makes the control function carry a different part of the fit
each time, and the prior on the effect moves with it.
[`coef()`](https://rdrr.io/r/stats/coef.html) returns the identified
contrasts, never the raw \\\tilde f\\.

It needs between two and twenty distinct values, and a continuous
covariate is refused rather than quietly binned.

**At two values it is free, and it is what
[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) uses.**
What it buys is that the answer stops depending on which level was
written as 1. Fitting the same data with the treatment coded `0`/`1` and
again `1`/`0` and adding the two effects, which is zero if the coding
does not matter, gives 0.0045 under a drawn coding against 0.0125 under
a fixed zero; on weak data, where the prior has more to say, 0.0389
against 0.0830.

**Above two values it stops being free, and the restriction is a real
one.** Every contrast is then the same \\\tilde f\\ times a scalar, so
the levels share one shape of heterogeneity where the symmetric coding
gives each its own. That is worth a great deal when it holds and
expensive when it does not: on a three-level covariate, `"estimate"`
takes about a quarter off the error when every level really does share a
shape, and nearly triples it when they do not.

So the two are what they look like: `"estimate"` is the parsimonious
model and the default is the general one. Reach for it when the levels
plausibly differ in degree rather than in kind.

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

## References

Hahn, P. R., Murray, J. S., & Carvalho, C. M. (2020). Bayesian
regression tree models for causal inference: regularization,
confounding, and heterogeneous effects. *Bayesian Analysis*, 15(3),
965–1056. [doi:10.1214/19-BA1195](https://doi.org/10.1214/19-BA1195)

## See also

- [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
  for the formula interface

- [`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) for
  the causal case, which is this term with the priors and the propensity
  score set up for it

- [`coef.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/coef.bartisan_fit.md)
  for reading the coefficients out

- [bartisan-families](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
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
