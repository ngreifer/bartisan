# How often each predictor is used

Reports, for every predictor, how many splitting rules the forest spends
on it and how often it is used at all. This is the quantity people mean
by "variable importance" for a BART model, and it is the same table
[`summary.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/print.bartisan_fit.md)
prints – this returns it as a data frame instead, ready to sort, filter
or plot.

## Usage

``` r
variable_importance(object, level = 0.95)
```

## Arguments

- object:

  a fit from
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).

- level:

  the width of the interval reported for `splits`. Default `0.95`.

## Value

A data frame, one row per predictor, sorted by `prop_used` and then
`splits`, both decreasing. Columns are `variable`, `splits`,
`splits_lower`, `splits_upper` and `prop_used`. A family with more than
one additive predictor has a forest for each, and gains a leading
`predictor` column naming which.

## Which column answers which question

`prop_used` – the proportion of posterior draws in which the predictor
received at least one splitting rule – is the one to read first. It
behaves like a posterior probability that the predictor belongs in the
model, and it separates signal from noise sharply once `sparsity = TRUE`
in
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md),
which puts a Dirichlet prior on how the rules are shared out and lets
unused predictors be dropped rather than merely used rarely.

`splits` – the mean number of rules per draw – says how much of the
forest's structure a predictor accounts for. It is the more familiar
number and the easier one to over-read.

## Three things this is not

**It is not an effect size.** A predictor can be split on constantly and
move the prediction very little, and the reverse happens too. If the
question is how much a predictor moves the outcome, that is a job for
[`marginaleffects::avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
on the fitted model, not for this table. See
[bartisan-marginaleffects](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md).

**It is not stable under correlated predictors.** When two predictors
carry the same information the trees split on whichever is convenient,
and the usage distributes between them more or less arbitrarily. A
predictor can matter and still show a low `prop_used` because a
collinear partner absorbed it. Treat a group of correlated predictors as
a group.

**It is not causal.** A ranking of predictors by usage is a description
of this fitted function, not of what would happen if any of them were
changed.

## Reading it as variable selection

With `sparsity = TRUE`, `prop_used` is usable as a selection rule:
predictors the forest genuinely needs sit near 1 and the rest fall near
0, usually with a wide gap rather than a continuum. There is no
threshold that is correct in general; look at the gap and check that
your conclusion does not depend on where in it you cut.

## See also

[`summary.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/print.bartisan_fit.md),
which prints the same table;
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
for `sparsity`;
[bartisan-marginaleffects](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
for effects rather than usage.

## Examples

``` r
set.seed(1)
n <- 300
d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n), x4 = runif(n))

# Only x1 and x2 are in the truth.
d$y <- 2 * d$x1 + sin(3 * d$x2) + rnorm(n, sd = 0.3)

fit <- bartisan(y ~ ., data = d, family = gaussian(),
               control = bartisan_control(sparsity = TRUE))

variable_importance(fit)
#>   variable splits splits_lower splits_upper prop_used
#> 1       x2 39.416           19       68.000     1.000
#> 2       x1 26.568            8       45.525     1.000
#> 3       x3  7.948            0       25.000     0.926
#> 4       x4  2.038            0        8.000     0.576
```
