# Fit a generalized BART model

The workhorse behind
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).
Not intended to be called directly: it assumes the design matrix has
already been mapped to the unit interval and the response already
coerced to the form the requested family expects.

## Usage

``` r
.bartisan_fit(
  X,
  has_na,
  y,
  weights,
  offset,
  group_probs,
  family_name,
  link,
  family_opts,
  control,
  random_spec,
  codes,
  cat_col,
  n_levels,
  vc_basis
)
```

## Arguments

- X:

  a design matrix with entries in `[0, 1]`, possibly with `NA`.

- has_na:

  `logical`; one entry per column of `X` saying whether that column
  contains a missing value. A rule on a column with none is not given a
  missing-value branch, so complete data reproduces the sampler exactly
  as it was.

- y:

  the response, coerced by the calling family.

- weights:

  `numeric`; the prior weights.

- offset:

  an `H` by `N` matrix of fixed contributions to the additive
  predictors.

- group_probs:

  a sparse matrix whose columns are predictor groups.

- family_name, link, family_opts:

  the family specification.

- control:

  a list of sampler and prior settings.

## Value

A list of posterior draws and the encoded forests.
