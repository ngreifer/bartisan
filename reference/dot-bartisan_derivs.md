# Score and information of a family, analytic or by differences

Exists so that the test suite can check each family's analytic
derivatives against central differences of its own log density.

## Usage

``` r
.bartisan_derivs(
  y,
  weights,
  eta_draws,
  family_name,
  link,
  family_opts,
  aux,
  component,
  by_difference,
  blocked = FALSE
)
```

## Arguments

- y, weights, eta_draws, family_name, link, family_opts, aux:

  as for
  [`.bartisan_logdens()`](https://ngreifer.github.io/bartisan/reference/dot-bartisan_logdens.md).

- component:

  which additive predictor to differentiate with respect to.

- by_difference:

  use central differences instead of the analytic form.

- blocked:

  evaluate a whole draw at once through the family's block methods
  rather than one observation at a time. The two paths should agree;
  they differ for a family whose per-observation route falls back on
  differences while its block route does not.

## Value

A list with matrices `d1` and `info`, draws by observations.
