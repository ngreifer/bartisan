# Evaluate stored forests at new data

Evaluate stored forests at new data

## Usage

``` r
.bartisan_predict(
  X,
  forest_flat,
  tree_start,
  bandwidth,
  num_forest,
  num_trees,
  num_draws,
  soft,
  gate,
  iterations,
  codes,
  tree_mask
)
```

## Arguments

- X:

  a design matrix with entries in `[0, 1]`.

- forest_flat, tree_start:

  the encoded forests returned by
  [`.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/dot-bartisan_fit.md).

- bandwidth:

  a matrix of per-tree bandwidths.

- num_forest, num_trees, num_draws:

  `integer`; the dimensions of the stored chain.

- soft:

  `logical`; whether the decision rules are soft.

- gate:

  `integer`; which gate the soft rules use; see `GateShape` in `node.h`.

- iterations:

  `integer`; the zero-based saved iterations to evaluate.

- codes:

  a matrix of level codes for the categorical rules.

- tree_mask:

  `logical`; one entry per stored tree, in the order
  [`.bartisan_tree_uses()`](https://ngreifer.github.io/bartisan/reference/dot-bartisan_tree_uses.md)
  reports, saying whether to evaluate it. A zero-length vector evaluates
  every tree, which is the ordinary case; a subset is what partial
  dependence uses to avoid re-evaluating the trees that cannot move
  across its grid.

## Value

A list of `num_forest` matrices of additive predictors.
