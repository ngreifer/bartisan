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
  codes
)
```

## Arguments

- X:

  design matrix with entries in `[0, 1]`.

- forest_flat, tree_start:

  the encoded forests returned by
  [`.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/dot-bartisan_fit.md).

- bandwidth:

  a matrix of per-tree bandwidths.

- num_forest, num_trees, num_draws:

  dimensions of the stored chain.

- soft:

  whether the decision rules are soft.

- gate:

  which gate the soft rules use; see `GateShape` in `node.h`.

- iterations:

  the zero-based saved iterations to evaluate.

## Value

A list of `num_forest` matrices of additive predictors.
