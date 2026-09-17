# Which stored trees split on a given set of predictor columns

A sum of trees is separable, so a tree that never splits on a column
contributes the same amount however that column is set. Partial
dependence evaluates one column over a grid with the rest of the data
held fixed, so the trees this marks out are the only ones that have to
be re-evaluated at each grid point; the rest are evaluated once.

## Usage

``` r
.bartisan_tree_uses(
  forest_flat,
  tree_start,
  num_forest,
  num_trees,
  num_draws,
  num_cols,
  cat_cols
)
```

## Arguments

- forest_flat, tree_start:

  the encoded forests returned by
  [`.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/dot-bartisan_fit.md).

- num_forest, num_trees, num_draws:

  `integer`; the dimensions of the stored chain.

- num_cols:

  `integer`; the zero-based columns of the design matrix to look for,
  for rules without a level mask.

- cat_cols:

  `integer`; the zero-based columns of the level-code matrix to look
  for, for rules with one.

## Value

A `logical` vector with one entry per stored tree, in the order the flat
encoding holds them, which is iteration-major and then forest and then
tree.
