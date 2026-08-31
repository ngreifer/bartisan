# Splitting rules on a factor. A rule takes a subset of the levels still
# available at the node, rather than a threshold on one indicator column, which
# is what lets a tree put two levels in the same leaf.

sim_grouped <- function(n = 500, seed = 1) {
  set.seed(seed)
  g <- factor(sample(letters[1:5], n, TRUE))
  d <- data.frame(g = g, x1 = stats::rnorm(n))
  # {a, b} share a mean and {c, d, e} share another, which is a partition with
  # two cells of size two or more. One-hot encoding cannot form one.
  d$y <- ifelse(g %in% c("a", "b"), 2, -2) + stats::rnorm(n, sd = 0.5)
  d
}

test_that("a factor's levels are pooled into groups the data support", {
  d <- sim_grouped()

  fit <- bartisan(y ~ g + x1, data = d, family = gaussian(),
                  control = bartisan_control(num_trees = 20, num_burn = 200,
                                             num_draws = 200, verbose = FALSE))

  by_level <- tapply(fitted(fit), d$g, mean)
  high <- unname(by_level[c("a", "b")])
  low <- unname(by_level[c("c", "d", "e")])

  # The claim is the grouping, not the level means, which the leaf prior shrinks
  # towards the overall mean: the levels that share a mean are estimated close
  # together, and the two groups are far apart.
  expect_lt(max(high) - min(high), 0.25)
  expect_lt(max(low) - min(low), 0.25)
  expect_gt(min(high) - max(low), 3)
})

test_that("the rule reaches partitions one-hot encoding cannot form", {
  skip_on_cran()

  # Read off a one-tree forest on pure noise, where the posterior is close to the
  # prior: two levels are in the same cell exactly when the tree predicts the
  # same value for them.
  set.seed(4)
  n <- 600
  K <- 5
  d <- data.frame(y = stats::rnorm(n),
                  g = factor(sample(letters[seq_len(K)], n, TRUE)))

  fit <- bartisan(y ~ g, data = d, family = gaussian(),
                  control = bartisan_control(num_trees = 1, num_burn = 300,
                                             num_draws = 1500, gate = "hard",
                                             sparsity = FALSE, verbose = FALSE))

  one <- d[match(levels(d$g), d$g), , drop = FALSE]
  draws <- predict(fit, newdata = one, type = "link", summary = FALSE,
                   draws = TRUE)

  # A partition one-hot cannot form is one with two or more cells holding two or
  # more levels each, since a threshold on an indicator peels off one level and
  # leaves the rest in a single cell.
  two_big_cells <- apply(draws, 1L, function(v) {
    sum(table(match(v, unique(v))) >= 2) >= 2
  })

  expect_gt(mean(two_big_cells), 0.1)
})

test_that("categorical = onehot restores the threshold-on-an-indicator rule", {
  set.seed(5)
  n <- 500
  K <- 5
  d <- data.frame(y = stats::rnorm(n),
                  g = factor(sample(letters[seq_len(K)], n, TRUE)))

  fit <- bartisan(y ~ g, data = d, family = gaussian(),
                  control = bartisan_control(num_trees = 1, num_burn = 200,
                                             num_draws = 1000, gate = "hard",
                                             sparsity = FALSE,
                                             categorical = "onehot",
                                             verbose = FALSE))

  one <- d[match(levels(d$g), d$g), , drop = FALSE]
  draws <- predict(fit, newdata = one, type = "link", summary = FALSE,
                   draws = TRUE)

  two_big_cells <- apply(draws, 1L, function(v) {
    sum(table(match(v, unique(v))) >= 2) >= 2
  })

  # Not merely rare: impossible.
  expect_identical(sum(two_big_cells), 0L)
})

test_that("prediction on new data reproduces the fit", {
  d <- sim_grouped(seed = 6)

  fit <- bartisan(y ~ g + x1, data = d, family = gaussian(),
                  control = bartisan_control(num_trees = 10, num_burn = 100,
                                             num_draws = 100, verbose = FALSE))

  # The level codes of `newdata` are built from the structure the fit recorded,
  # so a rule that names a set of levels finds the same levels again.
  expect_equal(as.numeric(predict(fit, newdata = d)), as.numeric(fitted(fit)),
               tolerance = 1e-10)

  # And a subset of the rows, where a level may be absent altogether.
  small <- d[d$g %in% c("a", "c"), ]
  expect_equal(as.numeric(predict(fit, newdata = small)),
               as.numeric(fitted(fit))[d$g %in% c("a", "c")],
               tolerance = 1e-10)
})

test_that("a rejected change of rule does not leave a stale level set", {
  # The change move restores the old rule when its proposal is rejected. The
  # level set has to be part of what it restores: without it a node could end up
  # with a numeric rule's `var`, which indexes X, and a categorical rule's level
  # set, which says to index the matrix of level codes instead. That read was
  # out of bounds on the first fit that mixed a factor with a numeric predictor.
  set.seed(7)
  n <- 300
  d <- data.frame(y = stats::rnorm(n),
                  g = factor(sample(letters[1:5], n, TRUE)),
                  x1 = stats::rnorm(n))

  expect_no_error(bartisan(y ~ g + x1, data = d, family = gaussian(),
                           control = bartisan_control(num_burn = 200,
                                                      num_draws = 200,
                                                      verbose = FALSE)))
})

test_that("a two-level factor and a missing level are both handled", {
  set.seed(8)
  n <- 300
  d <- data.frame(y = stats::rnorm(n),
                  b = factor(sample(c("no", "yes"), n, TRUE)),
                  g = factor(sample(letters[1:4], n, TRUE)))

  # Two levels lose nothing either way, since 2^2 - 2 is the second Bell number.
  expect_no_error(bartisan(y ~ b, data = d, family = gaussian(),
                           control = quick_control()))

  d$g[1:20] <- NA

  fit <- bartisan(y ~ b + g, data = d, family = gaussian(),
                  na.action = stats::na.pass, control = quick_control())

  expect_identical(nrow(fit[["model"]]), 300L)
  expect_no_error(predict(fit, newdata = d))
})

test_that("only groups whose columns are mutually exclusive indicators qualify", {
  d <- data.frame(g = factor(c("a", "b", "c", "a", "b")), z = c(1, 2, 3, 4, 5))
  mt <- stats::terms(~ g + z + g:z, data = d)
  mf <- stats::model.frame(mt, d)
  des <- build_design(mt, mf)
  info <- level_codes(des$x, des$assign)

  # The factor qualifies. The numeric does not, and neither does the factor
  # crossed with it, whose columns are not indicators.
  expect_identical(info[["n_levels"]], c(3L, 0L, 0L))
  expect_identical(info[["cat_col"]], c(0L, -1L, -1L))
  expect_identical(as.integer(info[["codes"]]), c(0L, 1L, 2L, 0L, 1L))
})
