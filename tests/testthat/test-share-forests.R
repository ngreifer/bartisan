# Sharing the tree topology across a family's forests. The defining property is
# checkable exactly rather than statistically: if the t-th tree of every forest
# is the same tree, then every forest takes the same rules, so the per-draw
# split counts and the per-tree bandwidths must agree to the last bit -- and
# must not agree when the forests are separate. Almost everything below is that
# one fact asked in a different place.
#
# The model is implemented but not offered, so there is no argument to set and
# every fit here goes through `with_shared_forests()`, which flips the internal
# flag. The reasoning for keeping it internal is recorded beside that flag in
# `R/control.R`.

sim_share_topology <- function(n = 300L, p = 8L, seed = 1L) {
  set.seed(seed)
  x <- matrix(stats::runif(n * p), n, p)
  colnames(x) <- paste0("x", seq_len(p))
  mu <- 3 * sin(pi * x[, 1L]) + 2 * x[, 2L]
  lsd <- -0.9 + 0.9 * x[, 1L] + 0.7 * x[, 2L]
  data.frame(y = stats::rnorm(n, mu, exp(lsd)), x)
}

# The bandwidth matrix is draws by (trees times forests), laid out a forest at a
# time, so this is one forest's block.
bandwidth_block <- function(fit, h) {
  k <- fit[["num_trees"]][[h]]
  fit[["bandwidth"]][, (h - 1L) * k + seq_len(k), drop = FALSE]
}

test_that("sharing is off, and there is no argument that turns it on", {
  # The whole of its public surface is that there isn't one, so that is what
  # gets asserted: a caller who has read about the model and reaches for the
  # obvious spelling of it gets told there is no such setting, in both of the
  # places a control setting can be given.
  expect_false(bartisan_control()[["share_forests"]])
  expect_false("share_forests" %in% names(formals(bartisan_control)))
  expect_error(bartisan_control(share_forests = TRUE), "unused argument")
  expect_error(
    bartisan(y ~ ., sim_share_topology(n = 60L, p = 3L),
             family = stats::gaussian(), share_forests = TRUE,
             control = quick_control(num_trees = 3L)),
    "not an argument")

  # And the internal flag is what does reach the engine's setting.
  expect_true(with_shared_forests(bartisan_control()[["share_forests"]]))
})

test_that("one forest is nothing to share, and is not an error", {
  d <- sim_share_topology(n = 120L, p = 4L, seed = 11L)

  fit <- with_shared_forests(
    bartisan(y ~ ., d, family = stats::gaussian(),
             control = quick_control(num_trees = 5L)))

  expect_s3_class(fit, "bartisan_fit")
  expect_length(fit[["counts"]], 1L)
})

test_that("sharing gives the forests the same trees", {
  d <- sim_share_topology(n = 200L, p = 6L, seed = 12L)
  ctrl <- function() {
    quick_control(num_trees = 10L, num_burn = 60L, num_draws = 80L)
  }

  set.seed(3L)
  apart <- bartisan(y ~ ., d, family = gaussian_ls(), control = ctrl())
  set.seed(3L)
  together <- with_shared_forests(
    bartisan(y ~ ., d, family = gaussian_ls(), control = ctrl()))

  # The same rules mean the same counts of them, in every draw and for every
  # predictor. Nothing weaker would do: two forests can agree on totals by
  # accident, but not on a whole matrix.
  expect_identical(together[["counts"]][[1L]], together[["counts"]][[2L]])
  expect_false(identical(apart[["counts"]][[1L]], apart[["counts"]][[2L]]))

  # And the same gates, since a gate's width decides how the shared support
  # divides and two widths would be two partitions.
  expect_identical(bandwidth_block(together, 1L),
                   bandwidth_block(together, 2L))
  expect_false(identical(bandwidth_block(apart, 1L),
                         bandwidth_block(apart, 2L)))
})

test_that("sharing leaves the forests their own leaf values", {
  d <- sim_share_topology(n = 250L, p = 6L, seed = 13L)

  fit <- with_shared_forests(
    bartisan(y ~ ., d, family = gaussian_ls(),
             control = quick_control(num_trees = 10L, num_burn = 100L,
                                     num_draws = 150L)))

  # This is the whole of the flexibility the model keeps: one partition, a
  # different value in each piece for each component. A mean and a log standard
  # deviation that came out equal would mean the leaves had been shared too.
  expect_false(isTRUE(all.equal(fit[["eta"]][[1L]], fit[["eta"]][[2L]])))
  expect_gt(stats::sd(fit[["sigma_mu"]][, 1L] - fit[["sigma_mu"]][, 2L]), 0)
})

test_that("a shared fit's predictor is the one its trees replay", {
  skip_on_cran()

  d <- sim_share_topology(n = 200L, p = 5L, seed = 14L)

  # Across the gates, because the shared moves copy the membership weights
  # rather than recomputing them and a compact gate drops observations that an
  # unbounded one keeps.
  for (gate in c("smoothstep", "smootherstep", "logistic", "hard")) {
    fit <- with_shared_forests(
      bartisan(y ~ ., d, family = gaussian_ls(),
               control = quick_control(num_trees = 8L, num_burn = 50L,
                                       num_draws = 60L, gate = gate)))
    expect_predictor_invariant(fit, d)
  }

  # And with the bandwidth held fixed, which takes the shared bandwidth move
  # out of the picture and leaves only the tree moves.
  fit <- with_shared_forests(
    bartisan(y ~ ., d, family = gaussian_ls(),
             control = quick_control(num_trees = 8L, num_burn = 50L,
                                     num_draws = 60L,
                                     update_bandwidth = FALSE)))
  expect_predictor_invariant(fit, d)
})

test_that("a zero-inflated fit shares its two forests", {
  skip_on_cran()

  d <- sim_x(n = 250L, p = 6L, seed = 15L)
  d$y <- ifelse(stats::runif(nrow(d)) < stats::plogis(-0.6 + 1.6 * d$x2),
                0L, stats::rpois(nrow(d), exp(0.4 + 1.2 * d$x1)))

  fit <- with_shared_forests(
    bartisan(y ~ ., d, family = zi_poisson(),
             control = quick_control(num_trees = 10L, num_burn = 60L,
                                     num_draws = 80L)))

  expect_identical(fit[["counts"]][[1L]], fit[["counts"]][[2L]])
  expect_predictor_invariant(fit, d)
})

test_that("varying coefficients share with the control function", {
  skip_on_cran()

  d <- sim_x(n = 250L, p = 6L, seed = 16L)
  d$z1 <- stats::rnorm(nrow(d))
  d$z2 <- stats::rnorm(nrow(d))
  d$y <- 3 * sin(pi * d$x1) + d$z1 * (1 + d$x2) + d$z2 * (1 - d$x3) +
    stats::rnorm(nrow(d))

  # Three forests, all feeding one additive predictor: a control function and
  # two coefficient surfaces. Sharing applies to every forest in the model, not
  # only to forests that belong to different parameters, so all three take the
  # same partition of the modifier space.
  fit <- with_shared_forests(
    bartisan(y ~ x1 + x2 + x3 + x4 + x5 + x6 + vc(z1) + vc(z2), d,
             family = stats::gaussian(),
             control = quick_control(num_trees = 10L, num_burn = 60L,
                                     num_draws = 80L)))

  expect_length(fit[["counts"]], 3L)
  expect_identical(fit[["counts"]][[1L]], fit[["counts"]][[2L]])
  expect_identical(fit[["counts"]][[1L]], fit[["counts"]][[3L]])
})

test_that("settings that describe the trees have to agree to share them", {
  skip_on_cran()

  d <- sim_share_topology(n = 150L, p = 5L, seed = 17L)

  # The trees are one object, so a per-forest count of them has no meaning.
  expect_error(
    with_shared_forests(
      bartisan(y ~ ., d, family = gaussian_ls(),
               control = quick_control(num_trees = c(8L, 4L)))),
    "same number of trees")

  # One topology has one prior over its shape.
  expect_error(
    with_shared_forests(
      bartisan(y ~ ., d, family = gaussian_ls(),
               control = quick_control(num_trees = 5L,
                                       gamma = c(0.95, 0.5)))),
    "same .gamma")

  # And one set of splitting proportions to draw its rules from.
  expect_error(
    with_shared_forests(
      bartisan(y ~ ., d, family = gaussian_ls(),
               control = quick_control(num_trees = 5L,
                                       sparsity = c(TRUE, FALSE)))),
    "same .sparsity")

  # A gate's width decides how the shared support divides.
  expect_error(
    with_shared_forests(
      bartisan(y ~ ., d, family = gaussian_ls(),
               control = quick_control(num_trees = 5L,
                                       bandwidth = c(0.1, 0.3)))),
    "same .bandwidth")
})

test_that("forests that may split on different predictors refuse to share", {
  skip_on_cran()

  d <- sim_x(n = 200L, p = 6L, seed = 18L)
  d$z <- stats::rbinom(nrow(d), 1L, 0.5)
  d$y <- d$x1 + d$z * d$x2 + stats::rnorm(nrow(d))

  # A rule names a predictor, so forests with different predictors available to
  # them have no common set of rules to hold. A `vc()` term restricted to its
  # own modifiers is the way to produce that.
  expect_error(
    with_shared_forests(
      bartisan(y ~ vc(z, ~ x1 + x2) + x1 + x2 + x3, d,
               family = stats::gaussian(),
               control = quick_control(num_trees = 5L))),
    "same predictors")
})

test_that("a shared fit recovers a truth the components really do share", {
  skip_on_cran()

  d <- sim_share_topology(n = 400L, p = 10L, seed = 19L)
  test <- sim_share_topology(n = 400L, p = 10L, seed = 20L)
  x <- as.matrix(test[, paste0("x", 1:10)])
  mu <- 3 * sin(pi * x[, 1L]) + 2 * x[, 2L]

  fit <- with_shared_forests(
    bartisan(y ~ ., d, family = gaussian_ls(),
             control = quick_control(num_trees = 20L, num_burn = 300L,
                                     num_draws = 500L)))

  eta <- predict(fit, newdata = test, type = "link")
  expect_lt(sqrt(mean((eta[, "mean"] - mu)^2)), 0.6 * stats::sd(mu))

  # The log standard deviation is driven by the same two predictors, so the
  # scale forest should be splitting on them too -- which is what it gets from
  # the shared topology rather than from its own residual signal.
  imp <- variable_importance(fit)
  scale_part <- imp[imp[["predictor"]] == "log_sd", ]
  top <- scale_part[["variable"]][order(-scale_part[["prop_splits"]])][1:3]
  expect_gte(length(intersect(top, c("x1", "x2"))), 1L)
})
