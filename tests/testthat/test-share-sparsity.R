# Pooling the splitting proportions across a family's forests. What makes this
# testable rather than a matter of taste is that the effect is visible in the
# splitting counts: a forest with no signal of its own splits on noise when it
# has its own Dirichlet, and on the other forest's predictors when they share
# one.

# Five predictors carry the mean and none carries the spread, so the scale
# forest of a `gaussian_ls()` fit has nothing of its own to find. Any predictor
# it concentrates on under a shared prior came from the mean forest.
sim_share <- function(n = 400L, p = 20L, seed = 1L) {
  set.seed(seed)
  x <- matrix(stats::runif(n * p), n, p)
  colnames(x) <- paste0("x", seq_len(p))
  mu <- 2 + 1.4 * sin(pi * x[, 1L] * x[, 2L]) + 1.6 * (x[, 3L] - 0.5)^2 +
    0.9 * x[, 4L] + 0.5 * x[, 5L]
  data.frame(y = stats::rnorm(n, mu, 0.5), x)
}

top_by_forest <- function(fit, k = 5L) {
  imp <- variable_importance(fit)
  lapply(split(imp, imp[["predictor"]]), function(part) {
    part[["variable"]][order(-part[["prop_splits"]])][seq_len(k)]
  })
}

test_that("share_sparsity is validated where it is set", {
  expect_false(bartisan_control()[["share_sparsity"]])
  expect_true(bartisan_control(share_sparsity = TRUE)[["share_sparsity"]])

  expect_error(bartisan_control(share_sparsity = 1), "logical")
  expect_error(bartisan_control(share_sparsity = c(TRUE, TRUE)), "logical")

  # Nothing to pool when the proportions are fixed rather than drawn, in either
  # of the two ways they can be fixed.
  expect_error(bartisan_control(sparsity = FALSE, share_sparsity = TRUE),
               "nothing to share")
  expect_error(bartisan_control(split_prior = c(x1 = 2), share_sparsity = TRUE),
               "nothing to share")

  # And it is allowed as long as at least one forest draws them, since the
  # engine shares among those and leaves the rest alone.
  expect_true(bartisan_control(sparsity = c(TRUE, FALSE),
                               share_sparsity = TRUE)[["share_sparsity"]])
})

test_that("one forest is nothing to share, and is not an error", {
  d <- sim_share(n = 150L, p = 5L, seed = 21L)

  fit <- bartisan(y ~ ., d, family = stats::gaussian(), share_sparsity = TRUE,
                  control = quick_control(num_trees = 5L))

  expect_s3_class(fit, "bartisan_fit")
  expect_length(fit[["counts"]], 1L)
})

test_that("a shared prior moves a signal-free forest onto the other's predictors", {
  skip_on_cran()

  d <- sim_share(n = 400L, p = 20L, seed = 22L)
  real <- paste0("x", 1:5)

  ctrl <- function(share) {
    quick_control(num_trees = 20L, num_burn = 300L, num_draws = 500L,
                  share_sparsity = share)
  }

  set.seed(2L)
  apart <- bartisan(y ~ ., d, family = gaussian_ls(), control = ctrl(FALSE))
  set.seed(2L)
  together <- bartisan(y ~ ., d, family = gaussian_ls(), control = ctrl(TRUE))

  a <- top_by_forest(apart)
  b <- top_by_forest(together)

  # Both fits should find the mean, which is where the signal is.
  expect_gte(length(intersect(a[["mean"]], real)), 4L)
  expect_gte(length(intersect(b[["mean"]], real)), 4L)

  # The scale forest is the test. Given its own prior it has no reason to prefer
  # the five that matter; sharing gives it one.
  expect_gt(length(intersect(b[["log_sd"]], real)),
            length(intersect(a[["log_sd"]], real)))
  expect_gte(length(intersect(b[["log_sd"]], real)), 4L)

  # Sharing a prior is not supposed to buy or cost fit on data this size, so a
  # large move in the log likelihood would mean something else changed.
  expect_equal(mean(together[["loglik"]]), mean(apart[["loglik"]]),
               tolerance = 0.05)
})

test_that("forests that may split on different predictors refuse to share", {
  skip_on_cran()

  d <- sim_share(n = 200L, p = 6L, seed = 23L)
  d$z <- stats::rbinom(nrow(d), 1L, 0.5)

  # A `vc()` term holds its coefficient forest to the modifiers, so with more
  # than one modifier the two forests draw their proportions over different
  # supports and a pooled Dirichlet over them would not be one distribution.
  expect_error(
    bartisan(y ~ vc(z, ~ x1 + x2) + x1 + x2 + x3, d, family = stats::gaussian(),
             control = quick_control(num_trees = 5L, share_sparsity = TRUE)),
    "same predictors")

  # One modifier is not that case, and is not an error: a forest held to a
  # single predictor has nothing to select between, so it does not draw its
  # proportions at all and there is nothing for it to share.
  expect_s3_class(
    bartisan(y ~ vc(z, ~ x1) + x1 + x2, d, family = stats::gaussian(),
             control = quick_control(num_trees = 5L, share_sparsity = TRUE)),
    "bartisan_fit")
})

test_that("sharing leaves the forests their own trees and leaf scales", {
  skip_on_cran()

  d <- sim_share(n = 300L, p = 10L, seed = 24L)

  fit <- bartisan(y ~ ., d, family = gaussian_ls(),
                  control = quick_control(num_trees = 20L, num_burn = 200L,
                                          num_draws = 300L,
                                          share_sparsity = TRUE))

  # Only the prior over predictors is common. The forests still differ in how
  # many rules they take and in their leaf scale, which is what distinguishes
  # this from sharing the topology.
  counts <- fit[["counts"]]
  expect_length(counts, 2L)
  expect_false(identical(counts[[1L]], counts[[2L]]))
  expect_gt(stats::sd(fit[["sigma_mu"]][, 1L] - fit[["sigma_mu"]][, 2L]), 0)

  expect_predictor_invariant(fit, d)
})
