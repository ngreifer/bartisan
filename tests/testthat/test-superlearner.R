# `SL.bartisan()` and its predict method, the learner SuperLearner calls by
# name. SuperLearner supplies a fixed set of arguments, so these check that each
# is used the way it is meant to be: the family, the weights, and the
# predictors under whatever names they arrive with.

sim_learner <- function(n = 80, seed = 1) {
  d <- sim_x(n, 3, seed)
  d$y <- 2 * d$x1 - d$x2 + stats::rnorm(n, sd = 0.5)
  d$b <- as.numeric(d$y > mean(d$y))
  d
}

test_that("`SL.bartisan()` returns predictions its predict method reproduces", {
  d <- sim_learner()
  X <- d[c("x1", "x2", "x3")]

  out <- SL.bartisan(d$b, X, newX = X[1:10, ], family = binomial(),
                     obsWeights = rep(1, nrow(d)), id = seq_len(nrow(d)),
                     num_trees = 5, num_burn = 30, num_draws = 30)

  expect_s3_class(out[["fit"]], "SL.bartisan")
  expect_s3_class(out[["fit"]][["object"]], "bartisan_fit")
  expect_identical(out[["fit"]][["object"]][["family"]][["family"]], "binomial")

  expect_length(out[["pred"]], 10L)
  expect_true(all(out[["pred"]] > 0 & out[["pred"]] < 1))

  # Posterior means from one set of stored draws, so the second computation is
  # the first one repeated.
  expect_equal(stats::predict(out[["fit"]], newdata = X[1:10, ]),
               out[["pred"]], tolerance = 1e-10)
})

test_that("`SL.bartisan()` passes weights only when they are not all 1", {
  d <- sim_learner()
  X <- d[c("x1", "x2", "x3")]
  w <- stats::runif(nrow(d), 0.5, 2)

  out <- SL.bartisan(d$y, X, newX = X, family = gaussian(), obsWeights = w,
                     id = NULL, num_trees = 5, num_burn = 30, num_draws = 30)
  expect_equal(out[["fit"]][["object"]][["prior_weights"]], w)

  # `dpm()` takes no weights, and would be an error if the unit weights
  # SuperLearner passes by default were forwarded.
  out <- SL.bartisan(d$y, X, newX = X, family = gaussian(),
                     obsWeights = rep(1, nrow(d)), id = NULL,
                     bartisan_family = dpm(),
                     num_trees = 5, num_burn = 30, num_draws = 30)
  expect_identical(out[["fit"]][["object"]][["family"]][["family"]], "dpm")
})

test_that("`SL.bartisan()` handles a predictor named `Y` and a matrix `X`", {
  d <- sim_learner()
  X <- as.matrix(d[c("x1", "x2")])
  colnames(X) <- c("Y", "x2")

  out <- SL.bartisan(d$y, X, newX = X[1:5, ], family = gaussian(),
                     obsWeights = rep(1, nrow(d)), id = NULL,
                     num_trees = 5, num_burn = 30, num_draws = 30)

  expect_setequal(out[["fit"]][["object"]][["term_labels"]], c("Y", "x2"))
  expect_length(out[["pred"]], 5L)
})

test_that("`SL.bartisan()` can be named in a SuperLearner library", {
  skip_if_not_installed("SuperLearner")

  d <- sim_learner(60)
  X <- d[c("x1", "x2", "x3")]

  # SuperLearner looks learners up by name in `env`, the caller's frame by
  # default. Its own namespace sees its learners and screens (`SL.mean`, `All`)
  # and, through the search path, this package's exports, as a session that
  # has attached both would.
  set.seed(1301)
  sl <- SuperLearner::SuperLearner(d$y, X, family = gaussian(),
                                   SL.library = c("SL.mean", "SL.bartisan"),
                                   cvControl = list(V = 2L),
                                   env = asNamespace("SuperLearner"))

  expect_named(sl[["coef"]], c("SL.mean_All", "SL.bartisan_All"))
  expect_s3_class(sl[["fitLibrary"]][["SL.bartisan_All"]], "SL.bartisan")

  pred <- stats::predict(sl, newdata = X[1:5, ], onlySL = TRUE)
  expect_length(pred[["pred"]], 5L)
})
