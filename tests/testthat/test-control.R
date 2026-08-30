test_that("bartisan_control validates its arguments", {
  expect_s3_class(bartisan_control(), "bartisan_control")

  expect_error(bartisan_control(num_trees = 0), "must be")
  expect_error(bartisan_control(num_trees = 2.5), "whole")
  expect_error(bartisan_control(gate = "yes"), "should be one of")
  expect_error(bartisan_control(sparsity = "yes"), "should be one of")

  # A vector of tree counts is accepted here; how many forests there are is not
  # known until the family is, so the length check belongs to `bartisan()`.
  expect_identical(bartisan_control(num_trees = c(10, 20))[["num_trees"]],
                   c(10L, 20L))
  expect_error(bartisan_control(bandwidth = -1), "must be")
  expect_error(bartisan_control(gamma = 1.5), "must be")
  expect_error(bartisan_control(sigma_mu = -1), "must be positive")
  expect_error(bartisan_control(sigma_mu_ramp = 2), "must be between")
  expect_error(bartisan_control(x_transform = "nope"), "should be one of")
})

test_that("bartisan rejects a control that did not come from bartisan_control", {
  d <- sim_x(seed = 31)
  d$y <- stats::rnorm(nrow(d))

  expect_error(bartisan(y ~ ., data = d, control = list(num_trees = 5)),
               "must be the result of")
})

test_that("the leaf scale defaults to the documented value", {
  d <- sim_x(seed = 32)
  d$y <- 5 * d$x1 + stats::rnorm(nrow(d))

  ctrl <- quick_control(num_trees = 8L, k = 2, update_sigma_mu = FALSE,
                        sigma_mu_ramp = 0)
  fit <- bartisan(y ~ ., data = d, control = ctrl)

  expected <- 3 * stats::sd(d$y) / (2 * sqrt(8))
  expect_equal(unname(fit[["sigma_mu"]][1L, 1L]), expected, tolerance = 1e-8)
})

test_that("both predictor transforms map into the unit interval", {
  d <- sim_x(seed = 33)
  d$x2 <- d$x2 * 1000 + 500
  d$y <- stats::rnorm(nrow(d))

  for (type in c("quantile", "range")) {
    fit <- bartisan(y ~ ., data = d, control = quick_control(x_transform = type))
    mapped <- vapply(seq_along(fit[["unit_maps"]]), function(j) {
      z <- fit[["unit_maps"]][[j]](d[[j]])
      all(z >= 0 & z <= 1)
    }, logical(1L))
    expect_true(all(mapped))
  }
})

test_that("a binary predictor maps to zero and one under either transform", {
  for (type in c("quantile", "range")) {
    f <- make_unit_map(c(0, 1, 1, 0), type)
    expect_identical(f(c(0, 1)), c(0, 1))
  }
})

test_that("a constant predictor maps to a single interior value", {
  f <- make_unit_map(rep(3, 10), "quantile")
  expect_identical(f(c(3, 3)), c(0.5, 0.5))
})

test_that("the sparsity prior is disabled when there is only one group", {
  d <- sim_x(p = 1, seed = 34)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  fit <- bartisan(y ~ x1, data = d, control = quick_control())

  expect_identical(ncol(fit[["counts"]][["eta"]]), 1L)
  expect_s3_class(fit, "bartisan_fit")
})

test_that("turning off the ramp leaves the leaf scale free to move", {
  d <- sim_x(seed = 35)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  fit <- bartisan(y ~ ., data = d,
                  control = quick_control(sigma_mu_ramp = 0,
                                          update_sigma_mu = TRUE))

  expect_true(stats::sd(fit[["sigma_mu"]][, 1L]) > 0)
})

test_that("control arguments can be passed through bartisan()'s dots", {
  d <- sim_x(seed = 41)
  d$y <- stats::rnorm(nrow(d))

  # A setting passed in `...` reaches the control list, and the ones passed to
  # bartisan_control() alongside it survive.
  fit <- bartisan(y ~ ., data = d, control = quick_control(),
                  num_trees = 3, gate = "hard")

  expect_identical(fit[["control"]][["num_trees"]], 3L)
  expect_false(fit[["control"]][["soft"]])
  expect_identical(fit[["control"]][["num_burn"]], 30L)

  # `...` wins over the same setting given in `control`.
  fit2 <- bartisan(y ~ ., data = d, control = quick_control(num_trees = 9),
                   num_trees = 4)

  expect_identical(fit2[["control"]][["num_trees"]], 4L)

  # Passing a setting through `...` is the same call as passing it through
  # `control`: same draws from the same seed.
  set.seed(1)
  a <- bartisan(y ~ ., data = d, control = quick_control(), num_trees = 3,
                gate = "hard")
  set.seed(1)
  b <- bartisan(y ~ ., data = d,
                control = quick_control(num_trees = 3, gate = "hard"))

  expect_equal(a[["eta"]], b[["eta"]])

  # An argument whose default is NULL is kept when passed explicitly, rather
  # than being dropped from the merged list.
  fit3 <- bartisan(y ~ ., data = d, control = quick_control(),
                   sigma_mu = NULL)
  expect_null(fit3[["control"]][["sigma_mu"]])
})

test_that("unknown dots arguments to bartisan() are an error", {
  d <- sim_x(seed = 42)
  d$y <- stats::rnorm(nrow(d))

  expect_error(bartisan(y ~ ., data = d, control = quick_control(),
                       num_tree = 3),
               "not an argument")
  expect_error(bartisan(y ~ ., data = d, control = quick_control(),
                       nonsense = 1, other = 2),
               "not arguments")

  # Validation still happens in bartisan_control(), rather than being skipped
  # for settings that arrive this way.
  expect_error(bartisan(y ~ ., data = d, control = quick_control(),
                       num_trees = -1),
               "num_trees")
})

test_that("`gate` covers both kinds of rule and `sparsity` sets the DART prior", {
  expect_true(bartisan_control()[["soft"]])
  expect_identical(bartisan_control()[["gate"]], "smoothstep")

  for (g in c("hard", "step")) {
    expect_false(bartisan_control(gate = g)[["soft"]])
  }

  for (g in c("logistic", "smoothstep", "smootherstep")) {
    expect_true(bartisan_control(gate = g)[["soft"]])
    expect_identical(bartisan_control(gate = g)[["gate"]], g)
  }

  # The negative binomial rewriting is a gain under hard rules only, so `gate`
  # has to reach `augment` for the default to mean what it says.
  expect_true("negbin" %in% bartisan_control(gate = "hard")[["augment"]])
  expect_false("negbin" %in% bartisan_control()[["augment"]])

  # The survival rewriting pays under both kinds of rule, so it is on either way.
  expect_true(all(c("aft") %in% bartisan_control()[["augment"]]))
  expect_identical(bartisan_control(augment = "aft")[["augment"]], "aft")

  on <- bartisan_control(sparsity = TRUE)
  off <- bartisan_control(sparsity = FALSE)

  expect_true(on[["update_s"]] && on[["update_alpha"]])
  expect_false(off[["update_s"]] || off[["update_alpha"]])
  expect_identical(off[c("update_s", "update_alpha")],
                   bartisan_control(sparsity = "none")[c("update_s",
                                                         "update_alpha")])
  expect_identical(on[["alpha_shape_2"]],
                   bartisan_control(sparsity = "moderate")[["alpha_shape_2"]])

  # Strength moves the prior on the concentration, and a knob given directly
  # still wins over the summary argument.
  expect_gt(bartisan_control(sparsity = "strong")[["alpha_shape_2"]],
            on[["alpha_shape_2"]])
  expect_lt(bartisan_control(sparsity = "weak")[["alpha_shape_1"]],
            bartisan_control(sparsity = "weak", alpha_shape_1 = 2)[["alpha_shape_1"]])
  expect_true(bartisan_control(sparsity = FALSE, update_s = TRUE)[["update_s"]])
})

test_that("a tree count per forest is honored, and the leaf scale follows it", {
  d <- sim_x(n = 200, seed = 61)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d), 0, exp(-1 + d$x2))

  fit <- bartisan(y ~ ., d, family = location_scale(),
                  control = quick_control(num_trees = c(20L, 5L)))

  expect_identical(fit[["num_trees"]], c(20L, 5L))

  # One bandwidth column per tree, the forests laid out back to back.
  expect_identical(ncol(fit[["bandwidth"]]), 25L)

  # One encoded tree per tree per draw.
  expect_identical(length(fit[["tree_start"]]),
                   25L * nrow(fit[["sigma_mu"]]) + 1L)

  # The stored forest replays to the reported predictor, which is what catches an
  # indexing error in the flat encoding.
  expect_predictor_invariant(fit, d)

  # A scalar is recycled, and recycling is the same call as saying it twice.
  set.seed(3)
  a <- bartisan(y ~ ., d, family = location_scale(),
                control = quick_control(num_trees = 8L))
  set.seed(3)
  b <- bartisan(y ~ ., d, family = location_scale(),
                control = quick_control(num_trees = c(8L, 8L)))

  expect_equal(a[["eta"]], b[["eta"]])

  # More values than forests is an error rather than a silent truncation.
  expect_error(bartisan(y ~ ., d, family = gaussian(),
                       control = quick_control(num_trees = c(5L, 5L))),
               "additive predictor")
})
