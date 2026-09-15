# `prior_summary()` reports the prior the engine was given rather than the one
# the caller wrote, which are not the same object: the per-forest arguments have
# been spread to one value each, `k` and `sigma_mu` are two ways of writing the
# same thing, and the default is neither of them. The fit carries a record of
# what the engine received, and these tests check the record against the formula
# it is supposed to be.

skip_if_no_rstantools <- function() {
  skip_if_not_installed("rstantools")
}

sim_prior <- function(n = 150L, seed = 1L) {
  set.seed(seed)
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  d$y <- stats::rnorm(n, d$x1)
  d
}

test_that("the leaf scale is the documented function of k and the tree count", {
  skip_if_no_rstantools()

  d <- sim_prior()
  fit <- bartisan(y ~ x1 + x2, d, family = stats::gaussian(),
                  control = quick_control(num_trees = 20L, k = 3,
                                          num_burn = 20L, num_draws = 30L))

  ps <- rstantools::prior_summary(fit)
  row <- ps[["forests"]]

  expect_s3_class(ps, "bartisan_prior_summary")
  expect_identical(nrow(row), 1L)

  # `?bartisan_control` documents the scale as 3 * s / (k * sqrt(num_trees)),
  # with `s` the response's scale on the link scale.
  expect_equal(row[["leaf_scale"]],
               3 * ps[["response"]][["scale"]] / (3 * sqrt(20)))
  expect_identical(row[["k"]], 3)
  expect_identical(row[["num_trees"]], 20L)

  # A Gaussian response is anchored at its own mean, and the scale it is
  # measured on is its standard deviation.
  expect_equal(ps[["response"]][["scale"]], stats::sd(d$y))
  expect_equal(ps[["response"]][["intercept"]], mean(d$y))
})

test_that("a per-forest argument gives one row per forest, and only it moves", {
  skip_if_no_rstantools()

  d <- sim_prior(seed = 2L)
  fit <- bartisan(y ~ x1 + x2, d, family = gaussian_ls(),
                  control = quick_control(num_trees = c(mean = 20L, log_sd = 6L),
                                          gamma = c(log_sd = 0.5),
                                          num_burn = 20L, num_draws = 30L))

  row <- rstantools::prior_summary(fit)[["forests"]]

  expect_identical(row[["forest"]], c("mean", "log_sd"))
  expect_identical(row[["num_trees"]], c(20L, 6L))

  # The forest the caller did not name keeps the argument's own default rather
  # than borrowing the one they did give, which is what `?bartisan_control`
  # promises and what the per-forest spreading got wrong once.
  expect_identical(row[["gamma"]], c(0.95, 0.5))
})

test_that("the concentration's scale is resolved the way the engine resolves it", {
  skip_if_no_rstantools()

  d <- sim_prior(seed = 3L)
  fit <- bartisan(y ~ x1 + x2, d, family = stats::gaussian(),
                  control = quick_control(num_trees = 5L, num_burn = 20L,
                                          num_draws = 30L))

  row <- rstantools::prior_summary(fit)[["forests"]]

  # `alpha_scale` is passed to the engine as zero to mean "the number of
  # predictors this forest may split on", and `Hypers()` substitutes that count;
  # reporting the zero would be reporting a sentinel rather than a prior.
  expect_identical(row[["candidates"]], 2L)
  expect_identical(row[["alpha_scale"]], 2L)

  # `sparsity = "moderate"`, the default, is Beta(0.5, 1).
  expect_identical(row[["shape_1"]], 0.5)
  expect_identical(row[["shape_2"]], 1)
})

test_that("the family's own parameters carry the prior each is given", {
  skip_if_no_rstantools()

  d <- sim_prior(seed = 4L)
  ctrl <- quick_control(num_trees = 5L, num_burn = 20L, num_draws = 30L)

  # A Gaussian scale is half-Cauchy at a scale read off the response.
  gaussian <- rstantools::prior_summary(
    bartisan(y ~ x1, d, family = stats::gaussian(), control = ctrl))[["family"]]

  expect_identical(gaussian[["name"]], "gaussian")
  expect_identical(gaussian[["parameters"]][["parameter"]], "sigma")
  expect_match(gaussian[["parameters"]][["prior"]], "HalfCauchy")

  # A location-scale family has no nuisance parameter at all: the scale is the
  # second forest, and is in the table rather than here.
  ls_family <- rstantools::prior_summary(
    bartisan(y ~ x1, d, family = gaussian_ls(), control = ctrl))[["family"]]

  expect_null(ls_family[["parameters"]])

  # The gamma-prior families are recognized from their option names, and
  # `shape_prior_shape` must not also read as a prior on something called
  # `shape_prior`.
  d$pos <- stats::rexp(nrow(d))
  gamma <- rstantools::prior_summary(
    bartisan(pos ~ x1, d, family = Gamma("log"), control = ctrl))[["family"]]

  expect_identical(gamma[["parameters"]][["parameter"]], "shape")
  expect_match(gamma[["parameters"]][["prior"]], "^Gamma\\(shape = ")
})

test_that("a cutpoint with no prior is reported as having none", {
  skip_if_no_rstantools()

  d <- sim_prior(seed = 5L)
  d$ord <- factor(sample(1:3, nrow(d), replace = TRUE), ordered = TRUE)
  ctrl <- quick_control(num_trees = 5L, num_burn = 20L, num_draws = 30L)

  params <- rstantools::prior_summary(
    bartisan(ord ~ x1, d, family = ordinal(),
             control = ctrl))[["family"]][["parameters"]]

  cuts <- params[params[["parameter"]] == "cutpoints", ]

  expect_identical(nrow(cuts), 1L)
  expect_match(cuts[["prior"]], "none")
})

test_that("the report prints every block it has, and none it has not", {
  skip_if_no_rstantools()

  d <- sim_prior(seed = 6L)
  d$g <- factor(sample(letters[1:4], nrow(d), replace = TRUE))

  # Soft rules and a random part, so both of the optional blocks are present.
  soft <- bartisan(y ~ x1 + x2 + (1 | g), d, family = stats::gaussian(),
                   control = quick_control(num_trees = 5L, num_burn = 20L,
                                           num_draws = 30L))

  out <- capture.output(print(rstantools::prior_summary(soft)))
  text <- paste(out, collapse = " ")

  for (heading in c("Trees", "Leaves", "Splitting variables",
                    "Decision rules", "Group intercepts", "Family")) {
    expect_match(text, heading, fixed = TRUE)
  }

  expect_match(text, "smoothstep", fixed = TRUE)

  # Hard rules have no bandwidth to report and no formula bar has no group
  # intercepts, so neither block should claim one.
  hard <- bartisan(y ~ x1 + x2, d, family = stats::gaussian(),
                   control = quick_control(num_trees = 5L, gate = "hard",
                                           num_burn = 20L, num_draws = 30L))

  bare <- paste(capture.output(print(rstantools::prior_summary(hard))),
                collapse = " ")

  expect_match(bare, "no bandwidth is used", fixed = TRUE)
  expect_no_match(bare, "Group intercepts", fixed = TRUE)

  # One forest, one set of numbers, and the table would say each of them twice.
  expect_no_match(bare, "num_trees", fixed = TRUE)
})

test_that("a fixed splitting prior is reported as fixed rather than drawn", {
  skip_if_no_rstantools()

  d <- sim_prior(seed = 7L)
  ctrl <- quick_control(num_trees = 5L, num_burn = 20L, num_draws = 30L)

  fixed <- bartisan(y ~ x1 + x2, d, family = stats::gaussian(),
                    control = ctrl, split_prior = c(x1 = 2, x2 = 1))

  ps <- rstantools::prior_summary(fixed)

  expect_true(ps[["split_prior"]])
  expect_match(paste(capture.output(print(ps)), collapse = " "),
               "Fixed by", fixed = TRUE)

  none <- bartisan(y ~ x1 + x2, d, family = stats::gaussian(),
                   control = ctrl, sparsity = FALSE)

  expect_false(rstantools::prior_summary(none)[["sparsity"]])
  expect_match(paste(capture.output(print(rstantools::prior_summary(none))),
                     collapse = " "),
               "equally likely", fixed = TRUE)
})

test_that("a prior-only fit says so, since its draws are what this describes", {
  skip_if_no_rstantools()

  d <- sim_prior(seed = 8L)
  fit <- bartisan(y ~ x1 + x2, d, family = stats::gaussian(),
                  prior_only = TRUE,
                  control = quick_control(num_trees = 5L, num_burn = 20L,
                                          num_draws = 30L))

  ps <- rstantools::prior_summary(fit)

  expect_true(ps[["prior_only"]])
  expect_match(paste(capture.output(print(ps)), collapse = " "),
               "itself a draw from the prior", fixed = TRUE)
})
