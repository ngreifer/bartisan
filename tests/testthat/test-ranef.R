# `ranef()` on a `(1 | group)` fit. The generic is nlme's, which lme4
# re-exports as the same function object, so one registration serves both names
# and both are checked here: the whole point of the method is that a
# mixed-model user's first guess reaches it.

skip_if_no_nlme <- function() {
  skip_if_not_installed("nlme")
}

sim_groups <- function(n = 400L, seed = 1L) {
  set.seed(seed)
  effect <- c(a = -1, b = -0.5, c = 0, d = 0.5, e = 1)
  d <- data.frame(x = stats::runif(n),
                  site = factor(sample(names(effect), n, replace = TRUE)))
  d$y <- stats::rnorm(n, 2 * d$x + effect[as.character(d$site)], 0.5)
  attr(d, "effect") <- effect
  d
}

group_control <- function(num_trees = 20L, num_burn = 100L,
                          num_draws = 200L, ...) {
  bartisan_control(num_trees = num_trees, num_burn = num_burn,
                   num_draws = num_draws, verbose = FALSE, ...)
}

test_that("the shape is the one lme4 returns, and both generics reach it", {
  skip_if_no_nlme()

  d <- sim_groups()
  fit <- bartisan(y ~ x + (1 | site), d, family = stats::gaussian(),
                  control = group_control())

  re <- nlme::ranef(fit)

  expect_type(re, "list")
  expect_named(re, "site")
  expect_s3_class(re[["site"]], "data.frame")

  # One row per level, in the order the levels are in, and one column for the
  # single additive predictor's intercept under lme4's name for it.
  expect_identical(rownames(re[["site"]]), levels(d$site))
  expect_named(re[["site"]], "(Intercept)")

  # lme4 re-exports nlme's generic rather than defining its own, so the two
  # names are the same function and one `@exportS3Method` covers both. If that
  # ever stops being true this is where it shows.
  skip_if_not_installed("lme4")

  expect_equal(lme4::ranef(fit), re)
})

test_that("the intercepts recover a known group effect", {
  skip_on_cran()
  skip_if_no_nlme()

  d <- sim_groups(n = 600L, seed = 2L)
  fit <- bartisan(y ~ x + (1 | site), d, family = stats::gaussian(),
                  control = group_control(num_trees = 50L, num_draws = 400L))

  got <- nlme::ranef(fit)[["site"]][["(Intercept)"]]

  # Tested by correlation and spread rather than by value: the intercepts are
  # deviations from an additive predictor that also fits `x`, and nothing
  # constrains them to sum to zero, so their level is not the thing to assert.
  expect_gt(stats::cor(got, attr(d, "effect")), 0.95)
  expect_lt(max(abs(got - mean(got) - attr(d, "effect"))), 0.4)
})

test_that("draws = TRUE gives the draws behind each posterior mean", {
  skip_if_no_nlme()

  d <- sim_groups(seed = 3L)
  fit <- bartisan(y ~ x + (1 | site), d, family = stats::gaussian(),
                  control = group_control())

  re <- nlme::ranef(fit)
  draws <- nlme::ranef(fit, draws = TRUE)

  expect_named(draws, "site")
  expect_named(draws[["site"]], "(Intercept)")

  m <- draws[["site"]][["(Intercept)"]]

  expect_true(is.matrix(m))
  expect_identical(ncol(m), length(levels(d$site)))
  expect_identical(colnames(m), levels(d$site))
  expect_identical(nrow(m), nrow(fit[["sigma_mu"]]))

  # The default is the mean of these, so the two must agree exactly.
  expect_equal(colMeans(m), re[["site"]][["(Intercept)"]],
               ignore_attr = TRUE)
})

test_that("each grouping factor gets its own entry, with only its own levels", {
  skip_if_no_nlme()

  set.seed(4); n <- 400L
  d <- data.frame(x = stats::runif(n),
                  school = factor(sample(letters[1:3], n, replace = TRUE)),
                  year = factor(sample(2001:2004, n, replace = TRUE)))
  d$y <- stats::rnorm(n, d$x)

  fit <- bartisan(y ~ x + (1 | school) + (1 | year), d,
                  family = stats::gaussian(), control = group_control())

  re <- nlme::ranef(fit)

  expect_named(re, c("school", "year"))
  expect_identical(rownames(re[["school"]]), levels(d$school))
  expect_identical(rownames(re[["year"]]), levels(d$year))
  expect_identical(nrow(re[["year"]]), 4L)
})

test_that("a family with several predictors gets a column per predictor", {
  skip_if_no_nlme()

  d <- sim_groups(seed = 5L)
  fit <- bartisan(y ~ x + (1 | site), d, family = gaussian_ls(),
                  control = group_control(num_trees = c(mean = 20L,
                                                        log_sd = 8L)))

  re <- nlme::ranef(fit)

  # The intercept is on each additive predictor and they are independent, so
  # `(Intercept)` would name two different things; the predictors name them.
  expect_named(re[["site"]], c("mean", "log_sd"))
  expect_identical(rownames(re[["site"]]), levels(d$site))

  expect_named(nlme::ranef(fit, draws = TRUE)[["site"]],
               c("mean", "log_sd"))
})

test_that("a fit with no bar in the formula says so rather than returning nothing", {
  skip_if_no_nlme()

  d <- sim_groups(n = 150L, seed = 6L)
  fit <- bartisan(y ~ x, d, family = stats::gaussian(),
                  control = group_control(num_trees = 5L, num_burn = 20L,
                                          num_draws = 30L))

  expect_error(nlme::ranef(fit), "no group intercepts")
})
