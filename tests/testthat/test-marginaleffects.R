# The point of these is not that marginaleffects works -- it is that the draws
# reach it correctly. Every estimand it computes is the same transformation
# applied to every draw, so if the "posterior_draws" attribute is transposed, or
# stacked in the wrong order for a categorical family, the numbers come out
# plausible and wrong.

skip_if_no_me <- function() {
  testthat::skip_if_not_installed("marginaleffects")
}

test_that("the base generics report what insight will ask them for", {
  d <- sim_x(n = 80, seed = 91)
  d$y <- d$x1 + stats::rnorm(nrow(d))

  fit <- bartisan(y ~ x1 + x2, d, control = quick_control())

  expect_s3_class(stats::formula(fit), "formula")
  expect_s3_class(stats::terms(fit), "terms")
  expect_identical(stats::nobs(fit), 80L)
  expect_s3_class(stats::model.frame(fit), "data.frame")
  # The model frame is the data the fit saw, which is what a counterfactual grid
  # is built from.
  expect_identical(nrow(stats::model.frame(fit)), 80L)
  expect_true(all(c("y", "x1", "x2") %in% names(stats::model.frame(fit))))
})

test_that("get_predict hands over draws in the orientation marginaleffects expects", {
  skip_if_no_me()

  d <- sim_x(n = 60, seed = 92)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))
  fit <- bartisan(y ~ x1 + x2, d, control = quick_control())

  out <- marginaleffects::get_predict(fit)

  expect_true(all(c("rowid", "group", "estimate") %in% names(out)))
  expect_identical(nrow(out), 60L)

  draws <- attr(out, "posterior_draws")
  # Observations by draws, not the other way round. Getting this backwards is
  # the failure that produces plausible nonsense.
  expect_identical(dim(draws), c(60L, 30L))

  # And the point estimate is the row mean of those draws, so the two agree
  # about which margin is which.
  expect_equal(out$estimate, rowMeans(draws), ignore_attr = TRUE)

  # The draws are the same numbers predict() gives, transposed.
  direct <- stats::predict(fit, newdata = d, type = "response", draws = TRUE)
  expect_equal(draws, t(direct), ignore_attr = TRUE)
})

test_that("a categorical family stacks rows and draws the same way", {
  skip_if_no_me()

  d <- sim_x(n = 60, seed = 93)
  set.seed(1093)
  z <- 2 * d$x1 - d$x2 + stats::rnorm(nrow(d))
  d$y <- ordered(rowSums(outer(z, stats::quantile(z, c(1, 2) / 3), ">")) + 1L)

  fit <- bartisan(y ~ x1 + x2, d, family = ordinal("probit"),
                  control = quick_control())

  out <- marginaleffects::get_predict(fit)
  draws <- attr(out, "posterior_draws")

  # Three categories, so three blocks of 60 rows.
  expect_identical(nrow(out), 180L)
  expect_identical(dim(draws), c(180L, 30L))
  expect_identical(as.character(out$group[c(1, 61, 121)]), c("1", "2", "3"))

  # Row order and draw order agree, which is the thing that would silently break.
  expect_equal(out$estimate, rowMeans(draws), ignore_attr = TRUE)

  # Probabilities over the categories sum to one for every row of every draw.
  by_cat <- array(draws, dim = c(60L, 3L, 30L))
  expect_equal(apply(by_cat, c(1L, 3L), sum), matrix(1, 60L, 30L),
               tolerance = 1e-10)
})

test_that("the estimands run, and their intervals come from the posterior", {
  skip_if_no_me()

  set.seed(94)
  n <- 300
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n),
                  z = factor(sample(c("a", "b"), n, TRUE)))
  d$y <- 2 * d$x1 - d$x2 + (d$z == "b") + stats::rnorm(n, sd = 0.5)

  fit <- bartisan(y ~ x1 + x2 + z, d, num_trees = 20L, num_burn = 200L,
                  num_draws = 200L)

  p <- marginaleffects::avg_predictions(fit)
  expect_identical(nrow(p), 1L)
  expect_true(p$conf.low < p$estimate && p$estimate < p$conf.high)

  cmp <- marginaleffects::avg_comparisons(fit, variables = "z")
  expect_identical(nrow(cmp), 1L)
  # The factor contrast is a whole unit in truth, and shrinkage pulls it in, so
  # this checks the sign and the order of magnitude rather than the value.
  expect_gt(cmp$estimate, 0.3)
  expect_lt(cmp$estimate, 1.5)

  # A grouped prediction splits by the factor and orders as the factor does.
  by_z <- marginaleffects::avg_predictions(fit, by = "z")
  expect_identical(nrow(by_z), 2L)
  expect_lt(by_z$estimate[1L], by_z$estimate[2L])
})

test_that("a multi-predictor family works on the response scale and not the link one", {
  skip_if_no_me()

  d <- sim_x(n = 80, seed = 95)
  d$y <- stats::rnorm(nrow(d))

  fit <- bartisan(y ~ x1 + x2, d, family = location_scale(),
                  control = quick_control())

  # The mean is one number per observation, so it goes through.
  out <- marginaleffects::get_predict(fit, type = "response")
  expect_identical(nrow(out), 80L)
  expect_identical(dim(attr(out, "posterior_draws")), c(80L, 30L))

  # The link is not, and saying so beats answering about the location and
  # leaving the scale unmentioned.
  expect_error(marginaleffects::get_predict(fit, type = "link"),
               "not a single quantity")

  # A zero-inflated family is the same shape of case.
  d$count <- stats::rpois(nrow(d), 2) * stats::rbinom(nrow(d), 1, 0.7)
  zi <- bartisan(count ~ x1 + x2, d, family = zi_poisson(),
                 control = quick_control())
  expect_identical(nrow(marginaleffects::get_predict(zi, type = "response")),
                   80L)
  expect_error(marginaleffects::get_predict(zi, type = "link"),
               "not a single quantity")
})

test_that("a multinomial fit gives one group per category", {
  skip_if_no_me()

  d <- sim_x(n = 60, seed = 97)
  set.seed(1097)
  d$y <- factor(sample(c("a", "b", "c"), nrow(d), TRUE))

  fit <- bartisan(y ~ x1 + x2, d, family = multinomial(),
                  control = quick_control())

  out <- marginaleffects::get_predict(fit)
  expect_identical(nrow(out), 180L)
  expect_identical(as.character(unique(out$group)), c("a", "b", "c"))
  expect_equal(out$estimate, rowMeans(attr(out, "posterior_draws")),
               ignore_attr = TRUE)
})

test_that("slopes are stable under a linear predictor transform and not under the quantile one", {
  skip_if_no_me()
  skip_on_cran()

  # The claim documented in ?bartisan-marginaleffects: the default transform is an
  # empirical distribution function, so the fit is a step function of the
  # original predictor and its difference quotient diverges as the step shrinks.
  set.seed(96)
  n <- 600
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  d$y <- 2 * d$x1 - d$x2 + stats::rnorm(n, sd = 0.5)

  slope_at <- function(transform, eps) {
    fit <- bartisan(y ~ x1 + x2, d, x_transform = transform, num_trees = 20L,
                    num_burn = 200L, num_draws = 200L)
    marginaleffects::avg_slopes(fit, variables = "x2", eps = eps)$estimate
  }

  linear_small <- slope_at("range", 1e-4)
  linear_big <- slope_at("range", 0.05)

  # Stable across a 500-fold change in the step, and of the right sign.
  expect_lt(abs(linear_small - linear_big), 0.25)
  expect_lt(linear_big, 0)

  # The quantile transform is not: at a small step the quotient is inflated by
  # roughly an order of magnitude.
  quantile_small <- slope_at("quantile", 1e-4)
  expect_gt(abs(quantile_small), 3 * abs(linear_small))
})

test_that("the ordinal mean and standardized latent scales are reachable", {
  skip_if_no_me()

  d <- sim_x(n = 150, seed = 101)
  set.seed(1101)
  latent <- 2 * d$x1 - d$x2 + stats::rlogis(nrow(d))
  d$y <- factor(findInterval(latent, c(-0.5, 0.8, 1.6)) + 1L,
                labels = c("1", "2", "3", "4"), ordered = TRUE)

  # Big enough to find x1: a contrast that is exactly zero because no tree ever
  # split on the predictor would pass the sign checks below without testing
  # anything.
  fit <- bartisan(y ~ x1 + x2, d, family = ordinal(),
                  control = quick_control(num_trees = 20L, num_burn = 100L,
                                          num_draws = 100L))

  # Both are one number per observation, so they come back ungrouped and with
  # the same draws predict() would give.
  for (type in c("mean", "stdlv")) {
    out <- marginaleffects::get_predict(fit, type = type)
    expect_identical(nrow(out), 150L)
    expect_identical(unique(out$group), "main_marginaleffect")
    expect_equal(out$estimate, stats::predict(fit, type = type))
    expect_identical(dim(attr(out, "posterior_draws")), c(150L, 100L))
  }

  # And they work through the estimands, which is the point of being reachable.
  # marginaleffects centers a posterior at its median, where predict() reports
  # its mean, so the two agree on the draws and not on the summary of them.
  average <- marginaleffects::avg_predictions(fit, type = "mean")
  drawn <- stats::predict(fit, type = "mean", draws = TRUE)
  expect_equal(average$estimate, stats::median(rowMeans(drawn)),
               tolerance = 1e-8)
  expect_true(average$conf.low < average$estimate)

  comparison <- marginaleffects::avg_comparisons(fit, type = "mean")
  expect_identical(nrow(comparison), 2L)
  expect_gt(comparison$estimate[comparison$term == "x1"], 0)

  # The standardized latent variable is a monotone transform of the predictor,
  # so a contrast on it has to agree in sign with one on the predictor.
  standardized <- marginaleffects::avg_comparisons(fit, type = "stdlv")
  on_link <- marginaleffects::avg_comparisons(fit, type = "link")
  expect_equal(sign(standardized$estimate), sign(on_link$estimate))
})

test_that("a binomial fit reaches the standardized latent scale too", {
  skip_if_no_me()

  d <- sim_x(n = 150, seed = 105)
  set.seed(1105)
  d$y <- stats::rbinom(nrow(d), 1L, stats::plogis(2 * d$x1 - d$x2))

  fit <- bartisan(y ~ x1 + x2, d,
                  control = quick_control(num_trees = 20L, num_burn = 100L,
                                          num_draws = 100L),
                 family = stats::binomial())

  out <- marginaleffects::get_predict(fit, type = "stdlv")
  expect_identical(nrow(out), 150L)
  expect_equal(out$estimate, stats::predict(fit, type = "stdlv"))
  expect_identical(dim(attr(out, "posterior_draws")), c(150L, 100L))

  # A contrast on the standardized scale is the contrast on the predictor divided
  # by the latent standard deviation, which for a logit link is at least
  # `pi / sqrt(3)`, so it is same-signed and no larger.
  standardized <- marginaleffects::avg_comparisons(fit, type = "stdlv")
  on_link <- marginaleffects::avg_comparisons(fit, type = "link")
  expect_equal(sign(standardized$estimate), sign(on_link$estimate))
  expect_true(all(abs(standardized$estimate) <= abs(on_link$estimate)))
  expect_gt(abs(standardized$estimate[standardized$term == "x1"]), 0)
})

test_that("values reaches predict, and the type vocabulary is checked", {
  skip_if_no_me()

  d <- sim_x(n = 120, seed = 103)
  set.seed(1103)
  latent <- 2 * d$x1 - d$x2 + stats::rlogis(nrow(d))
  d$y <- factor(findInterval(latent, c(0, 1)) + 1L,
                labels = c("1", "2", "3"), ordered = TRUE)

  fit <- bartisan(y ~ x1 + x2, d, family = ordinal(), control = quick_control())

  # `values` is this package's argument rather than one marginaleffects knows,
  # so it warns about not recognizing it -- and then passes it through, which is
  # what has to be true for the number to be right.
  worth <- c("1" = 0, "2" = 0, "3" = 100)
  out <- suppressWarnings(
    marginaleffects::avg_predictions(fit, type = "mean", values = worth))
  drawn <- stats::predict(fit, type = "mean", values = worth, draws = TRUE)
  expect_equal(out$estimate, stats::median(rowMeans(drawn)), tolerance = 1e-8)

  # And it is being used: without it the levels would be read as 1, 2 and 3.
  expect_gt(out$estimate, 10)

  # The aliases the same quantities go by elsewhere in marginaleffects.
  expect_equal(marginaleffects::get_predict(fit, type = "probs"),
               marginaleffects::get_predict(fit, type = "prob"))
  expect_equal(marginaleffects::get_predict(fit, type = "lp"),
               marginaleffects::get_predict(fit, type = "link"))

  # A type that is not a number per observation is refused by name rather than
  # failing somewhere downstream.
  expect_error(marginaleffects::get_predict(fit, type = "class"),
               "must be")
  expect_error(marginaleffects::get_predict(fit, type = "density"),
               "must be")
})

test_that("the estimand functions work for a survival response", {
  skip_on_cran()
  skip_if_not_installed("marginaleffects")

  # A survival response is a two-column matrix, which a model frame keeps as a
  # single matrix column. marginaleffects converts what it is handed to a
  # data.table, which cannot hold one, so every estimand failed for every
  # survival family -- whatever `type` was asked for -- with an error naming
  # data.table rather than the cause.
  # The interval checks below are coverage assertions at a single seed, so they
  # depend on how hard the prior shrinks this particular sample. Seed 43 sat on
  # the edge: the truth fell 0.004 inside the accelerated failure time interval,
  # against 0.03 to 0.10 at other seeds, so any change to the sampler flipped it.
  # Seed 3 has room on both families.
  set.seed(3)
  n <- 500
  d <- data.frame(x2 = stats::runif(n), trt = stats::rbinom(n, 1L, 0.5))
  event_time <- -log(stats::runif(n)) / exp(0.8 * d$trt)
  cens <- stats::runif(n, 0, 5)
  d$time <- pmin(event_time, cens)
  d$status <- as.numeric(event_time <= cens)

  ctrl <- bartisan_control(num_burn = 300, num_draws = 300, verbose = FALSE)

  for (family in list(ph(), weibull_aft())) {
    fit <- bartisan(cbind(time, status) ~ trt + x2, d, family = family,
                    control = ctrl)
    label <- family[["family"]]

    expect_no_error(marginaleffects::avg_comparisons(fit, variables = "trt"))

    # The estimand survival analysis is usually asked for: the difference in
    # t-year survival. The truth here is exp(-t e^0.8) - exp(-t).
    # marginaleffects checks the dots against a whitelist of its own, hardcoded
    # per model class, so it warns that it does not recognize `times` -- while
    # passing it through, which is what the warning says. There is no hook to
    # register an argument with it, so the warning is expected rather than a
    # fault, and is documented as such.
    at_one <- suppressWarnings(
      marginaleffects::avg_comparisons(fit, variables = "trt",
                                       type = "survival", times = 1))
    truth <- exp(-exp(0.8)) - exp(-1)

    # A point estimate shrunk toward zero is what BART's prior does, so the
    # check that matters is that the interval covers; the estimate only has to be
    # in the right neighbourhood. Measured errors here are 0.03 to 0.05.
    expect_lt(abs(at_one[["estimate"]] - truth), 0.08)
    expect_lt(at_one[["conf.low"]], truth)
    expect_gt(at_one[["conf.high"]], truth)

    # Survival falls further apart early and converges later, so the contrast at
    # a late time is smaller in magnitude than at an early one.
    at_two <- suppressWarnings(
      marginaleffects::avg_comparisons(fit, variables = "trt",
                                       type = "survival", times = 2))
    expect_lt(abs(at_two[["estimate"]]), abs(at_one[["estimate"]]))
  }
})

test_that("predictions() on newdata returns one row per row of newdata", {
  skip_if_not_installed("marginaleffects")

  # marginaleffects prepends rows of its own, marked `rowid = -1`, and drops them
  # again by that marker. Regenerating the column in get_predict() destroyed the
  # marker and let those rows into the result, which made `predictions()` return
  # more rows than it was given and biased every estimand built on it.
  d <- sim_x(n = 120, p = 3, seed = 991)
  set.seed(9911)
  d$g <- factor(sample(c("a", "b", "c"), nrow(d), replace = TRUE))
  d$y <- d$x1 + as.numeric(d$g) + stats::rnorm(nrow(d), sd = 0.3)

  fit <- bartisan(y ~ ., d, family = stats::gaussian(),
                  control = quick_control(num_draws = 60L))

  for (k in c(1L, 2L, 3L, 10L)) {
    got <- marginaleffects::predictions(fit, newdata = utils::head(d, k))
    expect_identical(nrow(as.data.frame(got)), k)
  }

  expect_identical(nrow(as.data.frame(marginaleffects::predictions(fit))),
                   nrow(d))

  grid <- marginaleffects::datagrid(model = fit, g = c("a", "b"))
  expect_identical(
    nrow(as.data.frame(marginaleffects::predictions(fit, newdata = grid))),
    nrow(grid))

  # The kept rows are the ones asked for, in order. marginaleffects summarizes
  # the draws with the median where predict() reports the mean, so the two agree
  # on the median rather than exactly.
  got <- as.data.frame(marginaleffects::predictions(fit, newdata = d[1:3, ]))
  drawn <- predict(fit, newdata = d[1:3, ], draws = TRUE)
  expect_equal(got$estimate, unname(apply(drawn, 2L, stats::median)),
               tolerance = 1e-8)
})
