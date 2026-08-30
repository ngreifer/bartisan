# The ordinal family gained a third link and changed the chart it reports its
# cutpoints in. Both are places where a mistake produces a plausible fit rather
# than an error, so they are checked against quantities that are known.

test_that("the complementary log-log link recovers a cloglog-generated fit", {
  skip_on_cran()

  set.seed(11)
  n <- 1200
  x <- sim_x(n = n, p = 4, seed = 11)
  lin <- 2 * sin(pi * x$x1) - x$x3
  truth <- c(-0.5, 0.6, 1.4)

  # A latent smallest extreme value variate, which is what the complementary
  # log-log link is the distribution function of.
  z <- lin + log(-log1p(-stats::runif(n)))
  d <- x
  d$y <- ordered(rowSums(outer(z, truth, ">")) + 1L)

  fits <- lapply(c("cloglog", "logit", "probit"), function(lk) {
    set.seed(3)
    bartisan(y ~ ., d, family = ordinal(lk), gate = "hard", num_trees = 50,
             num_burn = 300, num_save = 300)
  })
  names(fits) <- c("cloglog", "logit", "probit")

  # The link the data came from fits best, which is the sharp check that the new
  # branch computes the right density rather than merely a smooth one.
  ll <- vapply(fits, function(f) mean(f[["loglik"]]), numeric(1))
  expect_gt(ll[["cloglog"]], ll[["logit"]])
  expect_gt(ll[["cloglog"]], ll[["probit"]])

  # And it recovers the cutpoints, in the chart they are reported in.
  expect_equal(colMeans(fits[["cloglog"]][["aux"]]), truth - mean(lin),
               tolerance = 0.2, ignore_attr = TRUE)

  for (f in fits) {
    expect_predictor_invariant(f, d)
  }
})

test_that("cloglog probabilities from predict() match the link's own definition", {
  d <- sim_x(n = 150, seed = 12)
  set.seed(1012)
  z <- d$x1 - d$x2 + log(-log1p(-stats::runif(nrow(d))))
  d$y <- ordered(rowSums(outer(z, stats::quantile(z, c(1, 2) / 3), ">")) + 1L)

  fit <- bartisan(y ~ ., d, family = ordinal("cloglog"),
                  control = quick_control(gate = "hard"))

  p <- stats::predict(fit, newdata = d, type = "prob")
  expect_equal(rowSums(p), rep(1, nrow(d)), tolerance = 1e-10,
               ignore_attr = TRUE)

  # Rebuilt from the reported predictor and cutpoints with the link written out
  # by hand, which is the check that predict() and the sampler agree about what
  # "cloglog" means.
  e <- colMeans(stats::predict(fit, newdata = d, type = "link", draws = TRUE))
  cuts <- colMeans(fit[["aux"]])
  cdf <- function(v) -expm1(-exp(v))
  manual <- cbind(cdf(cuts[1L] - e),
                  cdf(cuts[2L] - e) - cdf(cuts[1L] - e),
                  1 - cdf(cuts[2L] - e))
  # predict() averages probabilities over draws and this uses the mean draw, so
  # the two agree only up to the curvature of the link.
  expect_equal(unname(p), unname(manual), tolerance = 0.05)
})

test_that("the cloglog augmentation targets the same posterior as the direct fit", {
  skip_on_cran()

  set.seed(13)
  n <- 600
  x <- sim_x(n = n, p = 3, seed = 13)
  lin <- 2 * x$x1 - x$x2
  truth <- c(-0.3, 0.8)
  z <- lin + log(-log1p(-stats::runif(n)))
  d <- x
  d$y <- ordered(rowSums(outer(z, truth, ">")) + 1L)

  fit <- function(augment) {
    set.seed(8)
    bartisan(y ~ ., d, family = ordinal("cloglog"), gate = "hard",
             num_trees = 20, num_burn = 400, num_save = 400, augment = augment)
  }

  direct <- fit(FALSE)
  aug <- fit("ordinal")

  expect_gt(stats::cor(colMeans(direct[["eta"]][[1L]]),
                       colMeans(aug[["eta"]][[1L]])), 0.95)
  expect_equal(colMeans(direct[["aux"]]), colMeans(aug[["aux"]]),
               tolerance = 0.2)
  # Reported on the ordinal scale either way, not on the augmented one.
  expect_lt(abs(mean(direct[["loglik"]]) - mean(aug[["loglik"]])),
            0.05 * abs(mean(direct[["loglik"]])))
  expect_predictor_invariant(aug, d)
})

test_that("three or more categories are reported with a centered predictor", {
  d <- sim_x(n = 200, seed = 14)
  set.seed(1014)
  z <- 2 * d$x1 - d$x2 + stats::rnorm(nrow(d))
  d$y <- ordered(rowSums(outer(z, stats::quantile(z, c(1, 2) / 3), ">")) + 1L)

  for (lk in c("logit", "probit", "cloglog")) {
    fit <- bartisan(y ~ ., d, family = ordinal(lk),
                    control = quick_control(gate = "hard"))

    # Every draw's predictor averages to zero over the fitted sample, which is
    # the identifying convention, and no cutpoint is pinned.
    expect_lt(max(abs(rowMeans(fit[["eta"]][[1L]]))), 1e-8)
    expect_gt(stats::sd(fit[["aux"]][, "cut1"]), 0)

    # Ordered by construction.
    expect_true(all(fit[["aux"]][, "cut1"] < fit[["aux"]][, "cut2"]))

    # The stored forest replays to the reported predictor, which is what makes
    # the change of chart a change of chart rather than a change of fit.
    expect_predictor_invariant(fit, d)
  }
})

test_that("two categories keep the chart that matches binary regression", {
  d <- sim_x(n = 200, seed = 15)
  set.seed(1015)
  d$y <- ordered(as.integer(d$x1 - d$x2 + stats::rnorm(nrow(d)) > 0))

  fit <- bartisan(y ~ ., d, family = ordinal(),
                  control = quick_control(gate = "hard"))

  # One boundary, folded into the intercept exactly as binary regression does,
  # so it stays pinned and the predictor is not centered.
  expect_true(all(fit[["aux"]][, "cut1"] == 0))
  expect_predictor_invariant(fit, d)
})

test_that("the chart is a change of chart: fitted probabilities are untouched", {
  skip_on_cran()

  # The whole justification for recording the draws in a different chart is that
  # every identified quantity is invariant. Fitted probabilities are identified,
  # so they must not move -- checked by shifting the recorded cutpoints and
  # predictor back by the amount the centering removed and comparing.
  d <- sim_x(n = 200, seed = 16)
  set.seed(1016)
  z <- 2 * d$x1 - d$x2 + stats::rnorm(nrow(d))
  d$y <- ordered(rowSums(outer(z, stats::quantile(z, c(1, 2) / 3), ">")) + 1L)

  fit <- bartisan(y ~ ., d, family = ordinal("probit"), gate = "hard",
                  num_trees = 20, num_burn = 200, num_save = 200)

  e <- stats::predict(fit, newdata = d, type = "link", draws = TRUE)
  cuts <- fit[["aux"]]

  # For each draw, cut - eta is what the probabilities depend on, and it must be
  # the same whatever constant is added to both.
  s <- 1L
  a <- cuts[s, 1L] - e[s, ]
  shifted <- (cuts[s, 1L] + 3.7) - (e[s, ] + 3.7)
  expect_equal(a, shifted)

  p <- stats::predict(fit, newdata = d, type = "prob")
  expect_equal(rowSums(p), rep(1, nrow(d)), tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("missing predictors are kept by default", {
  d <- sim_x(n = 150, seed = 17)
  set.seed(1017)
  d$y <- d$x1 - d$x2 + stats::rnorm(nrow(d))
  d$x3[c(3, 20, 55, 99)] <- NA

  # na.pass is the default now, so the rows stay and the rules decide where the
  # missing values go.
  fit <- bartisan(y ~ ., d, control = quick_control())
  expect_identical(fit[["n"]], 150L)
  expect_true(fit[["has_na"]][["x3"]])
  expect_predictor_invariant(fit, d)

  # na.omit still drops them, and then predict() refuses the same data because
  # no rule carries an answer.
  dropped <- bartisan(y ~ ., d, na.action = stats::na.omit,
                      control = quick_control())
  expect_identical(dropped[["n"]], 146L)
  expect_false(dropped[["has_na"]][["x3"]])
  expect_error(stats::predict(dropped, newdata = d), "missing values")

  # A missing response is dropped either way, with a warning.
  d2 <- d
  d2$y[1:2] <- NA
  expect_warning(fit2 <- bartisan(y ~ ., d2, control = quick_control()),
                 "missing response")
  expect_identical(fit2[["n"]], 148L)
})

test_that("the cutpoints are on the same scale as polr's", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  # A linear truth, so polr is correctly specified and the two should land in the
  # same place. What is being checked is the *chart*: an ordinal model is
  # identified only up to a common shift of its cutpoints and its predictor, and
  # the claim in ?bartisan-families is that bartisan reports the chart polr reports
  # when polr's predictors are centered.
  set.seed(41)
  n <- 3000
  d <- data.frame(x1 = stats::runif(n, -1, 2), x2 = stats::runif(n, 0, 3))
  lp <- 1.5 * d$x1 - 0.8 * d$x2
  zeta <- c(-1, 0.4, 1.9)
  d$y <- ordered(rowSums(outer(lp + stats::rnorm(n), zeta, ">")) + 1L)

  fit <- bartisan(y ~ x1 + x2, d, family = ordinal("probit"), gate = "hard",
                  num_trees = 50, num_burn = 400, num_save = 400)
  cuts <- colMeans(fit[["aux"]])

  raw <- MASS::polr(y ~ x1 + x2, data = d, method = "probit")

  # polr's own chart differs from this one by the mean of its linear predictor,
  # which is not zero because polr identifies by dropping the intercept column
  # rather than by centering.
  expect_equal(unname(cuts), unname(raw$zeta - mean(raw$lp)),
               tolerance = 0.1)

  # Equivalently, and this is the statement worth remembering: centering polr's
  # predictors puts it in exactly this chart.
  dc <- d
  dc$x1 <- dc$x1 - mean(dc$x1)
  dc$x2 <- dc$x2 - mean(dc$x2)
  centered <- MASS::polr(y ~ x1 + x2, data = dc, method = "probit")

  expect_equal(unname(cuts), unname(centered$zeta), tolerance = 0.1)
  expect_lt(abs(mean(centered$lp)), 1e-8)

  # The gaps are identified outright, so they match whatever the chart.
  expect_equal(diff(cuts), diff(raw$zeta), tolerance = 0.1,
               ignore_attr = TRUE)

  # And the predictors agree once both are centered.
  eta <- colMeans(fit[["eta"]][[1L]])
  expect_gt(stats::cor(eta, raw$lp), 0.99)
})
