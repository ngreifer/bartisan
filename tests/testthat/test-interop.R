# Interfaces to the packages that assess a fit: the posterior predictive
# sampler, the pointwise likelihood, and the methods built on them.

# The sampler and the log density are two independent statements about the same
# distribution, and this is what checks them against each other: draw many
# replicates from one fixed set of parameter values, and compare the empirical
# distribution against the density the C++ engine reports at the same values.
# `iterations` is allowed to repeat, so `rep(1L, R)` is R independent draws from
# draw one.
replicates_at_one_draw <- function(fit, row, reps = 4000L, ...) {
  as.vector(rstantools::posterior_predict(fit, newdata = row,
                                          iterations = rep(1L, reps), ...))
}

density_at_one_draw <- function(fit, row, response, values) {
  nd <- row[rep(1L, length(values)), , drop = FALSE]
  nd[[response]] <- values

  as.vector(stats::predict(fit, newdata = nd, type = "density",
                           iterations = 1L, draws = TRUE))
}

test_that("the posterior predictive sampler agrees with the log density", {
  skip_if_not_installed("rstantools")
  skip_on_cran()

  d <- sim_x(200, 2)
  set.seed(1201)
  signal <- 2 * d$x1 - d$x2

  # A discrete family: the empirical frequency of every value on a grid wide
  # enough to hold essentially all the mass, against its probability.
  d$y <- rpois(nrow(d), exp(0.5 + signal))
  fit <- bartisan(y ~ x1 + x2, data = d, family = poisson(),
                  control = quick_control())

  grid <- 0:40
  probability <- density_at_one_draw(fit, d[1L, ], "y", grid)
  expect_equal(sum(probability), 1, tolerance = 1e-6)

  set.seed(1202)
  drawn <- replicates_at_one_draw(fit, d[1L, ])
  empirical <- vapply(grid, function(g) mean(drawn == g), numeric(1L))
  expect_lt(max(abs(empirical - probability)), 0.02)

  # A continuous family: the quantiles of the draws against the quantiles of the
  # numerically integrated density.
  set.seed(1203)
  d$y <- signal + rnorm(nrow(d))
  fit <- bartisan(y ~ x1 + x2, data = d, family = gaussian(),
                  control = quick_control())

  grid <- seq(-8, 10, length.out = 1001L)
  density <- density_at_one_draw(fit, d[1L, ], "y", grid)
  cdf <- cumsum(density) * diff(grid)[1L]
  expect_equal(max(cdf), 1, tolerance = 0.01)

  set.seed(1204)
  drawn <- replicates_at_one_draw(fit, d[1L, ])
  probs <- c(0.1, 0.5, 0.9)
  target <- suppressWarnings(stats::approx(cdf, grid, xout = probs)$y)
  expect_lt(max(abs(stats::quantile(drawn, probs, names = FALSE) - target)),
            0.25)

  # And the mixture, whose sampler and density have to agree about the reporting
  # chart: the components are stored centred, the baseline a fresh component
  # comes from is not, and getting either shift wrong moves the draws off the
  # density by the amount the mixture was off centre. A skewed error makes that
  # amount large enough to see.
  set.seed(1205)
  d$y <- signal + 2 * (rgamma(nrow(d), 1.5, 1.5) - 1)
  fit <- bartisan(y ~ x1 + x2, data = d, family = dpm(),
                  control = quick_control(num_burn = 300L, num_draws = 300L))

  expect_gt(abs(mean(fit[["aux"]][, "center"])), 0.05)

  grid <- seq(-12, 16, length.out = 1401L)
  density <- density_at_one_draw(fit, d[1L, ], "y", grid)
  cdf <- cumsum(density) * diff(grid)[1L]
  expect_equal(max(cdf), 1, tolerance = 0.02)

  set.seed(1206)
  drawn <- replicates_at_one_draw(fit, d[1L, ])
  target <- suppressWarnings(stats::approx(cdf, grid, xout = probs)$y)
  expect_lt(max(abs(stats::quantile(drawn, probs, names = FALSE) - target)),
            0.5)
})

test_that("every category of a categorical family is drawn at its probability", {
  skip_if_not_installed("rstantools")
  skip_on_cran()

  d <- sim_x(200, 2)
  set.seed(1211)
  cut_points <- c(-0.5, 0.8, 1.6)
  d$y <- factor(findInterval(2 * d$x1 - d$x2 + rlogis(nrow(d)), cut_points) + 1L,
                labels = c("a", "b", "c", "d"), ordered = TRUE)

  for (link in c("logit", "probit", "cloglog")) {
    fit <- bartisan(y ~ x1 + x2, data = d, family = ordinal(link),
                    control = quick_control())

    nd <- d[rep(1L, 4L), ]
    nd$y <- factor(c("a", "b", "c", "d"), levels = c("a", "b", "c", "d"),
                   ordered = TRUE)
    probability <- as.vector(stats::predict(fit, newdata = nd,
                                            type = "density", iterations = 1L,
                                            draws = TRUE))
    expect_equal(sum(probability), 1, tolerance = 1e-6)

    set.seed(1212)
    drawn <- replicates_at_one_draw(fit, d[1L, ])

    # The codes index the categories, so they are within range by construction
    # and cover them all.
    expect_true(all(drawn %in% seq_along(fit[["levels"]])))

    empirical <- vapply(seq_along(fit[["levels"]]),
                        function(k) mean(drawn == k), numeric(1L))
    expect_lt(max(abs(empirical - probability)), 0.03)
  }
})

test_that("the zero-inflated and ordered beta samplers hit their point masses", {
  skip_if_not_installed("rstantools")
  skip_on_cran()

  d <- sim_x(200, 2)
  set.seed(1221)
  signal <- 2 * d$x1 - d$x2

  d$y <- ifelse(runif(nrow(d)) < 0.3, 0, rpois(nrow(d), exp(0.7 + signal)))
  fit <- bartisan(y ~ x1 + x2, data = d, family = zi_poisson(),
                  control = quick_control())

  grid <- 0:40
  probability <- density_at_one_draw(fit, d[1L, ], "y", grid)
  set.seed(1222)
  drawn <- replicates_at_one_draw(fit, d[1L, ])
  empirical <- vapply(grid, function(g) mean(drawn == g), numeric(1L))
  expect_lt(max(abs(empirical - probability)), 0.03)

  # An ordered beta response is two point masses and a density between them, and
  # the three have to partition the draws.
  mu <- stats::plogis(signal)
  y <- stats::rbeta(nrow(d), mu * 6, 6 - mu * 6)
  y[runif(nrow(d)) < 0.15] <- 0
  y[runif(nrow(d)) < 0.15] <- 1
  d$y <- y

  fit <- bartisan(y ~ x1 + x2, data = d, family = ordbeta(),
                  control = quick_control())

  ends <- density_at_one_draw(fit, d[1L, ], "y", c(0, 1))
  set.seed(1223)
  drawn <- replicates_at_one_draw(fit, d[1L, ])

  expect_equal(mean(drawn == 0), ends[1L], tolerance = 0.05)
  expect_equal(mean(drawn == 1), ends[2L], tolerance = 0.05)
  expect_true(all(drawn >= 0 & drawn <= 1))
})

# The families the tests above leave out, checked the same way: replicates drawn
# at one posterior draw against the distribution the log density gives at that
# draw. A count is compared value by value on a grid that holds essentially all
# its mass, and a continuous response by the largest gap between the empirical
# distribution function of the replicates and the integrated density. Over five
# seeds the largest gap was 0.021, with 4000 replicates.
test_that("each family's replicate draws follow its log density", {
  skip_if_not_installed("rstantools")
  skip_on_cran()
  skip_if_not_installed("patrick")

  patrick::with_parameters_test_that(
    "family:",
    {
      d <- sim_x(200, 2, seed = seed)
      d$y <- draw(2 * d$x1 - d$x2)

      fit <- bartisan(y ~ x1 + x2, data = d, family = family,
                      control = quick_control())

      probability <- density_at_one_draw(fit, d[1L, ], "y", grid)
      set.seed(seed + 1L)
      drawn <- replicates_at_one_draw(fit, d[1L, ])

      if (discrete) {
        expect_equal(sum(probability), 1, tolerance = 1e-6)
        empirical <- vapply(grid, function(g) mean(drawn == g), numeric(1L))
        expect_lt(max(abs(empirical - probability)), 0.03)
      }
      else {
        cdf <- c(0, cumsum(diff(grid) * (probability[-1L] +
                                           probability[-length(grid)]) / 2))
        expect_equal(max(cdf), 1, tolerance = 0.01)
        expect_lt(max(abs(stats::ecdf(drawn)(grid) - cdf)), 0.035)
      }
    },
    patrick::cases(
      Beta = list(
        family = Beta(), seed = 1281L, discrete = FALSE,
        draw = function(s) {
          mu <- stats::plogis(s)
          stats::rbeta(length(s), 20 * mu, 20 * (1 - mu))
        },
        grid = seq(0, 1, length.out = 4001L)[-c(1L, 4001L)]),
      negbin = list(
        family = negbin(), seed = 1282L, discrete = TRUE,
        draw = function(s) stats::rnbinom(length(s), size = 3, mu = exp(0.5 + s)),
        grid = 0:150),
      zi_negbin = list(
        family = zi_negbin(), seed = 1283L, discrete = TRUE,
        draw = function(s) {
          ifelse(stats::runif(length(s)) < 0.3, 0,
                 stats::rnbinom(length(s), size = 3, mu = exp(0.7 + s)))
        },
        grid = 0:150),
      Gamma = list(
        family = stats::Gamma("log"), seed = 1284L, discrete = FALSE,
        draw = function(s) stats::rgamma(length(s), shape = 3, rate = 3 / exp(s)),
        grid = seq(0, 60, length.out = 6001L)[-1L]),
      Gamma_ls = list(
        family = Gamma_ls(), seed = 1285L, discrete = FALSE,
        draw = function(s) stats::rgamma(length(s), shape = 3, rate = 3 / exp(s)),
        grid = seq(0, 60, length.out = 6001L)[-1L]),
      gaussian_ls = list(
        family = gaussian_ls(), seed = 1286L, discrete = FALSE,
        draw = function(s) stats::rnorm(length(s), s, exp(-0.5 + 0.5 * s)),
        grid = seq(-15, 15, length.out = 6001L))
    )
  )
})

test_that("a binomial replicate is a fraction of the trials", {
  skip_if_not_installed("rstantools")

  d <- sim_x(100, 2)
  set.seed(1234)
  d$y <- rbinom(nrow(d), 6L, stats::plogis(2 * d$x1 - d$x2)) / 6
  d$trials <- 6

  fit <- bartisan(y ~ x1 + x2, data = d, family = binomial(), weights = trials,
                  control = quick_control())

  # The number of trials is not a function of the predictors, so predicting for
  # new rows has to be told what it is rather than assuming one trial.
  expect_error(rstantools::posterior_predict(fit, newdata = d[1L, ]),
               "number of trials")

  set.seed(1235)
  drawn <- replicates_at_one_draw(fit, d[1L, ], weights = 6)
  expect_true(all(drawn %in% ((0:6) / 6)))

  # A fraction of a trial cannot be drawn.
  expect_error(rstantools::posterior_predict(fit, newdata = d[1L, ],
                                             weights = 2.5),
               "whole numbers of trials")

  # Binary data are the same statement with one trial, so they come back as
  # zeros and ones rather than as counts.
  set.seed(1236)
  d$y <- rbinom(nrow(d), 1L, stats::plogis(2 * d$x1 - d$x2))
  fit <- bartisan(y ~ x1 + x2, data = d, family = binomial(),
                  control = quick_control())
  expect_setequal(unique(as.vector(rstantools::posterior_predict(fit))), c(0, 1))
})

test_that("a binomial replicate is drawn at its probability", {
  skip_if_not_installed("rstantools")
  skip_on_cran()

  d <- sim_x(200, 2)
  set.seed(1231)
  d$y <- rbinom(nrow(d), 6L, stats::plogis(2 * d$x1 - d$x2)) / 6
  d$trials <- 6

  fit <- bartisan(y ~ x1 + x2, data = d, family = binomial(), weights = trials,
                  control = quick_control())

  set.seed(1232)
  drawn <- replicates_at_one_draw(fit, d[1L, ], weights = 6)

  nd <- d[rep(1L, 7L), ]
  nd$y <- (0:6) / 6
  probability <- as.vector(
    stats::predict(fit, newdata = nd, type = "density", iterations = 1L,
                   draws = TRUE, weights = rep(6, 7L)))
  empirical <- vapply((0:6) / 6, function(g) mean(abs(drawn - g) < 1e-8),
                      numeric(1L))
  expect_lt(max(abs(empirical - probability)), 0.03)
})

test_that("an accelerated failure time replicate is an event time", {
  skip_if_not_installed("rstantools")
  skip_if_not_installed("survival")

  d <- sim_x(200, 2)
  set.seed(1241)
  d$time <- exp(1 + 2 * d$x1 + 0.4 * rnorm(nrow(d)))
  d$event <- rbinom(nrow(d), 1L, 0.8)

  for (family in list(weibull_aft(), lognormal_aft(), loglogistic_aft())) {
    fit <- bartisan(survival::Surv(time, event) ~ x1 + x2, data = d,
                    family = family, control = quick_control())

    set.seed(1242)
    drawn <- replicates_at_one_draw(fit, d[1L, ])
    expect_true(all(drawn > 0))

    # The residual is on the log time scale the model is linear on.
    expect_equal(stats::residuals(fit),
                 log(d$time) - colMeans(fit[["eta"]][[1L]]))

    # The density is on the log time scale, so the comparison is too.
    grid <- seq(-4, 8, length.out = 1201L)
    nd <- d[rep(1L, length(grid)), ]
    nd$time <- exp(grid)
    nd$event <- 1L
    density <- as.vector(stats::predict(fit, newdata = nd, type = "density",
                                        iterations = 1L, draws = TRUE))
    cdf <- cumsum(density) * diff(grid)[1L]
    expect_equal(max(cdf), 1, tolerance = 0.01)

    target <- suppressWarnings(stats::approx(cdf, grid, xout = 0.5)$y)
    expect_equal(stats::median(log(drawn)), target, tolerance = 0.1)
  }
})

test_that("the predictive draws have the mean the response scale reports", {
  skip_if_not_installed("rstantools")

  d <- sim_x(150, 2)
  set.seed(1251)
  d$y <- rpois(nrow(d), exp(0.5 + 2 * d$x1))

  fit <- bartisan(y ~ x1 + x2, data = d, family = poisson(),
                  control = quick_control(num_draws = 400L))

  set.seed(1252)
  drawn <- rstantools::posterior_predict(fit)

  # Averaged over draws and observations the two have to agree; per observation
  # they would not, since one replicate per draw is a noisy estimate of a mean.
  expect_equal(mean(drawn), mean(stats::fitted(fit)), tolerance = 0.05)

  expect_identical(dim(rstantools::posterior_epred(fit)),
                   dim(fit[["eta"]][[1L]]))
  expect_equal(colMeans(rstantools::posterior_epred(fit)),
               stats::fitted(fit))
  expect_equal(rstantools::posterior_linpred(fit), fit[["eta"]][[1L]])
  expect_equal(rstantools::posterior_linpred(fit, transform = TRUE),
               rstantools::posterior_epred(fit))
})

test_that("newdata and iterations are respected", {
  skip_if_not_installed("rstantools")

  d <- sim_x(120, 2)
  set.seed(1261)
  d$y <- 2 * d$x1 + rnorm(nrow(d))

  fit <- bartisan(y ~ x1 + x2, data = d, control = quick_control())

  nd <- d[1:5, ]
  expect_identical(dim(rstantools::posterior_predict(fit, newdata = nd)),
                   c(30L, 5L))
  expect_identical(dim(rstantools::posterior_predict(fit, iterations = 1:4)),
                   c(4L, 120L))
  expect_identical(dim(rstantools::log_lik(fit, newdata = nd)), c(30L, 5L))
})

test_that("log_lik is the pointwise log density", {
  skip_if_not_installed("rstantools")

  d <- sim_x(120, 2)
  set.seed(1271)
  d$y <- rpois(nrow(d), exp(0.5 + 2 * d$x1))

  fit <- bartisan(y ~ x1 + x2, data = d, family = poisson(),
                  control = quick_control())

  expect_equal(rstantools::log_lik(fit),
               stats::predict(fit, type = "density", draws = TRUE, log = TRUE))

  # And it has to reconcile with the total the sampler recorded, which is the
  # weighted sum over observations.
  expect_equal(rowSums(rstantools::log_lik(fit)), fit[["loglik"]],
               tolerance = 1e-6)
})

test_that("simulate follows the stats contract", {
  skip_if_not_installed("rstantools")

  d <- sim_x(120, 2)
  set.seed(1281)
  d$y <- 2 * d$x1 + rnorm(nrow(d))

  fit <- bartisan(y ~ x1 + x2, data = d, control = quick_control())

  out <- stats::simulate(fit, nsim = 3L, seed = 42L)
  expect_s3_class(out, "data.frame")
  expect_identical(dim(out), c(120L, 3L))
  expect_named(out, c("sim_1", "sim_2", "sim_3"))
  expect_identical(attr(out, "seed"), 42L)

  # Same seed, same draws; and the caller's stream is left where it was found.
  set.seed(7)
  before <- .Random.seed
  again <- stats::simulate(fit, nsim = 3L, seed = 42L)
  expect_equal(out, again, ignore_attr = TRUE)
  expect_identical(.Random.seed, before)

  # With no stream at all yet, one is started so that there is one to restore,
  # and the seed still decides the draws.
  saved <- .Random.seed
  on.exit(assign(".Random.seed", saved, envir = globalenv()), add = TRUE)
  rm(".Random.seed", envir = globalenv())

  fresh <- stats::simulate(fit, nsim = 3L, seed = 42L)
  expect_equal(fresh, out, ignore_attr = TRUE)
  expect_true(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("simulate returns a factor when the response was one", {
  skip_if_not_installed("rstantools")

  d <- sim_x(150, 2)
  set.seed(1291)
  d$y <- factor(ifelse(rbinom(nrow(d), 1L, stats::plogis(2 * d$x1)) == 1L,
                       "yes", "no"))

  fit <- bartisan(y ~ x1 + x2, data = d, family = binomial(),
                  control = quick_control())
  out <- stats::simulate(fit, nsim = 2L)
  expect_s3_class(out$sim_1, "factor")
  expect_identical(levels(out$sim_1), levels(d$y))

  set.seed(1292)
  d$y <- factor(findInterval(2 * d$x1 + rlogis(nrow(d)), c(0, 1)) + 1L,
                labels = c("lo", "mid", "hi"), ordered = TRUE)
  fit <- bartisan(y ~ x1 + x2, data = d, family = ordinal(),
                  control = quick_control())
  out <- stats::simulate(fit, nsim = 2L)
  expect_s3_class(out$sim_1, "ordered")
  expect_identical(levels(out$sim_1), levels(d$y))
})

test_that("the accessors report what a glm's would", {
  d <- sim_x(120, 2)
  set.seed(1301)
  d$w <- rep(c(1, 2), length.out = nrow(d))
  d$y <- 2 * d$x1 + rnorm(nrow(d))

  # Named rather than inferred: a numeric response defaults to `dpm()`, which
  # refuses weights.
  fit <- bartisan(y ~ x1 + x2, data = d, family = gaussian(), weights = w,
                  control = quick_control())

  expect_equal(stats::fitted(fit), stats::predict(fit))
  expect_equal(stats::fitted(fit, type = "link"),
               stats::predict(fit, type = "link"))
  expect_equal(stats::residuals(fit), d$y - stats::fitted(fit))
  expect_equal(stats::weights(fit), d$w)
  expect_equal(stats::sigma(fit), mean(fit[["aux"]][, "sigma"]))

  # A family whose scale is not a single number has no sigma to report, and
  # inventing one would be worse than saying so.
  set.seed(1302)
  d$y <- rpois(nrow(d), exp(0.5 + 2 * d$x1))
  poisson_fit <- bartisan(y ~ x1 + x2, data = d, family = poisson(),
                          control = quick_control())
  expect_null(stats::sigma(poisson_fit))
})

test_that("a family with no mean and no sampler says so", {
  skip_if_not_installed("rstantools")

  d <- sim_x(120, 2)
  set.seed(1311)
  d$y <- factor(findInterval(2 * d$x1 + rlogis(nrow(d)), c(0, 1)) + 1L,
                labels = c("lo", "mid", "hi"), ordered = TRUE)

  fit <- bartisan(y ~ x1 + x2, data = d, family = ordinal(),
                  control = quick_control())
  expect_error(stats::residuals(fit), "no mean")

  set.seed(1312)
  d$y <- 2 * d$x1 + rnorm(nrow(d))
  custom <- custom_family(
    function(y, eta, ...) stats::dnorm(y, eta[, 1L], 1, log = TRUE),
    num_predictors = 1L)
  custom_fit <- bartisan(y ~ x1 + x2, data = d, family = custom,
                         control = quick_control())

  expect_error(rstantools::posterior_predict(custom_fit), "no way to draw")
  expect_error(stats::residuals(custom_fit), "no mean")

  # The log density is the one thing a custom family does supply, so the
  # pointwise likelihood still works and so does everything built on it.
  expect_identical(dim(rstantools::log_lik(custom_fit)), c(30L, 120L))

  # A proportional hazards fit has a survival curve and no sampler for event
  # times.
  d$time <- stats::rexp(nrow(d), exp(-d$x1))
  d$status <- stats::rbinom(nrow(d), 1L, 0.8)
  ph_fit <- bartisan(cbind(time, status) ~ x1 + x2, data = d, family = ph(),
                     control = quick_control())

  expect_error(rstantools::posterior_predict(ph_fit),
               "no posterior predictive sampler")
})

test_that("loo and waic run on the pointwise likelihood", {
  skip_if_not_installed("loo")
  skip_if_not_installed("rstantools")

  d <- sim_x(150, 2)
  set.seed(1321)
  d$y <- rpois(nrow(d), exp(0.5 + 2 * d$x1))

  fit <- bartisan(y ~ x1 + x2, data = d, family = poisson(), chains = 2L,
                  control = quick_control(num_draws = 100L))

  # The chain structure is what loo needs to judge the efficiency of the draws.
  expect_identical(chain_ids(fit), rep(1:2, each = 100L))

  result <- suppressWarnings(loo::loo(fit))
  expect_s3_class(result, "loo")
  expect_identical(dim(result[["pointwise"]])[1L], 150L)

  # The leave-one-out estimate is a penalized version of the in-sample one, so
  # it has to be the smaller of the two.
  lppd <- sum(stats::predict(fit, type = "density", log = TRUE))
  expect_lt(result[["estimates"]]["elpd_loo", "Estimate"], lppd)

  waic_result <- suppressWarnings(loo::waic(fit))
  expect_s3_class(waic_result, "waic")
  expect_equal(waic_result[["estimates"]]["elpd_waic", "Estimate"],
               result[["estimates"]]["elpd_loo", "Estimate"],
               tolerance = 0.05)
})

test_that("the Bayesian R-squared is a posterior of a ratio", {
  skip_if_not_installed("performance")
  skip_if_not_installed("rstantools")

  d <- sim_x(150, 2)
  set.seed(1333)
  d$y <- 3 * d$x1 - 2 * d$x2 + rnorm(nrow(d), sd = 0.5)

  fit <- bartisan(y ~ x1 + x2, data = d, control = quick_control())

  draws <- performance::r2_posterior(fit)[["R2_Bayes"]]
  expect_length(draws, 30L)
  expect_true(all(draws > 0 & draws < 1))

  summary_r2 <- performance::r2(fit)
  expect_equal(unname(summary_r2[["R2_Bayes"]]), stats::median(draws),
               tolerance = 1e-8)

  # A categorical response has no mean, so there is no such ratio to report.
  set.seed(1334)
  d$y <- factor(findInterval(3 * d$x1 + rlogis(nrow(d)), c(0, 1)) + 1L,
                labels = c("lo", "mid", "hi"), ordered = TRUE)
  ordinal_fit <- bartisan(y ~ x1 + x2, data = d, family = ordinal(),
                          control = quick_control())
  expect_warning(out <- performance::r2_posterior(ordinal_fit), "no mean")
  expect_null(out)
})

test_that("the Bayesian R-squared is well away from zero for a strong signal", {
  skip_if_not_installed("performance")
  skip_if_not_installed("rstantools")
  skip_on_cran()

  d <- sim_x(200, 2)
  set.seed(1331)
  d$y <- 3 * d$x1 - 2 * d$x2 + rnorm(nrow(d), sd = 0.5)

  fit <- bartisan(y ~ x1 + x2, data = d, control = quick_control(num_draws = 60L))

  draws <- performance::r2_posterior(fit)[["R2_Bayes"]]

  # A strong signal, so it should be well away from zero.
  expect_gt(mean(draws), 0.5)
})

test_that("model_performance collects the fit statistics", {
  skip_if_not_installed("performance")
  skip_if_not_installed("loo")

  d <- sim_x(150, 2)
  set.seed(1341)
  d$y <- 2 * d$x1 + rnorm(nrow(d))

  fit <- bartisan(y ~ x1 + x2, data = d, control = quick_control())

  out <- suppressWarnings(performance::model_performance(fit))
  expect_s3_class(out, "performance_model")
  expect_identical(nrow(out), 1L)
  expect_true(all(c("ELPD", "LOOIC", "WAIC", "R2", "RMSE", "Sigma") %in%
                    names(out)))
  expect_equal(out[["RMSE"]], sqrt(mean(stats::residuals(fit)^2)))
  expect_equal(out[["Sigma"]], stats::sigma(fit))

  # A subset is honored, and the metrics that were not asked for are absent.
  subset <- suppressWarnings(
    performance::model_performance(fit, metrics = c("R2", "RMSE")))
  expect_named(subset, c("R2", "RMSE"))
})

test_that("model_performance leaves out the residual metrics a response has none for", {
  skip_if_not_installed("performance")
  skip_if_not_installed("loo")

  d <- sim_x(80, 2, seed = 1342L)
  d$y <- factor(sample(1:3, nrow(d), replace = TRUE), ordered = TRUE)

  fit <- bartisan(y ~ x1 + x2, data = d, family = ordinal(),
                  control = quick_control(num_burn = 10L, num_draws = 10L))

  # Categories have no mean, so no residual, and neither metric built on one.
  out <- suppressWarnings(performance::model_performance(fit))
  expect_true(all(c("ELPD", "WAIC") %in% names(out)))
  expect_false(any(c("RMSE", "Sigma") %in% names(out)))
})

test_that("as_draws hands the scalar parameters over with their chain structure", {
  skip_if_not_installed("posterior")

  d <- sim_x(120, 2)
  set.seed(1351)
  d$y <- 2 * d$x1 + rnorm(nrow(d))

  # Named, because `aux.sigma` below is the Gaussian family's nuisance parameter
  # and the numeric default is now the mixture.
  fit <- bartisan(y ~ x1 + x2, data = d, family = gaussian(), chains = 3L,
                  control = quick_control(num_draws = 40L))

  draws <- posterior::as_draws(fit)
  expect_s3_class(draws, "draws_array")
  expect_identical(posterior::niterations(draws), 40L)
  expect_identical(posterior::nchains(draws), 3L)
  expect_true(all(c("loglik", "aux.sigma") %in% posterior::variables(draws)))

  # The draws have to be the same numbers in the same order the fit stores them,
  # chain by chain.
  expect_equal(as.vector(draws[, 1L, "loglik"]), fit[["loglik"]][1:40])
  expect_equal(as.vector(draws[, 3L, "loglik"]), fit[["loglik"]][81:120])

  # And the summary of them has to agree with the fit's own diagnostics, which
  # are computed from the same vectors by a different implementation of the same
  # estimators. Exactly, not approximately: at a tolerance of 1% this hid a
  # rank normalization that stretched the top of its scale and an effective
  # sample size that counted its last autocorrelations twice. The equality holds
  # for an even number of draws per chain; see `rank_normalize()` for odd ones.
  summary_draws <- posterior::summarise_draws(draws)
  own <- diagnose(fit)[["table"]]
  theirs <- summary_draws[summary_draws$variable == "loglik", ]
  ours <- own[own$quantity == "loglik", ]
  expect_equal(theirs$rhat, ours$rhat, tolerance = 1e-10)
  expect_equal(theirs$ess_bulk, ours$ess_bulk, tolerance = 1e-10)
  expect_equal(theirs$ess_tail, ours$ess_tail, tolerance = 1e-10)
})

test_that("pp_check runs a bayesplot check", {
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("rstantools")

  d <- sim_x(150, 2)
  set.seed(1361)
  d$y <- rpois(nrow(d), exp(0.5 + 2 * d$x1))

  fit <- bartisan(y ~ x1 + x2, data = d, family = poisson(),
                  control = quick_control())

  expect_s3_class(bayesplot::pp_check(fit, ndraws = 5L), "ggplot")
  expect_s3_class(bayesplot::pp_check(fit, "hist", ndraws = 3L), "ggplot")
  expect_error(bayesplot::pp_check(fit, "not_a_check"), "not a")
})

# A replicate is compared with the response on the scale it is drawn on: event
# times rather than the log times an accelerated failure time model stores, and
# category codes from one rather than the zero-based ones the sampler indexes.
test_that("pp_check compares replicates with the response on their own scale", {
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("survival")

  d <- sim_x(n = 80L, p = 2L, seed = 1362L)
  d$time <- stats::rexp(nrow(d), exp(-d$x1))
  d$status <- stats::rbinom(nrow(d), 1L, 0.7)
  d$grade <- factor(sample(1:3, nrow(d), replace = TRUE), ordered = TRUE)

  ctrl <- quick_control(num_burn = 10L, num_draws = 10L)

  aft <- bartisan(survival::Surv(time, status) ~ x1 + x2, data = d,
                  family = weibull_aft(), control = ctrl)

  expect_equal(observed_response(aft), d$time)

  # Censored times make any comparison of event times worth a warning.
  expect_warning(dens <- bayesplot::pp_check(aft, ndraws = 3L), "censored")
  expect_s3_class(dens, "ggplot")

  # The Kaplan-Meier check takes its events from the survival object. Drawing it
  # needs *ggfortify*, which *bayesplot* suggests rather than imports.
  expect_identical(survival_status(aft), as.numeric(d$status))

  if (rlang::is_installed("ggfortify")) {
    km <- suppressWarnings(bayesplot::pp_check(aft, "km_overlay", ndraws = 3L))
    expect_s3_class(km, "ggplot")
  }

  # A response given as a matrix carries no event indicator to take.
  matrix_fit <- bartisan(cbind(time, status) ~ x1 + x2, data = d,
                         family = weibull_aft(), control = ctrl)
  expect_error(suppressWarnings(bayesplot::pp_check(matrix_fit, "km_overlay")),
               "status_y")

  ordered_fit <- bartisan(grade ~ x1 + x2, data = d, family = ordinal(),
                          control = ctrl)

  expect_identical(observed_response(ordered_fit), as.integer(d$grade))
  expect_s3_class(bayesplot::pp_check(ordered_fit, "bars", ndraws = 3L),
                  "ggplot")
})

test_that("the easystats packages read the fit through the accessors", {
  skip_if_not_installed("insight")
  skip_if_not_installed("performance")

  d <- sim_x(200, 2)
  set.seed(1371)
  d$y <- 2 * d$x1 - d$x2 + rnorm(nrow(d))

  fit <- bartisan(y ~ x1 + x2, data = d, control = quick_control())

  expect_equal(as.numeric(insight::get_sigma(fit)), stats::sigma(fit))
  expect_equal(insight::get_residuals(fit), stats::residuals(fit),
               ignore_attr = TRUE)
  expect_equal(performance::performance_mse(fit),
               mean(stats::residuals(fit)^2))

  # The posterior predictive check runs through simulate(), so this is the whole
  # chain from the sampler to an easystats plot.
  checked <- performance::check_predictions(fit, iterations = 10L)
  expect_s3_class(checked, "performance_pp_check")
  expect_identical(ncol(checked), 11L)
})

test_that("as_draws() carries the additive predictor, which is what mixing is about", {
  skip_if_not_installed("posterior")

  d <- sim_x(n = 120, p = 3, seed = 781)
  set.seed(7811)
  d$y <- d$x1 + stats::rnorm(nrow(d), sd = 0.3)

  fit <- bartisan(y ~ ., d, family = stats::gaussian(), chains = 2L,
                  control = quick_control(num_draws = 60L))

  vars <- posterior::variables(posterior::as_draws(fit))

  # A representative spread of eta columns, not all of them: a draws array with
  # one column per observation is not something summarise_draws() can be used on.
  eta_vars <- grep("^eta\\[", vars, value = TRUE)
  expect_gt(length(eta_vars), 0L)
  expect_lte(length(eta_vars), 10L)
  expect_true(all(c("loglik", "sigma_mu.eta") %in% vars))

  # Opting out restores the scalars-only object.
  expect_false(any(grepl("^eta\\[", posterior::variables(
    posterior::as_draws(fit, eta = FALSE)))))

  # Naming observations takes exactly those.
  expect_true(all(c("eta[1]", "eta[2]") %in%
                    posterior::variables(posterior::as_draws(fit, eta = c(1, 2)))))

  expect_error(posterior::as_draws(fit, eta = c(1, 1e6)),
               "observation indices")

  # The values are the draws themselves, and summarise_draws() can read them.
  drawn <- posterior::as_draws(fit, eta = 3)
  expect_equal(as.vector(drawn[, , "eta[3]"]), unname(fit$eta$eta[, 3]),
               ignore_attr = TRUE)
  expect_true("eta[3]" %in% posterior::summarise_draws(drawn)$variable)
})

test_that("as_draws() names the predictor for each forest and keeps a small fit whole", {
  skip_if_not_installed("posterior")

  d <- sim_x(n = 40L, p = 2L, seed = 782L)
  d$y <- stats::rnorm(nrow(d), d$x1)
  ctrl <- quick_control(num_burn = 10L, num_draws = 10L)

  eta_names <- function(fit) {
    grep("^eta", posterior::variables(posterior::as_draws(fit)), value = TRUE)
  }

  # Ten or fewer observations are all taken, in order, rather than a spread.
  small <- bartisan(y ~ x1, d[1:8, ], family = stats::gaussian(),
                    control = ctrl)
  expect_identical(eta_names(small), sprintf("eta[%d]", 1:8))

  # With two predictors, each column says which one it belongs to.
  two <- bartisan(y ~ x1, d, family = gaussian_ls(), control = ctrl)
  names_two <- eta_names(two)
  expect_true(all(grepl("^eta[.](mean|log_sd)\\[[0-9]+\\]$", names_two)))
  expect_true(any(startsWith(names_two, "eta.mean[")))
  expect_true(any(startsWith(names_two, "eta.log_sd[")))
})

# The `ppc_loo_*` checks reweight the replicates towards the leave-one-out
# predictive, so they need the importance weights as well as the draws.
# `pp_check()` passed only `y` and `yrep`, and every one of them failed with
# "One of 'lw' and 'psis_object' must be specified."
test_that("pp_check() supplies the weights a leave-one-out check needs", {
  skip_on_cran()
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("loo")
  skip_if_not_installed("rstantools")

  d <- sim_x(n = 150L, p = 3L, seed = 71L)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d), 0, 0.5)

  fit <- bartisan(y ~ ., d, family = stats::gaussian(),
                  control = quick_control(num_trees = 5L, num_burn = 60L,
                                          num_draws = 100L))

  # `loo_calibration` is left out: it wants a binary response, which is
  # *bayesplot*'s own requirement and is checked below.
  drew <- vapply(c("loo_pit_ecdf", "loo_pit_overlay", "loo_pit_qq",
                   "loo_intervals", "loo_ribbon"),
                 function(type) {
                   p <- suppressMessages(suppressWarnings(
                     bayesplot::pp_check(fit, type = type)))
                   inherits(p, "ggplot")
                 }, logical(1L))

  expect_true(all(drew))
  expect_named(drew, c("loo_pit_ecdf", "loo_pit_overlay", "loo_pit_qq",
                       "loo_intervals", "loo_ribbon"))

  # Every draw is used, whatever `ndraws` says, because the weights and the
  # replicates have to be the same shape. Saying so beats ignoring it.
  #
  # Collected rather than matched with `expect_warning()`, because the Pareto
  # diagnostic may warn here too and which warnings fire is not the point.
  warns <- character()

  suppressMessages(withCallingHandlers(
    bayesplot::pp_check(fit, type = "loo_pit_ecdf", ndraws = 5L),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }))

  expect_true(any(grepl("does not apply", warns, fixed = TRUE)))

  # A check that does not reweight still subsamples, and still says nothing.
  expect_no_warning(
    suppressMessages(bayesplot::pp_check(fit, type = "dens_overlay",
                                         ndraws = 5L)))

  # Weights the caller computed themselves are used instead of ours.
  ll <- rstantools::log_lik(fit)
  psis <- suppressWarnings(loo::psis(-ll, r_eff = NA))

  expect_s3_class(
    suppressMessages(suppressWarnings(
      bayesplot::pp_check(fit, type = "loo_pit_ecdf", psis_object = psis))),
    "ggplot")
})

test_that("a binary response reaches the calibration check", {
  skip_on_cran()
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("loo")

  d <- sim_x(n = 200L, p = 3L, seed = 72L)
  d$y <- stats::rbinom(nrow(d), 1L, stats::plogis(-0.5 + 2 * d$x1))

  fit <- bartisan(y ~ ., d, family = stats::binomial(),
                  control = quick_control(num_trees = 5L, num_burn = 60L,
                                          num_draws = 100L))

  expect_s3_class(
    suppressMessages(suppressWarnings(
      bayesplot::pp_check(fit, type = "loo_calibration"))),
    "ggplot")

  expect_s3_class(
    suppressMessages(suppressWarnings(
      bayesplot::pp_check(fit, type = "calibration"))),
    "ggplot")
})

# A binned residual plot and a calibration plot bin their second argument and
# read the outcome within each bin, so replicate outcomes give two degenerate
# bins: every binned mean came back as exactly 0 or exactly 1 and the plot said
# nothing. They get the predictive mean instead.
test_that("the binned and calibration checks are passed probabilities", {
  skip_on_cran()
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("rstantools")

  d <- sim_x(n = 200L, p = 3L, seed = 73L)
  d$y <- stats::rbinom(nrow(d), 1L, stats::plogis(-0.5 + 2 * d$x1))

  fit <- bartisan(y ~ ., d, family = stats::binomial(),
                  control = quick_control(num_trees = 5L, num_burn = 60L,
                                          num_draws = 100L))

  binned <- suppressMessages(suppressWarnings(
    bayesplot::pp_check(fit, type = "error_binned", ndraws = 2L)))

  expect_s3_class(binned, "ggplot")

  # The binned means are strictly inside the unit interval, which is what
  # distinguishes probabilities from the zeros and ones that used to arrive.
  expect_true(all(binned[["data"]][["ey_bar"]] > 0))
  expect_true(all(binned[["data"]][["ey_bar"]] < 1))

  # `ndraws` still subsets, and does so on the draws of the mean.
  expect_identical(length(unique(binned[["data"]][["rep_id"]])), 2L)

  # A check that does compare replicates with the response keeps getting them.
  reps <- suppressMessages(suppressWarnings(
    bayesplot::pp_check(fit, type = "hist", ndraws = 2L)))

  expect_s3_class(reps, "ggplot")
  expect_true(all(reps[["data"]][["value"]] %in% c(0, 1)))
})

# An accelerated failure time family reports the density of log(T) and `ph()`
# the density of T, so a log score taken across that boundary is off by the
# Jacobian, by thousands of points, and the ordering it produces can be the
# wrong way round. `scale` puts both on one measure.
test_that("loo(scale=) puts the survival families on one measure", {
  skip_on_cran()
  skip_if_not_installed("loo")
  skip_if_not_installed("survival")

  d <- sim_x(n = 150L, p = 3L, seed = 81L)
  d$time <- stats::rexp(nrow(d), rate = exp(-1 - d$x1))
  d$status <- stats::rbinom(nrow(d), 1L, 0.7)

  f <- survival::Surv(time, status) ~ x1 + x2 + x3
  ctrl <- quick_control(num_trees = 5L, num_burn = 60L, num_draws = 100L)

  aft <- bartisan(f, d, family = weibull_aft(), control = ctrl)
  prop_haz <- bartisan(f, d, family = ph(), control = ctrl)

  elpd <- function(x) unname(x[["estimates"]]["elpd_loo", "Estimate"])
  quiet <- function(e) suppressWarnings(suppressMessages(e))

  jacobian <- sum(d$status * log(d$time))

  # The correction is the Jacobian, on events only.
  expect_equal(elpd(quiet(loo::loo(aft, scale = "time"))),
               elpd(quiet(loo::loo(aft))) - jacobian)

  expect_equal(elpd(quiet(loo::loo(prop_haz, scale = "log_time"))),
               elpd(quiet(loo::loo(prop_haz))) + jacobian)

  # A fit already on the scale named is returned untouched, which is what lets
  # one `scale` be named for every model in a comparison.
  expect_equal(elpd(quiet(loo::loo(aft, scale = "log_time"))),
               elpd(quiet(loo::loo(aft))))

  expect_equal(elpd(quiet(loo::loo(prop_haz, scale = "time"))),
               elpd(quiet(loo::loo(prop_haz))))

  # And the comparison does not depend on which of the two scales is chosen,
  # since a constant per observation cancels from the difference.
  on_time <- elpd(quiet(loo::loo(prop_haz, scale = "time"))) -
    elpd(quiet(loo::loo(aft, scale = "time")))

  on_log <- elpd(quiet(loo::loo(prop_haz, scale = "log_time"))) -
    elpd(quiet(loo::loo(aft, scale = "log_time")))

  expect_equal(on_time, on_log)

  # `waic()` takes it too, and a family with only one scale rejects it.
  expect_s3_class(quiet(loo::waic(aft, scale = "time")), "waic")

  gaussian_fit <- bartisan(y ~ ., sim_x(n = 120L, p = 2L, seed = 82L) |>
                             transform(y = stats::rnorm(120L)),
                           family = stats::gaussian(), control = ctrl)

  expect_error(loo::loo(gaussian_fit, scale = "time"), "names the measure")

  # A censored time is a bound rather than a value, so there is no residual on
  # the response scale to take, and no residual variance for an R2 to rest on.
  expect_error(stats::residuals(prop_haz), "censored")

  if (rlang::is_installed("performance")) {
    expect_warning(r2 <- performance::r2_posterior(prop_haz), "censored")
    expect_null(r2)
  }
})

# `loo()` estimates the leave-one-out density from one fit; `kfold()` refits and
# does not estimate it. The contract that matters is that `loo_compare()` takes
# the result beside a `<loo>` object, which needs three rows in `estimates`: with
# one, `loo_compare()` flattens a length-2 vector beside a length-6 one and
# `sapply()` returns a list rather than a matrix.
test_that("kfold() returns what loo_compare() accepts", {
  skip_on_cran()
  skip_if_not_installed("loo")
  skip_if_not_installed("rstantools")

  d <- sim_x(n = 150L, p = 3L, seed = 101L)
  d$y <- stats::rbinom(nrow(d), 1L, stats::plogis(d$x1))

  ctrl <- quick_control(num_trees = 5L, num_burn = 50L, num_draws = 80L)
  fit <- bartisan(y ~ ., data = d, family = stats::binomial(), control = ctrl)

  kf <- loo::kfold(fit, K = 3L)

  expect_s3_class(kf, "kfold")
  expect_true(loo::is.kfold(kf))
  expect_identical(attr(kf, "K"), 3L)

  expect_identical(rownames(kf[["estimates"]]),
                   c("elpd_kfold", "p_kfold", "kfoldic"))
  expect_identical(colnames(kf[["estimates"]]), c("Estimate", "SE"))
  expect_identical(nrow(kf[["pointwise"]]), nrow(d))

  # The estimate is the sum of the pointwise values, and `kfoldic` is -2 times
  # the elpd, which is what makes the three rows one object rather than three.
  expect_equal(kf[["estimates"]]["elpd_kfold", "Estimate"],
               sum(kf[["pointwise"]][, "elpd_kfold"]))
  expect_equal(kf[["estimates"]]["kfoldic", "Estimate"],
               -2 * kf[["estimates"]]["elpd_kfold", "Estimate"])

  # Every observation is scored exactly once, by a fit that did not see it.
  expect_setequal(kf[["folds"]], 1:3)
  expect_false(anyNA(kf[["pointwise"]]))

  # The contract itself, against a `<loo>` object rather than another `<kfold>`.
  cmp <- suppressWarnings(loo::loo_compare(list(kfold = kf,
                                                loo = loo::loo(fit))))
  expect_s3_class(cmp, "compare.loo")
  expect_true(all(c("elpd_diff", "se_diff") %in% colnames(cmp)))
  expect_identical(nrow(cmp), 2L)

  # Held out, a model predicts worse than it does in sample; that gap is what
  # `p_kfold` is.
  expect_gt(kf[["estimates"]]["p_kfold", "Estimate"], 0)

  # Folds may be supplied, which is how two models are scored on one split.
  same <- loo::kfold(fit, folds = kf[["folds"]])
  expect_identical(same[["folds"]], kf[["folds"]])

  expect_error(loo::kfold(fit, folds = c(1L, 2L)), "one fold per observation")
  expect_error(loo::kfold(fit, folds = rep(c(1L, 3L), length.out = nrow(d))),
               "every fold")

  # `save_fits` keeps the refits, and does not by default.
  expect_null(kf[["fits"]])
  expect_length(loo::kfold(fit, K = 2L, save_fits = TRUE)[["fits"]], 2L)
})

# `kfold()` refits from the recorded call, whose arguments were evaluated in the
# formula's environment. A fit made inside `lapply(trees, function(n) ...)` has
# `num_trees = n` in its call, and `n` lives only in the anonymous function's
# frame: with no `n` elsewhere the refit failed, and with one it silently refit
# a different model. `vignette("comparison")` hit the first; the second is the
# worse of the two and is the one set up here, with an `n` beside the formula.
test_that("kfold() refits a model made inside lapply() as it was made", {
  skip_on_cran()
  skip_if_not_installed("loo")

  d <- sim_x(n = 120L, p = 2L, seed = 103L)
  d$y <- stats::rnorm(nrow(d), sin(3 * d$x1))

  model <- y ~ x1 + x2
  n <- 3L
  folds <- rep(1:2, length.out = nrow(d))

  fit_one <- function(n) {
    bartisan(model, data = d, family = stats::gaussian(), num_trees = n,
             num_burn = 20L, num_draws = 20L, verbose = FALSE)
  }

  set.seed(1L)
  inside <- lapply(c(6L, 9L), fit_one)

  refits <- lapply(inside, loo::kfold, folds = folds, save_fits = TRUE)

  expect_identical(vapply(refits, function(k) k[["fits"]][[1L]][["num_trees"]],
                          integer(1L)),
                   c(6L, 9L))

  # And the same as the fit made outside a loop, draw for draw.
  set.seed(1L)
  outside <- bartisan(model, data = d, family = stats::gaussian(),
                      num_trees = 6L, num_burn = 20L, num_draws = 20L,
                      verbose = FALSE)
  expect_identical(outside[["eta"]], inside[[1L]][["eta"]])

  set.seed(2L)
  a <- loo::kfold(outside, folds = folds)
  set.seed(2L)
  b <- loo::kfold(inside[[1L]], folds = folds)
  expect_identical(a[["pointwise"]], b[["pointwise"]])

  # The data is the one argument not kept, so data that lived only in the loop
  # is an error saying how to make a fit that can be refit.
  local_data <- lapply(list(d), function(dd) {
    bartisan(model, data = dd, family = stats::gaussian(), num_trees = 5L,
             num_burn = 10L, num_draws = 10L, verbose = FALSE)
  })
  expect_error(loo::kfold(local_data[[1L]], folds = folds),
               "cannot be found from where its formula was written")
})

# A `bcf()` fit's forests split on a propensity score the caller never named, and
# the held-out density rebuilt the response from the full terms without it, so
# `predict(type = "density")` on new data, and `kfold()` with it, failed on
# `.propensity` for every such fit.
test_that("kfold() and held-out densities work on a bcf() fit", {
  skip_on_cran()
  skip_if_not_installed("loo")

  d <- sim_x(n = 120L, p = 2L, seed = 104L)
  d$z <- stats::rbinom(nrow(d), 1L, 0.5)
  d$y <- stats::rnorm(nrow(d), sin(3 * d$x1) + d$z)

  fit <- lapply(4L, function(m) {
    bcf(y ~ x1 + x2, treat = ~ z, data = d, family = stats::gaussian(),
        num_trees = c(m, 3L), num_burn = 20L, num_draws = 20L,
        verbose = FALSE,
        propensity_args = list(num_trees = 5L, num_burn = 10L,
                               num_draws = 10L))
  })[[1L]]

  dens <- stats::predict(fit, newdata = d[1:5, c("y", "x1", "x2", "z")],
                         type = "density", log = TRUE)
  expect_length(dens, 5L)
  expect_true(all(is.finite(dens)))

  kf <- loo::kfold(fit, folds = rep(1:2, length.out = nrow(d)),
                   save_fits = TRUE)
  expect_false(anyNA(kf[["pointwise"]]))
  expect_identical(kf[["fits"]][[1L]][["num_trees"]], c(4L, 3L))
})

# A score given to `bcf()` as numbers has a row for each row of the data. Each
# refit has to be given the training rows of it, and the held-out rows have to
# be scored with theirs, since a score that was supplied rather than fitted
# cannot be rebuilt for new rows. The refits were given all of it, and stopped.
test_that("kfold() takes a supplied propensity score fold by fold", {
  skip_if_not_installed("loo")
  skip_if_not_installed("patrick")

  d <- sim_x(n = 60L, p = 2L, seed = 109L)
  d$z <- stats::rbinom(nrow(d), 1L, stats::plogis(d$x1))
  d$g <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  d$y <- stats::rnorm(nrow(d), d$x1 + d$z)

  folds <- rep(1:2, length.out = nrow(d))
  held <- folds == 1L

  # One probability per level, each varying, since a constant column is dropped.
  level_scores <- exp(cbind(a = d$x1, b = d$x2, c = 0))
  level_scores <- level_scores / rowSums(level_scores)

  patrick::with_parameters_test_that(
    "score given as",
    {
      fit <- bcf(y ~ x1 + x2, treat = treat, data = d,
                 family = stats::gaussian(), propensity = score,
                 num_trees = 4L, num_burn = 10L, num_draws = 10L,
                 verbose = FALSE)

      kf <- loo::kfold(fit, folds = folds, save_fits = TRUE)
      expect_false(anyNA(kf[["pointwise"]]))

      # The first refit is trained on the second fold and scores the first.
      refit <- kf[["fits"]][[1L]]
      expect_equal(refit[["bcf"]][["propensity"]], ss(as.matrix(score), !held),
                   ignore_attr = TRUE)

      scored <- cbind(d[held, ], ss(fit[["bcf"]][["propensity"]], held))
      expect_equal(kf[["pointwise"]][held, "elpd_kfold"],
                   stats::predict(refit, newdata = scored, type = "density",
                                  log = TRUE),
                   ignore_attr = TRUE)
    },
    patrick::cases(
      vector = list(treat = ~ z, score = stats::plogis(d$x1)),
      matrix = list(treat = ~ g, score = level_scores)
    )
  )
})

# `predict(type = "density")` falls back to the fit's own prior weights when it
# is given none, so a weighted fit scored on held-out rows without them came
# back wrong rather than erroring. The refits and the scores both have to carry
# them.
test_that("kfold() carries prior weights into the refits and the scores", {
  skip_on_cran()
  skip_if_not_installed("loo")
  skip_if_not_installed("rstantools")

  d <- sim_x(n = 120L, p = 2L, seed = 102L)
  d$trials <- 6
  d$y <- stats::rbinom(nrow(d), 6L, stats::plogis(d$x1)) / 6

  ctrl <- quick_control(num_trees = 5L, num_burn = 50L, num_draws = 80L)
  fit <- bartisan(y ~ x1 + x2, data = d, family = stats::binomial(),
                  weights = trials, control = ctrl)

  kf <- loo::kfold(fit, K = 2L)

  expect_false(anyNA(kf[["pointwise"]]))

  # Six trials per row makes each contribution a sum of six Bernoulli terms, so
  # a score taken as though there were one trial is far too close to zero. The
  # in-sample log density is the scale to judge that against.
  in_sample <- sum(stats::predict(fit, type = "density", log = TRUE))
  held_out <- kf[["estimates"]]["elpd_kfold", "Estimate"]

  expect_lt(held_out, 0)
  expect_gt(held_out, 4 * in_sample)
  expect_lt(held_out, in_sample / 4)
})

# `kfold()` rebuilds each refit from the fit's call and the data that call names,
# so the checks here are about finding them: recorded or not, given as a list,
# changed since, or missing, and with a subset that would otherwise be applied
# twice. Every comparison is between two cross-validations on the same folds and
# the same seed, which is exact whatever the seed.
test_that("kfold() refits from the data and the call the fit was made from", {
  skip_if_not_installed("loo")

  d <- sim_x(n = 60L, p = 2L, seed = 105L)
  d$y <- stats::rnorm(nrow(d), d$x1)

  ctrl <- quick_control(num_burn = 10L, num_draws = 10L)
  folds <- rep(1:2, length.out = nrow(d))

  fit <- bartisan(y ~ x1 + x2, data = d, family = stats::gaussian(),
                  control = ctrl)

  expect_error(loo::kfold(fit, folds = replace(folds, 1L, 0L)), "from 1 up")
  expect_error(loo::kfold(fit, folds = replace(folds, 1L, NA)), "from 1 up")

  elpd <- function(x) {
    set.seed(1L)
    loo::kfold(x, folds = folds)[["pointwise"]][, "elpd_kfold"]
  }

  kept <- elpd(fit)

  # A fit that recorded no argument values, as one made before they were kept,
  # is refit from its call evaluated where its formula was written.
  bare <- fit
  bare[["call_values"]] <- NULL
  expect_identical(elpd(bare), kept)

  # Data given as a list is the data frame it describes.
  listed <- bartisan(y ~ x1 + x2, data = as.list(d), family = stats::gaussian(),
                     control = ctrl)
  expect_identical(elpd(listed), kept)

  # Data that has lost rows since the fit cannot be refit as it was made.
  changed <- local({
    fit <- bartisan(y ~ x1 + x2, data = d, family = stats::gaussian(),
                    control = ctrl)
    d <- d[-(1:5), ]
    fit
  })
  expect_error(loo::kfold(changed, folds = folds), "rows are gone")

  # And a fit whose call names no data has nothing to take folds from.
  loose <- local({
    y <- d$y
    x1 <- d$x1
    x2 <- d$x2
    bartisan(y ~ x1 + x2, family = stats::gaussian(), control = ctrl)
  })
  expect_error(loo::kfold(loose, folds = folds), "names no")

  # A subset given as a vector the length of the full data would be applied a
  # second time to training data that no longer has that many rows. The rows are
  # chosen by name instead, so each one is scored once and fit once.
  keep <- d$x1 > 0.3
  subsetted <- bartisan(y ~ x1 + x2, data = d, subset = keep,
                        family = stats::gaussian(), control = ctrl)
  kf <- loo::kfold(subsetted, K = 2L, save_fits = TRUE)

  expect_identical(nrow(kf[["pointwise"]]), sum(keep))
  expect_equal(sum(vapply(kf[["fits"]], stats::nobs, numeric(1L))), sum(keep))
})

# The offset is taken from the model frame rather than re-evaluated, so it has to
# reach each refit and each score, and once. Written in the formula it is already
# in the formula each refit is made from, and passing the whole of it as well
# counted that part twice; given as a matrix, one column per additive predictor,
# it has to be taken by rows.
test_that("kfold() carries an offset into the refits and the scores once", {
  skip_if_not_installed("loo")
  skip_if_not_installed("patrick")

  d <- sim_x(n = 60L, p = 2L, seed = 106L)
  d$exposure <- stats::runif(nrow(d), 1, 50)
  d$scale <- stats::runif(nrow(d), 0.5, 2)
  d$y <- stats::rpois(nrow(d), d$exposure * d$scale * exp(d$x1 - 2))
  d$zeros <- ifelse(stats::runif(nrow(d)) < 0.3, 0, d$y)

  ctrl <- quick_control(num_burn = 10L, num_draws = 10L)
  folds <- rep(1:2, length.out = nrow(d))
  held <- folds == 1L
  total <- log(d$exposure) + log(d$scale)

  rows_of <- function(v, i) {
    if (is.matrix(v)) v[i, , drop = FALSE] else v[i]
  }

  patrick::with_parameters_test_that(
    "offset written as",
    {
      fit <- eval(model)
      kf <- loo::kfold(fit, folds = folds, save_fits = TRUE)

      # The first refit is trained on the second fold and scores the first.
      refit <- kf[["fits"]][[1L]]

      expect_equal(stats::model.offset(refit[["model"]]),
                   rows_of(expected, !held), ignore_attr = TRUE)
      expect_equal(kf[["pointwise"]][held, "elpd_kfold"],
                   stats::predict(refit, newdata = d[held, ], type = "density",
                                  log = TRUE, offset = rows_of(expected, held)),
                   ignore_attr = TRUE)
    },
    patrick::cases(
      argument = list(
        model = quote(bartisan(y ~ x1 + x2, data = d, family = stats::poisson(),
                               offset = log(exposure) + log(scale),
                               control = ctrl)),
        expected = total),
      formula = list(
        model = quote(bartisan(y ~ x1 + x2 + offset(log(exposure) + log(scale)),
                               data = d, family = stats::poisson(),
                               control = ctrl)),
        expected = total),
      `formula and argument` = list(
        model = quote(bartisan(y ~ x1 + x2 + offset(log(exposure)), data = d,
                               family = stats::poisson(), offset = log(scale),
                               control = ctrl)),
        expected = total),
      matrix = list(
        model = quote(bartisan(zeros ~ x1 + x2, data = d, family = zi_poisson(),
                               offset = cbind(log(exposure) + log(scale), 0),
                               control = ctrl)),
        expected = cbind(total, 0))
    )
  )
})

test_that("kfold(scale=) moves the held-out and in-sample scores together", {
  skip_if_not_installed("loo")
  skip_if_not_installed("survival")

  d <- sim_x(n = 60L, p = 2L, seed = 107L)
  d$time <- stats::rexp(nrow(d), exp(-d$x1))
  d$status <- stats::rbinom(nrow(d), 1L, 0.7)
  folds <- rep(1:2, length.out = nrow(d))

  fit <- bartisan(survival::Surv(time, status) ~ x1 + x2, data = d,
                  family = weibull_aft(),
                  control = quick_control(num_burn = 10L, num_draws = 10L))

  set.seed(1L)
  log_time <- loo::kfold(fit, folds = folds)[["pointwise"]]
  set.seed(1L)
  on_time <- loo::kfold(fit, folds = folds, scale = "time")[["pointwise"]]

  # The Jacobian of the change of measure, on events only, as for `loo()`.
  expect_equal(on_time[, "elpd_kfold"],
               log_time[, "elpd_kfold"] - d$status * log(d$time))

  # Both sides move by it, so the gap between them does not.
  expect_equal(on_time[, "p_kfold"], log_time[, "p_kfold"])
})

# The folds go to workers when a plan has any, with the streams drawn before the
# branch, so the scores do not depend on whether a plan was set.
test_that("kfold() scores the same folds in parallel as in sequence", {
  skip_on_cran()
  skip_if_not_installed("loo")
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)

  d <- sim_x(n = 60L, p = 2L, seed = 108L)
  d$y <- stats::rnorm(nrow(d), d$x1)
  folds <- rep(1:2, length.out = nrow(d))

  fit <- bartisan(y ~ x1 + x2, data = d, family = stats::gaussian(),
                  control = quick_control(num_burn = 10L, num_draws = 10L))

  set.seed(1L)
  one <- loo::kfold(fit, folds = folds)

  # A machine that reports one core warns about the load two workers put on it,
  # which says nothing about the code under test.
  started <- tryCatch({
    suppressWarnings(future::plan(future::multisession, workers = 2L))
    isTRUE(future::nbrOfWorkers() >= 2L)
  }, error = function(e) FALSE)

  skip_if_not(started, "no second worker available")

  set.seed(1L)
  many <- loo::kfold(fit, folds = folds)

  expect_equal(many[["pointwise"]], one[["pointwise"]])
})
