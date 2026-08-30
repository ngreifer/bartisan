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
                  control = quick_control(num_burn = 300L, num_save = 300L))

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

test_that("a binomial replicate is a fraction of the trials", {
  skip_if_not_installed("rstantools")

  d <- sim_x(200, 2)
  set.seed(1231)
  d$y <- rbinom(nrow(d), 6L, stats::plogis(2 * d$x1 - d$x2)) / 6
  d$trials <- 6

  fit <- bartisan(y ~ x1 + x2, data = d, family = binomial(), weights = trials,
                  control = quick_control())

  # The number of trials is not a function of the predictors, so predicting for
  # new rows has to be told what it is rather than assuming one trial.
  expect_error(rstantools::posterior_predict(fit, newdata = d[1L, ]),
               "number of trials")

  set.seed(1232)
  drawn <- replicates_at_one_draw(fit, d[1L, ], weights = 6)
  expect_true(all(drawn %in% ((0:6) / 6)))

  nd <- d[rep(1L, 7L), ]
  nd$y <- (0:6) / 6
  probability <- as.vector(
    stats::predict(fit, newdata = nd, type = "density", iterations = 1L,
                   draws = TRUE, weights = rep(6, 7L)))
  empirical <- vapply((0:6) / 6, function(g) mean(abs(drawn - g) < 1e-8),
                      numeric(1L))
  expect_lt(max(abs(empirical - probability)), 0.03)

  # Binary data are the same statement with one trial, so they come back as
  # zeros and ones rather than as counts.
  set.seed(1233)
  d$y <- rbinom(nrow(d), 1L, stats::plogis(2 * d$x1 - d$x2))
  fit <- bartisan(y ~ x1 + x2, data = d, family = binomial(),
                  control = quick_control())
  expect_setequal(unique(as.vector(rstantools::posterior_predict(fit))), c(0, 1))
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
                  control = quick_control(num_save = 400L))

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
})

test_that("loo and waic run on the pointwise likelihood", {
  skip_if_not_installed("loo")
  skip_if_not_installed("rstantools")

  d <- sim_x(150, 2)
  set.seed(1321)
  d$y <- rpois(nrow(d), exp(0.5 + 2 * d$x1))

  fit <- bartisan(y ~ x1 + x2, data = d, family = poisson(), chains = 2L,
                  control = quick_control(num_save = 100L))

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

  d <- sim_x(200, 2)
  set.seed(1331)
  d$y <- 3 * d$x1 - 2 * d$x2 + rnorm(nrow(d), sd = 0.5)

  fit <- bartisan(y ~ x1 + x2, data = d, control = quick_control(num_save = 60L))

  draws <- performance::r2_posterior(fit)[["R2_Bayes"]]
  expect_length(draws, 60L)
  expect_true(all(draws > 0 & draws < 1))

  # A strong signal, so it should be well away from zero.
  expect_gt(mean(draws), 0.5)

  summary_r2 <- performance::r2(fit)
  expect_equal(unname(summary_r2[["R2_Bayes"]]), stats::median(draws),
               tolerance = 1e-8)

  # A categorical response has no mean, so there is no such ratio to report.
  set.seed(1332)
  d$y <- factor(findInterval(3 * d$x1 + rlogis(nrow(d)), c(0, 1)) + 1L,
                labels = c("lo", "mid", "hi"), ordered = TRUE)
  ordinal_fit <- bartisan(y ~ x1 + x2, data = d, family = ordinal(),
                          control = quick_control())
  expect_warning(out <- performance::r2_posterior(ordinal_fit), "no mean")
  expect_null(out)
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

test_that("as_draws hands the scalar parameters over with their chain structure", {
  skip_if_not_installed("posterior")

  d <- sim_x(120, 2)
  set.seed(1351)
  d$y <- 2 * d$x1 + rnorm(nrow(d))

  # Named, because `aux.sigma` below is the Gaussian family's nuisance parameter
  # and the numeric default is now the mixture.
  fit <- bartisan(y ~ x1 + x2, data = d, family = gaussian(), chains = 3L,
                  control = quick_control(num_save = 40L))

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
  # are computed from the same vectors by a different implementation.
  summary_draws <- posterior::summarise_draws(draws)
  own <- fit[["rhat"]]
  expect_equal(summary_draws$rhat[summary_draws$variable == "loglik"],
               own$rhat[own$quantity == "loglik"], tolerance = 0.01)
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
                  control = quick_control(num_save = 60L))

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
