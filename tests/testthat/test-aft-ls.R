# The two location-scale accelerated failure time families: `lognormal_aft_ls()`
# and `dpm_aft_ls()`, which is NFT BART. Both put a second forest on the log of
# the scale of log T, so survival curves for different covariates can cross.
# The tests are about the identities that tie what they report together: the
# densities against the sampler's log likelihood, the survival curve against
# its closed form, and the hazard against the survival curve.

sim_crossing <- function(n, seed) {
  set.seed(seed)
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  d$sd <- exp(-1 + 1.2 * d$x2)
  log_t <- 1 + 0.6 * d$x1 + d$sd * stats::rnorm(n)
  cens <- log(stats::runif(n, 1, 12))
  d$time <- exp(pmin(log_t, cens))
  d$status <- as.numeric(log_t <= cens)
  d
}

patrick::with_parameters_test_that(
  "a location-scale survival fit agrees with itself:",
  {
    d <- sim_crossing(150, 1501)

    fit <- bartisan(list(cbind(time, status) ~ x1 + x2, log_sd = ~ x1 + x2), d,
                    family = family,
                    control = quick_control(num_trees = c(10L, 5L),
                                            num_burn = 60L, num_draws = 60L))

    expect_identical(names(fit[["eta"]]), c("mean", "log_sd"))

    # The observed-data density, an event's density and a censoring's survival
    # probability, reproduces the log likelihood the sampler reported.
    dens <- stats::predict(fit, type = "density", draws = TRUE, log = TRUE)
    expect_equal(rowSums(dens), fit[["loglik"]], tolerance = 1e-8)

    eta <- stats::predict(fit, type = "link", draws = TRUE)
    replayed <- stats::predict(fit, newdata = d, type = "link", draws = TRUE)
    expect_equal(replayed, eta, tolerance = 1e-6, ignore_attr = TRUE)

    dens_new <- stats::predict(fit, newdata = d, type = "density",
                               draws = TRUE, log = TRUE)
    expect_equal(rowSums(dens_new), fit[["loglik"]], tolerance = 1e-6)

    # The survival curve decreases, the median survival time is where it
    # crosses one half, and the hazard is -d log S / dt.
    times <- c(1, 2, 4, 8)
    surv <- stats::predict(fit, type = "survival", times = times, draws = TRUE)
    expect_true(all(surv > 0 & surv < 1))
    expect_true(all(apply(surv, c(1L, 2L), function(z) all(diff(z) <= 1e-12))))

    median <- stats::predict(fit, type = "response", draws = TRUE)
    at_median <- vapply(seq_len(ncol(median)), function(i) {
      vapply(seq_len(nrow(median)), function(s) {
        stats::predict(fit, type = "survival", times = median[s, i],
                       draws = TRUE, iterations = s)[1L, i, 1L]
      }, numeric(1L))
    }, numeric(nrow(median)))[, 1:5]
    expect_equal(at_median, matrix(0.5, nrow(median), 5L), tolerance = 1e-6,
                 ignore_attr = TRUE)

    step <- 1e-4
    hazard <- stats::predict(fit, type = "hazard", times = 3, draws = TRUE)
    log_s <- function(t) {
      log(stats::predict(fit, type = "survival", times = t, draws = TRUE))
    }
    difference <- -(log_s(3 + step) - log_s(3 - step)) / (2 * step)
    expect_equal(hazard, difference, tolerance = 1e-4, ignore_attr = TRUE)
  },
  patrick::cases(
    lognormal = list(family = lognormal_aft_ls()),
    nft = list(family = dpm_aft_ls())
  )
)

test_that("`lognormal_aft_ls()` reports the closed-form survival curve", {
  d <- sim_crossing(120, 1502)

  fit <- bartisan(list(cbind(time, status) ~ x1 + x2, log_sd = ~ x2), d,
                  family = lognormal_aft_ls(),
                  control = quick_control(num_trees = c(10L, 5L)))

  eta <- stats::predict(fit, type = "link", draws = TRUE)
  surv <- stats::predict(fit, type = "survival", times = 2, draws = TRUE)
  expect_equal(surv[, , 1L],
               stats::pnorm((log(2) - eta[["mean"]]) / exp(eta[["log_sd"]]),
                            lower.tail = FALSE),
               ignore_attr = TRUE)

  # The median is exp(mean), whatever the scale.
  expect_equal(stats::predict(fit, type = "response", draws = TRUE),
               exp(eta[["mean"]]), ignore_attr = TRUE)
})

test_that("`lognormal_aft_ls()` takes prior weights and `dpm_aft_ls()` refuses them", {
  d <- sim_crossing(80, 1503)
  w <- stats::runif(nrow(d), 0.5, 1.5)
  ctrl <- quick_control(num_trees = c(5L, 5L))

  # With weights the censored times are not imputed, so the fit takes the
  # direct likelihood; its density is the same function either way, and at
  # the fitted data it carries the weights.
  fit <- bartisan(list(cbind(time, status) ~ x1, log_sd = ~ x2), d,
                  family = lognormal_aft_ls(), weights = w, control = ctrl)
  dens <- stats::predict(fit, type = "density", draws = TRUE, log = TRUE)
  expect_equal(rowSums(dens), fit[["loglik"]], tolerance = 1e-8)

  expect_error(bartisan(list(cbind(time, status) ~ x1, log_sd = ~ x2), d,
                        family = dpm_aft_ls(), weights = w, control = ctrl),
               "does not take prior weights")
})

test_that("the scale forest of `dpm_aft_ls()` recovers a spread that varies", {
  skip_on_cran()

  d <- sim_crossing(600, 1504)

  fit <- bartisan(list(cbind(time, status) ~ x1 + x2, log_sd = ~ x1 + x2), d,
                  family = dpm_aft_ls(),
                  control = quick_control(num_trees = c(50L, 10L),
                                          num_burn = 400L, num_draws = 400L))

  log_sd <- stats::predict(fit, type = "link")[, "log_sd"]
  expect_gt(stats::cor(log_sd, log(d$sd)), 0.85)
  expect_lt(abs(mean(exp(log_sd) / d$sd) - 1), 0.2)
})

patrick::with_parameters_test_that(
  "a Dirichlet process survival fit draws replicates from its own distribution:",
  {
    skip_if_not_installed("rstantools")

    d <- sim_crossing(150, 1505)
    formula <- if (scaled) {
      list(cbind(time, status) ~ x1 + x2, log_sd = ~ x1 + x2)
    } else {
      cbind(time, status) ~ x1 + x2
    }

    fit <- bartisan(formula, d, family = family,
                    control = quick_control(num_trees = if (scaled) c(10L, 5L)
                                                        else 10L,
                                            num_burn = 60L, num_draws = 200L))

    # Within a draw a replicate exceeds t with probability S(t | x) for that
    # draw, so the share of replicates above t, pooled over draws and units,
    # matches the mean of the survival curve up to Monte Carlo error.
    set.seed(15)
    reps <- rstantools::posterior_predict(fit)
    expect_identical(dim(reps), c(200L, nrow(d)))
    expect_true(all(reps > 0))

    times <- c(2, 4, 8)
    surv <- stats::predict(fit, type = "survival", times = times, draws = TRUE)

    for (k in seq_along(times)) {
      p <- surv[, , k]
      se <- sqrt(mean(p * (1 - p)) / length(p))
      expect_lt(abs(mean(reps > times[k]) - mean(p)), 4 * se)
    }

    # On new data too, where the predictor is rebuilt from the trees.
    new_reps <- rstantools::posterior_predict(fit, newdata = d[1:20, ])
    expect_identical(dim(new_reps), c(200L, 20L))
  },
  patrick::cases(
    dpm_aft = list(family = dpm_aft(), scaled = FALSE),
    dpm_aft_ls = list(family = dpm_aft_ls(), scaled = TRUE)
  )
)
