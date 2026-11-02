# The heteroskedastic mixture, `dpm_ls()`, which scales the whole error by
# exp(log_sd). The tests are about the chart it reports in, and about the one
# thing it adds that nothing else in the package does: a mean whose reporting
# shift differs by observation, which the engine writes into the recorded draws
# and R has to apply again to anything it rebuilds from the trees.

sim_fanning <- function(n, seed) {
  set.seed(seed)
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  d$sd <- exp(-1 + 1.5 * d$x2)
  # A skewed error with mean zero and variance one, so `sd` is the conditional
  # standard deviation and the shape is one a single normal cannot take.
  d$y <- sin(3 * d$x1) +
    d$sd * (stats::rgamma(n, shape = 1.5) - 1.5) / sqrt(1.5)
  d
}

test_that("a dpm_ls() fit reports in its chart and agrees with itself", {
  d <- sim_fanning(150, 1401)

  fit <- bartisan(list(y ~ x1 + x2, log_sd = ~ x1 + x2), d, family = dpm_ls(),
                  control = quick_control(num_trees = c(10L, 5L),
                                          num_burn = 60L, num_draws = 60L))

  expect_identical(colnames(fit[["aux"]]),
                   c("alpha", "clusters", "center", "spread", "scale",
                     "error_sd"))

  # The stored mixture is centered, with unit variance, which is the unit the
  # scale forest's shift took out.
  for (s in c(1L, 30L, 60L)) {
    components <- mixture_at(fit, s)
    weight <- components[, "weight"]
    expect_equal(sum(weight * components[, "mean"]), 0, tolerance = 1e-10)
    expect_equal(sum(weight * (components[, "sd"]^2 + components[, "mean"]^2)),
                 1, tolerance = 1e-10)
  }

  # The density at the recorded predictors reproduces the log likelihood the
  # sampler computed in its own chart, so the recorded predictors and the
  # stored mixture are in one chart together.
  dens <- stats::predict(fit, type = "density", draws = TRUE, log = TRUE)
  expect_equal(rowSums(dens), fit[["loglik"]], tolerance = 1e-8)

  # Rebuilt from the trees the predictors match the recorded ones, which is the
  # test of `report_chart()`, the part of the mean no leaf can carry.
  eta <- stats::predict(fit, type = "link", draws = TRUE)
  replayed <- stats::predict(fit, newdata = d, type = "link", draws = TRUE)
  expect_equal(replayed, eta, tolerance = 1e-6, ignore_attr = TRUE)

  dens_new <- stats::predict(fit, newdata = d, type = "density",
                             draws = TRUE, log = TRUE)
  expect_equal(rowSums(dens_new), fit[["loglik"]], tolerance = 1e-6)

  expect_equal(stats::predict(fit, type = "response"),
               stats::predict(fit, type = "link")[, "mean"],
               ignore_attr = TRUE)
})

test_that("the scale forest of dpm_ls() is the conditional standard deviation", {
  skip_on_cran()

  d <- sim_fanning(600, 1402)

  fit <- bartisan(list(y ~ x1 + x2, log_sd = ~ x1 + x2), d, family = dpm_ls(),
                  control = quick_control(num_trees = c(50L, 10L),
                                          num_burn = 400L, num_draws = 400L))

  log_sd <- stats::predict(fit, type = "link")[, "log_sd"]
  expect_gt(stats::cor(log_sd, log(d$sd)), 0.9)
  expect_lt(abs(mean(exp(log_sd) / d$sd) - 1), 0.15)

  # The error density is drawn at the scale the second predictor has on average,
  # where W, which has unit variance, is stretched by `scale`; it integrates to
  # one, is centered, and keeps the skew the errors were drawn with.
  grid <- seq(-6, 8, length.out = 1401L)
  density <- error_density(fit, at = grid)[, "mean"]
  step <- diff(grid)[1L]
  moment <- function(k) sum(grid^k * density) * step

  expect_equal(moment(0), 1, tolerance = 1e-3)
  expect_equal(moment(1), 0, tolerance = 0.05)
  expect_equal(sqrt(moment(2)), mean(fit[["aux"]][, "scale"]),
               tolerance = 0.05)
  expect_gt(moment(3), 0)
})

test_that("an empty scale forest makes dpm_ls() a dpm()", {
  skip_on_cran()

  # `?bartisan-families`: with `~ 1` on the scale, every tree in that forest is
  # a stump and exp(log_sd) one number per draw multiplying W, which is the model
  # `dpm()` fits. The two cannot agree exactly, because the scalar is drawn under
  # the leaf prior, but they fit the same mean and the same density.
  d <- sim_fanning(400, 1403)
  d$y <- sin(3 * d$x1) + 0.5 * (stats::rgamma(nrow(d), shape = 1.5) - 1.5) /
    sqrt(1.5)
  control <- quick_control(num_trees = c(30L, 5L), num_burn = 300L,
                           num_draws = 300L)

  two <- bartisan(list(y ~ x1 + x2, log_sd = ~ 1), d, family = dpm_ls(),
                  control = control)
  one <- bartisan(y ~ x1 + x2, d, family = dpm(),
                  control = quick_control(num_trees = 30L, num_burn = 300L,
                                          num_draws = 300L))

  scale <- stats::predict(two, type = "link", draws = TRUE)[["log_sd"]]
  expect_lt(max(apply(scale, 1L, stats::sd)), 1e-8)

  expect_gt(stats::cor(stats::predict(two, type = "response"),
                       stats::predict(one, type = "response")), 0.98)
  # Per observation, where the two came out 0.014 apart.
  expect_lt(abs(mean(two[["loglik"]]) - mean(one[["loglik"]])) / nrow(d),
            0.03)
})

test_that("dpm_ls() refuses prior weights, and says why", {
  d <- sim_fanning(40, 1404)

  expect_error(bartisan(y ~ x1, d, family = dpm_ls(),
                        weights = stats::runif(40),
                        control = quick_control()),
               "does not take prior weights")
})
