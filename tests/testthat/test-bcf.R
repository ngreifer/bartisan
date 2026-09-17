# `bcf()` is a wrapper over the varying-coefficient interface. Everything it does
# is expressible in `bartisan()` with a `vc()` term, which is what these check:
# the wrapper's decisions, not the sampler underneath it.

sim_causal <- function(n = 600, seed = 1, kind = "binary") {
  set.seed(seed)
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))

  d$z <- switch(kind,
                binary = stats::rbinom(n, 1L, stats::plogis(0.8 * d$x1)),
                categorical = factor(sample(c("a", "b", "c"), n, TRUE)),
                continuous = d$x1 + stats::rnorm(n))

  effect <- 1 + 0.8 * d$x2
  d$tau <- effect
  numeric_z <- if (is.factor(d$z)) as.numeric(d$z) - 1 else d$z
  d$y <- 2 * d$x1 - d$x2 + numeric_z * effect + stats::rnorm(n, sd = 0.5)
  d
}

bcf_args <- function(...) {
  list(num_burn = 200, num_draws = 200, verbose = FALSE, ...)
}

test_that("bcf sets up the model the way it says it does", {
  d <- sim_causal()

  fit <- do.call(bcf, c(list(y ~ x1 + x2, treat = ~ z, data = d,
                             family = gaussian()), bcf_args()))

  expect_named(fit[["eta"]], c("(Intercept)", "z"))

  # Fewer trees for the effect than the control function, since effect
  # heterogeneity is usually simpler than a prognostic surface. The sparsity
  # prior is left at its default: the treatment is the coefficient rather than a
  # predictor the forest splits on, so nothing can drop it.
  expect_identical(fit[["num_trees"]], c(50L, 25L))
  expect_true(fit[["control"]][["sparsity"]])

  # The propensity score is a predictor of the control function.
  expect_true(".propensity" %in% attr(stats::terms(fit), "term.labels"))
})

test_that("bcf's tree counts give way to either spelling of the caller's", {
  d <- sim_causal(seed = 11)

  trees <- function(...) {
    fit <- do.call(bcf, c(list(y ~ x1 + x2, treat = ~ z, data = d,
                               family = gaussian()), list(...)))
    fit[["num_trees"]]
  }

  # `bartisan()` documents `num_trees = 7` and
  # `control = bartisan_control(num_trees = 7)` as the same thing, so the
  # asymmetric default has to stand aside for both. It once read only `...`, and
  # the control spelling lost to it without a word.
  expect_identical(trees(num_trees = 7L, num_burn = 20L, num_draws = 20L,
                         verbose = FALSE),
                   c(7L, 7L))
  expect_identical(trees(control = bartisan_control(num_trees = 7L,
                                                    num_burn = 20L,
                                                    num_draws = 20L,
                                                    verbose = FALSE)),
                   c(7L, 7L))

  # And still applies when the caller said nothing either way.
  expect_identical(trees(control = bartisan_control(num_burn = 20L,
                                                    num_draws = 20L,
                                                    verbose = FALSE)),
                   c(50L, 25L))
})

test_that("the propensity score reaches the control function and not the effect", {
  d <- sim_causal(seed = 2)

  fit <- do.call(bcf, c(list(y ~ x1 + x2, treat = ~ z, data = d,
                             family = gaussian()), bcf_args()))

  # This is the whole point of estimating it: the control function absorbs the
  # selection without the effect being allowed to vary with it.
  counts <- fit[["counts"]]
  expect_gt(mean(counts[["(Intercept)"]][, ".propensity"]), 0)
  expect_identical(unname(mean(counts[["z"]][, ".propensity"])), 0)
})

test_that("bcf recovers a treatment effect under confounding", {
  d <- sim_causal(seed = 3)

  fit <- do.call(bcf, c(list(y ~ x1 + x2, treat = ~ z, data = d,
                             family = gaussian()), bcf_args()))

  expect_gt(cor(coef(fit)[, "z"], d$tau), 0.85)
  expect_equal(mean(coef(fit)[, "z"]), mean(d$tau), tolerance = 0.25)
})

test_that("a treatment among the covariates is removed rather than refused", {
  d <- sim_causal(seed = 4)

  # `y ~ .` is how most callers will write the covariates, and the treatment is
  # not one of them.
  expect_no_error(
    fit <- do.call(bcf, c(list(y ~ x1 + x2 + z, treat = ~ z, data = d,
                               family = gaussian()), bcf_args())))

  expect_named(fit[["eta"]], c("(Intercept)", "z"))
})

test_that("a categorical treatment gets a forest per level and a score per level", {
  d <- sim_causal(seed = 5, kind = "categorical")

  fit <- do.call(bcf, c(list(y ~ x1 + x2, treat = ~ z, data = d,
                             family = gaussian()), bcf_args()))

  # The balancing score for a multi-valued treatment is the whole vector of
  # assignment probabilities, so all of them go in.
  expect_named(fit[["eta"]], c("(Intercept)", "za", "zb", "zc"))
  expect_true(all(c(".propensity1", ".propensity2", ".propensity3") %in%
                    attr(stats::terms(fit), "term.labels")))
})

test_that("a continuous treatment refuses a propensity score, with the reason", {
  d <- sim_causal(seed = 6, kind = "continuous")

  expect_error(do.call(bcf, c(list(y ~ x1 + x2, treat = ~ z, data = d,
                                   family = gaussian()), bcf_args())),
               "no propensity score that is a probability")

  # Without one it fits, and a score supplied as a number is taken as given.
  expect_no_error(do.call(bcf, c(list(y ~ x1 + x2, treat = ~ z, data = d,
                                      family = gaussian(), propensity = FALSE),
                                 bcf_args())))
  expect_no_error(do.call(bcf, c(list(y ~ x1 + x2, treat = ~ z, data = d,
                                      family = gaussian(),
                                      propensity = stats::runif(nrow(d))),
                                 bcf_args())))
})

test_that("moderators restrict what the effect may vary with", {
  d <- sim_causal(seed = 7)

  fit <- do.call(bcf, c(list(y ~ x1 + x2, treat = ~ z, data = d,
                             family = gaussian(), moderators = ~ x2),
                        bcf_args()))

  counts <- fit[["counts"]][["z"]]
  expect_identical(unname(mean(counts[, "x1"])), 0)
  expect_gt(mean(counts[, "x2"]), 0)
})

test_that("predict works on the caller's own data", {
  d <- sim_causal(seed = 9)

  fit <- do.call(bcf, c(list(y ~ x1 + x2, treat = ~ z, data = d,
                             family = gaussian()), bcf_args()))

  # The fit uses a propensity score the caller never named, so their own frame
  # does not carry it. Rebuilding it from the model the fit kept is what makes
  # this work; without it every prediction failed on a missing variable whose
  # name they never chose.
  expect_no_error(p <- predict(fit, newdata = d))

  # Close, but not to the last bit. The score is a predictor, and a predictor
  # goes through the quantile transform, which is a step function -- so a score
  # that rebuilds to 1e-11 can still land on the other side of a step. Supplying
  # the score instead of rebuilding it removes the reconstruction and the
  # predictions match exactly, which is what pins the cause.
  expect_lt(max(abs(as.numeric(p) - as.numeric(fitted(fit)))),
            0.1 * stats::sd(d$y))

  supplied_score <- cbind(d, fit[["bcf"]][["propensity"]])
  expect_equal(as.numeric(predict(fit, newdata = supplied_score)),
               as.numeric(fitted(fit)), tolerance = 1e-8)

  # And a supplied score cannot be rebuilt, so that says so rather than failing
  # on the missing column.
  supplied <- do.call(bcf, c(list(y ~ x1 + x2, treat = ~ z, data = d,
                                  family = gaussian(),
                                  propensity = stats::runif(nrow(d))),
                             bcf_args()))

  expect_error(predict(supplied, newdata = d), "cannot be rebuilt")
})

test_that("the coefficient is the contrast on the link scale", {
  d <- sim_causal(seed = 10)

  fit <- do.call(bcf, c(list(y ~ x1 + x2, treat = ~ z, data = d,
                             family = gaussian()), bcf_args()))

  treated <- predict(fit, newdata = transform(d, z = 1), type = "link",
                     summary = FALSE, draws = TRUE)
  control <- predict(fit, newdata = transform(d, z = 0), type = "link",
                     summary = FALSE, draws = TRUE)

  expect_equal(mean(rowMeans(treated - control)),
               mean(coef(fit)[, "z"]), tolerance = 1e-10)
})

test_that("bcf validates its treatment argument", {
  d <- sim_causal(seed = 8)

  expect_error(bcf(y ~ x1 + x2, treat = ~ z + x1, data = d),
               "exactly one variable")
  expect_error(bcf(y ~ x1 + x2, treat = ~ nope, data = d),
               "no column")
  expect_error(bcf(y ~ z, treat = ~ z, data = d),
               "at least one covariate")
})

# `bcf()` chooses a drawn coding for its treatment-effect forest, and the engine
# refuses one for a family whose leaf target is not quadratic. The refusal is
# expected and handled -- the fixed coding is used instead -- but it is raised
# inside the sampler, so with more than one chain *future.apply* announces the
# failing attempt with an immediate warning of its own before re-raising it. A
# caller who asked for a Tweedie BCF fit got "Caught Rcpp::exception. Canceling
# all iterations ..." from a fit that then completed correctly.
test_that("a refused drawn coding does not leak the retry's noise", {
  skip_on_cran()
  skip_if_not_installed("future.apply")

  d <- sim_x(n = 200L, p = 4L, seed = 41L)
  d$z <- stats::rbinom(nrow(d), 1L, 0.4)
  # Non-negative with a point mass at zero, which is what a Tweedie is for.
  d$y <- ifelse(stats::runif(nrow(d)) < 0.25, 0,
                stats::rgamma(nrow(d), shape = 2, rate = 1 / (200 * (1 + d$z))))

  ctrl <- quick_control(num_burn = 10L, num_draws = 10L)
  form <- y ~ x1 + x2 + x3 + x4

  # Several chains is the case that warned; the refusal itself happens whatever
  # the chain count.
  expect_no_warning(
    fit <- bcf(form, treat = ~ z, data = d, family = tweedie(),
               chains = 4L, control = ctrl))

  # And the fit it fell back to is a real one, on the fixed coding.
  expect_s3_class(fit, "bartisan_fit")
  expect_length(fit[["counts"]], 2L)

  # The holding-back is conditional on the attempt having failed. A genuine
  # warning still reaches the caller, on the path where the trial is refused
  # and on the path where it succeeds.
  d$y[1:3] <- NA

  expect_warning(bcf(form, treat = ~ z, data = d, family = tweedie(),
                     chains = 4L, control = ctrl),
                 "missing response")
  expect_warning(bcf(form, treat = ~ z, data = d, family = dpm(),
                     chains = 4L, control = ctrl),
                 "missing response")
})
