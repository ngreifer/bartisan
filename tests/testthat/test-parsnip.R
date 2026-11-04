# The "bartisan" engine for `parsnip::bart()`. parsnip reaches the package only
# through the registry `register_parsnip_engine()` fills, so these check the
# registry entries by what comes out of them: the translated arguments, and
# predictions of every type against the same quantities computed from the
# underlying fit directly.

sim_parsnip <- function(n = 80, seed = 1) {
  d <- sim_x(n, 3, seed)
  d$y <- 2 * d$x1 - d$x2 + stats::rnorm(n, sd = 0.5)
  d$g <- factor(ifelse(d$y > mean(d$y), "high", "low"))
  d$k <- cut(d$y, 3, labels = c("a", "b", "c"))

  log_t <- 1 + d$x1 - 0.5 * d$x2 + 0.5 * stats::rnorm(n)
  cens <- log(stats::runif(n, 1, 10))
  d$time <- exp(pmin(log_t, cens))
  d$status <- as.numeric(log_t <= cens)
  d
}

test_that("the engine is registered for all three modes", {
  skip_if_not_installed("parsnip")

  engines <- parsnip::show_engines("bart")
  expect_setequal(engines$mode[engines$engine == "bartisan"],
                  c("regression", "classification", "censored regression"))

  # Registering a second time, as happens when parsnip is reloaded, leaves one
  # entry per mode.
  register_parsnip_engine()
  engines <- parsnip::show_engines("bart")
  expect_identical(sum(engines$engine == "bartisan"), 3L)
})

test_that("the model arguments reach `bartisan_control()`", {
  skip_if_not_installed("parsnip")

  d <- sim_parsnip()

  fit <- parsnip::bart(trees = 7, prior_terminal_node_coef = 0.9,
                       prior_terminal_node_expo = 1.5,
                       prior_outcome_range = 3, mode = "regression") |>
    parsnip::set_engine("bartisan", family = gaussian(), num_burn = 30,
                        num_draws = 30) |>
    parsnip::fit(y ~ x1 + x2 + x3, data = d)

  engine <- parsnip::extract_fit_engine(fit)
  expect_s3_class(engine, "bartisan_fit")

  control <- engine[["control"]]
  expect_equal(control[["num_trees"]], 7)
  expect_equal(control[["gamma"]], 0.9)
  expect_equal(control[["beta"]], 1.5)
  expect_equal(control[["k"]], 3)
  expect_equal(control[["num_draws"]], 30)
})

test_that("regression predictions match those of the underlying fit", {
  skip_if_not_installed("parsnip")

  d <- sim_parsnip()
  nd <- d[1:6, ]

  fit <- parsnip::bart(trees = 5, mode = "regression") |>
    parsnip::set_engine("bartisan", family = gaussian(), num_burn = 30,
                        num_draws = 40) |>
    parsnip::fit(y ~ x1 + x2 + x3, data = d)
  engine <- parsnip::extract_fit_engine(fit)

  expect_equal(stats::predict(fit, nd)$.pred,
               stats::predict(engine, newdata = nd), tolerance = 1e-10)

  raw <- stats::predict(fit, nd, type = "raw")
  expect_identical(dim(raw), c(40L, 6L))

  # Equal-tailed intervals: the quantiles of the posterior draws of the mean.
  ci <- stats::predict(fit, nd, type = "conf_int", level = 0.9,
                       std_error = TRUE)
  expect_named(ci, c(".pred_lower", ".pred_upper", ".std_error"))
  expect_equal(ci$.pred_lower, apply(raw, 2L, stats::quantile, 0.05,
                                     names = FALSE), tolerance = 1e-10)
  expect_equal(ci$.pred_upper, apply(raw, 2L, stats::quantile, 0.95,
                                     names = FALSE), tolerance = 1e-10)
  expect_equal(ci$.std_error, apply(raw, 2L, stats::sd), tolerance = 1e-10)

  skip_if_not_installed("rstantools")

  # A predictive interval adds the residual variation, so it is wider.
  pi <- stats::predict(fit, nd, type = "pred_int", level = 0.9)
  expect_named(pi, c(".pred_lower", ".pred_upper"))
  expect_true(all(pi$.pred_upper - pi$.pred_lower >
                    ci$.pred_upper - ci$.pred_lower))
})

patrick::with_parameters_test_that(
  "classification predictions match those of the underlying fit:",
  {
    skip_if_not_installed("parsnip")

    d <- sim_parsnip()
    nd <- d[1:6, ]
    levels <- levels(d[[outcome]])

    fit <- parsnip::bart(trees = 5, mode = "classification") |>
      parsnip::set_engine("bartisan", family = family, num_burn = 30,
                          num_draws = 40) |>
      parsnip::fit(stats::reformulate(c("x1", "x2", "x3"), outcome), data = d)
    engine <- parsnip::extract_fit_engine(fit)

    prob <- stats::predict(fit, nd, type = "prob")
    expect_named(prob, paste0(".pred_", levels))
    expect_equal(as.matrix(prob),
                 stats::predict(engine, newdata = nd, type = "prob"),
                 tolerance = 1e-10, ignore_attr = TRUE)

    class <- stats::predict(fit, nd, type = "class")$.pred_class
    expect_identical(levels(class), levels)
    expect_identical(as.character(class),
                     levels[max.col(as.matrix(prob), ties.method = "first")])

    raw <- stats::predict(fit, nd, type = "raw")
    expect_identical(dim(raw), c(40L, 6L, length(levels)))

    ci <- stats::predict(fit, nd, type = "conf_int", level = 0.8)
    expect_named(ci, c(paste0(".pred_lower_", levels),
                       paste0(".pred_upper_", levels)))

    for (j in seq_along(levels)) {
      expect_equal(ci[[paste0(".pred_lower_", levels[j])]],
                   apply(raw[, , j], 2L, stats::quantile, 0.1, names = FALSE),
                   tolerance = 1e-10)
      expect_equal(ci[[paste0(".pred_upper_", levels[j])]],
                   apply(raw[, , j], 2L, stats::quantile, 0.9, names = FALSE),
                   tolerance = 1e-10)
    }
  },
  patrick::cases(
    binary = list(outcome = "g", family = binomial()),
    multinomial = list(outcome = "k", family = multinomial())
  )
)

patrick::with_parameters_test_that(
  "censored regression predictions match those of the underlying fit:",
  {
    skip_if_not_installed("parsnip")
    skip_if_not_installed("survival")

    d <- sim_parsnip()
    nd <- d[1:6, ]
    times <- c(0, 1, 3)

    fit <- parsnip::bart(trees = 5, mode = "censored regression") |>
      parsnip::set_engine("bartisan", family = family, num_burn = 30,
                          num_draws = 40) |>
      parsnip::fit(survival::Surv(time, status) ~ x1 + x2 + x3, data = d)
    engine <- parsnip::extract_fit_engine(fit)
    expect_identical(engine[["family"]][["family"]],
                     sub("^(weibull|lognormal)_aft$", "aft", family))

    expect_equal(stats::predict(fit, nd, type = "time")$.pred_time,
                 stats::predict(engine, newdata = nd), tolerance = 1e-10)

    raw <- stats::predict(fit, nd, type = "raw")
    expect_identical(dim(raw), c(40L, 6L))

    # One tibble of evaluation times per row, with the survival function at
    # each and, since parsnip allows a time of 0, a survival of 1 there.
    surv <- stats::predict(fit, nd, type = "survival", eval_time = times,
                           interval = "confidence", level = 0.8)
    expect_named(surv, ".pred")
    expect_length(surv$.pred, 6L)
    expect_named(surv$.pred[[1L]], c(".eval_time", ".pred_survival",
                                     ".pred_lower", ".pred_upper"))

    curve <- do.call(rbind, lapply(surv$.pred, function(p) p$.pred_survival))
    lower <- do.call(rbind, lapply(surv$.pred, function(p) p$.pred_lower))
    expect_equal(curve[, 1L], rep(1, 6L))
    expect_equal(lower[, 1L], rep(1, 6L))

    draws <- stats::predict(engine, newdata = nd, type = "survival",
                            times = times[-1L], draws = TRUE)
    expect_equal(curve[, -1L], apply(draws, c(2L, 3L), mean),
                 tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(lower[, -1L],
                 apply(draws, c(2L, 3L), stats::quantile, 0.1, names = FALSE),
                 tolerance = 1e-10, ignore_attr = TRUE)

    # The hazard has no value reported at 0, which for most families is not
    # finite.
    haz <- stats::predict(fit, nd, type = "hazard", eval_time = times)
    expect_named(haz$.pred[[1L]], c(".eval_time", ".pred_hazard"))

    rate <- do.call(rbind, lapply(haz$.pred, function(p) p$.pred_hazard))
    expect_true(all(is.na(rate[, 1L])))
    expect_equal(rate[, -1L],
                 stats::predict(engine, newdata = nd, type = "hazard",
                                times = times[-1L]),
                 tolerance = 1e-10, ignore_attr = TRUE)

    # The linear predictor increases with the survival time, which for `ph()`
    # means the negative of its log hazard, and `increasing = FALSE` reverses
    # it.
    link <- stats::predict(engine, newdata = nd, type = "link")
    if (is.matrix(link)) {
      link <- link[, 1L]
    }
    sign <- {
      if (identical(family, "ph")) -1
      else 1
    }

    expect_equal(stats::predict(fit, nd, type = "linear_pred")$.pred_linear_pred,
                 sign * unname(link), tolerance = 1e-10)
    expect_equal(stats::predict(fit, nd, type = "linear_pred",
                                increasing = FALSE)$.pred_linear_pred,
                 -sign * unname(link), tolerance = 1e-10)
  },
  patrick::cases(
    weibull_aft = list(family = "weibull_aft"),
    ph = list(family = "ph"),
    dpm_aft = list(family = "dpm_aft"),
    lognormal_aft_ls = list(family = "lognormal_aft_ls")
  )
)

test_that("survival predictions carry the censoring weights tidymodels adds", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("survival")

  d <- sim_parsnip()
  d$surv <- survival::Surv(d$time, d$status)

  fit <- parsnip::bart(trees = 5, mode = "censored regression") |>
    parsnip::set_engine("bartisan", family = "lognormal_aft", num_burn = 30,
                        num_draws = 30) |>
    parsnip::fit(surv ~ x1 + x2 + x3, data = d)

  # `augment()` asks for the weights the survival metrics in yardstick need,
  # which parsnip computes from the `.pred` tibbles, so this checks their
  # layout as much as their contents.
  out <- parsnip::augment(fit, d[1:10, ], eval_time = c(1, 3))
  expect_true(all(c(".pred", ".pred_time") %in% names(out)))
  expect_true(all(c(".eval_time", ".pred_survival", ".weight_censored") %in%
                    names(out$.pred[[1L]])))
})

test_that("intervals can be computed for a single row of new data", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("survival")

  d <- sim_parsnip()

  fit <- parsnip::bart(trees = 5, mode = "classification") |>
    parsnip::set_engine("bartisan", family = "multinomial", num_burn = 30,
                        num_draws = 30) |>
    parsnip::fit(k ~ x1 + x2 + x3, data = d)

  ci <- stats::predict(fit, d[1L, ], type = "conf_int")
  expect_identical(dim(ci), c(1L, 6L))
  expect_true(all(ci[1:3] <= ci[4:6]))

  fit <- parsnip::bart(trees = 5, mode = "censored regression") |>
    parsnip::set_engine("bartisan", family = "lognormal_aft", num_burn = 30,
                        num_draws = 30) |>
    parsnip::fit(survival::Surv(time, status) ~ x1 + x2 + x3, data = d)

  surv <- stats::predict(fit, d[1L, ], type = "survival", eval_time = c(1, 3),
                         interval = "confidence")
  expect_length(surv$.pred, 1L)
  expect_identical(dim(surv$.pred[[1L]]), c(2L, 4L))
})
