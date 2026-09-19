# `diagnose()`. The statistics are tested against series whose answer is known
# rather than against a fit, because a fit's convergence is exactly the thing
# they are supposed to be measuring; the fits below test the wiring and the
# advice.

sim_diagnose <- function(n = 200L, seed = 1L) {
  set.seed(seed)
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  d$y <- d$x1 + stats::rnorm(n)
  d
}

diagnose_control <- function(...) {
  bartisan_control(num_trees = 5, num_burn = 100, num_draws = 100,
                   verbose = FALSE, ...)
}

# A stationary autoregressive series: no drift, whatever its autocorrelation.
ar_series <- function(n, rho) {
  as.vector(stats::filter(stats::rnorm(n, sd = sqrt(1 - rho^2)), rho,
                          method = "recursive"))
}

test_that("diagnose() refuses anything that is not a fit", {
  expect_error(diagnose(1), "must be a fit")
  expect_error(diagnose(list()), "must be a fit")
})

test_that("the report has the documented shape", {
  fit <- bartisan(y ~ ., data = sim_diagnose(), family = gaussian(),
                  control = diagnose_control(chains = 2))

  out <- diagnose(fit)

  expect_s3_class(out, "bartisan_diagnosis")
  expect_named(out, c("table", "checks", "advice", "chains", "draws",
                      "rhat_max", "ess_min"))
  expect_true(all(c("quantity", "rhat", "rhat_late", "ess_bulk", "ess_tail") %in%
                    names(out[["table"]])))
  expect_named(out[["checks"]], c("check", "status", "detail"))
  expect_identical(out[["chains"]], 2L)
  expect_identical(out[["draws"]], 200L)
  expect_type(out[["advice"]], "character")

  expect_output(print(out), "Convergence and mixing")
})

test_that("the leaf scale is left out and the forest's size is put in", {
  fit <- bartisan(y ~ ., data = sim_diagnose(seed = 2L), family = gaussian(),
                  control = diagnose_control(chains = 2))

  quantities <- diagnose(fit)[["table"]][["quantity"]]

  # Excluded for the reason in `?diagnose`: it mixes badly in every
  # implementation and nothing reported depends on it.
  expect_false(any(startsWith(quantities, "sigma_mu")))

  # And included, because a forest still growing is the clearest sign warmup
  # ended too early and no generic diagnostic can see it.
  expect_true(any(startsWith(quantities, "splits.")))
  expect_true(any(startsWith(quantities, "eta.")))
})

test_that("one chain is called out, since R-hat has nothing to compare with", {
  fit <- bartisan(y ~ ., data = sim_diagnose(seed = 3L), family = gaussian(),
                  control = diagnose_control(chains = 1))

  out <- diagnose(fit)
  chains <- out[["checks"]][out[["checks"]][["check"]] == "chains", ]

  expect_identical(chains[["status"]], "warn")
  expect_match(out[["advice"]][1L], "chains = 4", fixed = TRUE)

  # It is still computed, by splitting the one chain, rather than left blank.
  expect_true(any(is.finite(out[["table"]][["rhat"]])))
})

test_that("the thresholds are the caller's to move", {
  fit <- bartisan(y ~ ., data = sim_diagnose(seed = 4L), family = gaussian(),
                  control = diagnose_control(chains = 2))

  strict <- diagnose(fit, ess_min = 1e6)
  loose <- diagnose(fit, rhat_max = 10, ess_min = 1)

  expect_true(any(strict[["checks"]][["status"]] == "warn"))
  expect_length(loose[["advice"]], 0L)

  expect_error(diagnose(fit, rhat_max = 0.5), "rhat_max")
})

test_that("the late-draw R-hat is the same statistic on the second half", {
  set.seed(11)

  # A chain that starts badly and then settles: R-hat over everything is
  # elevated, and R-hat over the late draws alone is not. That gap is the whole
  # basis for telling "warmup was too short" from "the chains disagree".
  settles <- cbind(c(seq(6, 0, length.out = 200), ar_series(800L, 0.5)),
                   c(seq(-6, 0, length.out = 200), ar_series(800L, 0.5)))

  expect_gt(rhat_rank(settles), 1.05)
  expect_lt(rhat_late(settles), rhat_rank(settles))

  # Two chains that never agree stay bad on the late draws too.
  apart <- cbind(ar_series(1000L, 0.5) + 5, ar_series(1000L, 0.5) - 5)

  expect_gt(rhat_rank(apart), 1.05)
  expect_gt(rhat_late(apart), 1.05)

  # Too short to halve and still report something meaningful.
  expect_true(is.na(rhat_late(matrix(1:6, ncol = 2L))))
})

test_that("a single chain is folded so that R-hat has two halves to compare", {
  x <- matrix(seq_len(100), ncol = 1L)
  folded <- fold_halves(x)

  expect_identical(dim(folded), c(50L, 2L))
  expect_equal(folded[, 1L], 1:50)
  expect_equal(folded[, 2L], 51:100)

  # Too short to halve, so it comes back untouched rather than empty.
  expect_identical(fold_halves(matrix(1:3, ncol = 1L)), matrix(1:3, ncol = 1L))
})

test_that("the advice follows which statistic failed, not merely that one did", {
  checks <- function(...) {
    data.frame(check = c(...), status = "warn", detail = "")
  }

  # A short warmup gets the one fix, and not also a recommendation to collect
  # more draws from a distribution the sampler has not reached.
  warmup <- diagnosis_advice(checks("warmup"))
  expect_length(warmup, 1L)
  expect_match(warmup, "num_burn", fixed = TRUE)

  # R-hat elevated with warmup not implicated is the other case: chains that
  # have each settled somewhere different.
  disagree <- diagnosis_advice(checks("rhat"))
  expect_match(disagree[1L], "settled")
  expect_true(any(grepl("num_trees", disagree, fixed = TRUE)))

  # Both at once is warmup, and the disagreement advice is withheld, since more
  # burn-in is what to try first.
  both <- diagnosis_advice(checks("rhat", "warmup"))
  expect_false(any(grepl("num_trees", both, fixed = TRUE)))

  # Low effective sample size on its own is the benign case, and thinning is
  # named as the thing not to reach for.
  thin <- diagnosis_advice(checks("bulk ESS"))
  expect_match(thin[1L], "num_thin", fixed = TRUE)

  expect_length(diagnosis_advice(data.frame(check = "rhat", status = "ok",
                                            detail = "")), 0L)

  # The forest's own size takes advice of its own, which names what it does bind
  # on rather than sending the reader up the ladder above.
  forest <- diagnosis_advice(checks("forest size"))
  expect_length(forest, 1L)
  expect_match(forest, "variable_importance", fixed = TRUE)
  expect_false(any(grepl("num_trees", forest, fixed = TRUE)))
})

# A sum of trees reaches one function through many partitions, so the number of
# splitting rules is not pinned down the way a fitted value is. The checks on
# R-hat and effective sample size read the reported quantities, and that row
# gets its own check -- separated rather than suppressed, since it still binds on
# anything computed from the split counts.
test_that("the internal states are graded apart from the reported quantities", {
  row <- function(quantity, rhat, ess, bad) {
    data.frame(quantity = quantity, rhat = rhat, rhat_late = rhat,
               ess_bulk = ess, ess_tail = ess, ess_frac = 0.5,
               rhat_bad = bad, late_bad = bad)
  }

  # Everything reported is fine; only the forest disagrees.
  table <- rbind(row("loglik", 1.002, 900, 0),
                 row("splits.eta", 1.900, 6, 1),
                 row("eta.eta (average over observations)", 1.001, 950, 0))

  checks <- diagnosis_checks(table, chains = 4L, draws = 4000L,
                             rhat_max = 1.01, ess_min = 400)

  status <- function(name) checks[["status"]][checks[["check"]] == name]

  expect_identical(status("rhat"), "ok")
  expect_identical(status("bulk ESS"), "ok")
  expect_identical(status("tail ESS"), "ok")
  expect_identical(status("forest size"), "warn")

  # The reported rows are what the R-hat and ESS lines quote, so the nuisance
  # row's 1.900 and its 6 effective draws appear on the forest line and nowhere
  # else.
  expect_false(any(grepl("splits.eta",
                         checks[["detail"]][checks[["check"]] != "forest size"],
                         fixed = TRUE)))

  # And a reported quantity failing is still caught with the forest fine, so the
  # separation does not swallow the case it is not for.
  other <- rbind(row("loglik", 1.500, 8, 1),
                 row("splits.eta", 1.002, 900, 0))

  expect_identical(
    diagnosis_checks(other, chains = 4L, draws = 4000L, rhat_max = 1.01,
                     ess_min = 400) |>
      subset(check == "rhat", "status", drop = TRUE),
    "warn")
})

# An estimand is a contrast, and a contrast can mix badly where the function it
# is a contrast of mixes well. Under the default splitting prior a draw that
# gives the treatment no rule puts the contrast at exactly zero, and the sampler
# can stay there for a long run while the other predictors keep the fitted
# function moving, so the fit's own table shows nothing wrong.
test_that("diagnose() reports on the estimand, not only on the fit", {
  skip_on_cran()

  d <- sim_x(n = 250L, p = 3L, seed = 91L)
  d$z <- stats::rbinom(nrow(d), 1L, 0.5)
  d$y <- stats::rbinom(nrow(d), 1L, stats::plogis(d$x1 + 0.7 * d$z))

  fit <- bartisan(y ~ ., d, family = stats::binomial(), chains = 2L,
                  control = quick_control(num_trees = 10L, num_burn = 100L,
                                          num_draws = 200L))

  eff <- estimate_effect(fit, treat = "z")
  dg <- diagnose(eff)

  expect_s3_class(dg, "bartisan_diagnosis")
  expect_identical(dg[["chains"]], 2L)
  expect_identical(dg[["draws"]], 400L)

  # The contrast, and under it the two averages it is a contrast of: a contrast
  # that mixes badly can have one of them to blame rather than both.
  expect_identical(dg[["table"]][["quantity"]],
                   c(names(attr(eff, "draws")), names(attr(eff, "po_draws"))))
  expect_identical(nrow(dg[["table"]]), 3L)
  expect_true(all(is.finite(dg[["table"]][["ess_bulk"]])))
  expect_true(all(is.finite(dg[["table"]][["rhat"]])))

  fit_diag <- diagnose(fit)
  expect_false(identical(dg[["table"]][["ess_bulk"]],
                         fit_diag[["table"]][["ess_bulk"]]))

  # `by` gives one row per reported group, plus the potential outcomes, and
  # `CATE` folds the units the way the fit's table folds observations rather
  # than printing one row each.
  by_rows <- diagnose(estimate_effect(fit, treat = "z", by = ~ x3 > 0.5))
  expect_identical(nrow(by_rows[["table"]]), 4L)

  cate <- diagnose(estimate_effect(fit, treat = "z", estimand = "CATE"))
  expect_match(cate[["table"]][["quantity"]][1L], "worst 5% of 250 units",
               fixed = TRUE)

  # Every check reads as a sentence that names what it is about, rather than as
  # a fragment continuing a label the reader cannot see.
  expect_match(paste(dg[["checks"]][["detail"]], collapse = " "), "R-hat is")
  expect_match(paste(dg[["checks"]][["detail"]], collapse = " "), "ESS is")

  # The atom is the failure the fit cannot show, so it is named when present.
  stuck <- mean(attr(eff, "draws")[[1L]] == 0)

  if (stuck > 0.01) {
    expect_true("atom" %in% dg[["checks"]][["check"]])
    expect_match(paste(dg[["advice"]], collapse = " "), "atom at zero")
  }

  expect_error(diagnose(1), "must be a fit")
  expect_error(diagnose(list()), "must be a fit")
})

# The convergence pass splits its columns over a `future` plan. Written as a
# bare call, `diagnosis_block()` is read by future as belonging to this package
# and dropped from the globals it ships, leaving the worker to find an
# unexported function on a search path that carries only exports. Installed that
# resolved anyway; under `pkgload::load_all()`, which is how the package is run
# while being worked on, `diagnose()` died with "could not find function
# diagnosis_block" as soon as a plan had more than one worker and the pass was
# wide enough to split.
test_that("the parallel convergence pass ships what the worker needs", {
  skip_on_cran()
  skip_if_not_installed("future")

  # Two real workers, or there is nothing to test; a plan that cannot be set up
  # here (a single core, a sandbox) is a skip rather than a failure.
  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)

  started <- tryCatch({
    future::plan(future::multisession, workers = 2L)
    isTRUE(future::nbrOfWorkers() >= 2L)
  }, error = function(e) FALSE)

  skip_if_not(started, "no second worker available")

  d <- sim_x(n = 200L, p = 3L, seed = 95L)
  d$z <- stats::rbinom(nrow(d), 1L, 0.5)
  d$y <- stats::rbinom(nrow(d), 1L, stats::plogis(d$x1 + 0.7 * d$z))

  fit <- bartisan(y ~ ., d, family = stats::binomial(), chains = 2L,
                  control = quick_control(num_trees = 5L, num_burn = 60L,
                                          num_draws = 100L))

  # 200 observations is over the 100-column floor, so the pass really splits.
  parallel <- diagnose(fit)

  expect_s3_class(parallel, "bartisan_diagnosis")

  # Splitting is an implementation detail and must not change the answer.
  future::plan("sequential")
  sequential <- diagnose(fit)

  expect_equal(parallel[["table"]], sequential[["table"]])

  # And the same for an estimand, whose `CATE` path takes the same branch.
  future::plan(future::multisession, workers = 2L)
  cate <- estimate_effect(fit, treat = "z", estimand = "CATE")

  expect_s3_class(diagnose(cate), "bartisan_diagnosis")
})
