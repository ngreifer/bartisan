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
})
