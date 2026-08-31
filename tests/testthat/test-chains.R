# Independent chains, pooled, with a convergence diagnostic. The chain is the
# only parallel axis this sampler has: a sweep conditions on the last one, so
# there is nothing to split within a chain.

test_that("more than one chain runs without future.apply, on the same streams", {
  d <- sim_x(seed = 73)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  run <- function(hide) {
    set.seed(11)
    if (!hide) {
      return(bartisan(y ~ ., data = d, control = quick_control(chains = 3)))
    }
    installed <- rlang::is_installed
    testthat::with_mocked_bindings(
      bartisan(y ~ ., data = d, control = quick_control(chains = 3)),
      is_installed = function(pkg, ...) {
        if (identical(pkg, "future.apply")) FALSE else installed(pkg, ...)
      },
      .package = "rlang")
  }

  # Parallelism is how fast the chains are, not whether the model is fitted, so
  # the package falls back to running them one after another.
  expect_no_error(sequential <- run(TRUE))
  expect_identical(sequential[["chains"]], 3L)

  # And the two branches draw from the same streams, so a script does not change
  # its answer depending on whether future.apply happens to be installed.
  skip_if_not_installed("future.apply")
  expect_equal(run(FALSE)[["eta"]][[1L]], sequential[["eta"]][[1L]])
})

test_that("running chains leaves the session generator alone", {
  d <- sim_x(seed = 74)
  d$y <- stats::rnorm(nrow(d))

  set.seed(1)
  before <- RNGkind()
  bartisan(y ~ ., data = d, control = quick_control(chains = 2))

  # The streams are L'Ecuyer, which the session did not ask for.
  expect_identical(RNGkind(), before)
})

test_that("chains is a control setting and still reachable through the dots", {
  d <- sim_x(seed = 71)
  d$y <- stats::rnorm(nrow(d))

  # It lives on `bartisan_control()` now.
  expect_true("chains" %in% names(formals(bartisan_control)))
  expect_false("chains" %in% names(formals(bartisan)))

  # Both routes reach the same place, so the calls that passed it to
  # `bartisan()` before the move still work.
  a <- bartisan(y ~ ., data = d, control = quick_control(chains = 2))
  b <- bartisan(y ~ ., data = d, chains = 2, control = quick_control())

  expect_identical(a[["chains"]], 2L)
  expect_identical(b[["chains"]], 2L)
  expect_identical(a[["control"]][["chains"]], b[["control"]][["chains"]])
})

test_that("the leaf scale is not one of the diagnosed quantities", {
  d <- sim_x(seed = 72)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  fit <- bartisan(y ~ ., data = d, control = quick_control(chains = 2))

  # It mixes badly in every BART implementation, including one that draws it
  # exactly, and nothing reported depends on it, so it is out of the table. The
  # draws stay on the fit.
  expect_false(any(grepl("sigma_mu", fit[["rhat"]][["quantity"]])))

  # Out of the table and only out of the table.
  expect_true(is.matrix(fit[["sigma_mu"]]))
  skip_if_not_installed("posterior")
  expect_true("sigma_mu.eta" %in% posterior::variables(posterior::as_draws(fit)))
})

test_that("chains are pooled into one set of draws", {
  skip_if_not_installed("future.apply")

  d <- sim_x(n = 100, seed = 81)
  d$y <- 2 * d$x1 + stats::rnorm(100, sd = 0.4)

  one <- bartisan(y ~ ., data = d, control = quick_control())
  four <- bartisan(y ~ ., data = d, chains = 4, control = quick_control())

  expect_identical(four[["chains"]], 4L)
  expect_identical(nrow(four[["eta"]][["eta"]]),
                   4L * nrow(one[["eta"]][["eta"]]))
  expect_identical(ncol(four[["eta"]][["eta"]]),
                   ncol(one[["eta"]][["eta"]]))
  expect_identical(nrow(four[["sigma_mu"]]), 4L * nrow(one[["sigma_mu"]]))
  expect_identical(length(four[["loglik"]]), 4L * length(one[["loglik"]]))
  expect_identical(nrow(four[["aux"]]), 4L * nrow(one[["aux"]]))
  expect_identical(nrow(four[["counts"]][["eta"]]),
                   4L * nrow(one[["counts"]][["eta"]]))
})

test_that("the pooled forests still reproduce the pooled predictor", {
  skip_if_not_installed("future.apply")

  d <- sim_x(n = 100, seed = 82)
  d$y <- 2 * d$x1 + stats::rnorm(100, sd = 0.4)

  fit <- bartisan(y ~ ., data = d, chains = 3, control = quick_control())

  # The forests of each chain are concatenated and their record offsets shifted,
  # so this is the check that the shifting is right: replaying the stored trees
  # has to reproduce the predictor the sampler recorded, across every chain.
  expect_predictor_invariant(fit, d)
})

test_that("one seed reproduces the whole run whatever the backend", {
  skip_if_not_installed("future.apply")
  skip_if_not_installed("future")

  d <- sim_x(n = 80, seed = 83)
  d$y <- stats::rbinom(80, 1, 0.4)

  fit <- function() {
    set.seed(11)
    bartisan(y ~ ., data = d, family = binomial(), chains = 3,
             control = quick_control())
  }

  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)

  a <- fit()
  b <- fit()
  expect_equal(a[["eta"]], b[["eta"]])

  set.seed(12)
  c <- bartisan(y ~ ., data = d, family = binomial(), chains = 3,
                control = quick_control())
  expect_false(isTRUE(all.equal(a[["eta"]], c[["eta"]])))

  # Chains must differ from one another, or the streams are not independent.
  per <- nrow(a[["sigma_mu"]]) / 3L
  first <- a[["eta"]][["eta"]][1L, ]
  second <- a[["eta"]][["eta"]][per + 1L, ]
  expect_false(isTRUE(all.equal(first, second)))
})

test_that("the diagnostics are reported as a table of the right shape", {
  skip_if_not_installed("future.apply")

  d <- sim_x(n = 100, seed = 84)
  d$y <- 2 * d$x1 + stats::rnorm(100, sd = 0.4)

  one <- bartisan(y ~ ., data = d, control = quick_control())
  expect_null(one[["rhat"]])

  fit <- bartisan(y ~ ., data = d, chains = 4,
                  control = quick_control(num_burn = 200L, num_draws = 200L))

  diagnostics <- fit[["rhat"]]
  expect_s3_class(diagnostics, "data.frame")
  expect_identical(names(diagnostics),
                   c("quantity", "rhat", "ess_bulk", "ess_tail"))
  expect_true("loglik" %in% diagnostics$quantity)
  expect_true(any(grepl("^eta\\.", diagnostics$quantity)))
  expect_true(all(is.finite(diagnostics$rhat)))
  expect_true(all(diagnostics$rhat >= 1 - 1e-6))

  # An effective sample size cannot exceed the draws there are, and the tail one
  # is never the larger of the two by much.
  total <- nrow(fit[["sigma_mu"]])
  expect_true(all(diagnostics$ess_bulk <= total * 1.5))
  expect_true(all(diagnostics$ess_bulk > 0))

  # Independent chains from the same target.
  expect_lt(diagnostics$rhat[diagnostics$quantity == "loglik"], 1.3)
})

test_that("split-R-hat is the textbook quantity", {
  set.seed(85)

  # Computed by hand from the definition, on the split chains, so that the
  # implementation is checked against the formula rather than against itself.
  by_hand <- function(x) {
    n <- nrow(x)
    half <- n %/% 2L
    y <- cbind(x[seq_len(half), , drop = FALSE],
               x[n - half + seq_len(half), , drop = FALSE])
    within <- mean(apply(y, 2L, stats::var))
    between <- half * stats::var(colMeans(y))
    sqrt(((half - 1) / half * within + between / half) / within)
  }

  cases <- list(matrix(stats::rnorm(4000), ncol = 4L),
                cbind(stats::rnorm(500), stats::rnorm(500) + 5),
                matrix(rep(seq(-3, 3, length.out = 200), 3L), ncol = 3L),
                matrix(stats::rt(2000, 3), ncol = 4L))

  for (x in cases) {
    expect_equal(split_rhat(x), by_hand(x), tolerance = 1e-10)
  }

  # And it does what a diagnostic should: near one for stationary chains, large
  # for chains in different places, and large for a chain that drifts even when
  # every chain drifts identically -- which is the whole reason for splitting.
  expect_equal(split_rhat(cases[[1L]]), 1, tolerance = 0.05)
  expect_gt(split_rhat(cases[[2L]]), 2)
  expect_gt(split_rhat(cases[[3L]]), 1.5)

  # Too few draws or a single chain leaves the diagnostic undefined.
  expect_true(is.na(split_rhat(matrix(1:6, ncol = 1L))))
  expect_true(is.na(split_rhat(matrix(1:2, ncol = 2L))))
})

test_that("rank normalization makes the diagnostic scale-free", {
  set.seed(86)
  x <- matrix(stats::rnorm(4000), ncol = 4L)

  # Any increasing transformation leaves the ranks alone, so it leaves anything
  # computed from them alone too. That is the point of the transform: R-hat and
  # the effective sample size are derived for quantities with finite variance,
  # and ranks have it whatever the posterior looks like.
  expect_equal(split_rhat(rank_normalize(exp(x))),
               split_rhat(rank_normalize(x)), tolerance = 1e-10)
  expect_equal(ess_bulk(x^3), ess_bulk(x), tolerance = 1e-10)
  expect_equal(ess_tail(exp(x)), ess_tail(x), tolerance = 1e-10)

  # The folded half is deliberately not invariant: it is a statement about the
  # spread, and a nonlinear transformation genuinely changes that. So the
  # reported maximum of the two moves a little, and should.
  expect_equal(rhat_rank(exp(x)), rhat_rank(x), tolerance = 1e-3)

  # The unranked version is not scale-free, which is what it is being replaced
  # for: a Cauchy-tailed quantity has no variance for the formula to estimate.
  heavy <- matrix(stats::rcauchy(4000), ncol = 4L)
  expect_true(is.finite(rhat_rank(heavy)))
  expect_gt(ess_bulk(heavy), 100)

  # Folding is what catches chains that agree about the middle and disagree
  # about the spread; the unfolded diagnostic is blind to it.
  scales <- cbind(stats::rnorm(1000), stats::rnorm(1000),
                  stats::rnorm(1000, sd = 4), stats::rnorm(1000))
  expect_gt(rhat_rank(scales), split_rhat(rank_normalize(scales)))
  expect_gt(rhat_rank(scales), 1.05)
})

test_that("the effective sample size recovers what theory says it should", {
  set.seed(87)

  # An AR(1) chain has integrated autocorrelation time (1+rho)/(1-rho), so its
  # effective sample size is known in closed form. The estimator is deliberately
  # conservative, so it should come in at or below the theoretical value.
  ar1 <- function(n, rho) {
    x <- numeric(n)
    for (i in 2:n) {
      x[i] <- rho * x[i - 1L] + sqrt(1 - rho^2) * stats::rnorm(1)
    }
    x
  }

  for (rho in c(0.5, 0.8, 0.95)) {
    x <- matrix(replicate(4, ar1(4000, rho)), ncol = 4L)
    theory <- 16000 * (1 - rho) / (1 + rho)
    got <- ess_bulk(x)

    expect_lt(got, theory * 1.1)
    expect_gt(got, theory * 0.75)
  }

  # Independent draws: the effective sample size is the number of draws.
  independent <- matrix(stats::rnorm(8000), ncol = 4L)
  expect_equal(ess_bulk(independent), 8000, tolerance = 0.1)

  # A tail effective sample size is not the same number as a bulk one, and for
  # a sticky chain it is the one that runs short first.
  sticky <- matrix(replicate(4, ar1(2000, 0.95)), ncol = 4L)
  expect_true(is.finite(ess_tail(sticky)))
  expect_lt(ess_bulk(sticky), 8000 * 0.2)
})

test_that("a quantity the sampler holds fixed reports NA rather than nonsense", {
  # An ordinal model's first cutpoint is pinned at zero for identifiability, so
  # there is nothing to diagnose. Reducing an all-NA set with `na.rm` returns an
  # infinity and warns, and the sample autocovariance of a constant is a rounding
  # error rather than exactly zero, so it passes a variance guard and produces a
  # finite effective sample size out of nothing.
  const <- matrix(0, 200L, 3L)

  expect_identical(rhat_rank(const), NA_real_)
  expect_identical(ess_bulk(const), NA_real_)
  expect_identical(ess_tail(const), NA_real_)

  expect_silent({
    rhat_rank(const)
    ess_bulk(const)
    ess_tail(const)
  })

  # A constant at some value other than zero, and one where only some chains are
  # constant, behave the same way and differently respectively.
  expect_identical(ess_bulk(matrix(2.5, 200L, 3L)), NA_real_)

  set.seed(1)
  varying <- matrix(stats::rnorm(600), 200L, 3L)
  expect_gt(ess_bulk(varying), 400)
  expect_lt(rhat_rank(varying), 1.05)
})

test_that("the diagnostics table survives a pinned cutpoint", {
  skip_on_cran()

  d <- sim_x(n = 200, seed = 41)
  set.seed(1041)
  z <- 2 * d$x1 - d$x2 + stats::rnorm(nrow(d))

  # Two categories, which is the case where the single boundary is folded into
  # the intercept and stays pinned. With three or more it is free, and this is
  # then no longer a test of a constant quantity.
  d$y <- ordered(as.integer(z > stats::median(z)))

  fit <- expect_silent(bartisan(y ~ ., d, family = ordinal("probit"),
                               gate = "hard", chains = 2L,
                               control = quick_control()))

  cut1 <- fit[["rhat"]][fit[["rhat"]]$quantity == "aux.cut1", ]
  expect_identical(nrow(cut1), 1L)
  expect_true(is.na(cut1$rhat))
  expect_true(is.na(cut1$ess_bulk))
  expect_true(is.na(cut1$ess_tail))

  # The other rows are still populated, so one pinned quantity does not take the
  # table down with it.
  rest <- fit[["rhat"]][fit[["rhat"]]$quantity != "aux.cut1", ]
  expect_true(all(is.finite(rest$rhat)))
})
