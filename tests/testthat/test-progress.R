# Progress reporting. The sampler is one C++ call, so what is under test is that
# it calls back into R the promised number of times and that nothing about the
# fit depends on whether anyone was listening.

sim_progress <- function(n = 200L, seed = 1L) {
  set.seed(seed)
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  d$y <- d$x1 + stats::rnorm(n)
  d
}

progress_control <- function(...) {
  bartisan_control(num_trees = 5, num_burn = 100, num_draws = 100,
                   verbose = FALSE, ...)
}

# Every `progression` condition the fit signals, counted. `with_progress()` is
# what installs the relay, and without it a condition signalled in a worker never
# reaches this session.
count_progress <- function(expr) {
  seen <- 0L

  rlang::local_options(progressr.enable = TRUE)

  progressr::with_progress(
    withCallingHandlers(
      force(expr),
      progression = function(cnd) {
        seen <<- seen + 1L

        # A condition relayed from a worker is re-signalled with
        # `signalCondition()`, which establishes no restart, so muffling has to
        # be attempted rather than assumed.
        tryInvokeRestart("muffleProgression")
      }),
    handlers = progressr::handler_void())

  seen
}

test_that("the tick count is bounded and follows the run length", {
  expect_identical(progress_ticks(progress_control()), PROGRESS_TICKS)

  # A run shorter than the tick budget reports once per sweep rather than
  # promising more steps than it can take.
  expect_identical(progress_ticks(bartisan_control(num_burn = 3, num_draws = 4)),
                   7L)

  # Thinned sweeps count: the work is the sweeps, not the retained draws.
  expect_identical(
    progress_ticks(bartisan_control(num_burn = 0, num_draws = 5, num_thin = 4)),
    20L)

  # `bartisan_control()` will not build a run with nothing in it, so the
  # zero case is reachable only through the engine entry points, which take a
  # control list assembled by hand.
  expect_identical(
    progress_ticks(list(num_burn = 0L, num_draws = 0L, num_thin = 1L)), 0L)
})

test_that("no reporter is built when there is nothing to report", {
  reporter <- progress_reporter(1L, list(num_burn = 0L, num_draws = 0L,
                                         num_thin = 1L))

  expect_null(reporter[["report"]])
  expect_identical(reporter[["ticks"]], 0L)
})

test_that("one chain signals exactly the promised number of steps", {
  skip_if_not_installed("progressr")

  d <- sim_progress()

  # The sampler is asked for `PROGRESS_TICKS` reports and has to deliver that
  # many: too few and a bar never fills, too many and *progressr* warns about
  # overshooting the steps it was sized for. The two extra conditions are the
  # progressor's own start and finish.
  seen <- count_progress(
    bartisan(y ~ ., data = d, family = gaussian(),
             control = progress_control()))

  expect_identical(seen, PROGRESS_TICKS + 2L)
})

test_that("the count is shared across chains rather than restarting", {
  skip_if_not_installed("progressr")

  d <- sim_progress(seed = 2L)

  # One bar for the whole fit, so three chains signal three times as much
  # against a total sized for three. A per-chain progressor would fill and reset
  # three times instead.
  seen <- count_progress(
    bartisan(y ~ ., data = d, family = gaussian(),
             control = progress_control(chains = 3)))

  expect_identical(seen, 3L * PROGRESS_TICKS + 2L)
})

test_that("progress reaches the caller from parallel workers", {
  skip_if_not_installed("progressr")
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_on_cran()

  d <- sim_progress(seed = 3L)

  # Forked workers, which is the backend that needs no serialization of the
  # reporter. A socket backend relays too, but it reloads the package in each
  # worker, so it tests the installed build rather than this one.
  old <- future::plan(future::multicore, workers = 2)
  on.exit(future::plan(old), add = TRUE)

  seen <- count_progress(
    bartisan(y ~ ., data = d, family = gaussian(),
             control = progress_control(chains = 2)))

  expect_identical(seen, 2L * PROGRESS_TICKS + 2L)
})

test_that("the draws do not depend on whether progress was reported", {
  skip_if_not_installed("progressr")

  d <- sim_progress(seed = 4L)

  set.seed(7)
  quiet <- bartisan(y ~ ., data = d, family = gaussian(),
                    control = progress_control())

  rlang::local_options(progressr.enable = TRUE)
  set.seed(7)
  loud <- progressr::with_progress(
    bartisan(y ~ ., data = d, family = gaussian(),
             control = progress_control()),
    handlers = progressr::handler_void())

  # Reporting must not consume randomness or otherwise touch the chain.
  expect_equal(quiet[["eta"]], loud[["eta"]])
})
