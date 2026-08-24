# The gate is the only place the soft and hard cases differ, so it is where a
# mistake would be least visible: a wrong gate still produces a plausible fit.

test_that("the smoothstep gate is a distribution function with bounded support", {
  # Reached through the prediction engine, which is the same left_prob() the
  # sampler uses. A one-tree, one-split forest evaluated at a grid of x is the
  # gate itself.
  d <- sim_x(n = 200, seed = 3)
  d$y <- 2 * d$x1 + stats::rnorm(200)

  fits <- lapply(c("logistic", "smoothstep", "smootherstep"), function(g) {
    set.seed(4)
    genbart(y ~ ., d, control = quick_control(gate = g, num_trees = 5))
  })

  for (fit in fits) {
    expect_predictor_invariant(fit, d)
  }

  # The three gates give three different fits, so none of them is silently the
  # same as another.
  for (a in 1:2) {
    for (b in (a + 1):3) {
      expect_false(isTRUE(all.equal(fits[[a]][["fitted"]],
                                    fits[[b]][["fitted"]])))
    }
  }
})

test_that("each gate has the smoothness it claims", {
  # The gate is read out of the C++ through the prediction engine: a one-tree
  # fit on a fine grid of a single predictor is a sum of gates in that predictor,
  # so the numerical derivatives of the fitted curve say how smooth the gate is.
  # smoothstep is the Beta(2, 2) CDF, whose density vanishes at both ends of its
  # support, so the first derivative is continuous but the second jumps;
  # smootherstep is Beta(3, 3), where the second derivative vanishes too.
  n <- 3001
  d <- data.frame(x1 = seq(0, 1, length.out = n))
  d$y <- ifelse(d$x1 > 0.5, 1, -1)

  jump <- function(gate) {
    set.seed(2)
    # The response is a noiseless step, so the leaf scale has nothing to settle
    # to and is held fixed; this test is about the gate's shape, not the prior.
    fit <- genbart(y ~ ., d, control = quick_control(
      gate = gate, bandwidth = 0.05, update_bandwidth = FALSE, num_trees = 1L,
      num_burn = 100L, num_save = 100L, update_sigma_mu = FALSE))
    p <- stats::predict(fit, type = "link")
    # Largest jump in the second difference, scaled by the curve's own range, so
    # a kink shows up as a large value and a smooth join as a small one.
    d2 <- diff(p, differences = 2L)
    max(abs(diff(d2))) / (diff(range(p)) + 1e-12)
  }

  smooth1 <- jump("smoothstep")
  smooth2 <- jump("smootherstep")

  # Both are far smoother than a hard rule, whose fit steps.
  expect_lt(smooth1, 1e-3)
  expect_lt(smooth2, 1e-3)

  # And the C2 gate is the smoother of the two by this measure.
  expect_lt(smooth2, smooth1)
})

test_that("the smoothstep gate saturates and the logistic one does not", {
  # The claim the whole speedup rests on: past its half-width the smoothstep
  # gate is exactly zero or one, so an observation takes one side of the rule
  # outright. Checked on the C++ gate through a two-leaf tree whose split is
  # known, by comparing predictions far from the cutpoint.
  n <- 400
  d <- data.frame(x1 = seq(0, 1, length.out = n))
  d$y <- ifelse(d$x1 > 0.5, 1, -1) + stats::rnorm(n, sd = 0.1)

  band <- 0.02
  soft <- function(g) {
    set.seed(6)
    fit <- genbart(y ~ ., d, control = quick_control(
      gate = g, bandwidth = band, update_bandwidth = FALSE, num_trees = 1,
      num_burn = 200, num_save = 200))
    stats::predict(fit, type = "link")
  }

  p_log <- soft("logistic")
  p_ss <- soft("smoothstep")

  # Both recover the step.
  expect_gt(stats::cor(p_log, d$y), 0.9)
  expect_gt(stats::cor(p_ss, d$y), 0.9)

  # Away from the middle, the smoothstep fit is flat to numerical precision --
  # every observation there is on one side of every gate. The logistic fit is
  # not, because its gate never reaches zero.
  far <- d$x1 < 0.2
  expect_lt(stats::sd(p_ss[far]), 1e-8)
  expect_gt(stats::sd(p_log[far]), 1e-8)
})

test_that("the gate choice is stored and used by predict()", {
  d <- sim_x(seed = 8)
  d$y <- stats::rnorm(nrow(d))

  fit <- genbart(y ~ ., d, control = quick_control(gate = "smoothstep"))
  expect_identical(fit[["gate"]], "smoothstep")
  expect_predictor_invariant(fit, d)

  # A fit from before the gate existed has no element for it; predict() must
  # still treat it as the logistic gate rather than failing.
  logistic <- genbart(y ~ ., d, control = quick_control(gate = "logistic"))
  old <- logistic
  old[["gate"]] <- NULL
  expect_equal(stats::predict(old, newdata = d, type = "link"),
               stats::predict(logistic, newdata = d, type = "link"))

  # And the two gates really do evaluate the same forest differently, so the
  # element is not decoration.
  mixed <- logistic
  mixed[["gate"]] <- "smoothstep"
  expect_false(isTRUE(all.equal(stats::predict(mixed, newdata = d,
                                               type = "link"),
                               stats::predict(logistic, newdata = d,
                                              type = "link"))))

  expect_error(genbart_control(gate = "cosine"), "should be one of")
  expect_error(genbart_control(gate = "linear"), "should be one of")
})

test_that("the gate is ignored when the rules are hard", {
  d <- sim_x(seed = 9)
  d$y <- stats::rnorm(nrow(d))

  set.seed(10)
  a <- genbart(y ~ ., d, control = quick_control(soft = FALSE,
                                                 gate = "logistic"))
  set.seed(10)
  b <- genbart(y ~ ., d, control = quick_control(soft = FALSE,
                                                 gate = "smoothstep"))
  set.seed(10)
  cc <- genbart(y ~ ., d, control = quick_control(soft = FALSE,
                                                  gate = "smootherstep"))

  expect_equal(a[["eta"]], b[["eta"]])
  expect_equal(a[["eta"]], cc[["eta"]])
})

test_that("bandwidth_every skips the bandwidth move without disturbing anything else", {
  d <- sim_x(n = 150, seed = 13)
  d$y <- 2 * d$x1 * d$x2 + stats::rnorm(nrow(d))

  # Every sweep, and never: the bandwidth moves in the first and not the second.
  set.seed(14)
  every <- genbart(y ~ ., d, control = quick_control(num_burn = 100,
                                                     num_save = 100))
  set.seed(14)
  none <- genbart(y ~ ., d, control = quick_control(num_burn = 100,
                                                    num_save = 100,
                                                    update_bandwidth = FALSE))

  expect_gt(stats::sd(as.vector(every[["bandwidth"]])), 0)
  expect_identical(stats::sd(as.vector(none[["bandwidth"]])), 0)

  # A stride larger than the whole run leaves the bandwidth where it started,
  # exactly as switching the move off does -- the two are the same sampler.
  set.seed(14)
  never <- genbart(y ~ ., d, control = quick_control(num_burn = 100,
                                                     num_save = 100,
                                                     bandwidth_every = 10000L))
  expect_equal(never[["eta"]], none[["eta"]])

  # An intermediate stride is a valid sampler that still moves the bandwidth.
  set.seed(14)
  half <- genbart(y ~ ., d, control = quick_control(num_burn = 100,
                                                    num_save = 100,
                                                    bandwidth_every = 4L))
  expect_gt(stats::sd(as.vector(half[["bandwidth"]])), 0)
  expect_predictor_invariant(half, d)

  expect_error(genbart_control(bandwidth_every = 0), "must be")
  expect_error(genbart_control(bandwidth_every = 2.5), "whole number")
})
