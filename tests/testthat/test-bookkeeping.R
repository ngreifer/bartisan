# The sampler's tree bookkeeping was reworked for speed: nodes are recycled
# rather than allocated per proposal, the child weights come out of the split
# that produces them, the supports are snapshotted so a rejected bandwidth move
# can be undone by copying, and a tree with no splits draws its bandwidth from
# its prior directly. None of that is supposed to change what the sampler
# computes, so this file checks the things that would break if it did.

test_that("recycling nodes leaves the fit unchanged", {
  # Every field a recycled node starts with comes from the same function a fresh
  # one uses, so a stale field is the failure mode to look for. It would show up
  # as the stored forest disagreeing with the predictor the sampler recorded --
  # a wrong cutpoint, variable or leaf value replays differently.
  d <- sim_x(n = 200, seed = 51)
  set.seed(1051)
  d$y <- 2 * sin(pi * d$x1) - d$x2 + stats::rnorm(nrow(d))

  for (gate in c("smoothstep", "hard")) {
    fit <- bartisan(y ~ ., d, control = quick_control(
      gate = gate, num_trees = 25L, num_burn = 200L, num_draws = 200L))
    expect_predictor_invariant(fit, d)
  }

  # A long run with many births and deaths, so the pool is exercised heavily,
  # and with a deep tree prior so nodes are taken and given at several depths.
  deep <- bartisan(y ~ ., d, control = quick_control(
    gate = "hard", gamma = 0.99, beta = 0.5, num_trees = 10L,
    num_burn = 300L, num_draws = 300L))
  expect_predictor_invariant(deep, d)
})

test_that("the bandwidth of a tree with no splits is drawn from its prior", {
  skip_on_cran()

  # With the branching probability at zero the tree-shape prior forbids splits,
  # so every tree stays a single leaf, the bandwidth never enters the likelihood,
  # and its conditional is exactly the Exp(bandwidth) prior. That makes this the
  # one part of the sampler with a closed-form answer to check against.
  set.seed(61)
  d <- data.frame(x1 = stats::runif(300), x2 = stats::runif(300))
  d$y <- stats::rnorm(300)

  fit <- bartisan(y ~ ., d, gamma = 1e-10, bandwidth = 0.1, num_trees = 20L,
                  num_burn = 200L, num_draws = 2000L)

  b <- as.vector(fit[["bandwidth"]])
  expect_gt(length(b), 1000)

  # Mean and standard deviation of Exp(mean = 0.1) are both 0.1.
  expect_equal(mean(b), 0.1, tolerance = 0.03)
  expect_equal(stats::sd(b), 0.1, tolerance = 0.05)

  # And the whole distribution, not just its first two moments.
  expect_gt(suppressWarnings(stats::ks.test(b, "pexp", 10))$p.value, 0.01)

  # A different prior scale moves it, so the draw is using the prior rather than
  # anything hard-coded.
  wide <- bartisan(y ~ ., d, gamma = 1e-10, bandwidth = 0.4, num_trees = 20L,
                   num_burn = 200L, num_draws = 2000L)
  expect_equal(mean(as.vector(wide[["bandwidth"]])), 0.4, tolerance = 0.03)
})

test_that("a rejected bandwidth move leaves the tree exactly as it was", {
  # The snapshot has to be a faithful stand-in for recomputing every gate. If it
  # were not, the memberships and the predictor would drift apart, which is what
  # the predictor invariant detects.
  d <- sim_x(n = 250, seed = 71)
  set.seed(1071)
  d$y <- 2 * d$x1 * d$x2 + stats::rnorm(nrow(d))

  for (gate in c("logistic", "smoothstep", "smootherstep")) {
    fit <- bartisan(y ~ ., d, control = quick_control(
      gate = gate, num_trees = 20L, num_burn = 300L, num_draws = 300L))
    expect_predictor_invariant(fit, d)
    # The bandwidth actually moved, so rejections really happened.
    expect_gt(stats::sd(as.vector(fit[["bandwidth"]])), 0)
  }
})

test_that("the child weights the split records are the ones the target needs", {
  # Dividing a node's support and recording what each side got used to be two
  # passes over the same gates. Fusing them is only safe if the recorded weights
  # still sum to the parent's, which is what the whole membership scheme rests
  # on -- if they did not, the predictor and the stored forest would disagree.
  d <- sim_x(n = 200, seed = 81)
  set.seed(1081)
  d$y <- stats::rbinom(nrow(d), 1, stats::plogis(2 * d$x1 - d$x2))

  for (gate in c("smoothstep", "hard")) {
    fit <- bartisan(y ~ ., d, family = binomial("probit"),
                    control = quick_control(gate = gate, num_trees = 20L,
                                            num_burn = 200L, num_draws = 200L))
    expect_predictor_invariant(fit, d)
  }

  # Including with missing values, where a rule decides the side outright and
  # the weight is hard even when the rules are soft. `na.pass` is what keeps the
  # rows, rather than the default `na.omit` dropping them.
  d2 <- d
  d2$x3[c(4, 19, 60, 130)] <- NA
  fit <- bartisan(y ~ ., d2, family = binomial("probit"),
                  na.action = stats::na.pass,
                 control = quick_control(num_trees = 20L, num_burn = 200L,
                                         num_draws = 200L))
  expect_predictor_invariant(fit, d2)
})

# `get_varnames()` is `all.vars()` except that an extraction with `$`, `[[` or
# `[` is one variable rather than the pieces it is written from, which is how
# `model.frame()` names it.
test_that("get_varnames keeps an extraction whole where all.vars splits it", {
  expect_identical(get_varnames(quote(x)), "x")
  expect_identical(get_varnames(quote(log(x))), "x")
  expect_identical(get_varnames(quote(poly(x, 2))), "x")
  expect_identical(get_varnames(quote(cbind(a, b))), c("a", "b"))

  # The three extraction operators, each of which `all.vars()` mishandles: `$`
  # splits into the object and the column, and `[[` and `[` lose the column.
  expect_identical(get_varnames(quote(d$g)), "d$g")
  expect_identical(get_varnames(quote(d[["g"]])), 'd[["g"]]')
  expect_identical(get_varnames(quote(log(d$g + 1))), "d$g")
  expect_identical(get_varnames(quote(x + d$g)), c("x", "d$g"))

  expect_identical(all.vars(quote(d$g)), c("d", "g"))
  expect_identical(all.vars(quote(d[["g"]])), "d")
})

test_that("get_varnames walks a terms object, which as.list() would derail", {
  d <- data.frame(y = 1:3, x1 = 1:3, x2 = 1:3)
  tt <- stats::terms(y ~ x1 + d$x2, data = d)

  expect_identical(get_varnames(tt), c("y", "x1", "d$x2"))
  expect_identical(get_varnames(stats::delete.response(tt)), c("x1", "d$x2"))

  # A terms object is a call carrying attributes, and `as.list()` on one does
  # not drop the operator the way it does on a plain call, so a traversal built
  # on `as.list(e)[-1L]` walks the wrong tree. Indexing positions avoids it.
  expect_length(as.list(tt)[-1L], 3L)

  # The same formula without the attributes must give the same answer.
  expect_identical(get_varnames(y ~ x1 + d$x2), c("y", "x1", "d$x2"))
})

# The smoothed distribution function is swept rather than computed pairwise: an
# observation the kernel has entirely passed contributes a weight of exactly 1
# and one it has not reached contributes 0, so only those inside the window are
# evaluated. These check the sweep against the definition it is an optimization
# of, which is the only thing that makes the optimization safe.
test_that("the smooth CDF sweep matches evaluating every pair", {
  pairwise <- function(x, grid, h, kernel) {
    W <- if (kernel == 0L) {
      function(u) ifelse(u <= -1, 0, ifelse(u >= 1, 1, 0.75 * (u - u^3 / 3) + 0.5))
    }
    else {
      function(u) stats::pnorm(u)
    }
    t(vapply(grid, function(v) {
      w <- W((v - x) / h)
      c(mean(w), sum(w^2), sum((w - (x <= v))^2))
    }, numeric(3L)))
  }

  set.seed(4)
  # Skewed, with an outlier at each end, so the window moves unevenly.
  x <- sort(c(stats::rlnorm(198), -3, 40))
  grid <- sort(c(seq(min(x) - 1, max(x) + 1, length.out = 61L), x[c(1L, 200L)]))
  h <- min(stats::sd(x), stats::IQR(x) / 1.349) * 200^(-1 / 3)

  # Exact for a kernel with bounded support: nothing outside the window is
  # discarded, it is genuinely 0 or 1 there.
  epan <- .bartisan_smooth_cdf(x, grid, h, EPANECHNIKOV)
  expect_equal(epan, pairwise(x, grid, h, 0L), tolerance = 1e-12,
               ignore_attr = TRUE)

  # The Gaussian reaches everywhere and is cut at five bandwidths, so what it
  # drops is a tail rather than nothing.
  gauss <- .bartisan_smooth_cdf(x, grid, h, GAUSSIAN)
  expect_equal(gauss, pairwise(x, grid, h, 1L), tolerance = 1e-5,
               ignore_attr = TRUE)

  # The distribution function it reports is one, whatever the kernel.
  for (m in list(epan, gauss)) {
    expect_false(is.unsorted(m[, 1L]))
    expect_true(all(m[, 1L] >= 0 & m[, 1L] <= 1))
  }

  expect_error(.bartisan_smooth_cdf(x, grid, 0, EPANECHNIKOV), "bandwidth")
  expect_error(.bartisan_smooth_cdf(numeric(), grid, h, EPANECHNIKOV), "at least one")
})

test_that("the smoothcdf map is increasing and spans the unit interval", {
  set.seed(11)

  # The outliers are placed rather than drawn, so how extreme they are is not a
  # property of the seed, and there are three of them in 300 so that the central
  # 98% checked below is the bulk rather than the bulk plus most of the tail.
  shapes <- list(lognormal = stats::rlnorm(300),
                 outliers = c(stats::rnorm(297), c(200, 500, 1000)),
                 bimodal = c(stats::rnorm(150, 0, .3), stats::rnorm(150, 8, .3)),
                 ties = rep(seq_len(12), length.out = 300))

  for (nm in names(shapes)) {
    x <- shapes[[nm]]
    m <- make_unit_map(x, "smoothcdf")
    u <- m(sort(x))

    expect_false(is.unsorted(u))
    expect_equal(range(u), c(0, 1), tolerance = 1e-8)
    expect_true(is.na(m(NA_real_)))

    # Clamped outside the observed range, as the other maps are.
    expect_identical(m(min(x) - 1e6), 0)
    expect_identical(m(max(x) + 1e6), 1)
  }

  # The point of it: an outlier no longer compresses the bulk into a sliver the
  # cutpoint prior cannot reach, which is what `"range"` does here.
  x <- shapes[["outliers"]]
  bulk <- function(type) {
    diff(range(make_unit_map(x, type)(stats::quantile(x, c(.01, .99)))))
  }
  expect_gt(bulk("smoothcdf"), 0.9)
  expect_lt(bulk("range"), 0.05)
})
