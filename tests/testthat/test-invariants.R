# Invariants that hold across families and across the structural wrappers, which
# is where the composition bugs live. Each of these is cheap and none of them
# needs a known truth: they are statements the package already makes about
# itself, turned into assertions.
#
# The one that matters most is the augmentation agreement. `?bartisan_control`
# says of `augment` that "the posterior is the same either way", and for a year
# it was not: a `vc()` model never forwarded `before_forest()` to the family it
# wrapped, so an augmented family's augmentation was never refreshed and a logit
# binomial recovered a sixth of the effect it was given. Nothing in the code
# looked wrong, because the bug was an absence.

sim_invariant <- function(seed = 3L, n = 600L) {
  set.seed(seed)
  d <- as.data.frame(matrix(stats::runif(n * 3L), n, 3L))
  names(d) <- paste0("x", 1:3)
  f <- 1.2 * d$x1 - 1.2 * d$x2
  d$z <- stats::rbinom(n, 1L, 0.5)
  lp <- f - mean(f) + 0.8 * d$z

  d$ybin <- stats::rbinom(n, 1L, stats::plogis(lp))
  d$yord <- cut(lp + stats::rnorm(n), c(-Inf, -0.5, 0.5, Inf), labels = 1:3,
                ordered_result = TRUE)
  d$ynb <- stats::rnbinom(n, mu = exp(lp), size = 3)
  d
}

invariant_cases <- list(
  list(name = "binomial logit", family = binomial(), response = "ybin"),
  list(name = "binomial probit", family = binomial("probit"), response = "ybin"),
  list(name = "ordinal", family = ordinal(), response = "yord"),
  list(name = "negbin", family = negbin(), response = "ynb")
)

invariant_fit <- function(case, shape, augment, d) {
  form <- {
    if (identical(shape, "plain")) {
      stats::reformulate(c("z", "x1", "x2", "x3"), response = case$response)
    }
    else {
      stats::reformulate(c("x1", "x2", "x3", "vc(z)"), response = case$response)
    }
  }

  bartisan(form, data = d, family = case$family,
           control = bartisan_control(
             num_trees = if (identical(shape, "vc")) c(20L, 10L) else 20L,
             num_burn = 200L, num_draws = 300L, chains = 1L, gate = "hard",
             augment = augment))
}

test_that("a discrete family reports a log likelihood that is not positive", {
  skip_on_cran()

  d <- sim_invariant()

  for (case in invariant_cases) {
    for (shape in c("plain", "vc")) {
      fit <- invariant_fit(case, shape, TRUE, d)

      # Every one of these has a probability mass function, so its log
      # likelihood is a sum of logs of numbers no greater than one. An augmented
      # family's *target* is the augmented density and may be anything; what is
      # reported is meant to be the likelihood the augmentation stands in for,
      # and reporting the target instead gave +37 where -830 was right.
      expect_lt(mean(fit[["loglik"]]),
                0, label = paste(case$name, shape, "reported log likelihood"))
    }
  }
})

test_that("augmenting a family does not move its posterior", {
  skip_on_cran()

  d <- sim_invariant()

  for (case in invariant_cases) {
    for (shape in c("plain", "vc")) {
      set.seed(9)
      with_aug <- invariant_fit(case, shape, TRUE, d)
      set.seed(9)
      without <- invariant_fit(case, shape, FALSE, d)

      # Two samplers for one posterior, so they agree up to Monte Carlo error
      # and not exactly. The tolerance is loose enough that a short chain passes
      # and tight enough that the failure this was written for, a coefficient
      # off by a factor of five, does not.
      one <- mean(colMeans(with_aug[["eta"]][[1L]]))
      two <- mean(colMeans(without[["eta"]][[1L]]))

      expect_equal(one, two, tolerance = 0.15 * (abs(two) + 0.5),
                   label = paste(case$name, shape, "predictor"))

      expect_equal(mean(with_aug[["loglik"]]), mean(without[["loglik"]]),
                   tolerance = 0.02 * abs(mean(without[["loglik"]])) + 5,
                   label = paste(case$name, shape, "log likelihood"))
    }
  }
})

test_that("every hook a decorator inherits is one it means to inherit", {
  skip_on_cran()
  skip_if_not(dir.exists(test_path("..", "..", "src")))

  # The varying-coefficient bug was a decorator that forwarded some of the base
  # class's hooks and silently inherited the rest. This is the mechanical form
  # of that question: for each wrapper, which of the base's virtual hooks does
  # it neither override nor deliberately leave alone? The answer should be a
  # list someone has looked at, which is what the expectation below pins.
  header <- readLines(test_path("..", "..", "src", "family.h"))
  hooks <- unique(sub(".*\\b(\\w+)\\($", "\\1",
                      regmatches(header, regexpr("virtual [^;{]*\\w+\\(",
                                                 header))))
  hooks <- hooks[nzchar(hooks)]

  expect_gt(length(hooks), 10L)

  body <- readLines(test_path("..", "..", "src", "family.cpp"))

  for (wrapper in c("LinkedFamily", "VaryingCoefficientFamily")) {
    at <- grep(paste0("^struct ", wrapper), body)
    expect_length(at, 1L)

    ends <- grep("^};", body)
    text <- paste(body[at:min(ends[ends > at])], collapse = "\n")

    # These two are what the bug was. Whatever else a wrapper decides about the
    # rest of the hooks, a decorator that does not pass a sweep boundary or a
    # likelihood report through to what it wraps is broken.
    expect_true(grepl("before_forest", text),
                label = paste(wrapper, "forwards before_forest()"))
    expect_true(grepl("reported_loglik", text),
                label = paste(wrapper, "forwards reported_loglik()"))
  }
})

test_that("a binomial fit keeps its invariants under every structural wrapper", {
  skip_on_cran()

  # The binomial row of the matrix in `_dev/AUDITING.md`. The wrappers compose,
  # and a hook dropped by one of them shows up only in the cell where it meets
  # a family that needed it, so the cells are what get tested rather than the
  # wrappers one at a time.
  set.seed(5)
  n <- 700L
  d <- as.data.frame(matrix(stats::runif(n * 3L), n, 3L))
  names(d) <- paste0("x", 1:3)
  d$g <- factor(sample(letters[1:8], n, TRUE))
  d$cat <- factor(sample(c("a", "b", "c"), n, TRUE))
  u <- stats::setNames(stats::rnorm(8L, 0, 0.6), letters[1:8])
  d$z <- stats::rbinom(n, 1L, 0.5)
  d$off <- stats::rnorm(n, 0, 0.3)
  f <- 1.2 * d$x1 - 1.2 * d$x2
  lp <- f - mean(f) + 0.8 * d$z + u[as.character(d$g)] + d$off
  d$y <- stats::rbinom(n, 1L, stats::plogis(lp))
  d$tr <- rep(4L, n)
  d$ycnt <- stats::rbinom(n, 4L, stats::plogis(lp)) / 4
  d$na1 <- d$x1
  d$na1[sample(n, 30L)] <- NA

  cells <- list(
    list("ranef", y ~ z + x1 + x2 + (1 | g), 15L, NULL),
    list("vc + ranef", y ~ x1 + x2 + vc(z) + (1 | g), c(15L, 8L), NULL),
    list("offset", y ~ z + x1 + x2, 15L, "off"),
    list("vc + offset", y ~ x1 + x2 + vc(z), c(15L, 8L), "off"),
    list("vc + weights", ycnt ~ x1 + x2 + vc(z), c(15L, 8L), "tr"),
    list("vc + categorical", y ~ x1 + cat + vc(z), c(15L, 8L), NULL),
    list("vc + missing", y ~ na1 + x2 + vc(z), c(15L, 8L), NULL)
  )

  for (cell in cells) {
    fits <- lapply(c(TRUE, FALSE), function(aug) {
      args <- list(formula = cell[[2L]], data = d, family = binomial(),
                   control = bartisan_control(num_trees = cell[[3L]],
                                              num_burn = 150L, num_draws = 200L,
                                              chains = 1L, gate = "hard",
                                              augment = aug))

      if (identical(cell[[4L]], "off")) args[["offset"]] <- d[["off"]]
      if (identical(cell[[4L]], "tr")) args[["weights"]] <- d[["tr"]]

      set.seed(13)
      do.call(bartisan, args)
    })

    for (fit in fits) {
      expect_lt(mean(fit[["loglik"]]), 0,
                label = paste(cell[[1L]], "reported log likelihood"))
    }

    expect_equal(mean(fits[[1L]][["loglik"]]), mean(fits[[2L]][["loglik"]]),
                 tolerance = 0.05 * abs(mean(fits[[2L]][["loglik"]])) + 10,
                 label = paste(cell[[1L]], "log likelihood under augmentation"))

    one <- mean(colMeans(fits[[1L]][["eta"]][[1L]]))
    two <- mean(colMeans(fits[[2L]][["eta"]][[1L]]))
    expect_equal(one, two, tolerance = 0.25 * (abs(two) + 0.5),
                 label = paste(cell[[1L]], "predictor under augmentation"))
  }
})

# What a closure would cost to send, measured with its source references
# removed.
#
# `serialize()` on a function also serializes its `srcref`, and a `srcref`
# carries the entire text of the file the function was defined in. An installed
# package has none, so the raw size is the retained data; but
# `pkgload::load_all()` keeps them, which is what `devtools::test()` and
# `testthat::test_local()` use, and there the text of `R/utils.R` adds 286 KB to
# every closure defined in it. The two thresholds below are tight enough that
# this drowned out what they measure, and the failure appeared only under
# `devtools::test()` and never under `R CMD check`. Stripping the references
# leaves the environment chain, which is the thing at issue: a closure holding a
# 20,000-element vector still measures 320 KB through this, and one holding an
# unforced promise to the frame containing it still measures 160 KB.
closure_bytes <- function(f) {
  length(serialize(utils::removeSource(f), NULL))
}

# A closure is sent to every worker, and a closure carries the frame it was
# written in. Nothing about reading one says how large it is, which is what makes
# this worth asserting rather than reviewing: the reporter `diagnose()` hands to
# its workers was four lines long and 336 MB to send, because it kept a promise
# alive that reached back to the frame holding the fit. Both ways of leaking a
# frame are covered here -- writing the closure where the frame is, and leaving
# an argument unforced -- because either one alone brings the fit along.
test_that("a progress reporter carries the progressor and nothing else", {
  skip_if_not_installed("progressr")

  # Shaped like `diagnose()`: the frame that asks for a reporter is also the
  # frame holding the draws, so `envir` points straight at them.
  asks_for_one <- function() {
    ballast <- numeric(2e6)
    list(report = diagnosis_reporter(10L),
         stepper = progress_budget(diagnosis_reporter(10L), 100L, 10L)(50L),
         ballast_size = length(serialize(ballast, NULL)))
  }

  got <- progressr::with_progress(asks_for_one(),
                                  handlers = progressr::handler_void())

  expect_gt(got[["ballast_size"]], 1e7)
  expect_lt(closure_bytes(got[["report"]]), 1e5)
  expect_lt(closure_bytes(got[["stepper"]]), 1e5)
})

# The same invariant one layer down, and the reason the two are next to each
# other: a unit map is a closure stored in the fit, one per predictor column, so
# whatever it retains is paid for as long as the fit exists. Both hazards were
# live here at once. The closures were written in `make_unit_map()`'s frame, so
# they held `x` and `ux` whether they used them or not; and the two-value branch
# never looked at `type`, so that argument stayed an unforced promise and kept
# the frame of the *caller* alive, which is `unit_transform()`'s and has the
# design matrix in it twice. A binary predictor's map, whose whole content is two
# numbers, serialized to 148 KB on a 400-row fit.
test_that("a unit map carries only the numbers it uses", {
  # Shaped like `unit_transform()`: the frame asking for the maps is also the one
  # holding the design matrix, and it holds a second copy while mapping it.
  asks_for_maps <- function(x, type) {
    out <- x
    maps <- lapply(seq_len(ncol(x)), function(j) make_unit_map(x[, j], type))

    for (j in seq_len(ncol(x))) {
      out[, j] <- maps[[j]](x[, j])
    }

    maps
  }

  set.seed(4L)
  n <- 2000L
  x <- cbind(continuous = stats::rnorm(n),
             binary = stats::rbinom(n, 1L, 0.4),
             constant = rep(1, n))
  whole <- length(serialize(x, NULL))

  # Without this the thresholds below could be met by there being nothing to
  # leak in the first place.
  expect_gt(whole, 4e4)

  for (type in c("quantile", "range")) {
    sizes <- vapply(asks_for_maps(x, type), closure_bytes, numeric(1L))

    # Two numbers, or none, however long the column is.
    expect_lt(sizes[[2L]], 2e3)
    expect_lt(sizes[[3L]], 2e3)
  }

  # A quantile map holds an `ecdf` and so has to hold the values it was built
  # from, but its own column only; a range map is two numbers like the rest.
  expect_lt(closure_bytes(asks_for_maps(x, "quantile")[[1L]]), whole)
  expect_lt(closure_bytes(asks_for_maps(x, "range")[[1L]]), 2e3)
})
