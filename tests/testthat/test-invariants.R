# Invariants that hold across families and across the structural wrappers, which
# is where the composition bugs live. Each of these is cheap and none of them
# needs a known truth: they are statements the package already makes about
# itself, turned into assertions.
#
# The one that matters most is the augmentation agreement. `?bartisan_control`
# says of `augment` that "the posterior is the same either way", and for four days
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

  # The varying-coefficient bugs were a decorator that forwarded some of the
  # base class's hooks and silently inherited the rest: `before_forest()` and
  # `reported_loglik()` first, then `report_shift()`, which left a `bcf()` fit
  # with `family = dpm()` reporting the raw forest and its level unidentified.
  # Reading the wrapper finds none of these, because the bug is an absence.
  # This is the mechanical form of the question: for each wrapper, which of the
  # base's virtual hooks does it neither override nor deliberately leave alone?
  # The lists below are the ones someone has looked at, with the reason each
  # hook can be inherited. A hook added to `family.h` lands in the difference
  # until it is either forwarded or added to a list with its reason.
  # The header is searched as one string, since a declaration's return type
  # and its name can sit on different lines.
  header <- paste(readLines(test_path("..", "..", "src", "family.h")),
                  collapse = "\n")
  hooks <- regmatches(header,
                      gregexpr("(?m)^\\s*virtual\\s+[^;{(]*?\\b\\w+\\s*\\(",
                               header, perl = TRUE))[[1L]]
  hooks <- unique(sub(".*\\b(\\w+)\\s*\\($", "\\1", hooks))
  hooks <- setdiff(hooks[nzchar(hooks)], "Family")

  expect_gt(length(hooks), 10L)

  body <- readLines(test_path("..", "..", "src", "family.cpp"))
  ends <- grep("^};", body)

  # Declarations span lines, so the struct is searched as one string: a hook
  # is overridden when its name opens a parameter list that closes on
  # `override`, which a forwarding call to `inner->hook()` never does.
  inherited_by <- function(wrapper) {
    at <- grep(paste0("^struct ", wrapper, "\\b"), body)
    expect_length(at, 1L)
    text <- paste(body[at:min(ends[ends > at])], collapse = "\n")

    overridden <- vapply(hooks, function(hook) {
      grepl(sprintf("\\b%s\\s*\\([^;]*?\\)\\s*(const\\s*)?override", hook),
            text, perl = TRUE)
    }, logical(1L))

    sort(hooks[!overridden])
  }

  # Generic over the unit hooks every wrapper does override: the blocked sums
  # and the likelihood deltas call `logdens_block()`, `logdens_unit()` and the
  # score hooks, so a wrapper that forwards those gets these right for free.
  generic <- c("accumulate1", "accumulate1_at", "accumulate1_from",
               "accumulate2", "accumulate2_from", "loglik_delta",
               "loglik_delta_rows")

  expect_setequal(inherited_by("VaryingCoefficientFamily"), generic)

  # A link applied from R puts the predictor on the caller's scale, so a shift
  # of the wrapped family's predictor is not a shift of this one, and the
  # sampler's own chart is the one to report in: `report_shift()`,
  # `aux_values_shifted()` and `mixture_flat()` stay at the base. The
  # derivatives are central differences through the link by design, the target
  # form is general for the same reason, which makes `exp_rate()` and the
  # coding exactness moot, and the wrapper has one predictor, so no family with
  # pinned forests reaches it. None of the families `with_link()` accepts
  # shifts its chart or pins a forest; one that did would have to be decided
  # here.
  expect_setequal(inherited_by("LinkedFamily"),
                  c(generic, "aux_values_shifted", "coding_is_exact",
                    "coding_not_exact", "dlogdens_unit", "exp_rate",
                    "info_unit", "mixture_flat", "num_pinned", "report_shift",
                    "target_form"))

  # `RFamily` wraps an R function rather than another family, so what it does
  # not override is its implementation and not a forward it forgot: the base
  # defaults are the behavior of a likelihood that is a black box, and as the
  # innermost family its own hook count is the answer to `sweep_delivered()`.
  expect_setequal(inherited_by("RFamily"),
                  c(generic, "aux_values_shifted", "before_forest",
                    "coding_is_exact", "coding_not_exact", "compute_eta_free",
                    "dlogdens_unit", "exp_rate", "info_unit", "log_norm_const",
                    "logdens_extra_total", "mixture_flat", "report_shift",
                    "reported_loglik", "sweep_delivered", "target_form"))
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

# An S3 method registered on another package's generic is only reachable if
# A method on another package's generic is reached only if it is registered, and
# the registration comes from a roxygen tag whose attachment is positional:
# inserting a helper between `@exportS3Method` and the function moves the
# registration onto the helper, and the method is then never called.
#
# Asserted against the method registry rather than against `NAMESPACE`, which is
# what an earlier version of this test read. The registry is what actually
# governs dispatch, it exists in every context the tests run in, and it records
# the package each generic belongs to, so it can tell a method on
# `bayesplot::pp_check` from one on a local generic of the same name. Reading
# the file could do none of those: `../../NAMESPACE` resolves in
# `tests/testthat/` and not under `R CMD check`, where the tests run from
# `bartisan.Rcheck/tests/`, so the test written to catch what only breaks in an
# installed package was the one thing that broke against an installed package.
test_that("every S3 method on a foreign generic is registered", {
  methods <- getNamespaceInfo("bartisan", "S3methods")

  foreign <- c("bayesplot::pp_check", "loo::loo", "loo::waic", "loo::kfold",
               "posterior::as_draws", "performance::model_performance",
               "insight::get_data", "marginaleffects::get_predict",
               "rstantools::posterior_predict", "rstantools::log_lik",
               "rstantools::prior_summary")

  for (generic in foreign) {
    parts <- strsplit(generic, "::", fixed = TRUE)[[1L]]

    found <- methods[, 1L] == parts[[2L]] &
      methods[, 2L] == "bartisan_fit" &
      !is.na(methods[, 4L]) & methods[, 4L] == parts[[1L]]

    expect_true(any(found),
                info = paste(generic, "has no registered bartisan_fit method"))
  }
})

# A print method writes to stdout. `cli::cli_bullets()` writes to stderr, so
# every bullet it produced was missing from `capture.output()` and from a
# knitted document: `vignette("bartisan")` showed `diagnose()`'s "What to do"
# heading with nothing under it, and `estimate_effect()`'s legend naming the
# contrast's levels never appeared at all. The bullets now go through
# `cli_bullets_cat()`, and this asserts none of them found their way back.
test_that("print methods write nothing to stderr", {
  skip_on_cran()

  d <- sim_x(n = 200L, p = 3L, seed = 74L)
  d$z <- stats::rbinom(nrow(d), 1L, 0.5)
  d$y <- stats::rbinom(nrow(d), 1L, stats::plogis(d$x1 + 0.8 * d$z))

  fit <- bartisan(y ~ ., d, family = stats::binomial(), chains = 2L,
                  control = quick_control(num_trees = 5L, num_burn = 60L,
                                          num_draws = 100L))

  shown <- list(fit = fit,
                summary = summary(fit),
                diagnosis = diagnose(fit),
                importance = variable_importance(fit),
                effect = estimate_effect(fit, treat = "z"),
                lnor = estimate_effect(fit, treat = "z",
                                       comparison = "lnor"),
                cate = estimate_effect(fit, treat = "z",
                                       estimand = "CATE"),
                partial = partial_dependence(fit, ~ x1))

  for (what in names(shown)) {
    on_err <- utils::capture.output(type = "message",
                                    print(shown[[what]]))

    expect_identical(on_err, character(), info = what)
  }

  # And the content that used to be lost is on stdout, where a reader of a
  # vignette will see it.
  expect_match(printed_text(shown[["diagnosis"]]), "What to do")
  expect_match(printed_text(shown[["lnor"]]), "O(y) is the odds", fixed = TRUE)
})

# The family-by-structure matrix of `_dev/AUDITING.md` § 3. One data set with a
# predictor, a binary treatment and a grouping factor, and a response for every
# family drawn from a known linear predictor on that family's own scale. The
# survival responses each use their family's own error so that the predictor
# is the location of log T in every case.
sim_matrix <- function(n = 400L, seed = 31L) {
  set.seed(seed)
  d <- as.data.frame(matrix(stats::runif(n * 3L), n, 3L))
  names(d) <- paste0("x", 1:3)
  d$z <- stats::rbinom(n, 1L, 0.5)
  d$g <- factor(sample(letters[1:6], n, TRUE))
  d$off <- stats::rnorm(n, 0, 0.3)
  d$w <- sample(1:3, n, TRUE)
  u <- stats::setNames(stats::rnorm(6L, 0, 0.4), letters[1:6])
  f <- 1.2 * d$x1 - 1.2 * d$x2
  lp <- f - mean(f) + 0.8 * d$z + u[as.character(d$g)]

  d$ynorm <- lp + stats::rnorm(n, 0, 0.5)
  d$ybin <- stats::rbinom(n, 1L, stats::plogis(lp))
  d$yprob <- stats::rbinom(n, 1L, stats::pnorm(lp))
  d$ycnt <- stats::rpois(n, exp(lp - 0.5))
  d$ynb <- stats::rnbinom(n, mu = exp(lp - 0.5), size = 3)
  d$yzi <- d$ycnt * stats::rbinom(n, 1L, 0.75)
  d$ygam <- stats::rgamma(n, shape = 4, rate = 4 / exp(lp - 0.5))
  mu <- stats::plogis(lp)
  d$ybeta <- pmin(pmax(stats::rbeta(n, mu * 10, (1 - mu) * 10), 1e-4), 1 - 1e-4)
  d$yob <- d$ybeta
  d$yob[lp < stats::quantile(lp, 0.08)] <- 0
  d$yob[lp > stats::quantile(lp, 0.92)] <- 1
  events <- stats::rpois(n, exp(lp - 0.5))
  d$ytw <- ifelse(events > 0, stats::rgamma(n, shape = 2 * pmax(events, 1), rate = 2), 0)
  d$yord <- cut(lp + stats::rlogis(n), c(-Inf, -0.5, 0.5, Inf), labels = 1:3,
                ordered_result = TRUE)
  d$ymult <- factor(cut(lp + stats::rlogis(n), c(-Inf, -0.5, 0.5, Inf),
                        labels = c("a", "b", "c")), ordered = FALSE)
  d$ydpm <- lp + 2 * (stats::rgamma(n, 1.5, 1.5) - 1)

  cens <- lp + stats::rnorm(n, 0.6, 0.7)
  for (nm in c("weib", "llog", "lnorm", "dpma")) {
    e <- switch(nm,
                weib = log(stats::rexp(n)),
                llog = stats::rlogis(n),
                lnorm = stats::rnorm(n),
                dpma = 1.2 * (stats::rgamma(n, 1.5, 1.5) - 1))
    lt <- lp + 0.5 * e
    d[[paste0("t_", nm)]] <- exp(pmin(lt, cens))
    d[[paste0("e_", nm)]] <- as.integer(lt <= cens)
  }

  # Proportional hazards with a Weibull baseline.
  lt <- log((-log(stats::runif(n)) / exp(lp))^(1 / 1.5))
  d$t_ph <- exp(pmin(lt, cens))
  d$e_ph <- as.integer(lt <= cens)

  attr(d, "lp") <- lp
  d
}

# `structures` names the cells a family is fit in; `multinomial()` and `mnp()`
# refuse a varying coefficient by design. `mnp()`'s density is simulated rather
# than computed, so the density identity is not exact for it and is skipped;
# its replay, range and sign checks still run.
matrix_cases <- list(
  list(name = "gaussian", family = gaussian(), response = "ynorm"),
  list(name = "gaussian_ls", family = gaussian_ls(), response = "ynorm"),
  list(name = "binomial logit", family = binomial(), response = "ybin",
       discrete = TRUE),
  list(name = "binomial probit", family = binomial("probit"), response = "yprob",
       discrete = TRUE),
  list(name = "poisson", family = poisson(), response = "ycnt", discrete = TRUE),
  list(name = "negbin", family = negbin(), response = "ynb", discrete = TRUE),
  list(name = "zi_poisson", family = zi_poisson(), response = "yzi",
       discrete = TRUE),
  list(name = "zi_negbin", family = zi_negbin(), response = "yzi",
       discrete = TRUE),
  list(name = "Gamma", family = Gamma("log"), response = "ygam"),
  list(name = "Gamma_ls", family = Gamma_ls(), response = "ygam"),
  list(name = "Beta", family = Beta(), response = "ybeta"),
  list(name = "ordbeta", family = ordbeta(), response = "yob"),
  list(name = "tweedie", family = tweedie(), response = "ytw"),
  list(name = "ordinal logit", family = ordinal(), response = "yord",
       discrete = TRUE),
  list(name = "ordinal probit", family = ordinal("probit"), response = "yord",
       discrete = TRUE),
  list(name = "ordinal cloglog", family = ordinal("cloglog"), response = "yord",
       discrete = TRUE),
  list(name = "multinomial", family = multinomial(), response = "ymult",
       discrete = TRUE, structures = c("plain", "ranef", "offset", "weights")),
  list(name = "mnp", family = multinomial("probit", replicates = 20L),
       response = "ymult", discrete = TRUE, simulated_density = TRUE,
       structures = c("plain", "ranef", "offset", "weights")),
  list(name = "dpm", family = dpm(), response = "ydpm"),
  list(name = "weibull_aft", family = weibull_aft(),
       response = "cbind(t_weib, e_weib)"),
  list(name = "loglogistic_aft", family = loglogistic_aft(),
       response = "cbind(t_llog, e_llog)"),
  list(name = "lognormal_aft", family = lognormal_aft(),
       response = "cbind(t_lnorm, e_lnorm)"),
  list(name = "dpm_aft", family = dpm_aft(), response = "cbind(t_dpma, e_dpma)"),
  list(name = "ph", family = ph(), response = "cbind(t_ph, e_ph)"),
  list(name = "custom", response = "ynorm",
       family = custom_family(function(y, eta) {
         stats::dnorm(y, eta[, 1], 0.5, log = TRUE)
       }))
)

matrix_structures <- c("plain", "vc", "drawn", "ranef", "vc + ranef", "offset",
                       "weights")

matrix_formula <- function(response, structure) {
  terms <- switch(structure,
                  plain = ,
                  offset = ,
                  weights = c("z", "x1", "x2", "x3"),
                  vc = c("x1", "x2", "x3", "vc(z)"),
                  drawn = c("x1", "x2", "x3", 'vc(z, center = "estimate")'),
                  ranef = c("z", "x1", "x2", "x3", "(1 | g)"),
                  `vc + ranef` = c("x1", "x2", "x3", "vc(z)", "(1 | g)"))

  stats::as.formula(paste(response, "~", paste(terms, collapse = " + ")),
                    env = globalenv())
}

# NULL when the cell does not exist: a drawn coding needs a leaf target that is
# quadratic in the predictor, and the Dirichlet process families take no prior
# weights; both refusals say so. Any other error is a failure of the cell. The
# cells refused are collected so that the set can be pinned below.
matrix_skipped <- new.env()

matrix_fit <- function(case, structure, d, control) {
  form <- matrix_formula(case$response, structure)
  args <- list(form, data = d, family = case$family, control = control)

  if (identical(structure, "offset")) {
    args[["offset"]] <- d$off
  }

  if (identical(structure, "weights")) {
    args[["weights"]] <- d$w
  }

  tryCatch(
    do.call(bartisan, args),
    error = function(e) {
      msg <- conditionMessage(e)
      refused <- (identical(structure, "drawn") &&
                    grepl("leaf target is\\s+quadratic", msg)) ||
        (identical(structure, "weights") &&
           grepl("does not take prior weights", msg))

      if (refused) {
        assign(paste(case$name, structure), TRUE, envir = matrix_skipped)
        return(NULL)
      }

      stop(sprintf("%s under %s: %s", case$name, structure, msg), call. = FALSE)
    })
}

# Nuisance parameters that are scales, counts or rates and so cannot be
# negative, by the column names the families report them under.
positive_aux <- function(aux) {
  if (is.null(aux)) {
    return(numeric(0))
  }

  named <- intersect(colnames(aux),
                     c("sigma", "phi", "theta", "shape", "error_sd", "alpha",
                       "clusters", "power"))
  hazards <- grep("^lambda[0-9]+$", colnames(aux), value = TRUE)
  aux[, c(named, hazards), drop = FALSE]
}

test_that("every family keeps its chart and its density under every structural wrapper", {
  skip_on_cran()

  # Two identities that need no truth and hold draw by draw, so a short chain
  # is enough. The density route reproduces the recorded log likelihood only
  # when the recorded predictor and the family's recorded nuisance parameters
  # are in one chart, which is what a `bcf()` fit with `family = dpm()` lost.
  # Replaying the stored forests reproduces the recorded predictor only when
  # the leaves were written in the same chart as the predictor. The cells are
  # what get tested rather than the wrappers one at a time, because a hook
  # dropped by a wrapper shows up only where it meets a family that needed it.
  # Both gates, since the replay of a soft tree goes through the membership
  # weights and the bandwidth, which a hard tree never touches.
  d <- sim_matrix()
  rm(list = ls(matrix_skipped), envir = matrix_skipped)

  for (gate in c("hard", "smoothstep")) {
    control <- quick_control(num_trees = 5L, num_burn = 30L, num_draws = 30L,
                             gate = gate)

    for (case in matrix_cases) {
      for (structure in case$structures %||% matrix_structures) {
        fit <- matrix_fit(case, structure, d, control)

        if (is.null(fit)) {
          next
        }

        label <- paste(case$name, "under", structure, "with", gate, "rules")

        expect_true(all(is.finite(fit[["loglik"]])), label = label)

        if (isTRUE(case$discrete)) {
          expect_lt(mean(fit[["loglik"]]), 0,
                    label = paste(label, "reported log likelihood"))
        }

        if (!isTRUE(case$simulated_density)) {
          dens <- stats::predict(fit, type = "density", draws = TRUE,
                                 log = TRUE)
          expect_equal(rowSums(dens), fit[["loglik"]], tolerance = 1e-8,
                       label = paste(label, "density against log likelihood"))
        }

        eta <- stats::predict(fit, type = "link", draws = TRUE)
        # Values only: the replay carries the data's row names and the stored
        # predictor does not. An offset is part of the predictor the sampler
        # recorded, so the replay is handed it.
        replayed <- stats::predict(
          fit, newdata = d, type = "link", draws = TRUE,
          offset = if (identical(structure, "offset")) d$off)
        expect_equal(replayed, eta, tolerance = 1e-6, ignore_attr = TRUE,
                     label = paste(label, "replayed predictor"))

        # Range checks: nothing a fit reports should be outside its support.
        # A probability sits in [0, 1] and a row of them sums to one; a scale,
        # count or rate is positive. R-hat and the effective sample sizes are
        # finite and positive, which is as far as a range check can go on the
        # rank-normalized estimators: split R-hat can read a little below one
        # and bulk ESS can exceed the number of draws for antithetic chains, so
        # neither "at least one" nor "at most the draws" is a valid bound.
        response <- stats::predict(fit, type = "response")
        expect_true(all(is.finite(response)),
                    label = paste(label, "finite response"))

        if (fit[["family"]][["family"]] %in% c("binomial", "ordinal",
                                                "multinomial", "mnp")) {
          prob <- stats::predict(fit, type = "prob")
          expect_true(all(prob >= 0 & prob <= 1),
                      label = paste(label, "probabilities in [0, 1]"))

          if (is.matrix(prob) && ncol(prob) > 1L) {
            expect_equal(unname(rowSums(prob)), rep.int(1, nrow(prob)),
                         tolerance = 1e-8, label = paste(label, "rows sum to one"))
          }
        }

        positive <- positive_aux(fit[["aux"]])
        if (length(positive)) {
          expect_true(all(positive > 0),
                      label = paste(label, "positive nuisance parameters"))
        }

        table <- suppressWarnings(diagnose(fit))[["table"]]
        rhat <- table[["rhat"]][is.finite(table[["rhat"]])]
        expect_true(all(rhat > 0.9), label = paste(label, "R-hat in range"))
        for (col in c("ess_bulk", "ess_tail")) {
          ess <- table[[col]][!is.na(table[[col]])]
          expect_true(all(is.finite(ess) & ess > 0),
                      label = paste(label, col, "positive and finite"))
        }
      }
    }
  }

  # The cells the fitter refused, which is a list someone has looked at: a
  # drawn coding is refused wherever the family's leaf target is not quadratic
  # (the two `_ls` families through their scale predictor), and prior weights
  # by the two Dirichlet process families. A family that changes its target
  # form, or starts or stops taking weights, moves a name here.
  refused <- sort(ls(matrix_skipped))
  expected <- sort(c(
    paste(c("gaussian_ls", "poisson", "negbin", "zi_poisson", "zi_negbin",
            "Gamma", "Gamma_ls", "Beta", "ordbeta", "tweedie",
            "ordinal cloglog", "weibull_aft", "ph", "custom"), "drawn"),
    paste(c("dpm", "dpm_aft"), "weights")))
  expect_setequal(refused, expected)
})

test_that("a family with a reporting chart keeps it under every structural wrapper", {
  skip_on_cran()

  # Three families record their draws in a chart other than the sampler's:
  # `ordinal()` centers the predictor over the fitted sample, and `dpm()` and
  # `dpm_aft()` put the mixture at mean zero so that the predictor is the
  # conditional mean. The level of the recorded predictor is then an identified
  # quantity and moves little; in the sampler's own chart it is the coordinate
  # the likelihood does not pin, and it follows `center`, the raw mixture mean,
  # at a correlation near minus one with the spread of `center` itself. That
  # is what a `bcf()` fit reported before the wrapper forwarded the shift.
  d <- sim_matrix(n = 500L, seed = 37L)
  chains <- 2L
  draws <- 200L
  control <- bartisan_control(num_trees = 20L, num_burn = 200L,
                              num_draws = draws, chains = chains,
                              gate = "hard", verbose = FALSE)

  chart_cases <- list(
    list(name = "ordinal logit", family = ordinal(), response = "yord"),
    list(name = "ordinal probit", family = ordinal("probit"),
         response = "yord"),
    list(name = "dpm", family = dpm(), response = "ydpm"),
    list(name = "dpm_aft", family = dpm_aft(), response = "cbind(t_dpma, e_dpma)")
  )

  for (case in chart_cases) {
    for (structure in matrix_structures) {
      fit <- matrix_fit(case, structure, d, control)

      if (is.null(fit)) {
        next
      }

      label <- paste(case$name, "under", structure)
      eta <- stats::predict(fit, type = "link", draws = TRUE)
      level <- rowMeans(eta)

      if (startsWith(case$name, "ordinal")) {
        expect_lt(max(abs(level)), 1e-8, label = paste(label, "centered"))
      }
      else {
        center <- fit[["aux"]][, "center"]
        expect_lt(stats::sd(level), 0.5 * stats::sd(center),
                  label = paste(label, "level against the raw center"))
        expect_lt(abs(stats::cor(level, center)), 0.8,
                  label = paste(label, "level follows the raw center"))
      }
    }
  }
})

# Invariants the documentation states, `_dev/AUDITING.md` § 1: each is a pair of
# routes to one answer, needs no known truth, and is written down as a claim
# already. The augmentation claim is the model for these and sits above.

invariant_data <- function(n = 600L, seed = 43L) {
  set.seed(seed)
  d <- as.data.frame(matrix(stats::runif(n * 3L), n, 3L))
  names(d) <- paste0("x", 1:3)
  d$z <- stats::rbinom(n, 1L, 0.5)
  f <- 1.2 * d$x1 - 1.2 * d$x2
  d$y <- f - mean(f) + 0.8 * d$z + stats::rnorm(n, 0, 0.5)
  d
}

invariant_control <- function(...) {
  bartisan_control(num_trees = 20L, num_burn = 200L, num_draws = 300L,
                   chains = 1L, gate = "hard", verbose = FALSE, ...)
}

ate_of <- function(fit) {
  mean(attr(estimate_effect(fit, treat = "z", estimand = "ATE"), "draws")[[1L]])
}

test_that("a vc() estimand does not depend on the centering", {
  skip_on_cran()

  # `?vc` says every coefficient and every estimand is identical under any
  # choice of `center`. Two samplers for one posterior, so identical up to
  # Monte Carlo error: the control function moves by the effect times the
  # reference and the estimand does not move at all.
  d <- invariant_data()
  fits <- lapply(list("mean", "zero", 0.25), function(center) {
    set.seed(11)
    bartisan(y ~ x1 + x2 + x3 + vc(z, center = center), data = d,
             family = gaussian(), control = invariant_control())
  })

  ates <- vapply(fits, ate_of, numeric(1L))
  expect_lt(diff(range(ates)), 0.15)
  expect_lt(abs(ates[1L] - 0.8), 0.25)

  # The effect is the same for every unit here, so the per-unit coefficients
  # are compared in level rather than by correlation, which has nothing to
  # correlate when the truth is a constant.
  effects <- lapply(fits, function(fit) coef(fit)[, "z"])
  expect_lt(mean(abs(effects[[1L]] - effects[[2L]])), 0.15)
  expect_lt(mean(abs(effects[[1L]] - effects[[3L]])), 0.15)
})

test_that("an empty scale forest makes gaussian_ls() a Gaussian", {
  skip_on_cran()

  # `?bartisan-families`: `gaussian_ls()` with `~ 1` on its scale is
  # `gaussian()`, to within the difference in the prior on the scale, which is
  # a leaf prior on one constant against the Gaussian family's half-Cauchy.
  d <- invariant_data()
  set.seed(12)
  plain <- bartisan(y ~ x1 + x2 + x3 + z, data = d, family = gaussian(),
                    control = invariant_control())
  set.seed(12)
  empty <- bartisan(list(y ~ x1 + x2 + x3 + z, ~ 1), data = d,
                    family = gaussian_ls(), control = invariant_control())

  mean_plain <- colMeans(plain[["eta"]][[1L]])
  mean_empty <- colMeans(empty[["eta"]][[1L]])
  expect_gt(stats::cor(mean_plain, mean_empty), 0.98)
  expect_lt(mean(abs(mean_plain - mean_empty)), 0.1)

  # The scale forest is one constant per draw, and it is the residual sd.
  sd_empty <- exp(empty[["eta"]][[2L]])
  expect_lt(max(apply(sd_empty, 1L, function(x) diff(range(x)))), 1e-8)
  expect_equal(mean(sd_empty[, 1L]), mean(plain[["aux"]][, "sigma"]),
               tolerance = 0.2)
})

test_that("the drawn coding restricts nothing at two levels", {
  skip_on_cran()

  # `?vc` says the drawn coding restricts nothing with two levels, so the
  # estimand should not move when it is switched off.
  d <- invariant_data()
  set.seed(13)
  fixed <- bartisan(y ~ x1 + x2 + x3 + vc(z), data = d, family = gaussian(),
                    control = invariant_control())
  set.seed(13)
  drawn <- bartisan(y ~ x1 + x2 + x3 + vc(z, center = "estimate"), data = d,
                    family = gaussian(), control = invariant_control())

  expect_lt(abs(ate_of(fixed) - ate_of(drawn)), 0.15)
  expect_lt(mean(abs(coef(fixed)[, "z"] - coef(drawn)[, "z"])), 0.15)
})

test_that("sparsity = FALSE is the uniform split_prior", {
  skip_on_cran()

  # `?bartisan_control`: turning the sparsity prior off recovers a uniform
  # prior over predictors, which is a `split_prior` a caller can write out.
  # Both hold the splitting proportions fixed, so with one seed the two fits
  # are the same chain.
  d <- invariant_data(n = 300L)
  set.seed(14)
  off <- bartisan(y ~ x1 + x2 + x3 + z, data = d, family = gaussian(),
                  control = invariant_control(sparsity = FALSE))
  set.seed(14)
  uniform <- bartisan(y ~ x1 + x2 + x3 + z, data = d, family = gaussian(),
                      control = invariant_control(
                        split_prior = c(x1 = 1, x2 = 1, x3 = 1, z = 1)))

  expect_equal(off[["eta"]], uniform[["eta"]], tolerance = 1e-10)
  expect_equal(off[["counts"]], uniform[["counts"]], tolerance = 1e-10)
  expect_equal(off[["loglik"]], uniform[["loglik"]], tolerance = 1e-10)
})
