test_that("stats family objects, functions and names all resolve", {
  expect_identical(as_bartisan_family(gaussian())[["family"]], "gaussian")
  expect_identical(as_bartisan_family(binomial)[["family"]], "binomial")
  expect_identical(as_bartisan_family("poisson")[["family"]], "poisson")

  expect_identical(as_bartisan_family(binomial("probit"))[["link"]], "probit")
  expect_identical(as_bartisan_family(stats::Gamma("log"))[["family"]], "Gamma")
})

test_that("bartisan-specific families carry their own names and links", {
  expect_identical(as_bartisan_family(negbin())[["family"]], "negbin")
  expect_identical(as_bartisan_family(ordinal("probit"))[["link"]], "probit")
  expect_identical(as_bartisan_family(multinomial())[["family"]], "multinomial")
  expect_identical(as_bartisan_family(weibull_aft())[["link"]], "weibull")
  expect_identical(as_bartisan_family(lognormal_aft())[["family"]], "aft")
  expect_identical(as_bartisan_family(location_scale())[["family"]],
                   "location_scale")
})

test_that("unsupported families and links are rejected with a clear message", {
  # A link the engine does not compile is composed onto the family's native
  # scale instead of being rejected, but only for the families where that
  # composition is defined.
  # The family constructors police their own links, so reaching the check in
  # as_bartisan_family() means handing it a family object built elsewhere.
  bad <- structure(list(family = "ordinal", link = "cauchit"), class = "family")
  expect_error(as_bartisan_family(bad), "link is not supported")
  expect_error(as_bartisan_family(ordinal("cauchit")), "should be one of")
  expect_error(as_bartisan_family(stats::inverse.gaussian()),
               "not supported by")
  expect_error(as_bartisan_family(1), "must be a family name")
})

test_that("a link the engine does not compile is carried as a composition", {
  f <- as_bartisan_family(binomial("cauchit"))
  expect_identical(f[["link"]], "cauchit")
  expect_identical(f[["native_link"]], "logit")
  expect_true(is.function(f[["custom_link"]][["linkinv"]]))

  # The composed map takes the caller's predictor to the log odds, and its
  # derivative is the chain rule applied to the caller's mu.eta.
  parts <- compose_link(f[["custom_link"]], "logit")
  eta <- c(-2, -0.5, 0, 1.5)
  expect_equal(parts[["link_theta"]](eta),
               stats::qlogis(stats::pcauchy(eta)))
  expect_equal(parts[["link_dtheta"]](eta),
               stats::dcauchy(eta) /
                 (stats::pcauchy(eta) * (1 - stats::pcauchy(eta))))

  # The native links stay native, so nothing is composed for them.
  expect_null(as_bartisan_family(binomial("probit"))[["custom_link"]])
  expect_null(as_bartisan_family(poisson())[["custom_link"]])
})

test_that("a family-specific option survives resolution", {
  expect_identical(as_bartisan_family(negbin(theta = 3))[["theta"]], 3)
  expect_identical(as_bartisan_family(multinomial(reference = "b"))[["reference"]],
                   "b")
  expect_null(as_bartisan_family(multinomial())[["reference"]])
})

test_that("a two-level factor response works for binomial, more does not", {
  prepared <- prepare_binomial(factor(c("no", "yes", "no")), rep(1, 3), 3)
  expect_identical(prepared[["y"]], c(0, 1, 0))
  expect_identical(prepared[["levels"]], c("no", "yes"))

  expect_error(prepare_binomial(factor(c("a", "b", "c")), rep(1, 3), 3),
               "exactly two levels")
})

test_that("binomial counts become proportions with trials as weights", {
  y <- cbind(c(2, 5), c(8, 5))
  prepared <- prepare_binomial(y, c(1, 1), 2)

  expect_identical(prepared[["y"]], c(0.2, 0.5))
  expect_identical(prepared[["weights"]], c(10, 10))
})

test_that("an unordered factor for the ordinal family warns but proceeds", {
  expect_warning(out <- prepare_ordered(factor(c("a", "b", "c")), "ordinal"),
                 "unordered factor")
  expect_identical(out[["num_cat"]], 3L)
})

test_that("survival responses are validated", {
  ok <- cbind(c(1, 2, 3), c(1, 0, 1))
  expect_identical(prepare_surv(ok, 3)[["event"]], c(1, 0, 1))

  expect_error(prepare_surv(cbind(c(0, 1), c(1, 1)), 2), "strictly positive")
  expect_error(prepare_surv(cbind(c(1, 2), c(0, 0)), 2), "very observation is censored")
  expect_error(prepare_surv(cbind(c(1, 2), c(1, 3)), 2), "event indicator")
  expect_error(prepare_surv(c(1, 2, 3), 3), "two-column response")
})

test_that("unused response levels are dropped for multinomial", {
  y <- factor(c("a", "b", "a"), levels = c("a", "b", "c"))
  expect_warning(out <- prepare_unordered(y, "multinomial"), "unused")
  expect_identical(out[["num_cat"]], 2L)
  expect_identical(out[["levels"]], c("a", "b"))
})

test_that("`Beta()` is the interior of `ordbeta()` and says so about boundaries", {
  d <- sim_x(n = 200, seed = 72)
  set.seed(172)
  mu <- stats::plogis(0.5 * d$x1)
  d$y <- stats::rbeta(nrow(d), mu * 10, 10 - mu * 10)

  expect_identical(as_bartisan_family(Beta())[["family"]], "beta")
  expect_identical(as_bartisan_family(Beta())[["link"]], "logit")
  expect_identical(as_bartisan_family("Beta")[["link"]], "logit")

  # Only the logit link is compiled. The other two are reached by composition,
  # and this checks that they actually are: listing them as native links made
  # them silently fit the logit model and then back-transform as though they had
  # not, so the fitted means were wrong.
  expect_null(as_bartisan_family(Beta("logit"))[["custom_link"]])

  for (link in c("probit", "cloglog")) {
    family <- as_bartisan_family(Beta(link))
    expect_identical(family[["link"]], link, label = link)
    expect_false(is_null(family[["custom_link"]]), label = link)
    expect_identical(family[["native_link"]], "logit", label = link)
  }

  fit_link <- function(link) {
    set.seed(9)
    bartisan(y ~ ., d, family = Beta(link),
             control = quick_control(num_trees = 5L, num_burn = 10L,
                                     num_draws = 10L))
  }
  expect_false(identical(fit_link("logit")[["eta"]],
                         fit_link("probit")[["eta"]]))
  expect_false(identical(fit_link("logit")[["eta"]],
                         fit_link("cloglog")[["eta"]]))

  # A response at either endpoint has no beta density, so it is an error rather
  # than something to nudge inward, and the message names the family that does
  # model the endpoints.
  for (bad in c(0, 1)) {
    dd <- d
    dd$y[1L] <- bad
    expect_error(bartisan(y ~ ., dd, family = Beta(),
                         control = quick_control()),
                 "strictly between 0 and 1")
  }

  # A fixed precision is respected exactly and reported as a constant.
  fit <- bartisan(y ~ ., d, family = Beta(phi = 7),
                  control = quick_control(num_trees = 5L, num_burn = 10L,
                                          num_draws = 10L))
  expect_true(all(fit[["aux"]][, "phi"] == 7))

  expect_error(Beta(phi = -1))

  # Any other link is composed, as it is for `binomial()`, so `cauchit` is
  # accepted rather than refused; a name that is not a link at all is not.
  expect_false(is_null(as_bartisan_family(Beta("cauchit"))[["custom_link"]]))
  expect_error(as_bartisan_family(Beta("not-a-link")))
})

test_that("the gamma family overrules any link but the log one", {
  # The inverse link -- base R's default, because it is the canonical link for
  # the gamma -- needs a positive mean, and this sampler's additive predictor is
  # unconstrained, so a draw that wanders non-positive has no gamma density at
  # all. Rather than accept that and warn, the link is replaced. `stats::Gamma()`
  # itself is left exactly as base R defines it, so that attaching this package
  # cannot change what `glm()` does.
  expect_false("Gamma" %in% getNamespaceExports("bartisan"))
  expect_identical(stats::Gamma()[["link"]], "inverse")

  for (link in c("inverse", "identity", "sqrt")) {
    expect_message(family <- as_bartisan_family(stats::Gamma(link)),
                   "is ignored")
    expect_identical(family[["link"]], "log")

    # And nothing composed survives, since the caller's link is not honored.
    expect_null(family[["custom_link"]])
  }

  # Naming the log link, or naming the family as this package's own string, is
  # silent: there is nothing to report.
  expect_no_message(family <- as_bartisan_family(stats::Gamma("log")))
  expect_identical(family[["link"]], "log")

  expect_no_message(family <- as_bartisan_family("Gamma"))
  expect_identical(family[["link"]], "log")

  # And the fixed-shape family is gone, so nothing offers it.
  expect_false(exists("Gamma_shape", where = asNamespace("bartisan"),
                      inherits = FALSE))
})

test_that("a composed link whose inverse has a restricted range is reported", {
  d <- sim_x(n = 120, seed = 71)
  d$y <- stats::rgamma(nrow(d), 4, 4 / exp(1 + d$x1))
  ctrl <- quick_control(num_trees = 5L, num_burn = 20L, num_draws = 20L)

  # The gamma is no longer an example of this, because its link is replaced
  # rather than composed; the Poisson identity link is the remaining case.
  d$count <- stats::rpois(nrow(d), 3)
  expect_message(bartisan(count ~ x1 + x2 + x3, d, family = poisson("identity"),
                         control = ctrl),
                 "does not cover the whole additive predictor")

  # A compiled link has nothing composed, and a composed link whose inverse does
  # cover the line is fine, so neither says anything.
  expect_no_message(bartisan(y ~ ., d, family = stats::Gamma("log"),
                            control = ctrl))

  d$bin <- stats::rbinom(nrow(d), 1, 0.4)
  expect_no_message(bartisan(bin ~ x1 + x2 + x3, d, family = binomial("cauchit"),
                            control = ctrl))
})
