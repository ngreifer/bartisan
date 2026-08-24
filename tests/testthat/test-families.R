test_that("stats family objects, functions and names all resolve", {
  expect_identical(as_genbart_family(gaussian())[["family"]], "gaussian")
  expect_identical(as_genbart_family(binomial)[["family"]], "binomial")
  expect_identical(as_genbart_family("poisson")[["family"]], "poisson")

  expect_identical(as_genbart_family(binomial("probit"))[["link"]], "probit")
  expect_identical(as_genbart_family(Gamma("log"))[["family"]], "Gamma")
})

test_that("genbart-specific families carry their own names and links", {
  expect_identical(as_genbart_family(negbin())[["family"]], "negbin")
  expect_identical(as_genbart_family(ordinal("probit"))[["link"]], "probit")
  expect_identical(as_genbart_family(multinomial())[["family"]], "multinomial")
  expect_identical(as_genbart_family(weibull_aft())[["link"]], "weibull")
  expect_identical(as_genbart_family(lognormal_aft())[["family"]], "aft")
  expect_identical(as_genbart_family(location_scale())[["family"]],
                   "location_scale")
})

test_that("unsupported families and links are rejected with a clear message", {
  # A link the engine does not compile is composed onto the family's native
  # scale instead of being rejected, but only for the families where that
  # composition is defined.
  # The family constructors police their own links, so reaching the check in
  # as_genbart_family() means handing it a family object built elsewhere.
  bad <- structure(list(family = "ordinal", link = "cauchit"), class = "family")
  expect_error(as_genbart_family(bad), "link is not supported")
  expect_error(as_genbart_family(ordinal("cauchit")), "should be one of")
  expect_error(as_genbart_family(stats::inverse.gaussian()),
               "not supported by")
  expect_error(as_genbart_family(1), "must be a family name")
})

test_that("a link the engine does not compile is carried as a composition", {
  f <- as_genbart_family(binomial("cauchit"))
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
  expect_null(as_genbart_family(binomial("probit"))[["custom_link"]])
  expect_null(as_genbart_family(poisson())[["custom_link"]])
})

test_that("a family-specific option survives resolution", {
  expect_identical(as_genbart_family(negbin(theta = 3))[["theta"]], 3)
  expect_identical(as_genbart_family(multinomial(reference = "b"))[["reference"]],
                   "b")
  expect_null(as_genbart_family(multinomial())[["reference"]])
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
