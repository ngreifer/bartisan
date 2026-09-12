# Causal estimands by g-computation, and the two ways of getting them wrong that
# the design exists to prevent: reading a link-scale contrast as a marginal
# effect, and contrasting before averaging.

sim_effect <- function(n = 200L, seed = 1L, binary = FALSE, levels = 2L) {
  set.seed(seed)
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n),
                  g = factor(sample(c("a", "b"), n, replace = TRUE)))

  d$z <- {
    if (levels == 2L) {
      stats::rbinom(n, 1L, stats::plogis(-0.3 + 1.2 * d$x1))
    }
    else {
      factor(sample(c("low", "mid", "high"), n, replace = TRUE),
             levels = c("low", "mid", "high"))
    }
  }

  lin <- d$x1 + as.numeric(d$z != d$z[1L]) * 0

  eff <- if (levels == 2L) as.numeric(d$z) * (1 + d$x2) else
    (as.integer(d$z) - 1L) * (0.8 + d$x2)

  d$y <- {
    if (binary) stats::rbinom(n, 1L, stats::plogis(-0.5 + d$x1 + eff))
    else d$x1 + eff + stats::rnorm(n)
  }

  d
}

fit_effect <- function(d, binary = FALSE) {
  suppressMessages(suppressWarnings(
    bcf(y ~ x1 + x2 + g, treatment = ~ z, data = d,
        family = if (binary) stats::binomial() else stats::gaussian(),
        control = quick_control(num_trees = 10L, num_burn = 50L,
                                num_draws = 100L))))
}

test_that("the estimands differ only in which units they average over", {
  d <- sim_effect()
  fit <- fit_effect(d)

  ate <- estimate_effect(fit)
  att <- estimate_effect(fit, estimand = "ATT")
  atc <- estimate_effect(fit, estimand = "ATC")

  expect_s3_class(ate, "bartisan_effect")
  expect_identical(nrow(ate), 1L)
  expect_identical(attr(ate, "n_units"), nrow(d))
  expect_identical(attr(att, "n_units"), sum(d$z == 1))
  expect_identical(attr(atc, "n_units"), sum(d$z == 0))

  # The ATE is a weighted average of the other two, with the group sizes as
  # weights, which is the arithmetic the estimands are defined by.
  pooled <- (att[["estimate"]] * sum(d$z == 1) +
               atc[["estimate"]] * sum(d$z == 0)) / nrow(d)
  expect_equal(ate[["estimate"]], pooled, tolerance = 1e-8)
})

test_that("on an identity link the ATT is the mean of the treated units' CATEs", {
  d <- sim_effect(seed = 2L)
  fit <- fit_effect(d)

  att <- estimate_effect(fit, estimand = "ATT")
  cate <- estimate_effect(fit, estimand = "CATE")
  treated <- cate[d$z[cate[["unit"]]] == 1L, ]

  expect_equal(att[["estimate"]], mean(treated[["estimate"]]),
               tolerance = 1e-8)
})

test_that("a marginal odds ratio is not the average conditional odds ratio", {
  d <- sim_effect(seed = 3L, binary = TRUE)
  fit <- fit_effect(d, binary = TRUE)

  marginal <- estimate_effect(fit, comparison = "or")
  conditional <- estimate_effect(fit, estimand = "CATE", comparison = "or")

  # This is the whole reason the estimand is computed on the response scale and
  # averaged before it is contrasted. If these two ever come out equal, either
  # the data has no effect heterogeneity or the averaging order has been
  # silently swapped.
  expect_false(isTRUE(all.equal(marginal[["estimate"]],
                                mean(conditional[["estimate"]]),
                                tolerance = 1e-3)))

  # And the difference contrast is the one where the two orders do agree.
  md <- estimate_effect(fit, comparison = "difference")
  cd <- estimate_effect(fit, estimand = "CATE", comparison = "difference")
  expect_equal(md[["estimate"]], mean(cd[["estimate"]]), tolerance = 1e-8)
})

test_that("the estimands agree with marginaleffects", {
  skip_if_not_installed("marginaleffects")
  skip_on_cran()

  d <- sim_effect(seed = 4L, binary = TRUE)
  fit <- fit_effect(d, binary = TRUE)

  # Both packages summarize the same draws; they differ only in where they put
  # the center, and *marginaleffects* uses the median by default.
  old <- options(marginaleffects_posterior_center = mean)
  on.exit(options(old), add = TRUE)

  ours <- estimate_effect(fit)
  theirs <- marginaleffects::avg_comparisons(fit, variables = "z")
  expect_equal(ours[["estimate"]], theirs[["estimate"]], tolerance = 1e-10)

  ours_by <- estimate_effect(fit, by = ~ g)
  theirs_by <- marginaleffects::avg_comparisons(fit, variables = "z", by = "g")
  expect_equal(ours_by[["estimate"]],
               as.data.frame(theirs_by)[["estimate"]], tolerance = 1e-10)

  # `"or"` summarizes the odds ratio and `"lnor"` its logarithm, so it is the
  # log scale that lines up with what *marginaleffects* reports as `"lnor"`.
  ours_ln <- estimate_effect(fit, estimand = "CATE", comparison = "lnor")
  theirs_ln <- marginaleffects::comparisons(fit, variables = "z",
                                            comparison = "lnor")
  expect_equal(ours_ln[["estimate"]],
               as.data.frame(theirs_ln)[["estimate"]], tolerance = 1e-10)
})

test_that("the potential outcomes are reported for the units averaged over", {
  skip_if_not_installed("marginaleffects")

  d <- sim_effect(seed = 5L, binary = TRUE)
  fit <- fit_effect(d, binary = TRUE)
  old <- options(marginaleffects_posterior_center = mean)
  on.exit(options(old), add = TRUE)

  att <- estimate_effect(fit, estimand = "ATT")
  po <- attr(att, "potential_outcomes")

  expect_identical(nrow(po), 2L)
  # The difference of the two potential outcomes is the effect itself.
  expect_equal(diff(po[["estimate"]]), att[["estimate"]], tolerance = 1e-8)
})

test_that("more than two treatment levels gives every pairwise contrast", {
  d <- sim_effect(seed = 6L, levels = 3L)
  fit <- fit_effect(d)

  ate <- estimate_effect(fit)
  expect_identical(nrow(ate), 3L)
  expect_setequal(ate[["contrast"]],
                  c("Y[mid] - Y[low]", "Y[high] - Y[low]",
                    "Y[high] - Y[mid]"))

  # `focal` has no default with three levels, since which group is the treated
  # one is not something the data says.
  expect_error(estimate_effect(fit, estimand = "ATT"), "must name the focal")

  att <- estimate_effect(fit, estimand = "ATT", focal = "high")
  expect_identical(attr(att, "n_units"), sum(d$z == "high"))
  expect_identical(nrow(att), 3L)

  # Everything is computed; `contrasts` is display only.
  shown <- printed(att)
  expect_true(any(grepl("Y[high] - Y[low]", shown, fixed = TRUE)))
  expect_false(any(grepl("Y[mid] - Y[low]", shown, fixed = TRUE)))

  all_shown <- printed(att, contrasts = "all")
  expect_true(any(grepl("Y[mid] - Y[low]", all_shown, fixed = TRUE)))

  expect_error(print(att, contrasts = "nope - nope"), "does not have")
})

test_that("a continuous treatment is refused, and says what to use instead", {
  d <- sim_effect(seed = 7L)
  d$w <- stats::runif(nrow(d))
  d$y <- d$x1 + d$w + stats::rnorm(nrow(d))

  fit <- suppressMessages(suppressWarnings(
    bartisan(y ~ x1 + x2 + w, d, family = stats::gaussian(),
             control = quick_control())))

  expect_error(estimate_effect(fit, treatment = "w"), "is continuous")
})

test_that("a fit without a named treatment asks for one", {
  d <- sim_effect(seed = 8L)
  fit <- suppressMessages(suppressWarnings(
    bartisan(y ~ x1 + x2 + z, d, family = stats::gaussian(),
             control = quick_control())))

  expect_error(estimate_effect(fit), "must name the treatment")

  # And works once named, since g-computation needs only a column to intervene
  # on and not the `bcf()` structure.
  eff <- estimate_effect(fit, treatment = "z")
  expect_s3_class(eff, "bartisan_effect")
  expect_identical(nrow(eff), 1L)
})

test_that("the interval type is what it says it is", {
  d <- sim_effect(seed = 9L)
  fit <- fit_effect(d)

  eti <- estimate_effect(fit, interval = "eti")
  hpdi <- estimate_effect(fit, interval = "hpdi")

  # The HPD interval is the shortest one at that level, so it cannot be wider.
  expect_lte(hpdi[["upper"]] - hpdi[["lower"]],
             eti[["upper"]] - eti[["lower"]] + 1e-8)

  # Matched against the collapsed output, since cli wraps the note to the
  # console width and the phrase can land across two lines.
  expect_match(printed_text(eti), "equal-tailed credible interval")
  expect_match(printed_text(hpdi), "highest posterior density interval")
})

test_that("an odds ratio is refused when the response is not a probability", {
  d <- sim_effect(seed = 10L)
  fit <- fit_effect(d)

  expect_error(estimate_effect(fit, comparison = "or"),
               "needs the response to be a probability")
})

test_that("the plot argument and the plot method draw the same thing", {
  skip_if_not_installed("ggplot2")

  d <- sim_effect(seed = 11L)
  fit <- fit_effect(d)

  for (spec in list(list(estimand = "CATE"), list(by = ~ g),
                    list(estimand = "ATE"))) {
    eff <- do.call(estimate_effect, c(list(fit), spec))
    from_method <- plot(eff)
    from_arg <- do.call(estimate_effect, c(list(fit), spec, list(plot = TRUE)))

    expect_s3_class(from_arg, "ggplot")
    expect_equal(from_arg[["data"]], from_method[["data"]])
    expect_equal(from_arg[["labels"]], from_method[["labels"]])
  }
})

test_that("print() reports the effect, its potential outcomes and its spread", {
  skip_if_not_installed("ggplot2")

  d <- sim_effect(seed = 12L)
  fit <- fit_effect(d)

  eff <- estimate_effect(fit)
  out <- printed(eff)

  # The two quantities the contrast is a contrast of, printed rather than left
  # in an attribute for the reader to find.
  expect_true(any(grepl("Average potential outcomes", out, fixed = TRUE)))
  expect_true(any(grepl("Y[0]", out, fixed = TRUE)))
  expect_true(any(grepl("Y[1]", out, fixed = TRUE)))
  expect_false(any(grepl("Average potential outcomes",
                         printed(eff, potential_outcomes = FALSE),
                         fixed = TRUE)))

  # A CATE object has one row per unit, so printing it whole is unusable; the
  # spread is what it reports instead, and the rows are still in the object.
  cate <- estimate_effect(fit, estimand = "CATE")
  expect_identical(nrow(cate), nrow(d))
  spread <- printed(cate)
  expect_match(printed_text(cate), "Quartiles of the per-unit estimates")
  expect_true(any(grepl("median", spread, fixed = TRUE)))
  # One row printed, not one per unit, which is the whole point of the change.
  expect_length(grep("Y\\[1\\] - Y\\[0\\]", spread), 1L)

  # `summary()` on a bcf fit is the same summary of the forests it is on any
  # other fit, and says where the effect is reported instead.
  s <- summary(fit)
  expect_s3_class(s, "summary.bartisan_fit")
  expect_identical(s[["treatment"]], "z")
  expect_match(printed_text(s), "estimate_effect")

  # And on a fit with no treatment it says nothing of the kind.
  plain <- suppressMessages(suppressWarnings(
    bartisan(y ~ x1 + x2, d, family = stats::gaussian(),
             control = quick_control())))
  expect_null(summary(plain)[["treatment"]])
  expect_false(grepl("estimate_effect", printed_text(summary(plain))))

  # `plot()` on the fit is the CATE forest, which is `plot()` on that object.
  expect_equal(plot(fit)[["data"]],
               plot(estimate_effect(fit, estimand = "CATE"))[["data"]])
  expect_true(any(grepl("estimate_effect", printed(fit))))
})

test_that("the contrast label names the quantity, not just the levels", {
  d <- sim_effect(seed = 18L, binary = TRUE)
  fit <- fit_effect(d, binary = TRUE)

  # `1 / 0` cannot distinguish a ratio from a log odds ratio, and with numeric
  # levels it reads as arithmetic on the numbers. The label carries the
  # quantity, and the level sits inside a symbol so it cannot be misread.
  expect_identical(estimate_effect(fit)[["contrast"]], "Y[1] - Y[0]")
  expect_identical(estimate_effect(fit, comparison = "ratio")[["contrast"]],
                   "Y[1] / Y[0]")
  expect_identical(estimate_effect(fit, comparison = "lnratio")[["contrast"]],
                   "log(Y[1] / Y[0])")
  expect_identical(estimate_effect(fit, comparison = "or")[["contrast"]],
                   "O(Y[1]) / O(Y[0])")
  expect_identical(estimate_effect(fit, comparison = "lnor")[["contrast"]],
                   "log(O(Y[1]) / O(Y[0]))")

  # And the legend that makes the symbol readable is printed, with the odds
  # clause only where an odds appears.
  expect_match(printed_text(estimate_effect(fit)),
               "Y\\[a\\] is the average response")
  expect_match(printed_text(estimate_effect(fit, comparison = "lnor")),
               "is the odds")
  expect_false(grepl("is the odds", printed_text(estimate_effect(fit))))

  # A factor treatment puts its level names in, so no legend entry is needed
  # per contrast however many there are.
  d3 <- sim_effect(seed = 19L, levels = 3L)
  f3 <- fit_effect(d3)
  expect_setequal(estimate_effect(f3)[["contrast"]],
                  c("Y[mid] - Y[low]", "Y[high] - Y[low]",
                    "Y[high] - Y[mid]"))
})

test_that("the marginal effect is drawn beside the units, not behind them", {
  skip_if_not_installed("ggplot2")

  d <- sim_effect(seed = 20L)
  fit <- fit_effect(d)
  p <- plot(estimate_effect(fit, estimand = "CATE"))

  # A band across the panel put the marginal effect under every conditional
  # interval, which is the one thing a reader wants to locate. It is now its own
  # interval past the right edge, so there is no `geom_rect` and there is a
  # layer whose data holds exactly one row per contrast.
  geoms <- vapply(p[["layers"]], function(l) class(l[["geom"]])[1L],
                  character(1L))
  expect_false("GeomRect" %in% geoms)

  marginal_layer <- vapply(p[["layers"]], function(l) {
    is.data.frame(l[["data"]]) && nrow(l[["data"]]) == 1L &&
      "rank" %in% names(l[["data"]])
  }, logical(1L))
  expect_true(any(marginal_layer))

  # And it sits to the right of every unit.
  marg <- p[["layers"]][[which(marginal_layer)[1L]]][["data"]]
  expect_gt(marg[["rank"]][1L], nrow(d))
})

test_that("intervening on the treatment leaves the propensity score alone", {
  d <- sim_effect(seed = 13L)
  fit <- fit_effect(d)

  # G-computation is only the right intervention if the score, which is a
  # function of the covariates, is carried through rather than recomputed from
  # the assigned treatment.
  frame <- fit[["model"]]
  s0 <- bartisan:::bcf_newdata_score(fit, transform(frame, z = 0L))
  s1 <- bartisan:::bcf_newdata_score(fit, transform(frame, z = 1L))

  expect_identical(s0, s1)
})

test_that("which level is treated follows WeightIt's rules", {
  tl <- function(levs, z) bartisan:::treated_level(levs, z)
  deduced <- function(x) !is.null(attr(x, "deduced"))

  # A logical treatment is treated where it is true.
  expect_identical(as.character(tl(c(FALSE, TRUE), c(TRUE, FALSE))), "TRUE")

  # A numeric one is untreated at zero, whatever the other level is.
  expect_equal(tl(c(0, 1), c(0L, 1L)), 1, ignore_attr = TRUE)
  expect_equal(tl(c(0, 5), c(0L, 5L)), 5, ignore_attr = TRUE)
  expect_true(deduced(tl(c(0, 1), c(0L, 1L))))

  # Levels that only look numeric are read the same way.
  expect_identical(as.character(tl(c("0", "1"), factor(c("0", "1")))), "1")

  # Otherwise the conventional names decide it, from either side.
  expect_identical(as.character(tl(c("control", "trt"),
                                   factor(c("control", "trt")))), "trt")
  expect_identical(as.character(tl(c("a", "treated"),
                                   factor(c("a", "treated")))), "treated")

  # And when nothing in the values says which is which, the second level is
  # taken and the guess is flagged as a guess.
  expect_false(deduced(tl(c("hi", "lo"), factor(c("hi", "lo")))))
  expect_false(deduced(tl(c(1, 2), c(1L, 2L))))
})

test_that("the treated level is guessed quietly when the values settle it", {
  d <- sim_effect(seed = 14L)
  fit <- fit_effect(d)

  # 0 and 1, so there is nothing to guess and nothing to say.
  expect_silent(estimate_effect(fit, estimand = "ATT"))

  att <- estimate_effect(fit, estimand = "ATT")
  atc <- estimate_effect(fit, estimand = "ATC")

  expect_identical(as.character(attr(att, "focal")), "1")
  expect_identical(as.character(attr(atc, "focal")), "0")
  expect_identical(attr(att, "n_units"), sum(d$z == 1))
  expect_identical(attr(atc, "n_units"), sum(d$z == 0))
})

test_that("a guessed treated level is announced", {
  d <- sim_effect(seed = 15L)
  # Level names that say nothing, so the level order is all there is to go on.
  # A numeric treatment cannot make this case, since `bcf()` reads a numeric one
  # that is not 0/1 as continuous and refuses it.
  d$z <- factor(d$z, levels = c(0L, 1L), labels = c("alpha", "beta"))

  fit <- suppressMessages(suppressWarnings(
    bcf(y ~ x1 + x2 + g, treatment = ~ z, data = d,
        family = stats::gaussian(),
        control = quick_control(num_trees = 10L, num_burn = 50L,
                                num_draws = 100L))))

  expect_message(estimate_effect(fit, estimand = "ATT"),
                 "is the treated level")
  expect_message(estimate_effect(fit, estimand = "ATC"),
                 "is the control level")

  # The ATE needs no focal group to be right, so it says nothing.
  expect_silent(estimate_effect(fit))

  # And naming it silences the message rather than changing anything else.
  named <- suppressMessages(estimate_effect(fit, estimand = "ATT"))
  expect_silent(quiet <- estimate_effect(fit, estimand = "ATT",
                                         focal = "beta"))
  expect_equal(quiet[["estimate"]], named[["estimate"]], tolerance = 1e-10)
})

test_that("with more than two levels ATT and ATC are the same estimand", {
  d <- sim_effect(seed = 16L, levels = 3L)
  fit <- fit_effect(d)

  # Both required to name a focal group, and both then average over it.
  expect_error(estimate_effect(fit, estimand = "ATT"), "must name the focal")
  expect_error(estimate_effect(fit, estimand = "ATC"), "must name the focal")

  att <- estimate_effect(fit, estimand = "ATT", focal = "mid")
  atc <- estimate_effect(fit, estimand = "ATC", focal = "mid")

  expect_identical(attr(att, "n_units"), sum(d$z == "mid"))
  expect_identical(attr(atc, "n_units"), attr(att, "n_units"))
  expect_equal(att[["estimate"]], atc[["estimate"]], tolerance = 1e-12)
})

test_that("a newdata holding one arm is still a contrast", {
  d <- sim_effect(seed = 17L)
  fit <- fit_effect(d)

  # The levels to contrast are a property of the fit, not of the units being
  # averaged over. Reading them from `newdata` made the conditional effect
  # among the treated impossible to ask for, since in that subset the treatment
  # takes one value and a one-valued numeric column reads as continuous.
  cate <- estimate_effect(fit, estimand = "CATE",
                          newdata = subset(d, z == 1))

  expect_identical(nrow(cate), sum(d$z == 1))
  expect_identical(unique(cate[["contrast"]]), "Y[1] - Y[0]")

  # Which is the same set of units the ATT averages over, so on an identity
  # link the two must agree exactly.
  att <- estimate_effect(fit, estimand = "ATT")
  expect_equal(att[["estimate"]], mean(cate[["estimate"]]), tolerance = 1e-8)

  # And the same for a factor treatment, where the assigned column has to carry
  # the model frame's levels rather than the subset's.
  d$zf <- factor(d$z, levels = c(0L, 1L), labels = c("ctrl", "trt"))
  ff <- suppressMessages(suppressWarnings(
    bcf(y ~ x1 + x2 + g, treatment = ~ zf, data = d,
        family = stats::gaussian(),
        control = quick_control(num_trees = 10L, num_burn = 50L,
                                num_draws = 100L))))

  one_arm <- estimate_effect(ff, estimand = "CATE",
                             newdata = subset(d, zf == "trt"))

  expect_identical(nrow(one_arm), sum(d$zf == "trt"))
  expect_identical(unique(one_arm[["contrast"]]), "Y[trt] - Y[ctrl]")
})
