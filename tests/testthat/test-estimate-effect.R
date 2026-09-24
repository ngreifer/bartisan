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
    bcf(y ~ x1 + x2 + g, treat = ~ z, data = d,
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

  expect_error(estimate_effect(fit, treat = "w"), "is continuous")
})

test_that("a fit without a named treatment asks for one", {
  d <- sim_effect(seed = 8L)
  fit <- suppressMessages(suppressWarnings(
    bartisan(y ~ x1 + x2 + z, d, family = stats::gaussian(),
             control = quick_control())))

  expect_error(estimate_effect(fit), "must name the treatment")

  # And works once named, since g-computation needs only a column to intervene
  # on and not the `bcf()` structure.
  eff <- estimate_effect(fit, treat = "z")
  expect_s3_class(eff, "bartisan_effect")
  expect_identical(nrow(eff), 1L)
})

test_that("the interval type is what it says it is", {
  d <- sim_effect(seed = 9L)
  fit <- fit_effect(d)

  eti <- estimate_effect(fit, interval = "eti")
  hpdi <- estimate_effect(fit, interval = "hpdi")

  # The HPD interval is the shortest window of draws at that level, so it cannot
  # be wider than the equal-tailed window of the same number of them. Compared
  # against that rather than against `eti`'s own bounds, which come from
  # `quantile()` and interpolate between order statistics: the two then cover
  # slightly different fractions of the draws, `floor(.95 n)` gaps against
  # `.95 (n - 1)`, and the shortest-window guarantee does not span the
  # difference.
  draws <- sort(attr(eti, "draws")[[1L]])
  n <- length(draws)
  k <- floor(0.95 * n)
  start <- (n - k) %/% 2L + 1L

  expect_lte(hpdi[["upper"]] - hpdi[["lower"]],
             draws[start + k] - draws[start] + 1e-8)

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

test_that("plot() draws each kind of effect", {
  skip_if_not_installed("ggplot2")

  d <- sim_effect(seed = 11L)
  fit <- fit_effect(d)

  for (spec in list(list(estimand = "CATE"), list(by = ~ g),
                    list(estimand = "ATE"))) {
    eff <- do.call(estimate_effect, c(list(fit), spec))
    expect_s3_class(plot(eff), "ggplot")
  }

  # Drawing is the method's job alone; there is no argument that does it.
  expect_error(estimate_effect(fit, plot = TRUE), "unused argument")
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
      "pct" %in% names(l[["data"]])
  }, logical(1L))
  expect_true(any(marginal_layer))

  # And it sits to the right of every unit, which run from 0% to 100% of the
  # ranking.
  marg <- p[["layers"]][[which(marginal_layer)[1L]]][["data"]]
  expect_gt(marg[["pct"]][1L], 100)
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
    bcf(y ~ x1 + x2 + g, treat = ~ z, data = d,
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

test_that("a by formula is evaluated rather than read for the names it mentions", {
  d <- sim_effect(seed = 18L)
  fit <- fit_effect(d)

  # `~ x1 > 0.5` groups on the comparison. Reducing the formula to the name it
  # mentions would group on `x1` itself, one group per distinct value, and a
  # continuous predictor has as many of those as it has units.
  cut <- estimate_effect(fit, by = ~ x1 > 0.5)

  expect_identical(nrow(cut), 2L)
  expect_identical(names(cut)[1L], "x1 > 0.5")
  expect_identical(attr(cut, "by"), "x1 > 0.5")
  expect_identical(cut[["n"]], c(sum(d$x1 <= 0.5), sum(d$x1 > 0.5)))

  # The formula is evaluated in its own environment, so a threshold the data
  # does not carry is looked up where the call was written and not in the
  # frame that happens to be evaluating it.
  local({
    threshold <- 0.25
    at_local <- estimate_effect(fit, by = ~ x1 > threshold)
    expect_identical(at_local[["n"]], c(sum(d$x1 <= 0.25), sum(d$x1 > 0.25)))
  })

  # Naming a column and writing a formula that names it are one grouping.
  by_string <- estimate_effect(fit, by = "g")
  by_formula <- estimate_effect(fit, by = ~ g)

  expect_identical(by_string[["g"]], by_formula[["g"]])
  expect_identical(by_string[["estimate"]], by_formula[["estimate"]])

  # One term, because the result has one column to report it in.
  expect_error(estimate_effect(fit, by = ~ x1 + g), "must name one grouping")
  expect_error(estimate_effect(fit, by = ~ nope), "could not be evaluated")
  expect_error(estimate_effect(fit, by = "nope"), "not a column")
})

# The difference-in-differences model of Souto and Louzada Neto (2025) is a
# varying coefficient on a 0/1 treated-now indicator, so every estimand it
# reports is an average of that coefficient over a chosen set of rows. This
# pins that `estimate_effect()` reaches all of them, since the alternative is
# asking a user to average `coef()` draws themselves.
test_that("the difference-in-differences estimands are reachable", {
  skip_on_cran()

  set.seed(20L)
  n_unit <- 80L
  d <- expand.grid(t = 1:6, id = seq_len(n_unit))
  d$cohort <- factor(c("never", "3", "5")[(d$id %% 3L) + 1L],
                     levels = c("never", "3", "5"))
  g <- c(never = Inf, `3` = 3, `5` = 5)[as.character(d$cohort)]
  d$k <- d$t - g
  d$ever <- as.integer(is.finite(g))
  d$Dit <- as.integer(is.finite(g) & d$t >= g)
  d$x1 <- stats::runif(nrow(d))
  d$y <- 0.3 * d$t + d$x1 + 2 * d$Dit + stats::rnorm(nrow(d), sd = 0.3)

  fit <- suppressMessages(
    bartisan(y ~ cohort + t + x1 + vc(Dit, ~ x1 + t + cohort), data = d,
             family = stats::gaussian(),
             control = quick_control(num_trees = c(20L, 10L), num_burn = 200L,
                                     num_draws = 300L)))

  treated <- d[d$Dit == 1, , drop = FALSE]

  # The ATT is the coefficient averaged over the treated rows, and that is what
  # `estimate_effect()` on those rows returns; the two are the same numbers, so
  # a user need not reach for `coef()` at all.
  att <- estimate_effect(fit, treat = "Dit", newdata = treated)
  by_hand <- mean(colMeans(coef(fit, draws = TRUE)[[1L]])[d$Dit == 1])

  expect_identical(nrow(att), 1L)
  expect_equal(att[["estimate"]], by_hand, tolerance = 1e-8)
  expect_equal(att[["estimate"]], 2, tolerance = 0.25)

  # The group-time effects, and the event study. `k` is not a predictor of the
  # model at all; `by` is evaluated in `newdata`, which is what lets an
  # event-time profile be asked for without the model carrying event time.
  gatt <- estimate_effect(fit, treat = "Dit", newdata = treated, by = ~ cohort)
  expect_identical(nrow(gatt), 2L)
  expect_identical(names(gatt)[1L], "cohort")

  es <- estimate_effect(fit, treat = "Dit", newdata = treated, by = ~ k)
  expect_identical(names(es)[1L], "k")
  expect_setequal(es[["k"]], as.character(sort(unique(treated$k))))
  expect_true(all(abs(es[["estimate"]] - 2) < 0.5))

  # And the per-observation effect, which is the coefficient itself.
  catt <- estimate_effect(fit, treat = "Dit", newdata = treated,
                          estimand = "CATE")
  expect_identical(nrow(catt), nrow(treated))
  expect_equal(catt[["estimate"]],
               unname(colMeans(coef(fit, draws = TRUE)[[1L]])[d$Dit == 1]),
               tolerance = 1e-8)
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
    bcf(y ~ x1 + x2 + g, treat = ~ zf, data = d,
        family = stats::gaussian(),
        control = quick_control(num_trees = 10L, num_burn = 50L,
                                num_draws = 100L))))

  one_arm <- estimate_effect(ff, estimand = "CATE",
                             newdata = subset(d, zf == "trt"))

  expect_identical(nrow(one_arm), sum(d$zf == "trt"))
  expect_identical(unique(one_arm[["contrast"]]), "Y[trt] - Y[ctrl]")
})

# The per-level predictions go to workers when a plan has any. The streams are
# drawn before the branch so both paths use the same ones, which is what keeps a
# family whose prediction simulates from depending on whether a plan was set.
test_that("the per-level predictions are the same in parallel as in sequence", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)

  d <- sim_x(n = 250L, p = 3L, seed = 96L)
  d$g <- factor(sample(c("a", "b", "c"), nrow(d), replace = TRUE))
  d$y <- stats::rbinom(nrow(d), 1L, stats::plogis(d$x1 + (d$g == "c")))

  fit <- bartisan(y ~ ., d, family = stats::binomial(), chains = 2L,
                  control = quick_control(num_trees = 5L, num_burn = 60L,
                                          num_draws = 100L))

  effect_now <- function() {
    set.seed(7L)
    estimate_effect(fit, treat = "g")
  }

  future::plan(future::sequential)
  one <- effect_now()

  started <- tryCatch({
    future::plan(future::multisession, workers = 2L)
    isTRUE(future::nbrOfWorkers() >= 2L)
  }, error = function(e) FALSE)

  skip_if_not(started, "no second worker available")

  many <- effect_now()

  # Three levels means three predictions and three pairwise contrasts.
  expect_identical(nrow(one), 3L)
  expect_equal(as.data.frame(many), as.data.frame(one))
  expect_equal(attr(many, "draws"), attr(one, "draws"))
})

test_that("the conditional-effects plot is laid out on a percentile axis", {
  skip_if_not_installed("ggplot2")

  d <- sim_effect(seed = 11L, binary = TRUE)
  fit <- fit_effect(d, binary = TRUE)
  p <- plot(estimate_effect(fit, estimand = "CATE"))
  b <- ggplot2::ggplot_build(p)

  # The units run from 0% to 100% of the ranking, and the axis says so.
  expect_identical(p$labels$x, "Units, ordered by their conditional effect")
  units <- which(vapply(b$data, function(z) {
    "shape" %in% names(z) && nrow(z) == nrow(d)
  }, logical(1L)))[1L]
  expect_equal(range(b$data[[units]]$x), c(0, 100))

  # The average sits centered between the divider and the panel's right edge,
  # which the scale does not pad.
  panel <- b$layout$panel_params[[1L]]$x.range
  divider <- b$data[[which(vapply(b$data, function(z) {
    "xintercept" %in% names(z)
  }, logical(1L)))[1L]]]$xintercept
  average <- b$data[[length(b$data)]]$x[1L]
  expect_equal(average - divider, panel[2L] - average)
  expect_true("Marginal\neffect" %in% b$layout$panel_params[[1L]]$x$get_labels())
})

test_that("an effect axis names the scale a difference is on", {
  lab <- effect_axis_label

  expect_identical(lab("CATE", "difference", "binomial", "response"),
                   "Conditional effect (difference in probability)")
  expect_identical(lab("ATE", "difference", "binomial", "response"),
                   "Effect (difference in probability)")
  expect_identical(lab("CATE", "difference", "binomial", "link"),
                   "Conditional effect (difference on the link scale)")
  expect_identical(lab("CATE", "difference", "gaussian", "link"),
                   "Conditional effect (difference on the link scale)")

  # The response's own units need no naming, and a ratio already names itself.
  expect_identical(lab("CATE", "difference", "gaussian", "response"),
                   "Conditional effect")
  expect_identical(lab("CATE", "or", "binomial", "response"),
                   "Conditional odds ratio")
  expect_identical(lab("ATE", "ratio", "binomial", "response"), "Ratio")
})

test_that("`marginal = FALSE` leaves the marginal effect out of the plot", {
  skip_if_not_installed("ggplot2")

  d <- sim_effect(seed = 11L, binary = TRUE)
  fit <- fit_effect(d, binary = TRUE)
  eff <- estimate_effect(fit, estimand = "CATE")

  with_it <- ggplot2::ggplot_build(plot(eff))
  without <- ggplot2::ggplot_build(plot(eff, marginal = FALSE))

  # No divider, no interval past the units, and an axis that ends with them.
  has_vline <- function(b) any(vapply(b$data, function(z) {
    "xintercept" %in% names(z)
  }, logical(1L)))
  expect_true(has_vline(with_it))
  expect_false(has_vline(without))
  expect_identical(length(without$data), length(with_it$data) - 3L)
  expect_lt(without$layout$panel_params[[1L]]$x.range[2L], 104)
  expect_false("Marginal\neffect" %in%
                 without$layout$panel_params[[1L]]$x$get_labels())

  # `plot()` on a bcf() fit passes it through.
  from_fit <- ggplot2::ggplot_build(plot(fit, marginal = FALSE))
  expect_false(has_vline(from_fit))

  expect_error(plot(eff, marginal = "no"), "marginal")
})
