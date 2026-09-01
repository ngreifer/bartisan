# Varying coefficients, from the formula through the fit. The first half tests
# the specification layer -- what the marker means, which forests may split on
# what, and the basis each coefficient multiplies -- and the second half fits.

sim_vc <- function(n = 120, seed = 1) {
  set.seed(seed)
  data.frame(y = stats::rnorm(n), x1 = stats::rnorm(n), x2 = stats::rnorm(n),
             z = stats::rbinom(n, 1L, 0.5),
             g = factor(sample(letters[1:3], n, TRUE)))
}

# The pipeline `bartisan()` runs, up to the point the engine takes over. One
# additive predictor, with every predictor group available to it, which is what a
# single-parameter family gives `resolve_vc()`.
spec_of <- function(formula, data, n_param = 1L, labels = NULL) {
  split <- split_random(formula)
  vc_split <- split_vc_terms(split[["fixed"]])
  mf <- stats::model.frame(vc_to_names(formula), data)
  mt <- stats::terms(vc_split[["fixed"]], data = mf)
  design <- build_design(mt, mf)
  groups <- unique(design[["term_labels"]][design[["assign"]]])

  forest_vc <- rep(list(list(specs = vc_split[["vc"]],
                             dot = uses_dot(split[["fixed"]]))), n_param)

  resolve_vc(forest_vc, mf, design,
             matrix(TRUE, nrow = length(groups), ncol = n_param),
             labels = labels)
}

test_that("vc() is a formula marker and refuses to be called", {
  expect_error(vc(1), "cannot be\\s+called on its own")
})

test_that("a vc() covariate leaves the design and stays in the frame", {
  d <- sim_vc()

  split <- split_vc_terms(y ~ x1 + x2 + vc(z))

  # Out of the design: the covariate is what a coefficient multiplies, not
  # something the control function splits on.
  expect_setequal(all.vars(split[["fixed"]]), c("y", "x1", "x2"))
  expect_named(split[["vc"]], "z")

  # In the frame, so it gets the same missing-value handling as any predictor.
  expect_setequal(all.vars(vc_to_names(y ~ x1 + x2 + vc(z))),
                  c("y", "x1", "x2", "z"))
})

test_that("vc() reads its arguments the way any call would", {
  one <- split_vc_terms(y ~ x1 + vc(z))[["vc"]][["z"]]
  expect_null(one[["modifiers"]])
  expect_identical(one[["center"]], "auto")

  positional <- split_vc_terms(y ~ x1 + vc(z, ~ x1))[["vc"]][["z"]]
  named <- split_vc_terms(y ~ x1 + vc(z, modifiers = ~ x1))[["vc"]][["z"]]
  expect_equal(positional[["modifiers"]], named[["modifiers"]])

  expect_identical(
    split_vc_terms(y ~ x1 + vc(z, center = "zero"))[["vc"]][["z"]][["center"]],
    "zero")
  expect_identical(
    split_vc_terms(y ~ x1 + vc(z, center = 3.5))[["vc"]][["z"]][["center"]],
    3.5)
})

test_that("several coefficients can vary, each with its own modifiers", {
  d <- sim_vc()
  spec <- spec_of(y ~ x1 + x2 + vc(z, ~ x1) + vc(g), d)

  # One forest for z, three for g's levels, and the control function.
  expect_identical(spec[["slopes"]], 4L)
  expect_identical(ncol(spec[["masks"]]), 5L)

  # z's forest was restricted to x1; g's levels were not restricted.
  z_at <- 2L
  expect_identical(rownames(spec[["masks"]])[spec[["masks"]][, z_at]], "x1")
  expect_true(all(spec[["masks"]][, 3L]))
})

test_that("`.` silently excludes a varying covariate from the control function", {
  d <- sim_vc()

  # The caller wrote `.`, not `z`, so dropping it is a correction rather than a
  # change of model, and it happens without comment.
  expect_no_warning(spec <- spec_of(y ~ . + vc(z), d))

  control <- spec[["masks"]][, 1L]
  expect_false(control[["z"]])
  expect_true(control[["x1"]])

  # The coefficient's own forest may still split on it: z * f(z) is a
  # nonlinearity in z, and it is identified.
  expect_true(spec[["masks"]][, 2L][["z"]])
})

test_that("naming a varying covariate as a modifier warns and fits as asked", {
  d <- sim_vc()

  expect_warning(spec <- spec_of(y ~ x1 + z + vc(z), d),
                 "not separately identified")

  # Fitted as written, which is what the warning says.
  expect_true(spec[["masks"]][, 1L][["z"]])
})

test_that("a categorical covariate is removed from its own forest quietly", {
  d <- sim_vc()

  # A numeric covariate in its own modifiers is kept: z * f(z) is a nonlinearity
  # in z and is identified. It has to be a design column to be a modifier at
  # all, which is why it is named in the fixed part here, and that is also what
  # trips the control-function rule -- so the two are checked together below.
  expect_warning(numeric_own <- spec_of(y ~ x1 + z + vc(z, ~ z + x1), d))
  expect_true(numeric_own[["masks"]][, 2L][["z"]])

  # Naming the covariate's own name in its modifiers is allowed even when the
  # covariate is not a design column, since `vc()` is what took it out.
  expect_no_error(spec_of(y ~ x1 + x2 + vc(g, ~ g + x1), d))
})

test_that("a modifier that is not a predictor is a typo, not a restriction", {
  d <- sim_vc()

  expect_error(spec_of(y ~ x1 + x2 + vc(z, ~ nope), d),
               "not a predictor")
  expect_error(spec_of(y ~ x1 + x2 + vc(z, ~ x1 + typo), d),
               "not .*predictors")
})

test_that("the two overlap rules are separate", {
  d <- sim_vc()

  # Naming a factor in the fixed part is the control-function rule, which warns;
  # the categorical rule then also removes it from its own forests, quietly.
  expect_warning(spec <- spec_of(y ~ x1 + g + vc(g), d),
                 "not separately identified")

  expect_true(spec[["masks"]][, 1L][["g"]])

  for (j in 2:4) {
    expect_false(spec[["masks"]][, j][["g"]])
  }
})

test_that("the basis is centred, and centring says what it says", {
  d <- sim_vc()

  mean_centred <- spec_of(y ~ x1 + vc(z, center = "mean"), d)
  expect_equal(mean(mean_centred[["basis"]][, 1L]), 0, tolerance = 1e-12)

  at_zero <- spec_of(y ~ x1 + vc(z, center = "zero"), d)
  expect_identical(as.numeric(at_zero[["basis"]][, 1L]), as.numeric(d$z))

  # `"mid"` is the midpoint of the range, so 0.5 for a 0/1 covariate.
  at_mid <- spec_of(y ~ x1 + vc(z, center = "mid"), d)
  expect_setequal(unique(as.numeric(at_mid[["basis"]][, 1L])), c(-0.5, 0.5))
})

test_that("a factor gets one symmetric forest per level", {
  d <- sim_vc()
  spec <- spec_of(y ~ x1 + vc(g), d)

  expect_identical(spec[["slopes"]], 3L)
  expect_identical(colnames(spec[["basis"]]), c("ga", "gb", "gc"))

  # Mean-centred indicators, so the columns sum to zero within a row. That sum
  # is the one function-valued redundancy the coding carries, and it is the same
  # one symmetric multinomial coding has.
  expect_equal(unname(rowSums(spec[["basis"]])), rep(0, nrow(d)),
               tolerance = 1e-12)
  expect_equal(unname(colMeans(spec[["basis"]])), rep(0, 3),
               tolerance = 1e-12)
})

test_that("a factor's center names a level, and numeric references are refused", {
  d <- sim_vc()

  expect_no_error(spec_of(y ~ x1 + vc(g, center = "b"), d))
  expect_error(spec_of(y ~ x1 + vc(g, center = "zero"), d),
               "must be .*mean.* or")
  expect_error(spec_of(y ~ x1 + vc(g, center = 0.5), d),
               "must be .*mean.* or")
})

test_that("forest labels drop the parameter when the family has only one", {
  d <- sim_vc()
  split <- split_vc_terms(y ~ x1 + vc(z) + vc(g))
  mf <- stats::model.frame(vc_to_names(y ~ x1 + vc(z) + vc(g)), d)
  parts <- vc_basis(split[["vc"]], mf)[["parts"]]

  expect_identical(vc_forest_labels("eta", split[["vc"]], parts, TRUE),
                   c("(Intercept)", "z", "ga", "gb", "gc"))
  expect_identical(vc_forest_labels("mean", split[["vc"]], parts, FALSE),
                   c("mean", "mean:z", "mean:ga", "mean:gb", "mean:gc"))
  expect_identical(vc_forest_labels("eta", list(), list(), TRUE), "eta")
})

test_that("vc() must be added with + and must name a bare predictor", {
  expect_error(split_vc_terms(y ~ x1 - vc(z)), "added to the formula")
  expect_error(split_vc_terms(y ~ x1:vc(z)), "added to the formula")
  expect_error(split_vc_terms(y ~ x1 + vc(log(z))), "bare predictor")
  expect_error(split_vc_terms(y ~ x1 + vc(z, x1)), "one-sided formula")
  expect_error(split_vc_terms(y ~ x1 + vc(z, nope = 1)), "does not|unused")
})

test_that("formula shapes that are not vc() survive untouched", {
  # Only `+` is split, so everything else keeps its spelling and the design
  # formula can still be expanded against the data later.
  for (f in list(y ~ ., y ~ . - x1, y ~ x1:x2, y ~ I(x1^2), y ~ x1 * x2)) {
    expect_equal(split_vc_terms(f)[["fixed"]], f)
  }

  keeps <- split_vc_terms(y ~ . - x1 + vc(z))[["fixed"]]
  expect_identical(deparse1(keeps), "y ~ . - x1")
})

test_that("a model of nothing but vc() terms is refused", {
  d <- sim_vc()

  expect_error(bartisan(y ~ vc(z), data = d, family = gaussian(),
                        control = quick_control()),
               "at least one predictor")
})



# ---- fitting -------------------------------------------------------------

sim_effect <- function(n = 600, seed = 3) {
  set.seed(seed)
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                  z = stats::rbinom(n, 1L, 0.5))
  d$tau <- 1 + 0.8 * d$x2
  d$y <- 2 * d$x1 + d$z * d$tau + stats::rnorm(n, sd = 0.5)
  d
}

vc_control <- function(...) {
  bartisan_control(num_trees = 20, num_burn = 250, num_draws = 250,
                   verbose = FALSE, ...)
}

test_that("a varying coefficient is fitted and recovers its function", {
  d <- sim_effect()

  fit <- bartisan(y ~ x1 + x2 + vc(z), data = d, family = gaussian(),
                  control = vc_control())

  expect_identical(fit[["num_forest"]], 2L)
  expect_named(fit[["eta"]], c("(Intercept)", "z"))

  expect_gt(cor(colMeans(fit[["eta"]][["z"]]), d$tau), 0.9)
})

test_that("the estimand path and the coefficient path agree", {
  skip_if_not_installed("marginaleffects")

  d <- sim_effect(seed = 4)
  fit <- bartisan(y ~ x1 + x2 + vc(z), data = d, family = gaussian(),
                  control = vc_control())

  # Every estimand reaches this through `predict()` on modified `newdata`, which
  # recomputes the basis and contrasts the combination. For an identity link that
  # has to be *exactly* the coefficient, and if it is not, one of the two is
  # wrong.
  treated <- predict(fit, newdata = transform(d, z = 1), type = "link",
                     summary = FALSE, draws = TRUE)
  control <- predict(fit, newdata = transform(d, z = 0), type = "link",
                     summary = FALSE, draws = TRUE)

  expect_equal(mean(rowMeans(treated - control)),
               mean(colMeans(fit[["eta"]][["z"]])), tolerance = 1e-10)

  # And marginaleffects gets the same answer through its own grid, which builds
  # and averages the counterfactuals its own way, so it agrees to its own
  # precision rather than to the last bit.
  ate <- marginaleffects::avg_comparisons(fit, variables = "z")

  expect_equal(ate[["estimate"]], mean(colMeans(fit[["eta"]][["z"]])),
               tolerance = 1e-3)
})

test_that("predict rebuilds the basis and agrees with the fit", {
  d <- sim_effect(seed = 5)
  fit <- bartisan(y ~ x1 + x2 + vc(z), data = d, family = gaussian(),
                  control = vc_control())

  expect_equal(as.numeric(predict(fit, newdata = d)),
               as.numeric(fitted(fit)), tolerance = 1e-7)

  # The centring is the fit's, not the new data's. Re-centring on `newdata` would
  # move the control function's reference between fitting and prediction, so the
  # same covariate value would predict two different things -- which is exactly
  # what the counterfactual grids `avg_comparisons()` builds would trip over,
  # since every row of one holds the treatment at a single value.
  all_treated <- transform(d, z = 1)
  expect_equal(as.numeric(predict(fit, newdata = all_treated)),
               as.numeric(fitted(fit)) +
                 (1 - d$z) * coef(fit)[, "z"],
               tolerance = 1e-6)
})

test_that("a continuous covariate may modify its own coefficient", {
  # The effect is then no longer linear in the covariate: the slope moves across
  # its range. With the truth y = 2 x1 + z^2, the slope at z is z, which a model
  # whose coefficient cannot see z has no way to produce.
  set.seed(6)
  n <- 900
  d <- data.frame(x1 = stats::rnorm(n), z = stats::runif(n, -2, 2))
  d$truth <- 2 * d$x1 + d$z^2
  d$y <- d$truth + stats::rnorm(n, sd = 0.4)

  linear <- bartisan(y ~ x1 + vc(z), data = d, family = gaussian(),
                     control = vc_control())
  curved <- suppressWarnings(
    bartisan(y ~ x1 + z + vc(z, ~ z + x1), data = d, family = gaussian(),
             control = vc_control()))

  expect_lt(sqrt(mean((fitted(curved) - d$truth)^2)),
            sqrt(mean((fitted(linear) - d$truth)^2)) / 3)

  # And the coefficient traces the slope of z^2, which is z.
  at <- stats::approx(d$z[order(d$z)],
                      coef(curved)[order(d$z), "z"],
                      xout = c(-1, 0, 1))$y

  expect_equal(at, c(-1, 0, 1), tolerance = 0.35)
})

test_that("centring moves the control function and leaves the effect alone", {
  d <- sim_effect(seed = 7)

  at_mean <- bartisan(y ~ x1 + x2 + vc(z, center = "mean"), data = d,
                      family = gaussian(), control = vc_control())
  at_zero <- bartisan(y ~ x1 + x2 + vc(z, center = "zero"), data = d,
                      family = gaussian(), control = vc_control())

  # Both estimate the same coefficient function.
  expect_gt(cor(coef(at_mean)[, "z"], coef(at_zero)[, "z"]), 0.85)

  # The control function differs by the effect times the reference, since
  # f0_mean = f0_zero + mean(z) * f1.
  gap <- colMeans(at_mean[["eta"]][[1L]]) - colMeans(at_zero[["eta"]][[1L]])
  expect_gt(mean(gap), 0)
})

test_that("the default reference depends on the covariate", {
  d <- sim_effect(seed = 8)
  d$w <- d$z + 50

  binary <- spec_of(y ~ x1 + vc(z), d)
  shifted <- spec_of(y ~ x1 + vc(w), d)

  # Zero is a value a 0/1 covariate takes, and the control function there is the
  # surface among the untreated. Zero is nowhere near `w`, so its control
  # function there would be an extrapolation.
  expect_identical(binary[["parts"]][[1L]][["center"]], 0)
  expect_equal(shifted[["parts"]][[1L]][["center"]], mean(d$w))
})

test_that("a factor gets one forest per level and coef recenters them", {
  set.seed(9)
  n <- 700
  d <- data.frame(x1 = stats::rnorm(n),
                  g = factor(sample(c("a", "b", "c"), n, TRUE)))
  effect <- c(a = -1, b = 0, c = 1)[as.character(d$g)]
  d$y <- 2 * d$x1 + effect + stats::rnorm(n, sd = 0.4)

  fit <- bartisan(y ~ x1 + vc(g), data = d, family = gaussian(),
                  control = vc_control())

  expect_identical(fit[["num_forest"]], 4L)
  expect_named(fit[["eta"]], c("(Intercept)", "ga", "gb", "gc"))

  cf <- coef(fit)
  expect_identical(colnames(cf), c("ga", "gb", "gc"))

  # Recentred to sum to zero across levels, which is what makes each one the
  # deviation it is reported as.
  expect_equal(unname(rowSums(cf)), rep(0, n), tolerance = 1e-10)

  # And the deviations order the way the truth does.
  expect_lt(mean(cf[, "ga"]), mean(cf[, "gc"]))
})

test_that("a group intercept reaches the control function and not the effect", {
  set.seed(10)
  n <- 600
  d <- data.frame(x1 = stats::rnorm(n), z = stats::rbinom(n, 1L, 0.5),
                  g = factor(sample(letters[1:5], n, TRUE)))
  d$y <- d$x1 + as.numeric(d$g) * 0.5 + 1.5 * d$z + stats::rnorm(n, sd = 0.4)

  fit <- bartisan(y ~ x1 + (1 | g) + vc(z), data = d, family = gaussian(),
                  control = vc_control())

  # A group effect on a coefficient's forest is a random slope, and
  # `split_random()` refuses `(x | g)` in as many words -- so producing one here
  # without being asked would contradict that refusal silently.
  expect_false(all(fit[["ranef"]][[1L]] == 0))
  expect_true(all(fit[["ranef"]][[2L]] == 0))
})

test_that("an augmented family carries varying coefficients too", {
  # The augmented families are built a second time, from the name and the
  # options, by the rewriting that turns them into Gaussian ones. A wrapper put
  # on the first family is discarded there, so an augmented family reached the
  # sampler claiming one additive predictor while the rest of the fit expected
  # two. Every family that has an augmented counterpart is worth checking.
  set.seed(13)
  n <- 600
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                  z = stats::rbinom(n, 1L, 0.5))
  d$y <- stats::rbinom(n, 1L, stats::plogis(d$x1 + d$z * (1 + d$x2)))

  for (link in c("logit", "probit")) {
    fit <- bartisan(y ~ x1 + x2 + vc(z), data = d,
                    family = stats::binomial(link), control = vc_control())

    expect_identical(fit[["num_forest"]], 2L)
    expect_named(fit[["eta"]], c("(Intercept)", "z"))
    expect_gt(mean(coef(fit)[, "z"]), 0)
  }

  # And a family with a nuisance parameter, which the decorator has to pass
  # through to the wrapped family rather than claim as its own.
  d$count <- stats::rpois(n, exp(0.5 + 0.3 * d$z))
  expect_no_error(bartisan(count ~ x1 + vc(z), data = d,
                           family = stats::poisson(), control = vc_control()))
})

test_that("one formula gives every additive predictor the same coefficients", {
  d <- sim_effect(seed = 11)

  # A single formula applies to every forest, which is the recycling rule every
  # per-forest argument follows, so a `vc()` term in it reaches each parameter.
  fit <- bartisan(y ~ x1 + x2 + vc(z), data = d, family = location_scale(),
                  control = vc_control())

  expect_identical(fit[["num_forest"]], 4L)
  expect_named(fit[["eta"]], c("mean", "mean:z", "log_sd", "log_sd:z"))
  expect_identical(colnames(coef(fit)), c("mean:z", "log_sd:z"))
})

test_that("coef refuses a model that has no coefficients", {
  d <- sim_effect(seed = 12)
  fit <- bartisan(y ~ x1 + x2, data = d, family = gaussian(),
                  control = vc_control())

  expect_error(coef(fit), "no varying coefficients")
})

# `center = "estimate"`: the coding is drawn rather than fixed, which is the
# parameter expansion of Hahn, Murray and Carvalho (2020, section 5.3).

test_that("a drawn coding reports one coefficient and its own coding draws", {
  d <- sim_effect(seed = 20)

  fit <- bartisan(y ~ x1 + x2 + vc(z, center = "estimate"), data = d,
                  family = gaussian(), control = vc_control())

  # One column and one forest whatever the number of levels, and the coding
  # coefficients are reported as nuisance parameters of their own.
  expect_identical(fit[["num_forest"]], 2L)
  expect_named(fit[["eta"]], c("(Intercept)", "z"))
  expect_true(all(c("b.z.0", "b.z.1") %in% colnames(fit[["aux"]])))

  expect_gt(cor(coef(fit)[, "z"], d$tau), 0.85)
})

test_that("coef under a drawn coding is the identified contrast, not the forest", {
  d <- sim_effect(seed = 21)

  fit <- bartisan(y ~ x1 + x2 + vc(z, center = "estimate"), data = d,
                  family = gaussian(), control = vc_control())

  # The effect is (b_1 - b_0) * f, so the reported coefficient must be the raw
  # forest scaled by the drawn contrast. Getting this wrong is the trap the
  # accessor exists to close.
  raw <- fit[["eta"]][["z"]]
  contrast <- fit[["aux"]][, "b.z.1"] - fit[["aux"]][, "b.z.0"]

  expect_equal(coef(fit)[, "z"], colMeans(raw * contrast),
               ignore_attr = TRUE)
})

test_that("a drawn coding does not depend on which level was written as 1", {
  d <- sim_effect(seed = 22)
  flipped <- transform(d, z = 1 - z)

  effect <- function(form, data) {
    mean(coef(bartisan(form, data = data, family = gaussian(),
                       control = vc_control()))[, 1L])
  }

  # Relabelling the treatment is the same model, so the effect should come back
  # the same size with the sign flipped and the two should sum to zero. The
  # drawn coding is closer to that than a fixed one, which is the whole point of
  # the parameter expansion.
  drawn <- effect(y ~ x1 + x2 + vc(z, center = "estimate"), d) +
    effect(y ~ x1 + x2 + vc(z, center = "estimate"), flipped)
  fixed <- effect(y ~ x1 + x2 + vc(z, center = "zero"), d) +
    effect(y ~ x1 + x2 + vc(z, center = "zero"), flipped)

  expect_lt(abs(drawn), abs(fixed))
})

test_that("a factor's drawn coding shares one forest across levels", {
  d <- sim_vc(n = 300, seed = 23)
  d$y <- d$x1 + c(a = 0, b = 1, c = 2)[as.character(d$g)] * (1 + d$x2) +
    stats::rnorm(300)

  shared <- bartisan(y ~ x1 + x2 + vc(g, center = "estimate"), data = d,
                     family = gaussian(), control = vc_control())
  per_level <- bartisan(y ~ x1 + x2 + vc(g), data = d, family = gaussian(),
                        control = vc_control())

  # One shared forest against one per level: the restriction the drawn coding
  # imposes above two levels, and the reason it stays off by default.
  expect_identical(shared[["num_forest"]], 2L)
  expect_identical(per_level[["num_forest"]], 4L)

  # And the reported contrasts are against the first level, so there are K - 1.
  expect_identical(colnames(coef(shared)), c("gb", "gc"))
  expect_identical(colnames(coef(per_level)), c("ga", "gb", "gc"))
})

test_that("a drawn coding needs a covariate with few enough values", {
  d <- sim_effect(seed = 24)

  expect_error(
    bartisan(y ~ x2 + vc(x1, center = "estimate"), data = d,
             family = gaussian(), control = vc_control()),
    "too many to code")

  d$constant <- 1
  expect_error(
    bartisan(y ~ x1 + vc(constant, center = "estimate"), data = d,
             family = gaussian(), control = vc_control()),
    "at least two values")
})

test_that("an abbreviated center is matched before it is acted on", {
  d <- sim_effect(seed = 25)

  # `match_arg()` completes an abbreviation, so testing the unmatched string
  # would take `"est"` down the fixed-centring path and silently fit the wrong
  # model. It has to reach the drawn coding.
  fit <- bartisan(y ~ x1 + x2 + vc(z, center = "est"), data = d,
                  family = gaussian(), control = vc_control())

  expect_true("b.z.0" %in% colnames(fit[["aux"]]))
})

test_that("the density is right when the coding is drawn", {
  d <- sim_effect(seed = 26)

  fit <- bartisan(y ~ x1 + x2 + vc(z, center = "estimate"), data = d,
                  family = gaussian(), control = vc_control())

  # The stored basis for a drawn coding is a zero placeholder, so a density that
  # reached for it would price every row as though the covariate had no effect.
  # Against the plain normal density at the fitted predictor it must not.
  dens <- predict(fit, type = "density", log = TRUE)
  manual <- stats::dnorm(d$y, fitted(fit), mean(fit[["aux"]][, "sigma"]),
                         log = TRUE)

  expect_lt(abs(sum(dens) - sum(manual)) / abs(sum(manual)), 0.05)
})

test_that("bcf draws the coding for a binary treatment and not for more levels", {
  d <- sim_vc(n = 250, seed = 27)
  d$y <- d$x1 + d$z * (1 + d$x2) + stats::rnorm(250)

  binary <- bcf(y ~ x1 + x2, treatment = ~ z, data = d, family = gaussian(),
                propensity = FALSE, num_trees = 15L, control = vc_control())
  several <- bcf(y ~ x1 + x2, treatment = ~ g, data = d, family = gaussian(),
                 propensity = FALSE, num_trees = 15L, control = vc_control())

  expect_true("b.z.0" %in% colnames(binary[["aux"]]))
  expect_identical(binary[["num_forest"]], 2L)

  # A treatment with more levels keeps the symmetric coding, so there is no
  # coding to report and a forest per level.
  expect_false(any(startsWith(colnames(several[["aux"]]), "b.")))
  expect_identical(several[["num_forest"]], 4L)
})

test_that("a two-level factor treatment sizes its forests from the coding", {
  d <- sim_vc(n = 250, seed = 28)
  d$zf <- factor(ifelse(d$z == 1L, "yes", "no"))
  d$y <- d$x1 + d$z * (1 + d$x2) + stats::rnorm(250)

  # The drawn coding is one forest whatever the type, so the default tree count
  # has to follow the coding rather than the treatment being a factor.
  fit <- bcf(y ~ x1 + x2, treatment = ~ zf, data = d, family = gaussian(),
             propensity = FALSE, num_trees = 15L, control = vc_control())

  expect_identical(fit[["num_forest"]], 2L)
  expect_identical(colnames(coef(fit)), "zf")
})

test_that("an augmented family qualifies for a drawn coding", {
  d <- sim_effect(seed = 29)
  d$y <- stats::rbinom(nrow(d), 1L, stats::plogis(d$x1 + d$z))

  # The guard has to be asked of the family that reaches the sampler, not of the
  # one the formula named: Polya-Gamma and Albert-Chib augmentation both make a
  # binomial's leaf target quadratic, so a logit or probit binomial qualifies
  # even though the unaugmented family does not.
  for (link in c("logit", "probit")) {
    fit <- bartisan(y ~ x1 + x2 + vc(z, center = "estimate"), data = d,
                    family = binomial(link = link), control = vc_control())
    expect_true("b.z.0" %in% colnames(fit[["aux"]]))
  }

  # And a family whose target is not quadratic either way is still refused,
  # because the conjugate step would be a Laplace approximation with no
  # Metropolis correction.
  expect_error(
    bartisan(y ~ x1 + x2 + vc(z, center = "estimate"), data = d,
             family = binomial(link = "cloglog"), control = vc_control()),
    "leaf target is\\s+quadratic")
})

test_that("bcf falls back to a fixed coding where a drawn one is not exact", {
  d <- sim_vc(n = 250, seed = 30)
  d$count <- stats::rpois(250, exp(0.3 * d$x1 + 0.2 * d$z))

  # The coding is bcf's choice rather than the caller's, so a family that cannot
  # have it drawn gets the fixed default instead of an error.
  fit <- bcf(count ~ x1 + x2, treatment = ~ z, data = d, family = poisson(),
             propensity = FALSE, num_trees = 15L, control = vc_control())

  expect_false(any(startsWith(colnames(fit[["aux"]]) %or% character(), "b.")))
  expect_identical(colnames(coef(fit)), "z")

  # And the tree count follows the coding that was used, so a two-level factor
  # falling back to the symmetric coding gets a forest per level rather than the
  # single forest a drawn coding would have had.
  d$zf <- factor(ifelse(d$z == 1L, "yes", "no"))

  fell_back <- bcf(count ~ x1 + x2, treatment = ~ zf, data = d,
                   family = poisson(), propensity = FALSE,
                   control = vc_control())
  drawn <- bcf(y ~ x1 + x2, treatment = ~ zf, data = d, family = gaussian(),
               propensity = FALSE, control = vc_control())

  expect_identical(fell_back[["num_forest"]], 3L)
  expect_identical(drawn[["num_forest"]], 2L)
})

test_that("bcf reports its own call rather than the one do.call made", {
  d <- sim_vc(n = 200, seed = 31)

  fit <- bcf(y ~ x1 + x2, treatment = ~ z, data = d, family = gaussian(),
             propensity = FALSE, num_trees = 15L, control = vc_control())

  expect_match(deparse(fit[["call"]])[1L], "^bcf\\(")
})

# Varying coefficients across a family's additive predictors. The forest space is
# then two-dimensional -- parameter by coefficient -- and each parameter's own
# formula says which coefficients it has.

sim_two <- function(n = 900, seed = 40) {
  set.seed(seed)
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                  z = stats::rbinom(n, 1L, 0.5))
  d$tau_mean <- 1 + 0.8 * d$x2
  d$tau_sd <- 0.5 + 0.4 * d$x1
  d$y <- 2 * d$x1 + d$z * d$tau_mean +
    stats::rnorm(n, sd = exp(-0.7 + 0.3 * d$x1 + d$z * d$tau_sd))
  d
}

test_that("a vc() term reaches only the parameter whose formula names it", {
  d <- sim_two()

  on_mean <- bartisan(list(mean = y ~ x1 + x2 + vc(z), log_sd = ~ x1 + x2),
                      data = d, family = location_scale(),
                      control = vc_control())
  on_sd <- bartisan(list(mean = y ~ x1 + x2 + z, log_sd = ~ x1 + x2 + vc(z)),
                    data = d, family = location_scale(),
                    control = vc_control())

  expect_named(on_mean[["eta"]], c("mean", "mean:z", "log_sd"))
  expect_named(on_sd[["eta"]], c("mean", "log_sd", "log_sd:z"))

  expect_identical(colnames(coef(on_mean)), "mean:z")
  expect_identical(colnames(coef(on_sd)), "log_sd:z")
})

test_that("both parameters recover their own coefficient function", {
  d <- sim_two(seed = 41)

  fit <- bartisan(y ~ x1 + x2 + vc(z), data = d, family = location_scale(),
                  control = bartisan_control(num_trees = 25, num_burn = 400,
                                             num_draws = 400, verbose = FALSE))

  cf <- coef(fit)

  # Two coefficient functions of the same covariate, on different parameters and
  # varying with different predictors. Fitting them at once is the point.
  expect_gt(cor(cf[, "mean:z"], d$tau_mean), 0.9)
  expect_gt(cor(cf[, "log_sd:z"], d$tau_sd), 0.9)

  expect_lt(abs(mean(cf[, "mean:z"]) - mean(d$tau_mean)), 0.2)
  expect_lt(abs(mean(cf[, "log_sd:z"]) - mean(d$tau_sd)), 0.2)
})

test_that("each parameter's coefficient takes its own modifiers", {
  d <- sim_two(seed = 42)

  fit <- bartisan(list(mean = y ~ x1 + x2 + vc(z, ~ x2),
                       log_sd = ~ x1 + x2 + vc(z, ~ x1)),
                  data = d, family = location_scale(), control = vc_control())

  # A forest may only split on what its own modifiers allow, and the weights are
  # how that is enforced, so the mask is what to read.
  masks <- fit[["vc"]][["masks"]]

  expect_false(masks["x1", 2L])
  expect_true(masks["x2", 2L])
  expect_true(masks["x1", 4L])
  expect_false(masks["x2", 4L])
})

test_that("the same covariate may vary on two parameters but not twice on one", {
  d <- sim_two(seed = 43)

  # Once per parameter is the feature; twice in one formula is a mistake. The
  # union across the formulas has `z` wrapped twice either way, so the check has
  # to be per formula rather than on the union.
  expect_named(
    bartisan(list(mean = y ~ x1 + x2 + vc(z, ~ x2),
                  log_sd = ~ x1 + x2 + vc(z, ~ x1)),
             data = d, family = location_scale(),
             control = vc_control())[["eta"]],
    c("mean", "mean:z", "log_sd", "log_sd:z"))

  expect_error(bartisan(y ~ x1 + vc(z) + vc(z), data = d,
                        family = location_scale(), control = vc_control()),
               "a varying coefficient more than once")
})

test_that("a named formula list lines its vc() terms up with the parameter", {
  d <- sim_two(seed = 44)

  # Given out of order, so the `vc()` term has to follow the reordering rather
  # than the position it was written in.
  fit <- bartisan(list(log_sd = ~ x1 + x2, mean = y ~ x1 + x2 + vc(z)),
                  data = d, family = location_scale(), control = vc_control())

  expect_named(fit[["eta"]], c("mean", "mean:z", "log_sd"))
})

test_that("a drawn coding is judged per parameter and named per parameter", {
  d <- sim_two(seed = 45)

  # `location_scale()` is quadratic in the mean and not in the log standard
  # deviation, so the same request is exact on one and not on the other. The
  # guard has to be asked of the predictor the coding actually feeds.
  fit <- bartisan(list(mean = y ~ x1 + x2 + vc(z, center = "estimate"),
                       log_sd = ~ x1 + x2),
                  data = d, family = location_scale(), control = vc_control())

  expect_true(all(c("b.mean:z.0", "b.mean:z.1") %in% colnames(fit[["aux"]])))
  expect_identical(colnames(coef(fit)), "mean:z")

  expect_error(
    bartisan(list(mean = y ~ x1 + x2,
                  log_sd = ~ x1 + x2 + vc(z, center = "estimate")),
             data = d, family = location_scale(), control = vc_control()),
    "leaf target is\\s+quadratic")
})

test_that("group intercepts reach every control function and no coefficient", {
  d <- sim_two(seed = 46)
  d$grp <- factor(sample(1:10, nrow(d), TRUE))

  fit <- bartisan(y ~ x1 + x2 + vc(z) + (1 | grp), data = d,
                  family = location_scale(), control = vc_control())

  # One set per forest, and a coefficient's must be exactly zero: a group-varying
  # coefficient is a random slope, which `split_random()` refuses outright.
  expect_length(fit[["ranef"]], 4L)
  expect_true(any(fit[["ranef"]][[1L]] != 0))
  expect_true(all(fit[["ranef"]][[2L]] == 0))
  expect_true(any(fit[["ranef"]][[3L]] != 0))
  expect_true(all(fit[["ranef"]][[4L]] == 0))
})

test_that("per-forest settings are keyed by the two-part forest names", {
  d <- sim_two(seed = 47)

  fit <- bartisan(y ~ x1 + x2 + vc(z), data = d, family = location_scale(),
                  num_trees = c(mean = 20, `mean:z` = 8, log_sd = 10,
                                `log_sd:z` = 5),
                  control = vc_control())

  expect_identical(fit[["num_trees"]], c(20L, 8L, 10L, 5L))
})

test_that("a family whose forests are one parameter's levels refuses vc()", {
  d <- sim_two(seed = 48)
  d$g <- factor(sample(letters[1:3], nrow(d), TRUE))

  # Those forests are identified only up to a shared function, which reporting
  # removes; a coefficient forest per level would add one such direction per
  # coefficient and the reporting does not carry them.
  expect_error(bartisan(g ~ x1 + x2 + vc(z), data = d, family = multinomial(),
                        control = vc_control()),
               "not available for this family")
})

test_that("a pinned nuisance forest survives a varying coefficient", {
  d <- sim_effect(seed = 49)

  fam <- custom_family(
    logdens = function(y, eta, aux) {
      stats::dnorm(y, eta[, 1], exp(aux[1]), log = TRUE)
    },
    num_predictors = 1L, aux_names = "logsig", start = 0, aux_start = 0)

  fit <- bartisan(y ~ x1 + x2 + vc(z), data = d, family = fam,
                  control = vc_control())

  # The nuisance forest trails the additive predictors, carries no coefficient
  # and stays pinned -- without which the engine treats it as an ordinary forest
  # and the parameter is never drawn.
  expect_named(fit[["eta"]], c("(Intercept)", "z"))
  expect_gt(stats::sd(fit[["aux"]][, "logsig"]), 0)
  expect_lt(abs(exp(mean(fit[["aux"]][, "logsig"])) - 0.5), 0.15)
})

test_that("a two-part family with a discrete response takes coefficients too", {
  set.seed(50)
  n <- 700
  d <- data.frame(x1 = stats::rnorm(n), z = stats::rbinom(n, 1L, 0.5))
  d$count <- stats::rpois(n, exp(0.4 * d$x1 + 0.6 * d$z))

  fit <- bartisan(list(count = count ~ x1 + vc(z), zero = ~ x1), data = d,
                  family = zi_poisson(), control = vc_control())

  expect_named(fit[["eta"]], c("count", "count:z", "zero"))
  expect_identical(colnames(coef(fit)), "count:z")
})

test_that("the overlap warning says which formula it means", {
  d <- sim_two(seed = 51)

  expect_warning(bartisan(list(mean = y ~ x1 + x2 + z + vc(z),
                               log_sd = ~ x1 + x2),
                          data = d, family = location_scale(),
                          control = vc_control()),
                 "the mean formula")

  # With one additive predictor there is only one formula, so naming it would be
  # noise.
  expect_warning(bartisan(y ~ x1 + x2 + z + vc(z), data = d,
                          family = gaussian(), control = vc_control()),
                 "the model formula")
})

test_that("each parameter's estimand path and coefficient path agree", {
  d <- sim_two(seed = 52)

  fit <- bartisan(y ~ x1 + x2 + vc(z), data = d, family = location_scale(),
                  control = vc_control())

  # Every estimand reaches the model through `predict()` on modified `newdata`,
  # which rebuilds the basis and combines the forests. The mean's link is the
  # identity, so its contrast has to be the mean's coefficient and nothing else
  # -- and the log standard deviation's contrast has to be its own, on its own
  # scale. If either disagrees, one of the two paths is wrong.
  at <- function(value, type) {
    predict(fit, newdata = transform(d, z = value), type = type, draws = TRUE)
  }

  mean_gap <- colMeans(at(1, "response") - at(0, "response"))
  expect_equal(mean_gap, coef(fit)[, "mean:z"], tolerance = 1e-8,
               ignore_attr = TRUE)

  link <- at(1, "link")
  link0 <- at(0, "link")
  sd_gap <- colMeans(link[[2L]] - link0[[2L]])
  expect_equal(sd_gap, coef(fit)[, "log_sd:z"], tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that("bcf composes with a family that has two additive predictors", {
  d <- sim_two(seed = 53)

  # Newly reachable, and it exercises two rules at once: one formula reaches both
  # parameters, and the drawn coding a binary treatment would otherwise get is
  # refused on the log standard deviation, so the whole fit falls back to a fixed
  # coding rather than failing.
  fit <- bcf(y ~ x1 + x2, treatment = ~ z, data = d,
             family = location_scale(), propensity = FALSE, num_trees = 12L,
             control = vc_control())

  expect_named(fit[["eta"]], c("mean", "mean:z", "log_sd", "log_sd:z"))
  expect_identical(colnames(coef(fit)), c("mean:z", "log_sd:z"))
  expect_false(any(startsWith(colnames(fit[["aux"]]) %or% character(), "b.")))
})
