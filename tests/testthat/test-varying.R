# Reading `vc()` out of a formula. The sampler does not carry varying
# coefficients yet, so these test the specification layer: what the marker means,
# which forests may split on what, and the basis each coefficient multiplies.

sim_vc <- function(n = 120, seed = 1) {
  set.seed(seed)
  data.frame(y = stats::rnorm(n), x1 = stats::rnorm(n), x2 = stats::rnorm(n),
             z = stats::rbinom(n, 1L, 0.5),
             g = factor(sample(letters[1:3], n, TRUE)))
}

# The pipeline `bartisan()` runs, up to the point the engine takes over.
spec_of <- function(formula, data) {
  split <- split_random(formula)
  vc_split <- split_vc_terms(split[["fixed"]])
  mf <- stats::model.frame(vc_to_names(formula), data)
  mt <- stats::terms(vc_split[["fixed"]], data = mf)
  design <- build_design(mt, mf)

  resolve_vc(vc_split[["vc"]], mf, design,
             uses_dot(split[["fixed"]]))
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

test_that("a varying coefficient needs a single-predictor family", {
  d <- sim_effect(seed = 11)

  expect_error(bartisan(y ~ x1 + x2 + vc(z), data = d,
                        family = location_scale(), control = vc_control()),
               "one additive predictor")
})

test_that("coef refuses a model that has no coefficients", {
  d <- sim_effect(seed = 12)
  fit <- bartisan(y ~ x1 + x2, data = d, family = gaussian(),
                  control = vc_control())

  expect_error(coef(fit), "no varying coefficients")
})
