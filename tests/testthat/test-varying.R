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
  expect_identical(one[["center"]], "mean")

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

  mean_centred <- spec_of(y ~ x1 + vc(z), d)
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

test_that("the sampler refuses varying coefficients rather than ignoring them", {
  d <- sim_vc()

  # Until the engine carries them, a fit would silently be the model without
  # them, which is the failure this feature exists to avoid.
  expect_error(bartisan(y ~ x1 + x2 + vc(z), data = d, family = gaussian(),
                        control = quick_control()),
               "not fitted yet")
})
