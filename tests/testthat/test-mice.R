# `mice.impute.bartisan()`, the imputation method mice calls by name. mice hands
# it the incomplete variable, a numeric design matrix, and the rows to fit and
# to fill, and wants back one value of the variable's own type per row filled.

sim_incomplete <- function(n = 120, seed = 1) {
  d <- sim_x(n, 2, seed)
  x <- as.matrix(d)
  signal <- 2 * d$x1 - d$x2

  list(x = x,
       ry = seq_len(n) > 30,
       numeric = signal + stats::rnorm(n, sd = 0.5),
       binary = factor(ifelse(signal + stats::rnorm(n) > 0.5, "b", "a")),
       categorical = cut(signal + stats::rnorm(n), 3,
                         labels = c("lo", "mid", "hi")) |>
         factor(ordered = FALSE),
       ordered = cut(signal + stats::rnorm(n), 3,
                     labels = c("lo", "mid", "hi"), ordered_result = TRUE),
       logical = signal + stats::rnorm(n) > 0.5)
}

patrick::with_parameters_test_that(
  "`mice.impute.bartisan()` returns values of the variable's own type:",
  {
    s <- sim_incomplete()
    y <- s[[variable]]

    set.seed(1401)
    imputed <- mice.impute.bartisan(y, s$ry, s$x, num_trees = 5, num_burn = 30)

    expect_length(imputed, sum(!s$ry))
    expect_false(anyNA(imputed))
    expect_identical(class(imputed), class(y))

    if (is.factor(y)) {
      expect_identical(levels(imputed), levels(y))
    }
  },
  patrick::cases(
    numeric = list(variable = "numeric"),
    binary = list(variable = "binary"),
    categorical = list(variable = "categorical"),
    ordered = list(variable = "ordered"),
    logical = list(variable = "logical")
  )
)

test_that("`mice.impute.bartisan()` fills only the rows in `wy`", {
  s <- sim_incomplete()
  wy <- !s$ry
  wy[1:10] <- FALSE

  imputed <- mice.impute.bartisan(s$numeric, s$ry, s$x, wy = wy,
                                  family = gaussian(), num_trees = 5,
                                  num_burn = 30)
  expect_length(imputed, sum(wy))
})

test_that("`mice.impute.bartisan()` keeps the model's settings and drops mice's", {
  # One draw is all an imputation uses, unless the caller asks for more either
  # way it can.
  control <- imputation_control(list(num_trees = 7, donors = 5, eps = 1e-4,
                                     type = c(1, 1)))
  expect_s3_class(control, "bartisan_control")
  expect_identical(control[["num_draws"]], 1L)
  expect_equal(control[["num_trees"]], 7)

  control <- imputation_control(list(num_draws = 20))
  expect_equal(control[["num_draws"]], 20)

  control <- imputation_control(list(control = bartisan_control(num_draws = 15)))
  expect_equal(control[["num_draws"]], 15)

  control <- imputation_control(list(control = bartisan_control(num_trees = 9),
                                     num_burn = 40))
  expect_identical(control[["num_draws"]], 1L)
  expect_equal(control[["num_trees"]], 9)
  expect_equal(control[["num_burn"]], 40)

  expect_error(imputation_control(list(control = list(num_draws = 5))),
               "bartisan_control")
})

test_that("`mice.impute.bartisan()` refuses a variable with no predictors", {
  s <- sim_incomplete()

  expect_error(mice.impute.bartisan(s$numeric, s$ry, s$x[, 0L, drop = FALSE]),
               "predictorMatrix")
})

test_that("`mice.impute.bartisan()` draws from the predictive distribution", {
  skip_on_cran()

  # An imputation at the conditional mean would understate the spread; a draw
  # from the posterior predictive reproduces the residual standard deviation,
  # here 1, up to the uncertainty in the fit.
  set.seed(1402)
  n <- 600
  x <- matrix(stats::runif(n * 2), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
  mu <- 3 * x[, 1] - 2 * x[, 2]
  y <- mu + stats::rnorm(n)
  ry <- seq_len(n) > 200

  imputed <- mice.impute.bartisan(y, ry, x, family = gaussian(), num_trees = 20)

  expect_equal(stats::sd(imputed - mu[!ry]), 1, tolerance = 0.15)

  # Two calls use different random numbers, as consecutive calls inside mice()
  # must.
  again <- mice.impute.bartisan(y, ry, x, family = gaussian(), num_trees = 20)
  expect_false(isTRUE(all.equal(imputed, again)))
})

test_that("`mice()` finds and uses the method", {
  skip_if_not_installed("mice")

  s <- sim_incomplete(80)
  d <- data.frame(s$x, y = s$numeric, g = s$binary)
  is.na(d$y[1:15]) <- TRUE
  is.na(d$g[16:30]) <- TRUE

  set.seed(1403)
  imp <- mice::mice(d, method = "bartisan", m = 2L, maxit = 2L,
                    num_trees = 5, num_burn = 30, printFlag = FALSE)

  expect_identical(unname(imp[["method"]][c("y", "g")]),
                   c("bartisan", "bartisan"))

  completed <- mice::complete(imp, 2L)
  expect_false(anyNA(completed))
  expect_identical(levels(completed$g), levels(d$g))
})
