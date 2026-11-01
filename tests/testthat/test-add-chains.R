# Chains added to a fit after the fact. They run on the streams a fit with more
# chains would have used next, so the result has to be that fit exactly, not
# merely one like it.

local({
  d <- sim_x(n = 60L, p = 2L, seed = 121L)
  set.seed(122L)
  d$z <- stats::rbinom(nrow(d), 1L, stats::plogis(d$x1 - 0.5))
  d$g <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  d$w <- stats::runif(nrow(d), 0.5, 1.5)
  d$e <- stats::runif(nrow(d), 1, 3)
  d$y <- stats::rnorm(nrow(d), d$x1 + d$z)
  d$count <- stats::rpois(nrow(d), d$e * exp(d$x1))

  ctrl <- quick_control()

  patrick::with_parameters_test_that(
    "added chains reproduce a fit that ran them from the start:",
    {
      direct_call <- model
      direct_call[["chains"]] <- start + add

      set.seed(123L)
      direct <- eval(direct_call)

      start_call <- model
      start_call[["chains"]] <- start

      set.seed(123L)
      fit <- eval(start_call)

      # The seed in effect when chains are added plays no part.
      set.seed(999L)
      added <- add_chains(fit, add)

      expect_s3_class(added, class(direct), exact = TRUE)
      expect_identical(added[["chains"]], start + add)
      expect_equal(unclass(added), unclass(direct))

      # A `bcf()` fit's propensity score is a covariate of every chain, so the
      # new chains have to condition on the score the old ones did rather than
      # one estimated again.
      if (!is.null(fit[["bcf"]])) {
        score <- colnames(fit[["bcf"]][["propensity"]])
        expect_identical(added[["bcf"]][["propensity"]],
                         fit[["bcf"]][["propensity"]])
        expect_identical(added[["model"]][score], fit[["model"]][score])
      }
    },
    patrick::cases(
      `one chain plus three, default family` = list(
        model = quote(bartisan(y ~ x1 + x2, data = d, control = ctrl)),
        start = 1L, add = 3L),
      `two chains plus two` = list(
        model = quote(bartisan(y ~ x1 + x2, data = d, family = stats::gaussian(),
                               control = ctrl)),
        start = 2L, add = 2L),
      `random effects, weights and an offset` = list(
        model = quote(bartisan(count ~ x1 + x2 + (1 | g), data = d,
                               family = stats::poisson(), weights = w,
                               offset = log(e), control = ctrl)),
        start = 1L, add = 2L),
      `bcf() with a fitted propensity score` = list(
        model = quote(bcf(y ~ x1 + x2, treat = ~ z, data = d,
                          family = stats::gaussian(),
                          propensity_args = small_propensity(),
                          num_trees = 4L, num_burn = 10L, num_draws = 10L,
                          verbose = FALSE)),
        start = 1L, add = 1L)
    )
  )

  test_that("adding chains one at a time is the same as adding them together", {
    set.seed(124L)
    fit <- bartisan(y ~ x1 + x2, data = d, family = stats::gaussian(),
                    control = ctrl)

    once <- add_chains(fit, 2L)
    twice <- add_chains(add_chains(fit), 1L)

    expect_equal(unclass(twice), unclass(once))

    # And the result is diagnosed as the fit it now is.
    expect_identical(diagnose(once)[["chains"]], 3L)
  })

  test_that("add_chains() refuses what it cannot extend", {
    set.seed(125L)
    fit <- bartisan(y ~ x1 + x2, data = d, family = stats::gaussian(),
                    control = ctrl)

    expect_error(add_chains(list()), "bartisan_fit")
    expect_error(add_chains(fit, 0))
    expect_error(add_chains(fit, 1.5))

    old <- fit
    old[["next_stream"]] <- NULL
    expect_error(add_chains(old), "random number")

    # The new chains would be fitted to whatever `d` is now, so a changed
    # response is caught rather than pooled with the old chains.
    d$y <- d$y + 1
    expect_error(add_chains(fit), "not the data")
  })
})

test_that("consecutive fits draw different streams without a new seed", {
  d <- sim_x(n = 40L, p = 2L, seed = 126L)
  d$y <- stats::rnorm(nrow(d))

  set.seed(127L)

  # Restoring the session's stream after drawing the chains' streams from it,
  # without moving it on, handed the next fit the same streams.
  for (chains in 1:2) {
    a <- bartisan(y ~ ., data = d, family = stats::gaussian(),
                  control = quick_control(chains = chains))
    b <- bartisan(y ~ ., data = d, family = stats::gaussian(),
                  control = quick_control(chains = chains))

    expect_false(isTRUE(all.equal(a[["eta"]], b[["eta"]])))
  }
})
