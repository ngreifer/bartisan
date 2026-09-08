# Small, fast settings. These tests check that the plumbing is right, not that
# the fits are good; the statistical behavior is checked in test-recovery.R.
quick_control <- function(...) {
  args <- list(num_trees = 5L, num_burn = 30L, num_draws = 30L, verbose = FALSE)
  args[names(list(...))] <- list(...)
  do.call(bartisan_control, args)
}

sim_x <- function(n = 60, p = 3, seed = 1) {
  set.seed(seed)
  out <- as.data.frame(matrix(stats::runif(n * p), nrow = n))
  names(out) <- paste0("x", seq_len(p))
  out
}

# The invariant that ties the two halves of the package together: the additive
# predictor recorded while sampling must equal the one obtained by replaying the
# stored trees. A mismatch means the sampler's bookkeeping and the saved forests
# have drifted apart, which would silently corrupt every prediction.
#
# Two fits this does not apply to, both because `predict(type = "link")` is then
# a different quantity from `fit$eta` rather than the same one recomputed. A
# `vc()` fit combines its forests into one predictor on the way out while `eta`
# keeps them apart, so the two differ by the modifier and not by an error. And a
# fit with an offset needs the offset passed to `predict()`, which this does not
# do. Call it on the plain shapes.
expect_predictor_invariant <- function(fit, data, tolerance = 1e-6) {
  # Explicitly on the link scale: the predictor is what the trees encode, and
  # predict() defaults to the response scale.
  replayed <- stats::predict(fit, newdata = data, type = "link", draws = TRUE)

  if (!is.list(replayed)) {
    replayed <- list(replayed)
  }

  for (h in seq_along(replayed)) {
    testthat::expect_equal(as.vector(replayed[[h]]),
                           as.vector(fit[["eta"]][[h]]),
                           tolerance = tolerance)
  }

  invisible(fit)
}

# Sharing the tree topology across a family's forests is implemented but not
# offered: `bartisan_control()` has no argument for it and it is engaged by an
# internal flag, for the reasons recorded beside that flag in `R/control.R`.
# Tests reach it the only way anything can.
# The environment is taken by reference first: `bartisan:::the$x <- value` is
# not an assignment R will make, because the target of `<-` would have to be
# `:::`, and it fails with "there is no package called '*tmp*'".
with_shared_forests <- function(code) {
  flags <- bartisan:::the
  old <- flags[["share_forests"]]
  flags[["share_forests"]] <- TRUE

  on.exit({
    if (is.null(old)) {
      rm("share_forests", envir = flags)
    }
    else {
      flags[["share_forests"]] <- old
    }
  }, add = TRUE)

  force(code)
}
