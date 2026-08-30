# Small, fast settings. These tests check that the plumbing is right, not that
# the fits are good; the statistical behavior is checked in test-recovery.R.
quick_control <- function(...) {
  args <- list(num_trees = 5L, num_burn = 30L, num_save = 30L, verbose = FALSE)
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
