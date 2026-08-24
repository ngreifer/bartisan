# Package-level state, for the handful of things that are properties of the
# session rather than of a fit.
the <- new.env(parent = emptyenv())

is_null <- function(x) {
  isTRUE(length(x) == 0L)
}

`%or%` <- function(x, y) {
  if (is_null(x)) y else x
}

# Map a predictor to the unit interval. The cutpoint prior is uniform on a
# node's live range, so the transformation decides what "uniform" means: the
# quantile version makes the prior invariant to any monotone reparameterization
# of the predictor, while the range version keeps the original spacing. Both
# return a closure so that `predict()` can apply the identical map to new data.
make_unit_map <- function(x, type = "quantile") {
  ux <- unique(x[!is.na(x)])

  # Every map has to send a missing value to a missing value, so that the
  # splitting rules can decide what to do with it rather than being handed a
  # number that looks observed.
  if (length(ux) < 2L) {
    return(function(z) ifelse(is.na(z), NA_real_, 0.5))
  }

  if (length(ux) == 2L) {
    lo <- min(ux)
    hi <- max(ux)
    return(function(z) as.numeric(z > (lo + hi) / 2))
  }

  if (identical(type, "range")) {
    lo <- min(x, na.rm = TRUE)
    hi <- max(x, na.rm = TRUE)
    return(function(z) pmin(pmax((z - lo) / (hi - lo), 0), 1))
  }

  f <- stats::ecdf(x)

  function(z) f(z)
}

# Prior weight of each predictor group in the sparsity prior. One group per term
# in the formula, so the dummy columns of a factor are selected as a unit rather
# than competing with each other.
make_group_probs <- function(assign, term_labels) {
  groups <- unique(assign)
  n_col <- length(assign)

  i <- integer(0)
  j <- integer(0)
  x <- numeric(0)

  for (g in seq_along(groups)) {
    cols <- which(assign == groups[g])
    i <- c(i, cols)
    j <- c(j, rep(g, length(cols)))
    x <- c(x, rep(1 / length(cols), length(cols)))
  }

  probs <- Matrix::sparseMatrix(i = i, j = j, x = x,
                                dims = c(n_col, length(groups)))

  labels <- {
    if (is_null(term_labels)) as.character(groups)
    else term_labels[groups]
  }

  dimnames(probs) <- list(NULL, labels)

  probs
}

# Posterior summaries used by the print and summary methods.
post_summary <- function(x, level = 0.95) {
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  q <- stats::quantile(x, probs = probs, names = FALSE, na.rm = TRUE)

  c(mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    lower = q[1L],
    upper = q[2L])
}
