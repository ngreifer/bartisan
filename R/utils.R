# Package-level state, for the handful of things that are properties of the
# session rather than of a fit.
the <- new.env(parent = emptyenv())

is_null <- function(x) {
  isTRUE(length(x) == 0L)
}

`%or%` <- function(x, y) {
  if (is_null(x)) y else x
}

# The package's own inline markup. cli has no underline class, and a heading set
# in `{.strong}` comes out bold, which in a lot of terminals and in the fonts an
# editor pane uses sits close enough to the body text to read as an accident
# rather than as a heading. Adding a class is what `cli_div()` is for; see
# `?cli::"inline-markup"`. A container closes when the frame that opened it
# exits, so each helper below opens its own and the class is available in every
# string they format.
cli_theme <- list(
  span.underline = list("text-decoration" = "underline")
)

cli_cat <- function(..., .envir = parent.frame()) {
  cli::cli_div(theme = cli_theme)

  cli::format_inline(..., .envir = .envir) |>
    cli::cat_line()
}

# The bullet-list counterpart of `cli_cat()`. `cli::cli_bullets()` writes to
# stderr, which for a print method means the bullets are missing from
# `capture.output()` and from a knitted document: a reader of the vignette saw
# `diagnose()`'s "What to do" heading with nothing under it. `cli_fmt()` renders
# the same call to a character vector instead, which then goes to stdout with
# the rest of the output.
cli_bullets_cat <- function(text, .envir = parent.frame()) {
  cli::cli_div(theme = cli_theme)

  cli::cli_fmt(cli::cli_bullets(text, .envir = .envir)) |>
    cli::cat_line()
}

unrowname <- function(x) {
  rownames(x) <- NULL
  x
}

# Map a predictor to the unit interval. The cutpoint prior is uniform on a
# node's live range, so the transformation decides what "uniform" means: the
# quantile version makes the prior invariant to any monotone reparameterization
# of the predictor, while the range version keeps the original spacing. Both
# return a closure so that `predict()` can apply the identical map to new data.
#
# Every branch hands its work to one of the four constructors below rather than
# writing its closure here, for the reason `reporter_for()` in `R/progress.R`
# gives at length: a closure carries the frame it was written in, and one written
# here would carry `x` and `ux` whether it used them or not. The maps are stored
# in the fit, one per predictor column, so what they hold is paid for as long as
# the fit exists.
#
# The bigger leak was `type`, and it is why this wants to be strict rather than
# tidy. The two-value branch never looks at `type`, so its promise was never
# forced, and an unforced promise holds the frame it came from: here that is
# `unit_transform()`'s frame, which has the whole design matrix in it twice. A
# binary predictor's map, whose entire content is two numbers, serialized to
# 148 KB on a 400-row fit, and every indicator column of every factor had one.
make_unit_map <- function(x, type = "quantile") {
  ux <- unique(x[!is.na(x)])

  # Every map has to send a missing value to a missing value, so that the
  # splitting rules can decide what to do with it rather than being handed a
  # number that looks observed.
  if (length(ux) < 2L) {
    return(constant_map)
  }

  if (length(ux) == 2L) {
    return(midpoint_map(min(ux), max(ux)))
  }

  # Over every value rather than over the distinct ones, which is the same range
  # and is how it was written.
  if (identical(type, "range")) {
    return(range_map(min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
  }

  if (identical(type, "smoothcdf")) {
    return(smooth_cdf_map(x))
  }

  ecdf_map(stats::ecdf(x))
}

# The distribution function, estimated smoothly rather than empirically.
#
# This is the default because the two obvious maps each fail, and for opposite
# reasons. `"range"` rescales linearly, so a handful of outliers leaves the bulk
# of the data inside a sliver of `[0, 1]`; a cutpoint drawn uniformly on a node's
# live range then almost never lands where the structure is, and the failure
# deepens as the sample grows, because the extremes grow with it. `"quantile"`
# uses the empirical distribution function, which is a step, so the fitted
# function is a step in the predictor: it cannot be differentiated, and a
# relationship that is straight in `x` becomes a jump wherever the data are
# gappy. Smoothing the same distribution function keeps the cutpoints where the
# data are and leaves the coordinate strictly increasing and differentiable.
#
# Convolving with a normal is what does the smoothing, at the bandwidth rate
# Azzalini (1981) gives for a distribution function, which is `n^(-1/3)` rather
# than the `n^(-1/5)` that is right for a density. Tabulating on a grid and
# interpolating is what keeps it affordable: evaluating the sum costs `O(n)` a
# point, which is seconds per predictor on a large fit, where the table is built
# once and read in constant time. Linear interpolation of an increasing function
# is increasing, so the map keeps the one property everything here rests on.
#
# The scale is the smaller of the standard deviation and the interquartile range
# over 1.349, so that a long tail sets the bandwidth from the bulk rather than
# from itself.
smooth_cdf_map <- function(x, n_grid = 1024L) {
  x <- sort(x[!is.na(x)])
  n <- length(x)

  scale <- min(stats::sd(x), stats::IQR(x) / 1.349)

  if (!is.finite(scale) || scale <= 0) {
    scale <- stats::sd(x)
  }

  # No spread to smooth over, so the linear map is the honest answer.
  if (!is.finite(scale) || scale <= 0) {
    return(range_map(min(x), max(x)))
  }

  h <- scale * n^(-1 / 3)
  grid <- seq(x[1L], x[n], length.out = n_grid)
  at <- vapply(grid, function(v) mean(stats::pnorm((v - x) / h)), 0)

  # Stretched so that the observed range covers `[0, 1]` exactly, as the other
  # maps do. Normalizing on the smoothed values instead would leave a few
  # percent of the coordinate below the smallest observation and above the
  # largest, where the kernel puts mass and no cutpoint can ever be useful.
  lo <- at[1L]
  hi <- at[n_grid]

  if (!(hi > lo)) {
    return(range_map(min(x), max(x)))
  }

  f <- stats::approxfun(grid, (at - lo) / (hi - lo), rule = 2L)

  function(z) {
    pmin(pmax(f(z), 0), 1)
  }
}

# The four maps. The constant one needs nothing and so is the map itself; the
# other three are constructors, and `force()` on what they capture is not a
# detail of them but the whole reason they exist, since an argument left as a
# promise reaches back to the frame the call was made from.
constant_map <- function(z) {
  ifelse(is.na(z), NA_real_, 0.5)
}

midpoint_map <- function(lo, hi) {
  force(lo)
  force(hi)

  function(z) {
    as.numeric(z > (lo + hi) / 2)
  }
}

range_map <- function(lo, hi) {
  force(lo)
  force(hi)

  function(z) {
    pmin(pmax((z - lo) / (hi - lo), 0), 1)
  }
}

ecdf_map <- function(f) {
  force(f)

  function(z) {
    f(z)
  }
}

# Which predictor groups are a set of mutually exclusive indicators, and each
# observation's level within them.
#
# A rule on one indicator column of a factor can only peel one level off the
# rest, so the tree prior reaches 2^K - K of the B_K partitions of K levels:
# 27 of 52 at K = 5, and 1,014 of 115,975 at K = 10. Deshpande (2024) shows what
# that costs, which is partial pooling: the bulk is never divided, so a typical
# tree puts one level alone and the rest together whether the data want that or
# not. A rule that takes a *subset* of the levels reaches all of them, and it
# needs to know which level each observation is in.
#
# The test is on the columns rather than on the terms, because it is the columns
# that have to be mutually exclusive indicators for a level code to exist. That
# admits a factor's main effect and an interaction of factors, whose cells are
# also a partition, and correctly excludes a factor crossed with a numeric
# predictor, whose columns are not indicators.
level_codes <- function(x, assign) {
  groups <- sort(unique(assign))
  n <- nrow(x)

  n_levels <- integer(length(groups))
  cat_col <- rep.int(-1L, length(groups))
  codes <- vector("list", length(groups))

  for (g in seq_along(groups)) {
    cols <- which(assign == groups[g])

    if (length(cols) < 2L) {
      next
    }

    block <- x[, cols, drop = FALSE]
    observed <- !apply(is.na(block), 1L, any)

    if (!any(observed)) {
      next
    }

    seen <- block[observed, , drop = FALSE]

    # Indicators, and exactly one of them set per observation.
    if (!all(seen == 0 | seen == 1) || !all(rowSums(seen) == 1)) {
      next
    }

    code <- rep.int(-1L, n)
    code[observed] <- max.col(seen, ties.method = "first") - 1L

    n_levels[g] <- length(cols)
    codes[[g]] <- code
  }

  keep <- which(n_levels > 0L)

  if (is_null(keep)) {
    return(list(codes = matrix(-1L, n, 0L), cat_col = cat_col,
                n_levels = n_levels))
  }

  cat_col[keep] <- seq_along(keep) - 1L

  list(codes = matrix(unlist(codes[keep], use.names = FALSE), nrow = n),
       cat_col = cat_col,
       n_levels = n_levels)
}

# The level codes of new data, under the structure the fit recorded.
#
# Which groups are categorical, and how many levels each has, is a property of
# the fitted trees rather than of the data being predicted for, so it is read
# from the fit. Re-deriving it would let a `newdata` that happens not to satisfy
# the indicator test -- one row, or a group missing throughout -- disagree with
# the rules the trees actually carry, which would be wrong rather than merely
# unsupported.
apply_level_codes <- function(x, assign, info) {
  n_levels <- info[["n_levels"]]
  cat_col <- info[["cat_col"]]
  groups <- sort(unique(assign))
  n_cat <- sum(n_levels > 0L)

  if (n_cat == 0L) {
    return(matrix(-1L, nrow = nrow(x), ncol = 0L))
  }

  out <- matrix(-1L, nrow = nrow(x), ncol = n_cat)

  for (g in seq_along(groups)) {
    if (n_levels[g] == 0L) {
      next
    }

    cols <- which(assign == groups[g])
    block <- x[, cols, drop = FALSE]
    observed <- !apply(is.na(block), 1L, any)

    if (any(observed)) {
      seen <- block[observed, , drop = FALSE]
      out[observed, cat_col[g] + 1L] <- max.col(seen, ties.method = "first") - 1L
    }
  }

  out
}

# Prior weight of each predictor group in the sparsity prior. One group per term
# in the formula, so the dummy columns of a factor are selected as a unit rather
# than competing with each other.
make_group_probs <- function(assign, term_labels) {
  groups <- unique(assign)
  n_col <- length(assign)

  i <- integer()
  j <- integer()
  x <- numeric()

  for (g in seq_along(groups)) {
    cols <- which(assign == groups[g])
    i <- c(i, cols)
    j <- c(j, rep.int(g, length(cols)))
    x <- c(x, rep.int(1 / length(cols), length(cols)))
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

# The caller's relative weights, spread over every predictor group and
# normalized to a probability vector. Unnamed groups take a weight of 1, so
# `c(x1 = 3, x3 = 0.5)` on three predictors gives 3/4.5, 1/4.5 and 0.5/4.5.
#
# An unmatched name is an error rather than a silent no-op. A weight is a claim
# about a specific predictor, and a typo in it would otherwise change nothing
# and say nothing, which is the failure mode worth spending an error on.
# One column of splitting weights per forest the engine builds.
#
# Two things set them, and they multiply. `split_prior` is what the caller asked
# for, per forest or once for all of them. `masks` is what each forest's own
# formula allows: a term a forest's formula leaves out gets a weight of zero,
# which is what holds the forest to its own predictors without giving it a design
# matrix of its own.
#
# Returns NULL when there is nothing to say -- no weights asked for and every
# forest using every predictor -- so that the ordinary case reaches the engine
# exactly as it did before and the sparsity prior is left alone.
# Two matrices, because the engine has to tell two things apart. `prior` is
# weights the caller fixed, and fixing them is a statement that nothing should
# draw over them. `mask` is which predictors a forest's own formula lets it split
# on, which restricts what may be drawn without saying anything about how. Sent
# down one channel, as they were, the second silently cancelled the sparsity
# prior: a `vc()` term whose moderators are a strict subset of the covariates
# produced a matrix the caller never asked for, and the engine read any matrix as
# fixed weights. That is how `bcf()` with a propensity score lost its sparsity
# prior on both forests at once.
resolve_split_matrix <- function(split_prior, groups, labels, masks,
                                 n_forest) {
  restricted <- !all(masks)

  if (is_null(split_prior) && !restricted) {
    return(list(prior = NULL, mask = NULL))
  }

  # `split_prior` is one set of weights for every forest, unlike the other
  # per-forest settings. Its own names are predictors, so a named vector cannot
  # also be read as keyed by forest: `c(x3 = 0)` is a weight on `x3` everywhere,
  # and there is no spelling that would mean it for one forest only. A list was
  # once accepted here and never reachable, `bartisan_control()` having checked
  # `split_prior` with `arg_numeric` the whole time. A forest is held to its own
  # predictors by `masks` below rather than by a weight of the caller's.
  per_forest <- rep(list(split_prior), ncol(masks))

  uniform <- rep.int(1 / length(groups), length(groups))
  allow_all <- rep.int(1, length(groups))

  # An intercept-only forest never splits, so neither column is ever read; an
  # all-allowed mask keeps the arithmetic on the engine's side away from a
  # support of no groups, exactly as for the trailing pinned forests below.
  mask <- vapply(seq_len(ncol(masks)), function(h) {
    if (any(masks[, h])) as.numeric(masks[, h]) else allow_all
  }, numeric(length(groups)))

  prior <- if (!is_null(split_prior)) {
    vapply(seq_len(ncol(masks)), function(h) {
      if (!any(masks[, h])) {
        return(uniform)
      }

      weights <- resolve_split_weights(per_forest[[h]], groups) %or% uniform

      # A term this forest's formula does not name is not a term it may split
      # on, so the caller's weight on it is not one the engine should see.
      weights[!masks[, h]] <- 0

      if (sum(weights) == 0) {
        arg::err(c("The {.val {labels[h]}} forest has no predictor left to split on.",
                   i = "Its formula names only predictors its {.arg split_prior}
                          gives a weight of zero."))
      }

      weights / sum(weights)
    }, numeric(length(groups)))
  }

  # The trailing pinned forests are one leaf that never splits, so their columns
  # are never read; these keep the matrices rectangular.
  pad <- function(m, fill) {
    if (is_null(m) || ncol(m) >= n_forest) {
      return(m)
    }

    cbind(m, matrix(fill, nrow = length(groups), ncol = n_forest - ncol(m)))
  }

  list(prior = pad(prior, 1 / length(groups)), mask = pad(mask, 1))
}

resolve_split_weights <- function(split_prior, groups) {
  if (is_null(split_prior)) {
    return(NULL)
  }

  unknown <- setdiff(names(split_prior), groups)

  if (!is_null(unknown)) {
    arg::err(c("{.arg split_prior} names {?a predictor/predictors} the model
                does not have: {.val {unknown}}",
               i = "The model's predictors are {.val {groups}}."))
  }

  weights <- rep.int(1, length(groups))
  names(weights) <- groups
  weights[names(split_prior)] <- as.numeric(split_prior)

  # A zero weight holds one predictor out of every tree, which is a use for the
  # argument. Every weight zero holds all of them out, which leaves nothing to
  # split on and no model to fit. Checked here rather than in
  # `bartisan_control()`, which sees only the names the caller gave and not the
  # unnamed predictors that default to one.
  if (sum(weights) == 0) {
    arg::err(c("{.arg split_prior} gives every predictor a weight of zero",
               i = "At least one predictor must be available to split on."))
  }

  weights / sum(weights)
}

# The settings the engine holds one copy of per forest, with the value a forest
# takes when a named argument does not mention it. Each is a hyperparameter of
# one forest's prior, so there is a real question a caller might have about it
# per predictor -- how hard to shrink the scale forest against the mean forest,
# whether the sparsity prior should apply to one and not the other. `num_trees`,
# `sigma_mu` and `split_prior` are per-forest too and are resolved separately,
# because each needs something the others do not.
PER_FOREST_DEFAULTS <- list(
  gamma = 0.95,
  beta = 2,
  alpha = 1,
  alpha_scale = 0,
  alpha_shape_1 = 0.5,
  alpha_shape_2 = 1,
  update_s = TRUE,
  update_alpha = TRUE,
  update_sigma_mu = TRUE,
  bandwidth = 0.1,
  update_bandwidth = TRUE,
  bandwidth_every = 1L
)

# The variables a formula or expression names, as `model.frame()` will name
# them. Adapted from `WeightIt:::get_varnames()`.
#
# This is `all.vars()` with one difference, and the difference is the reason it
# exists: an extraction with `$`, `[[` or `[` is one variable, not the pieces it
# is written from. `all.vars(quote(d$g))` gives `c("d", "g")`, and neither is a
# column of the model frame, which holds one column named `d$g`. So a caller who
# writes `y ~ x + d$g` gets a name that matches nothing, and whatever the name
# was being used for silently does not happen to that variable.
#
# `[[` on the call rather than `as.list()`, because a terms object is a call
# carrying attributes and `as.list()` on one does not subset the way it does on
# a plain call: `as.list(tt)[-1L]` returns three elements for a three-element
# terms object rather than dropping the operator, which walks the wrong tree.
# Indexing positions directly is immune to that, so this takes a terms object
# and a bare expression alike, as `all.vars()` does.
get_varnames <- function(expr) {
  recurse <- function(e) {
    if (is.symbol(e)) {
      return(as.character(e))
    }

    if (!is.call(e)) {
      return(NULL)
    }

    fn <- e[[1L]]

    if (identical(fn, quote(`$`)) || identical(fn, quote(`[[`)) ||
        identical(fn, quote(`[`))) {
      return(deparse1(e))
    }

    seq_along(e)[-1L] |>
      lapply(function(i) recurse(e[[i]])) |>
      unlist()
  }

  recurse(expr)
}

# One formula per forest, from either a formula or a list of them.
#
# The first is the model for the main parameter and carries the response. The
# rest are the models for the other additive predictors, in the order the family
# documents, and need no response because it is the same response: a one-sided
# formula is the natural way to write them, and a two-sided one is accepted so
# long as it names the same outcome.
split_formula_list <- function(formula, arg = "formula") {
  if (rlang::is_formula(formula)) {
    return(list(formula))
  }

  if (is_null(formula) || !is.list(formula)) {
    arg::err("{.arg {arg}} must be a formula, or a list of formulas with one
              per forest")
  }

  bad <- which(!vapply(formula, rlang::is_formula, logical(1L)))

  if (!is_null(bad)) {
    arg::err("{.arg {arg}} must contain formulas; element {bad[1L]} is
              {.cls {class(formula[[bad[1L]]])[1L]}}")
  }

  two_sided <- vapply(formula, rlang::is_formula, logical(1L), lhs = TRUE)

  if (!any(two_sided)) {
    arg::err(c("One element of {.arg {arg}} must be two-sided.",
               i = "The model for the main parameter carries the response; the
                  rest need no response."))
  }

  # An unnamed list is positional, so the response has to be on the first entry.
  # A named one says which forest each formula is for, so the response can be
  # wherever the caller put that forest.
  if (is_null(names(formula)) && !two_sided[[1L]]) {
    arg::err(c("The first element of {.arg {arg}} must be two-sided.",
               i = "It is the model for the main parameter and carries the
                  response. Name the list to give the formulas in another
                  order."))
  }

  lhs <- formula[[which(two_sided)[[1L]]]][[2L]]

  out <- lapply(seq_along(formula), function(i) {
    f <- formula[[i]]

    if (!two_sided[[i]]) {
      return(f)
    }

    # A response on a later formula is allowed and has to be the same one, since
    # there is only one response being modeled.
    if (!identical(deparse(f[[2L]]), deparse(lhs))) {
      arg::err(c("Element {i} of {.arg {arg}} has a different response from the
                  first: {.code {deparse(f[[2L]])}} against {.code {deparse(lhs)}}.",
                 i = "Only the first formula needs a response, and every forest
                    models the same one."))
    }

    f
  })

  setNames(out, names(formula))
}

# The term labels of one forest's formula, with `.` expanded against the data the
# way it is for a single-formula fit.
#
# `response` is the first formula's left-hand side, put back on a one-sided
# formula before the terms are taken. Without it `~ .` expands over every column
# including the outcome, which would make the outcome a predictor of itself.
forest_terms <- function(f, data, response = NULL) {
  if (!is_null(response) && length(f) == 2L) {
    f <- stats::as.formula(call("~", response, f[[2L]]), env = environment(f))
  }

  tt <- {
    if (is_null(data)) stats::terms(f)
    else stats::terms(f, data = data)
  }

  attr(tt, "term.labels")
}

# A term label reduced to its set of variables, so that an interaction written
# `x1:x2` in one forest's formula and `x2:x1` in another's is recognized as the
# same predictor rather than becoming two.
term_key <- function(labels) {
  vapply(labels, function(l) {
    paste(sort(strsplit(l, ":", fixed = TRUE)[[1L]]), collapse = ":")
  }, character(1L), USE.NAMES = FALSE)
}

# The formula the model frame is built from: the response, and every predictor
# any forest uses. Each forest is then held to its own subset by its splitting
# weights rather than by a design matrix of its own, so a predictor absent from
# one forest's formula is present in that forest's data and never split on.
#
# Built from term labels rather than by joining the right-hand sides, so that a
# predictor two forests share appears once, an interaction written either way
# lands in one term, and `.` expands per formula.
union_formula <- function(formulas, data = NULL) {
  first <- formulas[[1L]]

  if (length(formulas) == 1L) {
    return(first)
  }

  # Not `formulas[[1L]]`: a named list may carry the response on any element,
  # because its names say which forest each formula is for, and this runs before
  # the family is known and so before the list can be put in forest order.
  response <- response_of(formulas)

  labels <- unlist(lapply(formulas, forest_terms, data = data,
                          response = response),
                   use.names = FALSE)

  # Deduplicated by variable set rather than by label, so the first spelling of
  # an interaction wins and the second does not become a separate predictor.
  labels <- labels[!duplicated(term_key(labels))]

  if (is_null(labels)) {
    arg::err("{.arg formula} must include at least one predictor")
  }

  stats::reformulate(labels, response = response,
                     env = environment(first))
}

# The response every forest models, from whichever formula carries it.
response_of <- function(formulas) {
  for (f in formulas) {
    if (rlang::is_formula(f, lhs = TRUE)) {
      return(f[[2L]])
    }
  }

  NULL
}

# Which of the model's predictor groups each forest may split on: one logical
# column per forest, over the terms of the union formula. A forest whose formula
# names every term gets a column of TRUE, which is the single-formula case and
# carries no restriction.
forest_masks <- function(formulas, groups, data, response) {
  keys <- term_key(groups)

  # `matrix()` rather than relying on `vapply()` to build one: with a single
  # predictor group it returns a bare vector, and everything downstream indexes
  # by column.
  vapply(formulas, function(f) {
    keys %in% term_key(forest_terms(f, data, response))
  }, logical(length(groups))) |>
    matrix(nrow = length(groups),
           ncol = length(formulas))
}

# Spread a per-forest argument over the reported forests.
#
# A setting that means something different for each additive predictor may be
# given once, to apply everywhere; positionally, with one entry per forest in the
# order the family documents; or keyed by the forest names, in which case a
# forest nobody named keeps the argument's own default. Lists and atomic vectors
# are both accepted, because some of these settings are strings and some are
# numbers and a caller should not have to remember which.
#
# `joint` is for the families whose forests are the levels of one vector-valued
# parameter rather than separate components of the response distribution. There
# is nothing a caller could mean by giving those forests different settings, so
# more than one value is an error rather than a silent recycling.
resolve_per_forest <- function(value, labels, arg, default = NULL,
                               joint = FALSE) {
  if (is_null(value)) {
    return(NULL)
  }

  n <- length(labels)
  value <- if (is.list(value)) value else as.list(value)
  nms <- names(value)

  if (length(value) == 1L && is_null(nms)) {
    return(rep(value, n))
  }

  if (joint) {
    arg::err(c("{.arg {arg}} must be one value for this family, not {length(value)}.",
               i = "Its {n} forests are the levels of one parameter and act
                  together, so they take the same setting."))
  }

  if (is_null(nms) || !any(nzchar(nms))) {
    if (length(value) != n) {
      arg::err(c("{.arg {arg}} has {length(value)} values but the model has
                  {n} forest{?s}",
                 i = "Give one value, or {n}, or name them:
                    {.val {labels}}."))
    }

    return(value)
  }

  named <- nzchar(nms)
  unknown <- setdiff(nms[named], labels)

  if (!is_null(unknown)) {
    arg::err(c("{.arg {arg}} names {?a forest/forests} this family does not
                have: {.val {unknown}}",
               i = "Its forests are {.val {labels}}."))
  }

  if (anyDuplicated(nms[named])) {
    dup <- unique(nms[named][duplicated(nms[named])])
    arg::err("{.arg {arg}} names {.val {dup}} more than once")
  }

  # Named entries take the forest they name, and the rest fill what is left in
  # order, which is how R matches arguments to a function's formals. It is what
  # makes `list(y ~ x1 + x2, log_sd = ~ x2)` mean the obvious thing.
  out <- rep(list(default), n)
  names(out) <- labels
  out[nms[named]] <- value[named]

  free <- setdiff(labels, nms[named])

  if (sum(!named) > length(free)) {
    arg::err(c("{.arg {arg}} has {sum(!named)} unnamed entr{?y/ies} but only
                {length(free)} forest{?s} left to take {?it/them}",
               i = "Its forests are {.val {labels}};
                  {.val {nms[named]}} {?is/are} already named."))
  }

  if (any(!named)) {
    out[free[seq_len(sum(!named))]] <- value[!named]
  }

  # A forest nobody named and nobody filled positionally keeps the default,
  # rather than borrowing a value the caller chose for a different forest.
  unname(out)
}

# The same, reduced to an atomic vector, for the settings the engine reads as
# one. `NULL` entries are what a named argument leaves behind where the caller
# named no forest, and they take the default that is passed in.
per_forest_vector <- function(value, labels, arg, default, joint = FALSE) {
  out <- resolve_per_forest(value, labels, arg, default = default, joint = joint)

  if (is_null(out)) {
    return(NULL)
  }

  if (!all(lengths(out) == 1L)) {
    arg::err("{.arg {arg}} must be one value per forest")
  }

  unlist(out, use.names = FALSE)
}

# One named field from every element of a list, as a vector of a known type.
#
# Replaces a bare `vapply()` over the extraction operator, which had drifted
# into two spellings, one of them wrapping an anonymous function around what the
# operator already does. `type` is the prototype `vapply()` checks each field
# against, so a field that is missing or of the wrong type is an error here
# rather than a surprise downstream. Where the fields are not all one type, or
# are not scalars, `lapply()` over the operator is still the way.
pluck <- function(x, field, type = character(1L)) {
  vapply(x, `[[`, type, field)
}

post_summary <- function(x, level = 0.95) {
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  q <- stats::quantile(x, probs = probs, names = FALSE, na.rm = TRUE)

  c(mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    lower = q[1L],
    upper = q[2L])
}

# Check if future and future.apply are installed and more than one
# worker is requested in the future plan.
use_future <- function() {
  rlang::is_installed(c("future", "future.apply")) &&
    future::nbrOfWorkers() > 1L
}
