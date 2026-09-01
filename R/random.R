# Random-effect terms in the formula, in lme4's notation: `y ~ x + (1 | g)`.
#
# Only random *intercepts* are supported, and the reason is worth stating rather
# than leaving as an omission. A random intercept is a scalar entering the
# predictor with weight one for the observations in its level, which is exactly
# what a leaf is once its gate is removed -- so the sampler's leaf machinery
# handles it, including the closed forms for a quadratic or exponential target.
# A random *slope* is a scalar multiplying a covariate, which is a different
# shape of parameter and would need its own update; and a correlated pair of
# them would need a joint one. Neither is here, and asking for either says so.
#
# Several grouping factors are fine, and nesting is fine because `(1 | a/b)`
# expands to two intercept terms before this code sees it.

# The bar terms and the fixed part, with the bars validated.
split_random <- function(formula) {
  bars <- reformulas::findbars(formula)

  if (is_null(bars)) {
    return(list(fixed = formula, bars = list()))
  }

  for (b in bars) {
    lhs <- b[[2L]]

    # `1` parses as a numeric literal, and `(1 | g)` is the only left-hand side
    # that names an intercept alone. Anything else is a slope.
    if (!is.numeric(lhs) || length(lhs) != 1L || lhs != 1) {
      arg::err("only random intercepts are supported, so the left of every
                {.code |} must be {.code 1}; {.code ({deparse(lhs)} |
                {deparse(b[[3L]])})} asks for a random slope. Put the variable in
                the fixed part of the formula instead, where a tree can use it")
    }
  }

  list(fixed = reformulas::nobars(formula), bars = bars)
}

# Each grouping factor as a vector of zero-based level codes, plus its levels.
# The grouping expression is evaluated in the model frame, which is where
# `subbars()` put it, so `(1 | a:b)` works for the same reason `(1 | g)` does.
random_terms <- function(bars, mf) {
  if (is_null(bars)) {
    return(list())
  }

  out <- lapply(bars, function(b) {
    expr <- b[[3L]]
    label <- deparse(expr)
    value <- eval(expr, mf, environment())

    if (is_null(value)) {
      arg::err("grouping variable {.val {label}} was not found in the data")
    }

    f <- as.factor(value)

    if (anyNA(f)) {
      arg::err("grouping variable {.val {label}} has missing values, and a row
                with no group cannot be given a group effect. Drop those rows,
                or give them a level of their own")
    }

    levs <- levels(f)

    if (length(levs) < 2L) {
      arg::err("grouping variable {.val {label}} has {length(levs)} level{?s},
                so there is nothing for a group effect to vary over")
    }

    list(label = label, expr = expr, levels = levs,
         num_levels = length(levs),
         codes = as.integer(f) - 1L)
  })

  labels <- vapply(out, function(z) z[["label"]], character(1L))

  if (anyDuplicated(labels) > 0L) {
    arg::err("the same grouping variable appears twice in the formula:
              {.val {labels[anyDuplicated(labels)]}}")
  }

  setNames(out, labels)
}

# What the engine needs: the label, the codes and the level count.
random_spec <- function(terms) {
  lapply(terms, function(z) {
    list(label = z[["label"]], levels = z[["codes"]],
         num_levels = z[["num_levels"]])
  })
}

# The random part of the predictor for new data: one column per additive
# predictor, or NULL when the model has no random effects.
#
# A level that was not present at fitting time has no drawn intercept, so it is
# given zero -- the prior mean, which is the only defensible answer and is what
# `lme4::predict(allow.new.levels = TRUE)` does. It is reported rather than done
# silently, because a whole column of unseen levels usually means the grouping
# variable was coded differently.
random_predict <- function(object, newdata, iterations) {
  terms <- object[["random"]]

  if (is_null(terms)) {
    return(NULL)
  }

  draws <- object[["ranef"]]
  n_new <- nrow(newdata)
  unseen <- character()

  out <- lapply(seq_along(draws), function(h) {
    total <- matrix(0, nrow = length(iterations), ncol = n_new)
    at <- 0L

    for (z in terms) {
      value <- eval(z[["expr"]], newdata, environment())
      code <- match(as.character(value), z[["levels"]])
      block <- draws[[h]][iterations, at + seq_len(z[["num_levels"]]),
                          drop = FALSE]

      if (anyNA(code)) {
        unseen <<- c(unseen, z[["label"]])
        # A column of zeros for the unseen levels, which is the prior mean.
        block <- cbind(block, 0)
        code[is.na(code)] <- z[["num_levels"]] + 1L
      }

      total <- total + block[, code, drop = FALSE]
      at <- at + z[["num_levels"]]
    }

    total
  })

  if (!is_null(unseen)) {
    bad <- unique(unseen)
    arg::wrn(c("{.arg newdata} has levels of {.val {bad}} that were not
                present when the model was fit.",
               i = "Those rows get the group effect's prior mean of zero,
                    which is what a group with no data can be given."))
  }

  out
}
