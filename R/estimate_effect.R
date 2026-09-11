#' Causal effects from a fitted model
#'
#' Estimates the average or conditional effect of a treatment by g-computation
#' over the posterior draws. `estimate_effect()` predicts every unit under every
#' level of the treatment, contrasts those potential outcomes within each draw,
#' and averages over the units the estimand asks for.
#'
#' @param object a `<bartisan_fit>` object; the output of a call to [bcf()] or
#'   [bartisan()].
#' @param estimand `string`; which units the effect is averaged over. `"ATE"`
#'   (the default) uses all of them, and `"ATT"` and `"ATC"` those in the focal
#'   group, which by default is the treatment's second level for the former and
#'   its first for the latter. With a binary treatment that makes them the
#'   treated and the untreated without `focal` being named. `"CATE"` does not
#'   average at all and returns one effect per unit. Abbreviations and lowercase
#'   spellings are allowed.
#' @param treatment `string`; the name of the treatment variable. A fit from
#'   [bcf()] carries its own and needs none, so this is for a fit from
#'   [bartisan()], where nothing marks one predictor as the treatment.
#' @param comparison `string`; how the two potential outcomes are contrasted.
#'   Allowable options include `"difference"` (the default), `"ratio"`,
#'   `"lnratio"`, `"or"`, and `"lnor"`. The last two are available only when the
#'   response is a probability.
#' @param by optional; a one-sided formula or a variable name naming a grouping
#'   variable, in which case the effect is averaged within each of its levels
#'   rather than over the whole sample.
#' @param newdata optional; a data frame of units to average over. Default is
#'   the data the model was fit to, which is what makes the default estimand the
#'   sample average effect.
#' @param level `numeric`; the level of the credible interval. Default is `.95`.
#' @param interval `string`; `"eti"` (the default) for an equal-tailed interval
#'   from the quantiles of the draws, or `"hpdi"` for the highest posterior
#'   density interval, which is the shortest interval containing `level` of the
#'   posterior mass.
#' @param focal optional; the treatment level whose units the effect is averaged
#'   over for `"ATT"` and `"ATC"`. With a binary treatment it defaults to the
#'   second level for `"ATT"` and the first for `"ATC"`, which makes them the
#'   effect among the treated and among the untreated. With a treatment of more
#'   than two levels it has no default and must be supplied, and then `"ATT"` and
#'   `"ATC"` differ only in the level named.
#' @param type `string`; the prediction scale the effect is computed on, passed
#'   to [predict.bartisan_fit()]. Default is `"response"`, which is the scale on
#'   which an average of unit-level differences is the marginal effect. See
#'   Details before changing it.
#' @param plot `logical`; whether to draw the result rather than return it.
#'   Default is `FALSE`. `plot = TRUE` calls [plot.bartisan_effect()].
#' @param x a `<bartisan_effect>` object; the output of a call to
#'   `estimate_effect()`.
#' @param digits `integer`; the number of significant digits to print.
#' @param contrasts `string`; which contrasts to display when the treatment has
#'   more than two levels. `"focal"` (the default when a focal group is known)
#'   shows only the contrasts involving it, `"all"` shows every pairwise
#'   contrast, and a character vector of contrast labels shows those. All of
#'   them are computed either way; this only decides what is printed.
#' @param ... ignored.
#'
#' @returns
#' A `<bartisan_effect>` object, a data frame with one row per reported effect
#' and columns `contrast`, `estimate`, `lower`, and `upper`, plus `unit` when
#' `estimand = "CATE"` and the `by` variable's name when `by` is used. The
#' posterior draws of every reported quantity are kept in the `"draws"`
#' attribute, and the marginal mean of each potential outcome in
#' `"potential_outcomes"`, so a caller can re-contrast or re-summarize without
#' refitting.
#'
#' `print()` returns its input invisibly. `plot()` returns a \pkg{ggplot2}
#' object.
#'
#' @details
#' ## Why the Response Scale
#'
#' A varying coefficient is a contrast on the *link* scale: on a
#' `binomial("logit")` fit, `coef(object)` is a per-unit difference in log odds.
#' The average of those is the average conditional log odds ratio, which is not
#' the marginal log odds ratio and is not usually the quantity a treatment
#' question asks for. So the default here is `type = "response"`, where every
#' unit's contrast is on the scale the response is measured on and averaging
#' them is the marginal effect.
#'
#' `type = "link"` is still available and is the right choice for looking at how
#' the effect *varies*, since that is the scale the forest models it on. It is
#' the wrong choice for reporting an average.
#'
#' ## Marginal Rather Than Average Conditional
#'
#' For `"ratio"`, `"lnratio"`, `"or"` and `"lnor"` the potential outcomes are
#' averaged over units first and contrasted afterward, which gives the marginal
#' ratio. With `estimand = "CATE"` there is no averaging to do, so a ratio reported there is
#' a conditional ratio and the `print()` method says so.
#'
#' For the same reason `"or"` and `"lnor"` are not each other's `exp()` and
#' `log()`. Each summarizes the posterior of the quantity it names, and a
#' posterior mean does not survive a nonlinear transformation: `"or"` reports the
#' mean of the odds ratio and `"lnor"` the mean of its logarithm, which
#' exponentiates to something smaller. Report whichever scale the interval
#' should be symmetric on, which for a ratio is usually the log. The same holds
#' of `"ratio"` against `"lnratio"`.
#'
#' ## Which Level Is Treated
#'
#' With a binary treatment and no `focal`, which level is the treated one is
#' worked out from the treatment's values by the rules
#' \pkgfun{WeightIt}{weightit} uses, so that the same variable is read the same
#' way by both packages. In order: a `logical` treatment is treated at `TRUE`; a
#' `numeric` one is untreated at 0; a `character` or `factor` one whose levels
#' parse as numbers is untreated at 0, or else treated at the largest; and
#' otherwise the conventional names are matched, `"t"`, `"tr"`, `"treat"`,
#' `"treated"` and `"exposed"` against `"c"`, `"co"`, `"ctrl"`, `"control"` and
#' `"unexposed"`. When none of those applies the second level is taken as the
#' treated one and a message says so, that being the case where the guess can be
#' wrong. `"ATC"` takes the other level as its focal group, so both estimands
#' average over the group `focal` names and only the default differs.
#'
#' ## Multi-category Treatments
#'
#' There is nothing to work out from the values, so `focal` is required for
#' `"ATT"` and `"ATC"`, and those two then name the same estimand: the effect
#' among the units in the level named.
#'
#' Every pairwise contrast is computed. Which ones are shown is a display
#' choice, made by `contrasts` in the `print()` method: with a focal group the
#' default is to show only the contrasts involving it, matching how
#' \pkgfun{WeightIt}{weightit} and \pkgfun{cobalt}{bal.tab} use `focal`.
#'
#' A continuous treatment is not currently supported, because the effect of a continuous
#' treatment is a slope or a dose-response curve rather than a contrast of
#' levels. See the \CRANpkg{adrftools} package for tools to visualize and summarize the effect of a continuous treatment.
#'
#' @seealso [bcf()], which fits the model this is usually called on;
#'   [summary.bcf_fit()], which calls this with defaults and prints a compact
#'   block; [bartisan-marginaleffects] for the same estimands through
#'   \pkg{marginaleffects}, which also covers the ones not offered here
#'
#' @examplesIf rlang::is_installed("ggplot2")
#' data("rhc")
#' set.seed(123)
#'
#' fit <- bcf(death ~ age + sex + meanbp + aps, treatment = ~ rhc,
#'            data = rhc, num_trees = 10, num_burn = 50, num_draws = 50,
#'            verbose = FALSE)
#'
#' # The risk difference, averaged over everyone
#' estimate_effect(fit)
#'
#' # Among the treated, and as a risk ratio rather than a difference
#' estimate_effect(fit, estimand = "ATT", comparison = "ratio")
#'
#' # One effect per unit, ordered, with the average behind them
#' cate <- estimate_effect(fit, estimand = "CATE")
#' plot(cate)
#'
#' # By a subgroup, which is where a forest plot earns its keep
#' estimate_effect(fit, by = ~ sex)
#'
#' @export
estimate_effect <- function(object, treatment = NULL, estimand = "ATE",
                            comparison = "difference", by = NULL,
                            newdata = NULL, level = 0.95, interval = "eti",
                            focal = NULL, type = "response", plot = FALSE) {

  arg::arg_is(object, "bartisan_fit")
  arg::arg_number(level)
  arg::arg_between(level, c(0, 1), inclusive = FALSE)
  arg::arg_flag(plot)

  estimand <- arg::match_arg(toupper(estimand),
                             c("ATE", "ATT", "ATC", "CATE"))

  comparison <- arg::match_arg(comparison,
                               c("difference", "ratio", "lnratio", "or",
                                 "lnor"))

  interval <- arg::match_arg(tolower(interval), c("eti", "hpdi"))

  treatment <- effect_treatment(object, treatment)
  newdata <- effect_newdata(object, newdata, treatment)

  # Which levels there are to contrast is a property of the fitted model, not of
  # whatever units the caller asks to average over. Reading them from `newdata`
  # would make the conditional effect among the treated impossible to ask for,
  # since in that subset the treatment takes one value and a one-valued numeric
  # column is indistinguishable from a continuous one.
  fitted_z <- {
    if (is_null(object[["model"]]) || is_null(object[["model"]][[treatment]])) {
      newdata[[treatment]]
    }
    else {
      object[["model"]][[treatment]]
    }
  }

  z <- newdata[[treatment]]
  kind <- treatment_kind(fitted_z)

  if (identical(kind, "continuous")) {
    arg::err(c("{.arg treatment} {.val {treatment}} is continuous, and the
                effect of a continuous treatment is a slope rather than a
                contrast of levels",
               i = "Use {.fn marginaleffects::avg_slopes} for an average slope
                    or {.fn marginaleffects::plot_predictions} for a
                    dose-response curve; see {.topic `bartisan-marginaleffects`}.",
               i = "A treatment of more than two levels is supported when it is
                    a {.cls factor}."))
  }

  levs <- effect_levels(fitted_z)
  focal <- effect_focal(focal, levs, estimand, treatment, fitted_z)
  by <- effect_by(by, newdata)

  # The potential outcomes: one draws-by-units matrix per treatment level, each
  # from a frame in which every unit is assigned that level. The propensity
  # score is a function of the covariates alone, so it is carried through
  # unchanged rather than recomputed, which is what makes this the right
  # intervention and not a different model.
  po <- lapply(levs, function(a) {
    d <- newdata
    d[[treatment]] <- effect_assign(z, a, levels(fitted_z))
    draws <- stats::predict(object, newdata = d, type = type, draws = TRUE)

    if (!is.matrix(draws)) {
      arg::err(c("{.code type = \"{type}\"} does not give one number per
                  observation for this family, so there is no contrast to take.",
                 i = "A multinomial or ordinal response needs a scale that is
                      one-dimensional; try {.code type = \"mean\"} or {.code type = \"stdlv\"}."))
    }

    draws
  })

  names(po) <- as.character(levs)

  if (comparison %in% c("or", "lnor")) {
    rng <- range(unlist(lapply(po, range)), na.rm = TRUE)

    if (rng[1L] < 0 || rng[2L] > 1) {
      arg::err(c("{.code comparison = \"{comparison}\"} needs the response to be
                  a probability, and these predictions run from
                  {.val {signif(rng[1L], 3)}} to {.val {signif(rng[2L], 3)}}.",
                 i = "Use {.val difference}, {.val ratio} or {.val lnratio}."))
    }
  }

  # Which units the average runs over. `focal` names the group in both cases,
  # and the estimands differ only in which level it defaults to: the second for
  # `ATT` and the first for `ATC`, so that a binary treatment gives the treated
  # and the untreated without anything being named.
  keep <- switch(estimand,
                 ATE = , CATE = rep(TRUE, nrow(newdata)),
                 ATT = , ATC = as.character(z) == as.character(focal))

  if (!any(keep)) {
    arg::err("no observation is in the group {.arg estimand} = {.val {estimand}}
              asks to average over")
  }

  pairs <- effect_pairs(levs)

  out <- {
    if (identical(estimand, "CATE")) {
      effect_cate(po, pairs, comparison, level, interval, keep, newdata)
    }
    else {
      effect_marginal(po, pairs, comparison, level, interval, keep, newdata, by)
    }
  }

  attr(out, "contrast_map") <- data.frame(
    contrast = vapply(pairs, contrast_label, character(1L),
                      comparison = comparison),
    hi = vapply(pairs, function(p) as.character(p[["hi"]]), character(1L)),
    lo = vapply(pairs, function(p) as.character(p[["lo"]]), character(1L)),
    stringsAsFactors = FALSE)

  attr(out, "potential_outcomes") <- effect_po_summary(po, keep, level,
                                                       interval)
  attr(out, "estimand") <- estimand
  attr(out, "comparison") <- comparison
  attr(out, "interval") <- interval
  attr(out, "level") <- level
  attr(out, "treatment") <- treatment
  attr(out, "focal") <- focal
  attr(out, "type") <- type
  attr(out, "n_units") <- sum(keep)
  attr(out, "by") <- by[["name"]]

  class(out) <- c("bartisan_effect", "data.frame")

  # The marginal effect is kept beside the conditional ones so that the forest
  # plot can draw the band without a second call, and so that a reader of the
  # object can see what the units average to.
  if (identical(estimand, "CATE")) {
    attr(out, "marginal") <- effect_marginal(po, pairs, comparison, level,
                                             interval, keep, newdata, NULL)
  }

  if (!plot) {
    return(out)
  }

  plot(out)
}

# The treatment's name: a `bcf()` fit knows it, a `bartisan()` fit cannot.
effect_treatment <- function(object, treatment) {
  if (!is_null(treatment)) {
    arg::arg_string(treatment)
    return(treatment)
  }

  spec <- object[["bcf"]]

  if (!is_null(spec) && !is_null(spec[["treatment"]])) {
    return(spec[["treatment"]])
  }

  arg::err(c("{.arg treatment} must name the treatment variable.",
             i = "A fit from {.fn bcf} carries the name and needs none, but
                  nothing in a fit from {.fn bartisan} marks one predictor as
                  the treatment.",
             i = "For example {.code treatment = \"z\"}."))
}

# The units to average over, and the check that the treatment is among them.
effect_newdata <- function(object, newdata, treatment) {
  if (is_null(newdata)) {
    newdata <- object[["model"]]

    if (is_null(newdata)) {
      arg::err("this fit kept no data to average over, so {.arg newdata} is
                required")
    }
  }

  if (!is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }

  if (!treatment %in% names(newdata)) {
    arg::err("{.arg newdata} has no column {.val {treatment}}, which is the
              treatment this effect is a contrast on")
  }

  newdata
}

effect_levels <- function(z) {
  if (is.factor(z)) {
    return(levels(droplevels(z)))
  }

  sort(unique(stats::na.omit(z)))
}

# Assignment has to keep the column's type, since a factor predictor's levels
# are part of the model frame and a numeric one's are not.
effect_assign <- function(z, a, levs = NULL) {
  if (is.factor(z)) {
    # The model frame's levels rather than the subset's, so a `newdata` holding
    # one arm still produces a column the fit's design matrix recognizes.
    return(factor(rep(as.character(a), length(z)),
                  levels = levs %or% levels(z)))
  }

  if (is.character(z)) {
    return(rep(as.character(a), length(z)))
  }

  rep(as.numeric(a), length(z))
}

effect_focal <- function(focal, levs, estimand, treatment, z) {
  if (!is_null(focal)) {
    if (length(focal) != 1L || !as.character(focal) %in% as.character(levs)) {
      arg::err("{.arg focal} must be one of {.var {treatment}}'s levels,
                {.val {levs}}")
    }

    return(levs[[match(as.character(focal), as.character(levs))]])
  }

  # With more than two levels there is nothing to guess: which group is the
  # focal one is a question about the estimand rather than about the data, and
  # `"ATT"` and `"ATC"` then name the same quantity, following how
  # \pkgfun{WeightIt}{weightit} and \pkgfun{cobalt}{bal.tab} use `focal`.
  if (length(levs) > 2L) {
    if (!estimand %in% c("ATT", "ATC")) {
      return(NULL)
    }

    arg::err(c("{.arg focal} must name the focal treatment level for
                {.code estimand = \"{estimand}\"} when {.var {treatment}} has
                {length(levs)} levels.",
               i = "Its levels are {.val {levs}}.",
               i = "With more than two levels {.val ATT} and {.val ATC} are the
                    same estimand, the effect among the units in the level
                    named."))
  }

  guess <- treated_level(levs, z)

  # `"ATC"` is the effect among the untreated, so its focal group is the other
  # level. Both estimands average over the focal group; only the default
  # differs.
  chosen <- if (identical(estimand, "ATC")) {
    levs[[match(as.character(guess), as.character(levs)) %% 2L + 1L]]
  } else {
    guess
  }

  # Said out loud only when the heuristic had to fall back on the level order,
  # since that is the case where it can be wrong.
  if (is.null(attr(guess, "deduced")) && estimand %in% c("ATT", "ATC")) {
    what <- if (identical(estimand, "ATC")) "control" else "treated"
    lab <- as.character(chosen)
    arg::msg("assuming {.val {lab}} is the {what} level of {.var {treatment}};
              supply {.arg focal} if not")
  }

  chosen
}

# Which of two levels is the treated one, by the rules
# \pkgfun{WeightIt}{weightit} uses, so that the same treatment read by both
# packages is read the same way. In order: a logical treatment is treated at
# `TRUE`; a numeric one is control at 0; a character one that parses as numbers
# is control at 0; otherwise the conventional names are matched. The `"deduced"`
# attribute records whether any of those applied, because the last resort is the
# level order and that is the one worth warning about.
treated_level <- function(levs, z) {
  deduced <- function(x) structure(x, deduced = TRUE)

  if (is.logical(z)) {
    return(deduced(levs[[match("TRUE", as.character(levs))]]))
  }

  chars <- as.character(levs)

  if (is.numeric(levs)) {
    if (any(levs == 0)) {
      return(deduced(levs[[which(levs != 0)[1L]]]))
    }
  }
  else if (!anyNA(suppressWarnings(as.numeric(chars)))) {
    nums <- as.numeric(chars)

    if (any(nums == 0)) {
      return(deduced(levs[[which(nums != 0)[1L]]]))
    }

    return(deduced(levs[[which.max(nums)]]))
  }
  else {
    treated_names <- c("t", "tr", "treat", "treated", "exposed")
    control_names <- c("c", "co", "ctrl", "control", "unexposed")

    hit <- which(tolower(chars) %in% treated_names)

    if (length(hit) == 1L) {
      return(deduced(levs[[hit]]))
    }

    hit <- which(tolower(chars) %in% control_names)

    if (length(hit) == 1L) {
      return(deduced(levs[[which(seq_along(levs) != hit)[1L]]]))
    }
  }

  # Nothing in the values says which is which, so the second level it is, and
  # the caller is told.
  levs[[2L]]
}

effect_by <- function(by, newdata) {
  arg::when_not_null(
    by,
    arg::arg_or(
      arg::arg_formula(one_sided = TRUE),
      arg::arg_string()
    ))

  if (is_null(by)) {
    return(NULL)
  }

  if (rlang::is_formula(by)) {
    by <- all.vars(by)

    if (length(by) != 1L) {
      arg::err("{.arg by} must name exactly one variable, as in
                  {.code by = ~ sex}")
    }
  }

  if (!by %in% names(newdata)) {
    arg::err("{.arg by} names {.var {by}}, which is not a column of the data
              the effect is averaged over")
  }

  list(name = by, value = newdata[[by]])
}

# Every ordered pair, reference second, so a label reads "treated - control".
effect_pairs <- function(levs) {
  p <- length(levs)

  if (p == 2L) {
    return(list(list(hi = levs[[2L]], lo = levs[[1L]])))
  }

  out <- vector("list", p * (p - 1) / 2)

  k <- 1L
  for (i in seq_len(p)) {
    for (j in seq_len(i - 1L)) {
      out[[k]] <- list(hi = levs[[i]], lo = levs[[j]])
      k <- k + 1L
    }
  }

  out
}

contrast_label <- function(pair, comparison) {
  sep <- if (comparison %in% c("difference")) " - " else " / "
  paste0(as.character(pair[["hi"]]), sep, as.character(pair[["lo"]]))
}

# The contrast itself, on two vectors of draws that are already whatever the
# estimand asked for: marginal means, or one unit's potential outcomes.
apply_comparison <- function(hi, lo, comparison) {
  switch(comparison,
         difference = hi - lo,
         ratio = hi / lo,
         lnratio = log(hi / lo),
         or = (hi / (1 - hi)) / (lo / (1 - lo)),
         lnor = log((hi / (1 - hi)) / (lo / (1 - lo))))
}

# The shortest interval holding `level` of the draws. Reported rather than the
# equal-tailed one only on request, since the two differ only for a skewed
# posterior and the equal-tailed one is what a quantile of the draws gives.
hpd_interval <- function(x, level) {
  x <- sort(x[is.finite(x)])
  n <- length(x)

  if (n == 0L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  k <- max(1L, floor(level * n))

  if (k >= n) {
    return(c(lower = x[1L], upper = x[n]))
  }

  starts <- seq_len(n - k)
  width <- x[starts + k] - x[starts]
  i <- which.min(width)

  c(lower = x[i], upper = x[i + k])
}

effect_summary <- function(draws, level, interval) {

  if (identical(interval, "hpdi")) {
    bounds <- hpd_interval(draws, level)
  }
  else {
    probs <- c((1 - level) / 2, 1 + (level - 1) / 2)
    q <- stats::quantile(draws, probs = probs, names = FALSE, na.rm = TRUE)
    bounds <- c(lower = q[1L], upper = q[2L])
  }

  c(estimate = mean(draws, na.rm = TRUE),
    lower = bounds[["lower"]],
    upper = bounds[["upper"]])
}

# The marginal estimand: average the potential outcomes over the units within
# each draw, and only then contrast them. The other order gives the average of
# the conditional contrasts, which for anything but a difference is a different
# quantity.
effect_marginal <- function(po, pairs, comparison, level, interval, keep,
                            newdata, by) {

  if (is_null(by)) {
    groups <- list(list(label = NULL, keep = keep))
  }
  else {
    g <- by[["value"]]
    levs <- if (is.factor(g)) levels(droplevels(g)) else sort(unique(g))

    groups <- lapply(levs, function(l) {
      list(label = as.character(l),
           keep = keep & !is.na(g) & as.character(g) == as.character(l))
    })
  }

  # Empty groups are dropped before anything is allocated, so the number of
  # rows is known: one per surviving group per contrast.
  groups <- groups[vapply(groups, function(g) any(g[["keep"]]), logical(1L))]

  rows <- vector("list", length(groups) * length(pairs))
  draws <- vector("list", length(rows))
  keys <- character(length(rows))
  k <- 0L

  for (g in groups) {
    means <- lapply(po, function(m) rowMeans(m[, g[["keep"]], drop = FALSE]))

    for (pair in pairs) {
      lab <- contrast_label(pair, comparison)
      d <- apply_comparison(means[[as.character(pair[["hi"]])]],
                            means[[as.character(pair[["lo"]])]],
                            comparison)

      k <- k + 1L
      keys[k] <- if (is_null(g[["label"]])) lab else paste(g[["label"]], lab)
      draws[[k]] <- d

      row <- data.frame(contrast = lab, stringsAsFactors = FALSE)

      if (!is_null(g[["label"]])) {
        row[[by[["name"]]]] <- g[["label"]]
        row <- row[c(by[["name"]], "contrast")]
      }

      s <- effect_summary(d, level, interval)
      rows[[k]] <- cbind(row, as.data.frame(as.list(s)),
                         n = sum(g[["keep"]]))
    }
  }

  names(draws) <- keys

  out <- do.call(rbind, rows) |>
    unrowname()

  attr(out, "draws") <- draws
  out
}

# One effect per unit, with no averaging, so the contrast is conditional.
effect_cate <- function(po, pairs, comparison, level, interval, keep, newdata) {
  idx <- which(keep)
  rows <- vector("list", length(pairs))
  draws <- vector("list", length(pairs))

  for (i in seq_along(pairs)) {
    pair <- pairs[[i]]
    lab <- contrast_label(pair, comparison)
    hi <- po[[as.character(pair[["hi"]])]][, idx, drop = FALSE]
    lo <- po[[as.character(pair[["lo"]])]][, idx, drop = FALSE]
    d <- apply_comparison(hi, lo, comparison)

    draws[[i]] <- d

    s <- t(apply(d, 2L, effect_summary, level = level, interval = interval))
    rows[[i]] <- data.frame(contrast = lab, unit = idx,
                            estimate = s[, "estimate"],
                            lower = s[, "lower"],
                            upper = s[, "upper"],
                            stringsAsFactors = FALSE)
  }

  names(draws) <- vapply(pairs, contrast_label, character(1L),
                         comparison = comparison)

  out <- do.call(rbind, rows) |>
    unrowname()

  attr(out, "draws") <- draws
  out
}

# The marginal mean of each potential outcome, which is worth reporting beside
# the difference: a few hundred dollars means one thing against a baseline of
# six thousand and another against six hundred.
effect_po_summary <- function(po, keep, level, interval) {
  rows <- lapply(names(po), function(nm) {
    s <- po[[nm]][, keep, drop = FALSE] |>
      rowMeans() |>
      effect_summary(level, interval)

    data.frame(level = nm,
               estimate = s[["estimate"]],
               lower = s[["lower"]],
               upper = s[["upper"]],
               stringsAsFactors = FALSE)
  })

  do.call(rbind, rows) |>
    unrowname()
}

#' @rdname estimate_effect
#' @export
print.bartisan_effect <- function(x, digits = 3L, contrasts = NULL, ...) {

  arg::arg_whole_number(digits)

  estimand <- attr(x, "estimand")
  comparison <- attr(x, "comparison")
  interval <- attr(x, "interval")
  level <- attr(x, "level")
  treatment <- attr(x, "treatment")
  focal <- attr(x, "focal")
  by <- attr(x, "by")

  cli_cat("{.strong {effect_title(estimand, comparison)}}")
  cli::cat_line()

  cli_cat("Treatment: {.val {treatment}}")

  units <- attr(x, "n_units")

  if (estimand %in% c("ATT", "ATC") && !is_null(focal)) {
    group <- as.character(focal)
    cli_cat("Averaged over the {units} unit{?s} in group {.val {group}}")
  }
  else if (!identical(estimand, "CATE")) {
    cli_cat("Averaged over {units} unit{?s}")
  }

  if (!is_null(by)) {
    cli_cat("Within levels of {.val {by}}")
  }

  show <- effect_display(x, contrasts, focal)

  cli::cat_line()
  print(effect_round(show, digits), row.names = FALSE)
  cli::cat_line()

  band <- switch(interval,
                 hpdi = "highest posterior density interval",
                 "equal-tailed credible interval")

  cli::cli_bullets(c(i = "{.field estimate} is the posterior mean;
                          {.field lower} and {.field upper} bound the
                          {100 * level}% {band}."))

  if (identical(estimand, "CATE") && !identical(comparison, "difference")) {
    cli::cli_bullets(c(i = "These are {.emph conditional} {comparison}s, and
                            their average is not the marginal {comparison}."))
  }

  if (!identical(attr(x, "type"), "response")) {
    scale <- attr(x, "type")
    cli::cli_bullets(c(i = "Computed on the {.val {scale}} scale, where an
                            average of unit-level contrasts need not be the
                            marginal effect."))
  }

  invisible(x)
}

effect_title <- function(estimand, comparison) {
  what <- switch(comparison,
                 difference = "difference",
                 ratio = "ratio",
                 lnratio = "log ratio",
                 or = "odds ratio",
                 lnor = "log odds ratio")

  label <- switch(estimand,
                  ATE = "Average treatment effect",
                  ATT = "Average treatment effect on the treated",
                  ATC = "Average treatment effect on the controls",
                  CATE = "Conditional average treatment effects")

  sprintf("%s (%s)", label, what)
}

# Which contrasts to print. Everything is computed; this is display only, which
# is why it lives here and not in `estimate_effect()`.
effect_display <- function(x, contrasts, focal) {
  show <- as.data.frame(x)
  map <- attr(x, "contrast_map")
  labels <- map[["contrast"]]

  if (length(labels) < 2L) {
    return(show)
  }

  arg::when_not_null(
    contrasts,
    arg::arg_character
  )

  if (is_null(contrasts)) {
    contrasts <- if (is_null(focal)) "all" else "focal"
  }

  if (identical(contrasts, "all")) {
    return(show)
  }

  if (identical(contrasts, "focal")) {
    if (is_null(focal)) {
      arg::err("{.code contrasts = \"focal\"} needs a focal group, which is
                  set by {.arg focal} in {.fn estimate_effect}")
    }

    # From the levels the contrast was built from rather than from its
    # printed label, since a level's name can contain the separator.
    key <- as.character(focal)
    wanted <- labels[map[["hi"]] == key | map[["lo"]] == key]
  }
  else {
    unknown <- setdiff(contrasts, labels)

    if (!is_null(unknown)) {
      arg::err(c("{.arg contrasts} names {length(unknown)} contrast{?s} this
                    object does not have: {.val {unknown}}.",
                 i = "Its contrasts are {.val {labels}}."))
    }

    wanted <- contrasts
  }

  show[show[["contrast"]] %in% wanted, , drop = FALSE]
}

effect_round <- function(show, digits) {
  for (nm in c("estimate", "lower", "upper")) {
    if (!is_null(show[[nm]])) {
      show[[nm]] <- signif(show[[nm]], digits)
    }
  }

  show
}

#' @rdname estimate_effect
#' @export
plot.bartisan_effect <- function(x, ...) {
  require_ggplot2("treatment effects")

  estimand <- attr(x, "estimand")
  comparison <- attr(x, "comparison")
  by <- attr(x, "by")
  ylab <- effect_axis_label(estimand, comparison)
  null_at <- if (comparison %in% c("difference", "lnratio", "lnor")) 0 else 1

  if (identical(estimand, "CATE")) {
    return(effect_forest_units(x, ylab, null_at))
  }

  if (!is_null(by)) {
    return(effect_forest_groups(x, by, ylab, null_at))
  }

  effect_density(x, ylab, null_at)
}

effect_axis_label <- function(estimand, comparison) {
  what <- switch(comparison,
                 difference = "Effect",
                 ratio = "Ratio",
                 lnratio = "Log ratio",
                 or = "Odds ratio",
                 lnor = "Log odds ratio")

  if (identical(estimand, "CATE")) sprintf("Conditional %s", tolower(what))
  else what
}

# Units ordered by their estimate, with the marginal effect as a band behind
# them, which is what makes the spread readable as heterogeneity rather than as
# a list of numbers.
effect_forest_units <- function(x, ylab, null_at) {
  d <- as.data.frame(x)
  marg <- attr(x, "marginal")

  d <- do.call(rbind, lapply(split(d, d[["contrast"]]), function(z) {
    z <- z[order(z[["estimate"]]), , drop = FALSE]
    z[["rank"]] <- seq_len(nrow(z))
    z
  }))

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$rank, y = .data$estimate))

  if (!is_null(marg)) {
    m <- as.data.frame(marg)
    m <- m[match(unique(d[["contrast"]]), m[["contrast"]]), , drop = FALSE]
    p <- p +
      ggplot2::geom_rect(data = m,
                         ggplot2::aes(xmin = -Inf, xmax = Inf,
                                      ymin = .data$lower, ymax = .data$upper),
                         inherit.aes = FALSE, fill = "steelblue", alpha = 0.15) +
      ggplot2::geom_hline(data = m,
                          ggplot2::aes(yintercept = .data$estimate),
                          color = "steelblue")
  }

  p <- p +
    ggplot2::geom_hline(yintercept = null_at, linetype = 2, color = "grey40") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                           width = 0, color = "grey60") +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::labs(x = NULL, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank())

  if (length(unique(d[["contrast"]])) > 1L) {
    p <- p + ggplot2::facet_wrap(~ .data$contrast)
  }

  p
}

effect_forest_groups <- function(x, by, ylab, null_at) {
  d <- as.data.frame(x)
  d[["group"]] <- factor(d[[by]], levels = rev(unique(d[[by]])))

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$estimate, y = .data$group)) +
    ggplot2::geom_vline(xintercept = null_at, linetype = 2, color = "grey40") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
                           orientation = "y", width = 0, color = "grey60") +
    ggplot2::geom_point() +
    ggplot2::labs(x = ylab, y = by) +
    ggplot2::theme_bw()

  if (length(unique(d[["contrast"]])) > 1L) {
    p <- p + ggplot2::facet_wrap(~ .data$contrast)
  }

  p
}

# A scalar estimand has one posterior, and the posterior is the answer, so it is
# drawn rather than reduced to three numbers.
effect_density <- function(x, ylab, null_at) {
  draws <- attr(x, "draws")

  d <- do.call(rbind, lapply(names(draws), function(nm) {
    data.frame(contrast = nm, value = draws[[nm]], stringsAsFactors = FALSE)
  }))

  bands <- as.data.frame(x)
  bands[["contrast"]] <- {
    if (nrow(bands) == length(draws)) names(draws) else bands[["contrast"]]
  }

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_density(fill = "steelblue", alpha = 0.3, color = "steelblue") +
    ggplot2::geom_vline(xintercept = null_at, linetype = 2, color = "grey40") +
    ggplot2::geom_vline(data = bands,
                        ggplot2::aes(xintercept = .data$estimate)) +
    ggplot2::geom_vline(data = bands, ggplot2::aes(xintercept = .data$lower),
                        linetype = 3) +
    ggplot2::geom_vline(data = bands, ggplot2::aes(xintercept = .data$upper),
                        linetype = 3) +
    ggplot2::labs(x = ylab, y = "Posterior density") +
    ggplot2::theme_bw()

  if (length(draws) > 1L) {
    p <- p + ggplot2::facet_wrap(~ .data$contrast, scales = "free")
  }

  p
}
