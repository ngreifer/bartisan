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
#'   average at all and returns the conditional effect at each unit's
#'   covariates, which is not that unit's own individual effect; see Details.
#'   Abbreviations and lowercase spellings are allowed.
#' @param treat `string`; the name of the treatment variable. A fit from
#'   [bcf()] carries its own and needs none, so this is for a fit from
#'   [bartisan()], where nothing marks one predictor as the treatment.
#' @param comparison `string`; how the two potential outcomes are contrasted.
#'   Allowable options include `"difference"` (the default), `"ratio"`,
#'   `"lnratio"`, `"or"`, and `"lnor"`. The last two are available only when the
#'   response is a probability.
#' @param by optional; a one-sided formula or a variable name naming a grouping
#'   variable, in which case the effect is averaged within each of its levels
#'   rather than over the whole sample. A formula is evaluated with
#'   [stats::model.frame()] rather than read for the names it mentions, so it
#'   can define a grouping the data has no column for, as in `by = ~ age > 50`
#'   or `by = ~ interaction(sex, region)`. It must give exactly one grouping
#'   variable, and the term as written names the column it occupies in the
#'   output.
#' @param newdata optional; a data frame of units to average over. Default is
#'   the data the model was fit to, so that the average runs over the covariates
#'   of the fitted sample; see Details for what that average is.
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
#' @param x a `<bartisan_effect>` object; the output of a call to
#'   `estimate_effect()`.
#' @param digits `integer`; the number of significant digits to print.
#' @param contrasts `string`; which contrasts to display when the treatment has
#'   more than two levels. `"focal"` (the default when a focal group is known)
#'   shows only the contrasts involving it, `"all"` shows every pairwise
#'   contrast, and a character vector of contrast labels shows those. All of
#'   them are computed either way; this only decides what is printed.
#' @param potential_outcomes `logical`; whether to print the average response
#'   under each treatment level below the contrasts, those being what the
#'   contrasts were computed from. Default is `TRUE`. They are in the result's
#'   `"potential_outcomes"` attribute either way.
#' @param marginal `logical`; for `plot()` on conditional effects (i.e.,
#'   `estimand = "CATE"`), whether to draw the marginal effect beside the units.
#'   Default is `TRUE`. Ignored by the other plots, which draw no marginal
#'   effect beside their estimates.
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
#' `print()` returns its input invisibly, and prints the potential outcomes
#' below the contrasts unless `potential_outcomes = FALSE`. `plot()` returns a
#' \pkg{ggplot2} object.
#'
#' @details
#' ## The Estimand
#'
#' The ATE, ATT and ATC are averages of the conditional effect over the
#' covariates of the units averaged over, computed within each posterior draw.
#' This is sometimes called the mixed average treatment effect. The covariates
#' are treated as fixed, so the interval reflects uncertainty about the outcome
#' model but not the further variation a population average would carry, and
#' the units' observed outcomes are not conditioned on, as they would be for a
#' sample average of individual effects. Likewise, `estimand = "CATE"` reports
#' the expected effect at each unit's covariates rather than the unit's own
#' effect, which depends on how its two potential outcomes are associated and is
#' not identified. See `vignette("causal")` for the distinction between these
#' estimands.
#'
#' ## Setting `type`
#'
#' A varying coefficient is a contrast on the link scale: on a
#' `binomial("logit")` fit, `coef(object)` is a per-unit difference in log odds,
#' and the average of those is the average conditional log odds ratio rather than
#' the marginal one a treatment question usually asks for. The default is
#' therefore `type = "response"`, where every unit's contrast is on the scale the
#' response is measured on and averaging them gives the marginal effect.
#'
#' `type = "link"` is the right choice for looking at how the effect varies,
#' since that is the scale the forest models it on, and the wrong one for
#' reporting an average.
#'
#' ## Setting `comparison`
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
#' ## Setting `focal`
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
#' The effect of a continuous treatment is a slope or a dose-response curve
#' rather than a contrast of levels, and the \CRANpkg{adrftools} package has
#' tools for summarizing and visualizing one.
#'
#' @seealso [bcf()], which fits the model this is usually called on;
#'   [print.bcf_fit()] for its other methods; [`bartisan-marginaleffects`] for the
#'   same estimands through
#'   \pkg{marginaleffects}, which also covers the ones not offered here
#'
#' @examples
#' data("rhc")
#' set.seed(123)
#'
#' fit <- bcf(death ~ age + sex + meanbp + aps, treat = ~ rhc,
#'            data = rhc, num_trees = 10, num_burn = 50, num_draws = 50,
#'            verbose = FALSE)
#'
#' # The risk difference, averaged over everyone
#' estimate_effect(fit)
#'
#' # Among the treated, and as a risk ratio rather than a difference
#' estimate_effect(fit, estimand = "ATT", comparison = "ratio")
#'
#' # The effect at each unit's covariates, ordered, with the marginal effect
#' # beside them
#' cate <- estimate_effect(fit, estimand = "CATE")
#' plot(cate)
#'
#' # By a subgroup, which is where a forest plot earns its keep
#' estimate_effect(fit, by = ~ sex)
#'
#' @export
estimate_effect <- function(object, treat = NULL, estimand = "ATE",
                            comparison = "difference", by = NULL,
                            newdata = NULL, level = 0.95, interval = "eti",
                            focal = NULL, type = "response") {

  arg::arg_is(object, "bartisan_fit")
  arg::arg_number(level)
  arg::arg_between(level, c(0, 1), inclusive = FALSE)

  estimand <- arg::match_arg(toupper(estimand),
                             c("ATE", "ATT", "ATC", "CATE"))

  comparison <- arg::match_arg(comparison,
                               c("difference", "ratio", "lnratio", "or",
                                 "lnor"))

  interval <- arg::match_arg(tolower(interval), c("eti", "hpdi"))

  treat <- effect_treatment(object, treat)
  newdata <- effect_newdata(object, newdata, treat)

  # Which levels there are to contrast is a property of the fitted model, not of
  # whatever units the caller asks to average over. Reading them from `newdata`
  # would make the conditional effect among the treated impossible to ask for,
  # since in that subset the treatment takes one value and a one-valued numeric
  # column is indistinguishable from a continuous one.
  fitted_z <- {
    if (is_null(object[["model"]]) || is_null(object[["model"]][[treat]])) {
      newdata[[treat]]
    }
    else {
      object[["model"]][[treat]]
    }
  }

  z <- newdata[[treat]]
  kind <- treatment_kind(fitted_z)

  if (identical(kind, "continuous")) {
    arg::err(c("{.arg treat} {.val {treat}} is continuous, and the
                effect of a continuous treatment is a slope rather than a
                contrast of levels",
               i = "Use {.fn marginaleffects::avg_slopes} for an average slope
                    or {.fn marginaleffects::plot_predictions} for a
                    dose-response curve; see {.topic `bartisan-marginaleffects`}.",
               i = "A treatment of more than two levels is supported when it is
                    a {.cls factor}."))
  }

  levs <- effect_levels(fitted_z)
  focal <- effect_focal(focal, levs, estimand, treat, fitted_z)
  by <- effect_by(by, newdata)

  # The potential outcomes: one draws-by-units matrix per treatment level, each
  # from a frame in which every unit is assigned that level. The propensity
  # score is a function of the covariates alone, so it is carried through
  # unchanged rather than recomputed, which is what makes this the right
  # intervention and not a different model.
  predict_at <- function(a) {
    d <- newdata
    d[[treat]] <- effect_assign(z, a, levels(fitted_z))
    draws <- stats::predict(object, newdata = d, type = type, draws = TRUE)

    if (!is.matrix(draws)) {
      arg::err(c("{.code type = \"{type}\"} does not give one number per
                  observation for this family, so there is no contrast to take.",
                 i = "A multinomial or ordinal response needs a scale that is
                      one-dimensional; try {.code type = \"mean\"} or {.code type = \"stdlv\"}."))
    }

    draws
  }

  # One prediction per level, and none of them depends on another, so they go to
  # workers when a plan has any. The streams are drawn here rather than left to
  # `future.seed = TRUE` so that both branches below use the same ones: a family
  # whose prediction simulates rather than evaluating in closed form, as a
  # multinomial probit's orthant probabilities do, would otherwise give
  # different draws depending on whether a plan happened to be set.
  seeds <- parallel_streams(length(levs))

  if (use_future()) {
    po <- future.apply::future_lapply(levs, predict_at, future.seed = seeds,
                                      future.packages = "bartisan")
  }
  else {
    restore <- restore_stream()
    on.exit(restore(), add = TRUE)

    po <- lapply(seq_along(levs), function(i) {
      assign(".Random.seed", seeds[[i]], envir = globalenv())
      predict_at(levs[[i]])
    })
  }

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
    arg::err("no observation is in the group {.code estimand = \"{estimand}\"}
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

  # The draws behind that table, kept so that `diagnose()` can report on the two
  # averages as well as on their contrast. They are what the contrast is built
  # from, and a contrast that mixes badly can have one of them to blame rather
  # than both.
  attr(out, "po_draws") <- lapply(po, function(m) {
    rowMeans(m[, keep, drop = FALSE])
  })

  names(attr(out, "po_draws")) <- sprintf("Y[%s]", names(po))
  attr(out, "estimand") <- estimand
  attr(out, "comparison") <- comparison
  attr(out, "interval") <- interval
  attr(out, "level") <- level
  attr(out, "treatment") <- treat
  attr(out, "focal") <- focal
  attr(out, "type") <- type
  attr(out, "family") <- object[["family"]][["family"]]
  attr(out, "n_units") <- sum(keep)

  # What `diagnose()` needs to fold the draws back into chains and to name the
  # settings its advice would change. The estimand is what gets reported, and
  # its own mixing is not implied by the fit's: a contrast can stick where the
  # fitted function does not.
  attr(out, "chains") <- object[["chains"]] %or% 1L
  attr(out, "control") <- object[["control"]]
  attr(out, "by") <- by[["name"]]

  class(out) <- c("bartisan_effect", "data.frame")

  # The marginal effect is kept beside the conditional ones so that the forest
  # plot can draw it without a second call, and so that a reader of the object
  # can see what the units average to.
  if (identical(estimand, "CATE")) {
    attr(out, "marginal") <- effect_marginal(po, pairs, comparison, level,
                                             interval, keep, newdata, NULL)
  }

  out
}

# The treatment's name: a `bcf()` fit knows it, a `bartisan()` fit cannot.
effect_treatment <- function(object, treat) {
  if (!is_null(treat)) {
    arg::arg_string(treat)
    return(treat)
  }

  spec <- object[["bcf"]]

  if (!is_null(spec) && !is_null(spec[["treatment"]])) {
    return(spec[["treatment"]])
  }

  arg::err(c("{.arg treat} must name the treatment variable.",
             i = "A fit from {.fn bcf} carries the name and needs none, but
                  nothing in a fit from {.fn bartisan} marks one predictor as
                  the treatment.",
             i = "For example {.code treat = \"z\"}."))
}

# The units to average over, and the check that the treatment is among them.
effect_newdata <- function(object, newdata, treat) {
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

  if (!treat %in% names(newdata)) {
    arg::err("{.arg newdata} has no column {.val {treat}}, which is the
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

effect_focal <- function(focal, levs, estimand, treat, z) {
  if (!is_null(focal)) {
    if (length(focal) != 1L || !as.character(focal) %in% as.character(levs)) {
      arg::err("{.arg focal} must be one of {.var {treat}}'s levels,
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
                {.code estimand = \"{estimand}\"} when {.var {treat}} has
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
  chosen <- {
    if (identical(estimand, "ATC"))
      levs[[match(as.character(guess), as.character(levs)) %% 2L + 1L]]
    else
      guess
  }

  # Said out loud only when the heuristic had to fall back on the level order,
  # since that is the case where it can be wrong.
  if (is_null(attr(guess, "deduced")) && estimand %in% c("ATT", "ATC")) {
    what <- if (identical(estimand, "ATC")) "control" else "treated"
    lab <- as.character(chosen)
    arg::msg("assuming {.val {lab}} is the {what} level of {.var {treat}};
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
      arg::arg_string
    ))

  if (is_null(by)) {
    return(NULL)
  }

  if (!rlang::is_formula(by)) {
    if (!by %in% names(newdata)) {
      arg::err("{.arg by} names {.var {by}}, which is not a column of the data
                the effect is averaged over")
    }

    return(list(name = by, value = newdata[[by]]))
  }

  # The formula is evaluated rather than reduced to the names it mentions.
  # `all.vars()` read `by = ~ x3 > 0` as `x3` and then grouped by a continuous
  # predictor, one group per distinct value, which is the wrong answer and
  # gives no sign of being one.
  #
  # `model.frame()` does the evaluating, the same way the predictors themselves
  # are built from the fit's formula in `predict()`. It takes the terms one at
  # a time, it resolves a name the data does not have in the formula's own
  # environment rather than in this frame -- so `by = ~ x3 > cutoff` finds a
  # `cutoff` local to the function the call was written in, which evaluating
  # the right-hand side here would not -- and it names the column the way the
  # term reads, which is the name the result carries. `na.pass` because the
  # default action would drop the units whose group is missing and leave a
  # value that no longer lines up with the units it groups.
  label <- rlang::as_label(rlang::f_rhs(by))

  mf <- tryCatch(
    stats::model.frame(by, data = newdata, na.action = stats::na.pass),
    error = function(e) {
      arg::err("{.arg by} could not be evaluated in the data the effect is
                averaged over: {.code {label}}: {conditionMessage(e)}")
    })

  # A formula separates terms on `+` and `:`, so `~ sex + region` asks for two
  # groupings and there is one column here for each. Averaging within one
  # variable is what the result has a column for.
  if (ncol(mf) != 1L) {
    arg::err(c("{.arg by} must name one grouping variable, and {.code {label}}
                names {ncol(mf)}.",
               i = "It takes one term: cross two variables with
                    {.fn interaction}, as in
                    {.code by = ~ interaction(sex, region)}."))
  }

  value <- mf[[1L]]

  # A variable taken from the formula's environment rather than from the data
  # is not length-checked by `model.frame()`: a short one gives a frame whose
  # column disagrees with its own row names, and only printing it complains.
  if (length(value) != nrow(newdata)) {
    arg::err("{.arg by} must give one value per unit, and {.code {label}}
              gave {length(value)} for {nrow(newdata)} units")
  }

  list(name = names(mf)[1L], value = value)
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

# A label that names the quantity rather than leaving it to the heading. `1 / 0`
# is ambiguous between a ratio and a log odds ratio, and with numeric levels it
# reads as arithmetic on the numbers themselves, so the level goes inside a
# symbol: `Y[1]` cannot be misread as the number 1. `print()` prints one legend
# line saying what `Y[]` and `O()` are, whatever the number of contrasts.
contrast_label <- function(pair, comparison) {
  hi <- sprintf("Y[%s]", as.character(pair[["hi"]]))
  lo <- sprintf("Y[%s]", as.character(pair[["lo"]]))

  switch(comparison,
         difference = sprintf("%s - %s", hi, lo),
         ratio = sprintf("%s / %s", hi, lo),
         lnratio = sprintf("log(%s / %s)", hi, lo),
         or = sprintf("O(%s) / O(%s)", hi, lo),
         lnor = sprintf("log(O(%s) / O(%s))", hi, lo))
}

# The legend the labels above need, which depends on the comparison and on
# whether an average or a single unit is being reported.
contrast_legend <- function(comparison, treat, estimand) {
  what <- if (identical(estimand, "CATE")) {
    sprintf("{.field Y[a]} is the predicted response for that unit with
             {.var %s} set to {.val a}", treat)
  } else {
    sprintf("{.field Y[a]} is the average response with {.var %s} set to
             {.val a}", treat)
  }

  if (comparison %in% c("or", "lnor")) {
    paste0(what, ", and {.field O(y)} is the odds {.code y/(1-y)}.")
  } else {
    paste0(what, ".")
  }
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

  out <- do_rbind(rows) |>
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

  out <- do_rbind(rows) |>
    unrowname()

  attr(out, "draws") <- draws
  out
}

# The marginal mean of each potential outcome, which is worth reporting beside
# the difference: a few hundred dollars means one thing against a baseline of
# six thousand and another against six hundred.
effect_po_summary <- function(po, keep, level, interval) {
  lapply(names(po), function(nm) {
    s <- po[[nm]][, keep, drop = FALSE] |>
      rowMeans() |>
      effect_summary(level, interval)

    data.frame(quantity = sprintf("Y[%s]", nm),
               estimate = s[["estimate"]],
               lower = s[["lower"]],
               upper = s[["upper"]],
               stringsAsFactors = FALSE)
  }) |>
    do_rbind() |>
    unrowname()
}

#' @rdname estimate_effect
#' @export
print.bartisan_effect <- function(x, digits = 3L, contrasts = NULL,
                                  potential_outcomes = TRUE, ...) {

  arg::arg_whole_number(digits)
  arg::arg_flag(potential_outcomes)

  estimand <- attr(x, "estimand")
  comparison <- attr(x, "comparison")
  interval <- attr(x, "interval")
  level <- attr(x, "level")
  treat <- attr(x, "treatment")
  focal <- attr(x, "focal")
  by <- attr(x, "by")

  cli_cat("{.underline {effect_title(estimand, comparison)}}")
  cli::cat_line()

  cli_cat("Treatment: {.var {treat}}")

  units <- attr(x, "n_units")

  if (estimand %in% c("ATT", "ATC") && !is_null(focal)) {
    group <- as.character(focal)
    cli_cat("Averaged over the {units} unit{?s} in group {.val {group}}")
  }
  else if (!identical(estimand, "CATE")) {
    cli_cat("Averaged over {units} unit{?s}")
  }

  if (!is_null(by)) {
    cli_cat("Within levels of {.var {by}}")
  }

  show <- effect_display(x, contrasts, focal)

  cli::cat_line()

  if (identical(estimand, "CATE")) {
    print(effect_round(cate_spread(show), digits), row.names = FALSE)
  }
  else {
    print(effect_round(show, digits), row.names = FALSE)
  }

  # The two quantities the contrast is a contrast of, printed rather than left in
  # an attribute: a difference of a few points means one thing against a
  # baseline of .6 and another against .05, and a reader should not have to
  # reach for `attr()` to see which.
  po <- attr(x, "potential_outcomes")

  if (potential_outcomes && !is_null(po)) {
    cli::cat_line()
    cli_cat("{.underline Average potential outcomes}")
    cli::cat_line()
    print(effect_round(as.data.frame(po), digits), row.names = FALSE)
  }

  cli::cat_line()

  band <- switch(interval,
                 hpdi = "highest posterior density interval",
                 "equal-tailed credible interval")

  cli_bullets_cat(c(i = "{.field estimate} is the posterior mean;
                        {.field lower} and {.field upper} bound the
                        {100 * level}% {band}.",
                    i = contrast_legend(comparison, treat, estimand)))

  if (identical(estimand, "CATE")) {
    cli_bullets_cat(c(i = "Quartiles of the per-unit estimates. The object
                          itself holds one row per unit, with an interval
                          each."))
  }

  if (identical(estimand, "CATE") && !identical(comparison, "difference")) {
    cli_bullets_cat(c(i = "These are {.emph conditional} {comparison}s, and
                          their average is not the marginal {comparison}."))
  }

  if (!identical(attr(x, "type"), "response")) {
    scale <- attr(x, "type")
    cli_bullets_cat(c(i = "Computed on the {.val {scale}} scale, where an
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

# The quartiles of the per-unit estimates, which say how much the effect varies
# rather than how well any one unit is estimated.
cate_spread <- function(show) {
  lapply(split(show, show[["contrast"]]), function(z) {
    q <- stats::quantile(z[["estimate"]], c(0, 0.25, 0.5, 0.75, 1),
                         names = FALSE)
    data.frame(contrast = z[["contrast"]][1L], units = nrow(z), min = q[1L],
               q25 = q[2L], median = q[3L], q75 = q[4L], max = q[5L],
               stringsAsFactors = FALSE)
  }) |>
    do_rbind() |>
    unrowname()
}

effect_round <- function(show, digits) {
  for (nm in c("estimate", "lower", "upper", "min", "q25", "median", "q75",
               "max")) {
    if (!is_null(show[[nm]])) {
      show[[nm]] <- signif(show[[nm]], digits)
    }
  }

  show
}

#' @rdname estimate_effect
#' @export
plot.bartisan_effect <- function(x, marginal = TRUE, ...) {
  arg::arg_flag(marginal)

  estimand <- attr(x, "estimand")
  comparison <- attr(x, "comparison")
  by <- attr(x, "by")
  ylab <- effect_axis_label(estimand, comparison, attr(x, "family"),
                            attr(x, "type"))
  null_at <- if (comparison %in% c("difference", "lnratio", "lnor")) 0 else 1

  if (identical(estimand, "CATE")) {
    return(effect_forest_units(x, ylab, null_at, marginal))
  }

  if (!is_null(by)) {
    return(effect_forest_groups(x, by, ylab, null_at))
  }

  effect_density(x, ylab, null_at)
}

# The ratios name their scale in their own words; a difference does not, and
# "Effect" alone leaves a reader guessing whether 0.06 is a probability, a
# log odds, or the response's own units. So a difference says what it is a
# difference in where that is known: a probability for a binomial response,
# and the link scale whatever the family when `type = "link"`. Anything else is
# on the response's own scale, which the reader already knows.
effect_axis_label <- function(estimand, comparison, family = NULL,
                              type = NULL) {
  what <- switch(comparison,
                 difference = "Effect",
                 ratio = "Ratio",
                 lnratio = "Log ratio",
                 or = "Odds ratio",
                 lnor = "Log odds ratio")

  label <- {
    if (identical(estimand, "CATE")) sprintf("Conditional %s", tolower(what))
    else what
  }

  if (identical(comparison, "difference")) {
    scale <- {
      if (identical(type, "link")) "difference on the link scale"
      else if (identical(family, "binomial") &&
               (is_null(type) || identical(type, "response")))
        "difference in probability"
      else NULL
    }

    if (!is_null(scale)) {
      label <- sprintf("%s (%s)", label, scale)
    }
  }

  label
}

# Units ordered by their estimate, with the marginal effect beside them, which
# is what makes the spread readable as heterogeneity rather than as a list of
# numbers.
effect_forest_units <- function(x, ylab, null_at, marginal = TRUE) {
  d <- as.data.frame(x)
  marg <- if (marginal) attr(x, "marginal")

  d <- lapply(split(d, d[["contrast"]]), function(z) {
    z <- z[order(z[["estimate"]]), , drop = FALSE]
    z[["rank"]] <- seq_len(nrow(z))
    z
  }) |>
    do_rbind() |>
    unrowname()

  n <- max(d[["rank"]])

  # The units run from 0% to 100% of the ranking rather than 1 to n, so a
  # reader can say where in the distribution an effect sits ("the median unit",
  # "the top quarter") without counting. The first and last units sit on 0 and
  # 100 exactly.
  d[["pct"]] <- 100 * (d[["rank"]] - 1) / max(n - 1, 1)

  # The marginal effect goes beside the units rather than behind them. A band
  # across the panel is the obvious way to draw it and the wrong one: it sits
  # under every conditional interval, so the thing the reader most wants to
  # locate is the thing hardest to see. One interval past the right edge, in its
  # own color, is comparable by eye against any of them. The divider sits a
  # short way past the last unit, and the marginal effect is centered between
  # the divider and the edge of the panel, which the scale below does not pad.
  #
  # It is labeled as the marginal effect rather than as an average because it
  # is one only for a difference. It is the contrast of the potential outcomes
  # averaged over the units, and for a ratio or an odds ratio that is not the
  # average of the units' own contrasts, which is the point of reporting it.
  divider <- 104
  edge <- 116
  at <- (divider + edge) / 2

  # The marks shrink as the units multiply. At a few hundred they merge into a
  # solid block at any fixed size, which loses the very thing the plot is for.
  dot <- max(0.15, min(1.2, 60 / n))
  bar_alpha <- max(0.35, min(0.9, 200 / n))
  bar_width <- max(0.25, dot / 2)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$pct, y = .data$estimate)) +
    ggplot2::geom_hline(yintercept = null_at, linetype = 2, color = "grey40") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                           width = 0, color = "grey60", alpha = bar_alpha,
                           linewidth = bar_width) +
    ggplot2::geom_point(size = dot)

  if (!is_null(marg)) {
    m <- as.data.frame(marg)
    m <- m[m[["contrast"]] %in% unique(d[["contrast"]]), , drop = FALSE]
    m[["pct"]] <- at

    p <- p +
      ggplot2::geom_vline(xintercept = divider, color = "grey85") +
      ggplot2::geom_errorbar(data = m,
                             ggplot2::aes(ymin = .data$lower,
                                          ymax = .data$upper),
                             width = 0, linewidth = 1, color = "firebrick") +
      ggplot2::geom_point(data = m, size = 2, color = "firebrick") +
      ggplot2::scale_x_continuous(
        breaks = c(seq(0, 100, by = 25), at),
        labels = c(sprintf("%d%%", seq(0L, 100L, by = 25L)), "Marginal\neffect"),
        limits = c(0, edge),
        expand = ggplot2::expansion(mult = c(0.02, 0)))
  }
  else {
    p <- p +
      ggplot2::scale_x_continuous(
        breaks = seq(0, 100, by = 25),
        labels = sprintf("%d%%", seq(0L, 100L, by = 25L)),
        limits = c(0, 100),
        expand = ggplot2::expansion(mult = 0.02))
  }

  p <- p +
    ggplot2::labs(x = "Units, ordered by their conditional effect", y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
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

  d <- lapply(names(draws), function(nm) {
    data.frame(contrast = nm, value = draws[[nm]], stringsAsFactors = FALSE)
  }) |>
    do_rbind()

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
