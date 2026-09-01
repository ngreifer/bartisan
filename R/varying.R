# Varying coefficients: the `vc()` marker and the machinery that reads it out of
# a formula.
#
# The model is
#
#   g(mu_i) = f_0(Z_i) + sum_j (X_ij - c_j) f_j(Z_i)
#
# with a forest for the control function f_0 and one for each coefficient f_j.
# Deshpande, Bai, Balocchi, Starling and Weiss (2026) call this VCBART; Hahn,
# Murray and Carvalho (2020) is the case of one binary X_j, and Woody, Carvalho,
# Hahn and Murray (2020) the case of one continuous one, whose vocabulary this
# file follows: f_0 is the *control function* over *control variables*, f_j is an
# *exposure moderating function* over *moderators*.

#' Give a predictor a varying coefficient
#'
#' Used only inside a [bartisan()] formula, where it says that the coefficient of
#' `x` is a function of other predictors rather than a constant. It is not meant
#' to be called directly and does nothing if it is.
#'
#' @param x the predictor whose coefficient varies. A numeric variable gets one
#'   forest; a factor gets one per level, coded symmetrically as
#'   [multinomial()] codes its predictors.
#' @param modifiers a one-sided formula naming the predictors this coefficient's
#'   forest may split on. The default, `NULL`, is every predictor in the model
#'   except `x` itself.
#' @param center the value of `x` at which the control function is read, which is
#'   both how the model is fitted and how it is reported. `"auto"`, the default,
#'   uses `"zero"` for a `0`/`1` covariate and `"mean"` for any other numeric
#'   one, for the reason in Details. `"mean"` centers `x`, so the control
#'   function is the surface at its average; `"zero"` leaves `x` alone, so the
#'   control function is the surface at `x = 0`; `"mid"` uses the midpoint of
#'   `x`'s range; a number uses that number. For a factor, `center` is `"mean"`
#'   or the name of a level to report against.
#'
#' @returns
#' Nothing. `vc()` is a marker read out of the formula and never evaluated.
#'
#' @details
#' The model is
#'
#' \deqn{g(\mu_i) = f_0(Z_i) + \sum_j (X_{ij} - c_j) f_j(Z_i)}
#'
#' with a forest for the control function \eqn{f_0} and one for each varying
#' coefficient \eqn{f_j}. Every forest is fitted at once, so the coefficient has
#' a prior of its own rather than being whatever difference a single forest with
#' `x` among its predictors happens to produce.
#'
#' # Which predictors a coefficient may vary with
#'
#' By default every predictor in the model except `x` itself. The `modifiers`
#' argument narrows that, and naming something that is not a predictor is an
#' error rather than a silent restriction.
#'
#' **A varying covariate may not modify the control function.** With `z` among
#' \eqn{f_0}'s predictors, \eqn{f_0(Z) + z f_1(Z)} is not identified: any
#' function of `z` moves between the two. Reached through `.`, the covariate is
#' dropped from the control function without comment, since `.` did not name it.
#' Named outright, the model is fitted as asked and a warning says why that is a
#' choice.
#'
#' **A numeric covariate may modify its own coefficient, and doing so is how the
#' effect stops being linear.** `vc(z, ~ z + x1)` fits \eqn{z f_1(z, x_1)}, so
#' the slope itself moves across `z`'s range and the dose response is a curve
#' rather than a line. This is worth reaching for whenever the effect of a
#' continuous predictor might not be proportional to it. On a simulation where
#' the truth is \eqn{y = 2x_1 + z^2}, letting \eqn{f_1} split on `z` recovers
#' the fitted surface to a root mean squared error of 0.099 where forbidding it
#' gives 1.187, and the fitted coefficient traces the truth: -1.15, 0.15 and 1.09
#' at `z` of -1, 0 and 1, where the slope of \eqn{z^2} is -1, 0 and 1.
#'
#' **A categorical covariate may not**, and is removed from its own forests
#' quietly. A level's indicator is nonzero only on the rows where that level
#' holds, and the variable is constant on exactly those rows, so such a split
#' separates rows that contribute from rows that contribute nothing. It is wasted
#' rather than unidentified.
#'
#' # Where the control function sits
#'
#' Centring is a reparameterization of \eqn{f_0} alone: every coefficient and
#' every estimand is identical under any choice, and what changes is what the
#' control function means. `"auto"` picks by the covariate, because neither
#' answer wins everywhere. For a `0`/`1` covariate it uses zero, so \eqn{f_0} is
#' the surface among the untreated -- a quantity with its own meaning, and the
#' one that recovers the coefficient best, at a correlation of 0.987 against
#' 0.975 for mean-centring on the simulation in `_dev/`. For any other numeric
#' covariate it uses the mean, because zero may be nowhere near the data: with a
#' covariate around 50 the control function at zero is an extrapolation and
#' recovery collapses to a correlation of 0.42.
#'
#' A factor is always fitted mean-centred and gets one forest per level, coded
#' symmetrically the way [multinomial()] codes its predictors rather than as
#' contrasts against a level that happened to sort first. That coding carries one
#' spare function-valued dimension, which is what makes the reference a reporting
#' choice: `center` names the level [coef()] reports against, and no refit is
#' needed to change it.
#'
#' @seealso [bartisan()] for the formula interface, [bcf()] for the causal case,
#'   [coef.bartisan_fit()] for reading the coefficients out, and
#'   [bartisan-families] for the order the forests come in.
#'
#' @examples
#' # The coefficient of `z` varies with `x1` and `x2`.
#' y ~ x1 + x2 + vc(z)
#'
#' # ... and with `x1` alone.
#' y ~ x1 + x2 + vc(z, ~ x1)
#'
#' # The effect of `z` varies across `z` itself, so the dose response is a
#' # curve rather than a line.
#' y ~ x1 + z + vc(z, ~ z + x1)
#'
#' @export
vc <- function(x, modifiers = NULL, center = "auto") {
  arg::err(c("{.fn vc} is a marker for a {.fn bartisan} formula and cannot be
              called on its own",
             i = "Write it in the formula, as in
                {.code y ~ x1 + x2 + vc(z)}."))
}

# Every `vc()` term of a formula, and the formula with them removed.
#
# The shape follows `split_random()`: the marked terms come out before anything
# else looks at the formula, and their variables are kept so the model frame
# carries them and they get the same missing-value handling as any predictor.
#
# Arguments are matched against `vc()`'s own formals with `match.call()`, so
# partial and positional matching behave the way they would in any other call,
# without the term being evaluated -- `vc(z)` names a variable and evaluating it
# would look for an object called `z`.
split_vc_terms <- function(formula) {
  rhs <- formula[[length(formula)]]
  parts <- formula_addends(rhs)
  marked <- vapply(parts, is_vc_call, logical(1L))

  if (!any(marked)) {
    # A `vc()` buried under something other than `+` is not a term this can
    # read, and silently treating it as a predictor would fit a different model.
    check_no_buried_vc(rhs)
    return(list(fixed = formula, vc = list()))
  }

  for (part in parts[!marked]) {
    check_no_buried_vc(part)
  }

  specs <- lapply(parts[marked], function(call) {
    label <- deparse1(call)

    # `match.call()` is what gives partial and positional matching without
    # evaluating the term. Its own message for a bad argument does not say which
    # term the argument was in, which is the one thing worth adding.
    matched <- tryCatch(match.call(vc, call), error = function(e) {
      arg::err(c("{.code {label}} is not a call {.fn vc} accepts:
                  {conditionMessage(e)}",
                 i = "{.fn vc} takes {.arg x}, {.arg modifiers} and
                    {.arg center}."))
    })

    if (is_null(matched[["x"]])) {
      arg::err("{.code {label}} must name the predictor whose coefficient
                varies, as in {.code vc(z)}")
    }

    covariate <- matched[["x"]]

    if (!is.symbol(covariate)) {
      arg::err(c("{.code {label}} must name a bare predictor, not
                  {.code {deparse(covariate)}}",
                 i = "Compute it in {.arg data} first, then name the column."))
    }

    modifiers <- matched[["modifiers"]]

    # Checked as an expression before it is evaluated. Evaluating first turns
    # `vc(z, x1)` -- a plausible slip, since the modifiers read like predictors
    # -- into `object 'x1' not found`, which says nothing about the real problem.
    if (!is_null(modifiers)) {
      if (!is.call(modifiers) || !identical(modifiers[[1L]], quote(`~`))) {
        arg::err("the modifiers of {.code {label}} must be a one-sided formula,
                  as in {.code vc({deparse(covariate)}, ~ x1 + x2)}")
      }

      modifiers <- eval(modifiers, envir = environment(formula))

      if (!rlang::is_formula(modifiers, lhs = FALSE)) {
        arg::err("the modifiers of {.code {label}} must be a one-sided formula,
                  as in {.code vc({deparse(covariate)}, ~ x1 + x2)}")
      }
    }

    center <- matched[["center"]] %or% "auto"

    if (!is.numeric(center)) {
      center <- eval(center, envir = environment(formula))
    }

    list(label = label,
         covariate = as.character(covariate),
         modifiers = modifiers,
         center = center)
  })

  names(specs) <- vapply(specs, `[[`, character(1L), "covariate")

  duplicated_at <- duplicated(names(specs))

  if (any(duplicated_at)) {
    arg::err("{.arg formula} gives {.val {unique(names(specs)[duplicated_at])}}
              a varying coefficient more than once")
  }

  # Rebuilt from the addends that are left, so an unmarked term keeps its
  # spelling -- `.`, `. - x1`, `a:b` and `I(x^2)` all survive untouched, which
  # matters because the design formula is expanded against the data later.
  keep <- parts[!marked]

  new_rhs <- {
    if (is_null(keep)) quote(1)
    else Reduce(function(a, b) call("+", a, b), keep)
  }

  fixed <- formula
  fixed[[length(fixed)]] <- new_rhs

  list(fixed = fixed, vc = specs)
}

is_vc_call <- function(expr) {
  is.call(expr) && identical(expr[[1L]], quote(vc))
}

# The top-level addends of a right-hand side. Only `+` is split, so `. - x1`,
# `a:b` and `I(x^2)` each stay whole, which is what lets them be put back
# verbatim.
formula_addends <- function(expr) {
  if (is.call(expr) && identical(expr[[1L]], quote(`+`)) && length(expr) == 3L) {
    return(c(formula_addends(expr[[2L]]), formula_addends(expr[[3L]])))
  }

  list(expr)
}

# `vc()` is a term of the formula, not a function that can appear inside one.
# `y ~ x1 - vc(z)` or `y ~ x1:vc(z)` would otherwise be read as a predictor
# called `vc(z)` and fitted as a different model without comment.
check_no_buried_vc <- function(expr) {
  if (!is.call(expr)) {
    return(invisible(NULL))
  }

  if (is_vc_call(expr)) {
    arg::err(c("{.code {deparse1(expr)}} must be added to the formula with
                {.code +}",
               i = "{.fn vc} names one of the model's terms; it cannot be
                  subtracted, crossed or nested inside another term."))
  }

  for (i in seq_along(expr)[-1L]) {
    if (!missing_arg_at(expr, i)) {
      check_no_buried_vc(expr[[i]])
    }
  }

  invisible(NULL)
}

# Whether a formula's right-hand side is, or contains, a bare `.`.
#
# This is what separates a slip from a choice in the overlap rule below: a
# caller who wrote `.` did not name the varying covariate, so removing it from
# the control function's modifiers is a correction rather than a change of
# model. A caller who wrote the name did name it.
uses_dot <- function(formula) {
  if (is_null(formula)) {
    return(FALSE)
  }

  any(all.vars(formula[[length(formula)]]) == ".") ||
    any(vapply(as.list(formula[[length(formula)]]), function(z) {
      identical(z, quote(.))
    }, logical(1L)))
}

# The modifiers of each forest, as term labels over the union design.
#
# `groups` is every predictor group the model has. The control function gets
# every group its own formula allows, less any varying covariate it may not see;
# each coefficient's forest gets what its own `modifiers` formula allows.
#
# Two rules, and they differ in whether the caller named the variable:
#
#   * The control function cannot split on a varying covariate. With `z` among
#     its modifiers, f_0(Z) + z f_1(Z) is not identified -- any function of `z`
#     moves between the two -- and no other implementation of this model checks.
#     Reached through `.`, the covariate is dropped without comment. Named
#     outright, the model is fitted as asked with a warning.
#   * A *categorical* covariate cannot split its own forest either. A level's
#     indicator is nonzero only on the rows where the level holds, and the
#     variable is constant on exactly those rows, so such a split separates rows
#     that contribute from rows that contribute nothing. Wasted rather than
#     unidentified, so it goes quietly.
vc_modifiers <- function(specs, groups, dot, categorical) {
  covariates <- vapply(specs, `[[`, character(1L), "covariate")
  present <- intersect(covariates, groups)

  control <- setdiff(groups, present)

  if (!is_null(present) && !dot) {
    arg::wrn(c("{.arg formula} lets the control function split on
                {length(present)} predictor{?s} whose coefficient varies:
                {.val {present}}",
               i = "The control function and those coefficients are then not
                  separately identified, because any function of a varying
                  covariate can move between them.",
               i = "Leave them out of the formula's fixed part to fit the
                  identified model."))
    control <- groups
  }

  slope <- lapply(specs, function(spec) {
    own <- spec[["covariate"]]

    allowed <- groups

    if (!is_null(spec[["modifiers"]])) {
      asked <- attr(stats::terms(spec[["modifiers"]]), "term.labels")

      # A name that is not a predictor is a typo, not a restriction, and
      # silently fitting a forest with fewer modifiers than were asked for
      # is the failure mode worth spending an error on. The covariate's own
      # name is the exception: it is legitimately absent from the design,
      # since `vc()` took it out.
      unknown <- setdiff(asked, c(groups, own))

      if (!is_null(unknown)) {
        arg::err(c("the modifiers of {.code {spec[['label']]}} name {length(unknown)} thing{?s} that {?is/are} not
                         {?a predictor/predictors} of this model: {.val {unknown}}",
                   i = "Its predictors are {.val {groups}}."))
      }

      allowed <- intersect(asked, groups)
    }

    # A numeric covariate modifying its own coefficient is a nonlinearity
    # in it, is identified, and is something a caller might mean.
    if (isTRUE(categorical[[own]])) {
      allowed <- setdiff(allowed, own)
    }

    allowed
  })

  list(control = control,
       covariates = covariates,
       slope = slope)
}

# The basis columns each varying coefficient multiplies, and the centring behind
# them.
#
# One column per forest, in the order the forests are built: a numeric covariate
# contributes one, a factor contributes one per level. `mf` is the model frame,
# so the covariate has already been through the same missing-value handling as
# any predictor.
#
# Centring is a reparameterization of the control function alone -- every
# coefficient and every estimand is identical under any choice -- so what it
# decides is what the control function *means* and how well the two forests mix.
# An uncentred binary covariate leaves the coefficient informed only by the rows
# where it is nonzero while the control function absorbs the rest, which is the
# correlation that makes the pair mix badly.
vc_basis <- function(specs, mf) {
  out <- lapply(specs, function(spec) {
    x <- mf[[spec[["covariate"]]]]

    if (is.factor(x) || is.character(x)) {
      return(vc_basis_factor(x, spec))
    }

    if (!is.numeric(x)) {
      arg::err("{.code {spec[['label']]}} needs a numeric or categorical
                predictor; {.val {spec[['covariate']]}} is
                {.cls {class(x)[1L]}}")
    }

    vc_basis_numeric(x, spec)
  })

  list(columns = do.call(cbind, lapply(out, `[[`, "columns")),
       parts = out)
}

vc_basis_numeric <- function(x, spec) {
  center <- spec[["center"]]

  arg::arg_or(
    center,
    arg::arg_string,
    arg::arg_number,
    .arg = sprintf("center of %s", spec[["label"]])
  )

  if (is.character(center)) {
    center <- switch(
      arg::match_arg(center, c("auto", "mean", "zero", "mid"),
                     .arg = sprintf("center of %s", spec[["label"]])),
      # Measured rather than assumed, because neither choice wins everywhere.
      # For a 0/1 covariate, zero is a value it actually takes and the control
      # function is then the surface among the untreated, which is the prognostic
      # score Hahn, Murray and Carvalho (2020) parameterize -- and it recovers
      # the coefficient better, at correlation 0.987 against 0.975 and RMSE 0.144
      # against 0.208 on the simulation in `_dev/varying-coefficients.md`. For
      # anything else zero may be nowhere near the data: with a covariate around
      # 50 the control function at zero is an extrapolation, and recovery
      # collapses to a correlation of 0.42 where mean-centring holds at 0.99.
      auto = if (all(stats::na.omit(x) %in% c(0, 1))) 0 else mean(x, na.rm = TRUE),
      mean = mean(x, na.rm = TRUE),
      zero = 0,
      mid = mean(range(x, na.rm = TRUE)))
  }

  columns <- matrix(x - center, ncol = 1L,
                    dimnames = list(NULL, spec[["covariate"]]))

  list(columns = columns, center = center, levels = NULL,
       kind = "numeric", scale = stats::sd(x, na.rm = TRUE))
}

# A factor gets one forest per level, symmetrically, the way `multinomial()`
# codes its predictors -- not K-1 contrasts against a reference level, which
# would shrink every level towards whichever one sorted first and give that one a
# different prior from the rest.
#
# The coding is over-parameterized by exactly one function. Mean-centred, the
# columns sum to zero across levels within a row, so adding any g(Z) to *every*
# level's forest changes the fit by g(Z) times zero -- the control function is
# not even involved. That is the same redundancy the symmetric multinomial coding
# carries, and it is handled the same way: a proper leaf prior, the 1/sqrt(2)
# scale correction, and recentring at reporting.
#
# The upshot is that the reference is a *reporting* choice here rather than a
# fitting one, which the numeric case does not get. So the fit is always
# mean-centred and `center` names the level to report against.
vc_basis_factor <- function(x, spec) {
  x <- as.factor(x)
  levels <- levels(x)

  if (length(levels) < 2L) {
    arg::err("{.code {spec[['label']]}} needs a predictor with at least two
              levels; {.val {spec[['covariate']]}} has
              {length(levels)}")
  }

  center <- spec[["center"]]

  # A factor has no numeric zero to sit at, so `"auto"` is the average level
  # composition.
  if (identical(center, "auto")) {
    center <- "mean"
  }

  if (is.numeric(center) || !center %in% c("mean", levels)) {
    arg::err(c("the center of {.code {spec[['label']]}} must be {.val mean} or
                one of {.val {levels}}",
               i = "{.val {spec[['covariate']]}} is categorical, so
                  {.val zero}, {.val mid} and a number do not name a value it
                  can take."))
  }

  indicators <- vapply(levels, function(l) as.numeric(!is.na(x) & x == l),
                       numeric(length(x)))
  indicators[is.na(x), ] <- NA_real_

  # Mean-centred, which puts the control function at the average composition of
  # the levels. Reporting against a particular level is exact from these draws.
  shares <- colMeans(indicators, na.rm = TRUE)
  columns <- sweep(indicators, 2L, shares, "-")
  colnames(columns) <- sprintf("%s%s", spec[["covariate"]], levels)

  list(columns = columns, center = center, levels = levels, shares = shares,
       kind = "factor", scale = rep.int(1, length(levels)))
}

# The formula the model frame is built from: every `vc(z, ...)` replaced by the
# bare `z`.
#
# The same division of labour `subbars()` and `nobars()` perform for random
# effects. The frame is built from this version, so the covariate comes along and
# gets the same missing-value handling as any predictor; the design is built from
# `split_vc_terms()$fixed`, so the covariate is not also a splitting predictor
# unless the caller asked for that.
vc_to_names <- function(expr) {
  if (!is.call(expr)) {
    return(expr)
  }

  if (identical(expr[[1L]], quote(vc))) {
    matched <- match.call(vc, expr)
    return(matched[["x"]])
  }

  for (i in seq_along(expr)[-1L]) {
    if (!missing_arg_at(expr, i)) {
      expr[[i]] <- vc_to_names(expr[[i]])
    }
  }

  expr
}

missing_arg_at <- function(expr, i) {
  identical(expr[[i]], quote(expr = ))
}

# The label of every forest, in the order they are built: the control function
# first, then each varying coefficient's, a factor contributing one per level.
#
# `parameter` is the distributional parameter these forests belong to, which is
# dropped from the labels when the family has only one -- so the common case is
# `(Intercept)` and the covariate's own name rather than `eta` and `eta:z`.
vc_forest_labels <- function(parameter, specs, parts, drop_parameter) {
  if (is_null(specs)) {
    return(parameter)
  }

  slopes <- unlist(lapply(seq_along(specs), function(j) {
    part <- parts[[j]]
    covariate <- specs[[j]][["covariate"]]

    if (identical(part[["kind"]], "factor")) {
      sprintf("%s%s", covariate, part[["levels"]])
    }
    else {
      covariate
    }
  }), use.names = FALSE)

  if (drop_parameter) {
    return(c("(Intercept)", slopes))
  }

  c(parameter, sprintf("%s:%s", parameter, slopes))
}

# Everything the engine and `predict()` need from the `vc()` terms of one
# formula: the centred basis, which predictor groups each forest may split on,
# and enough of the centring to rebuild the basis for new data.
resolve_vc <- function(specs, mf, design, dot) {
  groups <- unique(design[["term_labels"]][design[["assign"]]])

  if (length(specs) == 0L) {
    return(list(specs = list(), basis = NULL, parts = list(),
                slopes = 0L, masks = NULL, groups = groups))
  }

  missing_from_frame <- setdiff(vapply(specs, `[[`, character(1L), "covariate"),
                                names(mf))

  if (length(missing_from_frame) > 0L) {
    arg::err("{.arg data} has no column {.val {missing_from_frame}}")
  }

  basis <- vc_basis(specs, mf)

  categorical <- vapply(basis[["parts"]],
                        function(p) identical(p[["kind"]], "factor"),
                        logical(1L))
  names(categorical) <- vapply(specs, `[[`, character(1L), "covariate")

  modifiers <- vc_modifiers(specs, groups, dot, categorical)

  # One mask column per forest, in the order the forests are built: the control
  # function, then each coefficient's, a factor's levels consecutively. A
  # factor's levels share one modifier set, since they are one predictor.
  columns <- c(list(modifiers[["control"]]),
               unlist(lapply(seq_along(specs), function(j) {
                 rep(list(modifiers[["slope"]][[j]]),
                     length(basis[["parts"]][[j]][["scale"]]))
               }), recursive = FALSE))

  masks <- vapply(columns, function(allowed) groups %in% allowed,
                  logical(length(groups))) |>
    matrix(nrow = length(groups), ncol = length(columns),
           dimnames = list(groups, NULL))

  empty <- !apply(masks, 2L, any)

  if (any(empty)) {
    arg::err("{sum(empty)} of the model's forests {?has/have} no predictor left to split on")
  }

  list(specs = specs, basis = basis[["columns"]], parts = basis[["parts"]],
       slopes = ncol(basis[["columns"]]), masks = masks, groups = groups)
}

# The basis the fit was built with, for the entry points that evaluate the
# likelihood at stored draws. Empty when the model has no varying coefficient,
# which is what tells the engine to leave the family alone.
vc_stored_basis <- function(object) {
  basis <- object[["vc"]][["basis"]]

  if (is_null(basis)) matrix(0, 0L, 0L) else basis
}

# The basis for new data, under the centring the fit was built with.
#
# The centring values are the fit's, not the new data's: re-centring on
# `newdata` would move the control function's reference between the fit and the
# prediction, so the same covariate value would predict two different things.
# That matters most for the counterfactual grids `avg_comparisons()` builds,
# which hold every row at one covariate value.
vc_newdata_basis <- function(object, newdata) {
  vc <- object[["vc"]]

  if (is_null(newdata)) {
    return(vc[["basis"]])
  }

  columns <- lapply(seq_along(vc[["specs"]]), function(j) {
    spec <- vc[["specs"]][[j]]
    part <- vc[["parts"]][[j]]
    name <- spec[["covariate"]]

    if (!name %in% names(newdata)) {
      arg::err("{.arg newdata} has no column {.val {name}}, whose coefficient
                varies")
    }

    x <- newdata[[name]]

    if (identical(part[["kind"]], "factor")) {
      x <- factor(as.character(x), levels = part[["levels"]])

      if (anyNA(x) && !anyNA(newdata[[name]])) {
        arg::err("{.arg newdata} has {.val {name}} levels the model never saw:
                  {.val {setdiff(unique(as.character(newdata[[name]])),
                                 part[['levels']])}}")
      }

      indicators <- vapply(part[["levels"]],
                           function(l) as.numeric(!is.na(x) & x == l),
                           numeric(length(x)))
      indicators[is.na(x), ] <- NA_real_

      return(sweep(indicators, 2L, part[["shares"]], "-"))
    }

    matrix(as.numeric(x) - part[["center"]], ncol = 1L)
  })

  out <- do.call(cbind, columns)
  colnames(out) <- colnames(vc[["basis"]])
  out
}

# A factor's coefficients, recentered to sum to zero across its levels.
#
# The symmetric coding is over-parameterized by one function -- adding the same
# g(Z) to every level's forest changes nothing, because the centred indicators
# sum to zero within a row -- so the level forests are not individually
# identified and their raw draws are not meaningful on their own. Removing their
# mean across levels is the recentring that makes each one the deviation it is
# reported as, and it is exact rather than an approximation.
vc_recenter <- function(slopes, vc) {
  at <- 0L

  for (j in seq_along(vc[["specs"]])) {
    part <- vc[["parts"]][[j]]
    width <- length(part[["scale"]])
    columns <- at + seq_len(width)
    at <- at + width

    if (!identical(part[["kind"]], "factor")) {
      next
    }

    shared <- Reduce(`+`, slopes[columns]) / width

    for (k in columns) {
      slopes[[k]] <- slopes[[k]] - shared
    }
  }

  slopes
}
