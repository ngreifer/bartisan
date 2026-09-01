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
#'   both how the model is fitted and how it is reported. `"mean"`, the default,
#'   centers `x` so that the control function is the surface at its average;
#'   `"zero"` leaves `x` alone, so the control function is the surface at
#'   `x = 0`; `"mid"` uses the midpoint of `x`'s range, which is `0.5` for a
#'   binary predictor; a number uses that number. `"estimate"` draws the coding
#'   rather than fixing it, following Hahn et al. (2020) section 5.3. For a
#'   factor, `center` is `"mean"` or the name of a level.
#'
#' @returns
#' A list of class `bartisan_vc`, which is of no use outside a formula.
#'
#' @seealso [bartisan()] for the formula interface and [bartisan-families] for
#'   the order the forests come in.
#'
#' @examples
#' # The coefficient of `z` varies with `x1` and `x2`.
#' y ~ x1 + x2 + vc(z)
#'
#' # ... and with `x1` alone.
#' y ~ x1 + x2 + vc(z, ~ x1)
#'
#' @export
vc <- function(x, modifiers = NULL, center = "mean") {
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

    center <- matched[["center"]] %or% "mean"

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
    if (length(keep) == 0L) quote(1)
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

  if (length(present) > 0L && !dot) {
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

  list(control = control, covariates = covariates,
       slope = lapply(specs, function(spec) {
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

           if (length(unknown) > 0L) {
             arg::err(c("the modifiers of {.code {spec[['label']]}} name
                         {length(unknown)} thing{?s} that {?is/are} not
                         {?a predictor/predictors} of this model:
                         {.val {unknown}}",
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
       }))
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

  if (is.character(center)) {
    center <- switch(
      arg::match_arg(center, c("mean", "zero", "mid", "estimate"),
                     .arg = sprintf("center of %s", spec[["label"]])),
      # Drawn rather than fixed, so the column is the covariate as written and
      # the coding coefficients carry the reference; see `vc_coding()`.
      estimate = 0,
      mean = mean(x, na.rm = TRUE),
      zero = 0,
      mid = mean(range(x, na.rm = TRUE)))
  }

  arg::arg_number(center, .arg = sprintf("center of %s", spec[["label"]]))

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

  list(columns = columns, center = center, levels = levels,
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
  if (length(specs) == 0L) {
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
                  logical(length(groups)))
  masks <- matrix(masks, nrow = length(groups), ncol = length(columns),
                  dimnames = list(groups, NULL))

  empty <- which(!apply(masks, 2L, any))

  if (length(empty) > 0L) {
    arg::err("{length(empty)} of the model's forests {?has/have} no predictor
              left to split on")
  }

  list(specs = specs, basis = basis[["columns"]], parts = basis[["parts"]],
       slopes = ncol(basis[["columns"]]), masks = masks, groups = groups)
}
