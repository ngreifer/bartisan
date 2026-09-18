#' Fit a generalized Bayesian additive regression trees (BART) model
#'
#' @description
#' Fits a BART model in which the response distribution is arbitrary rather than
#' restricted to the conditionally conjugate cases, using the
#' Laplace-approximation reversible-jump sampler of Linero (2025). Decision rules
#' may be soft, as in Linero and Yang (2018), which gives smoother fits than the
#' step functions of standard BART. The interface mirrors that of [stats::glm()]:
#' a formula, a data frame, and a family, with the families that `glm()` has no
#' counterpart for documented at [bartisan-families].
#'
#' @param formula a model formula. The right-hand side lists candidate
#'   predictors; the model finds interactions and nonlinearity on its own, so
#'   `y ~ x1 + x2 + x3` is usually the right specification. Survival families take a
#'   \pkgfun{survival}{Surv} object on the left. A `(1 | group)` term adds a
#'   group-level random intercept, in the notation of \pkg{lme4}; see Details.
#'
#'   For a family with more than one additive predictor this may be a *list* of
#'   formulas, one per forest, to give each one its own predictors. The first is
#'   the model for the main parameter and carries the response; the rest need no
#'   response, and follow the order in [bartisan-families], under "Several
#'   additive predictors", which also gives the name of each forest so the list
#'   can be named instead of ordered:
#'
#'   ```r
#'   bartisan(list(y ~ x1 + x2, ~ x2 + x3), data = d,
#'            family = gaussian_ls())
#'   bartisan(list(mean = y ~ x1 + x2, log_sd = ~ x2), data = d,
#'            family = gaussian_ls())
#'   ```
#'
#'   One formula applies to every forest, which is the ordinary case. A predictor
#'   left out of one forest's formula is still in the data and is never split on
#'   by that forest.
#'
#'   A formula naming no predictor at all makes that parameter a constant.
#'   `~ 1` leaves its forest nothing to split on, so every tree in it is a stump
#'   and the forest is a single drawn scalar rather than a function of the
#'   predictors.
#'
#'   [vc()] terms are read out of each formula in turn, so a parameter has the
#'   varying coefficients its own formula asks for and no others, which makes the
#'   forests two-dimensional (one axis the parameter, the other the coefficient).
#'   [vc()] documents how they are then named and keyed.
#' @param data a data frame containing the variables named in `formula`.
#' @param family the response distribution, given as a [stats::family] object,
#'   as one of the families in [bartisan-families], or as the name of either. A
#'   `family` object is accepted when the distribution it names is one this
#'   package implements, since the likelihood is the package's rather than the
#'   object's: a `family` object carries a link and a variance function and not a
#'   density, so one naming anything else (e.g., [stats::inverse.gaussian], or a
#'   Tweedie from another package) is an error rather than something a likelihood
#'   can be built from, and [custom_family()] is the route for those. Links are
#'   used as supplied, and a link the package does not compile is composed onto
#'   the scale its family works on. Default is `NULL`, in which case the family
#'   is read off the response and a message reports the choice; see Details for
#'   the rules.
#' @param prior_only `logical`; whether to draw from the prior rather than the
#'   posterior, which is what a prior predictive check reads. Default is `FALSE`.
#'   Not available for every family; see Details.
#' @param weights optional; prior weights, one per observation. For a binomial
#'   response given as proportions, these are the numbers of trials, as in
#'   `glm()`.
#' @param offset optional; a known component of the additive predictor, on the
#'   link scale.
#' @param subset optional; a vector specifying the subset of rows to use.
#' @param na.action how missing values are handled. Default is [stats::na.pass],
#'   which keeps rows whose *predictors* are missing and lets the splitting rules
#'   decide where they go, which is something the trees can do and `lm()` and
#'   `glm()` cannot; see Details. Pass [stats::na.omit] to drop any row with a
#'   missing value anywhere instead. Note that rows with a missing response,
#'   weight, or offset are dropped either way, with a warning, since there is
#'   nothing to fit them to.
#' @param control a `<bartisan_control>` object; the output of a call to
#'   [bartisan_control()], containing the sampler and prior settings.
#' @param ... further arguments to [bartisan_control()], which are merged into
#'   `control` and override any value given there, so that
#'   `bartisan(..., num_trees = 20)` and
#'   `bartisan(..., control = bartisan_control(num_trees = 20))` are the same
#'   call. A name that is not an argument of [bartisan_control()] is an error
#'   rather than being silently ignored.
#'
#' @details
#' ## The Sampler
#'
#' Standard BART relies on the leaf parameters being integrable in closed form,
#' which restricts it to a Gaussian response, or to models that can be reduced to
#' one by data augmentation. Linero's (2025) algorithm removes that restriction. At each
#' candidate move it builds a Gaussian approximation to the conditional posterior
#' of the affected leaf parameters, by Fisher scoring, and uses that
#' approximation as the proposal in a reversible-jump Metropolis step. The
#' approximation only has to be good enough to be accepted often; the stationary
#' distribution is the exact posterior either way.
#'
#' What a new family therefore has to supply is only the log density of one
#' observation and its first two derivatives with respect to the additive
#' predictor. Families whose response has more than one unconstrained parameter,
#' such as `multinomial()` and `gaussian_ls()`, carry one forest per
#' parameter. Because that is the whole interface, it can be reached from R:
#' [custom_family()] takes the log density as an R function and differences it for
#' the derivatives.
#'
#' ## Inferring the Family
#'
#' `family` may be left alone, in which case it is read off the response:
#'
#' | Response | Family |
#' |---|---|
#' | [survival::Surv()] object, or a two-column matrix of times and events | `dpm_aft()` |
#' | ordered factor | `ordinal()` |
#' | logical, or two levels, or numeric zeros and ones | `binomial()` |
#' | factor or character with more than two levels | `multinomial()` |
#' | two-column matrix of successes and failures | `binomial()` |
#' | anything else | `dpm()` |
#'
#' A message reports the choice, and naming `family` is what silences it, which
#' is also what changes it.
#'
#' Two scenarios are worth noting. A count is read as a numeric variable and therefore has `dpm()` as its default. And a
#' numeric response with exactly two values other than zero and one (e.g.,
#' `c(1, 2)`) is also given `dpm()` as its default rather than `binomial()`, which of the two counts as the
#' success not being something to guess at.
#'
#' ## Soft Decision Rules
#'
#' By default a decision rule is a smooth gate rather than a step, so an
#' observation reaches every leaf with some weight and the fitted function is
#' smooth. `gate` in [bartisan_control()] chooses both whether the rules are soft
#' and, if they are, the gate's shape; the default is the bounded `"smoothstep"`,
#' and `"logistic"` is Linero and Yang's (2018) original. Soft rules cost more
#' per iteration, since a leaf now touches every observation rather than only the
#' ones inside its cell, and they make the leaf parameters of a tree dependent on
#' one another. Combining them with a non-conjugate likelihood is an extension of
#' Linero (2025), which leaves it as an open problem; it is handled here by
#' giving the reversible-jump move a bivariate Laplace proposal for the pair of
#' child leaves, which reduces to Linero's independent pair exactly when the
#' rules are hard. Set `gate = "hard"` in [bartisan_control()] for the faster
#' hard-rule sampler.
#'
#' ## Random Intercepts
#'
#' A `(1 | group)` term in the formula adds an intercept per level of `group`,
#' drawn from a common mean-zero normal whose standard deviation is itself drawn
#' under the same half-Cauchy prior the leaf scale uses. Several grouping factors
#' are allowed, and `(1 | a/b)` expands to nesting as it does in \pkg{lme4}:
#'
#' ```r
#' bartisan(y ~ x1 + x2 + (1 | school), data = d)
#' bartisan(y ~ x1 + (1 | school) + (1 | year), data = d)
#' ```
#'
#' The intercepts are in `fit$ranef` and their standard deviations in `fit$tau`,
#' one matrix per additive predictor. A family with several predictors gets a
#' separate set for each (i.e., a zero-inflated count model has a group effect on
#' the count part and another on the inflation part), and they are independent of
#' one another.
#'
#' Only random *intercepts* are supported, and a random slope is refused rather
#' than ignored. The reason is that a random intercept is a scalar entering the
#' predictor with weight one for the observations in its level, which is what a
#' leaf is once its gate is removed, so the sampler's leaf machinery handles it
#' exactly; a slope is a different shape of parameter. A variable whose effect
#' varies by group belongs in the fixed part of the formula, where a tree can
#' split on the group and on the variable together and get an interaction of any
#' shape.
#'
#' A grouping factor can also go in the fixed part, where a tree splits on it
#' like anything else, and with few large groups that is the better choice: the
#' group means are well determined without pooling and a split can interact the
#' group with the covariates. The random intercept wins where there are many
#' small groups, which is where partial pooling earns its keep.
#'
#' A level of `group` that was not present at fitting time is given the prior
#' mean of zero when predicting, with a warning.
#'
#' ## Missing Predictor Values
#'
#' A missing predictor is not imputed and its row is not dropped, which is the
#' default here because a tree can do something better with a missing value than
#' either. Instead each splitting rule carries the answer for itself. A
#' rule on a variable that has missing values is drawn as one of three, with
#' equal probability:
#'
#' - `x < c`, or missing, goes left;
#' - `x < c` goes left, missing goes right;
#' - missing goes left, present goes right.
#'
#' This is missingness incorporated in attributes, and the third rule is what
#' lets the model split on missingness itself, so a variable whose absence
#' carries the signal is usable even where its observed values say nothing.
#'
#' Two consequences are worth being clear about. `predict()` accepts missing
#' values in a column that had them at fitting time, those being the columns
#' whose rules carry an answer. And what the model estimates is the mean of the
#' response given the predictors and the pattern of missingness, which is the
#' quantity prediction calls for; where the estimand is a regression or causal
#' effect defined on complete data, multiple imputation is the right tool.
#'
#' ## Preprocessing
#'
#' Predictors are mapped to the unit interval, because the cutpoint prior is
#' uniform on a node's live range and the soft-rule bandwidth is measured on the
#' predictor scale. Factors share a
#' single weight in the sparsity prior, so that a factor is selected or not as a
#' whole rather than one level at a time. The additive predictor starts from an
#' intercept-only fit, so the leaf prior describes departures from that fit
#' rather than the absolute level of the response. That starting value is the
#' exact null-model estimate for most families; for the accelerated failure time
#' families, where censoring makes the sample mean of the log times biased, and
#' for the zero-inflated and ordered beta families, it is a moment
#' approximation, which the sampler then moves away from.
#'
#' ## Drawing From the Prior (`prior_only`)
#'
#' `prior_only = TRUE` fits the same model to no data. Every observation is given
#' a weight of zero, and since the weight multiplies that observation's log
#' density, its gradient, and its curvature, the likelihood is flat: each tree
#' move is accepted or rejected on the prior alone and each leaf is drawn from
#' its prior. A family that draws an auxiliary parameter from the response
#' directly rather than through the weighted density (the mixture atoms under
#' `dpm()`, the latent utilities under `multinomial(link = "probit")`) is told
#' separately that a weightless observation carries no information, and draws
#' that parameter from its own prior instead. The sampler is otherwise untouched,
#' so what comes back is an ordinary fit whose draws are prior draws, and
#' \pkgfun{rstantools}{posterior_predict} on it gives the prior predictive
#' distribution.
#'
#' It answers a question the priors themselves cannot. `k`, `gamma` and `beta`
#' are statements about trees and leaves, and what they
#' imply about an outcome is opaque. The replicates put that on the response's own scale,
#' where it can be judged: a prior predictive that puts its mass where the
#' outcome cannot go, or spread over an implausible range, is a prior worth
#' changing before the data are seen and before any of their information is
#' spent.
#'
#' The additive predictor is anchored at an intercept-only fit on the link scale
#' and the leaf scale is calibrated from the response, which is how a BART prior
#' is specified. The replicates therefore take their location and scale from the
#' response and everything else from the prior: which predictors are split on,
#' how deep, how far the fitted function departs from that anchor. Read them for
#' shape and spread rather than for level, and note that the wider a prior is the
#' wider its replicates, so `gaussian_ls()` and `Gamma_ls()`, which put a log
#' scale in a second forest, run far wider than the response ever does.
#'
#' `loo()`, `waic()` and `kfold()` score a fit against data, so they refuse a
#' prior-only fit, and \pkgfun{performance}{model_performance} leaves out the
#' three columns built on them.
#'
#' @returns
#' A `<bartisan_fit>` object, a list with the following components among others.
#'
#'   \item{`eta`}{a list with one matrix per additive predictor, each of
#'     posterior draws by observation, on the link scale.}
#'   \item{`fitted`}{fitted values on the response scale, averaged over draws.}
#'   \item{`counts`}{a list with one matrix per additive predictor, of the
#'     number of splitting rules using each predictor group in each draw. Useful
#'     for variable selection.}
#'   \item{`aux`}{draws of the nuisance parameters, such as the residual
#'     standard deviation or the ordinal cutpoints, when the family has any.}
#'   \item{`has_na`}{which predictor columns contained a missing value, which is
#'     what determines where `predict()` will accept one.}
#'   \item{`sigma_mu`, `bandwidth`}{draws of the leaf standard deviation and,
#'     for soft rules, the per-tree gate bandwidths.}
#'   \item{`loglik`}{the log likelihood at each draw.}
#'   \item{`control`}{the `<bartisan_control>` object the fit used, with any
#'     settings given in `...` merged in.}
#'
#' @references
#' Linero, A. R. (2025). Generalized Bayesian additive regression trees models:
#' beyond conditional conjugacy. *Journal of the American Statistical
#' Association*, 120(549), 356--369. \doi{10.1080/01621459.2024.2337156}
#'
#' Linero, A. R., & Yang, Y. (2018). Bayesian regression tree ensembles that
#' adapt to smoothness and sparsity. *Journal of the Royal Statistical Society
#' Series B*, 80(5), 1087--1110. \doi{10.1111/rssb.12293}
#'
#' @seealso
#' [bartisan_control()] for the sampler and prior settings;
#' [predict.bartisan_fit()] for prediction; [bartisan-families] for the
#' likelihoods, and `vignette("families")` for a family-by-family guide;
#' [bartisan-marginaleffects] for reading effects off a fit
#'
#' @examples
#' data("rhc")
#' set.seed(123)
#'
#' # Whether a patient died, with every other variable a candidate predictor
#' # and the family read off the response. `days` is the timing of the same
#' # event, so it is excluded rather than conditioned on
#' fit <- bartisan(death ~ . - days, data = rhc,
#'                 num_trees = 10, num_burn = 50, num_draws = 50,
#'                 verbose = FALSE)
#' fit
#'
#' # Fitted probabilities
#' head(predict(fit, type = "response"))
#'
#' # The forest has no coefficients, so an effect is a contrast of
#' # predictions, here of catheterization on the probability of death
#' if (rlang::is_installed("marginaleffects")) {
#'   marginaleffects::avg_comparisons(fit, variables = "rhc")
#' }
#'
#' @export
bartisan <- function(formula, data, family = NULL, weights = NULL,
                     offset = NULL, subset = NULL,
                     na.action = stats::na.pass,
                     control = bartisan_control(), prior_only = FALSE, ...) {

  cl <- match.call()

  arg::arg_flag(prior_only)

  warn_unoptimized()

  # One formula per forest. A single formula is the common case and comes back as
  # a list of one, so everything below is written once.
  #
  # `one_sided = FALSE` on the first rather than a plain formula check: a formula
  # with no response has nothing to fit, and catching it here is what keeps the
  # `.` expansion below from inventing one. `update(~ x1 + x2, . ~ .)` returns
  # `. ~ x1 + x2`, so the frame would go looking for a variable named `.` and the
  # caller would get `object '.' not found` instead of the real problem.
  # A bare formula is checked here; a list is checked inside
  # `split_formula_list()`, which knows that a named one may carry the response
  # on any element and so cannot just look at the first.
  if (rlang::is_formula(formula)) {
    arg::arg_formula(formula, one_sided = FALSE)
  }

  forest_formulas <- split_formula_list(formula)
  if (!inherits(control, "bartisan_control")) {
    arg::err("{.arg control} must be the result of {.fn bartisan_control}")
  }

  control <- merge_control(control, list(...))

  # Read after the merge, so that `chains` given in `...` reaches this the same
  # way every other setting does.
  chains <- control[["chains"]]

  # `family` is resolved after the model frame, not before, because the default
  # is read off the response and the response is not available until then.

  # Random-effect terms come out of the formula before anything else looks at
  # it. The model frame is built from the version with the bars replaced by
  # ordinary sums, so the grouping variables come along and get the same
  # missing-value handling as everything else; the design matrix is built from
  # the version with the bars removed, so they are not also predictors.
  # The frame is built from every predictor any forest uses. Each forest is held
  # to its own subset further down, by zeroing its splitting weights on the terms
  # its formula leaves out, rather than by carrying a design matrix of its own.
  if (!rlang::is_formula(formula)) {
    formula <- union_formula(forest_formulas,
                             if (!missing(data) && is.data.frame(data)) data)
  }

  split <- split_random(formula)

  # `vc()` terms come out of the barless fixed part, so `terms()` never has to
  # make sense of a `|`. What is left is the design: the covariate whose
  # coefficient varies is not also a splitting predictor unless the caller put it
  # there.
  #
  # This is the union across every forest's formula, and it settles the design
  # alone -- a covariate that one parameter wraps in `vc()` and another names
  # outright stays a design column, because the union names it both ways and only
  # the wrapped mention is removed. Which forest may split on what is settled
  # later, per parameter, from that parameter's own formula.
  vc_split <- split_vc_terms(split[["fixed"]], unique_covariates = FALSE)
  vc_any <- !is_null(vc_split[["vc"]])
  split[["fixed"]] <- vc_split[["fixed"]]

  mf <- match.call(expand.dots = FALSE)
  keep <- match(c("formula", "data", "subset", "weights", "offset"),
                names(mf), 0L)
  mf <- mf[c(1L, keep)]
  # The frame carries the varying covariates, so they get the same missing-value
  # handling as any predictor; the design does not.
  mf[["formula"]] <- reformulas::subbars(vc_to_names(formula))

  # Unused levels are dropped from the predictors below rather than here,
  # because `model.frame()` would drop them from the *response* too. An ordinal
  # model can estimate a threshold for a category nobody landed in -- the
  # categories either side of it inform the two thresholds that bound it -- and
  # a rating scale with an unselected point is the ordinary case, not a
  # degenerate one.
  mf[["drop.unused.levels"]] <- FALSE

  # Set from the formal rather than carried over from the call. `match.call()`
  # only records what the caller actually wrote, so an argument left at its
  # default is absent from the reconstructed call and `model.frame()` falls back
  # on its own default -- which is `getOption("na.action")`, usually `na.omit`.
  # That silently overrode this function's default of `na.pass`. `NULL` is left
  # out so that it still means "whatever the session's option says", which is
  # what `lm()` and `glm()` do with it.
  if (!is_null(na.action)) {
    mf[["na.action"]] <- na.action
  }

  # `model.frame()` keeps every variable named anywhere in the formula, so a term
  # removed with `-` survives in the frame even though the terms correctly drop
  # it. That leaves a column in the frame that is not a predictor, and the
  # packages that read the frame to find out what the model uses then treat it as
  # one: `avg_comparisons()` reported an effect for a variable the model had
  # never seen. Resolving the formula through its own terms first expands `.` and
  # carries out the subtraction, so what reaches `model.frame()` names exactly
  # the variables the model uses. `.` can only be expanded when there is a
  # data frame to expand it against.
  # Left alone when the formula has random-effect terms: `.` would then expand
  # over the grouping variables as well, and the fixed part is derived from the
  # frame further down on the assumption that it has not been rewritten.
  # The frame's formula and the design's are expanded separately, because they
  # differ: `y ~ . + vc(z)` puts `z` in the frame either way, and whether it is
  # also a splitting predictor is what `.` decides.
  if (!missing(data) && is.data.frame(data) && is_null(split[["bars"]])) {
    resolved <- stats::update(stats::terms(mf[["formula"]], data = data), . ~ .)
    environment(resolved) <- environment(mf[["formula"]])
    mf[["formula"]] <- resolved

    design_formula <- stats::update(
      stats::terms(split[["fixed"]], data = data), . ~ .)
    environment(design_formula) <- environment(split[["fixed"]])
    split[["fixed"]] <- design_formula

    # Stored as well as used. `insight::find_formula()` reads this field, and
    # from `death ~ . - days` it concluded that the model's one predictor was
    # `days`, which is the variable the formula removes. Everything built on
    # that -- `avg_comparisons()` most visibly -- then described the wrong model.
    # It is the *frame's* formula that is stored, so a varying covariate is a
    # variable of the model as far as the estimand packages are concerned, which
    # is what makes `avg_comparisons(fit, variables = "z")` work.
    formula <- resolved
  }

  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # An indicator column that can never fire is only a cost, so the predictors
  # lose their unused levels; the response keeps them.
  response_at <- attr(attr(mf, "terms"), "response")

  for (j in seq_along(mf)) {
    if (j != response_at && is.factor(mf[[j]])) {
      mf[[j]] <- droplevels(mf[[j]])
    }
  }

  # The terms of the fixed part, which is what the trees split on. The frame's
  # own terms include the grouping variables, because the frame was built from
  # the formula with the bars replaced by sums so that those variables would be
  # present and get the same missing-value handling; using them would make a
  # grouping factor a predictor as well.
  #
  # `.` expands against the columns of whatever is given as `data`, so the frame
  # is handed over without the two columns `model.frame()` adds for the weights
  # and the offset -- those are not predictors, and expanding over them would put
  # them in the terms and send `predict()` looking for them in `newdata`.
  mt <- stats::terms(split$fixed, data = mf[!startsWith(names(mf), "(")])

  # The classes come from the frame, which knows them, rather than from the
  # terms, which do not; `predict()` uses them to refuse `newdata` that has
  # turned a factor into something else.
  attr(mt, "dataClasses") <- attr(attr(mf, "terms"), "dataClasses")

  if (is_null(attr(mt, "term.labels"))) {
    arg::err(c("{.arg formula} must include at least one predictor a forest can
                split on",
               i = if (vc_any) {
                 "A {.fn vc} covariate is what a coefficient multiplies, not
                  something its own forest can split on, so a model of nothing
                  but {.fn vc} terms has no predictors left."
               }))
  }

  # With na.action = na.pass the caller is asking for missing predictors to be
  # kept and handled by the splitting rules. A missing response, weight or
  # offset is a different matter: there is nothing to fit those rows to.
  mf <- drop_unusable_rows(mf, mt)

  y <- stats::model.response(mf, "any")
  model_weights <- as.vector(stats::model.weights(mf))
  model_offset <- as.vector(stats::model.offset(mf))

  family <- as_bartisan_family(family %or% default_family(y, model_weights))

  design <- build_design(mt, mf)
  n <- nrow(design$x)

  if (n == 0L) {
    arg::err("no usable observations remain")
  }

  unit <- unit_transform(design$x, control[["x_transform"]])

  has_na <- vapply(seq_len(ncol(unit$x)), function(j) anyNA(unit$x[, j]),
                   logical(1L)) |>
    setNames(colnames(unit$x))

  random <- random_terms(split$bars, mf)

  response <- prepare_response(family, y, model_weights, model_offset,
                               unit$x, n)

  group_probs <- make_group_probs(design$assign, design$term_labels)

  # A predictor group whose columns are mutually exclusive indicators gets a
  # level code per observation, so a rule on it can name a subset of its levels
  # rather than a threshold on one indicator. Taken from the design matrix rather
  # than the unit-mapped one, since the mapping sends a two-valued column to
  # exactly 0 and 1 and leaves the indicators alone either way.
  levels_info <- level_codes(design$x, design$assign)

  # `categorical = "onehot"` is expressed by telling the engine that no group has
  # levels, which sends every rule down the threshold-on-one-column path. The
  # codes are still built and stored, so `predict()` needs no second case.
  if (identical(control[["categorical"]], "onehot")) {
    levels_info[["n_levels"]] <- integer(length(levels_info[["n_levels"]]))
  }

  engine_control <- as.list(control)

  # The names of the family's additive predictors, which is what one formula per
  # forest is keyed by. Varying coefficients then split each of these into a
  # control function and its coefficients, so there are two label sets from here
  # on: these, for the formulas, and the enlarged set below for everything that
  # is set per forest.
  #
  # The trailing pinned forests standing in for a custom family's nuisance
  # parameters are in neither: they are not additive predictors and nothing
  # about them is the caller's to set per forest.
  n_param <- response[["n_forest"]] - response[["n_aux"]]
  param_labels <- forest_labels(response[["family"]], response[["opts"]],
                                response[["levels"]], n_param)
  joint <- joint_forests(response[["family"]])

  if (length(forest_formulas) > 1L && joint) {
    arg::err(c("{.arg formula} must be a single formula for this family, not
                {length(forest_formulas)}",
               i = "Its {length(param_labels)} forests are the levels of one
                  parameter and act together, so they take the same
                  predictors."))
  }

  if (length(forest_formulas) > n_param) {
    arg::err(c("{.arg formula} has {length(forest_formulas)} formulas but this
                family has {n_param} forest{?s}",
               i = "Its forests are {.val {param_labels}}."))
  }

  # Names on the list of formulas are checked the same way any per-forest
  # argument's are, and reorder it, so that `list(log_sd = ~ x2, mean = y ~ x1)`
  # means what it says.
  if (!is_null(names(forest_formulas))) {
    forest_formulas <- resolve_per_forest(forest_formulas, param_labels,
                                          "formula",
                                          default = forest_formulas[[1L]],
                                          joint = joint)
  }

  # One formula given applies to every forest, which is the rule for every
  # per-forest argument and is what makes the ordinary single-formula call reach
  # the engine unchanged. Fewer formulas than forests, but more than one, is not
  # a recycling anyone would mean.
  if (length(forest_formulas) == 1L) {
    forest_formulas <- rep(forest_formulas, n_param)
  }
  else if (length(forest_formulas) < n_param) {
    arg::err(c("{.arg formula} has {length(forest_formulas)} formulas but this
                family has {n_param} forests",
               i = "Give one formula, or {n_param}, or name them:
                  {.val {param_labels}}."))
  }

  # `vc()` terms come out of each parameter's own formula, so the mean can have a
  # varying coefficient the log standard deviation does not -- and a single
  # formula gives every parameter the same ones, which is the recycling rule
  # every other per-forest argument follows. Read after the reordering above so
  # that a named list lines its `vc()` terms up with the right parameter.
  forest_vc <- lapply(forest_formulas, function(f) {
    fixed <- split_random(f)[["fixed"]]

    list(specs = split_vc_terms(fixed)[["vc"]],
         dot = uses_dot(fixed))
  })

  # The formulas the masks are built from carry no `vc()` terms: a covariate
  # whose coefficient varies is what a coefficient multiplies, not something a
  # forest splits on.
  forest_fixed <- lapply(forest_formulas, function(f) {
    parts <- split_random(f)
    out <- split_vc_terms(parts[["fixed"]])[["fixed"]]
    environment(out) <- environment(f)
    out
  })

  # The basis each varying coefficient multiplies, centered, plus which forests
  # may split on what. Built from the frame rather than the design, because the
  # covariate is deliberately not a design column unless the caller put it there.
  vc <- resolve_vc(forest_vc, mf, design,
                   forest_masks(forest_fixed, colnames(group_probs),
                                if (!missing(data) &&
                                    is.data.frame(data)) data,
                                response_of(forest_fixed)),
                   response[["n_aux"]],
                   if (n_param > 1L) param_labels)

  if (vc[["slopes"]] > 0L && joint) {
    arg::err(c("{.fn vc} is not available for this family.",
               i = "Its {n_param} forests are the levels of one parameter and are
                  identified only up to a shared function, which is removed when
                  they are reported. A coefficient forest per level would add one
                  such direction per coefficient, and the reporting does not
                  carry them."))
  }

  response <- expand_for_vc(response, vc, y, response[["intercept"]])
  response[["offset"]] <- build_offset(response[["intercept"]], model_offset,
                                       response[["n_forest"]], n)
  response[["intercept"]] <- NULL

  # The enlarged set: one label per forest the engine builds, which is what every
  # per-forest setting below is keyed by.
  labels <- forest_labels(response[["family"]], response[["opts"]],
                          response[["levels"]],
                          response[["n_forest"]] - response[["n_aux"]], vc)

  # Which additive predictor each forest feeds, and which basis column it is
  # multiplied by. Zero-based for the engine, with -1 for a control function --
  # and for a custom family's pinned nuisance forests, which carry no
  # coefficient. This is the whole of what the engine needs to know about the
  # shape: the family reads it to combine the forests, and the random-effect
  # builder reads it to give group intercepts to control functions only.
  if (vc[["slopes"]] > 0L) {
    response[["opts"]] <- c(response[["opts"]],
                            list(vc_param = as.integer(vc[["param"]]) - 1L,
                                 vc_column = as.integer(vc[["column"]]) - 1L,
                                 vc_labels = labels))
  }

  # The coding coefficients ride along in the family's options, since it is the
  # family that draws them: they are a nuisance parameter of the same kind as a
  # scale, drawn once a sweep on the same hook.
  coding <- vc_coding(vc)

  if (!is_null(coding)) {
    response[["opts"]] <- c(response[["opts"]],
                            list(vc_coding = coding[["codes"]],
                                 vc_coding_levels = coding[["levels"]],
                                 vc_coding_names = coding[["names"]]))
  }

  # Matched against the predictors here rather than in `bartisan_control()`,
  # which does not know them. Each is one value per predictor group, so a
  # factor's dummy columns share what the term was given, the way they already
  # share one entry of the sparsity prior, and one column per forest.
  #
  # `split_prior` is weights the caller fixed; `split_mask` is which predictors
  # each forest's own formula lets it split on. Separate, because the first says
  # nothing may be drawn and the second only says over what.
  split <- resolve_split_matrix(control[["split_prior"]], colnames(group_probs),
                                labels, vc[["masks"]],
                                response[["n_forest"]])

  engine_control[["split_prior"]] <- split[["prior"]]
  engine_control[["split_mask"]] <- split[["mask"]]

  engine_control[["gate"]] <- gate_code(control[["gate"]])

  # The rest of the settings the engine keeps one copy of per forest. Each is
  # spread to one value per forest here, so the engine never has to decide what a
  # scalar means, and each keeps its own default where a named argument left a
  # forest out. `k` is not among them: it is a way of writing `sigma_mu`, and
  # that is spread just below.
  # The default for a forest the caller did not name is the argument's own
  # default, not the first value they did give. Passing `control[[nm]][[1L]]`
  # here made `gamma = c(log_sd = 0.5)` give the mean forest 0.5 as well, which
  # is borrowing another forest's value -- the one thing `?bartisan_control`
  # promises this does not do. A scalar still spreads to every forest, and a
  # positional vector is still taken in order; only the partly-named case moved.
  for (nm in names(PER_FOREST_DEFAULTS)) {
    engine_control[[nm]] <- per_forest_vector(
      control[[nm]], labels, nm, PER_FOREST_DEFAULTS[[nm]], joint)
    engine_control[[nm]] <- rep(engine_control[[nm]],
                                length.out = response[["n_forest"]])
  }

  # A forest whose formula names no predictor is intercept-only, and a branching
  # probability of zero is how the engine holds every tree in a forest at a
  # single leaf. The forest is then one drawn scalar, which is what `~ 1` says.
  # This overrides whatever `gamma` the caller gave that forest, because there is
  # nothing for it to branch on either way.
  if (any(vc[["pinned"]])) {
    engine_control[["gamma"]][vc[["pinned"]]] <- 0
  }

  # One tree count per additive predictor. A scalar is recycled, so the common
  # case reads the same as before; a vector, or one keyed by the forest names,
  # lets a forest that needs less capacity be given less, which is most of what
  # makes `gaussian_ls()` affordable.
  engine_control[["num_trees"]] <- resolve_num_trees(
    per_forest_vector(control[["num_trees"]], labels, "num_trees", 50L, joint),
    response[["n_forest"]], response[["n_aux"]])

  # The leaf scale divides by the square root of that forest's *own* tree count,
  # so a forest with fewer trees gets a proportionally larger prior per leaf and
  # the prior on the sum is unchanged. `k` is the usual way to say it and is
  # per-forest for the same reason `sigma_mu` is.
  k <- rep(per_forest_vector(control[["k"]], labels, "k", 2, joint),
           length.out = response[["n_forest"]])

  engine_control[["sigma_mu"]] <-
    per_forest_vector(control[["sigma_mu"]], labels, "sigma_mu", NULL, joint) %or%
    (3 * response[["eta_scale"]] / (k * sqrt(engine_control[["num_trees"]])))

  if (length(engine_control[["sigma_mu"]]) != response[["n_forest"]]) {
    engine_control[["sigma_mu"]] <- rep(engine_control[["sigma_mu"]],
                                        length.out = response[["n_forest"]])
  }

  # Progress, if the caller has asked for any. Sized for the whole run across
  # every chain, so one bar fills once rather than one per chain restarting the
  # count. Built here, in the calling session, because that is what lets it work
  # under `future.apply`: the reporter is captured by the engine closure, sent to
  # each worker, and the conditions it signals are relayed back.
  progress <- progress_reporter(chains, control)
  engine_control[["progress"]] <- progress[["report"]]
  engine_control[["progress_ticks"]] <- progress[["ticks"]]

  # Everything above is a deterministic function of the data and is done once;
  # only the sampler itself is repeated per chain.
  # A prior-only fit hands the engine a zero weight for every observation. The
  # weight multiplies that observation's log density, its gradient and its
  # curvature in one place (`src/family.h`), so zeroing it flattens the
  # likelihood exactly: every tree move is then accepted on the prior alone and
  # every leaf is drawn from its prior. The sampler is not otherwise touched.
  engine_weights <- if (prior_only) {
    prior_only_check(response[["family"]])

    rep.int(0, length(response[["weights"]]))
  }
  else {
    response[["weights"]]
  }

  engine <- function(ignored) {
    .bartisan_fit(X = unit$x,
                  has_na = as.integer(has_na),
                  y = response[["y"]],
                  weights = engine_weights,
                  offset = response[["offset"]],
                  group_probs = group_probs,
                  family_name = response[["family"]],
                  link = response[["link"]],
                  family_opts = response[["opts"]],
                  control = engine_control,
                  random_spec = random_spec(random),
                  codes = levels_info[["codes"]],
                  cat_col = levels_info[["cat_col"]],
                  n_levels = levels_info[["n_levels"]],
                  vc_basis = vc[["basis"]] %or% matrix(0, 0L, 0L))
  }

  draws <- {
    if (chains == 1L) engine(1L)
    else combine_chains(run_chains(engine, chains))
  }

  # What the fit reports is whether a rewriting happened, not which families one
  # was permitted for. `control[["augment"]]` is a request on the way in -- a
  # flag, or the names of the families it may apply to -- and the answer on the
  # way out is a single yes or no about this fit. The request itself is still
  # recoverable from `attr(control, "supplied")` and from the call.
  control[["augment"]] <- isTRUE(draws[["augmented"]])

  out <- list(call = cl,
              formula = formula,
              terms = mt,
              family = family,
              control = control,
              n = n,
              chains = chains,
              num_forest = draws[["num_forest"]],
              num_trees = draws[["num_trees"]],
              soft = control[["soft"]],
              gate = control[["gate"]],
              # NULL rather than an empty list when the formula has no bars, so
              # that everything downstream can test for it with one idiom and
              # `print()` does not announce a random part that is not there.
              random = random %or% NULL,
              ranef = draws[["ranef"]],
              tau = draws[["tau"]],
              eta = draws[["eta"]],
              counts = draws[["counts"]],
              sigma_mu = draws[["sigma_mu"]],
              # Every prior setting as the engine received it, which is not
              # what the caller wrote: the per-forest arguments have been spread
              # to one value each, `k` and `sigma_mu` are two ways of saying the
              # same thing and the default is neither, and an intercept-only
              # forest has had its branching probability zeroed. Kept so that
              # `prior_summary()` reports the prior that was used rather than
              # rebuilding it, and small enough not to matter: a handful of
              # vectors as long as there are forests.
              prior = prior_record(engine_control, control, k,
                                   response[["eta_scale"]], split,
                                   ncol(group_probs),
                                   draws[["num_forest"]]),
              bandwidth = draws[["bandwidth"]],
              loglik = as.vector(draws[["loglik"]]),
              forest_flat = draws[["forest_flat"]],
              tree_start = draws[["tree_start"]],
              # The Dirichlet process mixture, when there is one: a flat vector
              # of (mean, standard deviation, weight) triples with one offset per
              # draw, because its component count changes from draw to draw.
              mixture_flat = draws[["mixture_flat"]],
              mixture_start = draws[["mixture_start"]],
              intercept = response[["offset"]][, 1L],
              has_offset = !is_null(model_offset),
              # Kept so that the conditional density of the training data can be
              # evaluated without asking the caller to hand the outcome back.
              y = response[["y"]],
              prior_only = prior_only,
              prior_weights = response[["weights"]],
              family_opts = response[["opts"]],
              levels = response[["levels"]],
              num_cat = response[["num_cat"]],
              # The model frame is kept because the packages that build
              # counterfactual grids -- marginaleffects through insight -- need
              # the data the model saw, not just its terms. `glm()` keeps it for
              # the same reason and by the same default.
              model = mf,
              xlevels = stats::.getXlevels(mt, mf),
              contrasts = design$contrasts,
              assign = design$assign,
              term_labels = design$term_labels,
              group_names = colnames(group_probs),
              unit_maps = unit$maps,
              level_codes = levels_info,
              vc = vc,
              has_na = has_na,
              x_transform = control[["x_transform"]])

  if (!is_null(draws[["aux"]])) {
    out[["aux"]] <- draws[["aux"]] |>
      setColnames(draws[["aux_names"]])
  }

  names(out[["eta"]]) <- predictor_names(out)

  # The random part is indexed by additive predictor too, so it takes the same
  # names -- which is what makes the diagnostics table readable when a family has
  # more than one.
  if (!is_null(out[["ranef"]])) {
    names(out[["ranef"]]) <- predictor_names(out)
    names(out[["tau"]]) <- predictor_names(out)

    labels <- names(random)
    levels_per <- pluck(random, "num_levels", integer(1L))

    for (h in seq_along(out[["ranef"]])) {
      colnames(out[["ranef"]][[h]]) <- unlist(lapply(random, function(z) {
        sprintf("%s:%s", z[["label"]], z[["levels"]])
      }), use.names = FALSE)
      colnames(out[["tau"]][[h]]) <- labels
    }
  }
  names(out[["counts"]]) <- predictor_names(out)

  for (h in seq_along(out[["counts"]])) {
    colnames(out[["counts"]][[h]]) <- colnames(group_probs)
  }
  colnames(out[["sigma_mu"]]) <- predictor_names(out)

  class(out) <- "bartisan_fit"

  out[["fitted"]] <- fitted_from_eta(out, vc_combine(out, out[["eta"]], NULL),
                                     average = TRUE)

  # Trimmed to the reported forests. The target carries one value per forest the
  # engine builds, which includes the depth-zero forests standing in for a custom
  # family's nuisance parameters; `out$sigma_mu` records only the forests that are
  # additive predictors. Passing the untrimmed target compared column 1 against
  # the target of whichever forest happened to line up under recycling, and with
  # two predictors and one nuisance parameter it also warned about the length.
  warn_runaway_scale(out, engine_control[["sigma_mu"]][seq_len(out[["num_forest"]])])

  out
}

# An unoptimized build of the compiled code runs five to twenty times slower and
# is otherwise indistinguishable, which makes it very easy to draw conclusions
# about the sampler's speed from the wrong numbers. `devtools::load_all()` and
# `devtools::install()` compile without optimization and leave the object files
# behind for a later `R CMD INSTALL` to reuse, so this is not a rare accident.
# Warned once per session, since it is a property of the installation.
warn_unoptimized <- function() {
  if (isTRUE(the$checked_optimized) || .bartisan_optimized()) {
    the$checked_optimized <- TRUE
    return(invisible(NULL))
  }

  the$checked_optimized <- TRUE

  arg::wrn(c(
    "{.pkg bartisan}'s compiled code was built without optimization, which makes
     it 5 to 20 times slower than it should be.",
    i = "Reinstall from a clean source directory:
         {.code pkgbuild::clean_dll(); R CMD INSTALL --preclean .}",
    i = "{.code bartisan:::.bartisan_optimized()} reports the state of the
         installed library."))

  invisible(NULL)
}

# Run the sampler `chains` times with independent random number streams.
#
# The parallel axis that fits this sampler is the chain: a single chain is
# sequential by construction, since each sweep conditions on the last, and the
# per-move work is too small for the synchronization a within-chain split would
# need. Chains are embarrassingly parallel and are also what makes a convergence
# diagnostic possible at all.
#
# The backend is whatever the caller has planned through the future framework,
# so `plan(multisession)`, `plan(multicore)`, a cluster, or mirai's
# `plan(mirai_multisession)` all work without this package choosing for them.
# `future.seed = TRUE` gives each chain its own L'Ecuyer stream, which is what
# makes the result reproducible from a single `set.seed()` regardless of how many
# workers happen to run it.
# Without future.apply the chains run one after another rather than refusing to
# run: parallelism is how fast the chains are, not whether the model is fitted,
# and several chains run sequentially is still what makes the convergence
# diagnostics available. The streams are drawn the same way in both branches, so
# the draws do not depend on which one ran.
run_chains <- function(engine, chains) {
  # Generated here rather than left to `future.seed = TRUE`, so that both
  # branches below draw from the same streams. Otherwise the same script would
  # give different draws depending on whether future.apply happened to be
  # installed, which is a worse failure than being slow.
  seeds <- parallel_streams(chains)

  if (use_future()) {
    return(future.apply::future_lapply(seq_len(chains), engine,
                                       future.seed = seeds,
                                       future.packages = "bartisan"))
  }

  restore <- restore_stream()
  on.exit(restore(), add = TRUE)

  lapply(seq_len(chains), function(i) {
    assign(".Random.seed", seeds[[i]], envir = globalenv())
    engine(i)
  })
}

# Which families a flattened likelihood actually leaves at the prior.
#
# Zeroing the weights switches off everything that reaches the likelihood
# through `Family::logdens()`, which is the additive predictor and the leaves.
# It does not switch off an update written against the response directly, and
# the weight is a convenience rather than the only way to flatten one: an
# update that reads the data can be told to read zeros instead. Two were, in
# `src/family.cpp`, by having them ask whether the observation carries any
# weight at all:
#
#   `dpm()` and `dpm_aft()` assign each observation to a mixture atom from the
#   residual `y - eta` and redraw the atoms from those residuals. At zero weight
#   the residual now drops out of both, so the label comes from the Chinese
#   restaurant prior and the atom from the base measure, which is the
#   conditional an observation contributing no likelihood leaves behind.
#
#   `multinomial(link = "probit")` drew each latent utility with variance
#   `1 / (w * prec)`, infinite at zero weight, and the covariance drawn from
#   those utilities was no longer symmetric. A weightless observation's
#   utilities now sit on the predictor, leaving the covariance at its prior.
#
# What is left are the cases where the weight is not the obstacle. Two families
# draw a cutpoint from a target that is the weighted log probability *with no
# prior term at all*, so at zero weight the target is not flattened but empty
# and the slice sampler walks a flat improper density out to whatever bound the
# code happens to put on it: 1e4 for `ordinal()`, and plus or minus 30 for
# `ordbeta()`. There is nothing to fall back on until the model specifies a
# prior over ordered cutpoints, which is a modeling decision and not a switch.
#
# Every other auxiliary parameter drawn by slice sampling carries its prior into
# the target, which was checked one at a time rather than assumed: the gamma
# prior on `negbin()`'s size, `Gamma()`'s shape, `Beta()`'s and `tweedie()`'s
# precision and `ph()`'s baseline hazard, the half-Cauchy on the AFT scale, and
# the uniform on `tweedie()`'s power. `multinomial(link = "probit")` draws its
# covariance from an inverse Wishart whose scatter is weighted, so at zero
# weight it falls back on the identity and the prior degrees of freedom.
#
# Refused rather than warned about, because the failure is silent: a fit comes
# back and the replicates look like replicates, all of them piled at one end.
# `ordinal()` and `ordbeta()` were both here. Their cutpoints now carry the
# induced-Dirichlet prior, so a flat likelihood leaves a proper density to draw
# from and every family supports `prior_only = TRUE`.
PRIOR_ONLY_REFUSED <- character()

# The prior as the engine received it, trimmed to the forests that are reported.
#
# `engine_control` holds the per-forest spread of everything in
# `PER_FOREST_DEFAULTS`, so these are one value per forest even where the caller
# wrote a scalar. `alpha_scale` is the one entry whose default the engine
# resolves rather than R: zero there means "the number of predictors this forest
# may split on", which is what `Hypers()` substitutes, so it is substituted here
# too and the recorded value is the one in force.
prior_record <- function(engine_control, control, k, eta_scale, split,
                         n_group, n_forest) {
  keep <- seq_len(n_forest)

  per_forest <- c(names(PER_FOREST_DEFAULTS), "sigma_mu")
  out <- lapply(engine_control[per_forest], function(x) {
    if (is_null(x)) NULL else x[keep]
  })

  mask <- split[["mask"]]

  out[["candidates"]] <- {
    if (is_null(mask)) rep.int(n_group, n_forest)
    else colSums(mask[, keep, drop = FALSE] != 0)
  }

  out[["alpha_scale"]] <- ifelse(out[["alpha_scale"]] > 0,
                                 out[["alpha_scale"]],
                                 out[["candidates"]])

  out[["k"]] <- rep(k, length.out = n_forest)[keep]
  out[["eta_scale"]] <- rep(eta_scale, length.out = n_forest)[keep]
  out[["sparsity"]] <- control[["sparsity"]]
  out[["split_prior"]] <- !is_null(split[["prior"]])
  out[["share_sparsity"]] <- isTRUE(control[["share_sparsity"]])
  out[["update_tau"]] <- isTRUE(control[["update_tau"]])

  out
}

prior_only_check <- function(family) {
  if (!family %in% names(PRIOR_ONLY_REFUSED)) {
    return(invisible(TRUE))
  }

  arg::err(c("{.code prior_only = TRUE} is not available for
              {.code {family}()}, because {PRIOR_ONLY_REFUSED[[family]]}.",
             i = "The cutpoints would walk out to the bound and the replicates
                  would pile at one end of the scale, with nothing in them to
                  say so.",
             i = "Every other family supports it. For an ordered outcome with a
                  modest number of categories, {.fn multinomial} is the nearest
                  thing that does."))
}

# One L'Ecuyer stream per chain, advanced from the current seed, which is what
# `future.seed = TRUE` does. Taking them from the session's own state is what
# makes a single `set.seed()` before the call reproduce the whole run.
parallel_streams <- function(chains) {
  old <- restore_stream()
  on.exit(old(), add = TRUE)

  RNGkind("L'Ecuyer-CMRG")
  seed <- get(".Random.seed", envir = globalenv())

  out <- vector("list", chains)

  for (i in seq_len(chains)) {
    out[[i]] <- seed
    seed <- parallel::nextRNGStream(seed)
  }

  out
}

# Captures the session's RNG state and returns the function that puts it back,
# including the kind, so that switching to L'Ecuyer for the streams does not
# leave the session on a generator it did not choose.
restore_stream <- function() {
  kind <- RNGkind()

  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    seed <- get(".Random.seed", envir = globalenv())
    function() {
      RNGkind(kind[1L], kind[2L], kind[3L])
      assign(".Random.seed", seed, envir = globalenv())
    }
  }
  else {
    function() {
      RNGkind(kind[1L], kind[2L], kind[3L])
    }
  }
}

# Stack the chains into one set of draws, in chain order. The stored forests are
# indexed by a flat position that runs iteration, then forest, then tree, so the
# record offsets of each chain after the first have to be shifted by the total
# length of the ones before it.
combine_chains <- function(fits) {
  first <- fits[[1L]]
  out <- first

  stack <- function(fits, name) {
    pluck(fits, name) |> do_rbind()
  }

  stack_list <- function(fits, name) {
    lapply(seq_along(first[[name]]), function(h) {
      lapply(fits, function(z) z[[name]][[h]]) |>
        do_rbind()
    })
  }

  out[["eta"]] <- stack_list(fits, "eta")
  out[["counts"]] <- stack_list(fits, "counts")
  out[["sigma_mu"]] <- stack(fits, "sigma_mu")
  out[["bandwidth"]] <- stack(fits, "bandwidth")
  out[["loglik"]] <- stack(fits, "loglik")

  if (!is_null(first[["aux"]])) {
    out[["aux"]] <- stack(fits, "aux")
  }

  if (!is_null(first[["ranef"]])) {
    out[["ranef"]] <- stack_list(fits, "ranef")
    out[["tau"]] <- stack_list(fits, "tau")
  }

  # A flat vector with per-draw offsets is concatenated by shifting every
  # chain's offsets past what came before it. Each chain's own offsets start
  # with a zero, which belongs only to the first.

  flat <- pluck(fits, "forest_flat")
  at <- utils::head(cumsum(c(0L, lengths(flat))), -1L)
  starts <- lapply(seq_along(fits), function(k) {
    fits[[k]][["tree_start"]][-1L] + at[k]
  })

  out[["forest_flat"]] <- unlist(flat, use.names = FALSE)
  out[["tree_start"]] <- c(0L, unlist(starts, use.names = FALSE))


  if (!is_null(first[["mixture_flat"]])) {
    flat <- pluck(fits, "mixture_flat")
    at <- utils::head(cumsum(c(0L, lengths(flat))), -1L)
    starts <- lapply(seq_along(fits), function(k) {
      fits[[k]][["mixture_start"]][-1L] + at[k]
    })

    out[["mixture_flat"]] <- unlist(flat, use.names = FALSE)
    out[["mixture_start"]] <- c(0L, unlist(starts, use.names = FALSE))
  }

  out
}

# Reduce over observations, returning NA rather than an infinity when every one
# of them is NA -- which happens when the quantity does not vary.
worst <- function(x, f) {
  if (!any(is.finite(x))) {
    return(NA_real_)
  }

  f(x, na.rm = TRUE)
}

# Split-R-hat (Gelman and Rubin, as revised in Gelman et al., BDA3): each chain
# is halved so that drift within a chain shows up as disagreement between the
# halves. `x` is draws by chains.
split_rhat <- function(x) {
  y <- split_chains(x)

  if (is_null(y)) {
    return(NA_real_)
  }

  half <- nrow(y)

  # A quantity the sampler holds fixed has nothing to diagnose, and the
  # arithmetic below cannot be trusted to say so. `var()` returns a clean zero
  # for a constant column; subtracting a column mean does not, because the mean
  # of many copies of a value need not be that value back -- exact for 0 and for
  # 2.5, not for `qnorm()` of an average rank, which is what this is handed. The
  # rounding error then passes the guard below and turns nothing into an R-hat of
  # 1. Asked here instead, the same way `ess_from()` asks it.
  if (isTRUE(diff(range(y)) == 0)) {
    return(NA_real_)
  }

  # The within-chain variances directly rather than through `apply()`, which
  # splits the matrix into a list and calls a closure per column. This runs four
  # times per column of draws (twice here, twice more on the late half), and the
  # arithmetic is the same two passes either way: 4.3x faster, agreeing to the
  # last bit of a variance.
  within <- mean(colSums((y - rep(colMeans(y), each = half))^2) / (half - 1))

  # `isTRUE()` rather than a `<=` comparison because the quantity being guarded
  # can be NaN as well as zero -- a chain of one draw, or a constant -- and
  # `NaN <= 0` is NA, which is not something `if` can act on.
  if (!isTRUE(within > 0)) {
    return(NA_real_)
  }

  between <- half * stats::var(colMeans(y))

  sqrt(((half - 1) / half * within + between / half) / within)
}

# Halve every chain, so that a chain that has drifted disagrees with itself.
split_chains <- function(x) {
  draws <- nrow(x)
  half <- draws %/% 2L

  if (half < 2L || ncol(x) < 2L) {
    return(NULL)
  }

  cbind(x[seq_len(half), , drop = FALSE],
        x[draws - half + seq_len(half), , drop = FALSE])
}

# Rank-normalization: replace the draws by the normal scores of their pooled
# ranks. The point is that R-hat and the effective sample size are derived for
# quantities with finite variance and behave badly without it, and a rank
# transform guarantees it whatever the posterior looks like -- which also makes
# the diagnostic invariant to any monotone reparameterization. Blom's offsets.
rank_normalize <- function(x) {
  n <- length(x)
  o <- order(x)
  sorted <- x[o]

  # With no ties, averaging them is a no-op and the ranks are just the inverse
  # of the ordering, which on its own is 1.67x faster than asking `rank()` for
  # average ties, and worth 3% to 5% of the convergence pass measured end to
  # end. Draws of a continuous quantity have no ties, and those are every column
  # the pass walks bar a handful of scalar rows, so the fast path is nearly all
  # of the work. Where there are ties, or a missing value, `rank()` answers as
  # before; that costs the ordering twice, about 20% on those rows.
  r <- {
    if (anyNA(x) || any(sorted[-1L] == sorted[-n], na.rm = TRUE)) {
      rank(x, ties.method = "average")
    }
    else {
      out <- numeric(n)
      out[o] <- seq_len(n)
      out
    }
  }

  stats::qnorm((r - 3 / 8) / (n - 1 / 4)) |>
    matrix(nrow = nrow(x), ncol = ncol(x))
}

# Rank-normalized, folded, split R-hat (Vehtari, Gelman, Simpson, Carpenter and
# Buerkner 2021). Two diagnostics, maximized: the rank-normalized one catches
# chains that disagree about the middle of the distribution, and the folded one
# -- the same computation applied to the distance from the median -- catches
# chains that agree about the middle and disagree about the spread, which the
# first is blind to.
rhat_rank <- function(x, normalized = NULL) {
  if (is_null(split_chains(x)) || !all(is.finite(x))) {
    return(NA_real_)
  }

  # `normalized` is `rank_normalize(x)` when a caller already has it;
  # `diagnosis_stats()` does, because the bulk effective sample size starts from
  # the same thing.
  bulk <- (normalized %or% rank_normalize(x)) |>
    split_rhat()

  folded <- abs(x - stats::median(x)) |>
    rank_normalize() |>
    split_rhat()

  # A quantity the sampler holds fixed -- an ordinal model's first cutpoint, say
  # -- has no between-chain variance to compare, so both are NA. Reducing that
  # with `na.rm` returns -Inf and warns; there is simply nothing to diagnose.
  if (is.na(bulk) && is.na(folded)) {
    return(NA_real_)
  }

  max(bulk, folded, na.rm = TRUE)
}

# Effective sample size, following the algorithm of Vehtari et al. (2021) as
# implemented in Stan. `y` is draws by chains, already split.
#
# The autocorrelations are pooled across chains in a way that borrows the
# between-chain variance: a chain sitting somewhere the others are not looks
# well mixed on its own, and dividing by the pooled variance rather than its own
# is what penalizes it.
# Biased autocovariance at every lag, one column per chain, by the Wiener-
# Khinchin route: the transform of the padded series times its own conjugate is
# the transform of its autocovariance. The padding is what makes the implied
# convolution linear rather than circular, so the result matches
# `acf(type = "covariance", demean = TRUE)` to the last bit rather than
# approximately.
#
# This used to be one `stats::acf()` call per chain, which is the same arithmetic
# but pays for building an `acf` object each time -- dimnames, an `outer()`, a
# `deparse1()` of the series name -- and that bookkeeping, not the arithmetic,
# was most of the cost of a diagnostics pass. Called once per observation per
# additive predictor, it added up: the swap is worth about nine times on the
# autocovariance and about four on the pass as a whole.
autocovariance <- function(y) {
  draws <- nrow(y)

  # Zero-padded to at least twice the length, at a power of two so the transform
  # takes its fast path.
  nfft <- as.integer(2^ceiling(log2(2 * draws)))

  centered <- y - rep(colMeans(y), each = draws)
  padded <- rbind(centered, matrix(0, nfft - draws, ncol(y)))

  transform <- stats::mvfft(padded)

  acov <- Re(stats::mvfft(transform * Conj(transform), inverse = TRUE)) / nfft

  acov[seq_len(draws), , drop = FALSE] / draws
}

ess_from_split <- function(y) {
  draws <- nrow(y)
  chains <- ncol(y)

  if (draws < 4L) {
    return(NA_real_)
  }

  pooled <- rowMeans(autocovariance(y))
  mean_var <- pooled[1L] * draws / (draws - 1)
  var_plus <- mean_var * (draws - 1) / draws

  if (chains > 1L) {
    var_plus <- var_plus + stats::var(colMeans(y))
  }

  # Guarded with `isTRUE()` for the reason given in `split_rhat()`.
  if (!isTRUE(var_plus > 0) || !isTRUE(mean_var > 0)) {
    return(NA_real_)
  }

  # Every autocorrelation at once, so the sequence below indexes a vector rather
  # than calling a closure per lag. `rho[k + 1L]` is the correlation at lag `k`.
  rho <- 1 - (mean_var - pooled) / var_plus

  # Geyer's initial positive sequence: walk the autocorrelations in adjacent
  # pairs and stop at the first pair whose sum goes negative, which is where the
  # estimates stop being informative.
  #
  # Preallocated rather than grown with `c()`. On a slowly mixing quantity the
  # sequence runs to the last available lag rather than stopping early -- the
  # fitted values of a forest do, at a median of 98 lags out of 98 on a
  # four-chain fit -- so growing the vector reallocated it on every pair, and
  # that copying was most of what a diagnostics pass cost once the
  # autocovariance itself was cheap. The zeros a negative pair contributes are
  # already there.
  kept <- numeric(draws + 2L)
  kept[1L] <- 1
  kept[2L] <- rho[2L]
  n_kept <- 2L
  t <- 1L

  while (t < draws - 4L && kept[n_kept - 1L] + kept[n_kept] > 0) {
    even <- rho[t + 2L]
    odd <- rho[t + 3L]

    if (even + odd >= 0) {
      kept[n_kept + 1L] <- even
      kept[n_kept + 2L] <- odd
    }

    n_kept <- n_kept + 2L
    t <- t + 2L
  }

  kept <- kept[seq_len(n_kept)]

  extra <- max(kept[n_kept - 1L], 0)

  # Force the paired sums to be non-increasing, which is what makes the
  # estimator conservative rather than merely unbiased. With too few kept lags
  # for a second pair there is nothing to compare against.
  last <- length(kept) - 3L
  pairs <- if (last >= 2L) seq(2L, last, by = 2L) else integer()

  for (k in pairs) {
    previous <- kept[k - 1L] + kept[k]

    if (kept[k + 1L] + kept[k + 2L] > previous) {
      kept[k + 1L] <- previous / 2
      kept[k + 2L] <- kept[k + 1L]
    }
  }

  tau <- max(-1 + 2 * sum(kept) + extra,
             1 / log10(draws * chains))

  draws * chains / tau
}

ess_from <- function(x) {
  y <- split_chains(x)

  if (is_null(y) || !all(is.finite(y))) {
    return(NA_real_)
  }

  # A quantity the sampler holds fixed has no autocorrelation to estimate. The
  # variance guard inside ess_from_split() does not catch it: the sample
  # autocovariance of a constant is a rounding error rather than exactly zero, so
  # it passes the guard and the ratios built on it are meaningless. Comparing the
  # values themselves is exact.
  if (length(unique(as.vector(y))) < 2L) {
    return(NA_real_)
  }

  ess_from_split(y)
}

ess_bulk <- function(x) {
  rank_normalize(x) |>
    ess_from()
}

# Tail effective sample size: the smaller of the two effective sample sizes for
# the indicator that a draw falls below the 5% and above the 95% quantile. It is
# reported separately because a chain can be perfectly adequate for a posterior
# mean and nowhere near adequate for an interval endpoint -- the mean is an
# average over every draw, and a tail quantile depends on the few draws out
# there.
ess_tail <- function(x) {
  q <- stats::quantile(x, c(0.05, 0.95), names = FALSE, na.rm = TRUE)

  # The indicator takes two values, so rank-normalizing it is an affine map, and
  # an effective sample size is built from ratios of autocovariances and so is
  # invariant to one. It used to be rank-normalized anyway, which cost two of
  # the seven rankings a column needs and changed nothing: measured over
  # randomized cases the two agree to 8e-16 or exactly. Dropping it is also what
  # the quantity is defined as, the effective sample size of the indicator
  # itself.
  worst(c(ess_from((x <= q[1L]) * 1),
          ess_from((x >= q[2L]) * 1)), min)
}

# The leaf scale is drawn under a half-Cauchy prior, which has no upper bound.
# Where the response is perfectly, or nearly, separated by the predictors the
# likelihood rewards an unbounded predictor and that prior is not enough to hold
# the scale down: it wanders instead of settling, and the additive predictor can
# reach values for which the fitted probabilities are numerically zero or one.
# The condition is worth naming, because the remedy is a setting the caller
# already has.
warn_runaway_scale <- function(object, target) {
  if (!isTRUE(object[["control"]][["update_sigma_mu"]])) {
    return(invisible(NULL))
  }

  ratio <- colMeans(object[["sigma_mu"]]) / target
  at <- which(ratio > 5)

  if (!is_null(at)) {
    arg::wrn(c(
      "The leaf scale settled {round(max(ratio[at]))} times above its prior
     median, which usually means the response is close to separable by the
     predictors.",
      i = "The additive predictor is then only weakly identified; fix the scale
         with {.code update_sigma_mu = FALSE} in {.fn bartisan_control} if the
         draws look unstable."))
  }
}

# Names for the additive predictors, which are the columns of most outputs.
predictor_names <- function(object) {
  forest_labels(object[["family"]][["family"]], object[["family_opts"]],
                object[["levels"]], object[["num_forest"]], object[["vc"]])
}

# Names of the additive predictors, in the order the engine builds them. The
# first is always the main parameter -- the one a single-forest family would
# have on its own -- and the rest follow in the order documented on
# [bartisan-families]. These are the names that label the columns of most
# outputs, and the names a per-forest argument may be keyed by.
forest_labels <- function(family, opts, levels, n_report, vc = NULL) {
  slopes <- vc[["slopes"]] %or% 0L

  # A varying-coefficient model's forests are, for each of the family's additive
  # predictors, its control function and one per coefficient. The parameter's own
  # name is dropped when the family has only one, so the common case is
  # `(Intercept)` and the covariate's name rather than `eta` and `eta:z`, which
  # is what per-forest arguments have to be typed as. With several parameters the
  # name has to stay, since `mean:z` and `log_sd:z` are different forests.
  if (slopes > 0L) {
    base <- forest_labels(family, opts, levels, n_report - slopes)
    drop <- length(base) == 1L && identical(base, "eta")

    return(unlist(lapply(seq_along(base), function(h) {
      keep <- which(pluck(vc[["specs"]], "param", integer(1L)) == h)
      vc_forest_labels(base[h], vc[["specs"]][keep], vc[["parts"]][keep], drop)
    }), use.names = FALSE))
  }

  if (identical(family, "multinomial")) {
    if (isTRUE(opts[["symmetric"]])) {
      return(levels)
    }

    return(levels[-1L])
  }

  # One latent variable per non-reference category, named for the contrast it
  # carries.
  if (identical(family, "mnp")) {
    return(sprintf("%s-%s", levels[-1L], levels[1L]))
  }

  if (identical(family, "gaussian_ls")) {
    return(c("mean", "log_sd"))
  }

  if (identical(family, "Gamma_ls")) {
    return(c("mean", "log_dispersion"))
  }

  if (family %in% c("zip", "zinb")) {
    return(c("count", "zero"))
  }

  if (identical(family, "custom") && n_report > 1L) {
    return(sprintf("eta%d", seq_len(n_report)))
  }

  "eta"
}

# Whether the family's forests are components of the response distribution or
# parts of one vector-valued parameter. The multinomial families are the second
# kind: their forests are the levels of one categorical parameter and act
# together rather than describing separate pieces of the distribution, so
# splitting a setting across them says nothing a caller would mean. Every
# per-forest argument therefore applies to all of their forests at once.
joint_forests <- function(family) {
  family %in% c("multinomial", "mnp")
}

# Rows the model cannot use. A missing predictor is handled by the splitting
# rules; a missing response, prior weight or offset is not something the model
# can work around, so those rows go. This only ever removes anything when the
# caller asked for missing values to be kept, since na.omit has already removed
# them otherwise.
drop_unusable_rows <- function(mf, mt) {
  columns <- c(attr(mt, "response"),
               match(c("(weights)", "(offset)"), names(mf), 0L))
  columns <- columns[columns > 0L]

  if (is_null(columns)) {
    return(mf)
  }

  keep <- stats::complete.cases(mf[columns])

  if (all(keep)) {
    return(mf)
  }

  arg::wrn("dropping {sum(!keep)} row{?s} with a missing response, weight
            or offset; rows missing only predictors are kept")

  mf[keep, , drop = FALSE]
}

# Design matrix with an indicator per factor level rather than contrast coding.
# A tree splits on a single column, so "is level j" should be available as a
# rule for every level; with contrast coding the reference level is only
# reachable as the conjunction of all the others.
build_design <- function(mt, mf) {
  predictors <- attr(mt, "term.labels")
  variables <- get_varnames(stats::delete.response(mt))

  categorical <- variables[vapply(variables, function(nm) {
    if (!nm %in% names(mf)) {
      return(FALSE)
    }
    z <- mf[[nm]]
    is.factor(z) || is.character(z)
  }, logical(1L))]

  contrasts <- NULL

  if (!is_null(categorical)) {
    contrasts <- lapply(categorical, function(nm) {
      stats::contrasts(as.factor(mf[[nm]]), contrasts = FALSE)
    }) |>
      setNames(categorical)
  }

  x <- stats::model.matrix(mt, mf, contrasts.arg = contrasts)
  assign <- attr(x, "assign")

  keep <- assign != 0L
  x <- x[, keep, drop = FALSE]
  assign <- assign[keep]

  if (ncol(x) == 0L) {
    arg::err("the model has no predictor columns")
  }

  # Columns that never vary cannot support a split. A column that is constant
  # where it is observed still varies in whether it is observed at all, which is
  # something a rule can split on, so it stays.
  varies <- apply(x, 2L, function(z) {
    length(unique(z[!is.na(z)])) > 1L || (anyNA(z) && !all(is.na(z)))
  })

  if (!any(varies)) {
    arg::err("no predictor varies across observations")
  }

  if (!all(varies)) {
    arg::wrn("dropping {sum(!varies)} constant predictor column{?s}:
              {.val {colnames(x)[!varies]}}")
    x <- x[, varies, drop = FALSE]
    assign <- assign[varies]
  }

  list(x = x, assign = assign, term_labels = predictors,
       contrasts = contrasts)
}

unit_transform <- function(x, type) {
  maps <- lapply(seq_len(ncol(x)), function(j) make_unit_map(x[, j], type)) |>
    setNames(colnames(x))

  out <- x

  for (j in seq_len(ncol(x))) {
    out[, j] <- maps[[j]](x[, j])
  }

  list(x = out, maps = maps)
}

apply_unit_maps <- function(x, maps) {
  out <- x

  for (j in seq_len(ncol(x))) {
    out[, j] <- maps[[j]](x[, j])
  }

  out
}

# Settings passed to `bartisan()` through `...` are control arguments, as in
# `glm()`. They are merged with whatever the caller supplied to
# `bartisan_control()` and the whole thing is re-validated, so a bad value passed
# this way fails the same way it would have failed there.
merge_control <- function(control, dots) {
  if (is_null(dots)) {
    return(control)
  }

  allowed <- names(formals(bartisan_control))
  nms <- names(dots)

  if (is_null(nms) || !all(nzchar(nms))) {
    arg::err("arguments passed to {.fn bartisan} in {.arg ...} must be named")
  }

  if (anyDuplicated(nms) > 0) {
    arg::err("arguments passed to {.fn bartisan} in {.arg ...} must have distinct names")
  }

  bad <- setdiff(nms, allowed)

  if (!is_null(bad)) {
    arg::err("{.val {bad}} {?is/are} not {?an argument/arguments} of {.fn bartisan_control}")
  }

  args <- attr(control, "supplied") %or% list()

  # Single-bracket assignment from a list, rather than `modifyList()`, so that an
  # explicit `NULL` (a legitimate value for `sigma_mu` and `alpha_scale`) is kept
  # rather than removing the entry.
  args[nms] <- dots

  do.call(bartisan_control, args)
}

# The engine takes the gate as an integer, matching GateShape in node.h. Kept in
# one place so the two orderings cannot drift apart. A hard rule has no gate
# shape to choose, so it takes the default code and the engine ignores it.
gate_code <- function(gate) {
  code <- match(gate, c("logistic", "smoothstep", "smootherstep")) - 1L

  if (is.na(code)) 1L else code
}

# One tree count per additive predictor. The default of 50 is set in
# `bartisan_control()`'s signature and does not depend on the rules, though the
# measurement suggested it might: on the Friedman function held-out error levels
# off by 20 trees with soft rules and keeps improving to 50 with hard ones. 50 is
# kept for both because a smaller forest mixes worse -- on `lalonde`, four chains
# at 20 soft trees disagreed by 35% on an average contrast against 9% at 50 --
# and the Friedman gain at 20 was 5%.
resolve_num_trees <- function(num_trees, n_forest, n_aux = 0L) {
  num_trees <- as.integer(num_trees)

  # `num_trees` is about the additive predictors, so the count the caller sees
  # excludes the trailing nuisance forests.
  n_report <- n_forest - n_aux

  if (length(num_trees) > n_report) {
    arg::err("{.arg num_trees} has {length(num_trees)} value{?s} but this
              family has {n_report} additive predictor{?s}")
  }

  # A nuisance parameter is one scalar, so its forest is one tree, and the engine
  # pins it so that the tree can never split.
  c(rep(num_trees, length.out = n_report), rep.int(1L, n_aux))
}
