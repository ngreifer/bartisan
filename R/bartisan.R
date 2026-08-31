#' Fit a generalized Bayesian additive regression trees model
#'
#' @description
#' Fits a BART model in which the response distribution is arbitrary rather than
#' restricted to the conditionally conjugate cases, using the
#' Laplace-approximation reversible-jump sampler of Linero (2025). Decision
#' rules may be soft, as in Linero and Yang (2018), which gives smoother fits
#' than the step functions of standard BART.
#'
#' The interface deliberately mirrors [stats::glm()]: a formula, a data frame
#' and a family. Ordinary [stats::family] objects work unchanged, including
#' their links, and the extra families that `glm()` has no counterpart for are
#' documented at [bartisan-families], along with `custom_family()` for a
#' likelihood of your own.
#'
#' @param formula a model formula. The right-hand side lists candidate
#'   predictors; the model finds interactions and nonlinearity on its own, so
#'   `y ~ .` is usually the right specification. Survival families take a
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
#'   bartisan(list(y ~ x1 + x2, ~ x2 + x3), data = d, family = location_scale())
#'   bartisan(list(mean = y ~ x1 + x2, log_sd = ~ x2), data = d,
#'            family = location_scale())
#'   ```
#'
#'   One formula applies to every forest, which is the ordinary case. A predictor
#'   left out of one forest's formula is still in the data and is never split on
#'   by that forest.
#' @param data a data frame containing the variables in `formula`.
#' @param family the response distribution, as a [stats::family] object, one of
#'   the families in [bartisan-families], or a name. The default, `NULL`, reads
#'   one off the response and says which it chose; see Details for the rules and
#'   for what is supported.
#' @param weights optional prior weights. For a binomial response given as
#'   proportions, these are the numbers of trials, as in `glm()`.
#' @param offset optional known component of the additive predictor, on the link
#'   scale.
#' @param subset optional vector specifying a subset of rows to use.
#' @param na.action how to handle missing values. The default,
#'   [stats::na.pass], keeps rows whose *predictors* are missing and lets the
#'   splitting rules decide where they go, which is what the trees are able to do
#'   and `lm()` and `glm()` are not; see Details. Pass [stats::na.omit] to drop
#'   any row with a missing value anywhere instead. Rows with a missing response,
#'   weight or offset are dropped either way, with a warning, since there is
#'   nothing to fit them to.
#' @param control a list of sampler and prior settings from
#'   [bartisan_control()].
#' @param ... further arguments to [bartisan_control()]. They are
#'   merged into `control`, overriding any value given there, so that
#'   `bartisan(..., num_trees = 20)` and
#'   `bartisan(..., control = bartisan_control(num_trees = 20))` are the same
#'   call. Names that are not arguments of [bartisan_control()] are an error
#'   rather than being silently ignored.
#'
#' @details
#' # What the sampler does
#'
#' Standard BART relies on the leaf parameters being integrable in closed form,
#' which restricts it to a Gaussian response, or to models that can be reduced
#' to one by data augmentation. Linero's algorithm removes that restriction. At
#' each candidate move it builds a Gaussian approximation to the conditional
#' posterior of the affected leaf parameters, by Fisher scoring, and uses that
#' approximation as the proposal in a reversible-jump Metropolis step. The
#' approximation only has to be good enough to be accepted often; the stationary
#' distribution is the exact posterior either way.
#'
#' What a new family therefore has to supply is only the log density of one
#' observation and its first two derivatives with respect to the additive
#' predictor. Families whose response has more than one unconstrained parameter,
#' such as `multinomial()` and `location_scale()`, carry one forest per
#' parameter. Because that is the whole interface, it can be reached from R:
#' [custom_family()] takes the log density as an R function and differences it
#' for the derivatives, and a link the package does not compile is composed onto
#' the scale its family works on the same way.
#'
#' # The family is inferred when you do not name one
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
#' A message reports the choice. Naming `family` yourself is what silences it,
#' which is the same thing you would do to change the choice.
#'
#' Two of these are worth saying out loud. A **count** is not inferred as
#' `poisson()`: a non-negative integer response is often Poisson and often not,
#' and the Poisson variance assumption is strong enough that making it silently
#' would be a modeling decision taken on the caller's behalf. Gaussian is the
#' weaker guess and the one whose failure is easy to see. And a numeric response
#' with exactly two values that are *not* zero and one -- `c(1, 2)`, say -- is
#' Gaussian rather than binomial, because which of the two counts as the success
#' is not something to guess at.
#'
#' # Soft decision rules
#'
#' By default a decision rule is a smooth gate rather than a step, so an
#' observation reaches every leaf with some weight and the fitted function is
#' smooth. `gate` in [bartisan_control()] chooses both whether the rules are soft
#' and, if they are, the gate's shape; the default is the bounded
#' `"smoothstep"`, and `"logistic"` is Linero and Yang's (2018) original. This
#' costs more per iteration, since a leaf now
#' touches every observation rather than the ones inside its cell, and it makes
#' the leaf parameters of a tree dependent on one another. Combining soft rules
#' with a non-conjugate likelihood is an extension of Linero (2025), which
#' leaves it as an open problem; it is handled here by giving the reversible-jump
#' move a bivariate Laplace proposal for the pair of child leaves, which reduces
#' to Linero's independent pair exactly when the rules are hard.
#'
#' Set `gate = "hard"` in [bartisan_control()] for the faster hard-rule sampler.
#'
#' # Random intercepts
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
#' separate set for each -- a zero-inflated count model has a group effect on the
#' count part and another on the inflation part -- and they are independent of one
#' another.
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
#' **When to reach for this rather than putting the group in as a predictor.** A
#' grouping factor can also go in the fixed part, where a tree splits on it like
#' anything else, and with few large groups that is the better choice --
#' measured, it beats a random intercept, because the group means are well
#' determined without pooling and a split can interact the group with the
#' covariates. The random intercept wins when there are many small groups, which
#' is where partial pooling earns its keep: at 250 groups of four observations it
#' cut held-out error by 30% against the factor route, and at five groups of a
#' hundred it lost to it.
#'
#' A level of `group` that was not present at fitting time is given the prior mean
#' of zero when predicting, with a warning.
#'
#' # Missing predictor values
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
#' This is *missingness incorporated in attributes* (Twala, Jones and Hand 2008;
#' for BART, Kapelner and Bleich 2015). The third rule is what lets the model
#' split on missingness itself, so a variable whose *absence* carries the signal
#' is usable even if its observed values say nothing. Since the choice is drawn
#' from its prior along with the variable and the cutpoint, it cancels from every
#' acceptance ratio, and a variable with no missing values is not given the extra
#' draw at all: complete data reproduces the sampler exactly as it was.
#'
#' A missing value takes a hard path through the tree even when the rules are
#' soft, which is the right thing -- there is nothing about being absent to
#' smooth over -- and it keeps the leaf weights summing to one.
#'
#' Two consequences to be clear about. `predict()` accepts missing values only in
#' columns that had them at fitting time, because only those columns' rules carry
#' an answer; elsewhere every rule would send the value the same arbitrary way,
#' so it is an error instead. And what this estimates is the mean of the response
#' given the predictors *and the pattern of missingness*. That is what you want
#' for prediction. If the estimand is a regression or causal effect defined on
#' complete data, multiple imputation is the right tool and this is not.
#'
#' # Preprocessing
#'
#' Predictors are mapped to the unit interval, because the cutpoint prior is
#' uniform on a node's live range and the soft-rule bandwidth is measured on the
#' predictor scale. Factors are expanded to an indicator per level and share a
#' single weight in the sparsity prior, so that a factor is selected or not as a
#' whole rather than one level at a time. The additive predictor starts from an
#' intercept-only fit, so the leaf prior describes departures from that fit
#' rather than the absolute level of the response. That starting value is the
#' exact null-model estimate for most families; for the accelerated failure time
#' families, where censoring makes the sample mean of the log times biased, and
#' for the zero-inflated and ordered beta families, it is a moment
#' approximation, which the sampler then moves away from.
#'
#' @returns
#' An object of class `bartisan`, a list with elements including:
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
#'   \item{`rhat`}{a data frame of convergence diagnostics, when more than one
#'     chain was run: rank-normalized folded split R-hat and the bulk and tail
#'     effective sample sizes (Vehtari et al. 2021) for the log likelihood, the
#'     leaf scales, the nuisance parameters and the additive predictor. R-hat
#'     above about 1.01 says the chains have not agreed; an effective sample size
#'     below about 400 says the run is too short for the quantity it belongs to,
#'     and the tail column is the one that governs interval endpoints.}
#'   \item{`sigma_mu`, `bandwidth`}{draws of the leaf standard deviation and,
#'     for soft rules, the per-tree gate bandwidths.}
#'   \item{`loglik`}{the log likelihood at each draw.}
#'
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
#' Albert, J. H., & Chib, S. (1993). Bayesian analysis of binary and
#' polychotomous response data. *Journal of the American Statistical
#' Association*, 88(422), 669--679. \doi{10.1080/01621459.1993.10476321}
#'
#' Polson, N. G., Scott, J. G., & Windle, J. (2013). Bayesian inference for
#' logistic models using Polya-Gamma latent variables. *Journal of the American
#' Statistical Association*, 108(504), 1339--1349.
#' \doi{10.1080/01621459.2013.829001}
#'
#' Kapelner, A., & Bleich, J. (2015). Prediction with missing data via Bayesian
#' additive regression trees. *Canadian Journal of Statistics*, 43(2), 224--239.
#' \doi{10.1002/cjs.11248}
#'
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Buerkner, P.-C. (2021).
#' Rank-normalization, folding, and localization: an improved \eqn{\widehat{R}}
#' for assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667--718.
#' \doi{10.1214/20-BA1221}
#'
#' Twala, B. E. T. H., Jones, M. C., & Hand, D. J. (2008). Good methods for
#' coping with missing data in decision trees. *Pattern Recognition Letters*,
#' 29(7), 950--956. \doi{10.1016/j.patrec.2008.01.010}
#'
#' @seealso
#' [predict.bartisan_fit()], [bartisan_control()], [bartisan-families], and
#' `vignette("families", package = "bartisan")` for a family-by-family guide.
#'
#' @examples
#' set.seed(1)
#'
#' n <- 200
#' d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
#' d$y <- rbinom(n, 1, plogis(3 * sin(pi * d$x1 * d$x2) - 1))
#'
#' fit <- bartisan(y ~ x1 + x2 + x3, data = d, family = binomial(),
#'                control = bartisan_control(num_trees = 10, num_burn = 50,
#'                                          num_draws = 50, verbose = FALSE))
#' fit
#'
#' head(predict(fit, type = "response"))
#'
#' @export
bartisan <- function(formula, data, family = NULL, weights = NULL,
                     offset = NULL, subset = NULL,
                     na.action = stats::na.pass,
                     control = bartisan_control(), ...) {

  cl <- match.call()

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

  mf <- match.call(expand.dots = FALSE)
  keep <- match(c("formula", "data", "subset", "weights", "offset"),
                names(mf), 0L)
  mf <- mf[c(1L, keep)]
  mf[["formula"]] <- reformulas::subbars(formula)
  mf[["drop.unused.levels"]] <- TRUE

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
  if (!missing(data) && is.data.frame(data) && identical(formula, split$fixed)) {
    resolved <- stats::update(stats::terms(mf[["formula"]], data = data), . ~ .)
    environment(resolved) <- environment(mf[["formula"]])
    mf[["formula"]] <- resolved
    split$fixed <- resolved

    # Stored as well as used. `insight::find_formula()` reads this field, and
    # from `death ~ . - days` it concluded that the model's one predictor was
    # `days`, which is the variable the formula removes. Everything built on
    # that -- `avg_comparisons()` most visibly -- then described the wrong model.
    formula <- resolved
  }

  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

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
    arg::err("{.arg formula} must include at least one predictor")
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
                   logical(1L))
  names(has_na) <- colnames(unit$x)

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

  # The names of the forests, which is what a per-forest argument may be keyed
  # by and the order a positional one is read in. The trailing pinned forests
  # standing in for a custom family's nuisance parameters are not among them:
  # they are not additive predictors and nothing about them is the caller's to
  # set per forest.
  n_report <- response[["n_forest"]] - response[["n_aux"]]
  labels <- forest_labels(response[["family"]], response[["opts"]],
                          response[["levels"]], n_report)
  joint <- joint_forests(response[["family"]])

  if (length(forest_formulas) > 1L && joint) {
    arg::err(c("{.arg formula} must be a single formula for this family, not
                {length(forest_formulas)}",
               i = "Its {length(labels)} forests are the levels of one parameter
                  and act together, so they take the same predictors."))
  }

  if (length(forest_formulas) > n_report) {
    arg::err(c("{.arg formula} has {length(forest_formulas)} formulas but this
                family has {n_report} forest{?s}",
               i = "Its forests are {.val {labels}}."))
  }

  # Names on the list of formulas are checked the same way any per-forest
  # argument's are, and reorder it, so that `list(log_sd = ~ x2, mean = y ~ x1)`
  # means what it says.
  if (!is_null(names(forest_formulas))) {
    forest_formulas <- resolve_per_forest(forest_formulas, labels, "formula",
                                          default = forest_formulas[[1L]],
                                          joint = joint)
  }

  # One formula given applies to every forest, which is the rule for every
  # per-forest argument and is what makes the ordinary single-formula call reach
  # the engine unchanged. Fewer formulas than forests, but more than one, is not
  # a recycling anyone would mean.
  if (length(forest_formulas) == 1L) {
    forest_formulas <- rep(forest_formulas, n_report)
  }
  else if (length(forest_formulas) < n_report) {
    arg::err(c("{.arg formula} has {length(forest_formulas)} formulas but this
                family has {n_report} forests",
               i = "Give one formula, or {n_report}, or name them:
                  {.val {labels}}."))
  }

  # Matched against the predictors here rather than in `bartisan_control()`,
  # which does not know them. The result is one weight per predictor group, so a
  # factor's dummy columns share the weight given to the term, the way they
  # already share one entry of the sparsity prior. One column per forest, since
  # the weights are also how a forest is held to its own formula.
  engine_control[["split_prior"]] <-
    resolve_split_matrix(control[["split_prior"]], colnames(group_probs),
                         labels, joint,
                         forest_masks(forest_formulas, colnames(group_probs),
                                      if (!missing(data) &&
                                            is.data.frame(data)) data,
                                      response_of(forest_formulas)),
                         response[["n_forest"]])

  engine_control[["gate"]] <- gate_code(control[["gate"]])

  # The rest of the settings the engine keeps one copy of per forest. Each is
  # spread to one value per forest here, so the engine never has to decide what a
  # scalar means, and each keeps its own default where a named argument left a
  # forest out. `k` is not among them: it is a way of writing `sigma_mu`, and
  # that is spread just below.
  for (nm in names(PER_FOREST_DEFAULTS)) {
    engine_control[[nm]] <- per_forest_vector(
      control[[nm]], labels, nm,
      control[[nm]][[1L]] %||% PER_FOREST_DEFAULTS[[nm]], joint)
    engine_control[[nm]] <- rep(engine_control[[nm]],
                                length.out = response[["n_forest"]])
  }

  # One tree count per additive predictor. A scalar is recycled, so the common
  # case reads the same as before; a vector, or one keyed by the forest names,
  # lets a forest that needs less capacity be given less, which is most of what
  # makes `location_scale()` affordable.
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

  # Everything above is a deterministic function of the data and is done once;
  # only the sampler itself is repeated per chain.
  engine <- function(ignored) {
    .bartisan_fit(X = unit$x,
                  has_na = as.integer(has_na),
                  y = response[["y"]],
                  weights = response[["weights"]],
                  offset = response[["offset"]],
                  group_probs = group_probs,
                  family_name = response[["family"]],
                  link = response[["link"]],
                  family_opts = response[["opts"]],
                  control = engine_control,
                  random_spec = random_spec(random),
                  codes = levels_info[["codes"]],
                  cat_col = levels_info[["cat_col"]],
                  n_levels = levels_info[["n_levels"]])
  }

  draws <- {
    if (chains == 1L) engine(1L)
    else combine_chains(run_chains(engine, chains))
  }

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
              has_na = has_na,
              x_transform = control[["x_transform"]])

  if (!is_null(draws[["aux"]])) {
    aux <- draws[["aux"]]
    colnames(aux) <- draws[["aux_names"]]
    out[["aux"]] <- aux
  }

  names(out[["eta"]]) <- predictor_names(out)

  # The random part is indexed by additive predictor too, so it takes the same
  # names -- which is what makes the diagnostics table readable when a family has
  # more than one.
  if (!is_null(out[["ranef"]])) {
    names(out[["ranef"]]) <- predictor_names(out)
    names(out[["tau"]]) <- predictor_names(out)

    labels <- names(random)
    levels_per <- vapply(random, function(z) z[["num_levels"]], integer(1L))

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

  out[["fitted"]] <- fitted_from_eta(out, out[["eta"]], average = TRUE)

  if (chains > 1L) {
    out[["rhat"]] <- chain_diagnostics(out)
  }

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

  if (rlang::is_installed("future.apply")) {
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
    return(function() {
      RNGkind(kind[1L], kind[2L], kind[3L])
      assign(".Random.seed", seed, envir = globalenv())
    })
  }

  function() RNGkind(kind[1L], kind[2L], kind[3L])
}

# Stack the chains into one set of draws, in chain order. The stored forests are
# indexed by a flat position that runs iteration, then forest, then tree, so the
# record offsets of each chain after the first have to be shifted by the total
# length of the ones before it.
combine_chains <- function(fits) {
  first <- fits[[1L]]
  out <- first

  stack <- function(fits, name) {
    do.call(rbind, lapply(fits, `[[`, name))
  }

  stack_list <- function(fits, name) {
    lapply(seq_along(first[[name]]), function(h) {
      do.call(rbind, lapply(fits, function(z) z[[name]][[h]]))
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

  flat <- lapply(fits, `[[`, "forest_flat")
  at <- utils::head(cumsum(c(0L, lengths(flat))), -1L)
  starts <- lapply(seq_along(fits), function(k) {
    fits[[k]][["tree_start"]][-1L] + at[k]
  })

  out[["forest_flat"]] <- unlist(flat, use.names = FALSE)
  out[["tree_start"]] <- c(0L, unlist(starts, use.names = FALSE))


  if (!is_null(first[["mixture_flat"]])) {
    flat <- lapply(fits, `[[`, "mixture_flat")
    at <- utils::head(cumsum(c(0L, lengths(flat))), -1L)
    starts <- lapply(seq_along(fits), function(k) {
      fits[[k]][["mixture_start"]][-1L] + at[k]
    })

    out[["mixture_flat"]] <- unlist(flat, use.names = FALSE)
    out[["mixture_start"]] <- c(0L, unlist(starts, use.names = FALSE))
  }

  out
}

# Convergence and precision diagnostics for the quantities a caller would look
# at: the log likelihood, the leaf scales, the nuisance parameters, and the
# additive predictor at every observation. The last is a vector of length n, so
# it is summarized by its worst value rather than reported in full -- for R-hat
# the largest, for the two effective sample sizes the smallest, since in both
# cases the worst case is what decides whether the run is usable.
#
# Three numbers per quantity, following Vehtari et al. (2021): rank-normalized
# folded split R-hat, and the bulk and tail effective sample sizes. The two
# effective sample sizes are reported separately because they answer different
# questions -- the bulk one governs a posterior mean, the tail one an interval
# endpoint, and a run can easily be adequate for the first and not the second.
# The scalar parameters of a fit, as a named list of draw vectors. These are the
# quantities that are one number per draw whatever the data are, which is what
# makes them the ones a convergence diagnostic or a draws object wants; the
# predictor and the group effects are one number per observation or per level and
# are handled separately.
scalar_draws <- function(object) {
  out <- list(loglik = object[["loglik"]])

  for (h in seq_len(ncol(object[["sigma_mu"]]))) {
    nm <- sprintf("sigma_mu.%s", colnames(object[["sigma_mu"]])[h])
    out[[nm]] <- object[["sigma_mu"]][, h]
  }

  if (!is_null(object[["aux"]])) {
    for (nm in colnames(object[["aux"]])) {
      out[[sprintf("aux.%s", nm)]] <- object[["aux"]][, nm]
    }
  }

  # The scale of each random-effect term is a scalar worth diagnosing; the
  # intercepts themselves are summarized like the predictor, over the worst
  # level, since there is one per level.
  for (h in seq_along(object[["tau"]])) {
    for (r in seq_len(ncol(object[["tau"]][[h]]))) {
      nm <- sprintf("tau.%s.%s", names(object[["tau"]])[h],
                    colnames(object[["tau"]][[h]])[r])
      out[[nm]] <- object[["tau"]][[h]][, r]
    }
  }

  out
}

chain_diagnostics <- function(object) {
  chains <- object[["chains"]]
  per <- nrow(object[["sigma_mu"]]) / chains
  index <- matrix(seq_len(per * chains), nrow = per, ncol = chains)

  shape <- function(x) matrix(x[index], nrow = per, ncol = chains)

  scalars <- scalar_draws(object)

  # The leaf scale is left out of the table, and only out of the table: it is
  # still in `fit$sigma_mu` and still reaches `as_draws()`, so anyone who wants
  # to diagnose it can.
  #
  # It mixes badly, and not for a reason this package can fix. On one dataset the
  # same quantity comes out at rhat 1.12 with 22 effective draws in dbarts (chi
  # hyperprior, slice sampler) and 1.16 with 17 in stochtree (inverse-gamma, an
  # exact Gibbs draw), against 1.19 and 15 here, and the BART package avoids the
  # question by never drawing it. Nothing beats an exact draw, so the sampler is
  # not the cause in any of them. Reported beside the additive predictors at
  # equal status it meant every fit on ordinary data showed a row above any
  # threshold a reader would apply, for a hyperparameter nobody reports and whose
  # disagreement between chains does not reach the fitted function -- on those
  # same fits `eta` has rhat 1.00 and thousands of effective draws.
  scalars <- scalars[!startsWith(names(scalars), "sigma_mu.")]

  rows <- lapply(names(scalars), function(nm) {
    x <- shape(scalars[[nm]])
    data.frame(quantity = nm, rhat = rhat_rank(x), ess_bulk = ess_bulk(x),
               ess_tail = ess_tail(x))
  })

  for (h in seq_along(object[["eta"]])) {
    label <- sprintf("eta.%s", names(object[["eta"]])[h])
    draws <- object[["eta"]][[h]]

    per_obs <- vapply(seq_len(ncol(draws)), function(j) {
      x <- shape(draws[, j])
      c(rhat_rank(x), ess_bulk(x), ess_tail(x))
    }, numeric(3L))

    rows[[length(rows) + 1L]] <- data.frame(
      quantity = sprintf("%s (worst over observations)", label),
      rhat = worst(per_obs[1L, ], max),
      ess_bulk = worst(per_obs[2L, ], min),
      ess_tail = worst(per_obs[3L, ], min))
  }

  for (h in seq_along(object[["ranef"]])) {
    label <- sprintf("ranef.%s", names(object[["ranef"]])[h])
    draws <- object[["ranef"]][[h]]

    per_level <- vapply(seq_len(ncol(draws)), function(j) {
      x <- shape(draws[, j])
      c(rhat_rank(x), ess_bulk(x), ess_tail(x))
    }, numeric(3L))

    rows[[length(rows) + 1L]] <- data.frame(
      quantity = sprintf("%s (worst over levels)", label),
      rhat = worst(per_level[1L, ], max),
      ess_bulk = worst(per_level[2L, ], min),
      ess_tail = worst(per_level[3L, ], min))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
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
  within <- mean(apply(y, 2L, stats::var))

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
  r <- rank(x, ties.method = "average")

  stats::qnorm((r - 3 / 8) / (length(r) - 1 / 4)) |>
    matrix(nrow = nrow(x), ncol = ncol(x))
}

# Rank-normalized, folded, split R-hat (Vehtari, Gelman, Simpson, Carpenter and
# Buerkner 2021). Two diagnostics, maximized: the rank-normalized one catches
# chains that disagree about the middle of the distribution, and the folded one
# -- the same computation applied to the distance from the median -- catches
# chains that agree about the middle and disagree about the spread, which the
# first is blind to.
rhat_rank <- function(x) {
  if (is_null(split_chains(x)) || !all(is.finite(x))) {
    return(NA_real_)
  }

  bulk <- rank_normalize(x) |>
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
ess_from_split <- function(y) {
  draws <- nrow(y)
  chains <- ncol(y)

  if (draws < 4L) {
    return(NA_real_)
  }

  # Biased autocovariance, one column per chain, lags 0 .. draws - 1.
  acov <- vapply(seq_len(chains), function(m) {
    stats::acf(y[, m], lag.max = draws - 1L, type = "covariance",
               plot = FALSE, demean = TRUE)$acf[, 1L, 1L]
  }, numeric(draws))

  if (!is.matrix(acov)) {
    acov <- matrix(acov, ncol = 1L)
  }

  pooled <- rowMeans(acov)
  mean_var <- pooled[1L] * draws / (draws - 1)
  var_plus <- mean_var * (draws - 1) / draws

  if (chains > 1L) {
    var_plus <- var_plus + stats::var(colMeans(y))
  }

  # Guarded with `isTRUE()` for the reason given in `split_rhat()`.
  if (!isTRUE(var_plus > 0) || !isTRUE(mean_var > 0)) {
    return(NA_real_)
  }

  rho <- function(t) {1 - (mean_var - pooled[t + 1L]) / var_plus}

  # Geyer's initial positive sequence: walk the autocorrelations in adjacent
  # pairs and stop at the first pair whose sum goes negative, which is where the
  # estimates stop being informative.
  kept <- c(1, rho(1L))
  t <- 1L

  while (t < draws - 4L && sum(utils::tail(kept, 2L)) > 0) {
    even <- rho(t + 1L)
    odd <- rho(t + 2L)

    if (even + odd >= 0) {
      kept <- c(kept, even, odd)
    }
    else {
      kept <- c(kept, 0, 0)
    }

    t <- t + 2L
  }

  extra <- max(utils::tail(kept, 2L)[1L], 0)

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

  worst(c(ess_from(rank_normalize((x <= q[1L]) * 1)),
          ess_from(rank_normalize((x >= q[2L]) * 1))), min)
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
                object[["levels"]], object[["num_forest"]])
}

# Names of the additive predictors, in the order the engine builds them. The
# first is always the main parameter -- the one a single-forest family would
# have on its own -- and the rest follow in the order documented on
# [bartisan-families]. These are the names that label the columns of most
# outputs, and the names a per-forest argument may be keyed by.
forest_labels <- function(family, opts, levels, n_report) {
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

  if (identical(family, "location_scale")) {
    return(c("mean", "log_sd"))
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
  variables <- all.vars(stats::delete.response(mt))

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

# One tree count per additive predictor. The default depends on the rules,
# because a soft rule makes a tree more expressive: measured on the Friedman
# function, held-out error levels off by 20 trees with soft rules and keeps
# improving to 50 with hard ones. 50 is kept for both because a smaller forest
# mixes worse -- on `lalonde`, four chains at 20 soft trees disagreed by 35% on
# an average contrast against 9% at 50 -- and the Friedman gain at 20 was 5%.
resolve_num_trees <- function(num_trees, n_forest, n_aux = 0L) {
  num_trees <- as.integer(num_trees %or% 50L)

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
