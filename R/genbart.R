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
#' documented at [genbart-families], along with `custom_family()` for a
#' likelihood of your own.
#'
#' @param formula a model formula. The right-hand side lists candidate
#'   predictors; the model finds interactions and nonlinearity on its own, so
#'   `y ~ .` is usually the right specification. Survival families take a
#'   [survival::Surv()] object on the left. A `(1 | group)` term adds a
#'   group-level random intercept, in the notation of \pkg{lme4}; see Details.
#' @param data a data frame containing the variables in `formula`.
#' @param family the response distribution, as a [stats::family] object, one of
#'   the families in [genbart-families], or a name. See Details for what is
#'   supported.
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
#' @param chains how many independent chains to run. More than one requires the
#'   \pkg{future.apply} package and runs them under whatever backend the caller
#'   has planned with [future::plan()] -- `multisession`, `multicore`, a cluster,
#'   or mirai's `mirai_multisession`. The draws are pooled and
#'   [split-R-hat][genbart()] is reported in the `rhat` element. One
#'   `set.seed()` before the call reproduces the whole run whatever the backend,
#'   because each chain is given its own L'Ecuyer stream.
#' @param control a list of sampler and prior settings from
#'   [genbart_control()].
#' @param ... further arguments to [genbart_control()]. They are
#'   merged into `control`, overriding any value given there, so that
#'   `genbart(..., num_trees = 20)` and
#'   `genbart(..., control = genbart_control(num_trees = 20))` are the same
#'   call. Names that are not arguments of [genbart_control()] are an error
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
#' # Soft decision rules
#'
#' With `soft = TRUE`, the default, a decision rule is a logistic gate rather
#' than a step, so an observation reaches every leaf with some weight and the
#' fitted function is smooth. This costs more per iteration, since a leaf now
#' touches every observation rather than the ones inside its cell, and it makes
#' the leaf parameters of a tree dependent on one another. Combining soft rules
#' with a non-conjugate likelihood is an extension of Linero (2025), which
#' leaves it as an open problem; it is handled here by giving the reversible-jump
#' move a bivariate Laplace proposal for the pair of child leaves, which reduces
#' to Linero's independent pair exactly when the rules are hard.
#'
#' Set `soft = FALSE` in [genbart_control()] for the faster hard-rule sampler.
#'
#' # Random intercepts
#'
#' A `(1 | group)` term in the formula adds an intercept per level of `group`,
#' drawn from a common mean-zero normal whose standard deviation is itself drawn
#' under the same half-Cauchy prior the leaf scale uses. Several grouping factors
#' are allowed, and `(1 | a/b)` expands to nesting as it does in \pkg{lme4}:
#'
#' ```r
#' genbart(y ~ x1 + x2 + (1 | school), data = d)
#' genbart(y ~ x1 + (1 | school) + (1 | year), data = d)
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
#' An object of class `genbart`, a list with elements including:
#' \describe{
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
#' }
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
#' @seealso [predict.genbart()], [genbart_control()], [genbart-families]
#'
#' @examples
#' set.seed(1)
#'
#' n <- 200
#' d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
#' d$y <- rbinom(n, 1, plogis(3 * sin(pi * d$x1 * d$x2) - 1))
#'
#' fit <- genbart(y ~ x1 + x2 + x3, data = d, family = binomial(),
#'                control = genbart_control(num_trees = 10, num_burn = 50,
#'                                          num_save = 50, verbose = FALSE))
#' fit
#'
#' head(predict(fit, type = "response"))
#'
#' @export
genbart <- function(formula, data, family = stats::gaussian(), weights = NULL,
                    offset = NULL, subset = NULL,
                    na.action = stats::na.pass,
                    chains = 1L, control = genbart_control(), ...) {

  cl <- match.call()

  warn_unoptimized()

  arg::arg_formula(formula)
  arg::arg_whole_number(chains)
  arg::arg_gte(chains, 1)
  chains <- as.integer(chains)

  if (!inherits(control, "genbart_control")) {
    arg::err("{.arg control} must be the result of {.fn genbart_control}")
  }

  control <- merge_control(control, list(...))

  family <- as_genbart_family(family)

  # Random-effect terms come out of the formula before anything else looks at
  # it. The model frame is built from the version with the bars replaced by
  # ordinary sums, so the grouping variables come along and get the same
  # missing-value handling as everything else; the design matrix is built from
  # the version with the bars removed, so they are not also predictors.
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

  engine_control <- as.list(control)
  engine_control[["gate"]] <- gate_code(control[["gate"]])
  engine_control[["sigma_mu"]] <- control[["sigma_mu"]] %or%
    (3 * response[["eta_scale"]] / (control[["k"]] *
                                      sqrt(control[["num_trees"]])))

  if (length(engine_control[["sigma_mu"]]) != response[["n_forest"]]) {
    engine_control[["sigma_mu"]] <- rep(engine_control[["sigma_mu"]],
                                        length.out = response[["n_forest"]])
  }

  # Everything above is a deterministic function of the data and is done once;
  # only the sampler itself is repeated per chain.
  engine <- function(ignored) {
    .genbart_fit(X = unit$x,
                 has_na = as.integer(has_na),
                 y = response[["y"]],
                 weights = response[["weights"]],
                 offset = response[["offset"]],
                 group_probs = group_probs,
                 family_name = response[["family"]],
                 link = response[["link"]],
                 family_opts = response[["opts"]],
                 control = engine_control,
                 random_spec = random_spec(random))
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
              random = if (length(random) == 0L) NULL else random,
              ranef = draws[["ranef"]],
              tau = draws[["tau"]],
              eta = draws[["eta"]],
              counts = draws[["counts"]],
              sigma_mu = draws[["sigma_mu"]],
              bandwidth = draws[["bandwidth"]],
              loglik = as.vector(draws[["loglik"]]),
              forest_flat = draws[["forest_flat"]],
              tree_start = draws[["tree_start"]],
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
        paste0(z[["label"]], ":", z[["levels"]])
      }), use.names = FALSE)
      colnames(out[["tau"]][[h]]) <- labels
    }
  }
  names(out[["counts"]]) <- predictor_names(out)

  for (h in seq_along(out[["counts"]])) {
    colnames(out[["counts"]][[h]]) <- colnames(group_probs)
  }
  colnames(out[["sigma_mu"]]) <- predictor_names(out)

  class(out) <- "genbart"

  out[["fitted"]] <- fitted_from_eta(out, out[["eta"]], average = TRUE)

  if (chains > 1L) {
    out[["rhat"]] <- chain_diagnostics(out)
  }

  warn_runaway_scale(out, engine_control[["sigma_mu"]])

  out
}

# An unoptimized build of the compiled code runs five to twenty times slower and
# is otherwise indistinguishable, which makes it very easy to draw conclusions
# about the sampler's speed from the wrong numbers. `devtools::load_all()` and
# `devtools::install()` compile without optimization and leave the object files
# behind for a later `R CMD INSTALL` to reuse, so this is not a rare accident.
# Warned once per session, since it is a property of the installation.
warn_unoptimized <- function() {
  if (isTRUE(the$checked_optimized) || .genbart_optimized()) {
    the$checked_optimized <- TRUE
    return(invisible(NULL))
  }

  the$checked_optimized <- TRUE

  cli::cli_warn(c(
    "{.pkg genbart}'s compiled code was built without optimization, which makes
     it 5 to 20 times slower than it should be",
    i = "reinstall from a clean source directory:
         {.code pkgbuild::clean_dll(); R CMD INSTALL --preclean .}",
    i = "{.code genbart:::.genbart_optimized()} reports the state of the
         installed library"))

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
run_chains <- function(engine, chains) {
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    arg::err("running more than one chain needs the {.pkg future.apply} package.
              Install it, or use {.code chains = 1}")
  }

  future.apply::future_lapply(seq_len(chains), engine, future.seed = TRUE,
                              future.packages = "genbart")
}

# Stack the chains into one set of draws, in chain order. The stored forests are
# indexed by a flat position that runs iteration, then forest, then tree, so the
# record offsets of each chain after the first have to be shifted by the total
# length of the ones before it.
combine_chains <- function(fits) {
  first <- fits[[1L]]
  out <- first

  stack <- function(name) {
    do.call(rbind, lapply(fits, `[[`, name))
  }

  stack_list <- function(name) {
    lapply(seq_along(first[[name]]), function(h) {
      do.call(rbind, lapply(fits, function(z) z[[name]][[h]]))
    })
  }

  out[["eta"]] <- stack_list("eta")
  out[["counts"]] <- stack_list("counts")
  out[["sigma_mu"]] <- stack("sigma_mu")
  out[["bandwidth"]] <- stack("bandwidth")
  out[["loglik"]] <- do.call(rbind, lapply(fits, function(z) z[["loglik"]]))

  if (!is_null(first[["aux"]])) {
    out[["aux"]] <- stack("aux")
  }

  if (!is_null(first[["ranef"]])) {
    out[["ranef"]] <- stack_list("ranef")
    out[["tau"]] <- stack_list("tau")
  }

  flat <- vector("list", length(fits))
  starts <- vector("list", length(fits))
  offset <- 0

  for (k in seq_along(fits)) {
    flat[[k]] <- fits[[k]][["forest_flat"]]
    # Each chain's tree_start begins with a zero, which is only the first
    # chain's; the rest are that chain's offsets shifted past what came before.
    starts[[k]] <- fits[[k]][["tree_start"]][-1L] + offset
    offset <- offset + length(flat[[k]])
  }

  out[["forest_flat"]] <- unlist(flat, use.names = FALSE)
  out[["tree_start"]] <- c(0L, unlist(starts, use.names = FALSE))

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
chain_diagnostics <- function(object) {
  chains <- object[["chains"]]
  per <- nrow(object[["sigma_mu"]]) / chains
  index <- matrix(seq_len(per * chains), nrow = per, ncol = chains)

  shape <- function(x) matrix(x[index], nrow = per, ncol = chains)

  scalars <- list(loglik = object[["loglik"]])

  for (h in seq_len(ncol(object[["sigma_mu"]]))) {
    nm <- paste0("sigma_mu.", colnames(object[["sigma_mu"]])[h])
    scalars[[nm]] <- object[["sigma_mu"]][, h]
  }

  if (!is_null(object[["aux"]])) {
    for (nm in colnames(object[["aux"]])) {
      scalars[[paste0("aux.", nm)]] <- object[["aux"]][, nm]
    }
  }

  # The scale of each random-effect term is a scalar worth diagnosing; the
  # intercepts themselves are summarized like the predictor, over the worst
  # level, since there is one per level.
  for (h in seq_along(object[["tau"]])) {
    for (r in seq_len(ncol(object[["tau"]][[h]]))) {
      nm <- paste0("tau.", names(object[["tau"]])[h], ".",
                   colnames(object[["tau"]][[h]])[r])
      scalars[[nm]] <- object[["tau"]][[h]][, r]
    }
  }

  rows <- lapply(names(scalars), function(nm) {
    x <- shape(scalars[[nm]])
    data.frame(quantity = nm, rhat = rhat_rank(x), ess_bulk = ess_bulk(x),
               ess_tail = ess_tail(x))
  })

  for (h in seq_along(object[["eta"]])) {
    label <- paste0("eta.", names(object[["eta"]])[h])
    draws <- object[["eta"]][[h]]

    per_obs <- vapply(seq_len(ncol(draws)), function(j) {
      x <- shape(draws[, j])
      c(rhat_rank(x), ess_bulk(x), ess_tail(x))
    }, numeric(3L))

    rows[[length(rows) + 1L]] <- data.frame(
      quantity = paste0(label, " (worst over observations)"),
      rhat = worst(per_obs[1L, ], max),
      ess_bulk = worst(per_obs[2L, ], min),
      ess_tail = worst(per_obs[3L, ], min))
  }

  for (h in seq_along(object[["ranef"]])) {
    label <- paste0("ranef.", names(object[["ranef"]])[h])
    draws <- object[["ranef"]][[h]]

    per_level <- vapply(seq_len(ncol(draws)), function(j) {
      x <- shape(draws[, j])
      c(rhat_rank(x), ess_bulk(x), ess_tail(x))
    }, numeric(3L))

    rows[[length(rows) + 1L]] <- data.frame(
      quantity = paste0(label, " (worst over levels)"),
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

  if (!(within > 0)) {
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
  z <- stats::qnorm((r - 3 / 8) / (length(r) - 1 / 4))
  matrix(z, nrow = nrow(x), ncol = ncol(x))
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

  bulk <- split_rhat(rank_normalize(x))
  folded <- split_rhat(rank_normalize(abs(x - stats::median(x))))

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

  if (!(var_plus > 0) || !(mean_var > 0)) {
    return(NA_real_)
  }

  rho <- function(t) 1 - (mean_var - pooled[t + 1L]) / var_plus

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

  extra <- if (utils::tail(kept, 2L)[1L] > 0) utils::tail(kept, 2L)[1L] else 0

  # Force the paired sums to be non-increasing, which is what makes the
  # estimator conservative rather than merely unbiased. With too few kept lags
  # for a second pair there is nothing to compare against.
  last <- length(kept) - 3L
  pairs <- if (last >= 2L) seq(2L, last, by = 2L) else integer(0)

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
  ess_from(rank_normalize(x))
}

# Tail effective sample size: the smaller of the two effective sample sizes for
# the indicator that a draw falls below the 5% and above the 95% quantile. It is
# reported separately because a chain can be perfectly adequate for a posterior
# mean and nowhere near adequate for an interval endpoint -- the mean is an
# average over every draw, and a tail quantile depends on the few draws out
# there.
ess_tail <- function(x) {
  lower <- stats::quantile(x, 0.05, names = FALSE, na.rm = TRUE)
  upper <- stats::quantile(x, 0.95, names = FALSE, na.rm = TRUE)

  worst(c(ess_from(rank_normalize((x <= lower) * 1)),
          ess_from(rank_normalize((x >= upper) * 1))), min)
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

  if (length(at) == 0L) {
    return(invisible(NULL))
  }

  cli::cli_warn(c(
    "the leaf scale settled {round(max(ratio[at]))} times above its prior
     median, which usually means the response is close to separable by the
     predictors",
    i = "the additive predictor is then only weakly identified; fix the scale
         with {.code update_sigma_mu = FALSE} in {.fn genbart_control} if the
         draws look unstable"))

  invisible(NULL)
}

# Names for the additive predictors, which are the columns of most outputs.
predictor_names <- function(object) {
  family <- object[["family"]][["family"]]

  if (identical(family, "multinomial")) {
    if (isTRUE(object[["family_opts"]][["symmetric"]])) {
      return(object[["levels"]])
    }

    return(object[["levels"]][-1L])
  }

  if (identical(family, "location_scale")) {
    return(c("mean", "log_sd"))
  }

  if (family %in% c("zip", "zinb")) {
    return(c("count", "zero"))
  }

  if (identical(family, "custom") && object[["num_forest"]] > 1L) {
    return(paste0("eta", seq_len(object[["num_forest"]])))
  }

  "eta"
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

  if (length(columns) == 0L) {
    return(mf)
  }

  keep <- stats::complete.cases(mf[columns])

  if (all(keep)) {
    return(mf)
  }

  cli::cli_warn("dropping {sum(!keep)} row{?s} with a missing response, weight
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

  if (length(categorical) > 0L) {
    contrasts <- lapply(categorical, function(nm) {
      stats::contrasts(as.factor(mf[[nm]]), contrasts = FALSE)
    })
    names(contrasts) <- categorical
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
    cli::cli_warn("dropping {sum(!varies)} constant predictor column{?s}:
                   {.val {colnames(x)[!varies]}}")
    x <- x[, varies, drop = FALSE]
    assign <- assign[varies]
  }

  list(x = x, assign = assign, term_labels = predictors,
       contrasts = contrasts)
}

unit_transform <- function(x, type) {
  maps <- lapply(seq_len(ncol(x)), function(j) make_unit_map(x[, j], type))
  names(maps) <- colnames(x)

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

# Settings passed to `genbart()` through `...` are control arguments, as in
# `glm()`. They are merged with whatever the caller supplied to
# `genbart_control()` and the whole thing is re-validated, so a bad value passed
# this way fails the same way it would have failed there.
merge_control <- function(control, dots) {
  if (length(dots) == 0) {
    return(control)
  }

  allowed <- names(formals(genbart_control))
  nms <- names(dots)

  if (is_null(nms) || !all(nzchar(nms))) {
    arg::err("arguments passed to {.fn genbart} in {.arg ...} must be named")
  }

  if (anyDuplicated(nms) > 0) {
    arg::err("arguments passed to {.fn genbart} in {.arg ...} must have distinct names")
  }

  bad <- setdiff(nms, allowed)

  if (length(bad) > 0) {
    arg::err(paste0("{.val {bad}} {?is/are} not {?an argument/arguments} of ",
                    "{.fn genbart_control}"))
  }

  args <- attr(control, "supplied") %or% list()

  # Single-bracket assignment from a list, rather than `modifyList()`, so that an
  # explicit `NULL` (a legitimate value for `sigma_mu` and `alpha_scale`) is kept
  # rather than removing the entry.
  args[nms] <- dots

  do.call(genbart_control, args)
}

# The engine takes the gate as an integer, matching GateShape in node.h. Kept in
# one place so the two orderings cannot drift apart.
gate_code <- function(gate) {
  match(gate, c("logistic", "smoothstep", "smootherstep")) - 1L
}
