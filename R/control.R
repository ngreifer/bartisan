#' Sampler and prior settings for `bartisan()`
#'
#' @description
#' Collects the tuning constants of the sampler and the hyperparameters of the
#' tree prior, for passing to the `control` argument of [bartisan()]. The
#' defaults follow Linero (2025) and, for the soft decision rules, Linero and
#' Yang (2018), and are meant to be usable without adjustment.
#'
#' @param chains `numeric`; how many independent chains to run. Default is 1.
#'   The draws are pooled, and [diagnose()] computes split-R-hat from them on
#'   request rather than the fit carrying it. With the \pkg{future.apply}
#'   package installed the chains run in parallel under whatever backend
#'   \pkgfun{future}{plan} has set. One `set.seed()` before the call reproduces
#'   the run either way, since each chain gets its own L'Ecuyer stream. See
#'   `vignette("diagnostics")`.
#' @param num_burn `numeric`; the number of warmup iterations to discard.
#'   Default is 200, which measurement finds ample for the samplers here. Warmup
#'   is where the trees grow into the data and the hyperparameters find their
#'   scale, so raising it buys convergence rather than precision: increase it
#'   when [diagnose()]'s `rhat` says the chains have not agreed, and for a
#'   `dpm()` fit, whose component-count state settles more slowly than a forest
#'   does. See Details.
#' @param num_draws `numeric`; the number of draws to keep. Default is 800.
#'   These are what every estimate and interval is computed from, so raising it
#'   narrows Monte Carlo error; increase it when `ess_bulk` or `ess_tail` is
#'   small relative to what the reported quantity needs.
#' @param num_thin `numeric`; keep one draw in every `num_thin` after warmup.
#'   Default is 1 to keep every draw. Thinning discards draws to make the kept
#'   ones less correlated, which costs information and is worth it only to hold
#'   down the memory a long chain would otherwise take: for a given amount of
#'   computing, more draws beat fewer less-correlated ones, so it is best to
#'   avoid setting `num_thin`.
#' @param num_trees `numeric`; the number of trees, as one number for every
#'   forest or one per additive predictor. Default is 50. A family with more
#'   than one additive predictor takes a vector, which can be worth using
#'   because the forests are neither equally expensive nor in need of equal
#'   capacity. See Details.
#' @param gate string; the shape of the decision rule at each node. Allowable
#'   options include `"smoothstep"` (the default), `"smootherstep"`,
#'   `"logistic"`, and `"hard"` (or equivalently `"step"`). `"hard"` gives the
#'   step functions of standard BART; the other three give soft rules, in which
#'   every observation reaches every leaf with some weight and the fit is
#'   smooth. Soft rules are slower but improve predictive performance and yield
#'   smooth fits. Which soft gate is chosen matters much less than that one is.
#'   See Details.
#' @param sparsity `logical` or string; the prior on which predictors are split
#'   on. `TRUE` (the default) is the Dirichlet sparsity prior of Linero (2018)
#'   (i.e., DART), which can drop a predictor from the forest entirely, and
#'   `FALSE` gives every predictor the same splitting probability, which is
#'   classic BART. The strings `"none"`, `"weak"`, `"moderate"`, and `"strong"`
#'   name four strengths, `"none"` equal to `FALSE` and `"moderate"` to `TRUE`.
#'   Note that this sets `update_s`, `update_alpha`, `alpha_shape_1` and
#'   `alpha_shape_2` together, and that supplying any of those overrides it.
#'   Read the trade-off in Details before turning it off.
#' @param k *Advanced.* `numeric`; controls the leaf prior. The prior standard
#'   deviation of a forest is `3 / k` times the natural scale of its additive
#'   predictor, so larger values shrink the fit harder toward the intercept-only
#'   model. Default is 2.
#' @param bandwidth *Advanced.* `numeric`; the prior mean of the gate bandwidth
#'   of a soft rule, on the scale of the transformed predictors, which lie in
#'   `[0, 1]`. Smaller values approach hard rules. Default is .1. Ignored when
#'   `gate = "hard"`. See also `update_bandwidth` and `bandwidth_every`.
#' @param split_prior `numeric`; a named vector of relative prior weights on the
#'   predictors, keyed by their names in the formula. A predictor not named gets
#'   a weight of 1, and the probability of splitting on one is its weight over
#'   the total. Weights must be finite and nonnegative; zero keeps the predictor
#'   out of every tree, and naming a predictor the model does not have is an
#'   error. Default is `NULL` to weight every predictor equally. Note that
#'   setting this overrides `sparsity`; see Details.
#' @param share_sparsity *Advanced.* `logical`; for a family with more than one
#'   additive predictor, whether the forests draw their splitting proportions
#'   from one pooled Dirichlet rather than one each. Default is `FALSE`. `TRUE`
#'   is an assumption about the data rather than a free improvement, and
#'   requires that the proportions be drawn at all and over the same predictors,
#'   so `sparsity` must not be `FALSE` for the forests that are to share.
#'   Ignored for a family with a single forest. See Details.
#' @param categorical *Advanced.* string; how a splitting rule divides the
#'   levels of a factor. Allowable options include `"subset"` (the default),
#'   which draws a subset of the levels still available at the node and sends
#'   those left, and `"onehot"`, which is what most BART implementations do: it
#'   splits on one indicator column, peeling a single level off the rest.
#' @param augment *Advanced.* `logical`; whether to rewrite the likelihood as the
#'   margin of a Gaussian or a Poisson one, which makes the target a shape the
#'   sampler can exploit. The posterior is the same either way, so this is a
#'   sampling setting rather than a modeling one. Default is `TRUE`, which
#'   rewrites wherever a rewriting exists: the binomial, ordinal, multinomial,
#'   negative binomial, zero-inflated, and survival families. A rewriting is
#'   faster and has not been measured to mix worse per second, so the default is
#'   rarely worth changing; see `vignette("implementation")`.
#' @param x_transform string; how numeric predictors are mapped to `[0, 1]`
#'   before any rule sees them. Allowable options include `"smoothcdf"` (the
#'   default), `"quantile"`, and `"range"`. The map decides where a cutpoint can
#'   fall, the cutpoint prior being uniform on a node's live range in the
#'   transformed coordinate, and for soft rules it is also the scale `bandwidth`
#'   is measured on. `"smoothcdf"` smooths the empirical distribution function
#'   with an Epanechnikov kernel, `"quantile"` uses the empirical distribution
#'   function itself, and `"range"` rescales linearly. Note that `"range"` is
#'   the one to use when a slope is the quantity of interest, being the only
#'   affine map of the three; see `vignette("effects")`. Where the other two
#'   fail, and why the default is the default, is in
#'   `vignette("implementation")`.
#' @param gamma,beta *Advanced.* `numeric`; the branching probability at depth
#'   `d` is `gamma * (1 + d)^(-beta)`. Defaults are .95 and 2.
#' @param sigma_mu *Advanced.* `numeric`; the prior median of the leaf standard
#'   deviation, one value per additive predictor. Default is `NULL` to derive it
#'   from `k` and that forest's own tree count.
#' @param update_sigma_mu *Advanced.* `logical`; whether to draw the leaf
#'   standard deviation under a half-Cauchy prior rather than fixing it. Default
#'   is `TRUE`. `FALSE` is worth reaching for when a binary fit mixes badly:
#'   where the predictors separate the response well the leaf scale is barely
#'   identified, and it wanders rather than settling, which drags the effective
#'   sample size of everything built on it down with it. Fixing it can raise the
#'   effective sample size of the additive predictor and the log likelihood
#'   several times over at little cost in fit; see `vignette("implementation")`.
#' @param sigma_mu_ramp *Advanced.* `numeric`; the fraction of warmup over which
#'   the leaf standard deviation is raised from near zero to its target. Default
#'   is .25; set to 0 to disable. Linero (2025) describes this as essential:
#'   started at its full value, the sampler can settle early into a poor
#'   configuration and fail to move.
#' @param update_tau *Advanced.* `logical`; whether to draw the standard
#'   deviation of each random-effect term under the same half-Cauchy prior the
#'   leaf scale uses, rather than fixing it at that prior's median. Default is
#'   `TRUE`. Only relevant when the formula has a `(1 | group)` term.
#' @param update_bandwidth *Advanced.* `logical`; whether to draw the bandwidth
#'   of each tree rather than holding it at `bandwidth`. Default is `TRUE`.
#' @param bandwidth_every *Advanced.* `numeric`; how many sweeps between
#'   bandwidth draws for a given tree. Default is 1. The bandwidth is one scalar
#'   per tree, drawn by an adaptive random walk, and every attempt costs a full
#'   rebuild of the tree's memberships, which is the largest single item in a
#'   soft-rule fit, so drawing it less often trades mixing in that one parameter
#'   for speed. Raising it is reasonable when the mean function is known to be
#'   smooth and the fit is compute-bound; see Details.
#' @param alpha *Advanced.* `numeric`; the concentration of the Dirichlet prior
#'   on the splitting proportions, where smaller values concentrate splits on
#'   fewer predictors. Default is `NULL`, which uses 1 as the starting value for
#'   a parameter that is then drawn.
#' @param alpha_scale,alpha_shape_1,alpha_shape_2 *Advanced.* `numeric`; the
#'   prior on `alpha`, in which `alpha / (alpha + alpha_scale)` is
#'   Beta(`alpha_shape_1`, `alpha_shape_2`). Default for `alpha_scale` is `NULL`
#'   to use the number of predictor groups; the two shapes default to whatever
#'   `sparsity` implies.
#' @param update_s,update_alpha *Advanced.* `logical`; whether to draw the
#'   splitting proportions and their concentration. Defaults are `NULL` to
#'   follow `sparsity`. Turning both off recovers a uniform prior over
#'   predictors, which is what `sparsity = FALSE` does.
#' @param verbose *Advanced.* `logical`; whether to print progress to the
#'   console while sampling. Default is `FALSE`. For a progress bar instead, see
#'   the Progress section below, which needs no argument here.
#' @param num_print *Advanced.* `numeric`; how many iterations between the
#'   reports `verbose` prints. Default is 100.
#' @param block_eval *Validation.* `logical`; whether to evaluate the likelihood
#'   one leaf at a time rather than one observation at a time. Default is
#'   `FALSE`. A family built by [custom_family()] does this regardless, since it
#'   must call back into R; setting it for a compiled family produces the same
#'   draws from the same seed at somewhat greater cost.
#' @param exact_quadratic *Validation.* `logical`; whether to use the closed
#'   forms that a target quadratic in the additive predictor allows, in which
#'   one pass over a node determines the log target everywhere, so that the
#'   Laplace approximation is the conditional posterior rather than an
#'   approximation to it. Default is `TRUE`, which is what makes a Gaussian
#'   response, or any of the rewritings in `augment`, cheap. `FALSE` falls back
#'   on the general path; the two agree, at greater cost.
#' @param generic_accumulate *Validation.* `logical`; whether to accumulate a
#'   leaf's sums through the family's virtual interface rather than through its
#'   own statically dispatched loop. Default is `FALSE`. The two compute the
#'   same thing; the second lets the compiler inline the family's arithmetic,
#'   which is most of the remaining per-observation cost.
#'
#' @returns
#' A `<bartisan_control>` object, a list containing the supplied settings and
#' the defaults for those not supplied, for passing to the `control` argument of
#' [bartisan()].
#'
#' @details
#' The arguments are roughly ordered by how often they are worth changing. Those
#' from `chains` to `x_transform` are modeling decisions, in that changing one
#' changes what is fitted or how long it is fitted for, and the chain lengths
#' and `num_trees` come first because they are adjusted most. Those marked
#' *Advanced.* are hyperparameters of a prior that one of the first group
#' summarizes, or switches whose default is almost always right. The three
#' marked *Validation.* compute the same posterior more slowly, and are
#' documented so that the checks using them can be read.
#'
#' ## Setting `num_burn`
#'
#' The default is comfortably longer than the transient in every family
#' measured, so raising `num_burn` buys convergence rather than precision: raise
#' it when [diagnose()] says the chains have not agreed. `dpm()` is the
#' exception, as its mixture carries a component-count state that settles more
#' slowly than a forest does, so a longer warmup is worth having there and
#' especially when the error distribution itself is the object of interest.
#'
#' ## Setting `num_trees`
#'
#' The performance of the BART model is relatively insensitive to the number of
#' trees. More trees create a more complicated model, which can mix poorly. A
#' model with soft trees needs fewer trees than one with hard trees, but the
#' default is the same for both. Generally, you don't need to change the number
#' of trees for a model with a single forest; the defaults tend to do well,
#' which is a strength of BART. The leaf prior scale divides by the square root
#' of a forest's own tree count, so shrinking one forest does not change the
#' prior on the sum.
#'
#' `num_trees` can be supplied as a vector when the family has more than one
#' forest, where each value corresponds to the number of trees for that forest.
#' The scale forest of location-scale models (e.g., `family = gaussian_ls()`)
#' dominates the run time, its target not being quadratic, and a variance
#' surface carries much less information than a mean surface, so giving it fewer
#' trees runs substantially faster at the same accuracy. It is not the default
#' because how much resolution a variance surface needs depends on the surface,
#' and quietly under-parameterizing it would show up as intervals that are
#' wrong.
#'
#' ## Choosing a gate type with `gate`
#'
#' A gate determines how much of each branch of a tree an observation
#' contributes to, based on its value of the numeric predictor the rule splits
#' on. For a predictor \eqn{x} and a threshold \eqn{c}, the contribution of the
#' observation to the right-hand branch is \eqn{w_R(x; c)}, and the contribution
#' to the left branch is \eqn{w_L(x; c) = 1 - w_R(x; c)}.
#'
#' ### Hard Gates
#'
#' For a hard gate (`gate = "hard"`), which is what is used in traditional BART
#' and most other tree-based models,
#'
#' \deqn{w_R(x; c) = \mathbb{I}(x > c)}
#'
#' where \eqn{\mathbb{I}(\cdot)} is the indicator function. An observation
#' contributes only to one branch. This is the fastest to compute but cannot
#' capture smooth functions well, yielding jagged fits.
#'
#' ### Soft Gates
#'
#' For a soft gate, used in the SoftBart model of Linero and Yang (2018),
#'
#' \deqn{w_R(x; c) = F\left(\frac{x - c}{h}\right)}
#'
#' for some S-shaped function \eqn{F(\cdot)} with a bandwidth \eqn{h}.
#'
#' The options for \eqn{F(\cdot)} are the different values `gate` can take.
#' Computing the contribution weights for soft gates is slower than for hard
#' gates, and the bandwidth is itself a parameter with a Metropolis step per
#' tree per sweep, each rebuilding every membership weight in the tree, which is
#' the single largest item in a soft-rule fit. These features make models with
#' soft gates slower to compute than models with hard gates. The difference
#' among the soft gates (`"smoothstep"`, `"smootherstep"`, and `"logistic"`) is
#' negligible in predictive performance. `"smoothstep"` is computationally
#' fastest and lets some observations have weights of exactly 0 or 1, which is
#' why it is the default; `"logistic"` is the one Linero and Yang (2018)
#' originally describe.
#'
#' `bandwidth_every` controls how often the bandwidth is updated from its
#' starting value specified by `bandwidth` for a soft gate when
#' `update_bandwidth = TRUE` (the default). Drawing the bandwidth is what lets a
#' rule sharpen toward a step, so setting `update_bandwidth = FALSE` is faster
#' and tends to be more accurate on smooth functions, but is much worse on
#' nonsmooth ones. Raising `bandwidth_every` is the middle course, recovering
#' some speed while keeping soft rules, at a real cost in mixing and a small one
#' in accuracy where the mean function jumps.
#'
#' ## Controlling Sparsity
#'
#' `sparsity = TRUE` (the default) is the Dirichlet prior of Linero (2018) on
#' the splitting proportions, and it is a genuine variable-selection prior: it
#' can and does drop a predictor from every tree at once, which is the point of
#' it in the high-dimensional problems it was built for.
#'
#' That has a consequence which is easy to misread. A contrast on a predictor
#' the prior has dropped is exactly zero because in that draw the fit does not
#' depend on that predictor at all, so the posterior of a contrast has an atom
#' at zero whose mass is one minus the predictor's inclusion probability; any
#' summary reporting a median lands on it once it holds half the mass. More
#' trees does not fix it, because a predictor whose splitting proportion has
#' gone small is rarely proposed and so is hard to get back in.
#'
#' How to specify `sparsity` follows from the estimand. For prediction or
#' variable selection, keep the default; reaching past `TRUE` is warranted only
#' when the predictors are many and nearly all are expected to be irrelevant.
#' For a contrast, a partial effect, or a treatment effect, set `sparsity =
#' FALSE`, or use `split_prior`, which cannot drop anything: on a weak signal
#' the sparsity prior attenuates the estimate substantially and its interval
#' covers below its nominal rate, while a strong effect is untouched.
#' `vignette("effects")` works this through and [`bartisan-marginaleffects`]
#' covers the atom.
#'
#' For a varying-coefficient model, the sparsity choice can differ by forest
#' (see more details below). What the prior can drop is a predictor a forest
#' splits on, and in a [vc()] model the treatment is the coefficient rather than
#' one of those, carried by a forest of its own, so no splitting proportion can
#' drop it; the prior on that forest selects among the moderators instead. So
#' `sparsity = c(FALSE, TRUE)`, which disables the sparsity prior for the
#' control function forest and enables it for the varying coefficient forest, is
#' coherent, and in [bcf()] it is the asymmetry worth considering.
#'
#' Sparsity can also be controlled by choosing the sparsity parameter `alpha` or
#' drawing it from a prior determined by the parameters `alpha_scale`,
#' `alpha_shape_1`, and `alpha_shape_2`. The `sparsity` argument is a shortcut
#' to setting these, since the relationship between them and the sparsity they
#' induce is a bit opaque.
#'
#' ### Telling the Prior What Is Already Known
#'
#' `sparsity` and `split_prior` answer different questions and cannot both be in
#' force, so giving `split_prior` turns `sparsity` off. `sparsity` is for when
#' which predictors matter is unknown and the prior is to work it out from the
#' data, and a predictor can be dropped entirely; `split_prior` is for when
#' something is known and is to be honored, with the proportions held at the
#' supplied values.
#'
#' A weight is a statement about relative attention, not about effect size. It
#' changes how often a split on a predictor is proposed, which is a prior, so
#' the data can still overrule it in either direction. Because the weights are
#' fixed, `split_prior` does not accumulate the atom-at-zero mass above, which
#' makes it a reasonable middle course when a particular contrast is the
#' estimand but the predictors are too many to treat alike. One weight per term
#' in the formula, not per column of the design matrix, so a factor is named
#' once and its levels share the weight.
#'
#' ## Arguments That Vary by Forest
#'
#' A family with several additive predictors has one forest per predictor, each
#' with its own prior. Every argument that could mean something different for
#' one of them may be given once, to apply to all, or once per forest, either
#' positionally or keyed by the forest names listed in [`bartisan-families`]. A
#' forest a named argument does not mention keeps that argument's default rather
#' than borrowing another forest's value. That covers `num_trees`, `k`,
#' `sigma_mu`, `sparsity`, `split_prior`, `bandwidth`, `gamma`, `beta`, the four
#' `alpha` arguments, and the three `update_` flags; `formula` works the same
#' way, as [bartisan()] describes.
#'
#' ```r
#' bartisan_control(num_trees = c(mean = 50, log_sd = 10),
#'                  sparsity = c(mean = TRUE, log_sd = FALSE))
#' ```
#'
#' The multinomial families are the exception, for the reason given in
#' [`bartisan-families`]: their forests act as one, so these arguments take a
#' single value.
#'
#' ## Progress
#'
#' `verbose = TRUE` prints a line every `num_print` iterations, which is the
#' whole of what this package decides about progress. A progress bar is
#' \CRANpkg{progressr}'s business, and the sampler reports to it
#' unconditionally: nothing is shown unless a handler is active, so there is no
#' argument to switch on. To request a progress bar, use one of the calls below:
#'
#' ```r
#' progressr::with_progress(
#'   bartisan(y ~ ., data = d, family = gaussian())
#' )
#'
#' # or once, for the session
#' progressr::handlers(global = TRUE)
#' ```
#'
#' The bar is sized for the whole fit, so `chains = 4` fills one bar once rather
#' than four in sequence, and chains running in parallel under \CRANpkg{future}
#' relay their progress back as it arrives.
#'
#' ## Parallelization with \pkg{future}
#'
#' To request parallel processing, specify a future plan, e.g., using
#'
#' ```r
#' future::plan("multisession")
#' ```
#' Parallelization only works with multiple chains (i.e., `chains` set to
#' greater than 1), and each chain is run in parallel. With large fits (many
#' draws or a large sample size), an error might be thrown by
#' \pkgfun{future}{plan}, which refuses a single export above
#' `future.globals.maxSize`, 500 MB by default. The error names the option;
#' raising it is the fix, and running sequentially is the alternative.
#'
#' Parallelization can be used with a \pkg{progressr} progress bar.
#'
#' @references
#' Linero, A. R. (2018). Bayesian regression trees for high-dimensional
#' prediction and variable selection. *Journal of the American Statistical
#' Association*, 113(522), 626--636. \doi{10.1080/01621459.2016.1264957}
#'
#' Linero, A. R. (2025). Generalized Bayesian additive regression trees models:
#' beyond conditional conjugacy. *Journal of the American Statistical
#' Association*, 120(549), 356--369. \doi{10.1080/01621459.2024.2337156}
#'
#' Linero, A. R., Sinha, D., & Lipsitz, S. R. (2020). Semiparametric mixed-scale
#' models using shared Bayesian forests. *Biometrics*, 76(1), 131--144.
#' \doi{10.1111/biom.13107}
#'
#' Linero, A. R., & Yang, Y. (2018). Bayesian regression tree ensembles that
#' adapt to smoothness and sparsity. *Journal of the Royal Statistical Society
#' Series B*, 80(5), 1087--1110. \doi{10.1111/rssb.12293}
#'
#' @seealso
#' [bartisan()], which takes the result as its `control` argument;
#' [`bartisan-families`] for the forest names the per-forest arguments are keyed
#' by
#'
#' `vignette("implementation")` for more details on the BART implementation.
#'
#' @examples
#' data("rhc")
#' set.seed(123)
#'
#' # Settings can be built up once and reused across fits
#' ctrl <- bartisan_control(num_trees = 20, gate = "hard", num_burn = 50,
#'                          num_draws = 50)
#'
#' fit <- bartisan(death ~ . - days, data = rhc, control = ctrl)
#'
#' # The same call, with the settings passed through `...` instead
#' fit2 <- bartisan(death ~ . - days, data = rhc, num_trees = 20,
#'                  gate = "hard", num_burn = 50, num_draws = 50)
#'
#' # A setting given once applies to every forest, and a vector gives each
#' # forest its own value. A variance surface needs less capacity than a
#' # mean surface
#' bartisan_control(num_trees = c(mean = 50, log_sd = 10))
#'
#' # Weighting the splitting prior toward the treatment, which the sparsity
#' # prior would otherwise be free to drop
#' bartisan_control(split_prior = c(rhc = 10))
#'
#' @export
bartisan_control <- function(chains = 1L,
                             num_burn = 200L, num_draws = 800L, num_thin = 1L,
                             num_trees = 50L,
                             gate = "smoothstep",
                             sparsity = TRUE,
                             split_prior = NULL,
                             x_transform = "smoothcdf",
                             k = 2,
                             bandwidth = 0.1,
                             share_sparsity = FALSE,
                             categorical = "subset",
                             augment = TRUE,
                             gamma = 0.95, beta = 2,
                             sigma_mu = NULL, update_sigma_mu = TRUE,
                             sigma_mu_ramp = 0.25,
                             update_tau = TRUE,
                             update_bandwidth = TRUE, bandwidth_every = 1L,
                             alpha = NULL, alpha_scale = NULL,
                             alpha_shape_1 = NULL, alpha_shape_2 = NULL,
                             update_s = NULL, update_alpha = NULL,
                             verbose = FALSE, num_print = 100L,
                             block_eval = FALSE,
                             exact_quadratic = TRUE,
                             generic_accumulate = FALSE) {

  # Record what the caller actually supplied, before any of the arguments below
  # are normalized. `bartisan()` uses this to rebuild the control list when extra
  # settings are passed through its `...`, the way `glm()` does.
  supplied <- mget(as.character(setdiff(names(match.call())[-1], "")),
                   environment())

  for (nm in c("num_burn", "num_draws", "num_thin", "num_print",
               "bandwidth_every")) {
    value <- get(nm)
    arg::arg_whole_number(value, .arg = nm)
    arg::arg_gte(value, 0, .arg = nm)
  }

  # One tree count per forest, or one for all of them. Resolving how many
  # forests there are needs the family, so that is left to `bartisan()`; all that
  # can be checked here is that the values themselves are usable. `NULL` is not
  # among them, and `arg_non_null()` rejects a zero-length vector along with it.
  arg::arg_non_null(num_trees)
  arg::arg_whole_numeric(num_trees)
  arg::arg_no_NA(num_trees)
  arg::arg_gte(num_trees, 1)

  arg::arg_gte(bandwidth_every, 1)
  arg::arg_gte(num_draws, 1)
  arg::arg_gte(num_thin, 1)

  arg::arg_whole_number(chains)
  arg::arg_gte(chains, 1)

  for (nm in c("update_tau", "verbose", "block_eval",
               "exact_quadratic", "generic_accumulate")) {
    value <- get(nm)
    arg::arg_flag(value, .arg = nm)
  }

  # The per-forest settings are checked elementwise, since each may be one value
  # or one per forest, named or not. Whether the length is right for the family
  # is `bartisan()`'s business, because the forests are not known until the
  # family is; all that can be checked here is that the values are usable.
  for (nm in c("update_bandwidth", "update_sigma_mu")) {
    for (value in as.list(get(nm))) {
      arg::arg_flag(value, .arg = nm)
    }
  }

  for (nm in c("update_s", "update_alpha")) {
    for (value in as.list(get(nm))) {
      arg::arg_flag(value, .arg = nm)
    }
  }

  for (nm in c("bandwidth", "k", "gamma", "beta")) {
    for (value in as.list(get(nm))) {
      arg::arg_number(value, .arg = nm)
      arg::arg_gt(value, 0, .arg = nm)
    }
  }

  for (nm in c("alpha", "alpha_scale", "alpha_shape_1", "alpha_shape_2")) {
    for (value in as.list(get(nm))) {
      arg::arg_number(value, .arg = nm)
      arg::arg_gte(value, 0, .arg = nm)
    }
  }

  for (value in as.list(gamma)) {
    arg::arg_lte(value, 1, .arg = "gamma")
  }

  arg::arg_number(sigma_mu_ramp)
  arg::arg_between(sigma_mu_ramp, c(0, 1))

  arg::when_not_null(
    sigma_mu,
    arg::arg_and(
      arg::arg_numeric,
      arg::arg_gt(0)
    )
  )

  x_transform <- arg::match_arg(x_transform, c("smoothcdf", "quantile", "range"))

  # One argument for both the shape of a soft rule's gate and the choice between
  # soft and hard rules, because they are one decision: a hard rule is the
  # limiting case of a soft one and there is no gate shape to pick for it.
  gate <- arg::match_arg(gate, c("smoothstep", "smootherstep", "logistic",
                                 "hard", "step"))
  soft <- !gate %in% c("hard", "step")

  split_prior <- resolve_split_prior(split_prior)
  categorical <- arg::match_arg(categorical, c("subset", "onehot"))

  # A named splitting prior replaces the Dirichlet one rather than seeding it.
  # The two say different things: `sparsity` says the caller does not know which
  # predictors matter and wants the prior to find out, and `split_prior` says
  # they do know and want it honored. Drawing `s` from a Dirichlet centered on
  # the supplied weights would answer neither question, so the weights are held
  # fixed and `sparsity` is ignored. Said out loud only when the caller asked
  # for both, since `sparsity = TRUE` is the default and is not a request.
  if (!is_null(split_prior)) {
    if (!missing(sparsity) && !isFALSE(sparsity)) {
      arg::wrn(c("{.arg split_prior} overrides {.arg sparsity}, which is ignored.",
                 i = "{.arg split_prior} fixes the splitting probabilities; {.arg sparsity} draws them."))
    }
    sparsity <- FALSE
  }

  arg::arg_flag(share_sparsity)

  # Nothing to pool when the proportions are not drawn at all, and saying so is
  # better than accepting a setting that does nothing.
  if (share_sparsity && !any(resolve_sparsity(sparsity)[["update_s"]])) {
    arg::err(c("{.arg share_sparsity} has nothing to share when the splitting
                proportions are not drawn.",
               i = "It pools the counts behind one Dirichlet draw, which
                    {.code sparsity = TRUE} is what asks for."))
  }

  sparse <- resolve_sparsity(sparsity)
  augment <- resolve_augment(augment)

  out <- list(num_trees = stats::setNames(as.integer(num_trees),
                                          names(num_trees)),
              gate = gate,
              soft = soft,
              sparsity = sparsity,
              share_sparsity = share_sparsity,
              # Sharing the tree topology across a family's forests: built,
              # measured, and deliberately not offered. It is the shared forest
              # model of Linero, Sinha and Lipsitz (2020) with independent leaf
              # priors, and the engine implements it in full -- see
              # `update_shared_forests()` in `src/mcmc.cpp` and the guards in
              # `src/model.cpp`.
              #
              # It is not an argument because it is not a setting: it is a
              # different model, and a strong one. Every forest is held to the
              # same partition of the covariate space, which asserts that the
              # components are functions of the same predictors. Measured over
              # 20 replicates, that assertion pays about 1.2x to 1.4x on the
              # median quantity when it is true and costs up to 1.8x when it is
              # false, and the time it saves is 1.00x to 1.12x -- flat in the
              # number of forests, because the leaf refresh dominates and
              # sharing leaves the number of leaf values alone. An asymmetric
              # gamble on an assumption the caller has probably not examined is
              # not something to put behind a flag.
              #
              # Reaching it takes a reference to that environment and an
              # assignment into it (`flags <- bartisan:::the` then
              # `flags$share_forests <- TRUE`; the one-line form is not an
              # assignment R will make), which is what
              # is what `tests/testthat/test-share-forests.R` and the
              # `_dev/shared-topology-*.R` scripts do. `_dev/TASKS.md` has the
              # full measurements and the acceptance-ratio derivation.
              share_forests = isTRUE(the$share_forests),
              split_prior = split_prior,
              categorical = categorical,
              k = k,
              bandwidth = bandwidth,
              chains = as.integer(chains),
              num_burn = as.integer(num_burn),
              num_draws = as.integer(num_draws),
              num_thin = as.integer(num_thin),
              augment = augment,
              x_transform = x_transform,
              gamma = gamma,
              beta = beta,
              sigma_mu = sigma_mu,
              update_sigma_mu = update_sigma_mu,
              sigma_mu_ramp = sigma_mu_ramp,
              update_tau = update_tau,
              update_bandwidth = update_bandwidth,
              bandwidth_every = as.integer(bandwidth_every),
              alpha = alpha %or% sparse[["alpha"]],
              alpha_scale = alpha_scale %or% 0,
              alpha_shape_1 = alpha_shape_1 %or% sparse[["alpha_shape_1"]],
              alpha_shape_2 = alpha_shape_2 %or% sparse[["alpha_shape_2"]],
              update_s = update_s %or% sparse[["update_s"]],
              update_alpha = update_alpha %or% sparse[["update_alpha"]],
              verbose = verbose,
              num_print = as.integer(num_print),
              block_eval = block_eval,
              exact_quadratic = exact_quadratic,
              generic_accumulate = generic_accumulate)

  class(out) <- "bartisan_control"
  attr(out, "supplied") <- supplied

  out
}

# One argument standing in for the four that actually parameterize the Dirichlet
# sparsity prior of Linero (2018). What a caller wants to say is how hard the
# prior should push splits onto a few predictors, and the four hyperparameters
# are a poor way to say it: `alpha / (alpha + P)` is Beta(a1, a2), so the mean
# selection pressure is a2 / (a1 + a2) and moving it means moving two numbers at
# once in opposite directions.
#
# `"none"` leaves the splitting proportions at their uniform prior, which is
# classic BART. The other three draw them, and differ in the prior on the
# concentration: Beta(1, 1) is uniform on the transformed concentration and
# selects gently, Beta(0.5, 1) is Linero's default, and Beta(0.5, 3) pushes
# harder toward a few predictors.
# Vectorized over forests: one value applies everywhere, and several are one per
# forest, positionally or keyed by the forest names. `bartisan()` does the
# spreading, because the forests are not known until the family is, so all this
# does is carry the names through and resolve each element on its own.
resolve_sparsity <- function(sparsity) {
  if (is_null(sparsity)) {
    arg::err("{.arg sparsity} must not be empty")
  }

  levels <- vapply(seq_along(sparsity), function(i) {
    one <- sparsity[[i]]

    arg::arg_or(one, arg::arg_flag, arg::arg_string, .arg = "sparsity")

    if (isTRUE(one)) "moderate"
    else if (isFALSE(one)) "none"
    else arg::match_arg(one, c("none", "weak", "moderate", "strong"),
                        .arg = "sparsity")
  }, character(1L))

  shapes <- vapply(levels, switch, numeric(2L),
                   none = c(0.5, 1),
                   weak = c(1, 1),
                   moderate = c(0.5, 1),
                   strong = c(0.5, 3))

  keep <- function(x) {
    setNames(x, names(sparsity))
  }

  list(update_s = keep(levels != "none"),
       update_alpha = keep(levels != "none"),
       alpha = keep(rep.int(1, length(levels))),
       alpha_shape_1 = keep(shapes[1L, ]),
       alpha_shape_2 = keep(shapes[2L, ]))
}

# Relative prior weights on the predictors, keyed by the name each has in the
# formula. Kept as the caller wrote them here, because the predictors are not
# known until the model frame is built; `bartisan()` matches the names against
# the terms and normalizes. What is checked here is everything that can be
# checked without the data: that it is a named numeric vector of finite positive
# numbers with no duplicate and no empty name.
resolve_split_prior <- function(split_prior) {
  if (is_null(split_prior)) {
    return(NULL)
  }

  arg::when_not_null(split_prior,
                     arg::arg_numeric)

  nm <- names(split_prior)

  if (is_null(nm) || any(!nzchar(nm))) {
    arg::err("{.arg split_prior} must be named, with one name per predictor
              given a weight")
  }

  if (anyDuplicated(nm)) {
    dup <- unique(nm[duplicated(nm)])
    arg::err("{.arg split_prior} names each predictor once; {.val {dup}} {?is/are} repeated")
  }

  if (anyNA(split_prior) || any(!is.finite(split_prior)) ||
      any(split_prior < 0)) {
    arg::err("all values in {.arg split_prior} must be finite and non-negative")
  }

  # Zero is allowed and means what it says: the predictor is never split on, so
  # it stays in the model frame and out of every tree. Whether that leaves
  # anything to split on depends on the predictors the model has, which are not
  # known here; `resolve_split_weights()` checks it.

  split_prior
}

# The engine family names whose likelihood should be rewritten as the margin of a
# Gaussian one. TRUE means the ones where that has been measured to pay, which
# for the negative binomial depends on the decision rules: its rewriting is a
# gain of 2x in effective samples per second with hard rules, where the target's
# exponential form collapses the leaf work to one pass, and a slight loss with
# soft rules, where it does not. Naming a family explicitly always honors the
# request.
# The engine is told which family names may be rewritten, and a fit has one
# family, so per-family control would only ever amount to whether this fit's own
# family is in the list. The argument is the flag; the list is how the engine
# asks the question.
resolve_augment <- function(augment) {
  arg::arg_flag(augment)

  if (!augment) {
    return(character())
  }

  c("binomial", "ordinal", "multinomial", "zip", "zinb", "aft", "negbin")
}
