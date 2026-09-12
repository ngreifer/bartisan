#' @keywords internal
#'
#' @details
#' Fits Bayesian additive regression trees (BART) the way [stats::glm()] fits a
#' generalized linear model: a formula, a data frame, and a family. The forest
#' replaces the linear predictor, so nothing has to be said about which terms
#' enter the model, which are curved, or which interact.
#'
#' ```r
#' fit <- bartisan(y ~ ., data = d)
#' ```
#'
#' ## What to Reach For
#'
#' Most of what a fitted model is used for is handled by a package that already
#' does that job well, and \pkg{bartisan} registers the methods those packages
#' need rather than reimplementing them. This table is the map.
#'
#' | Task | Use |
#' | --- | --- |
#' | fit a model | [bartisan()] |
#' | choose a likelihood | [bartisan-families], `vignette("families")` |
#' | change the sampler's settings | [bartisan_control()] |
#' | predict for new data | [predict.bartisan_fit()] |
#' | a prediction **with an interval** | `marginaleffects::predictions()` |
#' | an interval for a **new observation**, noise included | [posterior_predict()][bartisan-interop] |
#' | how much a predictor moves the outcome | `marginaleffects::avg_comparisons()` |
#' | a partial dependence plot | [partial_dependence()], or `marginaleffects::plot_predictions()` for more control over the grid |
#' | which predictors the forest uses | [variable_importance()] |
#' | to tell the prior which predictors matter | `split_prior` in [bartisan_control()] |
#' | to give one forest its own predictors or settings | a list of formulas, and per-forest arguments; see [bartisan-families] |
#' | a coefficient that varies with the other predictors | [vc()] in the formula, then [coef()] |
#' | a treatment effect with its own prior | [bcf()], `vignette("causal")` |
#' | the ATE, ATT or the effect for each unit | [estimate_effect()] |
#' | to check it converged and mixed | [diagnose()], then [as_draws()][bartisan-interop] with \pkg{bayesplot} |
#' | to check it fits | [pp_check()][bartisan-interop], [residuals()][bartisan-interop] |
#' | to compare two models | [loo()][bartisan-interop] |
#' | survival data | [ph()], [dpm_aft()], `vignette("survival")` |
#' | a likelihood of one's own | [custom_family()] |
#'
#' ## Before the First Fit
#'
#' Three things are worth knowing first, and none of them requires knowing
#' anything about Bayesian statistics.
#'
#' **There are no coefficients.** A forest has no slope to read off, so the
#' effect of a predictor is found by asking the fitted model what it predicts
#' under one value of that predictor and under another, and taking the
#' difference. `marginaleffects::avg_comparisons(fit, variables = "x")` does
#' exactly that, and reports an interval with it. This is a better habit than
#' reading coefficients even where coefficients exist (i.e., in a linear model),
#' and here it is the only habit available.
#'
#' **The intervals mean what they appear to mean.** A 95% interval from any of
#' the above is the range the model considers most plausible given the data and
#' the model. It already includes the uncertainty from not knowing the shape of
#' the relationship, which is the part a linear model leaves out by assuming it
#' away.
#'
#' **The defaults are meant to be used.** The settings in [bartisan_control()]
#' are there for the analyses that need them; the priors and the number of trees
#' are chosen to work across a wide range of problems, and tuning them is rarely
#' where the gains are. Choosing the right [family][bartisan-families] matters
#' much more.
#'
#' ## What the Fit Does Not Do
#'
#' The model is flexible about the shape of the relationship and about nothing
#' else. It will not establish that a predictor is a cause, that the sample
#' represents the population, or that the outcome was measured well. Note that a
#' forest fitted to confounded data returns a confounded answer with a tight
#' interval around it.
#'
#' @seealso
#' [bartisan()] to fit a model; [bartisan-families] to choose a likelihood;
#' [bartisan-marginaleffects] and [bartisan-interop] for the packages that read a
#' fit; `vignette("bartisan")` for how the sampler works.
#'
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats setNames
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang .data
#' @useDynLib bartisan, .registration = TRUE
## usethis namespace: end
NULL
