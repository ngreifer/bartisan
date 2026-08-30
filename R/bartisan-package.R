#' @keywords internal
#'
#' @details
#' Fits Bayesian additive regression trees the way [stats::glm()] fits a
#' generalized linear model: a formula, a data frame, and a family. The forest
#' replaces the linear predictor, so nothing has to be said about which terms
#' enter, which are curved, or which interact. Everything else about the workflow
#' stays where it was.
#'
#' ```r
#' fit <- bartisan(y ~ ., data = d)
#' ```
#'
#' # What to reach for
#'
#' Most of what you will want to do with a fitted model lives in a package that
#' already does that job well, and *bartisan* registers the methods those packages
#' need rather than reimplementing them. This table is the map.
#'
#' | You want to | Use |
#' | --- | --- |
#' | fit a model | [bartisan()] |
#' | choose a likelihood | [bartisan-families], `vignette("families")` |
#' | change the sampler's settings | [bartisan_control()] |
#' | predict for new data | [predict.bartisan_fit()] |
#' | a prediction **with an interval** | `marginaleffects::predictions()` |
#' | an interval for a **new observation**, noise included | [posterior_predict()][bartisan-interop] |
#' | how much a predictor moves the outcome | `marginaleffects::avg_comparisons()` |
#' | a partial dependence plot | `marginaleffects::plot_predictions()` |
#' | which predictors the forest uses | [variable_importance()] |
#' | to check it converged | `fit$rhat`, then [as_draws()][bartisan-interop] with \pkg{bayesplot} |
#' | to check it fits | [pp_check()][bartisan-interop], [residuals()][bartisan-interop] |
#' | to compare two models | [loo()][bartisan-interop] |
#' | survival data | [ph()], [dpm_aft()], `vignette("survival")` |
#' | a likelihood of your own | [custom_family()] |
#'
#' # If you are new to this
#'
#' Three things are worth knowing before the first fit, and none of them requires
#' knowing anything about Bayesian statistics.
#'
#' **There are no coefficients.** A forest has no slope to read off, so the
#' question "what is the effect of `x`" is answered by asking the fitted model
#' what it predicts under one value of `x` and under another, and taking the
#' difference. `marginaleffects::avg_comparisons(fit, variables = "x")` does
#' exactly that, and reports an interval with it. This is a better habit than
#' reading coefficients even when coefficients exist, and here it is the only
#' habit available.
#'
#' **The intervals mean what you would hope.** A 95% interval from any of the
#' above is the range the model considers most plausible, given the data and the
#' model. It already includes the uncertainty from not knowing the shape of the
#' relationship, which is the part a linear model leaves out by assuming it away.
#'
#' **The defaults are meant to be used.** The settings in [bartisan_control()] are
#' there for people who need them; the priors and the number of trees are chosen
#' to work across a wide range of problems, and tuning them is rarely where the
#' gains are. Choosing the right [family][bartisan-families] matters much more.
#'
#' # What the fit does not do for you
#'
#' It is flexible about the shape of the relationship, not about anything else.
#' It will not tell you that a predictor is a cause, that the sample represents
#' the population, or that the outcome was measured well. A forest fitted to
#' confounded data returns a confounded answer with a tight interval around it.
#'
#' @seealso
#' [bartisan()] to fit, [bartisan-families] to choose a likelihood,
#' [bartisan-marginaleffects] and [bartisan-interop] for the packages that read a
#' fit, and `vignette("bartisan")` for how the sampler works.
#'
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats setNames
#' @importFrom Rcpp sourceCpp
#' @useDynLib bartisan, .registration = TRUE
## usethis namespace: end
NULL
