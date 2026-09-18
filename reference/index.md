# Package index

## Overview

The package’s own page: what it fits, and a table pointing from the task
at hand to the function that does it.

- [`bartisan-package`](https://ngreifer.github.io/bartisan/reference/bartisan-package.md)
  : bartisan: Generalized Bayesian Additive Regression Trees

## Fitting a Model

[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
is the entry point and
[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) is its
counterpart for a treatment effect. Everything else in this section
changes what either one fits: the likelihood, the settings the sampler
runs under, and the one term that gives a predictor a coefficient of its
own.

- [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
  : Fit a generalized Bayesian additive regression trees (BART) model

- [`negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  [`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  [`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  [`dpm_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  [`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  [`weibull_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  [`loglogistic_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  [`lognormal_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  [`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  [`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  [`Gamma_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  [`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  [`zi_negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  [`Beta()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  [`ordbeta()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  [`tweedie()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  : Response families for generalized BART

- [`custom_family()`](https://ngreifer.github.io/bartisan/reference/custom_family.md)
  : Fit a model to a likelihood written in R

- [`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
  :

  Sampler and prior settings for
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)

- [`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) :
  Bayesian causal forests

- [`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) : Give a
  predictor a varying coefficient

## Predictions and Effects

What the fitted model says about units, real or counterfactual. Every
quantity here is computed from the posterior draws, so the summaries
carry credible intervals rather than standard errors, and
[`predict()`](https://rdrr.io/r/stats/predict.html) hands back the draws
themselves. Calling
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a fit draws
a partial dependence plot, and that method is documented on the last
page in this section.

- [`predict(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
  : Predictions from a generalized BART model
- [`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
  [`print(`*`<bartisan_effect>`*`)`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
  [`plot(`*`<bartisan_effect>`*`)`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
  : Causal effects from a fitted model
- [`formula(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
  [`terms(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
  [`model.frame(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
  [`nobs(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
  [`family(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
  [`get_predict(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
  [`get_group_names(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
  [`get_coef(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
  [`set_coef(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
  [`get_vcov(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
  [`get_data(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
  : Counterfactual estimands with marginaleffects
- [`partial_dependence()`](https://ngreifer.github.io/bartisan/reference/partial_dependence.md)
  [`print(`*`<bartisan_partial>`*`)`](https://ngreifer.github.io/bartisan/reference/partial_dependence.md)
  [`plot(`*`<bartisan_partial>`*`)`](https://ngreifer.github.io/bartisan/reference/partial_dependence.md)
  [`plot(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/partial_dependence.md)
  : Partial dependence on one or two predictors

## Describing a Fitted Model

A forest has no coefficient table, so what takes its place is the
printed summary, how the splitting rules were spent, and what the parts
of the model that do have parameters were drawn to be. The last page
holds the print and plot methods a
[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) fit adds
on top of these.

- [`print(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/print.bartisan_fit.md)
  [`summary(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/print.bartisan_fit.md)
  [`print(`*`<summary.bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/print.bartisan_fit.md)
  : Summarize a generalized BART model
- [`variable_importance()`](https://ngreifer.github.io/bartisan/reference/variable_importance.md)
  [`plot(`*`<bartisan_importance>`*`)`](https://ngreifer.github.io/bartisan/reference/variable_importance.md)
  : How often each predictor is used
- [`coef(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/coef.bartisan_fit.md)
  : Varying coefficients
- [`ranef(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/ranef.bartisan_fit.md)
  : Group intercepts from a random-effect term
- [`error_density()`](https://ngreifer.github.io/bartisan/reference/error_density.md)
  [`plot(`*`<bartisan_error_density>`*`)`](https://ngreifer.github.io/bartisan/reference/error_density.md)
  : Error distribution of a Dirichlet process mixture fit
- [`print(`*`<bcf_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/print.bcf_fit.md)
  [`plot(`*`<bcf_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/print.bcf_fit.md)
  : Methods for Bayesian causal forest fits

## Checking Convergence and Fit

Whether the sampler explored the posterior, and whether the model
describes the data, which are worth settling before anything above is
read off a fit. The second page collects the methods that *rstantools*,
*loo*, *bayesplot*, *posterior* and *performance* call on a fit, along
with the `stats` methods that report the same quantities, so posterior
predictive draws, cross-validation, the prior the fit was given, and the
fit statistics are all on it.

- [`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
  : Check whether a fit converged and mixed
- [`posterior_predict(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`posterior_epred(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`posterior_linpred(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`log_lik(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`simulate(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`fitted(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`residuals(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`weights(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`sigma(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`prior_summary(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`print(`*`<bartisan_prior_summary>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`loo(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`waic(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`kfold(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`pp_check(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`as_draws(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`r2_posterior(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`r2(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  [`model_performance(`*`<bartisan_fit>`*`)`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
  : Interfaces to other packages

## Data

The observational study the examples and most of the vignettes are
fitted to.

- [`rhc`](https://ngreifer.github.io/bartisan/reference/rhc.md) : Right
  heart catheterization in critically ill patients
