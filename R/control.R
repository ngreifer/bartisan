#' Sampler and prior settings for `genbart()`
#'
#' @description
#' Collects the tuning constants of the sampler and the hyperparameters of the
#' tree prior. The defaults follow Linero (2025) and, for the soft decision
#' rules, Linero and Yang (2018), and are intended to be usable without
#' adjustment.
#'
#' @param num_trees number of trees in each forest. Soft rules make each tree
#'   more expressive, so fewer are needed than for standard BART.
#' @param num_burn number of warmup iterations to discard.
#' @param num_save number of draws to keep.
#' @param num_thin keep one draw in every `num_thin` after warmup.
#' @param soft use soft decision rules, as in the SoftBart model, rather than
#'   the hard rules of standard BART. Soft rules give smoother fits at a higher
#'   cost per iteration, because every observation reaches every leaf with some
#'   weight.
#' @param bandwidth prior mean of the gate bandwidth of a soft rule, on the
#'   scale of the transformed predictors, which lie in `[0, 1]`. Smaller values
#'   approach hard rules. Ignored when `soft = FALSE`.
#' @param update_bandwidth draw the bandwidth of each tree rather than holding
#'   it at `bandwidth`.
#' @param bandwidth_every how many sweeps between bandwidth draws for a given
#'   tree. The bandwidth is one scalar per tree, drawn by an adaptive random
#'   walk, and every attempt costs a full rebuild of the tree's memberships --
#'   the single largest item in a soft-rule fit. Drawing it less often than the
#'   trees themselves trades mixing in that one parameter for time; see Details.
#' @param gate the shape of a soft rule's gate, as a function of the distance
#'   from the cutpoint. All three are cumulative distribution functions, so all
#'   are monotone and under all of them the two children's weights sum to the
#'   parent's. What separates them is how many derivatives the fitted function
#'   has: `"logistic"` is the logistic CDF of Linero and Yang (2018) and is
#'   infinitely differentiable; `"smoothstep"` is the Beta(2, 2) CDF, once
#'   differentiable and supported on a bounded interval; `"smootherstep"` is
#'   Beta(3, 3), twice differentiable and bounded. The two bounded gates are
#'   about 1.4 times faster than the logistic and equally accurate, and are
#'   within noise of each other, so choosing between them is a question of
#'   smoothness rather than speed; see Details. Ignored when `soft = FALSE`.
#' @param k controls the leaf prior. The prior standard deviation of a forest is
#'   `3 / k` times the natural scale of its additive predictor, so larger `k`
#'   shrinks the fit harder towards the intercept-only model.
#' @param sigma_mu prior median of the leaf standard deviation, one value per
#'   additive predictor. The default, `NULL`, derives it from `k` and
#'   `num_trees`.
#' @param update_sigma_mu draw the leaf standard deviation under a half-Cauchy
#'   prior rather than fixing it.
#' @param update_tau draw the standard deviation of each random-effect term under
#'   the same half-Cauchy prior the leaf scale uses, rather than fixing it at
#'   that prior's median. Only relevant when the formula has a `(1 | group)`
#'   term.
#' @param sigma_mu_ramp fraction of warmup over which the leaf standard
#'   deviation is raised from near zero to its target. Linero (2025) describes
#'   this as essential: started at its full value, the sampler can settle early
#'   into a poor configuration and fail to move. Set to `0` to disable.
#' @param gamma,beta the branching probability at depth `d` is
#'   `gamma * (1 + d)^(-beta)`.
#' @param alpha concentration of the Dirichlet prior on the splitting
#'   proportions. Smaller values concentrate splits on fewer predictors.
#' @param alpha_scale,alpha_shape_1,alpha_shape_2 the prior on `alpha`, in which
#'   `alpha / (alpha + alpha_scale)` is Beta(`alpha_shape_1`,
#'   `alpha_shape_2`). The default `alpha_scale` of `NULL` uses the number of
#'   predictor groups.
#' @param update_s,update_alpha draw the splitting proportions and their
#'   concentration. Turning both off recovers a uniform prior over predictors.
#' @param x_transform how numeric predictors are mapped to `[0, 1]`.
#'   `"quantile"` uses each predictor's empirical distribution function, which
#'   makes the cutpoint prior invariant to monotone reparameterization;
#'   `"range"` rescales linearly, preserving the original spacing.
#' @param verbose report progress while sampling.
#' @param num_print how many iterations between progress reports.
#' @param augment rewrite the likelihood as the margin of a Gaussian one, or as
#'   a Poisson one, which makes the target a shape the sampler can exploit and
#'   the Laplace approximation exact or nearly so. The default, `TRUE`, does it
#'   wherever it has been measured to pay: the binomial family always -- by
#'   Albert and Chib (1993) for a probit link and Polya-Gamma augmentation for a
#'   logit one -- and the negative binomial when the rules are hard. `FALSE`
#'   never does it, and a character vector of engine family names
#'   (`"binomial"`, `"multinomial"`, `"negbin"`) asks for exactly those. Every
#'   augmentation trades speed for mixing, so the measured effect on effective
#'   sample size per second is what matters and it differs by family; see
#'   Details.
#' @param exact_quadratic use the closed forms that a target quadratic in the
#'   additive predictor allows: one pass over a node determines the log target
#'   everywhere, so the Laplace approximation is the conditional posterior rather
#'   than an approximation to it. This is what a conjugate sampler does, and it is
#'   what makes a Gaussian response, or any of the rewritings in `augment`, cheap.
#'   Setting it `FALSE` falls back on the general path and is a diagnostic: the
#'   two agree, at somewhat greater cost.
#' @param generic_accumulate accumulate a leaf's sums through the family's
#'   virtual interface rather than through its own statically dispatched loop.
#'   The two compute the same thing; the second lets the compiler inline the
#'   family's arithmetic, which is most of the remaining per-observation cost.
#'   A diagnostic, for checking that they agree.
#' @param block_eval evaluate the likelihood one leaf at a time rather than one
#'   observation at a time. A family built by [custom_family()] does this
#'   regardless, since it must call back into R; setting it for a compiled
#'   family is a diagnostic, and produces the same draws from the same seed at
#'   somewhat greater cost.
#'
#' @returns A list of class `genbart_control`.
#'
#' @details
#' # Gaussian rewritings
#'
#' The expensive part of this sampler is not the non-conjugacy of the likelihood
#' but that the leaf-level target is not quadratic in the additive predictor.
#' Where a data augmentation makes it quadratic, Fisher scoring reaches the mode
#' in one step, the fitted normal *is* the conditional posterior, and the leaf
#' refresh becomes a Gibbs step with acceptance one.
#'
#' Every such rewriting costs mixing, because the chain now has to move the
#' latent variables and the predictor in alternation. Measured on 500
#' observations, 50 trees and 800 draws, with the ratio that matters being the
#' last column:
#'
#' | Family | Rules | Speed | Effective sample size | ESS per second |
#' |---|---|---|---|---|
#' | `binomial("probit")` | either | 5.7 to 7.2x | 0.66 to 0.75x | **3.8 to 5.3x** |
#' | `binomial("logit")` | either | 2.6 to 3.7x | 0.81 to 1.8x | **2.1 to 6.7x** |
#' | `ordinal("probit")` | soft | 14x | 0.79 to 0.95x | **11 to 13x** |
#' | `ordinal("probit")` | hard | 26 to 30x | 0.73 to 0.90x | **22 to 24x** |
#' | `ordinal("logit")` | soft | 7.0 to 7.2x | 0.70 to 0.82x | **5.0 to 5.9x** |
#' | `ordinal("logit")` | hard | 14.5 to 15.3x | 0.87 to 1.02x | **13 to 15x** |
#' | `ordinal("cloglog")` | soft | 2.8x | 0.57x | 1.6x |
#' | `ordinal("cloglog")` | hard | 5.1x | 1.05x | **5.2x** |
#' | `negbin()` | hard | 1.7 to 1.9x | 0.61 to 1.14x | 1.2 to 2.0x |
#' | `negbin()` | soft | 1.1x | 0.71x | 0.8x |
#' | `multinomial()` | either | 4.2x | 0.37x | 1.6x |
#'
#' The ranges are two problems of different size and shape, which is a fair
#' picture of how much this varies: what an augmentation costs in mixing depends
#' on the data, not only on the family. The negative binomial is the marginal
#' case -- a clear gain on one problem and a slight one on the other -- and worth
#' turning off if its diagnostics look poor.
#'
#' The binomial family uses Polya-Gamma augmentation for a logit link and Albert
#' and Chib's latent normal for a probit one. The negative binomial uses neither:
#' it is written as a Poisson whose rate is drawn from a gamma, which costs one
#' gamma draw per observation and leaves the target in the exponential form the
#' sampler can collapse to a single pass -- but only under hard rules, which is
#' why the gain is there and not under soft ones.
#'
#' `ordinal("probit")` is the largest gain of any of them, because the target it
#' replaces is the most expensive: two cumulative-normal evaluations per
#' observation per pass, and a target with no exploitable shape, so every value
#' of a leaf parameter costs its own pass. Conditional on the latent normal all
#' of that collapses to one pass of arithmetic. The cutpoints are *not* drawn
#' from their conditional given the latent variables, which would be uniform
#' between the two order statistics bracketing each one and would mix worse and
#' worse as the sample grows; they are drawn from the ordinal likelihood with the
#' latent variables integrated out, and the latent variables are redrawn
#' immediately afterwards. That is a partially collapsed Gibbs sampler (Van Dyk
#' and Park 2008) and is the standard remedy (Cowles 1996).
#'
#' `ordinal("logit")` gets the same treatment by a different route, because its
#' latent variable is logistic rather than normal. A logistic variate is a normal
#' whose precision is itself random, and the mixing distribution is usually
#' reached through the Kolmogorov-Smirnov density, which needs a sampler of its
#' own. It does not have to be: Polson, Scott and Windle's (2013) Theorem 1 at
#' `a = 1`, `b = 2` says the standard logistic density is
#' `(1/4) E[exp(-w x^2 / 2)]` with `w` Polya-Gamma(2, 0), so the conditional of
#' the precision given a residual `r` is exactly Polya-Gamma(2, |r|) -- an
#' integer-parameter draw, which the exact Devroye sampler already used for the
#' logistic family covers. Nothing approximate enters. It costs one extra draw
#' per observation per sweep relative to the probit, which is why its gain is
#' about half as large.
#'
#' `ordinal("cloglog")` has a latent variable of a different kind, and it is the
#' one the model is usually named for: the cumulative complementary log-log model
#' is the discrete proportional hazards model, so its survivor function
#' `exp(-exp(c_k - eta))` says exactly that a waiting time with exponential rate
#' `exp(-eta)` has passed `exp(c_k)`. Conditional on that time the log density is
#' `-eta - T * exp(-eta)`, which is the *exponential* form rather than the
#' quadratic one -- the same shape as the gamma family. That is worth 5.1x under
#' hard rules, where the exponential form applies; under soft rules, where it
#' does not, what is left is one `exp()` per observation instead of a difference
#' of two extreme-value distribution functions, which is a smaller gain and costs
#' more mixing.
#'
#' `multinomial()` is a real but modest gain bought with a severe loss of mixing,
#' only worth taking if the chain is lengthened to match, so it is available by
#' name rather than by default.
#'
#' # Soft rules and the cost of a gate
#'
#' A soft rule is what makes the fit smooth, and it is charged for in two places.
#' Every observation reaches more than one leaf, so a pass over a node covers
#' more than the sample -- measured at 2.5 times, at the default bandwidth, for
#' the logistic gate. And the bandwidth is itself a parameter with a Metropolis
#' step per tree per sweep, each of which rebuilds every membership weight in the
#' tree, twice when it is rejected. On a Gaussian response with a thousand
#' observations that one move is about half the total time.
#'
#' A bounded gate addresses the first: past its half-width from the cutpoint the
#' gate is exactly zero or one, so the observation takes one side outright, the
#' other subtree is never visited, and the gate is a polynomial rather than an
#' `exp()`. Measured at 1.4x with no loss of accuracy, on three test functions
#' and six replicates. The second is cheaper than it was -- a rejected proposal is
#' rolled back from a snapshot rather than by evaluating every gate again, and a
#' tree with no splits has no gate at all, so its bandwidth is drawn straight from
#' its prior, which is exactly its full conditional -- but it is still about half
#' of a soft-rule fit.
#'
#' The half-width is `pi * sqrt((2a + 1) / 3)` times `bandwidth` for the
#' Beta(a, a) gate -- 4.06 for `"smoothstep"`, 4.80 for `"smootherstep"` -- which
#' equates the gates' standard deviations, so `bandwidth` means the same amount
#' of smoothing whichever is chosen.
#'
#' *Which* bounded gate makes almost no difference, and not for the reason it
#' looks like it should. At a bandwidth wide enough that a bounded gate truncates
#' nothing at all it is still 1.45 times faster than the logistic: what the
#' bounded gates save is the `exp()`, not the work on the far side of the
#' cutpoint. So the two come out within noise of each other on time and on
#' accuracy, and the reason to prefer `"smootherstep"` is that it gives a
#' twice-differentiable fit.
#'
#' `bandwidth_every` addresses the second, and is a real trade rather than a free
#' one. Fixing the bandwidth entirely (`update_bandwidth = FALSE`) is 2.4x faster
#' again and *more* accurate on smooth functions -- but much worse on functions
#' with jumps, where what the update is for is letting the rules sharpen towards
#' hard ones. On a three-step function, measured RMSE was 0.42 fixed against 0.19
#' drawn. So the default draws it every sweep, and a larger `bandwidth_every`
#' buys time at the cost of how fast that adaptation happens.
#'
#' There is nothing to rewrite for the Poisson and gamma families -- their
#' targets are already in the exponential form -- and no known rewriting for the
#' accelerated failure time, ordered beta or location-scale families.
#'
#' @references
#' Albert, J. H., & Chib, S. (1993). Bayesian analysis of binary and
#' polychotomous response data. *Journal of the American Statistical
#' Association*, 88(422), 669--679.
#'
#' Cowles, M. K. (1996). Accelerating Monte Carlo Markov chain convergence for
#' cumulative-link generalized linear models. *Statistics and Computing*, 6(2),
#' 101--111.
#'
#' Linero, A. R. (2025). Generalized Bayesian additive regression trees models:
#' beyond conditional conjugacy. *Journal of the American Statistical
#' Association*, 120(549), 356--369.
#'
#' Linero, A. R., & Yang, Y. (2018). Bayesian regression tree ensembles that
#' adapt to smoothness and sparsity. *Journal of the Royal Statistical Society
#' Series B*, 80(5), 1087--1110.
#'
#' Van Dyk, D. A., & Park, T. (2008). Partially collapsed Gibbs samplers: theory
#' and methods. *Journal of the American Statistical Association*, 103(482),
#' 790--796.
#'
#' @seealso [genbart()]
#'
#' @examples
#' genbart_control(num_trees = 20, soft = FALSE)
#'
#' @export
genbart_control <- function(num_trees = 50L, num_burn = 500L, num_save = 500L,
                            num_thin = 1L, soft = TRUE, bandwidth = 0.1,
                            update_bandwidth = TRUE, bandwidth_every = 1L,
                            gate = "logistic",
                            k = 2, sigma_mu = NULL,
                            update_sigma_mu = TRUE, update_tau = TRUE,
                            sigma_mu_ramp = 0.25,
                            gamma = 0.95, beta = 2, alpha = 1,
                            alpha_scale = NULL, alpha_shape_1 = 0.5,
                            alpha_shape_2 = 1, update_s = TRUE,
                            update_alpha = TRUE, x_transform = "quantile",
                            verbose = FALSE, num_print = 100L,
                            augment = TRUE, block_eval = FALSE,
                            exact_quadratic = TRUE,
                            generic_accumulate = FALSE) {

  # Record what the caller actually supplied, before any of the arguments below
  # are normalized. `genbart()` uses this to rebuild the control list when extra
  # settings are passed through its `...`, the way `glm()` does.
  supplied <- mget(as.character(setdiff(names(match.call())[-1], "")),
                   environment())

  for (nm in c("num_trees", "num_burn", "num_save", "num_thin", "num_print",
               "bandwidth_every")) {
    value <- get(nm)
    arg::arg_whole_number(value, .arg = nm)
    arg::arg_gte(value, 0, .arg = nm)
  }

  arg::arg_gte(num_trees, 1)
  arg::arg_gte(bandwidth_every, 1)
  arg::arg_gte(num_save, 1)
  arg::arg_gte(num_thin, 1)

  for (nm in c("soft", "update_bandwidth", "update_sigma_mu", "update_tau",
               "update_s",
               "update_alpha", "verbose", "block_eval",
               "exact_quadratic", "generic_accumulate")) {
    value <- get(nm)
    arg::arg_flag(value, .arg = nm)
  }

  for (nm in c("bandwidth", "k", "gamma", "beta", "alpha", "alpha_shape_1",
               "alpha_shape_2")) {
    value <- get(nm)
    arg::arg_number(value, .arg = nm)
    arg::arg_gt(value, 0, .arg = nm)
  }

  arg::arg_lte(gamma, 1)
  arg::arg_number(sigma_mu_ramp)
  arg::arg_between(sigma_mu_ramp, c(0, 1))

  if (!is_null(sigma_mu)) {
    arg::arg_numeric(sigma_mu)

    if (any(sigma_mu <= 0)) {
      arg::err("{.arg sigma_mu} must be positive")
    }
  }

  if (!is_null(alpha_scale)) {
    arg::arg_number(alpha_scale)
    arg::arg_gt(alpha_scale, 0)
  }

  x_transform <- arg::match_arg(x_transform, c("quantile", "range"))
  gate <- arg::match_arg(gate, c("logistic", "smoothstep",
                                 "smootherstep"))
  augment <- resolve_augment(augment, soft)

  out <- list(num_trees = as.integer(num_trees),
              num_burn = as.integer(num_burn),
              num_save = as.integer(num_save),
              num_thin = as.integer(num_thin),
              soft = soft,
              bandwidth = bandwidth,
              update_bandwidth = update_bandwidth,
              bandwidth_every = as.integer(bandwidth_every),
              gate = gate,
              k = k,
              sigma_mu = sigma_mu,
              update_sigma_mu = update_sigma_mu,
              update_tau = update_tau,
              sigma_mu_ramp = sigma_mu_ramp,
              gamma = gamma,
              beta = beta,
              alpha = alpha,
              alpha_scale = alpha_scale %or% 0,
              alpha_shape_1 = alpha_shape_1,
              alpha_shape_2 = alpha_shape_2,
              update_s = update_s,
              update_alpha = update_alpha,
              x_transform = x_transform,
              verbose = verbose,
              num_print = as.integer(num_print),
              augment = augment,
              block_eval = block_eval,
              exact_quadratic = exact_quadratic,
              generic_accumulate = generic_accumulate)

  class(out) <- "genbart_control"
  attr(out, "supplied") <- supplied

  out
}

# The engine family names whose likelihood should be rewritten as the margin of a
# Gaussian one. TRUE means the ones where that has been measured to pay, which
# for the negative binomial depends on the decision rules: its rewriting is a
# gain of 2x in effective samples per second with hard rules, where the target's
# exponential form collapses the leaf work to one pass, and a slight loss with
# soft rules, where it does not. Naming a family explicitly always honors the
# request.
resolve_augment <- function(augment, soft) {
  if (is.logical(augment)) {
    arg::arg_flag(augment)

    if (!isTRUE(augment)) {
      return(character(0L))
    }

    if (isTRUE(soft)) {
      return(c("binomial", "ordinal"))
    }

    return(c("binomial", "ordinal", "negbin"))
  }

  arg::arg_character(augment)
  arg::arg_element(augment, c("binomial", "ordinal", "multinomial", "negbin"),
                   .arg = "augment")

  unique(augment)
}
