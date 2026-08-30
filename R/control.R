#' Sampler and prior settings for `bartisan()`
#'
#' @description
#' Collects the tuning constants of the sampler and the hyperparameters of the
#' tree prior. The defaults follow Linero (2025) and, for the soft decision
#' rules, Linero and Yang (2018), and are intended to be usable without
#' adjustment.
#'
#' The arguments fall into three groups, and the first group is the one worth
#' reading: `num_trees`, `gate`, `sparsity`, `k`, `bandwidth`, the three chain
#' lengths and `x_transform` are **modeling decisions**, in that changing one
#' changes what is being fitted or how long it is fitted for. Everything from
#' `augment` to `num_print` is an **advanced setting**: a hyperparameter of a
#' prior the first group summarizes, or a switch whose default is almost always
#' right. The last three -- `block_eval`, `exact_quadratic` and
#' `generic_accumulate` -- exist **for internal validation** and are documented
#' so that the checks that use them can be read; they compute the same posterior
#' more slowly.
#'
#' @param num_trees number of trees, either one number for every forest or one
#'   per additive predictor. The default, `NULL`, is 50 for both kinds of rule.
#'   The two do not want the same count -- a soft rule makes each tree more
#'   expressive, so soft rules reach their best held-out error at about 20 trees
#'   and get worse past that, while hard rules keep improving to 200 -- but a
#'   smaller forest mixes worse, which is why 50 is the compromise for both. A
#'   family with more than one additive predictor takes a vector, which is worth
#'   using: `location_scale()` spends about 90% of its time on the scale forest,
#'   and a variance surface needs less capacity than a mean surface, so
#'   `num_trees = c(50, 10)` runs about 2.5 times faster at the same accuracy.
#'   All of this is measured in Details.
#' @param gate the shape of a decision rule, which is also how hard and soft
#'   rules are chosen between. `"hard"` (or `"step"`) gives the step functions of
#'   standard BART. The other three give soft rules, as in the SoftBart model, in
#'   which every observation reaches every leaf with some weight and the fitted
#'   function is smooth: `"smoothstep"`, the default, is the Beta(2, 2)
#'   cumulative distribution function, once differentiable and supported on a
#'   bounded interval; `"smootherstep"` is Beta(3, 3), twice differentiable and
#'   bounded; `"logistic"` is the logistic function of Linero and Yang (2018) and
#'   is infinitely differentiable. Soft rules cost three to five times as much
#'   per iteration and cut held-out error by 35 to 40%, so they are the accuracy
#'   argument rather than a tax. The two bounded gates are about 1.4 times faster
#'   than the logistic and equally accurate, and are within noise of each other;
#'   see Details.
#' @param sparsity the prior on which predictors get split on. `TRUE`, the
#'   default, is the Dirichlet sparsity prior of Linero (2018) -- DART -- which
#'   concentrates splits on the predictors that earn them and can drop the rest
#'   from the forest entirely. `FALSE` gives every predictor the same splitting
#'   probability, which is classic BART. `"none"`, `"weak"`, `"moderate"` and
#'   `"strong"` name four strengths, with `"none"` equal to `FALSE`, `"moderate"`
#'   equal to `TRUE`, and the other two moving the prior on the concentration.
#'   This argument sets `update_s`, `update_alpha`, `alpha_shape_1` and
#'   `alpha_shape_2` together; supplying any of those directly overrides it. Read
#'   the trade-off in Details before turning it off or up.
#' @param k controls the leaf prior. The prior standard deviation of a forest is
#'   `3 / k` times the natural scale of its additive predictor, so larger `k`
#'   shrinks the fit harder towards the intercept-only model.
#' @param bandwidth prior mean of the gate bandwidth of a soft rule, on the
#'   scale of the transformed predictors, which lie in `[0, 1]`. Smaller values
#'   approach hard rules. Ignored for `gate = "hard"`.
#' @param num_burn number of warmup iterations to discard.
#' @param num_save number of draws to keep.
#' @param num_thin keep one draw in every `num_thin` after warmup.
#' @param x_transform how numeric predictors are mapped to `[0, 1]`.
#'   `"quantile"` uses each predictor's empirical distribution function, which
#'   makes the cutpoint prior invariant to monotone reparameterization;
#'   `"range"` rescales linearly, preserving the original spacing.
#' @param augment *Advanced.* Rewrite the likelihood as the margin of a Gaussian
#'   one, or of a Poisson one, which makes the target a shape the sampler can
#'   exploit and the Laplace approximation exact or nearly so. **The same
#'   posterior either way**, so this is a sampling setting rather than a modeling
#'   one. The default, `TRUE`, does it wherever it has been measured to pay: the
#'   binomial, ordinal, multinomial, zero-inflated and survival families always,
#'   and the negative binomial when the rules are hard. `FALSE` never does it,
#'   and a character vector of engine family names -- `"binomial"`, `"ordinal"`,
#'   `"multinomial"`, `"negbin"`, `"zip"`, `"zinb"`, `"aft"` -- asks for exactly
#'   those.
#'   Every rewriting trades speed for mixing, so the measured effect on effective
#'   sample size per second is what matters and it differs by family; see Details.
#' @param gamma,beta *Advanced.* The branching probability at depth `d` is
#'   `gamma * (1 + d)^(-beta)`.
#' @param sigma_mu *Advanced.* Prior median of the leaf standard deviation, one
#'   value per additive predictor. The default, `NULL`, derives it from `k` and
#'   that forest's own tree count.
#' @param update_sigma_mu *Advanced.* Draw the leaf standard deviation under a
#'   half-Cauchy prior rather than fixing it.
#' @param sigma_mu_ramp *Advanced.* Fraction of warmup over which the leaf
#'   standard deviation is raised from near zero to its target. Linero (2025)
#'   describes this as essential: started at its full value, the sampler can
#'   settle early into a poor configuration and fail to move. Set to `0` to
#'   disable.
#' @param update_tau *Advanced.* Draw the standard deviation of each
#'   random-effect term under the same half-Cauchy prior the leaf scale uses,
#'   rather than fixing it at that prior's median. Only relevant when the formula
#'   has a `(1 | group)` term.
#' @param update_bandwidth *Advanced.* Draw the bandwidth of each tree rather
#'   than holding it at `bandwidth`.
#' @param bandwidth_every *Advanced.* How many sweeps between bandwidth draws
#'   for a given tree. The bandwidth is one scalar per tree, drawn by an adaptive
#'   random walk, and every attempt costs a full rebuild of the tree's
#'   memberships -- the single largest item in a soft-rule fit. Drawing it less
#'   often than the trees themselves trades mixing in that one parameter for
#'   time; see Details.
#' @param alpha *Advanced.* Concentration of the Dirichlet prior on the
#'   splitting proportions. Smaller values concentrate splits on fewer
#'   predictors. The default, `NULL`, uses 1, which is the starting value for a
#'   parameter that is then drawn.
#' @param alpha_scale,alpha_shape_1,alpha_shape_2 *Advanced.* The prior on
#'   `alpha`, in which `alpha / (alpha + alpha_scale)` is
#'   Beta(`alpha_shape_1`, `alpha_shape_2`). The default `alpha_scale` of `NULL`
#'   uses the number of predictor groups; the two shapes default to whatever
#'   `sparsity` implies.
#' @param update_s,update_alpha *Advanced.* Draw the splitting proportions and
#'   their concentration. The defaults, `NULL`, follow `sparsity`. Turning both
#'   off recovers a uniform prior over predictors, which is what
#'   `sparsity = FALSE` does.
#' @param verbose *Advanced.* Report progress while sampling.
#' @param num_print *Advanced.* How many iterations between progress reports.
#' @param block_eval *Validation.* Evaluate the likelihood one leaf at a time
#'   rather than one observation at a time. A family built by [custom_family()]
#'   does this regardless, since it must call back into R; setting it for a
#'   compiled family produces the same draws from the same seed at somewhat
#'   greater cost.
#' @param exact_quadratic *Validation.* Use the closed forms that a target
#'   quadratic in the additive predictor allows: one pass over a node determines
#'   the log target everywhere, so the Laplace approximation is the conditional
#'   posterior rather than an approximation to it. This is what a conjugate
#'   sampler does, and it is what makes a Gaussian response, or any of the
#'   rewritings in `augment`, cheap. Setting it `FALSE` falls back on the general
#'   path; the two agree, at greater cost.
#' @param generic_accumulate *Validation.* Accumulate a leaf's sums through the
#'   family's virtual interface rather than through its own statically dispatched
#'   loop. The two compute the same thing; the second lets the compiler inline
#'   the family's arithmetic, which is most of the remaining per-observation
#'   cost.
#'
#' @returns A list of class `bartisan_control`.
#'
#' @details
#' # How many trees, and how many per forest
#'
#' The tree count is the setting most worth thinking about after the family, and
#' the default of 50 is a compromise rather than an optimum. Measured on the
#' Friedman function with 1000 training and 1000 test observations, four chains
#' of 500 draws after 500 warmup iterations, as held-out root mean squared error
#' against the true regression function:
#'
#' | Trees | Soft rules | Hard rules |
#' |---|---|---|
#' | 5 | 0.286 | 1.149 |
#' | 10 | 0.281 | 0.682 |
#' | 20 | **0.270** | 0.558 |
#' | 50 | 0.284 | **0.521** |
#' | 100 | 0.289 | 0.531 |
#' | 200 | 0.319 | 0.510 |
#'
#' Two things to read off it. **Soft rules need far fewer trees than hard ones**,
#' which is what makes 200 -- the default in most BART packages -- actively worse
#' here than 20. And **the two kinds of rule want different counts**, since hard
#' rules are still improving at 200 where soft rules peaked at 20.
#'
#' The default is 50 for both anyway, for a reason the table cannot show: a
#' smaller forest mixes worse. On `MatchIt::lalonde`, four chains at 20 soft
#' trees disagreed by 35% on an average contrast where 50 trees disagreed by 9%,
#' and the Friedman gain at 20 trees was 5%. So 50 buys reliable inference at a
#' small cost in point accuracy. Drop to 20 if prediction is the only goal and
#' the fit is soft; raise towards 200 with hard rules if it is not.
#'
#' **A vector is worth using when the family has more than one forest.** They are
#' not equally expensive and they do not need equal capacity. `location_scale()`
#' spends about 90% of its time on the scale forest, because that forest's target
#' is not quadratic, and a variance surface carries much less information than a
#' mean surface. Measured on 1000 observations with a smooth mean and a
#' log-linear standard deviation:
#'
#' | `num_trees` | Seconds | Mean RMSE | Log-SD RMSE | Log score |
#' |---|---|---|---|---|
#' | `c(50, 50)` | 14.5 | 0.092 | 0.050 | -1188 |
#' | `c(50, 20)` | 8.3 | 0.094 | 0.047 | -1188 |
#' | `c(50, 10)` | 5.9 | 0.093 | 0.051 | -1188 |
#' | `c(50, 5)` | **4.9** | 0.094 | 0.046 | -1187 |
#' | `c(20, 5)` | **2.9** | 0.084 | 0.041 | -1184 |
#'
#' A Gaussian fit on the same data takes 1.4 seconds, so `c(50, 50)` costs 10
#' times a Gaussian fit and `c(50, 5)` costs 3.5 times, at the same accuracy to
#' three decimal places. That is not the default, because how many trees a
#' variance surface needs depends on how complicated it is, and silently
#' under-parameterizing it would show up as intervals that are wrong -- which is
#' the thing `location_scale()` exists to get right. It is worth setting by hand.
#'
#' The leaf prior scale divides by the square root of each forest's own tree
#' count, so shrinking one forest does not change the prior on the sum it forms.
#'
#' # The sparsity prior, and what it costs
#'
#' `sparsity = TRUE` is the Dirichlet prior of Linero (2018) on the splitting
#' proportions. It is a genuine variable-selection prior: it can and does drop a
#' predictor from every tree in the forest at once, and that is the point of it in
#' the high-dimensional problems it was built for.
#'
#' It has a consequence that is easy to misread. **A contrast on a predictor the
#' prior has dropped is exactly zero, not nearly zero**, because in that draw the
#' fit does not depend on that predictor at all. So the posterior of a contrast
#' has an atom at zero whose mass is one minus the predictor's inclusion
#' probability, which `summary()` reports as `prop_used`. Any summary that
#' reports a median lands on that atom whenever it holds half the mass; see
#' [bartisan-marginaleffects].
#'
#' **More trees does not fix it**, which is worth stating because it looks as
#' though it should. On `MatchIt::lalonde`, four chains with the sparsity prior on
#' left `treat` out of every tree in 20% of draws at 50 trees and 18% at 200,
#' and disagreed with each other by more than 100% of the estimate at every count
#' from 20 trees to 200. With `sparsity = FALSE` the atom disappears by 50 trees
#' and the four chains agreed to within 9%. The variable-selection state mixes
#' slowly: a predictor whose splitting proportion has gone small is rarely
#' proposed, so it is hard to get back in.
#'
#' So: keep the default when prediction is the goal or the predictors are many
#' and mostly irrelevant, and use `sparsity = FALSE` when a contrast or a partial
#' effect on a particular predictor is the estimand. Run several chains either
#' way, because one chain can look far more settled than the posterior is.
#'
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
#' | `multinomial()` | soft | 9.3x | 1.09x | **10.1x** |
#' | `multinomial()` | hard | 14.5x | 0.66x | **9.6x** |
#' | `zi_poisson()` | soft | 4.6x | 0.85x | **3.9x** |
#' | `zi_poisson()` | hard | 7.1x | 1.42x | **10.1x** |
#' | `zi_negbin()` | soft | 5.9x | 0.98x | **5.8x** |
#' | `zi_negbin()` | hard | 9.4x | 0.84x | **7.9x** |
#' | `lognormal_aft()` | soft | 11.2 to 19.7x | 0.83 to 1.27x | **14 to 16x** |
#' | `lognormal_aft()` | hard | 19.7 to 29.3x | 0.74 to 1.01x | **20 to 22x** |
#' | `loglogistic_aft()` | soft | 8.4 to 9.3x | 0.58 to 0.84x | **4.9 to 7.7x** |
#' | `loglogistic_aft()` | hard | 11.3 to 12.3x | 0.87 to 0.89x | **9.8 to 11x** |
#'
#' The ranges are two problems of different size and shape, which is a fair
#' picture of how much this varies: what an augmentation costs in mixing depends
#' on the data, not only on the family. The negative binomial is the marginal
#' case -- a clear gain on one problem and a slight one on the other -- and worth
#' turning off if its diagnostics look poor.
#'
#' The two survival rows are the other end of the range from the negative
#' binomial, and for a reason worth stating: right-censoring is what makes the
#' direct likelihood expensive. An observed failure contributes a density in the
#' predictor and a censored observation contributes a survival function, and the
#' two have different shapes, so the target has no exploitable form and every
#' trial value of a leaf parameter costs its own pass over the data. Imputing
#' each censored failure time above its censoring time replaces the survival term
#' with a density, and then every observation contributes the same quadratic
#' shape. `weibull_aft()` is absent from the table because it needs none of this:
#' its likelihood is already of the exponential form, censoring included, so it
#' gets the single-pass treatment directly under hard rules.
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
#' # Rare events mix slowly, and it is the data rather than the augmentation
#'
#' With very few events the whole fit mixes slowly, augmented or not. Measured at
#' 2000 observations, 50 trees and 1000 draws, the effective sample size of the
#' level of the predictor falls from about 860 at half positives to 75 at 1.4%
#' and to 8 at 0.2%. That is worth knowing about but it is not a fault of the
#' rewritings: the *shape* of the fit mixes at the same rate as its level, moving
#' the anchor by three units on the probit scale does not change the figure, and
#' \pkg{dbarts} -- an independent implementation of the same latent normal, with
#' the same anchoring -- reproduces it to three digits. It is the information in a
#' handful of events. Lengthen the chain, and read the `rhat` element of the fit
#' rather than assuming the default draw count is enough.
#'
#' `multinomial()` uses the same Polya-Gamma identity as the binomial logit, and
#' can, because conditional on the other categories the likelihood of category
#' *j* is *exactly* binomial-logistic in `eta_j - log C_j`, where `C_j` sums
#' `exp(eta)` over the others. So no stick-breaking decomposition is needed: the
#' structure the sampler already has, one forest at a time, is the structure the
#' augmentation wants. The latent variables are redrawn per forest rather than
#' per sweep, since `C_j` moves during a sweep.
#'
#' The **zero-inflated** families need two latent variables, because what blocks
#' them is a mixture rather than a link. The zero contributes
#' `log[pi + (1 - pi) P_0]`, a log-sum-exp of the two components, so neither
#' predictor has a shape. Introducing the indicator the mixture is a mixture over
#' -- whether the observation is a structural zero -- separates them: conditional
#' on it the count forest sees a plain Poisson or negative binomial, and the
#' inflation forest sees a Bernoulli logistic likelihood, which Polya-Gamma
#' handles. The indicator is drawn from its exact conditional, which is zero
#' whenever the count is positive. For `zi_negbin()` the Poisson-gamma rate goes
#' on top of that, so the count forest gets the exponential form as well; the
#' indicator is drawn with the rate integrated out and the rate redrawn
#' afterwards, which is a valid partially collapsed Gibbs step in that order.
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
#' Linero, A. R. (2018). Bayesian regression trees for high-dimensional
#' prediction and variable selection. *Journal of the American Statistical
#' Association*, 113(522), 626--636. \doi{10.1080/01621459.2016.1264957}
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
#' @seealso [bartisan()]
#'
#' @examples
#' bartisan_control(num_trees = 20, gate = "hard")
#'
#' @export
bartisan_control <- function(num_trees = NULL,
                             gate = "smoothstep",
                            sparsity = TRUE,
                            k = 2,
                            bandwidth = 0.1,
                            num_burn = 500L, num_save = 500L, num_thin = 1L,
                            augment = TRUE,
                            x_transform = "quantile",
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

  for (nm in c("num_burn", "num_save", "num_thin", "num_print",
               "bandwidth_every")) {
    value <- get(nm)
    arg::arg_whole_number(value, .arg = nm)
    arg::arg_gte(value, 0, .arg = nm)
  }

  # One tree count per forest, or one for all of them. Resolving how many
  # forests there are needs the family, so that is left to `bartisan()`; all that
  # can be checked here is that the values themselves are usable.
  arg::when_not_null(
    num_trees,
    arg::arg_and(
      arg::arg_whole_numeric,
      arg::arg_gte(1)
    )
  )

  arg::arg_gte(bandwidth_every, 1)
  arg::arg_gte(num_save, 1)
  arg::arg_gte(num_thin, 1)

  for (nm in c("update_bandwidth", "update_sigma_mu", "update_tau",
               "verbose", "block_eval",
               "exact_quadratic", "generic_accumulate")) {
    value <- get(nm)
    arg::arg_flag(value, .arg = nm)
  }

  for (nm in c("update_s", "update_alpha")) {
    arg::when_not_null(get(nm), arg::arg_flag, .arg = nm)
  }

  for (nm in c("bandwidth", "k", "gamma", "beta")) {
    value <- get(nm)
    arg::arg_number(value, .arg = nm)
    arg::arg_gt(value, 0, .arg = nm)
  }

  for (nm in c("alpha", "alpha_scale", "alpha_shape_1", "alpha_shape_2")) {
    arg::when_not_null(
      get(nm),
      arg::arg_and(
        arg::arg_number,
        arg::arg_gt(0)
      ),
      .arg = nm
    )
  }

  arg::arg_lte(gamma, 1)
  arg::arg_number(sigma_mu_ramp)
  arg::arg_between(sigma_mu_ramp, c(0, 1))

  arg::when_not_null(
    sigma_mu,
    arg::arg_and(
      arg::arg_numeric,
      arg::arg_gt(0)
    )
  )

  x_transform <- arg::match_arg(x_transform, c("quantile", "range"))

  # One argument for both the shape of a soft rule's gate and the choice between
  # soft and hard rules, because they are one decision: a hard rule is the
  # limiting case of a soft one and there is no gate shape to pick for it.
  gate <- arg::match_arg(gate, c("smoothstep", "smootherstep", "logistic",
                                 "hard", "step"))
  soft <- !gate %in% c("hard", "step")

  sparse <- resolve_sparsity(sparsity)
  augment <- resolve_augment(augment, soft)

  out <- list(num_trees = if (is_null(num_trees)) NULL else as.integer(num_trees),
              gate = gate,
              soft = soft,
              sparsity = sparsity,
              k = k,
              bandwidth = bandwidth,
              num_burn = as.integer(num_burn),
              num_save = as.integer(num_save),
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
# harder towards a few predictors.
resolve_sparsity <- function(sparsity) {
  arg::arg_or(sparsity,
              arg::arg_flag,
              arg::arg_string)

  level <- {
    if (isTRUE(sparsity)) "moderate"
    else if (isFALSE(sparsity)) "none"
    else arg::match_arg(sparsity, c("none", "weak", "moderate", "strong"),
                        .arg = "sparsity")
  }

  shapes <- switch(level,
                   none = c(0.5, 1),
                   weak = c(1, 1),
                   moderate = c(0.5, 1),
                   strong = c(0.5, 3))

  list(update_s = !identical(level, "none"),
       update_alpha = !identical(level, "none"),
       alpha = 1,
       alpha_shape_1 = shapes[1L],
       alpha_shape_2 = shapes[2L])
}

# The engine family names whose likelihood should be rewritten as the margin of a
# Gaussian one. TRUE means the ones where that has been measured to pay, which
# for the negative binomial depends on the decision rules: its rewriting is a
# gain of 2x in effective samples per second with hard rules, where the target's
# exponential form collapses the leaf work to one pass, and a slight loss with
# soft rules, where it does not. Naming a family explicitly always honors the
# request.
resolve_augment <- function(augment, soft) {
  arg::arg_or(augment,
              arg::arg_flag,
              arg::arg_character)

  if (is.character(augment)) {
    arg::arg_element(augment,
                     c("binomial", "ordinal", "multinomial", "negbin", "zip",
                       "zinb", "aft"),
                     .arg = "augment")

    return(unique(augment))
  }

  if (!isTRUE(augment)) {
    return(character())
  }

  # The negative binomial only where it pays: with hard rules, where the target's
  # exponential form applies. Every other rewriting pays under both kinds of
  # rule.
  c("binomial", "ordinal", "multinomial", "zip", "zinb", "aft",
    "negbin"[!isTRUE(soft)])
}
