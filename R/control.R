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
#' @param categorical how a splitting rule divides the levels of a factor.
#'   `"subset"`, the default, draws a subset of the levels still available at the
#'   node and sends those left, which is the rule of Deshpande (2024).
#'   `"onehot"` is what most BART implementations do: it splits on one indicator
#'   column, which peels a single level off the rest. The choice matters because
#'   `"onehot"` reaches only `2^K - K` of the `B_K` partitions of `K` levels --
#'   27 of 52 at `K = 5`, and 1,014 of 115,975 at `K = 10` -- and the partitions
#'   it can form all have at most one cell with more than one level in it, so the
#'   bulk of the levels is never divided. See Details.
#' @param split_prior relative prior weight on each predictor, for when you
#'   expect some to matter more than others. A named numeric vector, keyed by the
#'   names the predictors have in the formula; every predictor not named gets a
#'   weight of 1. The prior probability of splitting on a predictor is its weight
#'   divided by the total, so on a three-predictor model
#'   `split_prior = c(x1 = 3, x3 = 0.5)` gives `x1` a probability of `3 / 4.5`,
#'   `x2` `1 / 4.5` and `x3` `0.5 / 4.5`. Weights must be finite and not
#'   than negative, and naming a predictor the model does not have is an error
#'   rather than silently ignored. A weight of zero is allowed and means the
#'   predictor is never split on: it stays in the model frame and out of every
#'   tree. The default, `NULL`, weights every predictor equally.
#'   Setting this overrides `sparsity`; see Details.
#' @param k controls the leaf prior. The prior standard deviation of a forest is
#'   `3 / k` times the natural scale of its additive predictor, so larger `k`
#'   shrinks the fit harder towards the intercept-only model.
#' @param bandwidth prior mean of the gate bandwidth of a soft rule, on the
#'   scale of the transformed predictors, which lie in `[0, 1]`. Smaller values
#'   approach hard rules. Ignored for `gate = "hard"`.
#' @param chains how many independent chains to run. The draws are pooled and
#'   [split-R-hat][bartisan()] is reported in the `rhat` element. With the
#'   \pkg{future.apply} package installed the chains run in parallel under
#'   whatever backend the caller has planned with \pkgfun{future}{plan} --
#'   `multisession`, `multicore`, a cluster, or \pkg{mirai}'s
#'   `mirai_multisession`; without it they run one after another, which is slower
#'   and otherwise identical. One `set.seed()` before the call reproduces the
#'   whole run either way, because each chain is given its own L'Ecuyer stream.
#' @param num_burn number of warmup iterations to discard. Warmup is where the
#'   trees grow into the data and the hyperparameters find their scale, so raising
#'   it buys convergence rather than precision: increase it when `rhat` says the
#'   chains have not agreed.
#' @param num_draws number of draws to keep. These are what every estimate and
#'   interval is computed from, so raising it narrows Monte Carlo error and does
#'   nothing about convergence: increase it when `ess_bulk` or `ess_tail` is small
#'   relative to what the reported quantity needs.
#' @param num_thin keep one draw in every `num_thin` after warmup. Thinning
#'   discards draws to make the kept ones less correlated, which costs
#'   information and is worth it only to hold down the memory a long chain would
#'   otherwise take: for a given amount of computing, more draws beat fewer
#'   less-correlated ones.
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
#' @param verbose *Advanced.* Print progress to the console while sampling. For
#'   a progress bar instead, see the Progress section below, which needs no
#'   argument here.
#' @param num_print *Advanced.* How many iterations between the reports
#'   `verbose` prints.
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
#' **What each setting is worth, measured.** The trade is real in both
#' directions, and the strength barely matters in either: what matters is whether
#' the prior is on at all.
#'
#' For prediction, any sparsity beats none. Scored against the true regression
#' function on held-out data, Friedman with five relevant predictors:
#'
#' | predictors | `"none"` | `"weak"` | `"moderate"` | `"strong"` |
#' | --- | --- | --- | --- | --- |
#' | 10 | 0.446 | 0.385 | 0.400 | 0.374 |
#' | 50 | 0.465 | 0.346 | 0.372 | 0.362 |
#'
#' For a contrast on a predictor whose signal is weak, any sparsity is actively
#' harmful, and not only in the atom-at-zero sense above. A binary treatment
#' among 20 predictors, continuous outcome, residual standard deviation 1,
#' n = 800, five replicates, with `covers` the share of replicates whose 95%
#' interval contains the truth:
#'
#' | true effect | setting | estimate | atom | covers |
#' | --- | --- | --- | --- | --- |
#' | 0.05 | `FALSE` | 0.031 | 0.08 | 0.80 |
#' | 0.05 | `TRUE` | 0.000 | 0.89 | 0.40 |
#' | 0.10 | `FALSE` | 0.131 | 0.03 | 1.00 |
#' | 0.10 | `TRUE` | 0.029 | 0.69 | 1.00 |
#' | 0.20 | `FALSE` | 0.161 | 0.05 | 1.00 |
#' | 0.20 | `TRUE` | 0.094 | 0.55 | 0.60 |
#' | 0.50 | `FALSE` | 0.475 | 0.00 | 1.00 |
#' | 0.50 | `TRUE` | 0.474 | 0.00 | 1.00 |
#'
#' The prior attenuates a weak effect by half or more and its interval covers
#' well below its nominal rate. A strong effect is untouched, because the prior
#' never has reason to drop a predictor that is earning its splits, so this is a
#' weak-signal failure rather than a general one.
#'
#' So the setting is close to a switch, and which way to throw it follows from
#' the estimand rather than from the data:
#'
#' - **Prediction, or variable selection**: keep the default. `"weak"` is already
#'   worth most of what `"strong"` is worth, so reach past `TRUE` only when the
#'   predictors are many and you expect nearly all of them to be irrelevant.
#' - **A contrast, a partial effect, or a treatment effect**: `sparsity = FALSE`,
#'   or `split_prior`, which fixes the weights and so cannot drop anything. Use
#'   `split_prior` in preference when the other predictors are numerous enough
#'   that weighting them all alike is wasteful.
#'
#' Run several chains either way, because one chain can look far more settled
#' than the posterior is.
#'
#' # Arguments that vary by forest
#'
#' A family with several additive predictors has one forest per predictor, each
#' with its own prior. Every argument that could mean something different for one
#' of them may be given once, to apply to all, or one per forest -- positionally,
#' or keyed by the forest names listed in [bartisan-families]. A forest that a
#' named argument does not mention keeps that argument's default rather than
#' borrowing the value chosen for another forest.
#'
#' That covers `num_trees`, `k`, `sigma_mu`, `sparsity`, `split_prior`,
#' `bandwidth`, `gamma`, `beta`, the four `alpha` arguments, and the three
#' `update_` flags for the leaf scale, the splitting proportions and the
#' bandwidth. `formula` works the same way; see [bartisan()].
#'
#' ```r
#' # A scale forest with less capacity than the mean forest, and no sparsity
#' # prior on it.
#' bartisan_control(num_trees = c(mean = 50, log_sd = 10),
#'                  sparsity = c(mean = TRUE, log_sd = FALSE))
#' ```
#'
#' The multinomial families are the exception, for the reason given in
#' [bartisan-families]: their forests act as one, so these arguments take a
#' single value.
#'
#' # Splitting a factor
#'
#' A rule on one indicator column of a factor can only separate one level from
#' the rest. Applied repeatedly down a path that produces a partition with some
#' number of singleton levels and one cell holding everything else, and those are
#' the only partitions available: `2^K - K` of the `B_K` partitions of `K`
#' levels, which is 27 of 52 at `K = 5` and under 1% at `K = 10`.
#'
#' What that costs is partial pooling. Simulating the tree prior directly at
#' `K = 10`, the probability that two levels land in the same leaf is 0.77 under
#' `"onehot"` against 0.46 under `"subset"`, and a typical `"onehot"` tree has
#' 1.18 singleton levels out of 2.18 leaves: one level alone, the rest together,
#' whether the data want that or not. A rule that takes a subset of the levels
#' reaches every partition, and the co-clustering probabilities come out far more
#' even.
#'
#' **What it is worth depends on the decision rules.** Twenty levels in four
#' groups of five sharing a mean, twelve replicates, paired within replicate,
#' RMSE against the true mean function:
#'
#' | per level | subset, hard | onehot, hard | subset, soft | onehot, soft |
#' | --- | --- | --- | --- | --- |
#' | 10 | 0.3705 | 0.3998 | 0.3435 | 0.3445 |
#' | 25 | 0.2617 | 0.2725 | 0.2265 | 0.2201 |
#' | 100 | 0.1421 | 0.1488 | 0.1098 | 0.1059 |
#'
#' Under hard rules `"subset"` is better at every size, clearly so at ten
#' observations per level, where the gap is 7% against a standard error of 2%,
#' and by 4% at the two larger sizes, where the standard error is around half
#' the gap. Under soft rules, which is the default, **the two are
#' indistinguishable**: the largest gap is 0.006 against a standard error of
#' 0.005. Soft rules are worth far more than either choice, which is the biggest
#' number in the table and the one to act on.
#'
#' So `"subset"` is the default because it is right about the prior and never
#' loses beyond noise, not because it will visibly improve a fit.
#'
#' Two consequences of `"subset"` worth knowing. A rule on a factor is always
#' hard, even in a soft tree: a gate is a smooth function of the distance from a
#' cutpoint and there is no distance between two levels. Under `"onehot"` the
#' distance-based gate did apply to the 0 and 1 of an indicator, and at the
#' default bandwidth it left one level with a fractional membership weight for
#' 81% of cutpoints, which is not a smoothing anyone asked for. And a rule can no
#' longer land on an indicator that a path has already used up, which `"onehot"`
#' did for between 1% and 7% of its draws on a factor depending on `K`.
#'
#' A two-level factor is unaffected either way, since `2^2 - 2` is `B_2`.
#'
#' # Telling the prior what you already know
#'
#' `sparsity` and `split_prior` answer different questions and cannot both be in
#' force, so giving `split_prior` turns `sparsity` off. `sparsity` is for when you
#' do not know which predictors matter and want the prior to work it out from the
#' data; the splitting proportions are drawn, and a predictor can be dropped from
#' the forest entirely. `split_prior` is for when you do know something and want
#' it honored; the proportions are held at the weights you gave and nothing draws
#' over them.
#'
#' A weight is a statement about relative attention, not about effect size. It
#' changes how often the sampler proposes a split on a predictor, which is a
#' prior, so the data can still overrule it: a predictor given a large weight
#' whose splits do not improve the fit will collect rules that go nowhere, and
#' one given a small weight that genuinely matters will still be found, more
#' slowly. On pure noise, where nothing in the data prefers any predictor, the
#' realized share of splitting rules matches the weights closely, which is the
#' cleanest way to see what the argument does and the only case where the weights
#' fully determine the outcome.
#'
#' Because the weights are fixed, `split_prior` does not accumulate the
#' atom-at-zero mass described above. A predictor can still miss out on a rule in
#' some draw of a small forest, but that is sampling variation at a fixed
#' probability rather than a prior state that persists, and unlike the sparsity
#' prior it does go away as trees are added: on four noise predictors with
#' `split_prior = c(x1 = 8)`, `prop_used` runs from 0.89 to 1 at 20 trees and is
#' 1 for every predictor at 50 and at 200. The sparsity prior on the same data at
#' 50 trees left the most heavily weighted predictor out of 41% of draws. So
#' `split_prior` is a reasonable middle course when a particular contrast is the
#' estimand but the predictors are still too many to treat alike.
#'
#' One weight per term in the formula, not per column of the design matrix, so a
#' factor is named once and its levels share the weight, the way they already
#' share one entry of the sparsity prior.
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
#' # Progress
#'
#' `verbose = TRUE` prints a line to the console every `num_print` iterations,
#' which is the whole of what this package decides about progress. A progress
#' bar is \CRANpkg{progressr}'s business, and the sampler reports to it
#' unconditionally: nothing is shown unless a handler is active, so there is no
#' argument to switch on.
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
#' than four in sequence, and progress from chains running in parallel under
#' \CRANpkg{future} is relayed back as it arrives. Warmup and sampling are one
#' run for this purpose, because they cost the same per iteration and a bar that
#' restarts halfway is a worse report than one that does not.
#'
#' Reporting is capped at 50 steps per chain rather than one per iteration. The
#' callback itself is nothing next to a sweep, but a handler that redraws a bar
#' is not, and it would otherwise be possible for the reporting to cost more than
#' the sampling. Progress does not touch the draws: the same seed gives the same
#' fit whether or not anything was listening.
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
                             split_prior = NULL,
                             categorical = "subset",
                             k = 2,
                             bandwidth = 0.1,
                             chains = 1L,
                             num_burn = 500L, num_draws = 500L, num_thin = 1L,
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

  for (nm in c("num_burn", "num_draws", "num_thin", "num_print",
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

  x_transform <- arg::match_arg(x_transform, c("quantile", "range"))

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

  sparse <- resolve_sparsity(sparsity)
  augment <- resolve_augment(augment, soft)

  out <- list(num_trees = if (!is_null(num_trees)) stats::setNames(as.integer(num_trees), names(num_trees)),
              gate = gate,
              soft = soft,
              sparsity = sparsity,
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
# harder towards a few predictors.
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
