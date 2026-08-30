#' Response families for generalized BART
#'
#' @description
#' `bartisan()` accepts the [stats::family] objects used by [stats::glm()], so
#' `gaussian()`, `binomial("probit")`, `poisson()` and [stats::Gamma()] all work
#' unchanged. The functions documented here supply the additional families that
#' have no `glm()` counterpart, in the same style, so that they can be passed to
#' the `family` argument the same way.
#'
#' One thing to know about the gamma family, because base R's default is the
#' wrong choice here: `stats::Gamma()` defaults to `link = "inverse"`, and the
#' inverse link is the worst link for this sampler. **Write `Gamma("log")`.**
#' Base R's function is left as base R defines it, so that attaching this package
#' cannot change what `glm()` does; naming the family as a string,
#' `family = "Gamma"`, gets the log link, since that spelling is this package's
#' own. `bartisan()` says so when a link's
#' inverse does not cover the additive predictor, which is the case that catches
#' it; see Details.
#'
#' @param link the link function. Each family compiles the links for which the
#'   additive predictor is the natural unconstrained scale; any other link is
#'   applied from R, for the families where that is well defined. See Details.
#' @param theta for `negbin()` and `zi_negbin()`, a fixed value for the
#'   dispersion parameter. The default, `NULL`, draws it along with everything
#'   else.
#' @param phi for `Beta()` and `ordbeta()`, a fixed value for the beta
#'   precision. The default, `NULL`, draws it.
#' @param reference for `multinomial()`, the response category to hold as the
#'   reference. With the logit link the default, `NULL`, fits one forest per
#'   category instead and leaves the model unidentified, which is what makes the
#'   prior symmetric in the categories; see Details. The probit link is always
#'   written as contrasts against a reference, so there the default is the first
#'   level.
#' @param replicates for `multinomial("probit")`, how many simulation draws to
#'   use for the category probabilities, which have no closed form. Larger is
#'   more accurate and slower.
#' @param logdens for `custom_family()`, the log density. A function of the
#'   response and the additive predictors, `function(y, eta)`, where `y` is a
#'   numeric vector of length `n` and `eta` an `n` by `num_predictors` matrix,
#'   returning a numeric vector of length `n`. With nuisance parameters it takes
#'   a third argument, `function(y, eta, aux)`, where `aux` is a numeric vector
#'   of their current values. It is the log density of *one unit of prior
#'   weight*, so that `weights` behave as they do elsewhere, and terms free of
#'   `eta` may be dropped.
#' @param num_predictors for `custom_family()`, how many additive predictors the
#'   density has, that is, how many forests to fit.
#' @param start for `custom_family()`, the value each additive predictor starts
#'   at, in place of the intercept-only fit the compiled families use. One value
#'   or one per predictor.
#' @param derivatives for `custom_family()`, an optional
#'   `function(y, eta, h)` returning a list with elements `score` and `info`,
#'   the first derivative of `logdens` with respect to the `h`th predictor and
#'   minus its second derivative, each a vector of length `n`. The default,
#'   `NULL`, takes central differences of `logdens`. It covers the additive
#'   predictors only: a nuisance parameter is always differenced, which costs
#'   three calls per sweep rather than three per leaf.
#' @param aux_names for `custom_family()`, the names of the nuisance parameters
#'   to draw. Naming them is what declares them, because the names label the
#'   columns of `fit$aux` and are what `summary()` and `fit$rhat` report them
#'   under. The default, `NULL`, means none, unless `aux_start` is given, in
#'   which case they are named positionally.
#' @param aux_start for `custom_family()`, the value each nuisance parameter
#'   starts at. One value or one per parameter. The sampler will walk to the
#'   posterior from a poor start, so this need only be the right order of
#'   magnitude.
#' @param num_bins *Advanced.* For `ph()`, how many pieces the baseline hazard
#'   has, with the edges at evenly spaced quantiles of the observed times. The
#'   default, `NULL`, uses about the cube root of the sample size, which is the
#'   order the Freedman-Diaconis rule gives for a histogram. **You should not need
#'   to set this**: the estimates are flat in it over a sixty-fold range, and it
#'   is here for checking that rather than for tuning. See Details.
#' @param lambda_shape for `ph()`, the shape of the gamma prior on each bin's
#'   baseline hazard. Its rate is drawn.
#' @param update_lambda for `ph()`, whether to draw the baseline hazards. `FALSE`
#'   holds them at their prior mean, which is for diagnosis rather than analysis.
#' @param nu,q for `dpm()`, the degrees of freedom of the baseline's
#'   inverse-chi-square prior on a component's variance and the quantile of that
#'   prior placed at a rough estimate of the residual standard deviation. The
#'   defaults, 10 and 0.95, are the paper's, and are tighter than BART's own 3
#'   and 0.90 because the mixture covers small errors with extra components
#'   rather than with one component's left tail.
#' @param k_s for `dpm()`, how many units of the baseline's own scale the
#'   component means are allowed to reach out to. The default, 10, places the
#'   marginal of a component mean so that it reaches the largest residual of a
#'   linear fit.
#' @param alpha for `dpm()`, a fixed Dirichlet process concentration. The
#'   default, `NULL`, draws it.
#' @param max_clusters,psi for `dpm()`, the largest number of mixture components
#'   thought plausible and the shape of the taper towards it, which together set
#'   the prior on `alpha`. The defaults are a tenth of the sample size and 0.5.
#' @param name for `custom_family()`, a label used when printing the fit.
#'
#' @details
#' Every family reduces to a scalar additive predictor, or to several of them,
#' together with the first two derivatives of the log density with respect to
#' each. That is the whole interface the sampler needs, which is why the set of
#' available families is not restricted to the conditionally conjugate ones.
#'
#' `vignette("families", package = "bartisan")` covers all of this at length: how
#' to choose a family, how to choose among the links a family offers, and what
#' each family is and is not for. What follows is the short version, and the
#' points that could lead to output being misread.
#'
#' The supported families and links are:
#'
#' | Family | Links | Additive predictors | Drawn nuisance parameters |
#' |---|---|---|---|
#' | `gaussian()` | `identity` | 1 | residual standard deviation |
#' | `binomial()` | `logit`, `probit`, `cloglog` | 1 | none |
#' | `poisson()` | `log` | 1 | none |
#' | `negbin()` | `log` | 1 | dispersion |
#' | `Gamma("log")` | `log` (also `inverse`, `identity`, any link) | 1 | shape |
#' | `ordinal()` | `logit`, `probit`, `cloglog` | 1 | cutpoints |
#' | `multinomial()` | `logit`, `probit` | one per category, or per non-reference level | latent covariance, for the probit link |
#' | `weibull_aft()`, `loglogistic_aft()`, `lognormal_aft()` | none | 1 | scale |
#' | `ph()` | none | 1 | baseline hazard per bin |
#' | `dpm_aft()` | none | 1 | error mixture, concentration |
#' | `location_scale()` | `identity` | 2 | none |
#' | `zi_poisson()` | `log` | 2 | none |
#' | `zi_negbin()` | `log` | 2 | dispersion |
#' | `Beta()` | `logit`, `probit`, `cloglog` | 1 | precision |
#' | `ordbeta()` | `logit` | 1 | 2 cutpoints, precision |
#' | `dpm()` | `identity` | 1 | error mixture, concentration |
#'
#' A family with more than one additive predictor fits one forest per predictor.
#' Nuisance parameters are drawn alongside the trees and reported in `fit$aux`.
#'
#' # Links the engine does not compile
#'
#' The links listed above are the ones the sampler evaluates in compiled code.
#' Any other link is accepted for `gaussian()`, `binomial()`, `poisson()`,
#' `negbin()` and `Gamma()`, and applied from R by composing the caller's inverse
#' link with the family's own, with the chain rule carrying the derivatives back.
#' So `binomial("cauchit")` works, as does any link object of the kind
#' [stats::make.link()] returns. It costs a call into R for every leaf the
#' sampler visits, and the leaf prior scale is calibrated for the compiled link.
#'
#' A link whose inverse has a restricted range -- `Gamma("inverse")`,
#' `Gamma("identity")`, `poisson("identity")` -- will give non-finite densities
#' for some predictors. Those proposals are rejected rather than breaking the
#' chain, but they are wasted work and the fit is worse for it, so `bartisan()`
#' says so when it starts. Prefer links whose inverse is defined on the whole
#' line.
#'
#' The families with more than one additive predictor, or whose link enters
#' somewhere other than a single mean -- `ordinal()`, `multinomial()`, the
#' accelerated failure time families, `location_scale()`, the zero-inflated
#' families and `ordbeta()` -- take only their listed links. `custom_family()` is
#' the way to reach anything else.
#'
#' # If no family is given
#'
#' `family` may be omitted, in which case it is read off the response: a `Surv`
#' object gets `dpm_aft()`, an ordered factor `ordinal()`, a two-valued
#' response `binomial()`, any other factor `multinomial()`, and any numeric
#' response `dpm()`. `dpm()` cannot take prior weights, so a weighted fit with no
#' family named is an error rather than a silent substitution. The choice is reported with a message, which setting `family`
#' silences. A count is *not* given `poisson()` and a numeric response taking two
#' values other than 0 and 1 is *not* given `binomial()`, since either would be a
#' modeling decision rather than a reading of the response's type.
#'
#' # What to know before reading the output
#'
#' `ordinal()` accepts a numeric response as well as an ordered factor, taking
#' its sorted unique values as the categories, and that is **a method rather than
#' a fallback**: the cutpoints absorb the marginal distribution of the response
#' and the forest explains only the ordering, so nothing is assumed about the
#' error distribution and the model for \eqn{P(Y \le y \mid x)} is invariant to
#' any monotone transformation of the response. Predict with
#' `type = "mean"`. Bin the response onto twenty-odd quantiles first: one
#' cutpoint per distinct value costs 73 seconds against 2.7 for `gaussian()` at
#' 1000 observations, and twenty-five bins was both sixteen times faster and
#' slightly more accurate. See the vignette.
#'
#' `ordinal()` uses the cumulative-link parameterization of \pkgfun{MASS}{polr},
#' in which \eqn{P(Y \le k) = F(c_k - \eta)}, so larger values of the additive
#' predictor shift mass towards higher categories. Only the differences
#' \eqn{c_k - \eta_i} are identified, so one location has to be pinned: with three
#' or more categories the draws are reported in the chart where **the additive
#' predictor has mean zero over the fitted sample and every cutpoint is free**,
#' which is the chart `polr()` reports in when its predictors are centered. With
#' exactly two categories the single boundary is folded into the intercept
#' instead, so a two-category response is exactly binary regression with the
#' matching link and on the same scale. `cut1` is therefore a free parameter
#' rather than a constant zero, which is a change from earlier versions.
#'
#' `multinomial()` by default fits one forest per category and leaves the model
#' unidentified, since adding any function of the predictors to every category's
#' forest leaves the probabilities alone. This is the parameterization of Murray
#' (2021), whose point is that the prior is then symmetric in the categories;
#' every identified quantity is recovered from the draws. Passing `reference`
#' instead pins that category at zero and fits one fewer forest, giving log odds
#' against it.
#'
#' `multinomial("probit")` lets the latent utilities correlate, which a
#' multinomial logit cannot express at all. \eqn{\Sigma} is normalized by the
#' trace constraint \eqn{\mathrm{tr}(\Sigma) = C} (Burgette and Nordheim 2012),
#' and its lower triangle appears in `fit$aux` as `sigma11`, `sigma21` and so on.
#' **Those correlations are weakly identified: read the fitted probabilities, not
#' the covariance.** They enter the likelihood only through orthant probabilities
#' of a distribution whose location is a sum of trees, so a flexible mean absorbs
#' much of the dependence they are meant to measure. At 900 observations a true
#' correlation of zero came back as -0.57, and posterior intervals ran up to 1.07
#' wide on a parameter confined to \eqn{(-1, 1)}; at 3000 observations the
#' posterior tracks the truth to within about 0.2. Two further consequences: the
#' likelihood has **no closed form**, so it and every category probability are
#' simulated with `replicates` draws, and `augment` does not apply, because the
#' latent variables are the model rather than a rewriting of it. The sampler is
#' Algorithm P2 of Xu et al. (2025).
#'
#' `dpm()` is not a distribution but a way of not choosing one. It is DPMBART
#' (George et al. 2019): a numeric response with the sum of trees for its mean, as
#' `gaussian()` has, and a Dirichlet process mixture of normals for its errors
#' instead of a single normal, so the error distribution comes out as whatever
#' mixture the data ask for. [error_density()] gives that density, which is the
#' object the method exists to produce.
#'
#' **It is the family to reach for by default on a numeric response**, because it
#' does not pay for its flexibility: on normal errors, where `gaussian()` is
#' exactly right, it came out slightly ahead on both held-out error and log score
#' at the same time to one decimal place, and on heavy-tailed, skewed and bimodal
#' errors it was ahead by a great deal -- on bimodal errors at a thousand
#' observations, 0.050 against 0.154 in held-out RMSE, a factor of three, at the
#' same time to a tenth of a second. So it is the family a numeric response gets
#' when none is named. The reasons to prefer `gaussian()` are not statistical:
#' it takes **prior weights**, which `dpm()` refuses, and it reports one
#' interpretable `sigma` where `dpm()` has a mixture. It is also faster, by 1.4
#' times at a thousand observations. The vignette has the comparison.
#'
#' Two things to know about `dpm()` itself. **It does not buy
#' heteroskedasticity** -- the error distribution is flexible but it is the same
#' distribution at every \eqn{x}, and `location_scale()` is the family for a
#' spread that depends on the predictors. And **the additive predictor is the
#' conditional mean**, as it is for `gaussian()`: nothing in the model forces the
#' mixture to be centered, so the sampler works in a chart where only the sum of
#' the predictor and the error mean is identified, but reporting is done in the
#' chart where the mixture has mean zero and the whole conditional mean sits on
#' the predictor. `type = "link"` and `type = "response"` therefore agree
#' exactly, and `fit$aux` reports the shift that was taken out as `center` rather
#' than an error mean, which is zero by construction. Prior weights are refused,
#' since a weight would have to be a multiplicity in the Dirichlet process.
#'
#' The gamma family puts the forest on the log mean and draws the shape, which
#' acts as the inverse dispersion; it does *not* regress the shape on the
#' predictors. `negbin()` and `ordbeta()` take `theta` and `phi` to fix their
#' equivalents, and the gamma shape has no such argument because a caller who
#' knows it is rare. **The link is where the care is needed.** Only `log` is
#' compiled, and the base R default of `inverse` is the worst case for this
#' sampler: its inverse maps a negative predictor to a negative mean, whose log
#' is not a number, so the proposal is rejected -- dozens of times per fit.
#' Measured on 600 observations and 50 trees, `stats::Gamma()` took 7.2 seconds
#' against 3.8 for `Gamma("log")` and fitted the mean slightly worse. So write
#' the link, or name the family as the string `"Gamma"`, which resolves to the
#' log link. Any composed link whose inverse has a restricted range is reported
#' when the fit starts.
#'
#' The accelerated failure time families expect a right-censored response,
#' supplied either as a \pkgfun{survival}{Surv} object or as a two-column matrix
#' of times and event indicators. They model \eqn{\log T = \eta + \sigma\epsilon}
#' with \eqn{\epsilon} standard Gumbel, logistic or normal respectively, giving
#' Weibull, log-logistic and log-normal survival times, so the predictor is a log
#' time ratio in each.
#'
#' `ph()` is the proportional hazards alternative, with a piecewise-constant
#' baseline: \eqn{\lambda(t \mid x) = \lambda_0(t)\exp(r(x))}, so its predictor
#' is a log *hazard* ratio and the baseline is free to take any shape rather than
#' the monotone one a Weibull imposes. `num_bins` sets how many pieces, with the
#' edges at evenly spaced quantiles of the observed times; the default is about
#' the cube root of the sample size. The bin hazards are drawn from their exact
#' gamma conditionals and reported as `lambda1`, `lambda2`, ... in `fit$aux`,
#' together with the rate of their own prior. The predictor and the baseline are
#' identified only jointly, so the baseline carries the level and the predictor
#' is reported centered on it.
#'
#' Cox's *partial* likelihood is what cannot be used here: it couples
#' observations through risk sets and so does not decompose into a sum over the
#' observations reaching a leaf. The full likelihood of the piecewise-exponential
#' model does decompose, and it approaches the partial likelihood as the bins
#' shrink, which is how `ph()` reaches proportional hazards without it.
#'
#' **`num_bins` is not a modeling decision, and its default should be left
#' alone.** It is exposed for checking that, not for tuning. Measured over three
#' replicates at 700 observations, sweeping it from 4 to 250 -- a sixty-fold
#' range, against a baseline hazard that turns over and against a Weibull one --
#' moved the error in the survival function between 0.035 and 0.048 and the error
#' in the log hazard ratio between 0.135 and 0.185, with no trend in either and
#' every difference inside the replicate-to-replicate spread. What the bin count
#' does change is the effective number of parameters, which grows with it: from 17
#' at four bins to 206 at 250. That is what makes the default matter for `loo()`
#' and `waic()` rather than for the estimates -- one parameter per event time
#' would leave each observation's density inflated by a parameter only it informs,
#' and leave-one-out unable to do its job.
#'
#' The three differ in cost, though not enough to decide a model on.
#' `lognormal_aft()` and `loglogistic_aft()` impute each censored failure time
#' above its censoring time, which makes their targets quadratic and is worth 8
#' to 30 times the speed; `weibull_aft()` needs no imputation because its
#' likelihood already has a form the sampler can collapse to a single pass, but
#' only under hard rules, which makes it the slowest of the three at the default
#' gate.
#'
#' `dpm_aft()` is the accelerated failure time model with the error distribution
#' estimated rather than assumed: \eqn{\log T = m(x) + W} with \eqn{W} a
#' Dirichlet process mixture of normals constrained to mean zero, and censored
#' log-times imputed. It is `dpm()`'s error model with censoring, so its predictor
#' is the conditional mean of \eqn{\log T}, `error_density()` reports the fitted
#' error density, and prior weights are refused for the same reason `dpm()`
#' refuses them. Following Henderson, Louis, Rosner and Varadhan (2020).
#'
#' Reach for it when the shape of the error is in doubt and you would rather not
#' assert one. Measured against a two-component error it was worth 210 held-out
#' log points and a third of the error in \eqn{S(t \mid x)} over the best
#' fixed-error family; against a log-normal error, where `lognormal_aft()` is
#' correctly specified, the two were within 0.1 log points of each other. So it
#' gains where the assumption would have been wrong and costs nothing where it
#' would have been right, which is the property `dpm()` has against `gaussian()`.
#'
#' `location_scale()` regresses the mean and the log standard deviation of a
#' normal response on separate forests, so the variance is an unrestricted
#' function of the predictors.
#'
#' `zi_poisson()` and `zi_negbin()` are zero-inflated counts. Both parts get their
#' own forest: the first predictor is the log mean of the count component and the
#' second the log odds that an observation is a structural zero, so the
#' excess-zero mechanism is free to depend on the predictors. The two are
#' reported as the `count` and `zero` predictors.
#'
#' `Beta()` is beta regression for a response strictly inside the unit interval:
#' a forest on the link of the mean, and a precision drawn alongside it. A
#' response *at* either endpoint has no beta density, so it is an error rather
#' than something to nudge inward.
#'
#' `ordbeta()` is the ordered beta regression of Kubinec (2023), for a response on
#' the closed unit interval with point masses at zero and one. One predictor
#' drives both the probability of landing on an endpoint, through a pair of
#' cutpoints as in an ordinal model, and the mean of the beta density in between.
#' Because the predictor also enters the beta mean it is identified, so unlike
#' `ordinal()` both cutpoints are drawn.
#'
#' Choose between them on whether the response can reach a boundary, not on
#' whether it happens to in the sample: the two ask different questions, and
#' `ordbeta()` fitted to a response with no boundary observations leaves its
#' cutpoints with nothing to identify them.
#'
#' # Supplying a likelihood
#'
#' `custom_family()` takes the log density itself, as an R function, and fits the
#' model that goes with it. Nothing else about the sampler changes: the
#' leaf-level Laplace proposal needs the first two derivatives of the log density
#' with respect to each additive predictor and nothing more, and central
#' differences of the supplied function produce both.
#'
#' ```r
#' # A Poisson model written out by hand. Terms free of eta may be dropped;
#' # they cancel from every acceptance ratio.
#' pois <- custom_family(function(y, eta) y * eta[, 1] - exp(eta[, 1]),
#'                       start = log(mean(d$y)))
#'
#' # Two predictors: a mean and a log standard deviation.
#' ls <- custom_family(function(y, eta) dnorm(y, eta[, 1], exp(eta[, 2]),
#'                                            log = TRUE),
#'                     num_predictors = 2, start = c(0, 0))
#' ```
#'
#' The function is called once per leaf per Fisher-scoring step with the
#' observations reaching that leaf, so it must be vectorized over `y` and the
#' rows of `eta`; it must not be vectorized *within* an observation, and it must
#' return exactly one value per row. Supplying `derivatives` cuts three calls to
#' one and removes the differencing error.
#'
#' Nuisance parameters are drawn alongside the trees when `aux_names` names
#' them, and `logdens` then takes a third argument holding their current values:
#'
#' ```r
#' # A Gaussian written out by hand, with its scale drawn rather than fixed.
#' by_hand <- custom_family(
#'   logdens = function(y, eta, aux) dnorm(y, eta[, 1], exp(aux[1]), log = TRUE),
#'   aux_names = "log_sigma", aux_start = 0)
#' ```
#'
#' They are reported in `fit$aux` under those names, and covered by `summary()`
#' and `fit$rhat` like any other family's. There is no prior argument and no
#' bounds argument, because a nuisance parameter here is carried as an additive
#' predictor whose forest is pinned at depth zero -- one tree that can never
#' split, so the forest is a single scalar -- and it is drawn by the same
#' Laplace-plus-Metropolis step as any leaf, under that step's Gaussian leaf
#' prior. So a parameter with a restricted range is handled the way it would be
#' for a real predictor, by writing the transform into `logdens`: the `exp()`
#' above is what keeps the scale positive.
#'
#' What `custom_family()` does not do: the response must be numeric, so a factor
#' has to be coded first; and since the package cannot know what the mean of the
#' density is, `predict(type = "response")` returns the additive predictors
#' rather than a fitted mean.
#'
#' @returns
#' A list of class `bartisan_family`, containing at least the elements `family`
#' and `link`. Objects of this class are recognized by `bartisan()` alongside
#' ordinary [stats::family] objects.
#'
#' @references
#' Burgette, L. F., & Nordheim, E. V. (2012). The trace restriction: an
#' alternative identification strategy for the Bayesian multinomial probit
#' model. *Journal of Business & Economic Statistics*, 30(3), 404--410.
#' \doi{10.1080/07350015.2012.680416}
#'
#' George, E., Laud, P., Logan, B., McCulloch, R., & Sparapani, R. (2019). Fully
#' nonparametric Bayesian additive regression trees. In *Topics in
#' Identification, Limited Dependent Variables, Partial Observability,
#' Experimentation, and Flexible Modeling: Part B* (Advances in Econometrics,
#' vol. 40B, pp. 89--110). Emerald Publishing.
#' \doi{10.1108/S0731-90532019000040B006}
#'
#' Kubinec, R. (2023). Ordered beta regression: a parsimonious, well-fitting
#' model for continuous data with lower and upper bounds. *Political Analysis*,
#' 31(4), 519--536. \doi{10.1017/pan.2022.20}
#'
#' Murray, J. S. (2021). Log-linear Bayesian additive regression trees for
#' multinomial logistic and count regression models. *Journal of the American
#' Statistical Association*, 116(534), 756--769.
#' \doi{10.1080/01621459.2020.1813587}
#'
#' Xu, Y., Hogan, J., Daniels, M., Kantor, R., & Mwangi, A. (2025). Augmentation
#' samplers for multinomial probit Bayesian additive regression trees. *Journal
#' of Computational and Graphical Statistics*, 34(2), 498--508.
#' \doi{10.1080/10618600.2024.2388605}
#'
#' @seealso
#' [bartisan()], [error_density()], and
#' `vignette("families", package = "bartisan")` for the long form.
#'
#' @examplesIf FALSE
#' bartisan(y ~ ., data = d, family = negbin())
#' bartisan(y ~ ., data = d, family = ordinal("probit"))
#' bartisan(survival::Surv(time, status) ~ ., data = d, family = weibull_aft())
#'
#' # The same response under proportional hazards, with a free baseline.
#' bartisan(survival::Surv(time, status) ~ ., data = d, family = ph())

#'
#' # An accelerated failure time model with the error distribution estimated.
#' bartisan(survival::Surv(time, status) ~ ., data = d, family = dpm_aft())
#' bartisan(count ~ ., data = d, family = zi_negbin())
#' bartisan(proportion ~ ., data = d, family = Beta())
#'
#' # The same response, when it can also sit exactly at 0 or 1.
#' bartisan(proportion ~ ., data = d, family = ordbeta())
#'
#' # The latent covariance is reported in `aux`, as its lower triangle.
#' fit <- bartisan(y ~ ., data = d, family = multinomial("probit"))
#' colMeans(fit$aux)
#'
#' # A numeric response whose error distribution is estimated rather than
#' # assumed. `error_density()` reports what shape it came out.
#' fit <- bartisan(y ~ ., data = d, family = dpm())
#' error_density(fit)
#'
#' # A link the engine does not compile, applied from R.
#' bartisan(y ~ ., data = d, family = binomial("cauchit"))
#'
#' # A likelihood supplied from R.
#' bartisan(y ~ ., data = d,
#'         family = custom_family(function(y, eta) y * eta[, 1] - exp(eta[, 1])))
#'
#' @name bartisan-families
NULL

#' @rdname bartisan-families
#' @export
negbin <- function(link = "log", theta = NULL) {
  link <- arg::match_arg(link, "log")

  arg::when_not_null(
    theta,
    arg::arg_and(
      arg::arg_number,
      arg::arg_gt(0)
    )
  )

  new_bartisan_family("negbin", link, theta = theta)
}

# `family = "Gamma"` resolves through here, which is why this exists and why it
# is not exported: `stats::Gamma()` defaults to the inverse link, which is the
# worst link for this sampler, and masking a base R family to change that default
# would change `glm()` too. Naming the family as a string is the one route this
# package owns outright, so that is the route that gets the better default.
# `stats::Gamma()` and `Gamma("inverse")` behave exactly as they always did, and
# the message from `warn_restricted_link()` says what the inverse link costs.
gamma_string_family <- function() {
  stats::Gamma(link = "log")
}

#' @rdname bartisan-families
#' @export
ordinal <- function(link = "logit") {
  link <- arg::match_arg(link, c("logit", "probit", "cloglog"))

  new_bartisan_family("ordinal", link)
}

#' @rdname bartisan-families
#' @export
multinomial <- function(link = "logit", reference = NULL,
                        replicates = 200L) {
  link <- arg::match_arg(link, c("logit", "probit"))

  if (!is_null(reference)) {
    if (length(reference) != 1L || is.na(reference)) {
      arg::err("{.arg reference} must be a single response category")
    }

    reference <- as.character(reference)
  }

  # The two links are different enough inside the engine to be separate
  # families -- the probit one carries a covariance matrix and has one fewer
  # forest under the same coding -- but they are one model to the caller, so
  # they are one function with a link argument, as `binomial()` and `ordinal()`
  # are.
  if (identical(link, "logit")) {
    return(new_bartisan_family("multinomial", link, reference = reference))
  }

  arg::arg_count(replicates)
  arg::arg_gte(replicates, 1)

  new_bartisan_family("mnp", link, reference = reference,
                      replicates = as.integer(replicates))
}

#' @rdname bartisan-families
#' @export
dpm_aft <- function(nu = 10, q = 0.95, k_s = 10, alpha = NULL,
                    max_clusters = NULL, psi = 0.5) {
  out <- dpm(nu = nu, q = q, k_s = k_s, alpha = alpha,
             max_clusters = max_clusters, psi = psi)
  out[["family"]] <- "dpm_aft"
  out
}

#' @rdname bartisan-families
#' @export
dpm <- function(nu = 10, q = 0.95, k_s = 10, alpha = NULL,
                max_clusters = NULL, psi = 0.5) {
  arg::arg_number(nu)
  arg::arg_gt(nu, 0)
  arg::arg_number(q)
  arg::arg_between(q, c(0, 1))
  arg::arg_number(k_s)
  arg::arg_gt(k_s, 0)
  arg::arg_number(psi)
  arg::arg_gte(psi, 0)

  arg::when_not_null(
    alpha,
    arg::arg_and(
      arg::arg_number,
      arg::arg_gt(0)
    )
  )

  arg::when_not_null(
    max_clusters,
    arg::arg_and(
      arg::arg_count,
      arg::arg_gte(2)
    )
  )

  new_bartisan_family("dpm", "identity", nu = nu, q = q, k_s = k_s,
                      alpha = alpha, max_clusters = max_clusters, psi = psi)
}

#' @rdname bartisan-families
#' @export
weibull_aft <- function() {
  new_bartisan_family("aft", "weibull")
}

#' @rdname bartisan-families
#' @export
loglogistic_aft <- function() {
  new_bartisan_family("aft", "loglogistic")
}

#' @rdname bartisan-families
#' @export
lognormal_aft <- function() {
  new_bartisan_family("aft", "lognormal")
}

#' @rdname bartisan-families
#' @export
ph <- function(num_bins = NULL, lambda_shape = 1, update_lambda = TRUE) {
  arg::when_not_null(
    num_bins,
    arg::arg_and(
      arg::arg_whole_number,
      arg::arg_gte(2)
    )
  )

  arg::arg_number(lambda_shape)
  arg::arg_gt(lambda_shape, 0)
  arg::arg_flag(update_lambda)

  new_bartisan_family("ph", "log",
                      num_bins = if (is_null(num_bins)) NULL
                                else as.integer(num_bins),
                     lambda_shape = lambda_shape,
                     update_lambda = update_lambda)
}

#' @rdname bartisan-families
#' @export
location_scale <- function(link = "identity") {
  link <- arg::match_arg(link, "identity")

  new_bartisan_family("location_scale", link)
}

#' @rdname bartisan-families
#' @export
zi_poisson <- function(link = "log") {
  link <- arg::match_arg(link, "log")

  new_bartisan_family("zip", link)
}

#' @rdname bartisan-families
#' @export
zi_negbin <- function(link = "log", theta = NULL) {
  link <- arg::match_arg(link, "log")

  arg::when_not_null(
    theta,
    arg::arg_and(
      arg::arg_number,
      arg::arg_gt(0)
    )
  )

  new_bartisan_family("zinb", link, theta = theta)
}

#' @rdname bartisan-families
#' @export
Beta <- function(link = "logit", phi = NULL) {
  # Only the logit link is compiled; `probit` and `cloglog` are reached by
  # composition, the way an uncompiled link is for any other family, so they are
  # not listed here and are validated by `as_bartisan_family()` instead.
  arg::arg_string(link)

  arg::when_not_null(
    phi,
    arg::arg_and(
      arg::arg_number,
      arg::arg_gt(0)
    )
  )

  new_bartisan_family("beta", link, phi = phi)
}

#' @rdname bartisan-families
#' @export
ordbeta <- function(link = "logit", phi = NULL) {
  link <- arg::match_arg(link, "logit")

  arg::when_not_null(
    phi,
    arg::arg_and(
      arg::arg_number,
      arg::arg_gt(0)
    )
  )

  new_bartisan_family("ordbeta", link, phi = phi)
}

#' @rdname bartisan-families
#' @export
custom_family <- function(logdens, num_predictors = 1L, start = 0,
                          derivatives = NULL, aux_names = NULL, aux_start = 0,
                          name = "custom") {
  if (!is.function(logdens)) {
    arg::err("{.arg logdens} must be a function of the response and the
              additive predictors")
  }

  arg::arg_whole_number(num_predictors)
  arg::arg_gte(num_predictors, 1)
  arg::arg_numeric(start)
  arg::arg_numeric(aux_start)
  arg::arg_string(name)

  arg::when_not_null(
    derivatives,
    arg::arg_function
  )

  arg::when_not_null(
    aux_names,
    arg::arg_character
  )

  num_predictors <- as.integer(num_predictors)

  if (length(start) != 1L && length(start) != num_predictors) {
    arg::err("{.arg start} must have one value, or one per additive predictor
              ({num_predictors})")
  }

  # The nuisance parameters are declared by naming them, because their names are
  # what labels the columns of `fit$aux` and what `summary()` and `fit$rhat`
  # report them under. Giving only starting values names them positionally.
  if (is_null(aux_names)) {
    aux_names <- {
      if (missing(aux_start)) character()
      else paste0("aux", seq_along(aux_start))
    }
  }

  num_aux <- length(aux_names)

  if (anyDuplicated(aux_names) > 0L || any(!nzchar(aux_names))) {
    arg::err("{.arg aux_names} must be distinct and non-empty")
  }

  if (num_aux > 0L && length(aux_start) != 1L &&
      length(aux_start) != num_aux) {
    arg::err("{.arg aux_start} must have one value, or one per nuisance
              parameter ({num_aux})")
  }

  if (num_aux > 0L && length(formals(logdens)) < 3L) {
    arg::err("{.arg logdens} must take a third argument for the nuisance
              parameters when there are any.",
             i = "It is called as {.code logdens(y, eta, aux)}, with {.arg aux} a
                  numeric vector of length {num_aux}.")
  }

  new_bartisan_family("custom", "identity", logdens = logdens,
                      num_predictors = num_predictors,
                     start = rep(start, length.out = num_predictors),
                     derivatives = derivatives,
                     num_aux = num_aux,
                     aux_names = aux_names,
                     aux_start = rep(aux_start, length.out = max(num_aux, 1L)),
                     name = name)
}

bartisan_family_names <- c("gaussian", "binomial", "poisson", "negbin", "Gamma",
                           "ordinal", "multinomial", "dpm",
                          "weibull_aft", "loglogistic_aft", "lognormal_aft",
                          "location_scale", "zi_poisson", "zi_negbin",
                          "Beta", "ordbeta", "ph", "dpm_aft")

new_bartisan_family <- function(family, link, ...) {
  structure(c(list(family = family, link = link), list(...)),
            class = c("bartisan_family", "family"))
}

# The scale each engine family's additive predictor natively lives on. A link
# the engine does not carry is handled by composing the caller's inverse link
# with the link named here, which is why only families with a single mean and a
# conventional link can take an arbitrary one.
native_links <- c(gaussian = "identity", binomial = "logit", poisson = "log",
                  negbin = "log", beta = "logit")

valid_links <- list(custom = "identity",
                    gaussian = "identity",
                    dpm = "identity",
                    dpm_aft = "identity",
                    binomial = c("logit", "probit", "cloglog"),
                    poisson = "log",
                    negbin = "log",
                    Gamma = "log",
                    ordinal = c("logit", "probit", "cloglog"),
                    multinomial = "logit",
                    mnp = "probit",
                    aft = c("weibull", "loglogistic", "lognormal"),
                    ph = "log",
                    location_scale = "identity",
                    zip = "log",
                    zinb = "log",
                    beta = "logit",
                    ordbeta = "logit")

# Normalize whatever the user passed to `family` into a bartisan family object.
# Accepts a string, a family-generating function, a stats::family object, or one
# of the bartisan families above, mirroring how glm() resolves the argument.
# The family to use when the caller does not name one, read off the shape of the
# response the way `glm()` reads off nothing and `lm()` assumes everything. The
# rules are the ones that have a single obvious answer; anything else is Gaussian,
# which is the assumption a caller who did not think about it is making anyway.
#
# Deliberately *not* inferred: a count. A non-negative integer response is often
# Poisson and often not, and `poisson()` carries a variance assumption strong
# enough that guessing it would be a substantive modeling decision made
# silently. Gaussian is the weaker guess and the one that is easy to see is
# wrong.
default_family <- function(y, weights = NULL) {
  chosen <- {
    if (inherits(y, "Surv") || is_survival_matrix(y)) {
      "dpm_aft"
    }
    else if (is.ordered(y)) {
      "ordinal"
    }
    else if (is_binary(y)) {
      "binomial"
    }
    else if (is.factor(y) || is.character(y)) {
      "multinomial"
    }
    else if (is.matrix(y) && ncol(y) == 2L) {
      # Successes and failures, which is one of the shapes glm() accepts.
      "binomial"
    }
    else {
      "dpm"
    }
  }

  # Neither `dpm()` nor `dpm_aft()` can take prior weights, and quietly dropping
  # them or quietly swapping the family would each be worse than saying so: the
  # caller asked for weights and has not asked for a family, and only they can
  # say which of the two they meant.
  if (!is_null(weights)) {
    if (identical(chosen, "dpm")) {
      arg::err("a numeric response defaults to {.fn dpm}, which does not take
                prior weights",
               i = "name a family: {.code family = gaussian()} keeps the weights,
                    and so do {.fn ordinal} and {.fn location_scale}")
    }

    if (identical(chosen, "dpm_aft")) {
      arg::err("a censored response defaults to {.fn dpm_aft}, which does not
                take prior weights",
               i = "name a family: {.code family = lognormal_aft()} keeps the
                    weights, and so do {.fn weibull_aft},
                    {.fn loglogistic_aft} and {.fn ph}")
    }
  }

  arg::msg(c(i = "Using {.code family = {chosen}()}.",
             i = "Set {.arg family} to choose another, which also silences
                  this message."))

  get(chosen, mode = "function", envir = asNamespace("bartisan"))()
}



# Two levels, or numeric zeros and ones -- the responses for which a binomial
# model is the only sensible reading. A numeric response with two values that
# are not zero and one is left to the Gaussian default, since treating `c(1, 2)`
# as a success indicator would be a guess about which value is the success.
is_binary <- function(y) {
  if (is.logical(y)) {
    return(TRUE)
  }

  if (is.factor(y) || is.character(y)) {
    return(length(unique(y[!is.na(y)])) == 2L)
  }

  # Exactly two values, not "at most two": a constant response is degenerate
  # whatever the family, and reading `rep(1, n)` as an all-successes binomial
  # would replace a clear complaint about no variation with a silent fit.
  if (is.numeric(y) && !is.matrix(y)) {
    values <- unique(y[!is.na(y)])
    return(length(values) == 2L && all(values %in% c(0, 1)))
  }

  FALSE
}

# A two-column matrix of times and event indicators, which is the form the
# survival families accept alongside a Surv object.
is_survival_matrix <- function(y) {
  if (!is.matrix(y) || ncol(y) != 2L || !is.numeric(y)) {
    return(FALSE)
  }

  events <- y[, 2L]
  events <- events[!is.na(events)]

  # Non-negative times and a 0/1 second column. A two-column binomial response
  # is counts of successes and failures, so it fails the second test unless
  # every count is zero or one -- and then it is ambiguous, and survival is the
  # reading a two-column *numeric* matrix more often has.
  all(y[, 1L] >= 0, na.rm = TRUE) && length(events) > 0L &&
    all(events %in% c(0, 1))
}

as_bartisan_family <- function(family) {
  if (is.character(family)) {
    arg::arg_string(family)
    arg::arg_element(family, bartisan_family_names)

    # "Gamma" is the one name whose function is not the exported family: see
    # gamma_string_family().
    family <- {
      if (identical(family, "Gamma")) gamma_string_family
      else get(family, mode = "function", envir = asNamespace("bartisan"))
    }
  }

  if (is.function(family)) {
    family <- family()
  }

  if (!inherits(family, "family")) {
    arg::err("{.arg family} must be a family name, a family function, or a
              family object, such as {.code binomial(\"logit\")}")
  }

  name <- family[["family"]]
  link <- family[["link"]]

  if (!name %in% names(valid_links)) {
    arg::err("family {.val {name}} is not supported by {.fn bartisan}. Supported
              families are {.val {names(valid_links)}}")
  }

  # The gamma family is the one place a supplied link is overruled rather than
  # honored. The inverse link -- base R's default, because it is the canonical
  # link for the gamma -- needs a positive mean, and this sampler's additive
  # predictor is unconstrained: a draw that wanders non-positive has no gamma
  # density, so it is rejected, and any prediction there is NaN. Measured over
  # eight replicates on heavy-tailed data, five of them had such draws. The log
  # link is defined on the whole line, was faster and more accurate in every
  # setting measured, and is what a gamma regression almost always means in
  # practice. So it is used regardless, and the caller is told once.
  if (identical(name, "Gamma") && !identical(link, "log")) {
    arg::msg(c(i = "The {.val {link}} link is ignored: {.fn Gamma} is fitted on
                    the log link here.",
               i = "The additive predictor is unconstrained, and only the log
                    link's inverse keeps the mean positive over the whole line;
                    the others give rejected draws and {.val NaN} predictions.",
               i = "Write {.code Gamma(\"log\")} to silence this."))
    link <- "log"
  }

  allowed <- valid_links[[name]]
  custom_link <- NULL

  if (!link %in% allowed) {
    if (!name %in% names(native_links)) {
      arg::err("the {.val {link}} link is not supported for the {.val {name}}
                family. Supported links are {.val {allowed}}")
    }

    # A link the engine does not carry is honored by composing it onto the
    # scale the engine's family works on. That needs the inverse link, which a
    # stats::family object carries; a bare name is resolved through make.link().
    custom_link <- link_functions(family, link)
  }

  # A bartisan family carries its own options -- a fixed dispersion, a reference
  # category, a supplied log density -- which have to survive; a stats::family
  # object carries link machinery that has already been read off above.
  extra <- {
    if (inherits(family, "bartisan_family")) {
      family[!names(family) %in% c("family", "link")]
    }
    else list()
  }

  family <- do.call(new_bartisan_family,
                    c(list(family = name, link = link), extra))

  if (!is_null(custom_link)) {
    family[["custom_link"]] <- custom_link
    family[["native_link"]] <- native_links[[name]]
    warn_restricted_link(custom_link, native_links[[name]], link)
  }

  family
}

# A composed link is only as good as its range. The engine's predictor is
# unconstrained, so if the caller's inverse link maps some of the line outside
# the domain the engine's own link needs -- a negative mean for `log`, something
# outside the unit interval for `logit` -- those predictors give a non-finite
# density. The sampler rejects them, which is correct but is wasted work, and it
# is almost never what the caller intended: `stats::Gamma()`'s default inverse
# link is the case that turns up in practice, and it cost 1.9 times the time of
# the log link and fitted slightly worse on 600 observations. So say so once,
# rather than leaving it to be discovered.
warn_restricted_link <- function(custom_link, native, link) {
  if (identical(native, "identity")) {
    return(invisible())
  }

  # A grid wide enough to cover any predictor the leaf prior would produce.
  mu <- suppressWarnings(custom_link[["linkinv"]](seq(-6, 6, length.out = 121L)))

  inside <- {
    if (identical(native, "log")) isTRUE(all(is.finite(mu) & mu > 0))
    else isTRUE(all(is.finite(mu) & mu > 0 & mu < 1))
  }

  if (inside) {
    return(invisible())
  }

  arg::msg(c(i = "The {.val {link}} link's inverse does not cover the whole
                  additive predictor, so some proposals will have a non-finite
                  density and be rejected.",
             i = "The fit is valid but slower and less accurate; a link whose
                  inverse is defined on the whole line, such as
                  {.val {native}}, avoids it."))

  invisible()
}

# The inverse link and its derivative, taken from the family object when it has
# them -- which is the case for anything stats::binomial() and friends return,
# including a link built by the caller and passed as a "link-glm" object -- and
# resolved from the name otherwise.
link_functions <- function(family, link) {
  out <- family[c("linkfun", "linkinv", "mu.eta")]

  if (!is.function(out[["linkinv"]]) || !is.function(out[["linkfun"]])) {
    resolved <- try(stats::make.link(link), silent = TRUE)

    if (inherits(resolved, "try-error")) {
      arg::err("the {.val {link}} link is not one {.fn stats::make.link} knows.
                Supply it as a link object, as in
                {.code binomial(link = my_link)}, so that the inverse link
                comes with it")
    }

    out <- resolved[c("linkfun", "linkinv", "mu.eta")]
  }

  if (!is.function(out[["mu.eta"]])) {
    out[["mu.eta"]] <- NULL
  }

  out
}

# Build the map from the caller's additive predictor to the scale the engine's
# family works on, together with its derivative when the link brought one.
compose_link <- function(custom_link, native) {
  linkinv <- custom_link[["linkinv"]]
  mu_eta <- custom_link[["mu.eta"]]

  # The warnings are suppressed rather than avoided, and the values are left
  # exactly as they were. A link the caller supplied maps an unconstrained
  # predictor to a mean that the engine's own scale constrains -- the inverse
  # link of `Gamma()` is the standard case, where a negative predictor gives a
  # negative mean and `log()` of it is NaN -- and the sampler already treats a
  # NaN density as an impossible proposal and rejects it. So the proposal is
  # handled correctly and the only thing the warning does is emit one line per
  # rejected proposal, which on a default `Gamma()` fit is dozens of them.
  theta <- switch(native,
    identity = function(eta) linkinv(eta),
    log = function(eta) suppressWarnings(log(linkinv(eta))),
    logit = function(eta) suppressWarnings(stats::qlogis(linkinv(eta))))

  dtheta <- {
    if (is_null(mu_eta)) NULL
    else switch(native,
      identity = function(eta) mu_eta(eta),
      log = function(eta) mu_eta(eta) / linkinv(eta),
      logit = function(eta) {
        mu <- linkinv(eta)
        mu_eta(eta) / (mu * (1 - mu))
      })
  }

  list(link_theta = theta, link_dtheta = dtheta)
}
