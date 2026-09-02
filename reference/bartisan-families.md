# Response families for generalized BART

[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
accepts the [stats::family](https://rdrr.io/r/stats/family.html) objects
used by [`stats::glm()`](https://rdrr.io/r/stats/glm.html), so
[`gaussian()`](https://rdrr.io/r/stats/family.html),
`binomial("probit")`, [`poisson()`](https://rdrr.io/r/stats/family.html)
and [`stats::Gamma()`](https://rdrr.io/r/stats/family.html) all work
unchanged. The functions documented here supply the additional families
that have no [`glm()`](https://rdrr.io/r/stats/glm.html) counterpart, in
the same style, so that they can be passed to the `family` argument the
same way.

One thing to know about the gamma family, because base R's default is
the wrong choice here:
[`stats::Gamma()`](https://rdrr.io/r/stats/family.html) defaults to
`link = "inverse"`, and the inverse link is the worst link for this
sampler. **Write `Gamma("log")`.** Base R's function is left as base R
defines it, so that attaching this package cannot change what
[`glm()`](https://rdrr.io/r/stats/glm.html) does; naming the family as a
string, `family = "Gamma"`, gets the log link, since that spelling is
this package's own.
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
says so when a link's inverse does not cover the additive predictor,
which is the case that catches it; see Details.

## Usage

``` r
negbin(link = "log", theta = NULL)

ordinal(link = "logit")

multinomial(link = "logit", reference = NULL, replicates = 200L)

dpm_aft(
  nu = 10,
  q = 0.95,
  k_s = 10,
  alpha = NULL,
  max_clusters = NULL,
  psi = 0.5
)

dpm(nu = 10, q = 0.95, k_s = 10, alpha = NULL, max_clusters = NULL, psi = 0.5)

weibull_aft()

loglogistic_aft()

lognormal_aft()

ph(num_bins = NULL, lambda_shape = 1, update_lambda = TRUE)

location_scale(link = "identity")

zi_poisson(link = "log")

zi_negbin(link = "log", theta = NULL)

Beta(link = "logit", phi = NULL)

ordbeta(link = "logit", phi = NULL)

custom_family(
  logdens,
  num_predictors = 1L,
  start = 0,
  derivatives = NULL,
  aux_names = NULL,
  aux_start = 0,
  name = "custom"
)
```

## Arguments

- link:

  the link function. Each family compiles the links for which the
  additive predictor is the natural unconstrained scale; any other link
  is applied from R, for the families where that is well defined. See
  Details.

- theta:

  for `negbin()` and `zi_negbin()`, a fixed value for the dispersion
  parameter. The default, `NULL`, draws it along with everything else.

- reference:

  for `multinomial()`, the response category to hold as the reference.
  With the logit link the default, `NULL`, fits one forest per category
  instead and leaves the model unidentified, which is what makes the
  prior symmetric in the categories; see Details. The probit link is
  always written as contrasts against a reference, so there the default
  is the first level.

- replicates:

  for `multinomial("probit")`, how many simulation draws to use for the
  category probabilities, which have no closed form. Larger is more
  accurate and slower.

- nu, q:

  for `dpm()`, the degrees of freedom of the baseline's
  inverse-chi-square prior on a component's variance and the quantile of
  that prior placed at a rough estimate of the residual standard
  deviation. The defaults, 10 and 0.95, are the paper's, and are tighter
  than BART's own 3 and 0.90 because the mixture covers small errors
  with extra components rather than with one component's left tail.

- k_s:

  for `dpm()`, how many units of the baseline's own scale the component
  means are allowed to reach out to. The default, 10, places the
  marginal of a component mean so that it reaches the largest residual
  of a linear fit.

- alpha:

  for `dpm()`, a fixed Dirichlet process concentration. The default,
  `NULL`, draws it.

- max_clusters, psi:

  for `dpm()`, the largest number of mixture components thought
  plausible and the shape of the taper towards it, which together set
  the prior on `alpha`. The defaults are a tenth of the sample size and
  0.5.

- num_bins:

  *Advanced.* For `ph()`, how many pieces the baseline hazard has, with
  the edges at evenly spaced quantiles of the observed times. The
  default, `NULL`, uses about the cube root of the sample size, which is
  the order the Freedman-Diaconis rule gives for a histogram. **You
  should not need to set this**: the estimates are flat in it over a
  sixty-fold range, and it is here for checking that rather than for
  tuning. See Details.

- lambda_shape:

  for `ph()`, the shape of the gamma prior on each bin's baseline
  hazard. Its rate is drawn.

- update_lambda:

  for `ph()`, whether to draw the baseline hazards. `FALSE` holds them
  at their prior mean, which is for diagnosis rather than analysis.

- phi:

  for `Beta()` and `ordbeta()`, a fixed value for the beta precision.
  The default, `NULL`, draws it.

- logdens:

  for `custom_family()`, the log density. A function of the response and
  the additive predictors, `function(y, eta)`, where `y` is a numeric
  vector of length `n` and `eta` an `n` by `num_predictors` matrix,
  returning a numeric vector of length `n`. With nuisance parameters it
  takes a third argument, `function(y, eta, aux)`, where `aux` is a
  numeric vector of their current values. It is the log density of *one
  unit of prior weight*, so that `weights` behave as they do elsewhere,
  and terms free of `eta` may be dropped.

- num_predictors:

  for `custom_family()`, how many additive predictors the density has,
  that is, how many forests to fit.

- start:

  for `custom_family()`, the value each additive predictor starts at, in
  place of the intercept-only fit the compiled families use. One value
  or one per predictor.

- derivatives:

  for `custom_family()`, an optional `function(y, eta, h)` returning a
  list with elements `score` and `info`, the first derivative of
  `logdens` with respect to the `h`th predictor and minus its second
  derivative, each a vector of length `n`. The default, `NULL`, takes
  central differences of `logdens`. It covers the additive predictors
  only: a nuisance parameter is always differenced, which costs three
  calls per sweep rather than three per leaf.

- aux_names:

  for `custom_family()`, the names of the nuisance parameters to draw.
  Naming them is what declares them, because the names label the columns
  of `fit$aux` and are what
  [`summary()`](https://rdrr.io/r/base/summary.html) and `fit$rhat`
  report them under. The default, `NULL`, means none, unless `aux_start`
  is given, in which case they are named positionally.

- aux_start:

  for `custom_family()`, the value each nuisance parameter starts at.
  One value or one per parameter. The sampler will walk to the posterior
  from a poor start, so this need only be the right order of magnitude.

- name:

  for `custom_family()`, a label used when printing the fit.

## Value

A list of class `bartisan_family`, containing at least the elements
`family` and `link`. Objects of this class are recognized by
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
alongside ordinary [stats::family](https://rdrr.io/r/stats/family.html)
objects.

## Details

Every family reduces to a scalar additive predictor, or to several of
them, together with the first two derivatives of the log density with
respect to each. That is the whole interface the sampler needs, which is
why the set of available families is not restricted to the conditionally
conjugate ones.

[`vignette("families", package = "bartisan")`](https://ngreifer.github.io/bartisan/articles/families.md)
covers all of this at length: how to choose a family, how to choose
among the links a family offers, and what each family is and is not for.
What follows is the short version, and the points that could lead to
output being misread.

The supported families and links are:

|  |  |  |  |
|----|----|----|----|
| Family | Links | Additive predictors | Drawn nuisance parameters |
| [`gaussian()`](https://rdrr.io/r/stats/family.html) | `identity` | 1 | residual standard deviation |
| [`binomial()`](https://rdrr.io/r/stats/family.html) | `logit`, `probit`, `cloglog` | 1 | none |
| [`poisson()`](https://rdrr.io/r/stats/family.html) | `log` | 1 | none |
| `negbin()` | `log` | 1 | dispersion |
| `Gamma("log")` | `log` (also `inverse`, `identity`, any link) | 1 | shape |
| `ordinal()` | `logit`, `probit`, `cloglog` | 1 | cutpoints |
| `multinomial()` | `logit`, `probit` | one per category, or per non-reference level | latent covariance, for the probit link |
| `weibull_aft()`, `loglogistic_aft()`, `lognormal_aft()` | none | 1 | scale |
| `ph()` | none | 1 | baseline hazard per bin |
| `dpm_aft()` | none | 1 | error mixture, concentration |
| `location_scale()` | `identity` | 2 | none |
| `zi_poisson()` | `log` | 2 | none |
| `zi_negbin()` | `log` | 2 | dispersion |
| `Beta()` | `logit`, `probit`, `cloglog` | 1 | precision |
| `ordbeta()` | `logit` | 1 | 2 cutpoints, precision |
| `dpm()` | `identity` | 1 | error mixture, concentration |

A family with more than one additive predictor fits one forest per
predictor. Nuisance parameters are drawn alongside the trees and
reported in `fit$aux`.

## Links the engine does not compile

The links listed above are the ones the sampler evaluates in compiled
code. Any other link is accepted for
[`gaussian()`](https://rdrr.io/r/stats/family.html),
[`binomial()`](https://rdrr.io/r/stats/family.html),
[`poisson()`](https://rdrr.io/r/stats/family.html), `negbin()` and
[`Gamma()`](https://rdrr.io/r/stats/family.html), and applied from R by
composing the caller's inverse link with the family's own, with the
chain rule carrying the derivatives back. So `binomial("cauchit")`
works, as does any link object of the kind
[`stats::make.link()`](https://rdrr.io/r/stats/make.link.html) returns.
It costs a call into R for every leaf the sampler visits, and the leaf
prior scale is calibrated for the compiled link.

A link whose inverse has a restricted range – `Gamma("inverse")`,
`Gamma("identity")`, `poisson("identity")` – will give non-finite
densities for some predictors. Those proposals are rejected rather than
breaking the chain, but they are wasted work and the fit is worse for
it, so
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
says so when it starts. Prefer links whose inverse is defined on the
whole line.

The families with more than one additive predictor, or whose link enters
somewhere other than a single mean – `ordinal()`, `multinomial()`, the
accelerated failure time families, `location_scale()`, the zero-inflated
families and `ordbeta()` – take only their listed links.
`custom_family()` is the way to reach anything else.

## If no family is given

`family` may be omitted, in which case it is read off the response: a
`Surv` object gets `dpm_aft()`, an ordered factor `ordinal()`, a
two-valued response [`binomial()`](https://rdrr.io/r/stats/family.html),
any other factor `multinomial()`, and any numeric response `dpm()`.
`dpm()` cannot take prior weights, so a weighted fit with no family
named is an error rather than a silent substitution. The choice is
reported with a message, which setting `family` silences. A count is
*not* given [`poisson()`](https://rdrr.io/r/stats/family.html) and a
numeric response taking two values other than 0 and 1 is *not* given
[`binomial()`](https://rdrr.io/r/stats/family.html), since either would
be a modeling decision rather than a reading of the response's type.

## What to know before reading the output

`ordinal()` accepts a numeric response as well as an ordered factor,
taking its sorted unique values as the categories, and that is **a
method rather than a fallback**: the cutpoints absorb the marginal
distribution of the response and the forest explains only the ordering,
so nothing is assumed about the error distribution and the model for
\\P(Y \le y \mid x)\\ is invariant to any monotone transformation of the
response. Predict with `type = "mean"`. Bin the response onto twenty-odd
quantiles first: one cutpoint per distinct value costs 73 seconds
against 2.7 for [`gaussian()`](https://rdrr.io/r/stats/family.html) at
1000 observations, and twenty-five bins was both sixteen times faster
and slightly more accurate. See the vignette.

`ordinal()` uses the cumulative-link parameterization of
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html) , in which
\\P(Y \le k) = F(c_k - \eta)\\, so larger values of the additive
predictor shift mass towards higher categories. Only the differences
\\c_k - \eta_i\\ are identified, so one location has to be pinned: with
three or more categories the draws are reported in the chart where **the
additive predictor has mean zero over the fitted sample and every
cutpoint is free**, which is the chart `polr()` reports in when its
predictors are centered. With exactly two categories the single boundary
is folded into the intercept instead, so a two-category response is
exactly binary regression with the matching link and on the same scale.
`cut1` is therefore a free parameter rather than a constant zero, which
is a change from earlier versions.

`multinomial()` by default fits one forest per category and leaves the
model unidentified, since adding any function of the predictors to every
category's forest leaves the probabilities alone. This is the
parameterization of Murray (2021), whose point is that the prior is then
symmetric in the categories; every identified quantity is recovered from
the draws. Passing `reference` instead pins that category at zero and
fits one fewer forest, giving log odds against it.

`multinomial("probit")` lets the latent utilities correlate, which a
multinomial logit cannot express at all. \\\Sigma\\ is normalized by the
trace constraint \\\mathrm{tr}(\Sigma) = C\\ (Burgette and Nordheim
2012), and its lower triangle appears in `fit$aux` as `sigma11`,
`sigma21` and so on. **Those correlations are weakly identified: read
the fitted probabilities, not the covariance.** They enter the
likelihood only through orthant probabilities of a distribution whose
location is a sum of trees, so a flexible mean absorbs much of the
dependence they are meant to measure. At 900 observations a true
correlation of zero came back as -0.57, and posterior intervals ran up
to 1.07 wide on a parameter confined to \\(-1, 1)\\; at 3000
observations the posterior tracks the truth to within about 0.2. Two
further consequences: the likelihood has **no closed form**, so it and
every category probability are simulated with `replicates` draws, and
`augment` does not apply, because the latent variables are the model
rather than a rewriting of it. The sampler is Algorithm P2 of Xu et al.
(2025).

`dpm()` is not a distribution but a way of not choosing one. It is
DPMBART (George et al. 2019): a numeric response with the sum of trees
for its mean, as [`gaussian()`](https://rdrr.io/r/stats/family.html)
has, and a Dirichlet process mixture of normals for its errors instead
of a single normal, so the error distribution comes out as whatever
mixture the data ask for.
[`error_density()`](https://ngreifer.github.io/bartisan/reference/error_density.md)
gives that density, which is the object the method exists to produce.

**It is the family to reach for by default on a numeric response**,
because it does not pay for its flexibility: on normal errors, where
[`gaussian()`](https://rdrr.io/r/stats/family.html) is exactly right, it
came out slightly ahead on both held-out error and log score at the same
time to one decimal place, and on heavy-tailed, skewed and bimodal
errors it was ahead by a great deal – on bimodal errors at a thousand
observations, 0.050 against 0.154 in held-out RMSE, a factor of three,
at the same time to a tenth of a second. So it is the family a numeric
response gets when none is named. The reasons to prefer
[`gaussian()`](https://rdrr.io/r/stats/family.html) are not statistical:
it takes **prior weights**, which `dpm()` refuses, and it reports one
interpretable `sigma` where `dpm()` has a mixture. It is also faster, by
1.4 times at a thousand observations. The vignette has the comparison.

Two things to know about `dpm()` itself. **It does not buy
heteroskedasticity** – the error distribution is flexible but it is the
same distribution at every \\x\\, and `location_scale()` is the family
for a spread that depends on the predictors. And **the additive
predictor is the conditional mean**, as it is for
[`gaussian()`](https://rdrr.io/r/stats/family.html): nothing in the
model forces the mixture to be centered, so the sampler works in a chart
where only the sum of the predictor and the error mean is identified,
but reporting is done in the chart where the mixture has mean zero and
the whole conditional mean sits on the predictor. `type = "link"` and
`type = "response"` therefore agree exactly, and `fit$aux` reports the
shift that was taken out as `center` rather than an error mean, which is
zero by construction. Prior weights are refused, since a weight would
have to be a multiplicity in the Dirichlet process.

The gamma family puts the forest on the log mean and draws the shape,
which acts as the inverse dispersion; it does *not* regress the shape on
the predictors. `negbin()` and `ordbeta()` take `theta` and `phi` to fix
their equivalents, and the gamma shape has no such argument because a
caller who knows it is rare. **The link is where the care is needed.**
Only `log` is compiled, and the base R default of `inverse` is the worst
case for this sampler: its inverse maps a negative predictor to a
negative mean, whose log is not a number, so the proposal is rejected –
dozens of times per fit. Measured on 600 observations and 50 trees,
[`stats::Gamma()`](https://rdrr.io/r/stats/family.html) took 7.2 seconds
against 3.8 for `Gamma("log")` and fitted the mean slightly worse. So
write the link, or name the family as the string `"Gamma"`, which
resolves to the log link. Any composed link whose inverse has a
restricted range is reported when the fit starts.

The accelerated failure time families expect a right-censored response,
supplied either as a
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) object
or as a two-column matrix of times and event indicators. They model
\\\log T = \eta + \sigma\epsilon\\ with \\\epsilon\\ standard Gumbel,
logistic or normal respectively, giving Weibull, log-logistic and
log-normal survival times, so the predictor is a log time ratio in each.

`ph()` is the proportional hazards alternative, with a
piecewise-constant baseline: \\\lambda(t \mid x) =
\lambda_0(t)\exp(r(x))\\, so its predictor is a log *hazard* ratio and
the baseline is free to take any shape rather than the monotone one a
Weibull imposes. `num_bins` sets how many pieces, with the edges at
evenly spaced quantiles of the observed times; the default is about the
cube root of the sample size. The bin hazards are drawn from their exact
gamma conditionals and reported as `lambda1`, `lambda2`, ... in
`fit$aux`, together with the rate of their own prior. The predictor and
the baseline are identified only jointly, so the baseline carries the
level and the predictor is reported centered on it.

Cox's *partial* likelihood is what cannot be used here: it couples
observations through risk sets and so does not decompose into a sum over
the observations reaching a leaf. The full likelihood of the
piecewise-exponential model does decompose, and it approaches the
partial likelihood as the bins shrink, which is how `ph()` reaches
proportional hazards without it.

**`num_bins` is not a modeling decision, and its default should be left
alone.** It is exposed for checking that, not for tuning. Measured over
three replicates at 700 observations, sweeping it from 4 to 250 – a
sixty-fold range, against a baseline hazard that turns over and against
a Weibull one – moved the error in the survival function between 0.035
and 0.048 and the error in the log hazard ratio between 0.135 and 0.185,
with no trend in either and every difference inside the
replicate-to-replicate spread. What the bin count does change is the
effective number of parameters, which grows with it: from 17 at four
bins to 206 at 250. That is what makes the default matter for `loo()`
and `waic()` rather than for the estimates – one parameter per event
time would leave each observation's density inflated by a parameter only
it informs, and leave-one-out unable to do its job.

The three differ in cost, though not enough to decide a model on.
`lognormal_aft()` and `loglogistic_aft()` impute each censored failure
time above its censoring time, which makes their targets quadratic and
is worth 8 to 30 times the speed; `weibull_aft()` needs no imputation
because its likelihood already has a form the sampler can collapse to a
single pass, but only under hard rules, which makes it the slowest of
the three at the default gate.

`dpm_aft()` is the accelerated failure time model with the error
distribution estimated rather than assumed: \\\log T = m(x) + W\\ with
\\W\\ a Dirichlet process mixture of normals constrained to mean zero,
and censored log-times imputed. It is `dpm()`'s error model with
censoring, so its predictor is the conditional mean of \\\log T\\,
[`error_density()`](https://ngreifer.github.io/bartisan/reference/error_density.md)
reports the fitted error density, and prior weights are refused for the
same reason `dpm()` refuses them. Following Henderson, Louis, Rosner and
Varadhan (2020).

Reach for it when the shape of the error is in doubt and you would
rather not assert one. Measured against a two-component error it was
worth 210 held-out log points and a third of the error in \\S(t \mid
x)\\ over the best fixed-error family; against a log-normal error, where
`lognormal_aft()` is correctly specified, the two were within 0.1 log
points of each other. So it gains where the assumption would have been
wrong and costs nothing where it would have been right, which is the
property `dpm()` has against
[`gaussian()`](https://rdrr.io/r/stats/family.html).

`location_scale()` regresses the mean and the log standard deviation of
a normal response on separate forests, so the variance is an
unrestricted function of the predictors.

`zi_poisson()` and `zi_negbin()` are zero-inflated counts. Both parts
get their own forest: the first predictor is the log mean of the count
component and the second the log odds that an observation is a
structural zero, so the excess-zero mechanism is free to depend on the
predictors. The two are reported as the `count` and `zero` predictors.

`Beta()` is beta regression for a response strictly inside the unit
interval: a forest on the link of the mean, and a precision drawn
alongside it. A response *at* either endpoint has no beta density, so it
is an error rather than something to nudge inward.

`ordbeta()` is the ordered beta regression of Kubinec (2023), for a
response on the closed unit interval with point masses at zero and one.
One predictor drives both the probability of landing on an endpoint,
through a pair of cutpoints as in an ordinal model, and the mean of the
beta density in between. Because the predictor also enters the beta mean
it is identified, so unlike `ordinal()` both cutpoints are drawn.

Choose between them on whether the response can reach a boundary, not on
whether it happens to in the sample: the two ask different questions,
and `ordbeta()` fitted to a response with no boundary observations
leaves its cutpoints with nothing to identify them.

## Several additive predictors

Most families model one parameter with one forest. Some model several,
and then every argument that could mean something different for each of
them may be given once, to apply to all, or one per forest:
positionally, or keyed by the names below. This includes `formula`, so a
forest can have predictors of its own; see
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).

The first forest is always the main parameter, the one a single-forest
family would have on its own. The order is:

|  |  |
|----|----|
| Family | Forests, in order |
| `location_scale()` | `mean`, `log_sd` |
| `zi_poisson()`, `zi_negbin()` | `count`, `zero` |
| `custom_family(num_predictors = k)` | `eta1` ... `etak` |
| `multinomial()` | one per level, or per non-reference level |
| `mnp()` | one per non-reference level, named for its contrast |
| everything else | `eta` |

`mean` is the mean and `log_sd` is the logarithm of the standard
deviation, which is the scale the forest works on. `count` is the linear
predictor of the count component and `zero` that of the inflation
component. A custom family's nuisance parameters are not on this list:
they are carried as trailing forests pinned to a single leaf, and
nothing about them is set per forest.

So, for a location-scale model with a smaller scale forest and a
restricted set of predictors for it:

    bartisan(list(y ~ x1 + x2 + x3, log_sd = ~ x1), data = d,
             family = location_scale(), num_trees = c(mean = 50, log_sd = 10))

**The multinomial families are the exception.** Their forests are the
levels of one categorical parameter and act together rather than
describing separate components of the response distribution, so there is
nothing a caller could mean by giving one level a different prior or a
different set of predictors from another. Every argument applies to all
of their forests at once, and more than one value is an error rather
than a silent recycling.

## Supplying a likelihood

`custom_family()` takes the log density itself, as an R function, and
fits the model that goes with it. Nothing else about the sampler
changes: the leaf-level Laplace proposal needs the first two derivatives
of the log density with respect to each additive predictor and nothing
more, and central differences of the supplied function produce both.

    # A Poisson model written out by hand. Terms free of eta may be dropped;
    # they cancel from every acceptance ratio.
    pois <- custom_family(function(y, eta) y * eta[, 1] - exp(eta[, 1]),
                          start = log(mean(d$y)))

    # Two predictors: a mean and a log standard deviation.
    ls <- custom_family(function(y, eta) dnorm(y, eta[, 1], exp(eta[, 2]),
                                               log = TRUE),
                        num_predictors = 2, start = c(0, 0))

The function is called once per leaf per Fisher-scoring step with the
observations reaching that leaf, so it must be vectorized over `y` and
the rows of `eta`; it must not be vectorized *within* an observation,
and it must return exactly one value per row. Supplying `derivatives`
cuts three calls to one and removes the differencing error.

Nuisance parameters are drawn alongside the trees when `aux_names` names
them, and `logdens` then takes a third argument holding their current
values:

    # A Gaussian written out by hand, with its scale drawn rather than fixed.
    by_hand <- custom_family(
      logdens = function(y, eta, aux) dnorm(y, eta[, 1], exp(aux[1]), log = TRUE),
      aux_names = "log_sigma", aux_start = 0)

They are reported in `fit$aux` under those names, and covered by
[`summary()`](https://rdrr.io/r/base/summary.html) and `fit$rhat` like
any other family's. There is no prior argument and no bounds argument,
because a nuisance parameter here is carried as an additive predictor
whose forest is pinned at depth zero – one tree that can never split, so
the forest is a single scalar – and it is drawn by the same
Laplace-plus-Metropolis step as any leaf, under that step's Gaussian
leaf prior. So a parameter with a restricted range is handled the way it
would be for a real predictor, by writing the transform into `logdens`:
the [`exp()`](https://rdrr.io/r/base/Log.html) above is what keeps the
scale positive.

What `custom_family()` does not do: the response must be numeric, so a
factor has to be coded first; and since the package cannot know what the
mean of the density is, `predict(type = "response")` returns the
additive predictors rather than a fitted mean.

## References

Burgette, L. F., & Nordheim, E. V. (2012). The trace restriction: an
alternative identification strategy for the Bayesian multinomial probit
model. *Journal of Business & Economic Statistics*, 30(3), 404–410.
[doi:10.1080/07350015.2012.680416](https://doi.org/10.1080/07350015.2012.680416)

George, E., Laud, P., Logan, B., McCulloch, R., & Sparapani, R. (2019).
Fully nonparametric Bayesian additive regression trees. In *Topics in
Identification, Limited Dependent Variables, Partial Observability,
Experimentation, and Flexible Modeling: Part B* (Advances in
Econometrics, vol. 40B, pp. 89–110). Emerald Publishing.
[doi:10.1108/S0731-90532019000040B006](https://doi.org/10.1108/S0731-90532019000040B006)

Kubinec, R. (2023). Ordered beta regression: a parsimonious,
well-fitting model for continuous data with lower and upper bounds.
*Political Analysis*, 31(4), 519–536.
[doi:10.1017/pan.2022.20](https://doi.org/10.1017/pan.2022.20)

Murray, J. S. (2021). Log-linear Bayesian additive regression trees for
multinomial logistic and count regression models. *Journal of the
American Statistical Association*, 116(534), 756–769.
[doi:10.1080/01621459.2020.1813587](https://doi.org/10.1080/01621459.2020.1813587)

Xu, Y., Hogan, J., Daniels, M., Kantor, R., & Mwangi, A. (2025).
Augmentation samplers for multinomial probit Bayesian additive
regression trees. *Journal of Computational and Graphical Statistics*,
34(2), 498–508.
[doi:10.1080/10618600.2024.2388605](https://doi.org/10.1080/10618600.2024.2388605)

## See also

[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md),
[`error_density()`](https://ngreifer.github.io/bartisan/reference/error_density.md),
and
[`vignette("families", package = "bartisan")`](https://ngreifer.github.io/bartisan/articles/families.md)
for the long form.

## Examples

``` r
if (FALSE) {
bartisan(y ~ ., data = d, family = negbin())
bartisan(y ~ ., data = d, family = ordinal("probit"))
bartisan(survival::Surv(time, status) ~ ., data = d, family = weibull_aft())

# The same response under proportional hazards, with a free baseline.
bartisan(survival::Surv(time, status) ~ ., data = d, family = ph())

# An accelerated failure time model with the error distribution estimated.
bartisan(survival::Surv(time, status) ~ ., data = d, family = dpm_aft())
bartisan(count ~ ., data = d, family = zi_negbin())
bartisan(proportion ~ ., data = d, family = Beta())

# The same response, when it can also sit exactly at 0 or 1.
bartisan(proportion ~ ., data = d, family = ordbeta())

# The latent covariance is reported in `aux`, as its lower triangle.
fit <- bartisan(y ~ ., data = d, family = multinomial("probit"))
colMeans(fit$aux)

# A numeric response whose error distribution is estimated rather than
# assumed. `error_density()` reports what shape it came out.
fit <- bartisan(y ~ ., data = d, family = dpm())
error_density(fit)

# A link the engine does not compile, applied from R.
bartisan(y ~ ., data = d, family = binomial("cauchit"))

# A likelihood supplied from R.
bartisan(y ~ ., data = d,
        family = custom_family(function(y, eta) y * eta[, 1] - exp(eta[, 1])))
}
```
