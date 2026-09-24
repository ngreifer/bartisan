# Response families for generalized BART

[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
accepts the [`stats::family`](https://rdrr.io/r/stats/family.html)
objects used by [`stats::glm()`](https://rdrr.io/r/stats/glm.html), so
[`gaussian()`](https://rdrr.io/r/stats/family.html),
`binomial("probit")`, [`poisson()`](https://rdrr.io/r/stats/family.html)
and [`stats::Gamma()`](https://rdrr.io/r/stats/family.html) all work
unchanged. The functions documented here supply the additional families
that have no [`glm()`](https://rdrr.io/r/stats/glm.html) counterpart, in
the same style, so that they can be passed to the `family` argument the
same way.

One caveat carries over from base R:
[`stats::Gamma()`](https://rdrr.io/r/stats/family.html) defaults to
`link = "inverse"`, which is the worst link for this sampler, so write
`Gamma("log")`, or name the family as the string `family = "Gamma"`,
which is this package's own spelling and resolves to the log link; see
Details.

## Usage

``` r
negbin(link = "log", theta = NULL)

ordinal(link = "logit", cut_alpha = 1)

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

gaussian_ls(link = "identity")

Gamma_ls(link = "log")

zi_poisson(link = "log")

zi_negbin(link = "log", theta = NULL)

Beta(link = "logit", phi = NULL)

ordbeta(link = "logit", phi = NULL, cut_alpha = 1)

tweedie(link = "log", power = 1.5, phi = NULL)
```

## Arguments

- link:

  string; the link function. Allowable options include `"logit"` (the
  default), `"probit"`, and `"cloglog"` for `ordinal()` and `Beta()`,
  and `"logit"` (the default) and `"probit"` for `multinomial()`. The
  remaining families take one link each, which is therefore the default:
  `"log"` for `negbin()`, `zi_poisson()`, `zi_negbin()`, and
  `Gamma_ls()`, `"logit"` for `ordbeta()`, and `"identity"` for
  `gaussian_ls()`. Each family compiles the links for which the additive
  predictor is the natural unconstrained scale; any other link is
  applied from R, for the families where that is well defined. See
  Details.

- theta:

  `numeric`; for `negbin()` and `zi_negbin()`, a fixed value for the
  dispersion parameter, which must be positive. Default is `NULL` to
  draw it along with everything else.

- cut_alpha:

  `numeric`; for `ordinal()` and `ordbeta()`, the concentration of the
  induced-Dirichlet prior on the cutpoints. Default is 1, which says the
  category probabilities they imply are uniform over the simplex; larger
  values pull them toward equal shares. See Details.

- reference:

  for `multinomial()`, the response category to hold as the reference,
  given as a single value naming one of the response's levels. Default
  is `NULL`, which with the logit link fits one forest per category
  instead and leaves the model unidentified, which is what makes the
  prior symmetric in the categories; see Details. The probit link is
  always written as contrasts against a reference, so there the default
  is the first level.

- replicates:

  `numeric`; for `multinomial("probit")`, how many simulation draws to
  use for the category probabilities, which have no closed form. Default
  is 200. More is always better but resulting calculations will take
  longer.

- nu, q:

  `numeric`; for `dpm()` and `dpm_aft()`, the degrees of freedom of the
  baseline's inverse-chi-square prior on a component's variance and the
  quantile of that prior placed at a rough estimate of the residual
  standard deviation. Defaults are 10 and .95, which are tighter than
  BART's own 3 and .90 because the mixture covers small errors with
  extra components rather than with one component's left tail.

- k_s:

  `numeric`; for `dpm()` and `dpm_aft()`, how many units of the
  baseline's own scale the component means are allowed to reach out to.
  Must be positive. Default is 10, which places the marginal of a
  component mean so that it reaches the largest residual of a linear
  fit.

- alpha:

  `numeric`; for `dpm()` and `dpm_aft()`, a fixed Dirichlet process
  concentration, which must be positive. Default is `NULL` to draw it.

- max_clusters, psi:

  `numeric`; for `dpm()` and `dpm_aft()`, the largest number of mixture
  components thought plausible (2 or greater) and the shape of the taper
  toward it, which together set the prior on `alpha`. Defaults are
  `NULL` to use a tenth of the sample size, and .5.

- num_bins:

  *Advanced.* `numeric`; for `ph()`, how many pieces the baseline hazard
  has, with the edges at evenly spaced quantiles of the observed times.
  Must be 2 or greater. Default is `NULL` to use about the cube root of
  the sample size, which is the order the Freedman-Diaconis rule gives
  for a histogram. This should not need to be set: the estimates are
  flat in it over a sixty-fold range, and it is here for checking that
  rather than for tuning. See Details.

- lambda_shape:

  `numeric`; for `ph()`, the shape of the gamma prior on each bin's
  baseline hazard, which must be positive. Its rate is drawn. Default is
  1.

- update_lambda:

  `logical`; for `ph()`, whether to draw the baseline hazards. Default
  is `TRUE`; `FALSE` holds them at their prior mean, which is for
  diagnosis rather than analysis.

- phi:

  `numeric`; for `Beta()` and `ordbeta()`, a fixed value for the beta
  precision, and for `tweedie()` a fixed value for the dispersion. Must
  be positive. Default is `NULL` to draw it.

- power:

  `numeric`; for `tweedie()`, the variance power \\p\\ in
  \\\mathrm{Var}(y) = \phi\mu^p\\, strictly between 1 and 2. Default is
  1.5, and the value is held fixed there rather than drawn, because the
  power is weakly identified from data of the sizes this package is used
  on and a badly determined power drags the dispersion around with it.
  Pass `NULL` to draw it, which is worth doing only with a large sample
  and a real interest in the shape rather than the mean.

## Value

A `<bartisan_family>` object, which is a list containing at least the
elements `family` and `link` and which inherits from `family`, so that
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
recognizes it wherever it recognizes an ordinary
[`stats::family`](https://rdrr.io/r/stats/family.html) object.

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
| `Gamma("log")` | `log`; any other link is ignored | 1 | shape |
| `ordinal()` | `logit`, `probit`, `cloglog` | 1 | cutpoints |
| `multinomial()` | `logit`, `probit` | one per category, or per non-reference level | latent covariance, for the probit link |
| `weibull_aft()`, `loglogistic_aft()`, `lognormal_aft()` | none | 1 | scale |
| `ph()` | none | 1 | baseline hazard per bin |
| `dpm_aft()` | none | 1 | error mixture, concentration |
| `gaussian_ls()` | `identity` | 2 | none |
| `Gamma_ls()` | `log` | 2 | none |
| `zi_poisson()` | `log` | 2 | none |
| `zi_negbin()` | `log` | 2 | dispersion |
| `Beta()` | `logit`, `probit`, `cloglog` | 1 | precision |
| `ordbeta()` | `logit` | 1 | 2 cutpoints, precision |
| `tweedie()` | `log` | 1 | dispersion, and the power if it is drawn |
| `dpm()` | `identity` | 1 | error mixture, concentration |

A family with more than one additive predictor fits one forest per
predictor. Nuisance parameters are drawn alongside the trees and
reported in `fit$aux`.

### Setting `link`

The links listed above are the ones the sampler evaluates in compiled
code. Any other link is accepted for
[`gaussian()`](https://rdrr.io/r/stats/family.html),
[`binomial()`](https://rdrr.io/r/stats/family.html),
[`poisson()`](https://rdrr.io/r/stats/family.html) and `Beta()`, and
applied from R by composing the caller's inverse link with the family's
own, with the chain rule carrying the derivatives back. So
`binomial("cauchit")` works, as does any link object of the kind
[`stats::make.link()`](https://rdrr.io/r/stats/make.link.html) returns.
It costs a call into R for every leaf the sampler visits, and the leaf
prior scale is calibrated for the compiled link.

Some single-predictor families are exceptions. `negbin()` takes `"log"`
alone, so a link given to it is an error rather than a composition.
[`stats::Gamma()`](https://rdrr.io/r/stats/family.html) accepts any link
and fits none of them but `"log"`: every other link is dropped with a
message, because the ones base R offers have inverses that go
non-positive over part of the line and the additive predictor is
unconstrained.
[`custom_family()`](https://ngreifer.github.io/bartisan/reference/custom_family.md)
is the route to a gamma response on another link.

A composed link whose inverse has a restricted range
(`poisson("identity")`, `poisson("sqrt")`) will give non-finite
densities for some predictors. Those proposals are rejected rather than
breaking the chain, but they are wasted work and the fit is worse for
it, so
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
says so when it starts. Prefer links whose inverse is defined on the
whole line.

The families with more than one additive predictor, or whose link enters
somewhere other than a single mean (`ordinal()`, `multinomial()`, the
accelerated failure time families, `gaussian_ls()`, the zero-inflated
families and `ordbeta()`), take only their listed links.
[`custom_family()`](https://ngreifer.github.io/bartisan/reference/custom_family.md)
is the way to reach anything else.

### Omitting `family`

`family` may be omitted, in which case it is read off the response's
type and the choice reported with a message;
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
tabulates the lookup, since `family` is its argument. A count is read as
[`gaussian()`](https://rdrr.io/r/stats/family.html) rather than
[`poisson()`](https://rdrr.io/r/stats/family.html), and a numeric
response taking two values other than 0 and 1 as
[`gaussian()`](https://rdrr.io/r/stats/family.html) rather than
[`binomial()`](https://rdrr.io/r/stats/family.html), both being modeling
decisions rather than readings of the response's type.

### Reading the Output

`ordinal()` accepts a numeric response as well as an ordered factor,
taking its sorted unique values as the categories, and that is a method
rather than a fallback: the cutpoints absorb the marginal distribution
of the response and the forest explains only the ordering, so nothing is
assumed about the error distribution and the model for \\P(Y \le y \mid
x)\\ is invariant to any monotone transformation of the response.
Predict with `type = "mean"`. Bin the response onto twenty-odd quantiles
first, since one cutpoint per distinct value is both far slower and
slightly less accurate than twenty-five bins;
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
works it through.

The cutpoints of `ordinal()` and `ordbeta()` carry an induced-Dirichlet
prior, which puts the prior on the category probabilities the cutpoints
imply rather than on the cutpoints themselves. `cut_alpha` is its
concentration, and it regularizes a thinly observed category without
disturbing a well observed one.

A level of an ordered factor that nobody selected is kept rather than
dropped, so a rating scale with an unused point is fitted on all of its
points. A numeric response is read as the distinct values it takes, so
such a scale has to be given as an ordered factor for the unused point
to be modeled, and a warning says so.
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
derives the prior and works both through.

`ordinal()` uses the cumulative-link parameterization of
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html) , in which
\\P(Y \le k) = F(c_k - \eta)\\, so larger values of the additive
predictor shift mass toward higher categories. Only the differences
\\c_k - \eta_i\\ are identified, so one location has to be pinned: with
three or more categories the draws are reported in the chart where the
additive predictor has mean zero over the fitted sample and every
cutpoint is free, which is the chart `polr()` reports in when its
predictors are centered. With exactly two categories the single boundary
is folded into the intercept instead, so a two-category response is
exactly binary regression with the matching link and on the same scale.
`cut1` is therefore a free parameter rather than a constant zero, which
is a change from earlier versions.

`multinomial()` by default fits one forest per category and leaves the
model unidentified, since adding any function of the predictors to every
category's forest leaves the probabilities alone. The point of that
parameterization is that the prior is then symmetric in the categories,
and every identified quantity is recovered from the draws. Passing
`reference` instead pins that category at zero and fits one fewer
forest, giving log odds against it.

`multinomial("probit")` lets the latent utilities correlate, which a
multinomial logit cannot express at all. \\\Sigma\\ is normalized by the
trace constraint \\\mathrm{tr}(\Sigma) = C\\, and its lower triangle
appears in `fit$aux` as `sigma11`, `sigma21` and so on. Read the fitted
probabilities rather than the covariance: the correlations enter the
likelihood only through orthant probabilities of a distribution whose
location is a sum of trees, so a flexible mean absorbs much of the
dependence they are meant to measure and they are weakly identified
until the sample is large. The likelihood has no closed form, so it and
every category probability are simulated with `replicates` draws, and
`augment` does not apply, the latent variables being the model rather
than a rewriting of it.

`dpm()` is not a distribution but a way of not choosing one. It is
DPMBART: a numeric response with the sum of trees for its mean, as
[`gaussian()`](https://rdrr.io/r/stats/family.html) has, and a Dirichlet
process mixture of normals for its errors instead of a single normal, so
the error distribution comes out as whatever mixture the data ask for.
[`error_density()`](https://ngreifer.github.io/bartisan/reference/error_density.md)
gives that density, which is the object the method exists to produce.

It is the family to reach for by default on a numeric response, because
it does not pay for its flexibility: on normal errors, where
[`gaussian()`](https://rdrr.io/r/stats/family.html) is exactly right, it
came out slightly ahead on both held-out error and log score at the same
time to one decimal place, and on heavy-tailed, skewed and bimodal
errors it was ahead by a great deal (on bimodal errors at a thousand
observations, 0.050 against 0.154 in held-out RMSE, a factor of three,
at the same time to a tenth of a second). So it is the family a numeric
response gets when none is named. The reasons to prefer
[`gaussian()`](https://rdrr.io/r/stats/family.html) are not statistical:
it takes prior weights, which `dpm()` refuses, and it reports one
interpretable `sigma` where `dpm()` has a mixture. It is also faster, by
1.4 times at a thousand observations. The vignette has the comparison.

Some properties of `dpm()` itself are worth knowing. It does not buy
heteroskedasticity: the error distribution is flexible but it is the
same distribution at every \\x\\, and `gaussian_ls()` is the family for
a spread that depends on the predictors. And the additive predictor is
the conditional mean, as it is for
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
caller who knows it is rare. The link is where the care is needed. Only
`log` is compiled, and the base R default of `inverse` is the worst case
for this sampler: its inverse maps a negative predictor to a negative
mean, whose log is not a number, so the proposal is rejected, dozens of
times per fit, and the fit is both slower and slightly worse for it. So
write the link, or name the family as the string `"Gamma"`, which
resolves to the log link; base R's own function is left as base R
defines it, so that attaching bartisan cannot change what
[`glm()`](https://rdrr.io/r/stats/glm.html) does. Any composed link
whose inverse has a restricted range is reported when the fit starts.

The accelerated failure time families expect a right-censored response,
supplied either as a
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) object
or as a two-column matrix of times and event indicators. They model
\\\log T = \eta + \sigma\epsilon\\ with \\\epsilon\\ a standard smallest
extreme value, logistic or normal variate respectively, giving Weibull,
log-logistic and log-normal survival times.

A contrast in the predictor is a log time ratio in all of them, and in
`dpm_aft()` too, because with \\\epsilon\\ independent of \\x\\ every
quantile and both means of \\T\\ scale by \\e^{\Delta\eta}\\ whatever
shape the error has. What differs is what \\e^{\eta}\\ is on its own,
each family pinning its error's location differently: the median of
\\T\\ for `loglogistic_aft()` and `lognormal_aft()`, the geometric mean
for `dpm_aft()`, and the Weibull scale for `weibull_aft()`. Contrasts
are unaffected by any of that;
[`vignette("survival")`](https://ngreifer.github.io/bartisan/articles/survival.md)
tabulates the levels and measures the difference between them.

`weibull_aft()` is also the one family whose predictor carries a log
*hazard* ratio, of \\-\Delta\eta/\sigma\\, alongside its log time ratio.
That is a property of the smallest extreme value error rather than of
the structural part it shares with the other three: it is the only error
making an accelerated failure time model a proportional hazards model as
well.

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

`num_bins` should be left at its default. The estimates are insensitive
to it, and what it does change is the effective number of parameters,
which grows with the bin count and so matters for
[`loo()`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
and `waic()`.
[`vignette("survival")`](https://ngreifer.github.io/bartisan/articles/survival.md)
sweeps it.

The three differ in cost, though not enough to decide a model on.
`lognormal_aft()` and `loglogistic_aft()` impute each censored failure
time above its censoring time, which makes their targets quadratic and
is worth a large multiple of the speed; `weibull_aft()` needs no
imputation because its likelihood already has a form the sampler can
collapse to a single pass, but only under hard rules, which makes it the
slowest of the three at the default gate.

`dpm_aft()` is the accelerated failure time model with the error
distribution estimated rather than assumed: \\\log T = m(x) + W\\ with
\\W\\ a Dirichlet process mixture of normals constrained to mean zero,
and censored log-times imputed. It is `dpm()`'s error model with
censoring, so its predictor is the conditional mean of \\\log T\\ and a
contrast in it is a log time ratio exactly as for the three parametric
families above,
[`error_density()`](https://ngreifer.github.io/bartisan/reference/error_density.md)
reports the fitted error density, and prior weights are refused for the
same reason `dpm()` refuses them.

Reach for it when the shape of the error is in doubt and there is no
reason to assert one. It gains substantially where a fixed-error family
would have been wrong and costs next to nothing where one would have
been right, which is the property `dpm()` has against
[`gaussian()`](https://rdrr.io/r/stats/family.html);
[`vignette("survival")`](https://ngreifer.github.io/bartisan/articles/survival.md)
measures both cases.

`gaussian_ls()` regresses the mean and the log standard deviation of a
normal response on separate forests, so the variance is an unrestricted
function of the predictors.

`Gamma_ls()` does the same for a gamma response: the first forest is the
log mean, exactly as `Gamma("log")`'s is, and the second is the log
dispersion, so the shape is `exp(-log_dispersion)` at each observation
rather than one value drawn for the whole sample. What it relaxes is the
assumption that the coefficient of variation is constant, which is what
a gamma with a single shape asserts, and giving its second forest an
intercept-only formula puts that assumption back; see "Several additive
predictors" above.

Note that the mean forest is quadratic in neither family, but it takes
the same cheap exponential form `Gamma("log")`'s does, while the
dispersion forest takes the general path, so the cost sits in the second
forest and it is worth giving that one fewer trees.

`zi_poisson()` and `zi_negbin()` are zero-inflated counts. Both parts
get their own forest: the first predictor is the log mean of the count
component and the second the log odds that an observation is a
structural zero, so the excess-zero mechanism is free to depend on the
predictors. The two are reported as the `count` and `zero` predictors.

`Beta()` is beta regression for a response strictly inside the unit
interval: a forest on the link of the mean, and a precision drawn
alongside it. A response *at* either endpoint has no beta density, so it
is an error rather than something to nudge inward.

`ordbeta()` is ordered beta regression, for a response on the closed
unit interval with point masses at zero and one. One predictor drives
both the probability of landing on an endpoint, through a pair of
cutpoints as in an ordinal model, and the mean of the beta density in
between. Because the predictor also enters the beta mean it is
identified, so unlike `ordinal()` both cutpoints are drawn.

Choose between them on whether the response can reach a boundary, not on
whether it happens to in the sample: the two ask different questions,
and `ordbeta()` fitted to a response with no boundary observations
leaves its cutpoints with nothing to identify them.

`tweedie()` is the compound Poisson-gamma, for a non-negative response
with a point mass at zero and a continuous positive part, which is the
shape of spending, rainfall, insurance claims and earnings. It is the
analogue of `ordbeta()` at the other end: one predictor again drives
both parts, but through the mean rather than through a cutpoint, since
\\\mu = \exp(\eta)\\ and \\\mathrm{Var}(y) = \phi\mu^p\\ together fix
the probability of a zero at \\\exp(-\mu^{2-p}/(\phi(2-p)))\\. That is
what makes it a single process and is also its restriction: the share of
zeros has no level of its own, so a response whose zeros are more or
less common than its mean implies wants a two-part model instead, which
`zi_poisson()` and `zi_negbin()` are for counts and which
[`custom_family()`](https://ngreifer.github.io/bartisan/reference/custom_family.md)
can supply for anything else.

Because the mean is exactly \\\exp(\eta)\\, a counterfactual mean
through
[marginaleffects](https://CRAN.R-project.org/package=marginaleffects)
needs nothing beyond the forest, unlike a two-part model where it has to
be recombined across predictors; and the fit is comparable with a
[`poisson()`](https://rdrr.io/r/stats/family.html) or `Gamma("log")` fit
of the same response, since all three put the same quantity on the same
scale.

### Several Additive Predictors

Most families model one parameter with one forest. Some model several,
and then every argument that could mean something different for each of
them may be given per forest, keyed by the names below or positionally;
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
states the recycling rule and lists which arguments it covers, and
`formula` is among them, so a forest can have predictors of its own.

The first forest is always the main parameter, the one a single-forest
family would have on its own. This table is the canonical list of the
names, which
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
reproduces:

|  |  |
|----|----|
| Family | Forests, in order |
| `gaussian_ls()` | `mean`, `log_sd` |
| `Gamma_ls()` | `mean`, `log_dispersion` |
| `zi_poisson()`, `zi_negbin()` | `count`, `zero` |
| `custom_family(num_predictors = k)` | `eta1` ... `etak` |
| `multinomial("logit")` | one per level, or per non-reference level |
| `multinomial("probit")` | one per non-reference level, named for its contrast |
| everything else | `eta` |

`mean` is the mean and `log_sd` is the logarithm of the standard
deviation, which is the scale the forest works on. `Gamma_ls()`'s `mean`
is the log mean, as `Gamma("log")`'s single forest is, and its
`log_dispersion` is the logarithm of the dispersion, so that in both
location-scale families a larger second predictor means more spread.
`count` is the linear predictor of the count component and `zero` that
of the inflation component. A custom family's nuisance parameters are
not on this list: they are carried as trailing forests pinned to a
single leaf, and nothing about them is set per forest.

A forest whose formula names no predictor is a constant. `~ 1` leaves
that forest nothing to split on, so every tree in it is a stump and the
parameter is one drawn scalar. Every family here that takes more than
one formula accepts that, which is what makes the distinction between a
nuisance parameter and an empty forest a thin one: `gaussian_ls()` with
`~ 1` on its scale is
[`gaussian()`](https://rdrr.io/r/stats/family.html), `Gamma_ls()` with
`~ 1` is `Gamma("log")`, and `zi_poisson()` with `~ 1` on its inflation
part is the zero-inflated Poisson with a single structural-zero
probability. Note that the scalar is drawn under the leaf prior rather
than under the prior the corresponding built-in family puts on its
nuisance parameter, so the two agree closely rather than exactly. The
multinomial families are the exception, for the reason given below.

So, for a location-scale model with a smaller scale forest and a
restricted set of predictors for it:

    bartisan(list(y ~ x1 + x2 + x3, log_sd = ~ x1), data = d,
             family = gaussian_ls(), num_trees = c(mean = 50, log_sd = 10))

The multinomial families are the exception. Their forests are the levels
of one categorical parameter and act together rather than describing
separate components of the response distribution, so there is nothing a
caller could mean by giving one level a different prior or a different
set of predictors from another. Every argument applies to all of their
forests at once, and more than one value is an error rather than a
silent recycling.

### Supplying a Likelihood

When none of the families above is the right one,
[`custom_family()`](https://ngreifer.github.io/bartisan/reference/custom_family.md)
takes the log density itself, as an R function, and fits the model that
goes with it. Its page has the details, including what a log density
cannot supply: a posterior predictive distribution and a fitted mean.

## See also

[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
for fitting a model with one of these families;
[`custom_family()`](https://ngreifer.github.io/bartisan/reference/custom_family.md)
for a likelihood none of them covers;
[`error_density()`](https://ngreifer.github.io/bartisan/reference/error_density.md)
for the error distribution a `dpm()` or `dpm_aft()` fit estimates;
[`vignette("families", package = "bartisan")`](https://ngreifer.github.io/bartisan/articles/families.md)
for the long form

## Examples

``` r
data("rhc")
set.seed(123)

# A right-censored response, given as the time and the event indicator,
# with the error distribution estimated rather than assumed
fit <- bartisan(cbind(days, death) ~ ., data = rhc, family = dpm_aft(),
                num_trees = 10, num_burn = 50, num_draws = 50)

# The shape the errors came out
head(error_density(fit))
#>          at         mean        lower        upper
#> 1 -8.149747 6.990040e-05 3.117357e-05 0.0001222741
#> 2 -8.068250 8.174272e-05 3.713237e-05 0.0001414876
#> 3 -7.986752 9.544615e-05 4.415258e-05 0.0001634804
#> 4 -7.905255 1.112774e-04 5.240790e-05 0.0001886150
#> 5 -7.823757 1.295371e-04 6.209755e-05 0.0002172952
#> 6 -7.742260 1.505634e-04 7.344957e-05 0.0002499697

# The same response under proportional hazards, whose predictor is a log
# hazard ratio and whose baseline is free to take any shape
bartisan(cbind(days, death) ~ ., data = rhc, family = ph(),
         num_trees = 10, num_burn = 50, num_draws = 50)
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = cbind(days, death) ~ ., data = rhc, family = ph(), 
#>     num_trees = 10, num_burn = 50, num_draws = 50)
#> 
#> Family: "ph" with the "log" link
#> Observations: 1500
#> Structure: 1 forest of 10 trees, soft decision rules
#> Draws: 50 kept after 50 warmup
#> 
#> Posterior means: lambda1 = 0.0105, lambda2 = 0.0158, lambda3 = 0.0107, lambda4 = 0.0071, lambda5 = 0.0029, lambda6 = 0.00137, lambda7 = 0.00101, lambda8 = 0.000722, lambda9 = 0.00135, lambda10 = 0.00203, lambda11 = 0.00164, lambda12 = 0.00185, lambda_rate = 233

# An unordered response, with one forest per category and a prior that is
# symmetric in them
bartisan(race ~ . - days - death, data = rhc, family = multinomial(),
         num_trees = 10, num_burn = 50, num_draws = 50)
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = race ~ . - days - death, data = rhc, family = multinomial(), 
#>     num_trees = 10, num_burn = 50, num_draws = 50)
#> 
#> Family: "multinomial" with the "logit" link
#> Observations: 1500
#> Structure: 3 forests of 10 trees, soft decision rules
#> Draws: 50 kept after 50 warmup

# A link the engine does not compile, applied from R
bartisan(death ~ . - days, data = rhc, family = binomial("cauchit"),
         num_trees = 10, num_burn = 50, num_draws = 50)
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = death ~ . - days, data = rhc, family = binomial("cauchit"), 
#>     num_trees = 10, num_burn = 50, num_draws = 50)
#> 
#> Family: "binomial" with the "cauchit" link (supplied from R)
#> Observations: 1500
#> Structure: 1 forest of 10 trees, soft decision rules
#> Draws: 50 kept after 50 warmup
```
