# Response families in bartisan

``` r

library(bartisan)

set.seed(2026)

# Small chains throughout, so that this vignette builds quickly. The defaults
# are 50 trees and 500 draws after 500 warmup iterations.
ctrl <- bartisan_control(num_trees = 10, num_burn = 150, num_draws = 150)
```

This vignette is for choosing a `family`. It covers what each one
assumes, when to prefer it over the alternatives, and the handful of
places where the choice has a consequence that is easy to miss.
[`?bartisan_control`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
covers how the sampler goes about fitting them.

`family` works the way it does in
[`glm()`](https://rdrr.io/r/stats/glm.html): the
[`stats::family`](https://rdrr.io/r/stats/family.html) objects work
unchanged, *bartisan* adds the families that have no
[`glm()`](https://rdrr.io/r/stats/glm.html) counterpart in the same
style, and
[`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
takes a log density written as an R function when nothing on the list
fits. Unlike ordinary BART, nothing here requires the response to be
conditionally conjugate ([Linero 2025](#ref-linero2025)), which is why
the list is as long as it is.

## The families

| Family | Links | Additive predictors | Drawn nuisance parameters |
|----|----|----|----|
| [`gaussian()`](https://rdrr.io/r/stats/family.html) | `identity` | 1 | residual standard deviation |
| [`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | `identity` | 1 | error mixture, concentration |
| [`location_scale()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | `identity` | 2 | none |
| `Gamma("log")` | `log` | 1 | shape |
| [`binomial()`](https://rdrr.io/r/stats/family.html) | `logit`, `probit`, `cloglog` | 1 | none |
| [`poisson()`](https://rdrr.io/r/stats/family.html) | `log` | 1 | none |
| [`negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | `log` | 1 | dispersion |
| [`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | `log` | 2 | none |
| [`zi_negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | `log` | 2 | dispersion |
| [`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | `logit`, `probit`, `cloglog` | 1 | cutpoints |
| [`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | `logit`, `probit` | one per category, or per non-reference level | latent covariance, for the probit link |
| [`Beta()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | `logit`, `probit`, `cloglog` | 1 | precision |
| [`ordbeta()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | `logit` | 1 | 2 cutpoints, precision |
| [`weibull_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), [`loglogistic_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), [`lognormal_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | none | 1 | scale |
| [`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | none | 1 | baseline hazard per bin |
| [`dpm_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | none | 1 | error mixture, concentration |
| [`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | your own | as many as you ask for | none |

A family with more than one additive predictor fits one forest per
predictor, and `predict(type = "link")` returns one column per forest.
Nuisance parameters are drawn alongside the trees and reported in
`fit$aux`.

## Choosing a family

Start from the shape of the response:

| The response is | Use |
|----|----|
| numeric, unbounded | [`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), the default; [`gaussian()`](https://rdrr.io/r/stats/family.html) if you need prior weights; [`location_scale()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) if the spread varies with the predictors |
| positive and continuous | [`gaussian()`](https://rdrr.io/r/stats/family.html) or [`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) on `log(y)`, or `Gamma("log")` if the mean on the original scale is the quantity of interest |
| a proportion strictly inside \\(0, 1)\\ | [`Beta()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) |
| a proportion that can equal 0 or 1 | [`ordbeta()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) |
| a 0/1 indicator, or a two-column matrix of successes and failures | [`binomial()`](https://rdrr.io/r/stats/family.html) |
| a count | [`poisson()`](https://rdrr.io/r/stats/family.html), [`negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) if overdispersed, or a zero-inflated family if there are more zeros than either can produce |
| an ordered factor | [`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) |
| an unordered factor | [`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) |
| a time with censoring | [`weibull_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md), [`loglogistic_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) or [`lognormal_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md); [`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) for proportional hazards with a free baseline; [`dpm_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) if the shape of the error is in doubt |
| something else | [`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) |

Two questions cut across that table.

**Does anything besides the mean vary with the predictors?** Every
family except
[`location_scale()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
and the zero-inflated pair puts a forest on one location parameter and
holds the rest of the distribution fixed across observations. If the
spread itself moves with \\x\\, that is the wrong assumption and
[`location_scale()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is what relaxes it.

**Is the extra structure real, or would a flexible mean absorb it?** A
nonparametric mean makes this sharper than it is in a GLM. A sum of
trees can produce excess zeros on its own, by driving a Poisson mean
very low where the zeros are, so a zero-inflated family is for when the
zero mechanism is a *separate process you want to model* — not merely
when a histogram spikes at zero. The same caution applies to the
multinomial probit’s latent correlations.

When two families are both defensible, compare them instead of arguing
about them. [`loo()`](https://mc-stan.org/loo/reference/loo.html) and
[`waic()`](https://mc-stan.org/loo/reference/waic.html) work on a fit:

``` r

n <- 300
d <- data.frame(x1 = runif(n), x2 = runif(n))
d$count <- rpois(n, exp(1.2 * sin(pi * d$x1) + 0.4))

loo::loo_compare(
  loo::loo(bartisan(count ~ ., d, family = poisson(), control = ctrl)),
  loo::loo(bartisan(count ~ ., d, family = negbin(), control = ctrl))
)
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
#>   model elpd_diff se_diff p_worse       diag_diff       diag_elpd
#>  model1       0.0     0.0      NA                 2 k_psis > 0.54
#>  model2      -0.3     0.8    0.63 |elpd_diff| < 4
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
```

## The default family

`family` may be omitted, in which case it is read off the response:

| The response is | Family used |
|----|----|
| a `Surv` object, or a numeric matrix of non-negative times and 0/1 events | [`dpm_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) |
| an ordered factor | [`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) |
| logical, a factor or character with two levels, or numeric taking only the values 0 and 1 | [`binomial()`](https://rdrr.io/r/stats/family.html) |
| any other factor or character | [`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) |
| a two-column numeric matrix | [`binomial()`](https://rdrr.io/r/stats/family.html), read as successes and failures |
| any other numeric | [`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) |

The choice is reported when it is made, and naming `family` silences the
message. That is the intended way to silence it: an inferred family is a
modeling decision taken on your behalf, so confirming it in the call is
a better outcome than suppressing the notice.

``` r

d$binary <- rbinom(n, 1, 0.4)

fit <- bartisan(binary ~ x1 + x2, data = d, control = ctrl)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` to choose another, which also silences this message.
```

Two boundaries are deliberate. A numeric response taking exactly two
values that are *not* 0 and 1 is not read as binomial, because deciding
that `c(1, 2)` means failure and success would be a guess about which
value is the success. And a count is not read as Poisson, because “the
variance equals the mean” is a substantive claim rather than a reading
of the response’s type.

The number of distinct values does not enter into it: every numeric
response gets
[`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md).
A response with only a handful of distinct values probably wants
[`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
but that is a modeling decision, and a threshold would make the default
arbitrary and hard to predict.

**Neither
[`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
nor
[`dpm_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
can take prior weights**, so a fit with weights and no family named is
an error rather than a silent substitution. Dropping the weights and
changing family are both defensible and only you can say which was
meant. [`gaussian()`](https://rdrr.io/r/stats/family.html),
[`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
and
[`location_scale()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
take weights for a numeric response; the three `*_aft()` families and
[`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
take them for a censored one.

## Numeric responses

Four families fit a numeric response by putting a forest on its mean.
They differ in what else they allow to vary:

- **[`gaussian()`](https://rdrr.io/r/stats/family.html)** assumes one
  normal error, with `sigma` drawn and reported.
- **[`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)**
  estimates the error distribution — a Dirichlet process mixture of
  normals ([George et al. 2019](#ref-george2019)) rather than a single
  normal — but keeps it the same at every \\x\\.
- **[`location_scale()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)**
  puts a second forest on the log standard deviation, so the spread is
  an unrestricted function of the predictors.
- **[`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)**
  assumes nothing about the error at all; see [a continuous outcome as
  ordinal](#a-continuous-outcome-as-ordinal).

Measured on 1000 training and 1000 test observations, 50 trees, 500
draws after 500 warmup, three replicates, with the error centered so
that \\E\[Y \mid x\]\\ is the same function in every row. RMSE and
coverage are for that regression function on held-out data; the log
score is the held-out predictive density.

| Errors | [`gaussian()`](https://rdrr.io/r/stats/family.html) | [`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | [`location_scale()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | `ordinal("probit")` |
|----|----|----|----|----|
| normal | 0.135 / -1439 | 0.140 / -1440 | 0.137 / -1442 | **0.129** |
| \\t_3\\ | 0.094 / -1412 | **0.080 / -1253** | 0.101 / -1400 | 0.105 |
| skewed | 0.078 / -1106 | **0.053 / -956** | 0.078 / -1115 | 0.089 |
| bimodal | 0.149 / -1657 | **0.060 / -1240** | 0.136 / -1656 | 0.166 |
| heteroskedastic | 0.168 / -1684 | 0.165 / -1653 | **0.136 / -1541** | 0.147 |
| seconds | **1.6** | 2.1 | 17.4 | 2.7 |

Coverage was between 0.94 and 1.00 everywhere, so it does not separate
them. Four things to take from this.

**On normal errors they all tie.** Within a hair on RMSE and within
three log points, even though
[`gaussian()`](https://rdrr.io/r/stats/family.html) is exactly right
there. This is what makes
[`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
a sensible default rather than a specialist tool: it costs nothing when
the simpler assumption holds.

**[`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is the family for a badly *shaped* error.** Heavy tails are worth 159
log points over [`gaussian()`](https://rdrr.io/r/stats/family.html),
skewness 150, and bimodality 417 with 40% of the RMSE. If you are unsure
what the errors look like, this is the safe choice, and it is why it is
the default.

**[`location_scale()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is the family for a *varying* error, and only that.** It is level with
[`gaussian()`](https://rdrr.io/r/stats/family.html) in the first four
rows and wins the fifth by 112 log points over
[`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md).
Reach for it when a residual plot fans out or when variability is itself
the question. It costs roughly ten times a Gaussian fit, almost all of
it in the second forest, so give that forest fewer trees —
`num_trees = c(50, 10)` was two and a half times faster than `c(50, 50)`
with the same accuracy on both surfaces.

``` r

d$het <- rnorm(n, 2 * d$x1, exp(-1 + 1.5 * d$x2))

fit_ls <- bartisan(het ~ x1 + x2, data = d, family = location_scale(),
                   num_trees = c(mean = 10, log_sd = 5), control = ctrl)

colnames(predict(fit_ls, type = "link"))
#> [1] "mean"   "log_sd"
```

### Giving one forest its own predictors and settings

The two forests of a location-scale model are separate models of
separate things, and every argument that could differ between them may.
The tree counts above are keyed by forest name rather than ordered,
which is the same thing written more legibly. The model formula works
this way too, so the scale forest can have predictors the mean forest
does not, or fewer:

``` r

fit_sub <- bartisan(list(het ~ x1 + x2, log_sd = ~ x2), data = d,
                    family = location_scale(),
                    num_trees = c(mean = 10, log_sd = 5),
                    sparsity = c(mean = TRUE, log_sd = FALSE),
                    control = ctrl)

colMeans(fit_sub$counts$log_sd)
#>   x1   x2 
#> 0.00 7.89
```

The scale forest never splits on `x1`, because its formula does not name
it. The predictor is still in the data, so nothing about
[`predict()`](https://rdrr.io/r/stats/predict.html) or `newdata`
changes; the forest simply never uses it. This is useful when you know
which covariates drive the spread, or when the second forest is spending
capacity on predictors that only matter for the mean.

`?bartisan-families` lists the forests of every multi-forest family in
order, with their names.

**[`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
never wins by much and never loses by much.** It has the lowest error on
normal errors and the second lowest on heteroskedastic ones, which is
what you would expect from a model that writes down no error
distribution. Prefer it when the outcome is bounded, heavily rounded, or
piles up at a floor or ceiling, where a continuous density smears mass
across values the outcome cannot take.

[`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)’s
additive predictor is the conditional mean, as
[`gaussian()`](https://rdrr.io/r/stats/family.html)’s is, so
`type = "link"` and `type = "response"` agree and either can be compared
against a truth or across families.
[`error_density()`](https://ngreifer.github.io/bartisan/reference/error_density.md)
returns the fitted error density with a pointwise interval, which is the
object the family exists to produce:

``` r

d$heavy <- 2 * d$x1 + rt(n, df = 3)

fit_dpm <- bartisan(heavy ~ x1 + x2, data = d, family = dpm(), control = ctrl)

head(error_density(fit_dpm, at = c(-2, 0, 2)))
#>   at   mean  lower  upper
#> 1 -2 0.0720 0.0536 0.0903
#> 2  0 0.3220 0.2828 0.3655
#> 3  2 0.0849 0.0643 0.1028
```

## Positive continuous responses

`Gamma("log")` puts the forest on the log mean and draws the shape,
which acts as the inverse dispersion. The shape is not regressed on the
predictors and cannot be fixed by hand.

**Write the link.**
[`stats::Gamma()`](https://rdrr.io/r/stats/family.html) defaults to
`link = "inverse"`, which is the canonical link for the gamma and the
wrong one here: the additive predictor is unconstrained, and only the
log link’s inverse keeps the mean positive over the whole line. A draw
that wanders non-positive has no gamma density at all, so it is rejected
during fitting and [`predict()`](https://rdrr.io/r/stats/predict.html)
returns `NaN` there. So *bartisan* ignores any other link and fits on
`log`, saying so once; writing `Gamma("log")` silences that.
[`stats::Gamma()`](https://rdrr.io/r/stats/family.html) itself is
untouched, so attaching the package cannot change what
[`glm()`](https://rdrr.io/r/stats/glm.html) does. If you genuinely want
another link for a gamma response,
[`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is the route.

`Gamma("log")` is not the only way to model a positive outcome and often
not the best. Measured on 800 training and 800 test observations, 50
trees, 500 draws after 500 warmup, four replicates, over four shapes of
error around the same mean function. `ordinal("probit")` bins the
outcome onto 25 quantiles and predicts with `type = "mean"`; the
log-Gaussian’s score carries the Jacobian, so the scores are comparable
across the rows that have one.

| Errors | `Gamma("log")` | [`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | [`gaussian()`](https://rdrr.io/r/stats/family.html) on `log(y)` | `ordinal("probit")` |
|----|----|----|----|----|
| gamma, RMSE | 0.401 | 0.856 | **0.388** | 0.414 |
| gamma, log score | **-2020** | -2096 | -2037 | – |
| lognormal, RMSE | 0.234 | 0.400 | **0.225** | 0.232 |
| lognormal, log score | -1526 | -1596 | **-1520** | – |
| heavy tail, RMSE | 2.454 | 2.215 | **2.110** | 2.269 |
| mixed, RMSE | 1.606 | 1.808 | **1.419** | 1.610 |
| mixed, log score | -1935 | **-1912** | -1989 | – |
| seconds | 7.5 | **1.5** | **1.3** | 2.2 |

**[`gaussian()`](https://rdrr.io/r/stats/family.html) or
[`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
on `log(y)` is the first thing to try.** It had the best RMSE for the
conditional mean in all four settings and is five times faster. It is a
different model rather than a reparameterization — the forest is on the
mean of the log rather than the log of the mean, and predictions come
back on the log scale — so it is the right choice when a multiplicative
error is what you believe in.

**`Gamma("log")` earns its place when the mean on the original scale is
the target.** It has the best predictive density where the gamma is
correctly specified, which no other family here manages. Choose it
because a constant coefficient of variation is what you believe, not
because the outcome is positive.

**`ordinal("probit")` on 25 bins tracks `Gamma("log")` closely** —
within about 3% of it in every row, and 8% better on heavy tails — at
under a third of its time and with no assumption about the error’s
shape. It cannot give a predictive density on the original scale, which
is the price, and it is behind the log-Gaussian on RMSE throughout.

**[`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
on the raw scale is the one clear mistake here**, at roughly twice the
RMSE of everything else under gamma or lognormal errors. A right-skewed
error on the raw scale is exactly what its flexible error distribution
absorbs signal from. Use it on `log(y)` instead.

The heavy-tail row has no log score because that column is not stable
enough to report: the score is dominated by the single worst test point,
and across eight replicates `Gamma("log")` ranged from -1745 to -557785
while
[`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
ranged from -1661 to -15602. The medians are close; the tail risk is
not, and the gamma carries much more of it. If the predictive density
far out in the tail matters, prefer the mixture.

## Binary responses

[`binomial()`](https://rdrr.io/r/stats/family.html) accepts a 0/1
numeric response, a logical, a two-level factor, or a two-column matrix
of successes and failures.

For prediction the link hardly matters — logit and probit differ by a
scale factor of about 1.6 and give fitted probabilities that are hard to
tell apart — so choose on grounds other than fit.

**`logit`** is the default. Its predictor is a log odds, which is the
scale most readers of a binary model expect, and it is the link under
which a contrast has an odds-ratio reading.

**`probit`** reads the model as a normal latent variable crossing a
threshold, which is the right link when that latent variable is the
object of interest, and it is the one `predict(type = "stdlv")` is most
natural for. It is also the fastest of the three.

**`cloglog`** is the one substantively different choice, because it is
not symmetric: swapping the labels of success and failure gives a
different model, where for logit and probit it gives the same model with
the predictor negated. Use it when that asymmetry is the point — when
the outcome is “at least one event occurred” and the underlying count is
Poisson, where the complementary log-log model is exactly right, or in
discrete time, where the cumulative version *is* the discrete
proportional hazards model and the predictor is a log hazard ratio.

With very few events the whole fit mixes slowly whatever the link,
because a handful of events carry little information: at 2000
observations the effective sample size of the level of the predictor
fell from about 860 at half positives to 8 at 0.2%. Lengthen the chain
and read `fit$rhat`.

## Counts

[`poisson()`](https://rdrr.io/r/stats/family.html) asserts that the
variance equals the mean.
[`negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
adds a dispersion parameter, drawn by default and reported as `theta`.
[`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
and
[`zi_negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
mix a point mass at zero with a count component and give each part its
own forest — the first predictor is the log mean of the count, the
second the log odds of being a structural zero — so the excess-zero
mechanism can depend on the predictors.

Two questions, in this order.

**Are there more zeros than the count component can produce?** Not “are
there many zeros”: a Poisson with a small mean produces plenty, and a
forest can drive the mean low exactly where the zeros are. Fit the plain
family and the zero-inflated one and compare with
[`loo()`](https://mc-stan.org/loo/reference/loo.html). Reach for zero
inflation when you want to *model* the zero mechanism, or when the two
processes have different predictors — not merely to improve a fit.

**Is the count component overdispersed once the zeros are accounted
for?** A spike at zero inflates the sample variance and looks like
dispersion.
[`zi_negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
separates the two at the cost of a parameter; if `theta` comes back
large with a tight posterior, the negative binomial is not buying
anything and
[`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is the better-conditioned fit.

**[`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is a fourth option worth knowing about for counts.** A count is ordered
and takes few distinct values, which is exactly what the cutpoint
structure is for: it makes no assumption about the count distribution,
needs no dispersion parameter, and handles excess zeros without a
mixture, since the zero category simply gets whatever probability the
cutpoints give it. Predict with `type = "mean"`. The limits are that it
cannot predict a count larger than the largest one observed, and it has
no rate interpretation — no log link, so no incidence-rate ratio. Bin
first if the counts range over more than a few dozen values, as below.

## Ordered categories

[`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
expects an ordered factor, though a numeric response is accepted and its
sorted unique values are taken as the categories. It uses the
cumulative-link parameterization of
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html), in which \\P(Y
\le k) = F(c_k - \eta)\\, so larger values of the additive predictor
shift mass towards higher categories.

Only the differences \\c_k - \eta_i\\ are identified, so one location
has to be pinned. With three or more categories the draws are reported
in the chart where **the additive predictor has mean zero over the
fitted sample and every cutpoint is free**, which is what `polr()`
reports when its predictors are centered and which makes the cutpoints
readable as category boundaries. With exactly two categories the single
boundary is folded into the intercept instead, so a two-category
response is exactly binary regression on the same scale.

The three links are read as for the binomial, plus one consideration
specific to the ordinal case: the link decides what is held constant
across categories. `logit` gives the proportional odds model, in which a
contrast is a log odds ratio that is the same at every cutpoint.
`probit` gives the ordered probit, and is the link to choose when the
latent variable is the object of interest — a trait, a utility, an
underlying measurement recorded in bins. `cloglog` gives the
proportional *hazards* model on the categories, so choose it when the
categories are ordered durations or stages. Accuracy is much the same
across the three; `probit` is the fastest by a wide margin.

### A continuous outcome as ordinal

[`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
on a numeric response is not a workaround. Every distinct value of \\Y\\
becomes a category, the cutpoints absorb the marginal distribution of
\\Y\\, and the forest is left to explain only the *ordering*. Nothing is
assumed about the error, and the model for \\P(Y \le y \mid x)\\ is
invariant to any monotone transformation of \\Y\\: fitting on \\Y\\,
\\\log Y\\ or \\\sqrt Y\\ gives the same model. This is the
semiparametric ordinal regression `rms::orm()` fits for continuous
outcomes, with a forest in place of the linear predictor. Predict with
`type = "mean"`, which reads the category labels as numbers and returns
\\\sum_k y\_{(k)} P(Y = y\_{(k)} \mid x)\\.

**Bin the outcome first.** One cutpoint per distinct value means \\n\\
cutpoints: at 1000 observations an unbinned fit took 73 seconds against
2.7 for [`gaussian()`](https://rdrr.io/r/stats/family.html). Collapsing
onto a grid of quantiles is faster *and* more accurate, because a
cutpoint vector with a thousand weakly-identified entries is worse
conditioned than one with twenty-five. Somewhere between 10 and 50
quantile bins, and the choice inside that range hardly matters; 25 bins
was sixteen times faster than no binning and slightly more accurate.
Keep each bin’s mean as the label so that `type = "mean"` still reports
on the outcome’s scale:

``` r

d$cont <- 2 * sin(pi * d$x1) + d$x2 + rt(n, df = 3)

bins <- quantile(d$cont, seq(0, 1, length.out = 26))
at <- cut(d$cont, unique(bins), include.lowest = TRUE, labels = FALSE)
d$binned <- ave(d$cont, at)

fit_oc <- bartisan(binned ~ x1 + x2, data = d, family = ordinal("probit"),
                   control = ctrl)

head(predict(fit_oc, type = "mean"))
#> [1] 2.237 2.333 1.243 2.238 2.704 0.203
```

Two limits. `type = "mean"` is a convex combination of observed outcome
values, so it can never predict outside the range of the training
outcome — a feature when the outcome has a hard floor or ceiling, a
liability when you need extrapolation. And the invariance is a property
of the model for \\P(Y \le y \mid x)\\, not of the mean read off it:
`type = "mean"` after fitting on \\\log Y\\ is not the log of
`type = "mean"` after fitting on \\Y\\.

## Unordered categories

[`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
expects an unordered factor. Start with `logit`: it is faster, its
predictor is a log odds, and it is the parameterization most readers
expect. By default it fits one forest per category and leaves the model
unidentified ([Murray 2021](#ref-murray2021)), which keeps the prior
symmetric in the categories; every identified quantity is still
recovered from the draws. Passing `reference` pins that category and
fits one fewer forest, giving log odds against it the way
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) does —
the right choice when a particular contrast is the quantity of interest.

Choose `probit` when the alternatives are plausibly *substitutes*: it
allows the latent utilities to be correlated, which a multinomial logit
cannot express at all, so it is the link for when the
independence-of-irrelevant-alternatives assumption is what you doubt.
Expect to pay in speed and in a simulated likelihood — the probabilities
are Gaussian orthant probabilities with no closed form, so they and the
reported log likelihood are computed by simulation.

**The latent correlations are weakly identified. Do not report them as
estimates.** They enter the likelihood only through orthant
probabilities of a distribution whose location is a sum of trees, and a
flexible mean absorbs much of the dependence they are meant to capture.
At 900 observations the posterior tracked a swept true correlation only
loosely and not monotonically, with 95% intervals from 0.38 to 1.07 wide
on a parameter confined to \\(-1, 1)\\; repeating the zero case over
eight draws of the data gave posterior means from \\-0.60\\ to \\0.33\\.
By 3000 observations it behaves. What the probit link is still good for
is the *fit*: allowing correlated errors changes the category
probabilities whether or not \\\Sigma\\ is pinned down. Read the
probabilities, not the covariance.

## Bounded responses

[`Beta()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is beta regression for a response strictly inside the unit interval: a
forest on the link of the mean and a precision \\\phi\\ drawn alongside
it, so the response is \\\mathrm{Beta}(\mu\phi, (1-\mu)\phi)\\ with
variance \\\mu(1-\mu)/(1+\phi)\\. The same three links as
[`binomial()`](https://rdrr.io/r/stats/family.html), read the same way.

``` r

d$rate <- rbeta(n, plogis(1.5 * sin(pi * d$x1)) * 12,
                12 - plogis(1.5 * sin(pi * d$x1)) * 12)

fit_beta <- bartisan(rate ~ x1 + x2, data = d, family = Beta(), control = ctrl)
mean(fit_beta$aux[, "phi"])
#> [1] 11
```

[`ordbeta()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is the ordered beta regression of Kubinec ([2023](#ref-kubinec2023)),
for a response on the *closed* unit interval with point masses at zero
and one — a percentage of a budget, a slider scale. One predictor drives
both the probability of landing on an endpoint, through a pair of
cutpoints as in an ordinal model, and the mean of the beta density in
between.

**Choose between them on whether the response can reach a boundary, not
on whether it happens to in your sample.** A response *at* zero or one
has no beta density, so
[`Beta()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
errors rather than nudging it inward. In the other direction,
[`ordbeta()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
fitted to a response with no boundary observations leaves its two
cutpoints with nothing to identify them; the fit is not much worse for
it, but the cutpoints are not interpretable and drift to large
magnitudes.

The other alternative to weigh is a plain
[`binomial()`](https://rdrr.io/r/stats/family.html) fit on the
proportion with `weights`, which is right when the response really is a
count of successes out of a known denominator and wrong when it is a
continuous proportion that happens to reach its bounds.

## Right-censored survival times

Five routes, all taking a
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) object
or a two-column matrix of times and 0/1 event indicators. They differ in
what the predictor means and in what the model leaves free:

| Family | Model | The predictor is | Left free |
|----|----|----|----|
| [`weibull_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | \\\log T = \eta + \sigma\epsilon\\, \\\epsilon\\ Gumbel | a log **time** ratio | \\\sigma\\ |
| [`loglogistic_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | the same with \\\epsilon\\ logistic | a log **time** ratio | \\\sigma\\ |
| [`lognormal_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | the same with \\\epsilon\\ normal | a log **time** ratio | \\\sigma\\ |
| [`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | \\\lambda(t \mid x) = \lambda_0(t)e^{r(x)}\\ | a log **hazard** ratio | the baseline \\\lambda_0\\ |
| [`dpm_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md) | \\\log T = m(x) + W\\, \\W\\ a mixture | \\E\[\log T \mid x\]\\ | the whole error density |

The three accelerated failure time families fix the shape of the error
and so of the hazard, and are quick.
[`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
frees the baseline hazard instead, so it fits a hazard of any shape, at
the price of asserting proportionality and of reading the predictor on
the hazard scale.
[`dpm_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
frees the error density ([Henderson et al. 2020](#ref-henderson2020)),
which costs little when a single normal was right and gains a great deal
when it was not.

[`dpm_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is the default for a `Surv` response, on the evidence in the survival
vignette and because it is cheaper to fit than most of the alternatives.
[`weibull_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is the one family that is both an accelerated failure time and a
proportional hazards model, so name it when a hazard-ratio reading is
what you want.

``` r

d$time <- rexp(n, exp(-(1 + d$x1)))
d$event <- rbinom(n, 1, 0.7)

fit_aft <- bartisan(survival::Surv(time, event) ~ x1 + x2, data = d,
                    family = weibull_aft(), control = ctrl)
colMeans(fit_aft$aux)
#> sigma 
#>  1.01
```

The survival function comes from `type = "survival"`, which takes the
times to report it at and works for every family in the table:

``` r

fit_ph <- bartisan(survival::Surv(time, event) ~ x1 + x2, data = d,
                   family = ph(), control = ctrl)

head(predict(fit_ph, type = "survival", times = c(1, 2, 5)), 3)
#>          1     2     5
#> [1,] 0.880 0.770 0.513
#> [2,] 0.861 0.736 0.458
#> [3,] 0.828 0.680 0.373
```

It is also the estimand you usually want. The question is rarely about
the predictor but about survival at a horizon — the difference in
one-year survival between two groups, say — and that is a contrast of
`type = "survival"` at one time:

``` r

avg_comparisons(fit_ph, variables = "trt", type = "survival", times = 1)
```

**A quick way to choose.** Take
[`dpm_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
as the default when you have no view about the shape of anything —
measured against five alternatives over six data-generating truths, it
was best or tied-best on four and never worse than third, and it matched
the correctly specified family on the two truths where one existed. Take
[`lognormal_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
or
[`loglogistic_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
when the fit has to be quick,
[`weibull_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
when a hazard-ratio reading is wanted from a parametric fit, and
[`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
when the shape of the baseline hazard is the point. If the covariate
effect may **move with time**, so that survival curves cross, none of
the five is right and the discrete-time route is.

That is the summary; the trade-offs behind it are measured in
[`vignette("survival", package = "bartisan")`](https://ngreifer.github.io/bartisan/articles/survival.md),
which covers what each estimand is, which hazard shapes each family can
and cannot represent, how they behave under censoring and under
misspecification, the discrete-time route for non-proportional hazards,
and one trap in comparing their log scores.

## Other links, and your own likelihood

Links beyond those in the table are accepted for
[`gaussian()`](https://rdrr.io/r/stats/family.html),
[`binomial()`](https://rdrr.io/r/stats/family.html),
[`poisson()`](https://rdrr.io/r/stats/family.html) and
[`negbin()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
and applied from R by composing your inverse link with the family’s own.
So `binomial("cauchit")` works, as does any link object of the kind
[`stats::make.link()`](https://rdrr.io/r/stats/make.link.html) returns,
including one written by hand. Two cautions: the fit is slower, because
each leaf costs a call into R; and the additive predictor is
unconstrained, so a link whose inverse has a restricted range
(`poisson("identity")`) gives non-finite densities for some predictors,
which are rejected rather than breaking the chain but are wasted work.
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
reports that when the fit starts. The families with more than one
predictor, or whose link enters somewhere other than a single mean, take
only their listed links.

[`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
takes the log density itself and fits the model that goes with it. The
sampler needs only the first two derivatives with respect to each
additive predictor, and central differences of your function produce
both. Terms free of `eta` may be dropped: they cancel from every
acceptance ratio.

``` r

pois_by_hand <- custom_family(function(y, eta) y * eta[, 1] - exp(eta[, 1]),
                              start = log(mean(d$count)),
                              name = "hand-rolled Poisson")

fit_custom <- bartisan(count ~ x1 + x2, data = d, family = pois_by_hand,
                       control = ctrl)

cor(predict(fit_custom, type = "link"),
    predict(bartisan(count ~ x1 + x2, data = d, family = poisson(),
                    control = ctrl), type = "link"))
#> [1] 0.997
```

The function is called once per leaf per Fisher-scoring step with the
observations reaching that leaf, so it must be vectorized over `y` and
the rows of `eta` and must return exactly one value per row. Ask for
several additive predictors with `num_predictors`; supply `derivatives`
when they are easy to write down, which cuts three calls to one and
removes the differencing error. See
[`?custom_family`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
for both.

Nuisance parameters are drawn too, if you name them. `logdens` then
takes a third argument holding their current values, and the draws come
back in `fit$aux` under the names you gave, covered by
[`summary()`](https://rdrr.io/r/base/summary.html) and `fit$rhat` like
any other family’s:

``` r

by_hand <- custom_family(
  logdens   = function(y, eta, aux) dnorm(y, eta[, 1], exp(aux[1]), log = TRUE),
  aux_names = "log_sigma",
  aux_start = 0,
  start     = mean(d$cont))

fit_aux <- bartisan(cont ~ x1 + x2, data = d, family = by_hand, control = ctrl)
exp(mean(fit_aux$aux[, "log_sigma"]))
#> [1] 2.11
```

There is no prior argument and no bounds argument. A parameter with a
restricted range is handled the way it would be for a real predictor —
by writing the transform into `logdens`, which is what the
[`exp()`](https://rdrr.io/r/base/Log.html) above is doing. `aux_start`
need only be the right order of magnitude; the sampler walks to the
posterior from a poor start.

What it cannot do: take a non-numeric response, so a factor has to be
coded first; or report a fitted mean, since the package cannot know what
the mean of your density is — `predict(type = "response")` returns the
additive predictors instead, though `type = "density"` works.

If a family you want is missing and
[`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is not enough, that is worth raising as an issue rather than working
around.

## References

George, Edward, Purushottam Laud, Brent Logan, Robert McCulloch, and
Rodney Sparapani. 2019. “Fully Nonparametric Bayesian Additive
Regression Trees.” In *Topics in Identification, Limited Dependent
Variables, Partial Observability, Experimentation, and Flexible
Modeling: Part b*, vol. 40B. Advances in Econometrics. Emerald
Publishing Limited. <https://doi.org/10.1108/S0731-90532019000040B006>.

Henderson, Nicholas C., Thomas A. Louis, Gary L. Rosner, and Ravi
Varadhan. 2020. “Individualized Treatment Effects with Censored Data via
Fully Nonparametric Bayesian Accelerated Failure Time Models.”
*Biostatistics* 21 (1): 50–68.
<https://doi.org/10.1093/biostatistics/kxy028>.

Kubinec, Robert. 2023. “Ordered Beta Regression: A Parsimonious,
Well-Fitting Model for Continuous Data with Lower and Upper Bounds.”
*Political Analysis* 31 (4): 519–36.
<https://doi.org/10.1017/pan.2022.20>.

Linero, Antonio R. 2025. “Generalized Bayesian Additive Regression Trees
Models: Beyond Conditional Conjugacy.” *Journal of the American
Statistical Association* 120 (549): 356–69.
<https://doi.org/10.1080/01621459.2024.2337156>.

Murray, Jared S. 2021. “Log-Linear Bayesian Additive Regression Trees
for Multinomial Logistic and Count Regression Models.” *Journal of the
American Statistical Association* 116 (534): 756–69.
<https://doi.org/10.1080/01621459.2020.1813587>.
