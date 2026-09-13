# Has it converged, and does it fit?

## Introduction

Two questions get asked together and are worth separating. Has the
sampler converged, meaning has it explored the posterior properly? And
does the model fit, meaning does it describe the data?

The first is about the algorithm and the second is about the model; a
fit can converge beautifully on a badly chosen family, and a well chosen
family can be fitted by a chain that has not run long enough.

In this guide, we will take the two questions in that order. First we’ll
fit a model with four chains, hand it to
[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md),
and then unpack what that summary reports and why, including the one row
(the additive predictor `eta`) whose conventional threshold does not
apply to a forest and the one quantity that is deliberately left out of
the table. Next we’ll ask whether the model describes the data, using
posterior predictive checks, a calibration plot for the binary case in
which those checks are weak, residuals, and a Bayesian \\R^2\\. A short
checklist at the end collects what we recommend running before reporting
anything from a fit.

``` r

library(bartisan)

# For parallelization; optional
if (rlang::is_installed("future")) {
  future::plan(future::multisession)
}

data("rhc")

set.seed(2026)

fit <- bartisan(death ~ rhc + age + sex + race + edu + aps + meanbp + resp +
                  hema + pafi + paco2 + crea + surv2m + card,
                data = rhc, family = binomial(), chains = 4)
```

## Convergence

### Every Statistic in One Call (`diagnose()`)

[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
computes every convergence and mixing statistic in this section, says
which of them fall short, and says what to do about each. Everything it
reports comes from the stored draws, so it needs no other package.

``` r

diagnose(fit)
#> Convergence and mixing
#> 
#>                             quantity rhat rhat_late ess_bulk ess_tail
#>                               loglik 1.12     1.036       23      234
#>                           splits.eta 1.01     1.039      264      731
#>  eta.eta (average over observations) 1.00     0.999     1610     2278
#>   eta.eta (worst 5% of observations) 1.04     1.073       80      305
#> 
#> ✔ 4 chains, 3200 draws kept in total
#> ✖ above 1.01 for loglik
#> ✖ that R-hat rests on 23 effective draws, where 4 chains average 1.171 even
#>   when they agree
#> ℹ not the fix: R-hat stays high on the second half of the draws alone as well,
#>   so a longer warmup is not what is missing
#> ✖ the chains disagree about how many splitting rules the forest has (R-hat
#>   1.01)
#> ✖ 23 for loglik, below 400
#> ✖ 234 for loglik, below 400
#> ℹ the chains disagree about individual observations and agree about their
#>   average (R-hat 1.00, 1610 effective draws)
#> ℹ loglik carries 0.7 effective draws per hundred kept
#> 
#> What to do
#> 
#> Raise `num_draws`, which was `800`. R-hat is above the threshold for a quantity
#> that carries too few effective draws for the threshold to mean anything: with
#> this many chains it would sit about where it does even if the chains agreed
#> exactly, as the check above reports. Effective sample size is what makes it
#> readable, and that grows with the total number of draws; using fewer chains
#> lowers the bar as well, since R-hat's null rises with the number of chains
#> being compared.
#> If that does not settle it, reduce `num_trees`. A smaller forest has fewer ways
#> to represent the same fit, so the sampler has less room to move between them.
#> Then check the family. A likelihood that fits the data badly can give a
#> posterior with no single place to be; `bayesplot::pp_check()` is the
#> diagnostic.
#> Note that the chains disagree about the fitted values of individual
#> observations and not about their average, which is the usual shape of this in a
#> forest. An estimand averaged over observations therefore carries far more
#> effective draws than the table's worst row does, and R-hat for that estimand is
#> worth computing rather than inferring; `posterior::as_draws()` hands the draws
#> over for it.
```

The rest of this section is what it is reporting and why, which is worth
reading once; when the summary itself is enough,
[`vignette("bartisan")`](https://ngreifer.github.io/bartisan/articles/bartisan.md)
is the shorter tour of the whole workflow.

### Running More Than One Chain (`chains`)

`chains = 4` above is doing the work. The default is one chain, which
produces estimates but no way to check most of what matters, because
R-hat compares chains to each other; with one chain it is computed by
splitting that chain, which catches drift but cannot catch two chains
settling in different places. Running four costs four times as much
sampling, which for most fits is a few seconds, and less than that under
a `future` plan.

More chains is not, however, a way to improve `rhat`, and it is worth
knowing which way it cuts before reaching for it. Effective sample size
depends on the total number of draws and not on how they are divided, so
twice as many chains and twice as long a chain buy the same amount of
it. `rhat` is not symmetric that way: because it compares chains against
each other, its value under a sampler with nothing wrong is about
`1 + chains / ess`, so at 400 effective draws four chains sit at 1.010
and sixteen chains at 1.040. Adding chains raises the bar it has to
clear. When a fit fails on `rhat`, longer chains are the fix; when it
fails only on effective sample size, either will do.

[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
is where the diagnostics are computed.

``` r

diagnose(fit)
#> Convergence and mixing
#> 
#>                             quantity rhat rhat_late ess_bulk ess_tail
#>                               loglik 1.12     1.036       23      234
#>                           splits.eta 1.01     1.039      264      731
#>  eta.eta (average over observations) 1.00     0.999     1610     2278
#>   eta.eta (worst 5% of observations) 1.04     1.073       80      305
#> 
#> ✔ 4 chains, 3200 draws kept in total
#> ✖ above 1.01 for loglik
#> ✖ that R-hat rests on 23 effective draws, where 4 chains average 1.171 even
#>   when they agree
#> ℹ not the fix: R-hat stays high on the second half of the draws alone as well,
#>   so a longer warmup is not what is missing
#> ✖ the chains disagree about how many splitting rules the forest has (R-hat
#>   1.01)
#> ✖ 23 for loglik, below 400
#> ✖ 234 for loglik, below 400
#> ℹ the chains disagree about individual observations and agree about their
#>   average (R-hat 1.00, 1610 effective draws)
#> ℹ loglik carries 0.7 effective draws per hundred kept
#> 
#> What to do
#> 
#> Raise `num_draws`, which was `800`. R-hat is above the threshold for a quantity
#> that carries too few effective draws for the threshold to mean anything: with
#> this many chains it would sit about where it does even if the chains agreed
#> exactly, as the check above reports. Effective sample size is what makes it
#> readable, and that grows with the total number of draws; using fewer chains
#> lowers the bar as well, since R-hat's null rises with the number of chains
#> being compared.
#> If that does not settle it, reduce `num_trees`. A smaller forest has fewer ways
#> to represent the same fit, so the sampler has less room to move between them.
#> Then check the family. A likelihood that fits the data badly can give a
#> posterior with no single place to be; `bayesplot::pp_check()` is the
#> diagnostic.
#> Note that the chains disagree about the fitted values of individual
#> observations and not about their average, which is the usual shape of this in a
#> forest. An estimand averaged over observations therefore carries far more
#> effective draws than the table's worst row does, and R-hat for that estimand is
#> worth computing rather than inferring; `posterior::as_draws()` hands the draws
#> over for it.
```

The table is in `diagnose(fit)$table` when it is wanted as a data frame,
and the checks and the advice are in `$checks` and `$advice`.

`rhat` compares the variance between chains to the variance within them
([Vehtari et al. 2021](#ref-vehtari2021)); if the chains have found the
same posterior the two agree and the ratio is near 1. `ess_bulk` and
`ess_tail` are effective sample sizes, meaning how many independent
draws the correlated ones are worth, in the middle of the distribution
(i.e., the bulk) and in the tails respectively.

The conventional thresholds are `rhat` below 1.01 and effective sample
sizes above about 400 for a quantity we intend to report precisely.
Those come from the general MCMC literature and are a reasonable default
here, with one exception described below.
[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
applies them, and takes both as arguments so they can be moved.

Two of the columns are worth naming separately. `rhat_late` is the same
statistic computed on the second half of the retained draws alone, which
is what tells a warmup that ended too early (i.e., too small a
`num_burn`) from chains that have each settled somewhere different:
throwing away the early draws is exactly what more `num_burn` would have
done, so if that fixes `rhat`, warmup was the problem. The `splits.*`
row is the total number of splitting rules in the forest at each draw,
the one quantity here that is about the trees rather than about the
fitted values, and which no general-purpose MCMC diagnostic would think
to look at.

### The Rows of the Table

`loglik` is the log likelihood of the whole dataset at each draw. It is
a useful scalar summary of the fit and mixes reasonably.

`aux.*` are the nuisance parameters of the family (i.e., the parameters
that are not part of the additive predictor), and are usually the best
behaved rows in the table. A binomial likelihood has none, so they do
not appear here; a Gaussian fit would show `aux.sigma`, and a survival
fit its scale or baseline hazard.

`eta.*` is the additive predictor, and it gets two rows. One summarizes
the worst 5% of observations rather than the single worst, because the
worst of a thousand values is extreme even when every chain has
converged; the other is the average over observations. Reading them
against each other is the point of having both, because they routinely
differ by orders of magnitude: a forest settles the level of the fitted
function within a sweep or two and takes far longer to settle which
observation gets which share of it. Which row binds depends on what we
mean to report. An average or a contrast of averages, which is what
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
and `avg_comparisons()` give, is governed by the average row; a
prediction for one observation is governed by the other.

One quantity is deliberately not in the table. The leaf prior scale, in
`fit$sigma_mu`, mixes badly in every BART implementation: on one dataset
its counterparts in *dbarts* and in *stochtree* come out at `rhat` 1.12
and 1.16 with effective sample sizes of 22 and 17 out of 3000 draws, and
*stochtree* draws it exactly from its full conditional, so no sampler
can do better. It is a hyperparameter whose disagreement between chains
does not reach the fitted function, which on those same fits has `rhat`
1.00 and thousands of effective draws, so reporting it beside `eta` only
produces an alarming row with nothing to act on. It is out of this table
and only out of this table: the draws are in `fit$sigma_mu` and reach
[`as_draws()`](https://mc-stan.org/posterior/reference/draws.html), so
it can still be diagnosed by anyone who wants to.

### The Exception for `eta`

The `eta` row is a high percentile over every observation, so it is
conservative by construction, and forests mix slowly on their fitted
values. Two chains can visit quite different collections of trees that
imply nearly identical predictions, which inflates a between-chain
statistic without meaning the two chains disagree about anything we
would report, so whether two chains found the same trees is not a
question worth asking of a forest.

On harder problems this shows up clearly. Below we fit the Friedman
function, first with a chain far too short and then at the defaults.

``` r

set.seed(3)
fr <- as.data.frame(matrix(runif(400 * 10), 400, 10))
names(fr) <- paste0("x", 1:10)
fr$y <- 10 * sin(pi * fr$x1 * fr$x2) + 20 * (fr$x3 - 0.5)^2 +
  10 * fr$x4 + 5 * fr$x5 + rnorm(400)

too_short <- bartisan(y ~ ., fr, family = gaussian(), chains = 4,
                      control = bartisan_control(num_trees = 20, num_burn = 50,
                                                 num_draws = 50))

diagnose(too_short)$table
#>                              quantity rhat rhat_late ess_bulk ess_tail ess_frac
#> 1                              loglik 1.72      1.79     6.85     12.3   0.0342
#> 2                           aux.sigma 1.26      1.10    12.33     29.6   0.0617
#> 3                          splits.eta 1.57      1.96     7.73     17.8   0.0387
#> 4 eta.eta (average over observations) 1.02      1.02   203.44    191.7   1.0172
#> 5  eta.eta (worst 5% of observations) 1.86      1.97     6.46     14.6   0.0323
#>   rhat_bad late_bad
#> 1    1.000        1
#> 2    1.000        1
#> 3    1.000        1
#> 4    1.000        1
#> 5    0.998        1
```

Everything is bad. Every `rhat` is far above the 1.01 threshold, and the
effective sample sizes are all under 30 out of the 200 draws kept (4
chains x 50 kept draws per chain). `rhat_late` is no better than `rhat`
on any row, which is the informative part: discarding the early draws
does not fix it, so this is not a warmup that ended too soon but chains
that have each settled somewhere different. Nothing from this fit should
be used.

``` r

long_enough <- bartisan(y ~ ., fr, family = gaussian(), chains = 4)

diagnose(long_enough)$table
#>                              quantity rhat rhat_late ess_bulk ess_tail ess_frac
#> 1                              loglik 1.12      1.14     26.9    127.8  0.00839
#> 2                           aux.sigma 1.03      1.04    130.5    550.1  0.04077
#> 3                          splits.eta 1.08      1.11     37.3    164.5  0.01166
#> 4 eta.eta (average over observations) 1.00      1.00   3322.6   3158.2  1.03830
#> 5  eta.eta (worst 5% of observations) 1.14      1.23     20.5     70.5  0.00641
#>   rhat_bad late_bad
#> 1    1.000    1.000
#> 2    1.000    1.000
#> 3    1.000    1.000
#> 4    0.000    0.000
#> 5    0.993    0.995
```

Everything has improved by a large factor, and the residual standard
deviation is nearly at the threshold. But the log likelihood and `eta`
are both still above 1.01, and `eta`’s effective sample size is still in
the tens. That residual is the phenomenon described above rather than a
chain that merely needs to be longer, and adding draws does not reliably
remove it.

When the `eta` row is elevated, the distribution behind it is more
informative than the maximum:

``` r

per_observation_rhat <- function(fit) {
  draws <- fit$eta$eta
  per <- nrow(draws) / fit$chains
  index <- matrix(seq_len(per * fit$chains), per, fit$chains)
  apply(draws, 2, function(column) {
    posterior::rhat(matrix(column[index], per, fit$chains))
  })
}

round(quantile(per_observation_rhat(fit), c(0.5, 0.9, 0.99, 1)), 3)
#>  50%  90%  99% 100% 
#> 1.02 1.04 1.06 1.07
```

The median is a little above 1 and so is most of the distribution, so
the maximum in the table is not one badly behaved observation: the
fitted function as a whole mixes slowly here. That is the phenomenon
described above rather than a broken chain, and the two are told apart
by what happens to the quantities we report. If the log likelihood has
converged, along with any nuisance parameters the family has, and an
`avg_comparisons()` estimate is stable across chains and across a longer
run, the fit is usable.

If instead the log likelihood were badly elevated too, or the estimate
moved when the chain was lengthened, that would be a chain that has not
converged, and the answer is more draws or a simpler model.

### Looking at the Chains (`as_draws()`)

[`as_draws()`](https://mc-stan.org/posterior/reference/draws.html) hands
the fit to *posterior* and *bayesplot*. It carries the scalar parameters
and a representative set of `eta` columns.

``` r

library(posterior)

as_draws(fit) |>
  summarise_draws()
#> # A tibble: 12 × 10
#>    variable         mean   median     sd    mad       q5      q95  rhat ess_bulk
#>    <chr>           <dbl>    <dbl>  <dbl>  <dbl>    <dbl>    <dbl> <dbl>    <dbl>
#>  1 loglik       -8.32e+2 -8.32e+2 6.26   6.11   -841.    -821.     1.13    23.5 
#>  2 sigma_mu.eta  2.50e-1  2.36e-1 0.0650 0.0577    0.165    0.372  1.33     9.85
#>  3 eta[119]     -1.52e+0 -1.50e+0 0.363  0.348    -2.13    -0.945  1.04   102.  
#>  4 eta[466]     -4.19e-1 -4.06e-1 0.336  0.322    -0.976    0.127  1.02   294.  
#>  5 eta[541]     -2.67e-5  5.60e-3 0.289  0.282    -0.482    0.462  1.01   374.  
#>  6 eta[673]      2.97e-1  2.92e-1 0.311  0.301    -0.219    0.818  1.01   281.  
#>  7 eta[763]      5.78e-1  5.64e-1 0.296  0.289     0.123    1.08   1.01   294.  
#>  8 eta[1261]     8.69e-1  8.63e-1 0.336  0.338     0.319    1.42   1.02   337.  
#>  9 eta[1228]     1.28e+0  1.28e+0 0.341  0.333     0.717    1.83   1.01   269.  
#> 10 eta[137]      1.66e+0  1.65e+0 0.388  0.391     1.04     2.31   1.02   157.  
#> 11 eta[1374]     2.01e+0  2.01e+0 0.485  0.473     1.24     2.80   1.01   288.  
#> 12 eta[1135]     3.16e+0  3.12e+0 0.512  0.491     2.39     4.08   1.06    48.8 
#> # ℹ 1 more variable: ess_tail <dbl>
```

``` r

library(bayesplot)

as_draws(fit, eta = 1) |>
  mcmc_trace(pars = c("loglik", "eta[1]"))
```

![](diagnostics_files/figure-html/trace-1.png)

What we want to see is four chains overlapping, wandering around the
same level, with no drift and no long excursions. `eta = 1` selects the
first observation; `eta = TRUE` gives a spread of ten and `eta = FALSE`
gives none.

### Remedies for a Chain That Has Not Converged

Three things are worth trying, in the order in which they usually help.

**More draws.** Increasing `num_burn` and `num_draws` is the first thing
to try, and it fixes most cases.

**A smaller forest.** Reducing `num_trees` leaves fewer ways to
represent the same function, so the sampler mixes faster.

**A different family.** A likelihood that fits the data badly can
produce a posterior that is hard to explore, so the family is worth
checking when more draws and a smaller forest have not helped;
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
covers the alternatives.

One thing to rule out first. `rhat` above 1.01 on a quantity carrying
only a handful of effective draws is not yet evidence that the chains
disagree, since that is roughly where `rhat` sits when nothing is wrong,
and
[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
says so rather than reporting a disagreement it cannot support. In that
case the effective sample size is what to fix, and more draws is the
whole of the remedy.

## Fit

Convergence says the sampler did its job; it says nothing at all about
whether the model is right.

### Posterior Predictive Checks (`pp_check()`)

A posterior predictive check simulates outcomes from the fitted model
and compares their distribution to the observed one.

``` r

pp_check(fit)
```

![](diagnostics_files/figure-html/ppc-1.png)

For a continuous outcome this is the workhorse check, and systematic
differences are what to look for: replicates that are too narrow, that
miss a second mode, or that put mass where the outcome cannot go.
Simulating negative values for an outcome that cannot be negative says
the family is wrong, and
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
covers the alternatives.

For a binary outcome it is a weak check, which is worth knowing before
reading too much into it. There are only two values the replicates can
take, so they will match the observed proportion unless the model has
gone badly wrong; passing this check tells us almost nothing. The
calibration check below is what to read instead.

#### Checks Worth Reaching For

`type` names any of *bayesplot*’s checks without its `ppc_` prefix, and
four of them earn their place for a forest.

The default compares whole distributions, which is a coarse question. A
test statistic is a sharper one: `type = "stat"` compares any summary of
the replicates with the same summary of the outcome, so the spread, a
tail quantile, or the share of zeros can each be asked about on its own.
The share of zeros is the one that finds an unmodeled point mass, which
is what
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
recommends
[`zi_poisson()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
and
[`tweedie()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
for.

``` r

pp_check(fit, type = "stat", stat = "sd")
```

![](diagnostics_files/figure-html/ppcstat-1.png)

**Calibration of the whole predictive distribution** is what
`type = "loo_pit_ecdf"` asks. Each observation is transformed through
its own leave-one-out predictive distribution, which should leave a
uniform if the model is calibrated, and the plot compares the result
with one. A curve that leaves the band says the predictive spread is
wrong even where the mean is right, which the default check cannot see.

``` r

pp_check(fit, type = "loo_pit_ecdf")
```

![](diagnostics_files/figure-html/ppcloo-1.png)

**Where a residual is left**, rather than how big it is, is what
`type = "error_scatter_avg_vs_x"` shows: the average residual against a
predictor, which is where a missing interaction appears as structure.
And `type = "intervals"` draws a predictive interval per observation, so
systematic under-coverage is visible as a run of points outside their
intervals.

The rest apply where the response allows, and say so when it does not: a
rootogram wants counts, a bar plot wants discrete values, a calibration
plot wants a binary outcome. `type = "km_overlay"` is the one written
for censored data, overlaying the replicate survival curves on the
observed Kaplan-Meier curve; it takes the censoring indicator from the
fit’s own response, and needs the *ggfortify* package installed
alongside *bayesplot*.

### Calibration

The useful question for a binary outcome is whether the predicted
probabilities mean what they say: among the patients the model gave a
30% chance of dying, did about 30% die?

``` r

pp_check(fit, type = "loo_calibration")
```

![](diagnostics_files/figure-html/calibration-1.png)

A line on the diagonal means the probabilities are calibrated, and this
one sits close to it across the whole range. The dots along the bottom
show where the predictions fall, so a departure over a stretch holding
few patients can be discounted.

A line above the diagonal at the left and below it at the right means
the predictions are too extreme, which is the usual failure of an
overfitted model; the opposite pattern means they are too timid, which
is what heavy shrinkage produces.

Each patient is judged here against a probability estimated without
them, which is what separates this from the in-sample version,
`type = "calibration"`. That one reads optimistically, by however much
the model has fitted the individual patients rather than the pattern.
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns the
probabilities themselves, for binning them some other way.

### Residuals

For a continuous outcome, residuals plotted against fitted values read
the way they would for a linear model, with one difference. A forest
shrinks its predictions toward the overall mean, so a mild negative
trend is expected even when the model is correct: the highest fitted
values are pulled down and the lowest pulled up. A strong slope is not
expected, and a fan shape means the spread of the outcome depends on the
predictors, which
[`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
models directly.

For a binary outcome, the residuals take two values for any given fitted
probability and the plot is not informative; the calibration check above
is what to use instead.

### How Much the Model Explains (`r2()`)

``` r

performance::r2(fit)
#> # Bayesian R2 with Compatibility Interval
#> 
#>   Conditional R2: 0.172 (95% CI [0.138, 0.204])
```

This is the Bayesian \\R^2\\ of Gelman et al. ([2019](#ref-gelman2019)),
computed from the posterior rather than from a single fit, so it comes
with an interval and cannot exceed 1. For a binary outcome it is bounded
well below 1 by the outcome’s own randomness, so it is best read as a
relative measure rather than against any absolute standard.

## A Checklist

Before reporting anything from a fit, we run it with `chains = 4` and
pass it to
[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md),
reading what that says to change: it applies the thresholds, names the
rows that are exempt from them, and separates a warmup that ended too
early from chains that have each genuinely settled somewhere different.
We then run
[`pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)
and look for systematic differences, checking in particular that the
replicates respect any bound the outcome has (e.g., a count that cannot
go below zero), and for a binary outcome we read a calibration plot
instead, since the predictive check has little to say there.

## Where to Go Next

[`vignette("comparison")`](https://ngreifer.github.io/bartisan/articles/comparison.md)
covers choosing between models once each of them fits, and
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
covers what to do when the posterior predictive check says the family is
wrong.

## References

Gelman, Andrew, Ben Goodrich, Jonah Gabry, and Aki Vehtari. 2019.
“R-Squared for Bayesian Regression Models.” *The American Statistician*
73 (3): 307–9. <https://doi.org/10.1080/00031305.2018.1549100>.

Vehtari, Aki, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
Paul-Christian Bürkner. 2021. “Rank-Normalization, Folding, and
Localization: An Improved r-Hat for Assessing Convergence of MCMC.”
*Bayesian Analysis* 16 (2): 667–718.
<https://doi.org/10.1214/20-BA1221>.
