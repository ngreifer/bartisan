# Check whether a fit converged and mixed

Reports the convergence and mixing diagnostics worth looking at before
anything is read off a fit, and says what to do about whichever of them
fall short. Everything is computed from the stored draws, so no other
package is needed.

## Usage

``` r
diagnose(object, rhat_max = 1.01, ess_min = 400)
```

## Arguments

- object:

  a fit from
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).

- rhat_max:

  the largest R-hat treated as acceptable. The default, `1.01`, is the
  threshold of Vehtari et al. (2021); `1.1` was the older convention and
  is now considered too permissive.

- ess_min:

  the smallest effective sample size treated as acceptable, for the bulk
  and the tail alike. The default, `400`, is Vehtari et al.'s
  recommendation of 100 per chain at four chains, which is about what it
  takes for the Monte Carlo error of an interval endpoint to be small
  next to the posterior's own width.

## Value

An object of class `bartisan_diagnosis`, with a
[`print()`](https://rdrr.io/r/base/print.html) method that shows the
table, the checks and the advice. Its components are

- `table`:

  a data frame with one row per quantity: `rhat`, `rhat_late` (the same
  statistic on the second half of the draws alone), `ess_bulk`,
  `ess_tail`, and `ess_frac`, the bulk effective sample size as a
  fraction of the draws kept.

- `checks`:

  a data frame of `check`, `status` (`"ok"`, `"warn"` or `"note"`) and
  `detail`.

- `advice`:

  a character vector, most important first, empty when everything
  passed.

- `chains`,`draws`:

  how many chains, and how many draws were kept in total.

## What is reported

One row per scalar the sampler draws – the log likelihood, the nuisance
parameters of the family, the scale of each random-effect term – plus
one row for the additive predictor and one for each set of group
intercepts, summarized over their worst 5% of observations or levels
rather than averaged, since an average over a thousand observations
hides the ones that have not converged.

`rhat` is split-R-hat (Gelman and Rubin, as revised in Gelman et al.
2013): every chain is halved and the halves are compared, so drift
inside a chain counts as disagreement rather than hiding inside a chain
mean. With one chain it is computed by splitting that chain into
segments, which detects drift but cannot detect two chains settling in
different places – which is why one chain draws a warning of its own.

`rhat_late` is that same statistic computed on the second half of the
retained draws alone, and it is what separates the two reasons chains
disagree – by running the experiment rather than by testing for it.
Discarding the early retained draws is exactly what a longer warmup
would have done, so if R-hat is high overall and acceptable late, warmup
ended too early. If it stays high late, the chains have each settled
somewhere different and a longer warmup will not help.

A within-chain drift statistic would answer that question more directly
and cannot be made to work at BART's autocorrelation. Three versions
were calibrated against stationary autoregressive series, where by
construction there is nothing to find: taking each half's Monte Carlo
error from that half alone fires 31% of the time at an autocorrelation
of 0.995 against a nominal 5%; taking it from the whole chain holds
specificity under 3% but then misses a linear trend of six standard
deviations three times in four; batch means catch everything and fire
92% of the time on a chain that has converged. A forest is sticky enough
to sit where all three fail, so there is no threshold to pick and the
statistic is not offered.

The rows summarized over observations or levels report the **worst 5%**
boundary rather than the single worst column, and the checks are keyed
to the *share* of columns that failed rather than to that boundary: the
worst of a thousand values is extreme even when every chain has
converged, so a threshold applied to a maximum would condemn every fit.

`ess_bulk` and `ess_tail` are the rank-normalized effective sample sizes
of Vehtari et al. (2021). The tail one is reported separately because a
chain can be ample for a posterior mean and nowhere near enough for an
interval endpoint.

The forest itself is checked too, through the total number of splitting
rules in it at each draw. Chains that disagree about how large the
forest is are exploring different tree structures, and no generic MCMC
diagnostic can see that, because none of them looks at the forest.

## The leaf scale is left out

`sigma_mu` is deliberately absent, as it is from `fit$rhat`. It mixes
badly and not for a reason this package can fix: on the same data the
same quantity comes out at R-hat 1.12 in
[dbarts](https://CRAN.R-project.org/package=dbarts) and 1.16 in
`stochtree`, both of which draw it a different way, against 1.19 here.
It is a hyperparameter nobody reports, and its disagreement between
chains does not reach the fitted function – on those same fits the
additive predictor has R-hat 1.00 and thousands of effective draws. It
is still in `fit$sigma_mu` and still reaches
[as_draws()](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
for anyone who wants to look.

## What to do about poor mixing

The advice the print method gives follows from which statistic failed,
and the order matters because the fixes are not interchangeable.

- **One chain.** Nothing else can be diagnosed properly. `chains = 4` is
  the first thing to set, and with
  [future](https://CRAN.R-project.org/package=future) installed the
  chains run in parallel, so it usually costs little wall clock.

- **R-hat elevated, acceptable on the late draws.** Warmup ended too
  early: raise `num_burn`. Raising `num_draws` instead adds draws from a
  distribution the sampler has not reached yet.

- **R-hat elevated on the late draws too.** The chains have each settled
  somewhere different. Raise `num_burn` and `num_draws` together, and if
  that does not settle it, reduce `num_trees` – a smaller forest has
  fewer ways to represent the same fit, so the sampler has less room to
  wander between them – and check the family, because a likelihood that
  fits badly can produce a posterior with no single place to be.

- **Effective sample size low, R-hat fine.** The benign case. Raise
  `num_draws`. Thinning does not help: `num_thin` discards draws that
  were already paid for, so it lowers the effective sample size per unit
  of time and is worth it only when storing the draws is the binding
  constraint.

- **Tail effective sample size low, bulk fine.** The posterior mean is
  fine and the interval endpoints are not. Raise `num_draws` if
  intervals are what gets reported.

## References

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., &
Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). Chapman and
Hall/CRC.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C.
(2021). Rank-normalization, folding, and localization: an improved
\\\hat{R}\\ for assessing convergence of MCMC. *Bayesian Analysis*,
16(2), 667–718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

## See also

[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
for the settings the advice names,
[as_draws()](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
for handing the draws to
[bayesplot](https://CRAN.R-project.org/package=bayesplot) or
[posterior](https://CRAN.R-project.org/package=posterior), and
[`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md)
for the fuller treatment including posterior predictive checks.

## Examples

``` r
set.seed(1)
n <- 200
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
d$y <- d$x1 + rnorm(n)

fit <- bartisan(y ~ x1 + x2, data = d, family = gaussian(),
                control = bartisan_control(chains = 2, num_trees = 10,
                                           num_burn = 100, num_draws = 100,
                                           verbose = FALSE))

diagnose(fit)
#> Convergence and mixing
#> 
#>                            quantity  rhat rhat_late ess_bulk ess_tail
#>                              loglik 1.047     1.046       49       65
#>                           aux.sigma 1.040     1.051      148      154
#>                          splits.eta 1.117     1.091       25       48
#>  eta.eta (worst 5% of observations) 1.203     1.198        8       20
#> 
#> ✔ 2 chains, 200 draws kept in total
#> ✖ above 1.01 for loglik
#> ℹ not the whole story: R-hat stays high on the second half of the draws alone,
#>   so the chains disagree rather than merely start badly
#> ✖ the chains disagree about how many splitting rules the forest has (R-hat
#>   1.12)
#> ✖ 8 for eta.eta (worst 5% of observations), below 400
#> ✖ 20 for eta.eta (worst 5% of observations), below 400
#> ℹ eta.eta (worst 5% of observations) carries 3.8 effective draws per hundred
#>   kept
#> 
#> What to do
#> 
#> Raise `num_burn` and `num_draws` together. R-hat stays high even on the second
#> half of the draws alone, so the chains have each settled somewhere different
#> rather than merely started badly.
#> If that does not settle it, reduce `num_trees`. A smaller forest has fewer ways
#> to represent the same fit, so the sampler has less room to move between them.
#> Then check the family. A likelihood that fits the data badly can give a
#> posterior with no single place to be; `pp_check()` is the diagnostic.
```
