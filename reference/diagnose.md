# Check whether a fit converged and mixed

Reports the convergence and mixing diagnostics worth looking at before
anything is read off a fit, and says what to do about whichever of them
fall short. Everything is computed from the stored draws, so no other
package is needed.

## Usage

``` r
diagnose(object, rhat_max = 1.01, ess_min = 400, ...)

# Default S3 method
diagnose(object, rhat_max = 1.01, ess_min = 400, ...)

# S3 method for class 'bartisan_fit'
diagnose(object, rhat_max = 1.01, ess_min = 400, ...)

# S3 method for class 'bartisan_effect'
diagnose(object, rhat_max = 1.01, ess_min = 400, ...)
```

## Arguments

- object:

  a `<bartisan_fit>` object, the output of a call to
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md),
  or a `<bartisan_effect>` object, the output of a call to
  [`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md).

- rhat_max:

  `numeric`; the largest R-hat treated as acceptable. Default is 1.01,
  the threshold of Vehtari et al. (2021); 1.1 was the older convention
  and is now considered too permissive.

- ess_min:

  `numeric`; the smallest effective sample size treated as acceptable,
  for the bulk and the tail alike. Default is 400, Vehtari et al.'s
  recommendation of 100 per chain at four chains, which is about what it
  takes for the Monte Carlo error of an interval endpoint to be small
  next to the posterior's own width.

- ...:

  ignored.

## Value

A `<bartisan_diagnosis>` object, a list with a
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

### Diagnosing an Estimand

Called on the output of
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md),
the table has one row per reported quantity rather than one per sampled
parameter, and everything else reads the same way.

It is worth doing rather than inferred from the fit's own table, because
the two can disagree. An estimand is a contrast, and a contrast can mix
badly where the function it is a contrast of mixes well: under the
default splitting prior a draw that gives the treatment no rule puts the
contrast at exactly zero, and the sampler can stay there for a long run
while the other predictors keep the fitted function moving. Nothing in
the fit's table shows that, since no parameter the sampler draws is
stuck. The check reports the share of draws sitting at the atom when
there is one.

This complements the fit's diagnosis rather than replacing it. Chains
that have settled on different fitted functions can still agree about an
average over them, so an estimand that looks converged is not on its own
evidence that the sampler did its job.

## Details

This is where the diagnostics are computed, and the only place:
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
does not carry a table of its own, because the per-observation
statistics cost more than the sampling and would then be recomputed
here. With a [future](https://CRAN.R-project.org/package=future) plan in
place the pass is spread over the workers, and it reports progress
through [progressr](https://CRAN.R-project.org/package=progressr) the
way the sampler does.

### What Is Reported

One row per scalar the sampler draws (the log likelihood, the nuisance
parameters of the family, and the scale of each random-effect term),
plus two rows for the additive predictor and two for each set of group
intercepts: one summarizing the worst 5% of observations or levels, and
one for their average.

Both are reported because they routinely disagree, and neither on its
own is the fit. An average over a thousand observations hides the ones
that have not converged, which is what the worst 5% is there to show;
but a forest settles the level of the fitted function within a sweep or
two and takes much longer to settle which observation gets which share
of it, so the worst 5% overstates the trouble for anything averaged. The
average can carry close to one effective draw for every draw kept while
individual observations carry a handful. Which row binds depends on what
is being reported: an average or a contrast of averages is governed by
the average row, and a prediction for one observation by the other.

`rhat` is split-R-hat (Gelman and Rubin, as revised in Gelman et al.
2013): every chain is halved and the halves are compared, so drift
inside a chain counts as disagreement rather than hiding inside a chain
mean. With one chain it is computed by splitting that chain into
segments, which detects drift but cannot detect two chains settling in
different places; that is why one chain draws a warning of its own.

`rhat_late` is that same statistic computed on the second half of the
retained draws alone, and it separates the two reasons chains disagree
by running the experiment rather than by testing for it. Discarding the
early retained draws is exactly what a longer warmup would have done, so
if R-hat is high overall and acceptable late, warmup ended too early. If
it stays high late, the chains have each settled somewhere different and
a longer warmup will not help.

A within-chain drift statistic would answer that question more directly
and cannot be made to work at BART's autocorrelation. Three versions
were calibrated against stationary series where by construction there is
nothing to find, and each failed in one of the two available directions:
too many false alarms on a converged chain, or too little power against
a trend that is really there. A forest is sticky enough to sit where all
three fail, so there is no threshold to pick and the statistic is not
offered.
[`vignette("implementation")`](https://ngreifer.github.io/bartisan/articles/implementation.md)
has the calibrations.

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

### R-hat Needs Effective Draws

R-hat is a ratio of two variance estimates taken from the same draws, so
with few effective draws it sits above 1 whether or not anything is
wrong, and how far above depends on how many chains are being compared.
Against a stationary autoregressive series, where every chain has the
same distribution by construction and there is nothing at all to find,
R-hat averages \\1 + m/S\\ for \\m\\ chains carrying \\S\\ effective
draws between them: 1.026, 1.050, 1.103 and 1.205 at 80 effective draws
over two, four, eight and sixteen chains, against 1.025, 1.050, 1.100
and 1.200 from the formula.

`rhat_max` is therefore not a threshold a quantity can be held to at any
effective sample size. Four chains need 400 effective draws before 1.01
is even the average of R-hat's null, which is where the pairing of the
two defaults comes from and why Vehtari et al. (2021) give them
together, and sixteen chains need 1600 for the same 1.01. When a
quantity fails R-hat while carrying fewer than that, the checks report
that the number cannot be read yet rather than a disagreement it is not
entitled to claim, and the advice sends the reader to the effective
sample size instead.

### The Leaf Scale Is Left Out

`sigma_mu` is deliberately absent from the table. It mixes badly and not
for a reason this package can fix. It is still in `fit$sigma_mu` and
still reaches
[as_draws()](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
for anyone who wants to look.

### What to Do About Poor Mixing

The advice the print method gives follows from which statistic failed,
and the order matters because the fixes are not interchangeable.

**One chain** comes first, because nothing else can be diagnosed
properly until there are several. `chains = 4` is the setting to reach
for, and with [future](https://CRAN.R-project.org/package=future)
installed and a parallel backend used, the chains run in parallel, so it
usually costs little wall clock.

**R-hat elevated but acceptable on the late draws** says warmup ended
too early, so `num_burn` is the one to raise; raising `num_draws`
instead adds draws from a distribution the sampler has not reached yet.
**R-hat elevated on the late draws too** says the chains have each
settled somewhere different. Raise `num_burn` and `num_draws` together,
and if that does not settle it, reduce `num_trees` (a smaller forest has
fewer ways to represent the same fit, so the sampler has less room to
wander between them) and check the family, because a likelihood that
fits badly can produce a posterior with no single place to be.

**More chains and longer chains are not interchangeable.** Effective
sample size depends on the total number of draws and not on how they are
divided between chains, so twice as many chains and twice as long a
chain buy the same amount of it; with a parallel backend the chains are
the cheaper of the two up to the number of workers, though each chain
pays its own warmup. R-hat is not symmetric in the same way, because it
compares chains against each other: lengthening a chain drives it toward
1, while adding chains at a fixed total leaves each chain with less to
say and drives it up. A fit failing on R-hat therefore wants longer
chains rather than more of them, and a fit failing only on effective
sample size can have either.

**An effective sample size that is low while R-hat is fine** is the
benign case, and it wants only more draws. Note that thinning does not
help: `num_thin` discards draws that were already paid for, so it lowers
the effective sample size per unit of time and is worth it only when
storing the draws is the binding constraint. **A low tail effective
sample size with the bulk fine** says the posterior mean is fine and the
interval endpoints are not, so raising `num_draws` is warranted when
intervals are what gets reported.

## References

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., &
Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). Chapman and
Hall/CRC.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Buerkner, P.-C.
(2021). Rank-normalization, folding, and localization: an improved
\\\hat{R}\\ for assessing convergence of MCMC. *Bayesian Analysis*,
16(2), 667–718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

## See also

[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
for the settings the advice names;
[as_draws()](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
for handing the draws to
[bayesplot](https://CRAN.R-project.org/package=bayesplot) or
[posterior](https://CRAN.R-project.org/package=posterior);
[`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md)
for the fuller treatment, including posterior predictive checks

## Examples

``` r
data("rhc")
set.seed(123)

# Two chains, both deliberately short, so that there is something to report
fit <- bartisan(death ~ . - days, data = rhc, num_trees = 10, chains = 2,
                num_burn = 50, num_draws = 50, verbose = FALSE)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` to choose another, which also silences this message.

# The table, the checks, and what to do about whichever of them failed
diagnose(fit)
#> Convergence and mixing
#> 
#>                             quantity  rhat rhat_late ess_bulk ess_tail
#>                               loglik 1.281     1.791        6       45
#>                           splits.eta 1.008     1.153       31       48
#>  eta.eta (average over observations) 1.004     1.161       77       64
#>   eta.eta (worst 5% of observations) 1.488     1.728        4       16
#> 
#> ✔ 2 chains, 100 draws kept in total
#> ✖ R-hat is above 1.01 for loglik
#> ✖ That R-hat rests on only 6 effective draws, where 2 chains average 1.346 even
#>   when they agree
#> ℹ A longer warmup is not the fix: R-hat stays high on the second half of the
#>   draws alone as well
#> ✔ The chains agree about the size of the forest
#> ✖ Bulk ESS is 4 for eta.eta (worst 5% of observations), below 400
#> ✖ Tail ESS is 16 for eta.eta (worst 5% of observations), below 400
#> ℹ The chains disagree about individual observations and agree about their
#>   average (R-hat 1.00, 77 effective draws)
#> ℹ Per-draw efficiency is lowest for eta.eta (worst 5% of observations), which
#>   carries 4.3 effective draws per hundred kept
#> 
#> What to do
#> 
#> • Raise `num_draws`, which was `50`. R-hat is above the threshold for a
#>   quantity that carries too few effective draws for the threshold to mean
#>   anything: with this many chains it would sit about where it does even if the
#>   chains agreed exactly, as the check above reports. Effective sample size is
#>   what makes it readable, and that grows with the total number of draws; using
#>   fewer chains lowers the bar as well, since R-hat's null rises with the number
#>   of chains being compared.
#> • If that does not settle it, reduce `num_trees`, which was `10`. A smaller
#>   forest has fewer ways to represent the same fit, so the sampler has less room
#>   to move between them.
#> • Then check the family. A likelihood that fits the data badly can give a
#>   posterior with no single place to be; `bayesplot::pp_check()` is the
#>   diagnostic.
#> • Note that the chains disagree about the fitted values of individual
#>   observations and not about their average, which is the usual shape of this in
#>   a forest. What that means for an estimand cannot be read off this table
#>   either way, since an estimand is a contrast and a contrast can mix badly
#>   where the function it contrasts mixes well. Compute it: `diagnose()` takes
#>   the output of `estimate_effect()`, and `posterior::as_draws()` hands the
#>   draws to `posterior::summarise_draws()` for anything else.

# A stricter effective sample size, which is what an interval endpoint needs
# and a posterior mean does not
diagnose(fit, ess_min = 1000)
#> Convergence and mixing
#> 
#>                             quantity  rhat rhat_late ess_bulk ess_tail
#>                               loglik 1.281     1.791        6       45
#>                           splits.eta 1.008     1.153       31       48
#>  eta.eta (average over observations) 1.004     1.161       77       64
#>   eta.eta (worst 5% of observations) 1.488     1.728        4       16
#> 
#> ✔ 2 chains, 100 draws kept in total
#> ✖ R-hat is above 1.01 for loglik
#> ✖ That R-hat rests on only 6 effective draws, where 2 chains average 1.346 even
#>   when they agree
#> ℹ A longer warmup is not the fix: R-hat stays high on the second half of the
#>   draws alone as well
#> ✔ The chains agree about the size of the forest
#> ✖ Bulk ESS is 4 for eta.eta (worst 5% of observations), below 1000
#> ✖ Tail ESS is 16 for eta.eta (worst 5% of observations), below 1000
#> ℹ The chains disagree about individual observations and agree about their
#>   average (R-hat 1.00, 77 effective draws)
#> ℹ Per-draw efficiency is lowest for eta.eta (worst 5% of observations), which
#>   carries 4.3 effective draws per hundred kept
#> 
#> What to do
#> 
#> • Raise `num_draws`, which was `50`. R-hat is above the threshold for a
#>   quantity that carries too few effective draws for the threshold to mean
#>   anything: with this many chains it would sit about where it does even if the
#>   chains agreed exactly, as the check above reports. Effective sample size is
#>   what makes it readable, and that grows with the total number of draws; using
#>   fewer chains lowers the bar as well, since R-hat's null rises with the number
#>   of chains being compared.
#> • If that does not settle it, reduce `num_trees`, which was `10`. A smaller
#>   forest has fewer ways to represent the same fit, so the sampler has less room
#>   to move between them.
#> • Then check the family. A likelihood that fits the data badly can give a
#>   posterior with no single place to be; `bayesplot::pp_check()` is the
#>   diagnostic.
#> • Note that the chains disagree about the fitted values of individual
#>   observations and not about their average, which is the usual shape of this in
#>   a forest. What that means for an estimand cannot be read off this table
#>   either way, since an estimand is a contrast and a contrast can mix badly
#>   where the function it contrasts mixes well. Compute it: `diagnose()` takes
#>   the output of `estimate_effect()`, and `posterior::as_draws()` hands the
#>   draws to `posterior::summarise_draws()` for anything else.
```
