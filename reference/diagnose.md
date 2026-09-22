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
while the other predictors keep the fitted function moving. The check
reports the share of draws sitting at the atom when there is one. Chains
that have settled on different fitted functions can still agree about an
average over them, so this complements the fit's own diagnosis rather
than replacing it.

## Details

This is where the diagnostics are computed, and the only place:
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
does not carry a table of its own, because the per-observation
statistics cost more than the sampling and would then be recomputed
here. With a [future](https://CRAN.R-project.org/package=future) plan in
place the pass is spread over the workers, and it reports progress
through [progressr](https://CRAN.R-project.org/package=progressr) the
way the sampler does.

### The Rows of the Table

One row per scalar the sampler draws (the log likelihood, the nuisance
parameters of the family, and the scale of each random-effect term),
plus two rows for the additive predictor and two for each set of group
intercepts: one summarizing the worst 5% of observations or levels, and
one for their average.

Both are reported because they routinely disagree, and which one binds
depends on what is being reported. A forest settles the level of the
fitted function within a sweep or two and takes much longer to settle
which observation gets which share of it, so an average or a contrast of
averages is governed by the average row and a prediction for one
observation by the worst 5%.

An
[`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
fit with three or more categories is the exception, and which of the two
rows to read is different there. Only the differences between the
cutpoints and the additive predictor are identified, so the draws are
recorded in the chart where the predictor has mean zero over the fitted
sample and every cutpoint is free (see
[`bartisan-families`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)).
The average over observations is then zero in every draw, which leaves
it nothing to diagnose, and it is reported as `NA`. The level of the
fitted function has not gone anywhere: the sampler pins the first
threshold, so `aux.cut1` is that level rather than a cutpoint, and it is
the row to read wherever the level is what matters, as it is for a
probability in the lowest categories. It is usually the slowest row in
such a fit, and the most pessimistic one, since it carries the level on
its own where every quantity computed from the draws mixes the level
with faster-moving ones.

`rhat` is split-R-hat, so drift inside a chain counts as disagreement
rather than hiding inside a chain mean. `rhat_late` is that same
statistic on the second half of the retained draws alone, which
separates the two reasons chains disagree: if R-hat is high overall and
acceptable late, warmup ended too early, and if it stays high late, the
chains have each settled somewhere different. `ess_bulk` and `ess_tail`
are rank-normalized effective sample sizes, reported separately because
a chain can be ample for a posterior mean and nowhere near enough for an
interval endpoint. The forest itself is checked through the total number
of splitting rules at each draw, since chains that disagree about how
large the forest is are exploring different tree structures.

That last row is graded apart from the others, and the distinction is
worth understanding before acting on either. A sum of trees represents
one function through many different partitions, so the number of
splitting rules is not pinned down by the fit the way a fitted value is:
two chains can agree to three figures about every value of the additive
predictor while using forests of different sizes. The quantities a fit
reports are integrals over the tree structure, so their convergence is a
separate question from its convergence. The checks on R-hat and
effective sample size therefore read the reported quantities, and the
splitting rules get a check and a remedy of their own. It is separated
rather than suppressed, because it does bind on anything computed from
the split counts themselves, which is
[`variable_importance()`](https://ngreifer.github.io/bartisan/reference/variable_importance.md)
and
[`vignette("importance")`](https://ngreifer.github.io/bartisan/articles/importance.md).

The grading rests on how the model is parameterized and not on that row
being the worst one, which it usually is not. Measured over 144 fits
spanning three families, hard and soft rules and sample sizes from 500
to 8000, the splitting rules carried the highest R-hat in 8 of them; a
reported quantity carried it in the other 136. Its R-hat runs above the
averaged predictor's by 0.13 on a Gaussian fit, 0.05 on a probit one and
0.03 on an ordinal one, which is real and small. **The wide gap in a fit
is not between the reported rows and this one, but within the reported
rows**: the predictor averaged over observations carried a median of 93
times the effective sample size of its own worst 5%, 7815 against 52.
Which of those two governs a given summary is the question the table is
for, and the note about individual observations against their average is
the line to read.

The leaf scale `sigma_mu` is left out of the table, and is in
`fit$sigma_mu` and
[`as_draws()`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
for anyone who wants to look.

### Setting `rhat_max` and `ess_min`

R-hat is a ratio of two variance estimates taken from the same draws, so
with few effective draws it sits above 1 whether or not anything is
wrong, and how far above depends on how many chains are being compared.
`rhat_max` is therefore not a threshold a quantity can be held to at any
effective sample size: four chains need about 400 effective draws before
1.01 is even the average of R-hat's null, which is where the pairing of
the two defaults comes from, and sixteen chains need about 1600 for the
same 1.01. Where a quantity fails R-hat while carrying fewer than that,
the checks report that the number cannot be read yet and send the reader
to the effective sample size instead.

### Remedies for Poor Mixing

The advice the print method gives follows from which statistic failed,
and the order matters because the fixes are not interchangeable.

One chain comes first, since nothing else can be diagnosed properly
until there are several; `chains = 4` is the setting to reach for, and
with [future](https://CRAN.R-project.org/package=future) installed and a
parallel backend in use it usually costs little wall clock. R-hat
elevated but acceptable on the late draws says warmup ended too early,
so `num_burn` is the one to raise. R-hat elevated on the late draws too
says the chains have each settled somewhere different: raise `num_burn`
and `num_draws` together, and failing that reduce `num_trees` and check
the family, since a likelihood that fits badly can produce a posterior
with no single place to be.

Effective sample size depends on the total number of draws rather than
on how they are divided between chains, where R-hat compares chains
against each other, so a fit failing on R-hat wants longer chains rather
than more of them and a fit failing only on effective sample size can
have either. A low effective sample size with R-hat fine is the benign
case and wants only more draws; a low tail effective sample size with
the bulk fine says the posterior mean is sound and the interval
endpoints are not.

[`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md)
works all of this through on a fit.

## See also

[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
for the settings the advice names;
[`as_draws()`](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
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
#> ℹ Set `family` explicitly to silence this message.

# The table, the checks, and what to do about whichever of them failed
diagnose(fit)
#> Convergence and mixing
#> 
#>                             quantity  rhat rhat_late ess_bulk ess_tail
#>                               loglik 1.456     1.161        5       24
#>                           splits.eta 1.141     1.834       11       13
#>  eta.eta (average over observations) 0.990     1.000       79      117
#>   eta.eta (worst 5% of observations) 1.386     1.468        5       17
#> 
#> ✔ 2 chains, 100 draws kept in total
#> ✖ R-hat is above 1.01 for loglik
#> ✖ That R-hat rests on only 5 effective draws, where 2 chains average 1.441 even
#>   when they agree
#> ℹ A longer warmup is not the fix: R-hat stays high on the second half of the
#>   draws alone as well
#> ✖ The chains disagree about how many splitting rules the forest has (R-hat
#>   1.14)
#> ✖ Bulk ESS is 5 for loglik, below 400
#> ✖ Tail ESS is 17 for eta.eta (worst 5% of observations), below 400
#> ℹ The chains disagree about individual observations and agree about their
#>   average (R-hat 0.99, 79 effective draws)
#> ℹ Per-draw efficiency is lowest for loglik, which carries 4.5 effective draws
#>   per hundred kept
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
#> • The forest's own size is a different kind of failure from the others and does
#>   not take the same advice. A sum of trees represents one function through many
#>   different partitions, so two chains can agree about every fitted value while
#>   disagreeing about how many rules they used to get there, and the quantities a
#>   fit reports are integrals over that structure. Raising `num_draws` moves this
#>   row slowly and may not clear the threshold at any affordable length. Act on
#>   it when split counts are themselves what gets reported --
#>   `variable_importance()` and `vignette("importance")` -- and not when fitted
#>   values, predictions or effects are.
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
#>                               loglik 1.456     1.161        5       24
#>                           splits.eta 1.141     1.834       11       13
#>  eta.eta (average over observations) 0.990     1.000       79      117
#>   eta.eta (worst 5% of observations) 1.386     1.468        5       17
#> 
#> ✔ 2 chains, 100 draws kept in total
#> ✖ R-hat is above 1.01 for loglik
#> ✖ That R-hat rests on only 5 effective draws, where 2 chains average 1.441 even
#>   when they agree
#> ℹ A longer warmup is not the fix: R-hat stays high on the second half of the
#>   draws alone as well
#> ✖ The chains disagree about how many splitting rules the forest has (R-hat
#>   1.14)
#> ✖ Bulk ESS is 5 for loglik, below 1000
#> ✖ Tail ESS is 17 for eta.eta (worst 5% of observations), below 1000
#> ℹ The chains disagree about individual observations and agree about their
#>   average (R-hat 0.99, 79 effective draws)
#> ℹ Per-draw efficiency is lowest for loglik, which carries 4.5 effective draws
#>   per hundred kept
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
#> • The forest's own size is a different kind of failure from the others and does
#>   not take the same advice. A sum of trees represents one function through many
#>   different partitions, so two chains can agree about every fitted value while
#>   disagreeing about how many rules they used to get there, and the quantities a
#>   fit reports are integrals over that structure. Raising `num_draws` moves this
#>   row slowly and may not clear the threshold at any affordable length. Act on
#>   it when split counts are themselves what gets reported --
#>   `variable_importance()` and `vignette("importance")` -- and not when fitted
#>   values, predictions or effects are.
#> • Note that the chains disagree about the fitted values of individual
#>   observations and not about their average, which is the usual shape of this in
#>   a forest. What that means for an estimand cannot be read off this table
#>   either way, since an estimand is a contrast and a contrast can mix badly
#>   where the function it contrasts mixes well. Compute it: `diagnose()` takes
#>   the output of `estimate_effect()`, and `posterior::as_draws()` hands the
#>   draws to `posterior::summarise_draws()` for anything else.
```
