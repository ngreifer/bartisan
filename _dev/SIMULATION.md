# How the simulation benches here are built

Written after `_dev/propensity-settings.R`, which worked, so that the next one
does not have to rediscover the same things. `_dev/acic2016.R` follows it.

These are benches for **choosing a setting**, not for validating an estimator.
The question is always "which of these configurations should the package
default to, or the documentation recommend", and everything below follows from
that question rather than from simulation practice in general.

## Score the estimand, not the nuisance

The measurement that started this was of a propensity model's held-out AUC, and
it was the wrong measurement. A propensity score in `bcf()` is a covariate in
the control function: nothing is reported about it, and its job is to let the
control function absorb selection. A setting can predict treatment slightly
worse and leave less confounding in the effect, and only the second matters.

So the metric is bias, root mean squared error, coverage and interval width **on
the effect**, and the fit quality of the nuisance is at most a diagnostic
alongside. The same applies to anything else fitted as a means to an end: score
what the user reads, not what the component does.

## Real covariates, simulated treatment and outcome

A truth has to come from somewhere and real data does not have one. Simulating
everything gives a truth and loses the correlations, the skew, the discreteness
and the near-collinearity that make a real design hard. Simulating only the
treatment and the outcome, on covariates taken from a real dataset, gives both.

`_dev/propensity-settings.R` uses `rhc` (n = 1500, 13 covariates) and cobalt's
`lalonde` (n = 614, 7), which differ in size and in aspect ratio, so a result
that only holds at one shape shows itself. `_dev/acic2016.R` uses the 4802 x 58
covariates from the Collaborative Perinatal Project, which is the same idea at a
scale that exposes different problems.

## Anchor the range before comparing within it

Two configurations that bracket the comparison say whether the differences
inside it are worth anything:

- an **oracle** that is handed the true nuisance, which is the best any estimate
  of it could do; and
- a **baseline** that does without it entirely, which is what the nuisance is
  worth at all.

If the settings under test all sit between two anchors that are themselves close
together, the choice does not matter and no number of replicates will make it
matter. This is the cheapest way to find that out, and it costs two more cells.

## Pair everything

Every configuration sees the **same simulated dataset and the same random
stream**:

```r
for (r in seq_len(reps)) {
  d <- simulate(x, ..., seed = 10000 * design + r)
  for (nm in names(settings)) {
    set.seed(97)
    fit <- fit_one(settings[[nm]], d, covs)
  }
}
```

The differences then carry far less noise than the levels do, and the standard
error to report is of the paired difference against a reference configuration,
not of each mean separately. On this package's own measurements the paired
standard error has repeatedly been three to ten times smaller than the unpaired
one, which is the difference between a conclusion and a shrug.

## Fix what is not being tested

The outcome model is held at one setting throughout while the propensity model
varies. Crossing both would multiply the cells and answer a question nobody
asked. When a second factor does matter, cross it deliberately and say so.

## Scale the signals, not the shapes

`simulate()` scales the treatment and outcome signals to a fixed strength and
solves for the intercept that gives a target prevalence:

```r
shift <- uniroot(function(a) mean(plogis(a + et)) - prevalence, c(-20, 20))$root
```

Without that, "linear" and "nonlinear" selection differ in *how much* there is
to find as well as in shape, and the comparison confounds the two. Scaling first
means a difference between them is a difference of shape.

## Budget before running

A cell is one fit and fits are seconds. Multiply early: configurations times
designs times replicates times seconds, and know the answer before starting. The
propensity bench is 8 x 4 x 25 at 10 to 25 seconds, which is about three hours,
so it goes to `pueue` and the session does something else. `_dev/acic2016.R`
takes the number of settings and sims as arguments for the same reason, and
spreads a small budget across the 77 settings rather than taking the first few,
so that all six factors are represented.

**One timing job at a time.** pueue's default group runs four, and this machine
has four performance cores, so two timing runs at once contend and both are
wrong. Anything measuring seconds waits for the queue.

## Report the breakdown, not only the average

Averaging 77 settings into one number hides the thing the 77 exist for.
`_dev/acic2016.R` reports root mean squared error by overlap, by effect
heterogeneity and by response surface as well as overall, because a
configuration that wins on average by losing badly where overlap is poor is not
a configuration to default to.

## Keep the script, keep the seeds, keep the output

Each bench writes an `.rds` next to itself and is allowlisted in `.gitignore`,
so a claim in the documentation can be traced to the run that produced it. Seeds
are literals in the script. Re-running it reproduces the table.
