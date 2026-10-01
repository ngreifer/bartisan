# Finding the bugs that reading the code does not find

Written after the `before_forest()` bug, which survived several deliberate
audits of this codebase and was found instead by a user fitting a model by hand
and noticing the answer was implausible. That is a bad way to find a bug and it
is worth understanding why nothing else did.

## Why review missed it

**The bug was an absence.** `VaryingCoefficientFamily` forwarded eighteen of the
base class's hooks to the family it wrapped and silently inherited the rest.
Every line that was there was correct. There is no line to point at, no wrong
expression, no bad name, nothing a reviewer reads and doubts. Reading code finds
code that is wrong; it does not find code that is not there.

**The output was plausible.** A treatment effect shrunk toward zero is what a
regularized model is *supposed* to do. There was no crash, no `NaN`, no warning,
no failing test. A shrunken effect from a BART model looks like the literature's
main complaint about BART, which is exactly the shape a reader's prior expects.

**The tests covered every line and no composition.** There are tests for the
augmented families. There are tests for `vc()`. There were none for an augmented
family *under* `vc()`, which is where the two met and one of them was dropped.
Line coverage was not the thing that was missing.

The lesson generalizes: **the bugs that survive review are the ones whose output
is believable.** Finding them takes an oracle, not another reading.

## What would have caught it, in order of cost

### 1. Invariants the package already states

`?bartisan_control` says of `augment` that **"the posterior is the same either
way"**. That is a claim about two samplers for one posterior, it is written
down, and for the four days between the wrapper's first commit (2026-09-01)
and the fix (2026-09-05) it was false by a factor of five. A test that fits with
`augment = TRUE` and `augment = FALSE` and compares would have failed the day
the wrapper was written.

This is the cheapest and best technique available here because the package is
unusually rich in such statements. Others already written down and now worth
testing:

- **Centering**: "every coefficient and every estimand is identical under any
  choice" of `center` in `vc()`.
- **The empty forest**: `gaussian_ls()` with `~ 1` on its scale "is `gaussian()`
  with its drawn `sigma`", to within the difference in the prior.
- **The drawn coding**: at two levels "it restricts nothing", so the estimand
  should not move when it is switched off.
- **Sparsity off**: `sparsity = FALSE` "recovers a uniform prior over
  predictors", which is a `split_prior` a caller can supply explicitly.
- **Chains**: results do not depend on how many workers ran them. This one is
  already tested, and it is the model for the rest.

Every one is a pair of routes to a single answer, needs no known truth, and runs
in seconds. As of 2026-10-01 all five are tests in
`tests/testthat/test-invariants.R`: the centering and the drawn coding are
compared on the ATE and on the per-unit effects, the empty scale forest on the
predictor and on the residual sd, and `sparsity = FALSE` against the written-out
uniform `split_prior` with one seed, where the two fits are the same chain to
1e-10. The chains invariant was already in `test-chains.R`.

### 2. Range and sign checks on what is reported

A discrete family's log likelihood is a sum of logs of probabilities and cannot
be positive. The broken fits reported **+37.5** where every comparable fit
reported about **-830**. That single check, applied to every family the package
has, is a dozen lines and would have caught this without anyone understanding
the cause.

The same shape applies elsewhere: a fitted probability in [0, 1] with a row of
them summing to one, a scale or count or rate that is positive, and finite,
positive effective sample sizes and R-hats. None of these needs a truth and all
of them fail loudly. Two bounds this file once proposed are wrong for the
rank-normalized estimators `diagnose()` uses and were dropped when the checks
became tests on 2026-10-01: bulk ESS can exceed the number of draws for
antithetic chains (3597 on 3200 draws in a `dpm()` fit that day), and split
R-hat can read a little below one (0.997 on a converged average). The checks
run in every cell of the matrix test in `test-invariants.R`, over the fitted
response, the probabilities of the categorical families, the named nuisance
columns and the baseline hazards, and every row of `diagnose()`.

### 3. Recovery of a known truth, across the *cross-product*

Simulate from the model with a known parameter, fit, and check it comes back.
The decisive measurement here was a log-odds effect of exactly 1 with no
confounding: `vc()` with a logit binomial returned 0.17 where every other
configuration returned 0.88 to 0.91.

The important word is cross-product. Families were tested. Structures were
tested. The failure was in a cell of families x structures that nothing visited.
Worth keeping an explicit matrix of which cells are covered:

The binomial row is filled in. Twelve structures crossed with both links, each
fitted with the augmentation on and off, and each checked for a log likelihood
that is not positive and for the two fits agreeing:

| structure | logit | probit |
|---|---|---|
| plain | yes | yes |
| `vc()` | yes | yes |
| random effects | yes | yes |
| `vc()` + random effects | yes | yes |
| offset | yes | yes |
| `vc()` + offset | yes | yes |
| weights, four trials | yes | yes, declines to augment |
| `vc()` + weights | yes | yes, declines to augment |
| categorical predictor | yes | yes |
| `vc()` + categorical | yes | yes |
| missing predictor | yes | yes |
| `vc()` + missing | yes | yes |

All 24 pass. The probit weight cells agree exactly rather than approximately,
which is right: a probit augmentation needs a Bernoulli response and declines
four trials, so both fits are the same fit.

The same twelve, scored against a known log-odds effect of 0.8 and a
random-effect standard deviation of 0.6:

| structure | effect | recovered | cor(u, uhat) | tau |
|---|---|---|---|---|
| plain | 0.731 | 91% | | |
| `vc()` | 0.734 | 92% | | |
| random effects | 0.816 | 102% | 0.964 | 0.681 |
| `vc()` + random effects | 0.794 | 99% | 0.967 | 0.680 |
| offset | 0.746 | 93% | | |
| `vc()` + offset | 0.785 | 98% | | |
| weights | 0.700 | 88% | | |
| `vc()` + weights | 0.703 | 88% | | |

Nothing here is broken. The 88% on the weighted cells is the same shrinkage the
unstructured fits show and is not specific to any wrapper.

**The rest of the matrix** is now a test rather than a sweep:
`test-invariants.R`, "every family keeps its chart and its density under every
structural wrapper", fits 25 families under seven structures and both gates
and asserts two identities that hold draw by draw and need no truth. The density route must
reproduce the recorded log likelihood, which it does only when the recorded
predictor and the recorded nuisance parameters are in one chart; and the
stored forests must replay to the recorded predictor, which they do only when
the leaves were written in that same chart. The structures are plain, `vc()`
with the fixed coding, `vc()` with the drawn coding `bcf()` uses (skipped,
with the refusal named, where the family's leaf target is not quadratic),
random intercepts, `vc()` with random intercepts, an offset, and frequency
weights (refused, with the refusal named, by the two Dirichlet process
families). Every cell runs under hard rules and under the default smoothstep
gate, since the replay of a soft tree goes through membership weights and a
bandwidth the hard arm never touches. Discrete families also get the sign check
of § 2 and every cell the range checks. `mnp()` is in for everything but the
density identity, which is not exact for a simulated density. As of 2026-10-01
every cell passes, 2738 expectations in about two minutes:

| | plain | `vc()` | drawn coding | random effects | `vc()` + ranef |
|---|---|---|---|---|---|
| gaussian, `gaussian_ls` | yes | yes | yes | yes | yes |
| binomial, both links | yes | yes | yes | yes | yes |
| poisson, negbin, `zi_poisson`, `zi_negbin` | yes | yes | yes | yes | yes |
| `Gamma`, `Gamma_ls`, `tweedie` | yes | yes | yes | yes | yes |
| `Beta`, `ordbeta` | yes | yes | yes | yes | yes |
| ordinal, three links | yes | yes | yes | yes | yes |
| `dpm` | yes | yes | yes | yes | yes |
| `weibull_aft`, `loglogistic_aft`, `lognormal_aft`, `dpm_aft`, `ph` | yes | yes | yes | yes | yes |
| `custom_family()` | yes | yes | yes | yes | yes |
| multinomial | yes | refused, by design | refused | yes | refused |

"yes" in the drawn-coding column means the cell is fit when the family's
target is quadratic and refused with the right message when it is not; the
test pins the set of refusals, so a family that changes its target form moves
a name. The offset and weights columns are not shown: every family passes both,
the Dirichlet process families refusing weights as documented.

**Three bugs the matrix found on the day it was written.** The offset arm
failed the replay identity for every family alike: `bartisan()` stored the
intercept as the first column of the combined offset, which with a user offset
is the intercept plus the first observation's offset, so every prediction on
new data from a fit with an offset was shifted by that one value while every
training-data quantity was right. Fixed in `R/bartisan.R`, with a test in
`test-bartisan.R`. The other two were in the first version of the matrix: The `dpm` cell under
the drawn coding, which is what `bcf(family = dpm())` fits, failed the density
identity by 1668 log points: `VaryingCoefficientFamily` never forwarded
`report_shift()`, so the recorded predictor was the raw forest and its level
the coordinate the likelihood does not identify (the `TASKS.md` entry of
2026-10-01 has the measurements). And every `Beta` cell, the plain one
included, failed it by about 2900: `BetaFamily::set_aux()` set the precision
without refreshing the eta-free terms, which carry `lgamma(phi)`, so
`predict(type = "density")` and `loo()` priced every draw with the
constructor's precision. The sampler was right in both cases and every
reported number looked plausible, which is the shape of bug this file is
about.

A third test, "a family with a reporting chart keeps it under every
structural wrapper", covers the three families that record their draws in a
chart other than the sampler's. For `ordinal()` the recorded predictor must
average to zero in every draw; for `dpm()` and `dpm_aft()` the level of the
recorded predictor must have less than half the spread of the raw mixture
mean and must not follow it, which in the sampler's chart it does at a
correlation near minus one.

**Recovery across the matrix** (`_dev/recovery-matrix.R`, 2026-10-01, pueue
task 235, 115 cells, 17 minutes of fitting). Twenty-three families under the
five structures, n = 1000, 50 trees, hard rules, 300 + 500 draws, a link-scale
effect of 0.8 and a known level. The reading rule was written before the run:
a cell off from the rest of its row is the finding, and a row shrunk alike in
every cell is the prior. `multinomial()` is left out for want of a clean
truth on its scale, and `mnp()` with it.

Effect on the link scale, truth 0.8:

| family | plain | vc | drawn | ranef | vc + ranef |
|---|---|---|---|---|---|
| `gaussian` | 0.79 | 0.78 | 0.78 | 0.79 | 0.78 |
| `gaussian_ls` | 0.79 | 0.78 | refused | 0.78 | 0.78 |
| `binomial logit` | 0.89 | 0.87 | 0.84 | 0.89 | 0.88 |
| `binomial probit` | 0.88 | 0.90 | 0.89 | 0.88 | 0.94 |
| `binomial cloglog` | 0.68 | 0.73 | refused | 0.75 | 0.78 |
| `poisson` | 0.78 | 0.75 | refused | 0.76 | 0.73 |
| `negbin` | 0.83 | 0.83 | refused | 0.81 | 0.81 |
| `zi_poisson` | 0.75 | 0.73 | refused | 0.73 | 0.67 |
| `zi_negbin` | 0.82 | 0.71 | refused | 0.84 | 0.67 |
| `Gamma` | 0.78 | 0.77 | refused | 0.75 | 0.74 |
| `Gamma_ls` | 0.78 | 0.78 | refused | 0.76 | 0.75 |
| `Beta` | 0.75 | 0.75 | refused | 0.77 | 0.77 |
| `ordbeta` | 0.79 | 0.86 | refused | 1.00 | 1.00 |
| `tweedie` | 0.86 | 0.81 | refused | 0.84 | 0.80 |
| `ordinal logit` | 0.83 | 0.88 | 0.82 | 0.84 | 0.89 |
| `ordinal probit` | 0.50 | 0.50 | 0.52 | 0.53 | 0.51 |
| `ordinal cloglog` | 0.61 | 0.65 | refused | 0.60 | 0.66 |
| `dpm` | 0.86 | 0.85 | 0.85 | 0.88 | 0.90 |
| `weibull_aft` | 0.80 | 0.79 | refused | 0.81 | 0.80 |
| `loglogistic_aft` | 0.78 | 0.77 | 0.77 | 0.77 | 0.77 |
| `lognormal_aft` | 0.81 | 0.80 | 0.81 | 0.80 | 0.80 |
| `dpm_aft` | 0.78 | 0.79 | 0.79 | 0.80 | 0.81 |
| `ph` | 0.74 | 0.74 | refused | 0.85 | 0.86 |

Level error, mean fitted link minus the true mean (n/a where the level is not identified):

| family | plain | vc | drawn | ranef | vc + ranef |
|---|---|---|---|---|---|
| `gaussian` | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 |
| `gaussian_ls` | 0.01 | 0.01 | refused | 0.01 | 0.01 |
| `binomial logit` | 0.02 | 0.03 | 0.02 | 0.04 | 0.06 |
| `binomial probit` | -0.11 | -0.09 | -0.10 | -0.10 | -0.06 |
| `binomial cloglog` | -0.06 | -0.02 | refused | -0.00 | 0.03 |
| `poisson` | 0.09 | 0.09 | refused | 0.05 | 0.04 |
| `negbin` | 0.05 | 0.04 | refused | -0.00 | -0.01 |
| `zi_poisson` | 0.19 | 0.14 | refused | 0.08 | 0.10 |
| `zi_negbin` | 0.02 | -0.00 | refused | -0.09 | -0.02 |
| `Gamma` | 0.05 | 0.05 | refused | -0.02 | -0.02 |
| `Gamma_ls` | 0.05 | 0.05 | refused | -0.02 | -0.02 |
| `Beta` | -0.00 | 0.00 | refused | 0.02 | 0.03 |
| `ordbeta` | 0.04 | 0.08 | refused | 0.08 | 0.10 |
| `tweedie` | 0.07 | 0.08 | refused | 0.00 | 0.01 |
| `ordinal logit` | n/a | n/a | n/a | n/a | n/a |
| `ordinal probit` | n/a | n/a | n/a | n/a | n/a |
| `ordinal cloglog` | n/a | n/a | refused | n/a | n/a |
| `dpm` | 0.03 | 0.03 | 0.03 | 0.03 | 0.03 |
| `weibull_aft` | 0.09 | 0.09 | refused | 0.02 | 0.02 |
| `loglogistic_aft` | 0.02 | 0.02 | 0.02 | -0.00 | 0.00 |
| `lognormal_aft` | 0.05 | 0.05 | 0.05 | 0.02 | 0.02 |
| `dpm_aft` | -0.01 | 0.00 | -0.00 | -0.03 | -0.03 |
| `ph` | n/a | n/a | refused | n/a | n/a |

The effect's spread across structures within a row is at most 0.07 for 19 of
the 23 families. Of the other four, `binomial("cloglog")` and `ph()` rise
under random intercepts (0.68 to 0.78, 0.74 to 0.86), which is
non-collapsibility: the cells without `(1 | g)` estimate a marginal effect
attenuated by the unmodeled intercepts, and the truth of 0.8 is conditional.
`zi_negbin()` dips under `vc()` (0.71 and 0.67 against 0.82 and 0.84) and
`ordbeta()` reads 1.00 in both random-intercept cells against 0.79 and 0.86.
Both were rerun with three seeds and 800 draws (`_dev/recovery-replicates.R`)
and then followed up (`_dev/recovery-followup.R`; both in `TASKS.md` under the
same date). The `ordbeta` overshoot was the generator's: the row above set the
response to 0 or 1 by a deterministic cut on the predictor, a step the model
represents with logistic boundary probabilities. Drawn from the model's own
mechanism the row reads 0.75, 0.75, 0.78 and 0.78 across the four structures,
and the matrix script now draws it that way. The `zi_negbin` dip is
identification between the count and zero processes: at n = 4000 the gap is
0.03, half a posterior sd, against 0.11 at n = 1000. The ordinal probit row reads 0.50 and the cloglog rows 0.61 to 0.78
because the latent truth is logistic, and 0.8 on that scale is about 0.47 on
the probit's. The level error is under 0.1 in every identified cell except
`binomial("probit")` at about -0.1 in every structure and `zi_poisson()` at
0.19 and 0.14 without random intercepts; both are the same across their row
and so belong to the model, not to a wrapper. The `dpm` cell under the drawn
coding, which is the `bcf(family = dpm())` composition that was wrong, has a
level error of 0.026 against 0.026 in the plain cell.

**The soft-rule arm** (`RECOVERY_GATE=smoothstep`, pueue task 240, 42
minutes, the `ordbeta` row drawn from the model's own boundary mechanism):
every row's spread across structures is at most 0.08 except `ordinal()`
logit at 0.13 and `ph()` at 0.11, the latter the same non-collapsibility as
under hard rules. The two cells flagged under hard rules read normally here,
`ordbeta` under random intercepts at 0.84 and `zi_negbin` under `vc()` with
random intercepts at 0.85, which is what the follow-ups had concluded. Every
family's soft row is within about 0.05 of its hard row except the three ordinal
rows, which read about 30% lower in every structure alike (logit 0.54 to 0.67
against 0.82 to 0.89, probit 0.34 to 0.37 against 0.50 to 0.53). That was not
the gate. The hard arm ran before the `ordbeta` generator was corrected, and
the correction added a `runif()` call that shifted the random stream for every
response drawn after it, so the two arms' ordinal, tweedie and survival rows
were fit to different data. `_dev/ordinal-soft-check.R` reproduced the soft
arm's data step for step and fit both gates on it at 800 and 3200 draws: logit
0.595 hard against 0.606 soft, probit 0.351 against 0.354, with no movement at
the longer chain, and the same agreement on two fresh data sets. The ordinal
effect's posterior sd is 0.12 and the two data sets happened to sit on either
side. Two lessons for the generator: a change to one response's draw moves
every response after it, so a cross-arm comparison wants a seed per response;
and a row that moves alike in every structure between two runs is a data
question before it is a gate question.

Not covered: the recovery matrix's weights and offset arms, which the identity
matrix has and the recovery matrix does not.

### 4. The decorator audit, which is mechanical

The bug has a shape: a class that wraps another and forwards part of its
interface. That shape can be checked without understanding any of the
statistics. Enumerate the base class's virtual hooks, enumerate what each
wrapper overrides, and look at the difference deliberately.

Run on this codebase it takes a minute and it found a second instance
immediately: `LinkedFamily` was missing the same two forwards. That one is not
reachable today, because a family only reaches that wrapper through a link the
engine does not carry natively and an augmented family is only built for the
native links, so nothing is currently wrong. It is fixed anyway, because "not
reachable" is a property of today's composition rules and not of the code.

The first version of the test pinned only those two hooks, and `report_shift()`
went through it. Since 2026-10-01 "every hook a decorator inherits is one it
means to inherit" computes the full difference for `VaryingCoefficientFamily`,
`LinkedFamily` and `RFamily` from the sources and asserts the exact set, with
the reason each inherited hook may be inherited written beside it. A hook
added to `family.h` fails the test until it is forwarded or listed with a
reason. The lists as they stand: the varying-coefficient wrapper inherits only
the blocked sums and likelihood deltas, which are generic over the unit hooks
it forwards; the linked wrapper also inherits `report_shift()`,
`aux_values_shifted()`, `mixture_flat()`, `num_pinned()`, the target form and
the analytic derivatives, each for a stated reason that holds for the six
families `with_link()` accepts; and `RFamily` inherits the base defaults
because it wraps an R function, not a family.

### 5. Instrumented counters, for the hooks that are about timing

`before_forest()` is a notification, not a computation: nothing checks its
result, so nothing notices when it does not happen. A debug build that counts
calls per sweep and asserts that an augmented family's refresh happened once per
forest would have failed immediately and pointed straight at the cause rather
than the symptom.

This is the general remedy for any hook whose contract is "you will be told when
X happens": the only way to test it is to count.

**Done 2026-10-01, and always on.** `Family` counts calls to `before_forest()`
and `update_aux()` through two entry points, `run_before_forest()` and
`run_update_aux()`, which the engine and every wrapper call instead of the hooks
themselves, so the count lands on the innermost family. At the end of every
sweep the engine asks `sweep_delivered()`, which is true when `before_forest()`
ran at least once and `update_aux()` exactly once, and the wrappers answer for
what they wrap. A wrapper that drops either hook now stops the fit on its first
sweep with a message naming the two hooks, instead of shrinking an effect by a
factor of five with every test passing. The cost is two integer increments and
one virtual call per sweep. It fires nowhere in the matrix test, which is every
family under every wrapper; it cannot be shown firing without breaking a
wrapper on purpose, so the test of the counters is that the matrix runs.

### 6. Simulation-based calibration, for the whole posterior

Everything above tests a point estimate. The interval is a separate claim and
takes its own instrument: draw a parameter from the prior, simulate data from
it, fit, and record where the truth falls in the posterior. Over many draws the
ranks are uniform if and only if the sampler targets the right posterior.

**Done, in `_dev/sbc.R`, and it passes**: 300 replicates, chi-square 4.2 on 9
degrees of freedom for uniformity of the ranks, p = 0.90, and the 95% interval
covering at 0.957. That was the logit binomial with its augmentation off, so
the Laplace leaf sampler. On 2026-10-01 the script gained a family argument
and two more arms, each a likelihood with no parameter beyond the predictor so
that the prior is the tree prior alone: `binomial("probit")` with its
augmentation on, which calibrates the latent-variable sampler, and
`poisson()`, which calibrates the exponential-form target. Probit, 200
replicates at n = 400 under hard rules: chi-square 2.5 on 9 df (p = 0.98), 48%
of ranks in the middle half, mean rank 49.1 of 100, coverage 0.960 (SE 0.014).
Poisson, the same design: chi-square 6.0 on 9 df (p = 0.74), 43% of ranks in
the middle half, mean rank 48.8, coverage 0.955 (SE 0.015). The 43% is two
standard errors under 50 on 200 replicates and reads as a posterior a shade
narrow if it reads as anything, with coverage at nominal saying it does not;
a longer run would settle it and nothing here depends on it. All three arms are
calibrated, so the Laplace leaf sampler, the latent-variable augmentation and
the exponential-form target each draw from the posterior they define. The sub-nominal coverage measured on the benches is therefore
the prior doing its job against a truth it does not favor, not a sampler that is
wrong, which is the difference between a credible interval and a confidence
interval.

The obstacle worth knowing about is that this model's prior is **empirical**,
which is standard for BART and makes SBC not quite well-posed: the leaf scale
comes from the response's spread, the residual scale from a regression of the
response on the predictors, and the predictor is centered at an empirical value.
`sigma_mu` can be fixed through `bartisan_control()`, a `binomial()` avoids the
residual scale entirely, and calibrating a *contrast* rather than a level makes
the empirical centering cancel. What is left, the tree prior, is data-free and
replicates exactly. A family with a nuisance parameter, a soft gate and a drawn
sparsity prior are each a further piece of prior to replicate and are not
covered yet.

### 7. Differential testing against another implementation

For the cases another package can fit, its answer is an oracle. `dbarts` and
`stochtree` fit Gaussian and probit BART; where the models agree the posteriors
should agree within Monte Carlo error. The benchmarking in `_dev/` already runs
these side by side for *speed*, which is most of the work; comparing the answers
costs almost nothing more.

### 8. The one that actually worked: a plausibility check on real data

A domain expert fitted a model to real data and knew the answer was wrong. That
is not a technique that scales, but it is not nothing either, and it can be
partly captured: the estimate from `vignette("causal")` on `rhc` has a published
literature around it, and pinning it to a plausible range in a test turns one
person's judgment into something that runs every time.

## What is now in place

- `tests/testthat/test-invariants.R`: discrete log likelihoods are not positive;
  augmenting does not move the posterior; the exact set of hooks each of the
  three wrappers inherits, with reasons; the density and replay identities, the
  sign and range checks, for 25 families under seven structures and both gates;
  the reporting chart of the three families that have one; and the five
  documented invariants of § 1.
- `src/family.h` and `src/model.cpp`: the hook delivery counters of § 5, on in
  every fit.
- `tests/testthat/test-varying.R`: recovery of a known effect through `vc()`
  with an augmented family.
- `tests/testthat/test-dpm.R` and `test-ordinal.R`: the reporting chart
  survives the varying-coefficient wrapper, for `dpm()`, `dpm_aft()` and
  `ordinal()`.

## What to do next, in order

Everything on the list of 2026-10-01 is done: the matrix in both of its forms
and both gates, the documented invariants, the range checks, SBC for three
families, and the hook counters. What is left is smaller and none of it stands
between the package and a submission.

1. Seed each response of `_dev/recovery-matrix.R` separately, so that editing
   one response's generator stops moving the data of every response after it;
   then the hard and soft arms are comparable row for row, which on 2026-10-01
   they were only for the rows drawn before `ordbeta`.
2. Weights and offset arms for the recovery matrix, which the identity matrix
   has and the recovery matrix does not.
3. Differential testing against another implementation (§ 7). `dbarts` fits
   Gaussian and probit BART; where the models agree the posteriors should
   agree within Monte Carlo error, and the speed benches already run them side
   by side.
4. The plausibility pin of § 8: the `rhc` estimate in `vignette("causal")`
   against its literature, as a test.
5. A longer Poisson SBC run, to settle the 43% middle-half share, and the SBC
   arms under soft rules, which `_dev/sbc.R` takes as an argument.
