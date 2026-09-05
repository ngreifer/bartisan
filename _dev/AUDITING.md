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
down, and for a year it was false by a factor of five. A test that fits with
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
in seconds. `tests/testthat/test-invariants.R` now holds the first two of these.

### 2. Range and sign checks on what is reported

A discrete family's log likelihood is a sum of logs of probabilities and cannot
be positive. The broken fits reported **+37.5** where every comparable fit
reported about **-830**. That single check, applied to every family the package
has, is a dozen lines and would have caught this without anyone understanding
the cause.

The same shape applies elsewhere: a fitted probability in [0, 1], an effective
sample size no larger than the number of draws, a variance not negative, an
R-hat not below 1. None of these needs a truth and all of them fail loudly.

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

**The rest of the matrix**, after a second sweep of 28 cells over seven
families crossed with plain, `vc()`, random effects, and `vc()` with random
effects. Every one passes the two invariants.

| | plain | `vc()` | random effects | `vc()` + ranef |
|---|---|---|---|---|
| gaussian | yes | yes | yes | yes |
| binomial, both links | **yes** | **yes** | **yes** | **yes** |
| ordinal, logit and cloglog | **yes** | **yes** | **yes** | **yes** |
| negbin | **yes** | **yes** | **yes** | **yes** |
| poisson | **yes** | **yes** | **yes** | **yes** |
| zi_poisson | **yes** | **yes** | **yes** | **yes** |
| zi_negbin | **yes** | **yes** | **yes** | **yes** |
| `gaussian_ls` | **yes** | **yes** | **yes** | **yes** |
| multinomial | yes | refused, by design | ? | n/a |
| the survival families | yes | ? | ? | ? |
| `Gamma`, `Gamma_ls`, `Beta`, `ordbeta`, `dpm` | yes | ? | ? | ? |

Only the binomial row carries known-truth recovery as well; the rest are the two
invariants, which is what caught the original bug and is much the cheaper half.

One cell flagged on the first pass and did not survive scrutiny: `zi_negbin`
under `vc()` with random effects disagreed between the augmented and
unaugmented fits at 250 draws on one chain, and agreed at 800 draws on two. Four
additive predictors and a short chain is where Monte Carlo error looks like a
bug, and the remedy is to re-run the cell rather than to widen the tolerance
until it passes.

The remaining question marks are where a hook can still be dropped with nothing
to say so. `poisson` and `gaussian_ls` are worth noting as the control: neither
is on the augmentation list, and their two fits agree *exactly* rather than
approximately, which is what a no-op should look like.

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

`test-invariants.R` now pins the two hooks that matter for any wrapper.

### 5. Instrumented counters, for the hooks that are about timing

`before_forest()` is a notification, not a computation: nothing checks its
result, so nothing notices when it does not happen. A debug build that counts
calls per sweep and asserts that an augmented family's refresh happened once per
forest would have failed immediately and pointed straight at the cause rather
than the symptom.

This is the general remedy for any hook whose contract is "you will be told when
X happens": the only way to test it is to count.

### 6. Simulation-based calibration, for the whole posterior

Everything above tests a point estimate. The interval is a separate claim and
takes its own instrument: draw a parameter from the prior, simulate data from
it, fit, and record where the truth falls in the posterior. Over many draws the
ranks are uniform if and only if the sampler targets the right posterior.

**Done, in `_dev/sbc.R`, and it passes**: 300 replicates, chi-square 4.2 on 9
degrees of freedom for uniformity of the ranks, p = 0.90, and the 95% interval
covering at 0.957. The sub-nominal coverage measured on the benches is therefore
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
  augmenting does not move the posterior; both wrappers forward the two hooks
  whose absence caused this. Four families, both structures.
- `tests/testthat/test-varying.R`: recovery of a known effect through `vc()`
  with an augmented family.

## What to do next, in order

1. Fill in the family x structure matrix above with recovery tests. Each cell is
   a place a hook can be dropped and nothing will say so.
2. Turn the remaining documented invariants in § 1 into tests. They are already
   written down as claims; they only need asserting.
3. Add the range checks of § 2 for every family, which is a loop over the family
   list and one comparison.
4. Simulation-based calibration for two or three families, to settle the
   interval question rather than measuring it obliquely.
