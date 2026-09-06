# What three mixing papers say, and what to do about it here

Three papers, read 2026-09-05:

- Kim and Ročková (2025), *On mixing rates for Bayesian CART*, Electronic
  Journal of Statistics 19(2). `kimMixingRatesBayesian2025`
- Ronen, Saarinen, Tan, Duncan and Yu (2022), *A Mixing Time Lower Bound for a
  Simplified Version of BART*, arXiv:2210.09352. `ronenMixingTimeLower2022`
- Tan, Ronen, Saarinen and Yu (2026), *On the Computational Efficiency of
  Bayesian Additive Regression Trees: An Asymptotic Analysis*,
  arXiv:2406.19958. `tanComputationalEfficiencyBayesian2026`

They are one research programme by an overlapping group, and they agree.

## What they establish

**Mixing gets worse as the sample grows.** Ronen et al. give the first lower
bound, for a simplified BART with one tree and a subset of the moves, and it
grows *exponentially* in the number of observations. Tan et al. give a lower
bound for the hitting time of the high-posterior-density region under a
near-default sampler, growing at least like the square root of the sample size,
and attribute it to a posterior whose modes get harder to travel between as the
data grow.

**The consequence they name is ours.** Tan et al.: if default iteration counts
are used, "the posterior approximation of the BART sampler may deviate from the
true posterior at large n, despite posterior contraction, potentially resulting
in undesirable properties such as poor calibration of credible intervals."

**Local moves are the mechanism.** Kim and Ročková show that grow-and-prune
cannot reach a deep isolated signal in faster than superpolynomial time, and
mixes well only when the signal is "connected on a tree", meaning every ancestor
of a node that matters also matters.

**Four things fix it, with proofs.**

| remedy | source | effect on mixing time |
|---|---|---|
| more trees | Tan (i) | constant in n |
| global moves on tree structure | Tan (ii), Kim and Ročková's Twiggy | constant in n |
| tempering the tree-structure marginal | Tan (iii) | slower than any power of n |
| informed, thresholded proposals | Kim and Ročková's LIT-MH | linear |
| more chains as n grows | Ronen | not a fix, a mitigation |

Tan et al. add the point that matters most for whether any of this is
acceptable: the first three "do not affect the HDPR set, which means that their
target posterior has essentially the same inferential properties as the original
BART model."

Their empirical numbers. Tan et al., across six datasets, find R-hat rising with
n at the default and the rise "consistently dampened" by more trees or by
temperature 2 or 3; coverage effects "highly ambiguous". Kim and Ročková, on
call-centre data, report a minimum local Gelman-Rubin of 3.58 for plain Bayesian
CART, 1.21 informed, **1.01 Twiggy**, and 9.35 for informed Twiggy, which is a
warning that the two remedies do not compose.

## Why none of this transfers automatically to this package

Three differences, and each has to be measured before anything is built.

**The sampler is not Chipman's.** This package uses Linero's
Laplace-approximation reversible jump: a birth draws the two child leaves from a
Laplace fit rather than integrating them out in closed form. The acceptance
ratio, and therefore the geometry the chain moves through, is not the one any of
the three papers analyzed.

**The rules are soft by default.** All three analyze hard splits on discrete
covariates. A soft gate makes the likelihood a smooth function of the cutpoint,
which is exactly the structure that creates the barriers they describe: a hard
split's contribution is a step, a soft one's is a ramp. The multimodality may
already be much milder here, and the bandwidth move is a further global
rearrangement they have no analogue for.

**The measurements here do not obviously show the problem yet.** SBC passes at
n = 400, chi-square 4.2 on 9 degrees of freedom for rank uniformity, p = 0.90,
with 95% intervals covering at 0.957. The transient is 34 to 70 sweeps, and
burn-in length past that buys nothing, which is consistent with Tan et al.
finding the trend "robust to changes in ... number of burn-in iterations".

So the first question is not which remedy to implement. It is whether this
package has the disease.

## The plan

### Stage 1: does the n-trend exist here (no code changes)

1. **SBC against n.** Re-run `_dev/sbc.R` at n = 250, 500, 1000, 2000, 4000,
   hard rules and soft, and look at whether rank uniformity degrades. This is
   the direct test of Tan et al.'s central claim in this package's own terms,
   and it is the experiment that connects their theory to the under-coverage
   measured on the benches. **Do this first**: if uniformity holds to n = 4000,
   the rest of the programme is speculative here.
2. **R-hat against n**, their Experiment 1 shape: fixed iterations, n from 500
   to 8000, R-hat on held-out fitted values at several quantiles, soft against
   hard. Tests the soft-gate hypothesis directly. Reuse `diagnose()`.
3. **Run-to-run variation**, which is Ronen et al.'s practical claim: fit the
   same data with different seeds and measure the spread of the effect estimate
   against n.

Cost: hours of compute, no new code. Stage 2 is conditional on stage 1.

### Stage 2: the remedies already reachable from the R side

4. **Trees against n.** `num_trees` is already a control. Their claim is that
   more trees flattens the R-hat trend; the counter-consideration is this
   package's own measurement that 200 trees costs 5.8x the time for a third of
   the mixing benefit on a propensity model. Measure the trade in ESS per second
   rather than in R-hat alone, which is the metric the papers do not use.
5. **Chains against n.** Ronen et al.'s recommendation. `chains` is a control
   and the fits are already parallel, so this is a documentation question:
   whether `?bartisan_control` should tell a user to raise `chains` with n.

### Stage 3: new sampler code, cheapest first

6. **Warmup-only tempering.** Tan et al.'s (iii) raises the tree-structure
   marginal likelihood to the power 1/T. In this package's acceptance ratios
   that is the `log_f_after - log_f_before` term in `birth()`, `death()` and
   `change_rule()` in `src/mcmc.cpp`, and multiplying that difference by 1/T is
   a small, local change.

   Doing it **during warmup only, with T returning to 1 for the draws**, keeps
   the sampled target exact and uses tempering purely to find the mode, which
   sidesteps the question of whether a tempered posterior is the one to report
   from. It also composes with the existing `sigma_mu_ramp`, which is already a
   warmup-only schedule, so the machinery for a schedule is there.

   Control: `temperature`, a numeric of length one or two, defaulting to 1.
   Verify with: SBC unchanged (the sampled target must not move), the R-hat
   trend from stage 1, ESS per second, and the audit matrix invariants.

7. **A twig move.** Kim and Ročková's proposal that attaches or detaches an
   entire subtree, and Tan et al.'s (ii). The reversible-jump ratio needs the
   prior probability of the drawn subtree, which this package already computes
   piecewise in `grow_prob()`, so the arithmetic is available. This is real work
   in `src/mcmc.cpp` and `src/node.cpp` and it is where the largest proven gain
   is, constant mixing time in n and a Gelman-Rubin of 1.01 against 3.58 in
   their experiment.

   Do it only if stage 1 shows the disease, and after tempering, because
   tempering is a tenth of the code and buys a bound of its own.

8. **Informed proposals, probably not.** LIT-MH weights each candidate by a
   thresholded posterior ratio, which means evaluating many candidates per move
   where the current sampler evaluates one. The Laplace fit that makes each
   evaluation cheap here is still the dominant cost of a move, so the cost
   multiplies by the number of candidates. Kim and Ročková's own table has
   informed Twiggy at 9.35 against Twiggy's 1.01, so the two do not compose and
   the informed half is not the half doing the work.

## How each will be judged

Nothing here gets adopted on a mixing statistic alone. The bar is:

- **SBC unchanged**, which is the one that matters if the target is touched. A
  tempered sampler that reports from the tempered posterior has to show it.
- **Effective sample size per second**, not R-hat. R-hat is what the papers
  measure and it says nothing about cost; a remedy that doubles mixing and
  triples runtime is a loss, and this package has already found two of those.
- **The audit matrix invariants**, since a new move is a new composition and
  `test-invariants.R` is what catches a hook that stops being called.
- **The estimand**, on `_dev/propensity-settings.R` and `_dev/acic2016.R`. A
  sampler change that moves an average treatment effect is a finding in itself.

## Stage 1 result: the disease is not here

Run 2026-09-05, six cells, and recorded in full in `_dev/TASKS.md`. Rank
uniformity holds at every sample size under both gates, the smallest p over six
cells being 0.37; the pooled trend in mean rank against log n is null under both
(t = 1.30 hard, t = 1.13 soft); and there is no U-shaped histogram anywhere,
which is the signature the theory predicts. At n = 4000 under hard rules, the
case the papers analyze, the fit is the cleanest of the six.

So **stages 2 and 3 below are not started**. They stay written down because the
argument for them is sound and the measurement is bounded: it reaches n = 4000,
and it generates from the prior, which is right for calibration and wrong for
finding a worst case. The soft-gate hypothesis in the section above is also
dead: hard rules are equally clean, so smoothness is not what is saving this
sampler. The remaining candidate is the Laplace reversible jump drawing leaves
at the proposal instead of integrating them out.

## The one thing to do first

Stage 1, item 1. Everything else is conditional on it, it needs no new code, and
it is the only experiment among these whose answer is interesting either way: if
SBC degrades with n, the package has the problem the theory predicts and there
is a proven fix to implement; if it does not, the soft gate is doing something
the theory does not cover and that is worth knowing and writing down.
