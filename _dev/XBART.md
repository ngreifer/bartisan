# XBART, and why none of it is here

XBART (He and Hahn 2023), XBCF (Krantsevich, He and Hahn 2023) and ASBART (Ran
and Bai 2023) have come up repeatedly, from the literature and from readers of
the benchmark table, always in the same form: *this is twenty times faster than
BART, why are you not using it?*

This document is the consolidated answer. It has three parts, and they are
independent, so a reader who rejects one is not obliged to reject the others.

1. **What XBART is**, stated precisely enough that the rest follows.
2. **What was measured here**, which refutes three specific claims in this
   literature and confirms a fourth.
3. **Why it is not implemented**, which rests on structural facts about this
   package rather than on the measurements.

The short version: XBART's speedup is real and is not a measurement artifact,
but it is the speed of computing *a different object*. Three of the claims built
on top of it do not survive measurement here. One derivative paper's headline
number is an iteration-count artifact and nothing more.

**A caution against overreading this document.** It is tempting to summarize
this as "the XBART literature is unfair comparisons all the way down." That
overstates what is here and is not the conclusion. §2.4 says plainly which
claims were tested and which were not, and XBART's own core timing claim is in
the second list.

## 1. What XBART is

Standard BART updates a tree with a reversible-jump Metropolis-Hastings step:
propose a birth, a death or a rule change, accept with a ratio involving the
integrated likelihood. Each proposal is local and most are rejected, so the
structure moves slowly and a run needs thousands of sweeps.

XBART replaces the whole move set with a **recursive grow-from-root sweep**. At
each node it evaluates the marginal likelihood of every candidate cutpoint of
every predictor, plus a null cutpoint that stops the recursion, and samples one
proportional to that likelihood. Then it recurses into both children. A tree is
rebuilt from scratch every sweep rather than perturbed.

Two things make this fast, and both are worth stating because they are the
substance of the speed claim:

- **Presorting.** Each predictor is sorted once. Within a node, the candidate
  cutpoints for a variable are visited in sorted order, so the sufficient
  statistics for every candidate split come from a single cumulative-sum pass,
  and each child's statistics follow from the parent's by subtraction. Evaluating
  all $n_b$ candidates for a variable costs the same as evaluating one.
- **Far fewer sweeps.** XBART is typically run at 40 sweeps with 15 discarded,
  against BART's 1000 warmup plus 1000 draws. This is not an oversight; it
  follows from what the algorithm is (below).

The presorting is a genuine algorithmic win and this document does not dispute
it. The sweep count is where the comparison stops being like-for-like.

### It is not a posterior sampler, and its authors say so

This is the fact that governs everything else. He and Hahn describe the
grow-from-root step as "not a proper full conditional" and the resulting
estimator as "a greedy stochastic approximation." Their only theorem establishes
that *a* stationary distribution exists for a modified version of the algorithm,
not that the stationary distribution is the BART posterior.

Their own tables show the consequence directly: 95% intervals covering as little
as 0.50.

So XBART at 40 sweeps and BART at 2000 are not two routes to the same answer at
different speeds. They compute different things. XBART produces a fast, accurate
*point* estimate with an ensemble of plausible trees around it; BART produces
draws from a posterior. A package whose stated selling point is exactness cannot
substitute the first for the second and describe the result the same way.

That is a statement about scope, not a criticism. XBART is good at what it is
for.

## 2. What was measured here

Four claims from this literature were tested against this package. Three do not
survive. The fourth reproduced and then dissolved on closer measurement, and the
way it dissolved is not the way anyone predicted.

### 2.1 "Burn-in length is the binding constraint" — refuted

Both warm-start items on the To Do list existed to shorten burn-in, on the
premise that warmup was where the time went. It is not.

Fitted with `num_burn = 0` so every sweep is retained, the log likelihood reaches
within two standard deviations of its eventual level at sweep 36 (Gaussian,
soft), 55 (Gaussian, hard) and 61 (`ordinal()`, hard). On the case that should be
worst — 30 predictors with 25 irrelevant — the share of splits landing on the
five that matter plateaus by sweep 34 under soft rules and 48 under hard.

**A warm start can therefore remove at most 34 to 70 sweeps of a 1000-sweep run,
about 5%**, and it cannot reach even that, because it would have to be better
than what the sampler does in those sweeps. The McCartan replication measures a
prior-drawn basis as *worse* on exactly the sparse case where warmup takes
longest, 0.874 against DART's 0.958. The sampler would spend its early sweeps
undoing the warm start.

### 2.2 "The transient grows with $n$, so the regimes differ" — refuted

The obvious objection to §2.1 is that it was measured at package scale, while He
and Hahn's warm start moved coverage from 0.74 to 0.96 at $n = 10000$ and
Krantsevich's at $n = 5000$. There was also a mechanism predicting growth:
`P_BIRTH_DEATH = 0.7` gives 0.35 birth attempts per tree per sweep, birth
acceptance is near a third, so a tree gains about 0.12 internal nodes per sweep,
and the number of internal nodes the data supports grows with $n$.

The trees do grow. Internal nodes per tree go 2.1, 3.9, 7.4 under hard rules at
$n = 1500$, 10000, 50000. **The transient does not follow.**
`_dev/transient-scaling.R` and `_dev/transient-trace.R`, log likelihood as a
fraction of its total rise:

| sweep | 1500 hard | 10000 hard | 50000 hard | 1500 soft | 10000 soft | 50000 soft |
|---|---|---|---|---|---|---|
| 25 | 0.78 | 0.75 | 0.75 | 0.94 | 0.96 | 0.88 |
| 36 | 0.86 | 0.81 | 0.84 | 0.98 | 0.98 | 0.97 |
| 100 | 0.95 | 0.95 | 0.96 | 0.99 | 1.00 | 1.00 |

Reading across a row, $n$ changes nothing across a 33-fold range. The mechanism
is wrong because the log likelihood is dominated by the first few splits of each
tree: going from two internal nodes to seven buys accuracy slowly, and 50 trees
reach their first few splits at the same rate whatever $n$ is.

**So the regime gap closes in the 5% figure's favor.** Whatever moved coverage in
those two papers, burn-in length is not it here.

### 2.3 "Separately initialized chains are what widen the intervals" — reproduced, then dissolved

This is the part of Krantsevich's protocol that survived §2.1 and §2.2. They run
$s - b$ separately initialized chains and pool them, against one chain for BCF,
and their own 20000-after-20000 run still does not match their warm start. So the
hypothesis was that **chain diversity**, which `chains` already provides and which
needs no grow-from-root anywhere, is the active ingredient.

`_dev/chains-vs-length.R`, Friedman at $n = 4000$ with 25 of 30 predictors
irrelevant, 20 replicates, coverage of the true regression function at 500
held-out points against a fixed 4800-sweep budget:

| arm | chains | warmup | draws | pooled | coverage | width |
|---|---|---|---|---|---|---|
| A | 1 | 200 | 4600 | 4600 | 0.949 | 0.538 |
| B | 4 | 200 | 1000 | 4000 | 0.975 | 0.636 |
| D | 16 | 200 | 100 | 1600 | 0.986 | 0.742 |
| E | 16 | 50 | 250 | 4000 | 0.994 | 0.734 |

The effect is large and highly significant: coverage rises 0.026, 0.037, 0.045
with $t$ of 8.0, 9.3, 10.3. It reproduced.

**And then it dissolved.** `_dev/chains-convergence.R` raises total sweeps at a
*constant* 800 stored draws per chain by thinning at 1, 4, 16, 64, so Monte Carlo
error from the draw count is held fixed and a change in width is a change in the
posterior being sampled:

| sweeps past warmup | one-chain-minus-four-chain width gap |
|---|---|
| 800 | +10.8% ($t$ = 3.4) |
| 3200 | +5.9% ($t$ = 3.0) |
| 12800 | +3.2% ($t$ = 2.4) |

Regressing the gap on sweeps gives an exponent of **−0.44**, against the −0.5 that
Monte Carlo error predicts. A real posterior feature would plateau. This decays
at the textbook rate, and the widths converge from both sides to the same place.

**So the apparent multi-chain gain was short runs overstating interval width, the
multi-chain arms more than the single-chain one.** Not a property of chain count.

### 2.4 What this does *not* show, and it matters

Two overreadings to guard against, because both are natural and both are wrong.

**It does not show that chains get stuck in local modes.** The opposite
hypothesis was the one tested. The "real" reading of §2.3 was that a single chain
settles into one variable-selection state and reports the spread within it while
several chains find several. That names a mechanism, so it was tested directly
by turning the Dirichlet prior off, which removes the multimodality in question.
**The hypothesis was refuted, and backwards**: `sparsity = FALSE` makes the gap
two to three times *larger* at every rung, with worse absolute bias (0.132
against 0.108) and worse R-hat. Nothing here found a chain stuck anywhere.

The "initialization determines the answer" finding is real but belongs to a
different model. It is the **Infinite BART** prototype, where the same prior and
the same data gave 161 active trees from a dense start and 79 from a prior start,
with non-overlapping ranges over four chains. That is recorded under Infinite
BART in `_dev/TASKS.md` and says nothing about XBART.

**XBART's own timing claim was never measured here.** The 20–28x has not been
reproduced or disputed. The structural objections in §3 make it moot rather than
false: the mechanism it rests on is real, and does not survive the configuration
this package defaults to. Nothing in this document licenses the claim that
XBART's speedup is fake.

### 2.5 ASBART: the one headline that is an arithmetic artifact

ASBART (Ran and Bai 2023, arXiv 2310.13975) claims "about 10 times faster than
Soft BART with comparable accuracy," which would matter a great deal, since soft
rules are this package's default and its cost center.

It does not hold up, for two reasons that compound.

**The method is XBART with a smoothing pass bolted on.** The paper's own summary
of its contribution is "Combine the SBART and XBART together", and the procedure
is "initially obtaining the tree structures through the XBART method and then
pursuing the optimal window width." So the structure search is hard-rule
grow-from-root, and the bandwidth is fitted afterward on fixed structures. It
inherits both objections in §3.1 and §1 wholesale.

**The 10x is an iteration count.** The baseline is ASBART at **40 iterations**.
Against it, SBART with 2000 warmup and 2000 draws is 276.7x, and SBART at 500
iterations is still 54x. Nearly the entire ratio is 40 draws against 4000. The
paper reports this openly in its own efficiency table; it is the framing that
turns it into a 10x speedup claim.

This is the one place where "artifact and unfair comparison" is exactly the right
description, and it is a derivative paper rather than He and Hahn.

## 3. Why it is not implemented

The measurements above remove the *motivation*. These are the reasons the feature
would not be built even if the motivation returned.

### 3.1 It is a hard-rule technique, and soft rules are the default

The speed comes entirely from presorted index vectors: candidate cutpoints share
one cumulative-sum pass, and each child's statistics follow from the parent's by
subtraction.

**Under soft rules there is no partition to sort.** Every observation reaches
every leaf with a weight that itself depends on the candidate cutpoint, so there
is no prefix to sum and no subtraction that yields a child's statistics. The
entire mechanism dies. Since soft rules are the default and are the package's
accuracy argument — 35–40% lower held-out error, at 3–5x the time — a hard-rules-only
speedup does nothing for the configuration users actually run.

### 3.2 The conjugacy substitution reaches further than expected, and still only serves hard rules

This deserves stating because the first version of this assessment got it wrong
in the conservative direction.

`Target1` in `src/mcmc.cpp` determines the whole log target over a leaf from
three sums at a reference point: value, score and information. All three are
additive over observations, so all three prefix-sum exactly as a residual sum
does. Three sums instead of two is the only difference. Where the target is
`TARGET_EXP_UP` or `TARGET_EXP_DOWN` the mode then comes from
`exponential_mode()`, a scalar Newton loop touching no data, so iterating it per
candidate costs scalar arithmetic rather than a factor of the grid size.

That splits the family list three ways:

- **Exact and closed form** (target quadratic): `gaussian()`, `binomial()` under
  both links with `augment = TRUE`, multinomial logit and probit, ordinal probit
  and logit, the zero-inflated families, lognormal and log-logistic AFT, `dpm()`.
- **Exact up to the same Laplace the leaf updates already take**: `poisson()`,
  `Gamma()`, negative binomial, ordinal cloglog, `ph()`, and the log-scale
  predictor of `gaussian_ls()`.
- **One-step only**: `Beta()`, `ordbeta()`, `tweedie()`, `custom_family()`,
  `Gamma_ls()`'s shape, and every family under `augment = FALSE`.

So the criterion is further out than "conjugate families only." It still serves
only the hard-rule path, which returns to §3.1.

### 3.3 There is no cutpoint grid to normalize over

`Node::draw_rule()` in `src/node.cpp` ends
`val = lower + (upper - lower) * unif_rand()`: the cutpoint prior is **continuous**
on the node's live range.

Grow-from-root needs a finite candidate set, and its null-cutpoint weight
$|\mathcal{C}|((1 + d)^\beta / \alpha - 1)$ is calibrated to that set's size,
which is how it reproduces BART's $\alpha (1 + d)^{-\beta}$ branching
probability. Introducing a grid means grow-from-root draws from a **different tree
prior** than the sampler targets, so a warm start would place the chain at a draw
from the wrong model.

The default `x_transform = "smoothcdf"` makes the mismatch mild, since uniform
cutpoints in a smoothed-CDF coordinate are already close to quantile-spaced in
$x$. Mild is not none, and it would have to be measured rather than assumed.

### 3.4 `encode_tree()` is one-way

There is no decoder in `src/model.cpp`, so a warm start cannot go through the
saved `forest_flat`. It has to build `Node`s in process and hand them to
`update_forest()`. That direction is easy given `birth_leaves()`,
`split_support()` and the node pool, but it puts the whole feature in C++ with no
R prototype, which is the wrong shape for something speculative.

### 3.5 XBCF adds a protocol, not an algorithm

Krantsevich, He and Hahn's grow-from-root is He and Hahn's applied tree by tree
to two forests. Their parameter updates (the $a$, $b_0$, $b_1$ rescaling and
separate $\sigma_0$, $\sigma_1$) are not this package's parameterization, since
`bcf()` puts the effect in a coefficient forest inside the H-component framework.

What is genuinely new is that the $s - b$ post-burn-in XBCF forests initialize
$s - b$ *independent* MCMC chains run in parallel, so the warm start buys chain
diversity rather than one starting point. Parallel chains are already here — and
§2.3 measures that ingredient and finds it is Monte Carlo error.

## 4. The one variant that was never assessed away

Both options above either change the starting point (warm start) or replace the
kernel (grow-from-root as the sampler). The coverage work finds that what binds
is **effective sample size per sweep**, which neither touches.

A third option does: **use the grow-from-root criterion as an informed proposal
for the existing birth move** rather than as a replacement for it. Draw the
cutpoint proportional to $L(c)$ over the candidate set and correct with the
Hastings ratio, whose normalizing constant $\sum_c L(c)$ falls out of the same
pass. That stays an exact sampler, raises birth acceptance from about a third
toward one, and puts the rule where the likelihood wants it rather than where the
prior put it.

Costs and caveats: $O(n_b p)$ per birth against $O(n_b)$ now, so it wins only if
acceptance and placement beat a factor of $p$; it needs the same presorting; and
it is still hard-rules-only, so it does nothing for the default. **Not measured.**

The recent literature has caught up with this idea and is worth reading before
anyone builds it. Locally-balanced Markov processes (Livingstone et al. 2025)
give the general theory for proposals on discrete spaces; the Taxicab Sampler
(Geels et al. 2021) is built explicitly with tree models as the application;
Multiple Jump MCMC (Vogels et al. 2026) and similarity-driven proposals (Aiello
et al. 2026) are the same idea from two other directions.

One result there is a direct warning about this variant. Sasli et al. (2026),
*Learned proposals in trans-dimensional inference are optimal at equilibrium, not
during assembly*, show that the optimal birth proposal differs between the two
regimes: **residual-matching during assembly, marginal-matching at equilibrium.**
In BART backfitting every tree is fitted to a residual, so the assembly regime is
not a transient the chain passes through — it is the normal state of every tree
update. A proposal tuned for equilibrium may therefore hurt exactly where this
one is meant to help. That is the thing to measure first.

## 5. Summary of the verdict

| Claim | Status here |
|---|---|
| XBART's presort-and-prefix-sum speedup is real | **Granted.** Not disputed, not measured. |
| XBART draws from the BART posterior | **False**, per its own authors. |
| Burn-in length binds, so a warm start pays | **Refuted.** Transient is 34–70 sweeps of 1000. |
| The transient grows with $n$ | **Refuted.** Flat across a 33-fold range. |
| Pooling separately initialized chains widens intervals | **Reproduced, then shown to be Monte Carlo error** (exponent −0.44). |
| A single chain gets stuck in one sparsity mode | **Refuted, backwards.** `sparsity = FALSE` makes it worse. |
| ASBART is 10x faster than SBART | **Artifact.** 40 iterations against 4000. |

**What would change the decision.** A soft-rule analogue of the presorting trick,
which nobody has proposed and which looks impossible rather than merely hard; or
a measurement showing the informed-birth variant of §4 beats a factor of $p$ on
effective sample size per second under hard rules, with Sasli et al.'s
assembly-regime warning accounted for.

## References

- He, J. and Hahn, P. R. (2023). Stochastic tree ensembles for regularized
  nonlinear regression. *JASA* 118(541), 551–570.
- Krantsevich, N., He, J. and Hahn, P. R. (2023). Stochastic tree ensembles for
  estimating heterogeneous effects. *AISTATS*.
- Ran, H. and Bai, Y. (2023). ASBART: Accelerated Soft Bayes Additive Regression
  Trees. arXiv:2310.13975.
- Geels, V., Pratola, M. T. and Herbei, R. (2021). The Taxicab Sampler: MCMC for
  Discrete Spaces with Application to Tree Models. arXiv:2107.07313.
- Livingstone, S. et al. (2025). Foundations of locally-balanced Markov
  processes. arXiv:2504.13322.
- Sasli, A. et al. (2026). Learned proposals in trans-dimensional inference are
  optimal at equilibrium, not during assembly. arXiv:2608.12392.
- Vogels, L. et al. (2026). Multiple Jump MCMC. arXiv:2603.22573.

Measurements live in `_dev/transient-scaling.R`, `_dev/transient-trace.R`,
`_dev/chains-vs-length.R` and `_dev/chains-convergence.R`. The narrative record
is under "Assessed and not adopted", "Assessment: the two warm-start items" and
"The chain-count effect is Monte Carlo error, measured" in `_dev/TASKS.md`.
