# Sampler and prior settings for `bartisan()`

Collects the tuning constants of the sampler and the hyperparameters of
the tree prior. The defaults follow Linero (2025) and, for the soft
decision rules, Linero and Yang (2018), and are intended to be usable
without adjustment.

The arguments fall into three groups and are ordered by how often they
are worth changing. Everything from `chains` to `x_transform` is a
**modeling decision**, in that changing one changes what is fitted or
how long it is fitted for; the chain lengths and `num_trees` come first
because they are adjusted most. Everything from `gamma` to `num_print`
is an **advanced setting**: a hyperparameter of a prior the first group
summarizes, or a switch whose default is almost always right. The last
three (`block_eval`, `exact_quadratic`, and `generic_accumulate`) exist
**for internal validation**; they compute the same posterior more slowly
and are documented so that the checks using them can be read.

## Usage

``` r
bartisan_control(
  chains = 1L,
  num_burn = 200L,
  num_draws = 800L,
  num_thin = 1L,
  num_trees = NULL,
  gate = "smoothstep",
  sparsity = TRUE,
  k = 2,
  bandwidth = 0.1,
  split_prior = NULL,
  share_sparsity = FALSE,
  categorical = "subset",
  augment = TRUE,
  x_transform = "smoothcdf",
  gamma = 0.95,
  beta = 2,
  sigma_mu = NULL,
  update_sigma_mu = TRUE,
  sigma_mu_ramp = 0.25,
  update_tau = TRUE,
  update_bandwidth = TRUE,
  bandwidth_every = 1L,
  alpha = NULL,
  alpha_scale = NULL,
  alpha_shape_1 = NULL,
  alpha_shape_2 = NULL,
  update_s = NULL,
  update_alpha = NULL,
  verbose = FALSE,
  num_print = 100L,
  block_eval = FALSE,
  exact_quadratic = TRUE,
  generic_accumulate = FALSE
)
```

## Arguments

- chains:

  `numeric`; how many independent chains to run. Default is 1. The draws
  are pooled, and
  [`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
  computes split-R-hat from them on request rather than the fit carrying
  it. With the future.apply package installed the chains run in parallel
  under whatever backend
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  has set. One [`set.seed()`](https://rdrr.io/r/base/Random.html) before
  the call reproduces the run either way, since each chain gets its own
  L'Ecuyer stream. See
  [`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md).

- num_burn:

  `numeric`; the number of warmup iterations to discard. Default is 200,
  which measurement finds ample for the samplers here. Warmup is where
  the trees grow into the data and the hyperparameters find their scale,
  so raising it buys convergence rather than precision: increase it when
  [`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)'s
  `rhat` says the chains have not agreed, and for a
  [`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  fit, whose component-count state settles more slowly than a forest
  does. See Details.

- num_draws:

  `numeric`; the number of draws to keep. Default is 800. These are what
  every estimate and interval is computed from, so raising it narrows
  Monte Carlo error and does nothing about convergence: increase it when
  `ess_bulk` or `ess_tail` is small relative to what the reported
  quantity needs.

- num_thin:

  `numeric`; keep one draw in every `num_thin` after warmup. Default is
  1 to keep every draw. Thinning discards draws to make the kept ones
  less correlated, which costs information and is worth it only to hold
  down the memory a long chain would otherwise take: for a given amount
  of computing, more draws beat fewer less-correlated ones.

- num_trees:

  `numeric`; the number of trees, as one number for every forest or one
  per additive predictor. Default is `NULL`, which is 50 whatever the
  decision rules. A family with more than one additive predictor takes a
  vector, which is worth using because the forests are neither equally
  expensive nor in need of equal capacity. See Details.

- gate:

  string; the shape of a decision rule, which is also how hard and soft
  rules are chosen between. Allowable options include `"smoothstep"`
  (the default), `"smootherstep"`, `"logistic"`, and `"hard"` (or
  equivalently `"step"`). `"hard"` gives the step functions of standard
  BART; the other three give soft rules, in which every observation
  reaches every leaf with some weight and the fit is smooth. Soft rules
  cost several times as much per iteration and cut held-out error
  substantially, so they are the accuracy argument rather than a tax;
  which soft gate is chosen matters much less than that one is. See
  Details.

- sparsity:

  `logical` or string; the prior on which predictors are split on.
  `TRUE` (the default) is the Dirichlet sparsity prior of Linero (2018)
  (i.e., DART), which can drop a predictor from the forest entirely, and
  `FALSE` gives every predictor the same splitting probability, which is
  classic BART. The strings `"none"`, `"weak"`, `"moderate"`, and
  `"strong"` name four strengths, `"none"` equal to `FALSE` and
  `"moderate"` to `TRUE`. Note that this sets `update_s`,
  `update_alpha`, `alpha_shape_1` and `alpha_shape_2` together, and that
  supplying any of those overrides it. Read the trade-off in Details
  before turning it off.

- k:

  `numeric`; controls the leaf prior. The prior standard deviation of a
  forest is `3 / k` times the natural scale of its additive predictor,
  so larger values shrink the fit harder toward the intercept-only
  model. Default is 2.

- bandwidth:

  `numeric`; the prior mean of the gate bandwidth of a soft rule, on the
  scale of the transformed predictors, which lie in `[0, 1]`. Smaller
  values approach hard rules. Default is .1. Ignored when
  `gate = "hard"`.

- split_prior:

  `numeric`; a named vector of relative prior weights on the predictors,
  keyed by their names in the formula. A predictor not named gets a
  weight of 1, and the probability of splitting on one is its weight
  over the total. Weights must be finite and nonnegative; zero keeps the
  predictor out of every tree, and naming a predictor the model does not
  have is an error. Default is `NULL` to weight every predictor equally.
  Note that setting this overrides `sparsity`; see Details.

- share_sparsity:

  `logical`; for a family with more than one additive predictor, whether
  the forests draw their splitting proportions from one pooled Dirichlet
  rather than one each. Default is `FALSE`. `TRUE` is an assumption
  about the data rather than a free improvement, and requires that the
  proportions be drawn at all and over the same predictors, so
  `sparsity` must not be `FALSE` for the forests that are to share.
  Ignored for a family with a single forest. See Details.

- categorical:

  string; how a splitting rule divides the levels of a factor. Allowable
  options include `"subset"` (the default), which draws a subset of the
  levels still available at the node and sends those left, and
  `"onehot"`, which is what most BART implementations do: it splits on
  one indicator column, peeling a single level off the rest. `"subset"`
  is right about the prior and never loses beyond noise; see Details.

- augment:

  *Advanced.* `logical` or `character`; whether to rewrite the
  likelihood as the margin of a Gaussian or a Poisson one, which makes
  the target a shape the sampler can exploit. **The posterior is the
  same either way**, so this is a sampling setting and not a modeling
  one. Default is `TRUE`, which rewrites wherever a rewriting exists:
  the binomial, ordinal, multinomial, negative binomial, zero-inflated,
  and survival families. `FALSE` never does, and a character vector of
  engine family names (`"binomial"`, `"ordinal"`, `"multinomial"`,
  `"negbin"`, `"zip"`, `"zinb"`, `"aft"`) asks for exactly those. A
  rewriting is always faster and has not been measured to mix worse per
  second, so the default is rarely worth changing; see
  [`vignette("implementation")`](https://ngreifer.github.io/bartisan/articles/implementation.md).

- x_transform:

  string; how numeric predictors are mapped to `[0, 1]`. Allowable
  options include `"smoothcdf"` (the default), `"quantile"`, and
  `"range"`. The choice places the cutpoint prior and, for soft rules,
  sets the scale the gate's width is measured on, and each has a failure
  it is worth knowing about; see Details.

  `"smoothcdf"` uses a kernel-smoothed estimate of each predictor's
  distribution function, with an Epanechnikov kernel and a plug-in
  bandwidth at the `n^(-1/3)` rate that is right for a distribution
  function rather than the `n^(-1/5)` that is right for a density
  (Tenreiro, 2006; the rate was pointed out for a second-order
  approximation by Azzalini, 1981). Cutpoints land where the data are,
  as under `"quantile"`, and the map is strictly increasing and
  differentiable, as under `"range"`.

  `"quantile"` uses the empirical distribution function, which is what
  `SoftBart::softbart()` does. It is a step function, so the fit is a
  step function of the predictor: it has no derivative, and a
  relationship that is straight in the predictor becomes a jump wherever
  the data are gappy.

  `"range"` rescales linearly and preserves the original spacing. It is
  the most accurate of the three on well-behaved predictors and the
  right choice when a derivative is the quantity of interest, and it
  fails where a few extreme values leave the bulk of a predictor inside
  a sliver of its range.

- gamma, beta:

  *Advanced.* `numeric`; the branching probability at depth `d` is
  `gamma * (1 + d)^(-beta)`. Defaults are .95 and 2.

- sigma_mu:

  *Advanced.* `numeric`; the prior median of the leaf standard
  deviation, one value per additive predictor. Default is `NULL` to
  derive it from `k` and that forest's own tree count.

- update_sigma_mu:

  *Advanced.* `logical`; whether to draw the leaf standard deviation
  under a half-Cauchy prior rather than fixing it. Default is `TRUE`.
  `FALSE` is worth reaching for when a binary fit mixes badly: where the
  predictors separate the response well the leaf scale is barely
  identified, and it wanders rather than settling, which drags the
  effective sample size of everything built on it down with it. Fixing
  it can raise the effective sample size of the additive predictor and
  the log likelihood several times over at little cost in fit; see
  [`vignette("implementation")`](https://ngreifer.github.io/bartisan/articles/implementation.md).

- sigma_mu_ramp:

  *Advanced.* `numeric`; the fraction of warmup over which the leaf
  standard deviation is raised from near zero to its target. Default is
  .25; set to 0 to disable. Linero (2025) describes this as essential:
  started at its full value, the sampler can settle early into a poor
  configuration and fail to move.

- update_tau:

  *Advanced.* `logical`; whether to draw the standard deviation of each
  random-effect term under the same half-Cauchy prior the leaf scale
  uses, rather than fixing it at that prior's median. Default is `TRUE`.
  Only relevant when the formula has a `(1 | group)` term.

- update_bandwidth:

  *Advanced.* `logical`; whether to draw the bandwidth of each tree
  rather than holding it at `bandwidth`. Default is `TRUE`.

- bandwidth_every:

  *Advanced.* `numeric`; how many sweeps between bandwidth draws for a
  given tree. Default is 1. The bandwidth is one scalar per tree, drawn
  by an adaptive random walk, and every attempt costs a full rebuild of
  the tree's memberships, which is the largest single item in a
  soft-rule fit, so drawing it less often trades mixing in that one
  parameter for speed. Raising it is reasonable when the mean function
  is known to be smooth and the fit is compute-bound; see Details.

- alpha:

  *Advanced.* `numeric`; the concentration of the Dirichlet prior on the
  splitting proportions, where smaller values concentrate splits on
  fewer predictors. Default is `NULL`, which uses 1 as the starting
  value for a parameter that is then drawn.

- alpha_scale, alpha_shape_1, alpha_shape_2:

  *Advanced.* `numeric`; the prior on `alpha`, in which
  `alpha / (alpha + alpha_scale)` is Beta(`alpha_shape_1`,
  `alpha_shape_2`). Default for `alpha_scale` is `NULL` to use the
  number of predictor groups; the two shapes default to whatever
  `sparsity` implies.

- update_s, update_alpha:

  *Advanced.* `logical`; whether to draw the splitting proportions and
  their concentration. Defaults are `NULL` to follow `sparsity`. Turning
  both off recovers a uniform prior over predictors, which is what
  `sparsity = FALSE` does.

- verbose:

  *Advanced.* `logical`; whether to print progress to the console while
  sampling. Default is `FALSE`. For a progress bar instead, see the
  Progress section below, which needs no argument here.

- num_print:

  *Advanced.* `numeric`; how many iterations between the reports
  `verbose` prints. Default is 100.

- block_eval:

  *Validation.* `logical`; whether to evaluate the likelihood one leaf
  at a time rather than one observation at a time. Default is `FALSE`. A
  family built by
  [`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  does this regardless, since it must call back into R; setting it for a
  compiled family produces the same draws from the same seed at somewhat
  greater cost.

- exact_quadratic:

  *Validation.* `logical`; whether to use the closed forms that a target
  quadratic in the additive predictor allows, in which one pass over a
  node determines the log target everywhere, so that the Laplace
  approximation is the conditional posterior rather than an
  approximation to it. Default is `TRUE`, which is what makes a Gaussian
  response, or any of the rewritings in `augment`, cheap. `FALSE` falls
  back on the general path; the two agree, at greater cost.

- generic_accumulate:

  *Validation.* `logical`; whether to accumulate a leaf's sums through
  the family's virtual interface rather than through its own statically
  dispatched loop. Default is `FALSE`. The two compute the same thing;
  the second lets the compiler inline the family's arithmetic, which is
  most of the remaining per-observation cost.

## Value

A `<bartisan_control>` object, a list containing the supplied settings
and the defaults for those not supplied, for passing to the `control`
argument of
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).

## Details

### How Long Warmup Needs to Be

The default is comfortably longer than the transient in every family
measured, so raising `num_burn` buys convergence rather than precision:
raise it when
[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
says the chains have not agreed.
[`dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is the exception, its mixture carrying a component-count state that
settles more slowly than a forest does, so a longer warmup is worth
having there and especially when the error distribution itself is the
object of interest.

### How Many Trees, and How Many per Forest

**Soft rules need far fewer trees than hard ones**, and can get worse as
trees are added, which makes the count most BART packages default to
actively wrong here. The default is the same for both anyway, because a
smaller forest mixes worse: chains can disagree markedly on an average
contrast at a small forest where they agree at the default, for only a
slight gain in point accuracy. Drop below it if prediction is the only
goal and the fit is soft; raise it with hard rules if it is not.

**A vector is worth using when the family has more than one forest.**
The scale forest of
[`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
dominates the run time, its target not being quadratic, and a variance
surface carries much less information than a mean surface, so giving it
fewer trees runs substantially faster at the same accuracy. It is not
the default because how much resolution a variance surface needs depends
on the surface, and quietly under-parameterizing it would show up as
intervals that are wrong, which is what
[`gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
exists to get right. The leaf prior scale divides by the square root of
each forest's own tree count, so shrinking one forest does not change
the prior on the sum.

### The Sparsity Prior, and What It Costs

`sparsity = TRUE` is the Dirichlet prior of Linero (2018) on the
splitting proportions, and it is a genuine variable-selection prior: it
can and does drop a predictor from every tree at once, which is the
point of it in the high-dimensional problems it was built for.

That has a consequence which is easy to misread. **A contrast on a
predictor the prior has dropped is exactly zero, not nearly zero**,
because in that draw the fit does not depend on that predictor at all,
so the posterior of a contrast has an atom at zero whose mass is one
minus the predictor's inclusion probability; any summary reporting a
median lands on it once it holds half the mass. More trees does not fix
it, because a predictor whose splitting proportion has gone small is
rarely proposed and so is hard to get back in.

Which way to throw the switch follows from the estimand, and the
strength barely matters either way. **For prediction or variable
selection**, keep the default; reaching past `TRUE` is warranted only
when the predictors are many and nearly all are expected to be
irrelevant. **For a contrast, a partial effect, or a treatment effect**,
set `sparsity = FALSE`, or use `split_prior`, which cannot drop
anything: on a weak signal the prior attenuates the estimate
substantially and its interval covers below its nominal rate, while a
strong effect is untouched.
[`vignette("effects")`](https://ngreifer.github.io/bartisan/articles/effects.md)
works this through and
[bartisan-marginaleffects](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
covers the atom.

**For a varying-coefficient model the two answers can differ by
forest.** What the prior can drop is a predictor a forest splits on, and
in a [`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) model
the treatment is the coefficient rather than one of those, carried by a
forest of its own, so no splitting proportion can drop it; the prior on
that forest selects among the moderators instead. So
`sparsity = c(FALSE, TRUE)` is coherent, and in
[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) it is
the asymmetry worth considering.

### Sharing the Sparsity Prior Across Forests

The trade is asymmetric, which is what makes `share_sparsity = TRUE`
worth taking: nothing to gain when the predictors are few and there is
no selection problem to transfer, a real gain for the weaker component
when they are many, and a comparatively small cost when the components
turn out to depend on different predictors.
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
shows it in use and
[`variable_importance()`](https://ngreifer.github.io/bartisan/reference/variable_importance.md)
is where to check it. This is the variable-selection content of the
shared forests of Linero et al. (2020) and not their model, which shares
the tree *topology* and so fixes the cut points too;
[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) sits at
the other extreme, sharing an entire function rather than a prior.

### Arguments That Vary by Forest

A family with several additive predictors has one forest per predictor,
each with its own prior. Every argument that could mean something
different for one of them may be given once, to apply to all, or once
per forest, either positionally or keyed by the forest names listed in
[bartisan-families](https://ngreifer.github.io/bartisan/reference/bartisan-families.md).
A forest a named argument does not mention keeps that argument's default
rather than borrowing another forest's value. That covers `num_trees`,
`k`, `sigma_mu`, `sparsity`, `split_prior`, `bandwidth`, `gamma`,
`beta`, the four `alpha` arguments, and the three `update_` flags;
`formula` works the same way, as
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
describes.

    bartisan_control(num_trees = c(mean = 50, log_sd = 10),
                     sparsity = c(mean = TRUE, log_sd = FALSE))

The multinomial families are the exception, for the reason given in
[bartisan-families](https://ngreifer.github.io/bartisan/reference/bartisan-families.md):
their forests act as one, so these arguments take a single value.

### The Predictor Transform (`x_transform`)

Numeric predictors are mapped to `[0, 1]` before any rule sees them, and
the map does two jobs at once: cutpoints are uniform on a node's live
range in that coordinate, and for soft rules the gate's bandwidth is
measured there too. The three options fail in different places, which is
the whole of the choice between them.

`"range"` is the most accurate on well-behaved predictors and the only
one that supports a derivative, since it is the only affine map of the
three. Its failure is a predictor whose bulk sits inside a sliver of its
range, which a few extreme values are enough to produce: a cutpoint
drawn uniformly on a node's live range then almost never lands where the
structure is. Measured over a simulation with 1% of a predictor at 300
times the scale and a truth that oscillates within the bulk, it was 4
times worse than the alternatives at `n = 500` and 6 times worse at
`n = 2000`, with 95% intervals covering .24 rather than .95. The
deterioration with `n` is the signature: the extremes grow with the
sample, so the bulk occupies an ever smaller share.

`"quantile"` cannot fail that way, since it places cutpoints by rank. It
fails instead by being a step function, so the fit is a step function of
the predictor: there is no derivative to take, and a relationship that
is straight in the predictor becomes a jump wherever the data are gappy.
On two tight clusters with a linear truth it was nearly three times
worse than the alternatives.

`"smoothcdf"`, the default, smooths the same distribution function and
so does neither. Across eight data-generating processes it was never the
most accurate and never far from it, where each of the other two was
badly wrong somewhere. It is not a free lunch: a derivative taken
through it is the derivative of the fitted function times an estimated
density, so `"range"` remains the better choice when a slope is the
quantity of interest rather than a prediction.

The bandwidth is a plug-in rather than one chosen from the data. The
data-based selectors of Bergmann and Zaehle (2026) were measured here
and are not used, because they optimize a different thing: they minimize
the error in the estimated distribution function, where this map exists
to place cutpoints and keep the coordinate smooth. The two come apart
exactly where it matters. Given two tight clusters with a gap between
them, they correctly choose a much narrower bandwidth, the distribution
function really being flat in the gap; that returns the coordinate to
nearly the step `"quantile"` would have given, and doubled the error of
the fit.

### Splitting a Factor

What `"onehot"` costs is partial pooling. A rule on one indicator column
can only peel a single level off the rest, so the partitions it reaches
are a small fraction of those available and none divides the bulk of the
levels; a typical tree leaves one level alone and the rest together,
whether the data want that or not. A rule on a subset of the levels
reaches every partition.

Under hard rules `"subset"` is better at every sample size, clearly so
when levels are thinly observed; under soft rules, which is the default,
**the two are indistinguishable**. So `"subset"` is the default because
it is right about the prior and never loses beyond noise, not because it
will visibly improve a fit. Note that a rule on a factor is always hard
even in a soft tree, since a gate is a function of the distance from a
cutpoint and there is no distance between two levels; and that a
two-level factor is unaffected either way.

### Telling the Prior What Is Already Known

`sparsity` and `split_prior` answer different questions and cannot both
be in force, so giving `split_prior` turns `sparsity` off. `sparsity` is
for when which predictors matter is unknown and the prior is to work it
out from the data, and a predictor can be dropped entirely;
`split_prior` is for when something is known and is to be honored, with
the proportions held at the supplied weights. Neither is about a forest
being held to the predictors its own formula names: a
[`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) term whose
moderators are only some of the covariates still draws its splitting
proportions, over those moderators.

A weight is a statement about relative attention, not about effect size.
It changes how often a split on a predictor is proposed, which is a
prior, so the data can still overrule it in either direction. Because
the weights are fixed, `split_prior` does not accumulate the
atom-at-zero mass above, which makes it a reasonable middle course when
a particular contrast is the estimand but the predictors are too many to
treat alike. One weight per term in the formula, not per column of the
design matrix, so a factor is named once and its levels share the
weight.

### What `augment` Does Not Cover

[`Gamma()`](https://rdrr.io/r/stats/family.html) and
[`poisson()`](https://rdrr.io/r/stats/family.html) need no rewriting,
their targets being already in the exponential form the sampler
exploits, and none is known for the accelerated failure time, ordered
beta, or location-scale families, so `augment` is silent for all of
them. Which families it does cover, what each rewriting is, and what
each buys are in
[`vignette("implementation")`](https://ngreifer.github.io/bartisan/articles/implementation.md).

### Rare Events Mix Slowly

With very few events the whole fit mixes slowly, augmented or not, and
the effective sample size of the level of the predictor can fall to a
small fraction of what a balanced response gives. That is the
information in a handful of events rather than a fault of the
rewritings, and dbarts, an independent implementation of the same latent
normal, reproduces it. Lengthen the chain, and read
[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)'s
`rhat`.

### Progress

`verbose = TRUE` prints a line every `num_print` iterations, which is
the whole of what this package decides about progress. A progress bar is
[progressr](https://CRAN.R-project.org/package=progressr)'s business,
and the sampler reports to it unconditionally: nothing is shown unless a
handler is active, so there is no argument to switch on.

    progressr::with_progress(
      bartisan(y ~ ., data = d, family = gaussian())
    )

    # or once, for the session
    progressr::handlers(global = TRUE)

The bar is sized for the whole fit, so `chains = 4` fills one bar once
rather than four in sequence, and chains running in parallel under
[future](https://CRAN.R-project.org/package=future) relay their progress
back as it arrives. Convergence diagnostics are not included because
they are not part of a fit; they run in
[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md),
which reports its own progress the same way. Progress does not touch the
draws.

### What Else Runs in Parallel, and One Limit on It

A `future` plan is used by three things, and only the first is the
sampler: the chains of a fit, the per-observation pass in
[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md),
and the repeated predictions in
[`partial_dependence()`](https://ngreifer.github.io/bartisan/reference/partial_dependence.md)
and
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md).
The sampler has no other axis, since a sweep conditions on the one
before it. Measured on 1500 observations with a 25-point grid, the
predictions run about three times faster on four workers.

The limit worth knowing is future's, not this package's. Anything that
predicts on a worker has to be sent the fit, and a fit is mostly its
stored predictor and flattened forests, which grow with the draws and
the sample: at 1500 observations and 3200 draws one serializes to about
70 MB, so four workers move 280 MB.
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
refuses a single export above `future.globals.maxSize`, 500 MB by
default, and a fit large enough on both counts will trip it. The error
names the option; raising it is the fix, and running sequentially is the
alternative.

### Soft Rules and the Cost of a Gate

A soft rule is charged for in two places: every observation reaches more
than one leaf, so a pass over a node covers several times the sample,
and the bandwidth is itself a parameter with a Metropolis step per tree
per sweep, each rebuilding every membership weight in the tree, which is
the single largest item in a soft-rule fit.

A bounded gate addresses the first: past its half-width from the
cutpoint the gate is exactly zero or one, so the observation takes one
side outright and the gate is a polynomial rather than an
[`exp()`](https://rdrr.io/r/base/Log.html). Most of what that saves is
the [`exp()`](https://rdrr.io/r/base/Log.html) rather than the work past
the cutpoint, which is why *which* bounded gate is chosen makes almost
no difference; prefer `"smootherstep"` because it gives a
twice-differentiable fit. Its half-width is `pi * sqrt((2a + 1) / 3)`
times `bandwidth` for the Beta(a, a) gate, which equates the gates'
standard deviations so that `bandwidth` means the same amount of
smoothing whichever is chosen.

`bandwidth_every` addresses the second. Drawing the bandwidth is what
lets a rule sharpen toward a step, so fixing it
(`update_bandwidth = FALSE`) is faster still and *more* accurate on
smooth functions and much worse on nonsmooth ones. Raising
`bandwidth_every` is the middle course, recovering some speed while
keeping soft rules, at a real cost in mixing and a small one in accuracy
where the mean function jumps.

## References

Azzalini, A. (1981). A note on the estimation of a distribution function
and quantiles by a kernel method. *Biometrika*, 68(1), 326–328.
[doi:10.1093/biomet/68.1.326](https://doi.org/10.1093/biomet/68.1.326)

Bergmann, T., & Zaehle, H. (2026). Data-based bandwidth selection for
kernel smoothing of empirical distribution functions. *Metrika*.
[doi:10.1007/s00184-026-01025-6](https://doi.org/10.1007/s00184-026-01025-6)

Tenreiro, C. (2006). Asymptotic behaviour of multistage plug-in
bandwidth selections for kernel distribution function estimators.
*Journal of Nonparametric Statistics*, 18(1), 101–116.
[doi:10.1080/10485250600578334](https://doi.org/10.1080/10485250600578334)

Linero, A. R. (2018). Bayesian regression trees for high-dimensional
prediction and variable selection. *Journal of the American Statistical
Association*, 113(522), 626–636.
[doi:10.1080/01621459.2016.1264957](https://doi.org/10.1080/01621459.2016.1264957)

Linero, A. R. (2025). Generalized Bayesian additive regression trees
models: beyond conditional conjugacy. *Journal of the American
Statistical Association*, 120(549), 356–369.

Linero, A. R., & Yang, Y. (2018). Bayesian regression tree ensembles
that adapt to smoothness and sparsity. *Journal of the Royal Statistical
Society Series B*, 80(5), 1087–1110.

Linero, A. R., Sinha, D., & Lipsitz, S. R. (2020). Semiparametric
mixed-scale models using shared Bayesian forests. *Biometrics*, 76(1),
131–144. [doi:10.1111/biom.13107](https://doi.org/10.1111/biom.13107)

## See also

[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md),
which takes the result as its `control` argument;
[bartisan-families](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
for the forest names the per-forest arguments are keyed by

## Examples

``` r
data("rhc")

# Settings can be built up once and reused across fits
ctrl <- bartisan_control(num_trees = 20, gate = "hard", num_burn = 50,
                         num_draws = 50)

fit <- bartisan(death ~ . - days, data = rhc, control = ctrl)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` to choose another, which also silences this message.

# The same call, with the settings passed through `...` instead
fit2 <- bartisan(death ~ . - days, data = rhc, num_trees = 20,
                 gate = "hard", num_burn = 50, num_draws = 50)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` to choose another, which also silences this message.

# A setting given once applies to every forest, and a vector gives each
# forest its own value. A variance surface needs less capacity than a
# mean surface
bartisan_control(num_trees = c(mean = 50, log_sd = 10))
#> $num_trees
#>   mean log_sd 
#>     50     10 
#> 
#> $gate
#> [1] "smoothstep"
#> 
#> $soft
#> [1] TRUE
#> 
#> $sparsity
#> [1] TRUE
#> 
#> $share_sparsity
#> [1] FALSE
#> 
#> $share_forests
#> [1] FALSE
#> 
#> $split_prior
#> NULL
#> 
#> $categorical
#> [1] "subset"
#> 
#> $k
#> [1] 2
#> 
#> $bandwidth
#> [1] 0.1
#> 
#> $chains
#> [1] 1
#> 
#> $num_burn
#> [1] 200
#> 
#> $num_draws
#> [1] 800
#> 
#> $num_thin
#> [1] 1
#> 
#> $augment
#> [1] "binomial"    "ordinal"     "multinomial" "zip"         "zinb"       
#> [6] "aft"         "negbin"     
#> 
#> $x_transform
#> [1] "smoothcdf"
#> 
#> $gamma
#> [1] 0.95
#> 
#> $beta
#> [1] 2
#> 
#> $sigma_mu
#> NULL
#> 
#> $update_sigma_mu
#> [1] TRUE
#> 
#> $sigma_mu_ramp
#> [1] 0.25
#> 
#> $update_tau
#> [1] TRUE
#> 
#> $update_bandwidth
#> [1] TRUE
#> 
#> $bandwidth_every
#> [1] 1
#> 
#> $alpha
#> [1] 1
#> 
#> $alpha_scale
#> [1] 0
#> 
#> $alpha_shape_1
#> [1] 0.5
#> 
#> $alpha_shape_2
#> [1] 1
#> 
#> $update_s
#> [1] TRUE
#> 
#> $update_alpha
#> [1] TRUE
#> 
#> $verbose
#> [1] FALSE
#> 
#> $num_print
#> [1] 100
#> 
#> $block_eval
#> [1] FALSE
#> 
#> $exact_quadratic
#> [1] TRUE
#> 
#> $generic_accumulate
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "bartisan_control"
#> attr(,"supplied")
#> attr(,"supplied")$num_trees
#>   mean log_sd 
#>     50     10 
#> 

# Weighting the splitting prior toward the treatment, which the sparsity
# prior would otherwise be free to drop
bartisan_control(split_prior = c(rhc = 10))
#> $num_trees
#> NULL
#> 
#> $gate
#> [1] "smoothstep"
#> 
#> $soft
#> [1] TRUE
#> 
#> $sparsity
#> [1] FALSE
#> 
#> $share_sparsity
#> [1] FALSE
#> 
#> $share_forests
#> [1] FALSE
#> 
#> $split_prior
#> rhc 
#>  10 
#> 
#> $categorical
#> [1] "subset"
#> 
#> $k
#> [1] 2
#> 
#> $bandwidth
#> [1] 0.1
#> 
#> $chains
#> [1] 1
#> 
#> $num_burn
#> [1] 200
#> 
#> $num_draws
#> [1] 800
#> 
#> $num_thin
#> [1] 1
#> 
#> $augment
#> [1] "binomial"    "ordinal"     "multinomial" "zip"         "zinb"       
#> [6] "aft"         "negbin"     
#> 
#> $x_transform
#> [1] "smoothcdf"
#> 
#> $gamma
#> [1] 0.95
#> 
#> $beta
#> [1] 2
#> 
#> $sigma_mu
#> NULL
#> 
#> $update_sigma_mu
#> [1] TRUE
#> 
#> $sigma_mu_ramp
#> [1] 0.25
#> 
#> $update_tau
#> [1] TRUE
#> 
#> $update_bandwidth
#> [1] TRUE
#> 
#> $bandwidth_every
#> [1] 1
#> 
#> $alpha
#> [1] 1
#> 
#> $alpha_scale
#> [1] 0
#> 
#> $alpha_shape_1
#> [1] 0.5
#> 
#> $alpha_shape_2
#> [1] 1
#> 
#> $update_s
#> [1] FALSE
#> 
#> $update_alpha
#> [1] FALSE
#> 
#> $verbose
#> [1] FALSE
#> 
#> $num_print
#> [1] 100
#> 
#> $block_eval
#> [1] FALSE
#> 
#> $exact_quadratic
#> [1] TRUE
#> 
#> $generic_accumulate
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "bartisan_control"
#> attr(,"supplied")
#> attr(,"supplied")$split_prior
#> rhc 
#>  10 
#> 
```
