# What is left before shipping, and the vignette series

Two questions, answered together because they turn out to be the same question:
sketching the workflow vignettes is what exposes which gaps are real. A gap that
no vignette needs to work around is not a gap worth holding the release for.

The standard applied throughout: **does this seriously affect a workflow the
package claims to support?** Speedups, polish and features that another package
already covers are excluded, however tempting.

## Part 1: what is actually missing

### Blocking

**1. The package name.** `bartisan` collides case-insensitively with the archived
CRAN package `genBart`. This is a hard block on submission and nothing else on
this list matters until it is settled. Candidates are in `TASKS.md`.

**2. Variable importance has no accessor.** *(done -- `variable_importance()`)* This is the largest genuine gap
against the stated workflow. The information exists and is correct --
`summary()` prints a "Predictor usage" block, and `fit$counts$eta` is a
draws-by-predictor matrix of split counts -- but nothing is exported. Every user
doing variable selection will reinvent

```r
colMeans(fit$counts$eta > 0)          # posterior probability of being used
colMeans(fit$counts$eta)              # mean splits per draw
```

and the second is the wrong one to reach for first. With `sparsity = TRUE` the
first separates cleanly (measured: 1.00, 1.00 for two real predictors against
0.005, 0.010, 0.095 for three noise ones), which is a genuinely useful selection
procedure that is currently undiscoverable.

*Needs:* one exported accessor returning a tidy data frame -- predictor, mean
splits, posterior probability of use, and the DART inclusion probability when
`sparsity = TRUE` -- plus documentation of which column answers which question.
Small work, high leverage, and vignette 3 cannot be written without it.

**3. `as_draws()` omits the fitted values.** *(done -- `eta` argument)* It exposes six variables: `loglik`,
`sigma_mu.eta`, and the nuisance parameters. Not `eta`. Convergence on `eta` is
*already diagnosed* -- `fit$rhat` carries `eta.eta (worst over observations)` --
so the package knows this is the quantity that matters and then declines to hand
it to `bayesplot`. A user wanting a trace plot of a fitted value has to build the
draws array themselves from `posterior_epred()`.

*Needs:* a `variables` argument on `as_draws.bartisan()`, or including a small
representative set of `eta` columns by default. Small work; vignette 2 needs it.

**4. `custom_family()` has no posterior predictive draws.** A log density supplies
no way to draw from it, so `simulate()`, `pp_check()`, `r2()` and everything
built on them are unavailable for a user-written likelihood. Already in `TASKS.md`.
Medium work: an optional `rng` argument alongside the density.

**5. `predict(type = "density")` returns NaN silently** *(done -- it warns; the value stays NaN)* when a composed link's
inverse sends the predictor outside the family's support. Already in `TASKS.md`,
deliberately left alone because it changes `predict()`'s output contract. It is a
silent-wrong-answer class of bug and should be settled before release rather
than after -- a warning is enough.

### Not blocking, and worth saying so explicitly

These are real work, but none of them blocks a workflow the package claims:

- **Partial dependence** -- `marginaleffects::plot_predictions()` covers it.
- **Interactions** -- `comparisons(by = )` covers the useful cases; Friedman's H
  is a dozen lines in a vignette.
- **Model selection** -- `loo()` and `waic()` are present and work.
- **Mixing diagnostics** -- `rhat`, `ess_bulk`, `ess_tail` are present and
  already cover `eta`; only the plotting hand-off is missing (item 3).
- **A separate treatment forest (BCF)** -- a real modeling gap and the right
  post-1.0 feature, but `marginaleffects` covers the estimand side of causal
  inference for a model fitted by hand, and regularization-induced confounding is
  a documentable caveat rather than a broken feature. **The documentation of how
  to do causal inference with what exists is blocking; the feature is not.**
- Correlated random effects across predictors, `quasi()`, the joint ordinal
  cutpoint update, grow-from-root warm start, categorical splits on level
  subsets, soft random tree features -- all post-1.0.

### Deliberately not on the list

Further speedups. `weibull_aft()` at ~8 seconds for 700 observations is the
slowest family and could plausibly reach ~1 second with a Frühwirth-Schnatter
mixture-of-normals augmentation of its Gumbel error -- but that augmentation is
an *approximation*, and adopting it would weaken the exactness claim the package
currently makes. That is a design decision, not a bug, and not a release blocker
now that `weibull_aft()` is no longer the default.

### Measured after the fact: what the other BART packages export

Surveyed `dbarts`, `BART`, `bartMachine` and `SoftBart` by their actual export
lists rather than from memory. Every recurring non-fitting feature, and where it
lives here:

| They export | Here |
| --- | --- |
| partial dependence (`pdbart`, `pd_plot`, `pdsoftbart`) | `marginaleffects::plot_predictions()` |
| variable importance (`investigate_var_importance`, `get_var_props_over_chain`, `posterior_probs`) | `variable_importance()` |
| variable selection (`var_selection_by_permute`, `mc.wbart.gse`) | `prop_used` with `sparsity = TRUE` -- a posterior probability rather than a permutation null |
| credible intervals (`calc_credible_intervals`) | `marginaleffects::predictions()` |
| prediction intervals (`calc_prediction_intervals`) | `posterior_predict()` |
| convergence plots (`plot_convergence_diagnostics`) | `fit$rhat`, `as_draws()` + bayesplot |
| fit checks (`plot_y_vs_yhat`, `check_bart_error_assumptions`) | `pp_check()`, `residuals()` |
| interaction detection (`interaction_investigator`) | `marginaleffects::comparisons(by = )`, partially |
| k-fold CV (`k_fold_cv`, `xbart`, `bartMachineCV`) | **absent**; `loo()` covers model comparison |
| model matrix helpers (`makeModelMatrixFromDataFrame`, `dummify_data`, `preprocess_df`, `bartModelMatrix`) | **not needed** -- the formula interface does it |

Two conclusions.

**There is no remaining feature gap.** After `variable_importance()`, everything
those four packages export is either present, covered by a better-maintained
general package, or unnecessary here. The only genuine absence is k-fold
cross-validation, and it is absent for a reason: `xbart` and `bartMachineCV`
exist to tune `k`, `q` and the tree count, which the defaults here are meant to
settle, and `loo()` covers model comparison without refitting.

**What is left is discoverability, not function.** A newcomer cannot guess that
intervals come from `marginaleffects::predictions()`, or that a prediction
interval needs `posterior_predict()`. Four of those packages ship one function
per task with an obvious name; this one ships a map. The map has to be good --
hence the rewritten `?bartisan-package` and the vignette series below.

## Part 2: the vignette series (written)

Three vignettes exist: `bartisan` (the model and the sampler), `families`
(choosing a likelihood), `survival` (censored responses). Those are
*reference* documents, organized by the package's own structure.

The series below is organized by **the analyst's workflow instead** -- the arc a
person actually walks from a data frame to a defensible claim. Each answers a
question someone would ask about a regression, and answers it with BART.

A design rule for all of them: they must be honest about the places where the
BART answer is *worse* or *harder* than the `lm()` answer, not only where it is
better. A reader who finds out later that BART's variable importance is unstable
under correlated predictors, having not been told, will not trust the rest.

---

### 1. `workflow` (done) -- "A regression workflow with BART"

**The orientation piece, and the one to write first.** One dataset, one question,
start to finish, with every step deliberately shallow and a pointer to the
vignette that goes deep. Someone who reads only this should be able to do a
competent analysis.

The arc: fit → did it converge → does it fit → what matters → what is the effect
→ is this model better than that one → report it.

Deliberately establishes the two habits the rest of the series depends on: **read
counterfactual estimands through `marginaleffects`, never coefficients** (BART
has none to read), and **check convergence on the fitted values, not just on
`sigma`**.

Ends with the comparison to `lm()`/`glm()` stated plainly: what you give up is
the coefficient table and the ability to state the functional form; what you get
is not having to state the functional form. Nothing about a BART fit is a
substitute for knowing what question you are asking.

### 2. `diagnostics` (done) -- "Has it converged, and does it fit?"

Two questions usually run together and worth separating.

*Converged*: what R-hat and ESS mean when the parameter is a function rather than
a number, why `eta`'s worst-over-observations R-hat is the one to watch, running
`chains > 1` and why the default of 1 is not an endorsement, trace and rank plots
through `posterior` and `bayesplot`, and what to do when it has not converged
(more draws, more trees, harder rules).

*Fits*: `pp_check()`, residuals for the families that have meaningful ones, and
the specific failure BART shows when it is under-fitting -- shrinkage toward the
mean at the edges of the predictor space, which looks like bias and is.

Also the honest part: a BART fit can have excellent R-hat on `sigma` and a badly
mixed forest, and the tree structure itself is not identified, so "did the chains
find the same trees" is not a question worth asking. Only functionals of the fit
are.

*Needs item 3 above.*

### 3. `importance` (done) -- "Which variables matter?"

Split counts, `prop_used`, and the DART posterior inclusion probability, with the
distinction between them made sharp: **how often a variable is split on is not
how much it matters**, and neither is a test.

Turning `sparsity = TRUE` into an actual selection procedure, including the
measured separation. Then the three cautions that make this vignette worth
writing rather than a paragraph:

- Under **correlated predictors** the split counts distribute arbitrarily among
  the correlated set; a variable can be genuinely important and rarely split on
  because a collinear partner absorbed it.
- Importance is **not an effect**: a variable can be split on constantly and move
  the prediction very little. `marginaleffects` answers "how much does it move
  the prediction", and that is usually the real question.
- Importance is **not causal**, and the vignette should say so before someone
  reads a ranking as a set of causes.

*Needs item 2 above.*

### 4. `effects` (done) -- "Reading the fitted function"

The vignette that replaces the coefficient table. Predictions, comparisons and
slopes through `marginaleffects`; partial dependence via `plot_predictions()`;
average versus conditional effects and why the average is usually the one to
report.

Then interactions, which is the part with no ready-made tool: `comparisons(by =)`
to get an effect at levels of a moderator, and how to tell a real interaction
from the shrinkage artifact that produces a similar picture -- the honest answer
being that the posterior interval on the *difference of differences* is the test,
not the eyeball.

Nonlinearity: how to show a fitted curve with an interval, and why a soft-rule
fit is the one to show.

### 5. `comparison` (done) -- "Choosing between models"

`loo()` and `waic()`, what they compare and what they cannot. The measure trap
from the survival vignette generalizes: **a log score is only comparable between
models that put their density on the same scale**, and the package has at least
one pair that does not.

What model selection *means* here, which differs from the `lm()` case and is the
reason this is a vignette rather than a paragraph: you are not selecting a mean
structure, because the forest selects that. You are selecting the **likelihood** --
the family, the link, whether the variance is constant, whether the error shape
is assumed. That reframing is the content.

Also: when `loo()` breaks and what it means (the `ph()` bin-count result is the
worked example), and why comparing a BART fit to a linear model by `loo` is a
fair and useful thing to do.

### 6. `causal` (done) -- "Causal inference with BART"

The one where being clear about what the package does *not* do matters most.

G-computation through `avg_comparisons()`: the ATE, the ATT, and effects on the
scale the question is asked in. Overlap and positivity as things to check before
the estimate, not after. Weighting through `WeightIt` when the design calls for
it, and the propensity score as a covariate (Hahn, Murray and Carvalho).

Then the caveat that has to be in this vignette rather than a footnote:
**regularization-induced confounding**. BART's prior shrinks the treatment effect
toward zero along with everything else, and when treatment assignment is strongly
predicted by the covariates that shrinkage biases the effect. BCF's separate
treatment forest is the standard answer and `bartisan` does not have one. Say what
the bias looks like, when it bites, and what to do meanwhile.

Uncertainty: the posterior interval on an average effect is a credible interval
for the estimand under the model, and is not a confidence interval for a causal
effect unless the identification assumptions hold. Worth one paragraph, plainly.

---

### Sequencing

1 first: it is the entry point and the others are its expansions. Then 4, which
is the highest-value single document (most people's real question is "what is the
effect") and needs no new package code. Then 3 and 2, each after its blocking
item lands. Then 5. Then 6, last, because it is the one where being wrong is most
costly and it benefits from the other five being settled.

`families` and `survival` stay as reference documents and get cross-links from 1
and 5. `bartisan` stays as the "how it works" document; it is the only one
organized around the method rather than the workflow, which is correct for it.
