# Categories B and D in depth, and a BART specification for staggered DiD

Notation as before: $i$ indexes units, $t = 1,\dots,T$ calendar periods, $G_i$ is unit $i$'s cohort
(first treated period, $G_i = \infty$ for never-treated), $D_{it} = \mathbf{1}\{t \ge G_i\}$,
$E_i = \mathbf{1}\{G_i < \infty\}$, $K_{it} = t - G_i$, $X_i$ baseline covariates,
$\mathcal{G}$ the set of finite cohorts.

---

## 1. Category B — one fit on all cells, treatment saturated in cohort $\times$ period

### 1.1 The defining algebra

Category B models the *observed* outcome on every cell of the panel, and relies on the treatment
block carrying enough free parameters that treated observations cannot influence the baseline. The
generic form is

$$
Y_{it} \;=\; \underbrace{b(i,t,X_i)}_{\text{baseline}} \;+\; \sum_{g \in \mathcal{G}} \mathbf{1}\{G_i = g\}\,\mathbf{1}\{t \ge g\}\, \underbrace{\big[\tau_{gt} + \rho_{gt}(X_i)\big]}_{\text{treatment block}} \;+\; \varepsilon_{it}.
$$

The claim that makes B legitimate is that the $\hat\tau_{gt}$ from this pooled fit equal the
Category A imputation estimates — fit the baseline on $\{(i,t) : D_{it} = 0\}$ only, then average
the gaps $Y_{it} - \hat b(i,t,X_i)$ within each treated cell. Wooldridge proves this for the
canonical link. I checked the algebra directly on simulated panels ($T=8$, cohorts $\{4,6,7\}$),
comparing pooled $\hat\tau_{gt}$ against the imputation estimate cell by cell:

| baseline $b(i,t,X_i)$ | treatment block | $\max_{g,t}\lvert \hat\tau^{\text{pooled}}_{gt} - \hat\tau^{\text{imp}}_{gt}\rvert$ |
|---|---|---|
| $\nu_g + \lambda_t$ | cell dummies | $9 \times 10^{-14}$ — **identical** |
| $\alpha_i + \lambda_t$ | cell dummies | $5 \times 10^{-14}$ — **identical** |
| $\alpha_i + \lambda_t + X_i'(\kappa + \xi_t)$ | cell dummies | $1.4 \times 10^{-2}$ — **differ** |
| $\alpha_i + \lambda_t + X_i'(\kappa + \xi_t)$ | cells $+$ *uncentred* $X_i \times$ cell | $8.3 \times 10^{-2}$ — **differ** |
| $\alpha_i + \lambda_t + X_i'(\kappa + \xi_t)$ | cells $+$ *within-cell-centred* $X_i \times$ cell | $4 \times 10^{-14}$ — **identical** |

Two distinct lessons, and they are the whole of Category B:

**(i) Saturation must cover every direction in which the baseline is flexible.** With a baseline of
only unit/cohort and period effects, a free scalar per treated cell is enough: the cell dummies span
the baseline restricted to each treated cell, so the normal equations for the baseline reduce to the
untreated-cell equations. Add $X_i'\xi_t$ to the baseline and that is no longer true — treated
observations now carry information about $\xi_t$ for $t \ge g$, they leak into the baseline, and the
equivalence breaks (row 3). Restoring it requires the treatment block to carry a covariate
interaction for every treated cell too (row 5).

**(ii) The interaction must be centred within the cell, or $\tau_{gt}$ stops being the ATT.** Row 4
adds the covariate interactions but uncentred, and the discrepancy gets *worse*, not better — because
$\tau_{gt}$ is then the effect at $X_i = 0$, a different estimand from the cell-average effect. This
is exactly why Wooldridge writes his treatment interactions in cohort-centred form,
$\dot X_{i,g} = X_i - \mathbb{E}[X_i \mid G_i = g]$, giving

$$
\mathrm{ATT}(g,t) = \tau_{gt}, \qquad \mathrm{CATT}(g,t,x) = \tau_{gt} + (x - \mathbb{E}[X_i \mid G_i = g])'\rho_{gt}.
$$

Centring is not cosmetic. It is what buys the coefficient its interpretation.

### 1.2 This is precisely the diagnosis of DiD-BCF, stated generally

Lesson (i) gives the general principle that DiD-BCF violates. Its baseline $\mu(E_i, t, X_{it})$ is
flexible in the $E_i \times t$ direction, and its treatment block $\tau(X_{it}, K_{it}) \cdot D_{it}$
does not — cannot — saturate that direction, because $\tau$ never sees $G_i$ and there is no free
parameter per $(g,t)$ cell. The baseline is flexible in a direction the treatment block does not
cover, so the effect leaks into the baseline. In the single-adoption-date case the leak is total
(the $D_{it}$ column lies exactly in the baseline's span); under staggered adoption it is partial,
and the prior decides how much leaks.

### 1.3 The cost of saturation, and the two responses

ETWFE's full specification is

$$
Y_{it} = \eta + \lambda_t + X_i'(\kappa + \xi_t) + \sum_{g \in \mathcal{G}} \mathbf{1}\{G_i = g\}\Big[\nu_g + X_i'\zeta_g + \mathbf{1}\{t \ge g\}\big(\tau_{gt} + \dot X_{i,g}'\rho_{gt}\big)\Big] + \varepsilon_{it}.
$$

With $d$ covariates, $\lvert\mathcal{G}\rvert$ cohorts and $T$ periods, the treatment block alone
carries roughly $\tfrac{1}{2}\lvert\mathcal{G}\rvert T (1 + d)$ parameters. Saturation buys freedom
from contamination and pays for it in variance. The literature offers two ways to buy some
efficiency back without reintroducing bias by hand:

- **Fusion penalty (FETWFE, Faletto).** Penalise differences between adjacent coefficients,
  $\lvert \tau_{gt} - \tau_{g,t-1}\rvert$ and the covariate analogues, with a single tuning
  parameter, so equality restrictions are selected from the data rather than imposed. Faletto's
  motivation is explicitly that ETWFE "adds many parameters" and that *ad hoc* restrictions "may
  reintroduce bias".
- **Partition prior (Arora & Wagle).** Put a Dirichlet-process mixture prior on
  $\{\tau_{gt}\}$, which clusters cells with equal effects without fixing the number of clusters,
  returns co-clustering probabilities for every pair, and yields intervals that marginalise over the
  unknown partition. Their MAP partition under a pairwise penalty reduces to an
  $\ell_0$-penalised regression, so it is the Bayesian counterpart of the fusion penalty.

This is the natural home for a Bayesian package: the frequentist version needs a tuning parameter
chosen by cross-validation, whereas a prior over the partition of the $\tau$ array does the same job
and propagates the uncertainty into the intervals automatically.

### 1.4 The nonlinear-link version

Wooldridge (2023) carries B to non-Gaussian outcomes. Parallel trends is assumed on the **index
scale**: with link $h(\cdot)$,

$$
h\big(\mathbb{E}[Y_{it}(0) \mid X_i, G_i = g]\big) = \eta + \lambda_t + X_i'(\kappa + \xi_t) + \nu_g + X_i'\zeta_g,
$$

estimated by pooled quasi-maximum likelihood in the linear exponential family. His result: using the
conditional mean associated with the **canonical** link, imputation and pooling over the whole sample
give identical estimates. Leading cases are a logit mean with the Bernoulli quasi-log-likelihood for
binary and fractional outcomes, and an exponential mean with the Poisson QLL for nonnegative
outcomes. Two consequences: the equivalence in §1.1 is a canonical-link result, so with a
non-canonical link you must choose imputation deliberately; and the ATT must be formed on the outcome
scale by averaging imputed counterfactuals, never by back-transforming an average of index-scale
predictions.

### 1.5 Summary of what B requires

| requirement | why |
|---|---|
| baseline additive (or otherwise parametric) in $(i,t)$ | so a free scalar per cell can saturate it |
| a free parameter per treated $(g,t)$ cell | prevents the cell mean leaking into the baseline |
| an interaction per treated cell for every flexible covariate direction | prevents covariate-direction leakage |
| interactions centred within cohort/cell | keeps $\tau_{gt}$ equal to $\mathrm{ATT}(g,t)$ |
| parallel trends on the link scale; canonical link for pooling | equivalence and outcome-scale ATT |

What you get in return: all estimands are functionals of one fit, and $\mathrm{CATT}$ is available
analytically from the interaction coefficients.

---

## 2. Category D — saturate the group $\times$ time surface, then difference the fitted intercepts

### 2.1 The mechanics

Karim & Webb's DID-INT runs, with no constant,

$$
Y_{ist} = \sum_{s}\sum_{t} \lambda_{st}\, I(s,t) \;+\; \sum_{k=1}^{K} f_k\big(X^k_{ist}\big) \;+\; \varepsilon_{ist},
$$

where $I(s,t)$ is the group-by-period intersection indicator. Note what is absent: **there is no
treatment-effect parameter anywhere in the regression.** The $\hat\lambda_{st}$ are simply
covariate-adjusted group-by-period cell means. Effects are constructed afterwards, from the long
difference against each cohort's last pre-treatment period $t^{-s}$ and a comparison group $s'$:

$$
\hat\theta_{s,t} = \big(\hat\lambda_{s,t} - \hat\lambda_{s,t^{-s}}\big) - \big(\hat\lambda_{s',t} - \hat\lambda_{s',t^{-s}}\big),
\qquad
\hat\theta = \sum_{s}\sum_{t} \mathbf{1}\{t^s \le t\}\, w_{s,t}\, \hat\theta_{s,t}.
$$

Their underlying potential-outcome model is

$$
Y_{ist}(0) = \sum_{k} \gamma^k_{s,t}\, f\big(X^k_{ist}\big) + \alpha_i + \delta_t + \varepsilon_{ist},
$$

so covariate effects $\gamma^k_{s,t}$ are indexed by group *and* period — their relaxation of the
two-way common-causal-covariates condition — while $f$ is a flexible transform applied to each
covariate **separately**, additive across $k$, not a joint surface. Parallel trends is assumed on the
covariate-adjusted residuals.

### 2.2 Why an arbitrarily flexible surface is safe here

In Category B, identification lives in the *parameterisation*: $\tau_{gt}$ is a coefficient, and the
restriction that makes it meaningful is the additivity of the baseline. In Category D, identification
lives in the *contrast you form afterwards*. The estimation step is deliberately assumption-free —
a saturated cell-mean model cannot be misspecified in $(s,t)$, because it has a parameter for every
cell. Parallel trends enters only when you decide that

$$
\mathrm{ATT}(s,t) = \big(\lambda_{s,t} - \lambda_{s,t^{-s}}\big) - \big(\lambda_{s',t} - \lambda_{s',t^{-s}}\big)
$$

is the right contrast — that is, that absent treatment cohort $s$ would have moved like $s'$.

Because there is no competing treatment parameter, there is nothing for the surface to steal from.
You may make the $(s,t)$ surface as flexible as you like. This is the structural reason D is the only
family in which a *flexible* group-by-time baseline is legitimate, and it is what DiD-BCF should have
done.

Callaway–Sant'Anna belongs to the same logic: estimate conditional cell means (or long differences)
on the comparison group, then contrast. D is that idea with the cells pooled into one regression.

### 2.3 B and D are the same projection in two coordinate systems

In the linear unregularised case with matching covariate specifications, B and D are
reparameterisations of each other. The saturated cell means $\lambda_{st}$ and the ETWFE coefficients
$(\nu_g, \lambda_t, \tau_{gt})$ span the same column space; B writes the DiD contrast *into* the
design so that a coefficient equals the estimand, D leaves the design descriptive and applies the
contrast to fitted values. That is why the numbers in §1.1 come out identical to machine precision:
both are the same orthogonal projection, read off differently.

They diverge the moment anything is regularised or made nonparametric, and *where* you put the
flexibility determines which form survives:

| where flexibility / shrinkage goes | B valid? | D valid? |
|---|---|---|
| covariates only, $(i,t)$ additive | yes | yes |
| the group $\times$ time surface itself | **no** — the surface competes with $\tau_{gt}$ | yes — no competitor exists |
| the effects $\{\tau_{gt}\}$ (fusion, DP partition) | yes — they are parameters to shrink | awkward — effects are not parameters |
| both surface and effects | no | no |

So the design decision is not "B or D is better" but: *do you want a prior on the effects, or a
flexible baseline in $(i,t)$?* You can have either, not both. For a package whose selling point is a
prior over the $\mathrm{ATT}(g,t)$ array, B is the right frame, with the $(i,t)$ part kept
parametric.

---

## 3. A BART specification with parametric components

### 3.1 Precedent: BART models that carry non-tree parameters

This is an established architecture, not an improvisation. Four relevant precedents, all retrieved
and checked:

- **semi-BART** — Zeldow, Lo Re & Roy (arXiv 1806.04200). Confounders "not of scientific interest"
  are given unspecified functional form via BART; treatment and covariates of substantive interest
  get "the usual linear form from parametric regression", and the posterior of the linear part "can
  be interpreted as in parametric Bayesian regression". Exactly the split needed here.
- **General BART** — Tan & Roy (arXiv 1901.07504). A framework unifying semiparametric BART,
  correlated outcomes, and weaker distributional assumptions — i.e. the machinery for adding
  parametric and random-effect blocks to a tree ensemble.
- **stan4bart** — two separate records. The R package is *stan4bart: Bayesian Additive Regression
  Trees with Stan-Sampled Parametric Extensions* (Vincent Dorie, 2021, single-authored); the
  accompanying methods paper is *Stan and BART for Causal Inference: Estimating Heterogeneous
  Treatment Effects* (Dorie, Perrett & Hill, 2022). BART for the nonparametric part, Stan for fixed
  and multilevel parametric terms, sampled jointly — the closest existing engine to what is
  described below.
- **VCBART** — Deshpande, Bai, Balocchi, Starling & Weiss (arXiv 2003.06416). A linear
  varying-coefficient model $Y = \sum_j \beta_j(Z)\, W_j$ in which each coefficient function is a
  BART ensemble over effect modifiers. The staggered-DiD effect block below is a varying-coefficient
  model of exactly this shape.

Two cautions from the same literature. Prado, Parnell & Murphy ("Accounting for shared covariates in
semi-parametric BART") address the case where a covariate appears in *both* the linear and the tree
part — which will happen here, since $X_i$ enters both the baseline and the effect block; this is a
mixing and identification hazard worth reading before implementing. And Prevot, Häring, Nichols,
Holmes & Ganjgahi (arXiv 2508.08418) note that BART and BCF "assume independence across observations
and are fundamentally limited in their ability to model within-individual correlation over time",
which is the defect I measured in DiD-BCF (residual SD inflated from $0.50$ to $1.11$, i.i.d.
standard error for $\tau$ understating the unit-clustered one by $1.47\times$).

### 3.2 The proposed model

$$
Y_{it} \;=\; \underbrace{\alpha_i + \lambda_t}_{\text{parametric}} \;+\; \underbrace{f_0\big(X_i,\, t\big)}_{\text{BART}} \;+\; D_{it}\Big[\; \underbrace{\tau_{G_i,\,t}}_{\text{parametric}} \;+\; \underbrace{f_\tau\big(\dot X_{i,G_i},\, G_i,\, K_{it}\big)}_{\text{BART, centred}} \Big] \;+\; \varepsilon_{it},
$$

$$
\dot X_{i,g} = X_i - \bar X_g, \qquad \bar X_g = \frac{1}{n_g}\sum_{i : G_i = g} X_i .
$$

**Parametric blocks (not in any tree):**

$$
\alpha_i \sim \mathcal{N}(0, \sigma_\alpha^2), \qquad \sigma_\alpha \sim \text{half-}\mathcal{N}, \qquad \lambda_t \ \text{free}, \qquad \varepsilon_{it} \sim \mathcal{N}(0, \sigma^2).
$$

$\alpha_i$ as a hierarchical random effect rather than a fixed effect is the deliberate choice: it
supplies the within-unit correlation whose absence made DiD-BCF's intervals too narrow, and a treated
unit's $\alpha_i$ is then shrunk toward the population mean using only its pre-treatment periods.
$\{\tau_{gt}\}$ are free scalars — one per treated cell, which is the saturation requirement from
§1.5 — carrying a shrinkage prior across the array (§3.5).

**Nonparametric blocks:**

$f_0$ is the $Y(0)$ covariate surface. It receives $X_i$ and $t$ and therefore admits the
covariate-by-time interaction that ETWFE writes as $X_i'\xi_t$ — but nonparametrically, which is the
actual gain over ETWFE. **It must never receive $E_i$, $G_i$, or $D_{it}$.** That single restriction
is the difference between this model and DiD-BCF.

$f_\tau$ is the effect-heterogeneity surface, a BCF-style second ensemble with tighter
regularisation, over cohort-centred covariates, cohort, and event time. It replaces ETWFE's
$\dot X_{i,g}'\rho_{gt}$ with a nonparametric function, and lets effects borrow strength across
cohorts and event times instead of having a free parameter per cell.

### 3.3 Why this is identified where DiD-BCF is not

The baseline $\alpha_i + \lambda_t + f_0(X_i,t)$ has no access to treatment status, cohort, or
ever-treated status. It cannot represent a group-by-time interaction *except* through $X_i$ — which
is the legitimate conditional-parallel-trends channel, not a leak. The treatment block supplies a
free scalar per $(g,t)$ cell, so the cell means of the treated block are saturated exactly as
§1.5 requires and cannot be absorbed by the baseline.

The residual risk is the covariate-direction leakage of §1.1 row 3, and it is now an *overlap*
condition rather than a structural defect. $f_0$ can only absorb effect variation in a region of
$(X_i, t)$ space where treated units have no untreated counterparts. So the assumption to state and
check is

$$
0 < \Pr\big(G_i = g \mid X_i = x\big) < 1 \quad \text{for all } x \text{ in the support of interest},
$$

the same overlap condition Callaway–Sant'Anna require. This is what went wrong in DiD-BCF in the
sharpest possible way: giving $\mu$ the indicator $E_i$ makes that propensity exactly $0$ or $1$
everywhere, so there is no overlap at all and nothing is identified. With $f_0(X_i, t)$, overlap is
an empirical property you can diagnose.

Two practical safeguards follow. Enforce the centring of $f_\tau$ numerically — subtract the
within-cell mean of $f_\tau$ from the ensemble and add it to $\tau_{gt}$ at each iteration — so that
$\tau_{gt}$ retains its $\mathrm{ATT}(g,t)$ meaning (§1.1 lesson ii). And because the baseline is
conditionally uninformed by treated cells given the $\tau_{gt}$, the Gibbs update for
$(\alpha, \lambda, f_0)$ may be run on untreated cells alone — not as an approximation but as the
exact conditional, which also makes the Category A behaviour available as an option.

### 3.4 Reading the estimands off one fit

Do not read the ATT off the coefficient; average the fitted effect. For each posterior draw $m$:

$$
\mathrm{ATT}^{(m)}(g,t) = \frac{1}{n_g}\sum_{i : G_i = g} \Big[\tau^{(m)}_{gt} + f^{(m)}_\tau\big(\dot X_{i,g},\, g,\, t-g\big)\Big],
$$

$$
\mathrm{CATT}^{(m)}(g,t,x) = \tau^{(m)}_{gt} + f^{(m)}_\tau\big(x - \bar X_g,\, g,\, t-g\big),
$$

and any aggregate is a weighted average of the cell-level draws,

$$
\theta^{(m)} = \sum_{g \in \mathcal{G}} \sum_{t \ge g} w(g,t)\, \mathrm{ATT}^{(m)}(g,t),
$$

with $w(g,t) \propto \mathbf{1}\{t - g = e\}\Pr(G_i = g)$ for event-time $e$, $w \propto \mathbf{1}\{g = g_0\}$
for a cohort effect, and $w \propto \Pr(G_i = g)$ over all treated cells for the overall ATT. Because
these are functionals of the same posterior draws, credible intervals for every aggregate are
coherent and come free — the one place where this design beats the doubly-robust machinery outright.
Pre-treatment gaps $Y_{it} - \hat Y_{it}(0)$ for $K_{it} < 0$ are the built-in placebo.

### 3.5 Prior on the effect array

The $\{\tau_{gt}\}$ are the natural place for the shrinkage that §1.3 motivates. Options in
increasing sophistication: independent $\mathcal{N}(0, \sigma_\tau^2)$ with $\sigma_\tau$ learned;
a random walk in event time within cohort, $\tau_{g,t} = \tau_{g,t-1} + \omega_{gt}$, which is the
Bayesian analogue of FETWFE's fusion penalty; or a Dirichlet-process partition prior over the cells
following Arora & Wagle, which additionally returns co-clustering probabilities for pairs of cells.
Keep the prior on $f_\tau$ tighter than on $f_0$ as BCF does — that asymmetry is *correct* here,
because with the saturated $\tau_{gt}$ and no treatment indicator in $f_0$ there is no longer a
collinear direction for the effect to leak into.

### 3.6 Generalized outcomes

For bartisan's non-Gaussian links, state parallel trends on the link scale,

$$
h\big(\mathbb{E}[Y_{it}(0) \mid X_i]\big) = \alpha_i + \lambda_t + f_0(X_i, t),
$$

and form effects on the outcome scale as $h^{-1}(\cdot)$ evaluated with and without the treatment
block, averaged over units:

$$
\mathrm{ATT}^{(m)}(g,t) = \frac{1}{n_g}\sum_{i:G_i=g}\Big[ h^{-1}\big(\hat\eta^{(m)}_{it} + \tau^{(m)}_{gt} + f^{(m)}_\tau\big) - h^{-1}\big(\hat\eta^{(m)}_{it}\big)\Big],
$$

with $\hat\eta_{it}$ the baseline index. Never back-transform an average of index-scale predictions.
Note that under a nonlinear link the $\tau_{gt}$ scalars are index-scale parameters, so the
centring of $f_\tau$ no longer makes $\tau_{gt}$ equal the ATT — the averaged-contrast formula above
is then the only correct route, and the canonical link is worth preferring because it is where
Wooldridge's pooling/imputation equivalence holds.

### 3.7 Sampler notes for the Linero (2025) engine

The additions to a generalized-BART reversible-jump sampler are three conditionally conjugate or
Laplace-approximated blocks interleaved with the tree updates: $(\lambda_t, \tau_{gt})$ jointly,
which is a linear-Gaussian block update in the identity-link case and a Laplace step otherwise;
$\alpha_i$ with its variance component, a standard hierarchical update; and the two ensembles updated
by backfitting against the current partial residual (the `stan4bart` package implements exactly this
interleaving for Gaussian and binary outcomes and is worth reading as a reference implementation).
The $\tau_{gt}$ block and the $f_\tau$ ensemble
are the pair most at risk of poor mixing, since both describe the treated cells — which is where the
shared-covariate problem flagged by Prado et al. bites, and where the within-cell recentring step of
§3.3 also helps by making the two blocks approximately orthogonal.

---

## 4. One-line summary

Category B writes the DiD contrast into the design and requires the baseline to be parametric in
$(i,t)$ and the treatment block to saturate every flexible baseline direction, centred within cell.
Category D leaves the design descriptive, saturates the group-by-time surface, and imposes parallel
trends only when contrasting fitted values — which is why it alone tolerates a flexible baseline. A
BART model for staggered DiD should adopt B's frame: parametric $\alpha_i + \lambda_t$ and
saturated $\tau_{gt}$ carrying a shrinkage prior, a BART $Y(0)$ surface over covariates and time that
never sees treatment or cohort, and a second, more tightly regularised BART ensemble for effect
heterogeneity over cohort-centred covariates, cohort, and event time.
