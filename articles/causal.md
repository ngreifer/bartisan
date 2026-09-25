# Causal Inference with BART

## Introduction

BART is a popular choice for causal inference because it allows one to
fit the nuisance functionals required for effect estimation (the
relationships among the outcome, the treatment, and the potential
confounders) in a flexible manner ([Hill 2011](#ref-hill2011)). BART has
repeatedly been shown to outperform other effect estimation methods in
competitions ([Dorie et al. 2019](#ref-dorie2019)). In this vignette, we
demonstrate how to use BART implemented in *bartisan* to estimate
treatment effects under the assumption of no unmeasured confounding.
This includes both standard BART as well as Bayesian causal forests
(BCF), which is an improvement on traditional BART that works by fitting
separate BART models for the outcome under control and for the treatment
effect itself as a varying coefficient model.

However, it’s important to remember that BART is not a “causal inference
method”; it is just a method of estimating certain quantities, which
often are interpretable as associations. The assumptions that turn an
association into a causal effect are assumptions about the design, and
no model supplies them.

In addition to its use here, BART can be used with instrumental
variables analysis ([McCulloch et al.,
n.d.](#ref-mccullochCausalInferenceInstrumental2021)) and regression
discontinuity ([Alcantara et al.
2024](#ref-alcantaraModifiedBARTLearning2024)).

The example we use here is a dataset used to answer whether right heart
catheterization helps or harms critically ill patients ([Connors et al.
1996](#ref-connors1996)).

``` r

library(bartisan)

data(rhc)

set.seed(2026)
```

## Confounding

Catheterization was not randomized. Doctors chose it, and they chose it
more often for patients who were already doing badly. Those patients
were also more likely to die. We can use
[`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html)
to see how the patients differ between groups:

``` r

bal.tab(rhc ~ age + sex + race + edu + aps + meanbp + resp +
          hema + pafi + paco2 + crea + surv2m + card,
        data = rhc, stats = c("m", "ovl"), disp = "m")
#> Balance Measures
#>               Type  M.0.Un  M.1.Un Diff.Un OVL.Un
#> age        Contin.  61.816  60.756  -0.065  0.101
#> sex_male    Binary   0.524   0.591   0.067  0.067
#> race_white  Binary   0.793   0.768  -0.024  0.024
#> race_black  Binary   0.157   0.168   0.011  0.011
#> race_other  Binary   0.050   0.064   0.013  0.013
#> edu        Contin.  11.568  11.767   0.061  0.048
#> aps        Contin.  51.369  62.099   0.530  0.223
#> meanbp     Contin.  85.302  66.871  -0.505  0.231
#> resp       Contin.  28.841  27.465  -0.097  0.055
#> hema       Contin.  32.545  30.237  -0.290  0.159
#> pafi       Contin. 238.232 182.995  -0.504  0.208
#> paco2      Contin.  40.045  36.871  -0.262  0.096
#> crea       Contin.   1.905   2.506   0.292  0.179
#> surv2m     Contin.   0.604   0.559  -0.229  0.105
#> card_yes    Binary   0.299   0.405   0.106  0.106
#> 
#> Sample sizes
#>     Control Treated
#> All     935     565
```

The `M.0.Un` column indicates the mean for each variable in the control
group and the `M.1.Un` column indicates the mean for each variable in
the treated group. `Diff.Un` and `OVL.Un` are measures of the
distributional difference between the groups for each covariate; values
far from 0 indicate imbalance due to differential selection into
treatment. In particular, we can see that patients with higher values of
`aps` and `crea` and lower values of `meanbp`, `hema`, `pafi`, `paco2`,
and `surv2m` are overrepresented among treated units.

Given that sicker patients are both more likely to die and more likely
to receive RHC, it would not be unexpected to see that patients who
receive RHC are more likely to die, and we do see that:

``` r

with(rhc, tapply(death, rhc, mean))
#>      0      1 
#> 0.6193 0.7115
```

However, that doesn’t mean RHC causes death; to disentangle the effects
of RHC from the confounding effects of patients’ characteristics, we
need to adjust for these characteristics.

## Assumptions for Causal Inference

Several assumptions are required to interpret an adjusted effect
estimate as causal, none of which the fit can check:

1.  **No unmeasured confounding.** Every common cause of catheterization
    and death is in the model. Here that is the crux: the covariates
    include a physiological profile and the study’s own prognostic
    score, which is a serious attempt, but a doctor’s judgment at the
    bedside may not be fully captured by these thirteen variables.

2.  **Positivity.** Every kind of patient could have received the
    procedure or not. This one is partly checkable and is checked below.

3.  **Consistency.** “Catheterization” names a single well defined
    intervention.

4.  **No interference.** One patient’s treatment does not affect
    another’s outcome.

The first is the one that usually fails and the one to be explicit
about. State it as an assumption in what you write, rather than letting
the interval imply it has been handled.

### Checking Positivity

Positivity fails when some combination of covariates makes the treatment
nearly certain. To look for that, model the treatment and inspect the
fitted probabilities in both groups. BART is a good choice for this
model too, for the same reason it is a good choice for the outcome: the
functional form is not known.

``` r

ps_fit <- bartisan(
  rhc ~ age + sex + race + edu + aps + meanbp + resp + hema + pafi +
    paco2 + crea + surv2m + card,
  data = rhc, family = binomial()
)

prop_score <- fitted(ps_fit)
```

We can use
[`cobalt::bal.plot()`](https://ngreifer.github.io/cobalt/reference/bal.plot.html)
to examine the overlap of the propensity score distribution between the
groups:

``` r

bal.plot(rhc ~ prop_score, data = rhc, type = "hist", mirror = TRUE)
```

![](causal_files/figure-html/balplot-1.png)

What matters for positivity is that the two groups overlap over most of
their range, and they do. Distributions pushed against zero and one with
little overlap indicate a failure of positivity, and the honest response
then is to restrict the analysis to the region of overlap rather than
let the model extrapolate into a part of the covariate space where one
treatment was never observed.

This check belongs before the outcome model, not after. Where one
treatment was never observed, any estimate the outcome model produces
there is an extrapolation, and whether its interval shows that depends
on the model. On the example Li et al.
([2023](#ref-liBayesianCausalInference2023)) use to illustrate the
problem, a BART model fit separately to each treatment group kept the
same interval width where the groups did not overlap, and missed the
true effect there. In our replication of that example, the models used
in this vignette widened their intervals four- to ninefold in those
regions, but the wider intervals did not remove the bias:
[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md)’s
estimates there were still off by more than 1 on an effect of 5. So a
wide interval is worth noticing, but it does not replace checking
overlap before fitting the outcome model.

## The Estimator: G-Computation

To estimate the average treatment effect \\\tau\_{\text{ATE}} =
E\[Y(1)\] - E\[Y(0)\]\\, which is a function of the unobserved potential
outcomes \\Y(1)\\ and \\Y(0)\\, we can use the assumptions above to
express it as a function of observed quantities:

\\E\[Y(1)\] - E\[Y(0)\] = E \left\[ E\[Y \| X, A = 1\] \right\] - E
\left\[ E\[Y \| X, A = 0\] \right\]\\ To estimate \\E \left\[ E \left\[
Y \| X, A = a \right\] \right\] = \theta_a\\, we use g-computation
([Snowden et al.
2011](#ref-snowdenImplementationGComputationSimulated2011)):

\\\hat{\theta}\_a = \frac{1}{n}\sum\_{i=1}^n {\mu(a, x_i)}\\

where \\\mu(a, x_i)\\ is the predicted value from a regression of \\Y\\
on \\A\\ and \\X\\ for a unit with covariate profile \\X=x_i\\ and
treatment \\A\\ set to \\a\\. This estimator is sometimes also known as
the “regression estimator” or the “plug-in estimator”.

Here, we use traditional BART and BCF to model \\\mu(A, X)\\. The rest
of the analysis comes from the definitions above.

### Conditional, Sample, and Population Averages

Because the estimator averages the model’s predictions over the
covariates of the units in the sample, what it estimates is the average
of the conditional effect \\\mu(1, x) - \mu(0, x)\\ over those covariate
values, computed within each posterior draw. Li et al.
([2023](#ref-liBayesianCausalInference2023)) call this the “mixed
average treatment effect”, and note, as Ding and Li
([2018](#ref-dingCausalInferenceMissing2018)) did, that most Bayesian
causal analyses in fact report it.

It differs from the sample average treatment effect and the population
average treatment effect, the two averages usually defined, in what it
holds fixed. The sample average treatment effect is the average of each
unit’s own \\Y_i(1) - Y_i(0)\\, which conditions on the units’ observed
outcomes as well as their covariates. Inference for it requires imputing
each unit’s missing potential outcome, and that imputation depends on
how strongly a unit’s two potential outcomes are associated, about which
the data carry no information. The population average treatment effect
is the average over the population the sample was drawn from, which
treats the covariates as random as well and is the most uncertain of the
three ([Li et al. 2023](#ref-liBayesianCausalInference2023)). The mixed
average needs neither a model for the covariates nor an assumption about
how the potential outcomes are associated, which is why it is the usual
target.
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
estimates it for the ATE, the ATT, and the ATC alike, averaging over the
treated or untreated units’ covariates for the latter two.

### Controlling Sparsity

Before we fit the outcome model, we need to change one setting to make
traditional BART suitable for estimating the ATE. The default splitting
prior is a variable-selection prior. It can drop a predictor from every
tree in the forest at once, which makes it worth having when the goal is
prediction, but not when the estimand is a contrast on one particular
predictor. The traditional BART outcome model below is fitted with
`sparsity = FALSE`, which weights the predictors equally and cannot drop
any of them.

The propensity score model above keeps the default, and should:
predicting who was treated is a prediction problem, and no contrast is
read off it.
[`?bartisan_control`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
explains both halves of this, and using the `split_prior` argument can
be an alternative when there are enough covariates that weighting them
all alike is wasteful.

## The Outcome Model

We can fit the outcome model as a binary BART regression of the outcome
on the treatment and covariates, as in Hill ([2011](#ref-hill2011)),
with the propensity score estimated above added as one more covariate:

``` r

rhc$prop_score <- prop_score

fit <- bartisan(death ~ rhc + age + sex + race + edu + aps + meanbp + resp +
                  hema + pafi + paco2 + crea + surv2m + card + prop_score,
                data = rhc, family = binomial(), sparsity = FALSE)
```

We include the treatment, the covariates, and the propensity score in
the formula’s right-hand side, set `family = binomial()` to model the
binary outcome with logistic regression, and `sparsity = FALSE` to
remove the sparsity-inducing prior. Including the propensity score is
recommended by Carnegie ([2019](#ref-carnegie2019)), though Souto and
Louzada ([2024](#ref-soutoAblationStudiesNovel2024)) report that it
often makes little difference to the estimate. Linero
([2024](#ref-lineroNonparametricHighDimensionalModels2024)) gives the
Bayesian reason for including it: when the outcome model is flexible and
its prior is independent of the model for the treatment, the prior it
implies for the amount of confounding bias is concentrated near zero,
strongly enough that the data may not overcome it, and adjusting for the
propensity score is one way to relax that. Because the score is a
function of the covariates, adding it does not change what the model
conditions on, and the g-computation below holds it at each patient’s
own value.

The fits in this vignette use the default single chain for brevity.
Normally, we would run several chains, possibly change the number of
burn-in and retained draws, and examine convergence diagnostics to make
sure the model was fit correctly; see
[`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md)
for how to set these fit controls and check them.

## Potential Outcomes

The quantities underlying the effect estimate are the two average
potential outcomes: the proportion who would die if every patient were
catheterized, and the proportion who would die if none were.
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
computes them on the way to the effect and prints them beneath it. On a
fit from
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md),
the treatment has to be named, since nothing in the formula marks one
predictor as the treatment:

``` r

ate <- estimate_effect(fit, treat = "rhc")

ate
#> Average treatment effect (difference)
#> 
#> Treatment: `rhc`
#> Averaged over 1500 units
#> 
#>     contrast estimate   lower upper    n
#>  Y[1] - Y[0]   0.0638 0.00742  0.12 1500
#> 
#> Average potential outcomes
#> 
#>  quantity estimate lower upper
#>      Y[0]    0.630 0.599 0.660
#>      Y[1]    0.694 0.654 0.733
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ Y[a] is the average response with `rhc` set to "a".
```

Under the assumptions above this is the average treatment effect:
catheterization raises the probability of death by about 6 percentage
points, with an interval running from roughly 0.7 to 12.

The two rows below the contrast are the estimates of \\E\[Y(0)\]\\ and
\\E\[Y(1)\]\\, each averaged over the observed covariate distribution.
Reporting both is often more informative than reporting their difference
alone, because a difference of a few percentage points means something
different against a baseline of 63% than it would against 5%. Setting
`potential_outcomes = FALSE` in the
[`print()`](https://rdrr.io/r/base/print.html) call drops them where the
difference is all that is wanted.

### Computation and Scale of the Estimate

[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
performs the g-computation above literally. Every unit is predicted
twice, once with the treatment set to each of its levels; the two sets
of predictions are averaged over the units the estimand asks for,
*within each posterior draw*; and only then are the two averages
contrasted. What comes back is a posterior for the estimand, summarized
by its mean and a credible interval.

That order makes a ratio here the marginal ratio rather than the average
of the conditional ones, which is a different quantity. For a difference
the two orders agree, so it only shows up once a ratio is asked for.

The scale matters, and it is the reason the default is
`type = "response"`. On that scale, the average of the unit-level
differences *is* the marginal effect. A contrast read off the link scale
is not: on a logistic fit, the average of the conditional log odds
ratios is not the marginal log odds ratio, and the two can differ by a
good deal. The `comparison` argument asks for the contrast rather than
the scale, so a risk ratio or an odds ratio is available without leaving
the response scale:

``` r

estimate_effect(fit, treat = "rhc", comparison = "lnor")
#> Average treatment effect (log odds ratio)
#> 
#> Treatment: `rhc`
#> Averaged over 1500 units
#> 
#>                contrast estimate  lower upper    n
#>  log(O(Y[1]) / O(Y[0]))    0.287 0.0329 0.544 1500
#> 
#> Average potential outcomes
#> 
#>  quantity estimate lower upper
#>      Y[0]    0.630 0.599 0.660
#>      Y[1]    0.694 0.654 0.733
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ Y[a] is the average response with `rhc` set to "a", and O(y) is the odds
#>   `y/(1-y)`.
```

### The Same Answer Through *marginaleffects*

A fit also works with *marginaleffects*, and for the average effect the
two routes compute the same thing from the same draws. The one thing to
set is the posterior summary: *marginaleffects* reports the median by
default and
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
reports the mean, so an unadjusted comparison shows a difference that is
not there.

``` r

options(marginaleffects_posterior_center = mean)

avg_comparisons(fit, variables = "rhc")
#> 
#>  Estimate   2.5 % 97.5 %
#>    0.0638 0.00742   0.12
#> 
#> Term: rhc
#> Type: response
#> Comparison: 1 - 0
```

The point estimate and the interval match
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
above to machine precision. Which to reach for is a question of what
else is wanted:
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
covers the estimands a treatment question asks for and needs no extra
package, while *marginaleffects* covers a much wider class of
quantities.
[`vignette("effects")`](https://ngreifer.github.io/bartisan/articles/effects.md)
is about the second.

## Using Bayesian Causal Forests

Traditional BART shrinks the fitted function toward a constant; that
shrinkage applies to everything the forest fits, including the part of
the outcome that depends on the exposure. When the exposure is strongly
predicted by the covariates, the forest can explain the outcome using
the covariates alone, leaving little for the exposure to explain, and
shrink the estimated effect toward zero. The interval shrinks with it,
so the result is a confident estimate biased toward no effect.

Hahn et al. ([2020](#ref-hahn2020)) identified this mechanism; their
remedy is to give the treatment effect its own forest with its own
prior, so that shrinking the confounding part does not shrink the
effect. This is the Bayesian causal forest (BCF) model, a special case
of the varying coefficients BART model;
[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) fits it.

One specifies the control function in the model formula and identifies
the treatment in the `treat` argument.
[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) then
fits a varying coefficient BART model, the BCF. By default,
[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md)
estimates a propensity score using a logistic BART model and includes
that as a covariate in the control function, as recommended by Hahn et
al. ([2020](#ref-hahn2020)) and for the reason given above for the plain
BART model.

``` r

fit_bcf <- bcf(death ~ age + sex + race + edu + aps + meanbp + resp + hema +
                 pafi + paco2 + crea + surv2m + card,
               treat = ~ rhc, data = rhc,
               family = binomial())

fit_bcf
#> Generalized BART
#> 
#> Call:
#> bcf(formula = death ~ age + sex + race + edu + aps + meanbp + 
#>     resp + hema + pafi + paco2 + crea + surv2m + card, treat = ~rhc, 
#>     data = rhc, family = binomial())
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 1500
#> Structure: 2 forests of 50 and 25 trees, soft decision rules
#> Draws: 800 kept after 200 warmup
#> 
#> Posterior means: b.rhc.0 = -0.357, b.rhc.1 = -0.396
#> 
#> Treatment: "rhc"
#> Effect moderators: "age", "sex", "race", "edu", "aps", "meanbp", "resp", "hema", "pafi", "paco2", "crea", "surv2m", and "card"
#> ℹ `estimate_effect()` reports the treatment effect, with the average potential
#>   outcomes beside it; `plot()` draws the conditional ones.
```

Using [`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md)
rather than writing the varying coefficient BART model by hand sets four
settings to improve effect estimation:

1.  the effect gets a forest of its own
2.  the estimated propensity score goes into the control function and
    not into the effect forest
3.  the effect forest gets fewer trees because effect heterogeneity is
    usually simpler than a prognostic surface
4.  a binary treatment’s coding is drawn rather than fixed, so the
    answer does not depend on which arm was written as 1.

Note that `sparsity = FALSE` is not among them: the reason for it in the
section above is that the prior can drop the predictor whose contrast we
want, and here the treatment is the coefficient rather than a predictor
the forest splits on, so nothing can drop it.

Because the treatment is named in the call,
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
needs nothing else:

``` r

estimate_effect(fit_bcf)
#> Average treatment effect (difference)
#> 
#> Treatment: `rhc`
#> Averaged over 1500 units
#> 
#>     contrast estimate    lower  upper    n
#>  Y[1] - Y[0]   0.0385 -0.00195 0.0953 1500
#> 
#> Average potential outcomes
#> 
#>  quantity estimate lower upper
#>      Y[0]    0.639 0.607  0.67
#>      Y[1]    0.678 0.641  0.72
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ Y[a] is the average response with `rhc` set to "a".
```

The two potential outcomes are printed beneath the contrast, since a
difference of a few percentage points means one thing against a baseline
of 64% and another against 5%; setting `potential_outcomes = FALSE` in
[`print()`](https://rdrr.io/r/base/print.html) suppresses them. The
contrast’s label names the quantity rather than leaving it to the
heading, which matters once a ratio is asked for: `Y[1] - Y[0]` is a
difference of average responses where `log(O(Y[1]) / O(Y[0]))` is a log
odds ratio, and the note beneath the table says what `Y[a]` is, adding
what `O(y)` is once an odds ratio is asked for.

The output of [`summary()`](https://rdrr.io/r/base/summary.html) on a
[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) fit is
the same summary of the forests it is on any other fit, and says at the
end where the effect is reported. That way the same call means the same
thing whichever way the model was written.

### Conditional Effects

The effect forest gives one coefficient per patient, and
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
turns those coefficients into the conditional effect at each patient’s
covariates with `estimand = "CATE"`, reported on the response scale
rather than on the forest’s own link scale. A conditional effect is the
average effect among patients who share those covariates, not the
individual effect for that patient, which depends on how the patient’s
two potential outcomes are associated and which the data cannot identify
([Li et al. 2023](#ref-liBayesianCausalInference2023)). Here we ask for
them as odds ratios, which for a conditional effect is a conditional
odds ratio:

``` r

cate <- estimate_effect(fit_bcf, estimand = "CATE", comparison = "or")

quantile(cate$estimate, probs = c(0, .25, .5, .75, 1))
#>    0%   25%   50%   75%  100% 
#> 1.092 1.223 1.274 1.326 1.469
```

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
conditional effects, as differences unless the `comparison` argument
says otherwise: patients ordered by their estimate, with a credible
interval each, and the marginal effect as a single interval past the
right edge in its own color. Setting `marginal = FALSE` leaves it out.

``` r

plot(fit_bcf)
```

![](causal_files/figure-html/bcfplot-1.png)

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the fit is
the conditional effects on the response scale; passing a
`<bartisan_effect>` object to
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws whatever
that object holds, so `plot(estimate_effect(fit_bcf, by = ~ sex))` is a
subgroup forest plot instead.

Note that a conditional odds ratio is not the marginal one, and their
average is not it either. That is the scale caution from earlier in a
second place: an average of conditional contrasts equals the marginal
contrast only when the contrast is a difference. For an identity link,
which the next example uses, the ATE *is* the average of the conditional
effects.

## A Second Example: A Continuous Outcome, and the ATT

The catheterization question is about a whole population, so the ATE is
the estimand. Many questions are not. When a program is offered to a
particular group and the question is whether it helped *them*, the ATT
is the estimand to report, and the outcome is often continuous rather
than binary.

`lalonde`, from *cobalt*, is the standard example: a job training
program, with earnings in 1978 as the outcome.

``` r

data("lalonde", package = "cobalt")

with(lalonde, tapply(re78, treat, mean))
#>    0    1 
#> 6984 6349
```

Participants earned less than non-participants. Taken at face value, the
program looks harmful, but that doesn’t mean it actually is:
participants were selected for being out of work, so they would have
earned less anyway. This is confounding in the opposite direction to the
catheterization example, and it is why the raw comparison is worth
showing before the adjusted one.

Below, we fit the BCF model, this time setting `family = tweedie()`,
which is useful for continuous outcomes with a point mass at 0, as we
have here. See
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
for other families that could be used here.

``` r

fit_earn_bcf <- bcf(
  re78 ~ age + educ + race + married + nodegree + re74 + re75,
  treat = ~ treat,
  data = lalonde, family = tweedie()
)

estimate_effect(fit_earn_bcf, estimand = "ATT")
#> Average treatment effect on the treated (difference)
#> 
#> Treatment: `treat`
#> Averaged over the 185 units in group "1"
#> 
#>     contrast estimate lower upper   n
#>  Y[1] - Y[0]     1010  -316  2590 185
#> 
#> Average potential outcomes
#> 
#>  quantity estimate lower upper
#>      Y[0]     5370  4280  6500
#>      Y[1]     6390  5400  7380
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ Y[a] is the average response with `treat` set to "a".
```

Setting `estimand = "ATT"` averages over the treated rather than over
everyone.[^1]

The potential outcomes earn their place here: a difference of about a
thousand dollars means something different against a baseline of about
five thousand than it would against five hundred.

Because the link is the identity, the ATT is exactly the average of the
conditional effects among the treated, which is worth checking once to
see that the two agree:

``` r

att <- estimate_effect(fit_earn_bcf, estimand = "ATT")

cate_att <- estimate_effect(fit_earn_bcf, estimand = "CATE",
                            newdata = subset(lalonde, treat == 1))

c(ATT = att$estimate, mean_CATE = mean(cate_att$estimate))
#>       ATT mean_CATE 
#>      1015      1015
```

And the conditional effects, drawn:

``` r

plot(cate_att)
```

![](causal_files/figure-html/lalondeplot-1.png)

Most of the conditional effects are positive, as the ATT is, and only
two of their 185 intervals exclude zero. That is the usual picture: the
effect at one unit’s covariates is estimated from far less information
than an average over all of them, so its interval is much wider.

## Interpreting the Credible Interval

The posterior interval is a credible interval for the estimand under the
model and under the identification assumptions. It covers uncertainty
about the outcome surface at the covariate values in the sample, which
comes from having a finite sample and from not knowing the surface’s
shape. As described in the section on conditional, sample, and
population averages, it treats those covariates as fixed, so it does not
include the additional uncertainty an interval for the population
average would carry.

For a [`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md)
fit, the interval also takes the propensity score as known. The score is
the posterior mean of a separate model for the treatment, estimated once
before the outcome model is fit, so uncertainty in it does not reach the
interval. Estimating the two models in stages like this is the usual
practice, and it keeps the outcome from informing the propensity score,
which a joint model of the two would allow; Li et al.
([2023](#ref-liBayesianCausalInference2023)) discuss both approaches and
the alternatives between them.

The interval does not cover uncertainty about whether the assumptions
hold.

An unmeasured confounder does not widen the interval. It moves the
estimate and leaves the interval where it was. This is why a sensitivity
analysis, asking how strong a confounder would have to be to explain the
result away, is a more informative addition than any refinement of the
model.

## Further Reading

[`?estimate_effect`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
is the reference for the estimands used here, including the subgroup
form and the contrast types.
[`vignette("effects")`](https://ngreifer.github.io/bartisan/articles/effects.md)
covers what *marginaleffects* adds beyond them, which is a much wider
class of quantities and the route to take when the question is not a
treatment contrast.
[`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md)
covers checking the fit, which should happen before any of this is
interpreted.

## References

Alcantara, Rafael, Meijia Wang, P. Richard Hahn, and Hedibert Lopes.
2024. *Modified BART for Learning Heterogeneous Effects in Regression
Discontinuity Designs*. <https://doi.org/10.48550/arXiv.2407.14365>.

Carnegie, Nicole Bohme. 2019. “Comment: Contributions of Model Features
to BART Causal Inference Performance Using ACIC 2016 Competition Data.”
*Statistical Science* 34 (1): 90–93.
<https://doi.org/10.1214/18-STS682>.

Connors, Alfred F., Theodore Speroff, Neal V. Dawson, et al. 1996. “The
Effectiveness of Right Heart Catheterization in the Initial Care of
Critically Ill Patients.” *JAMA* 276 (11): 889–97.
<https://doi.org/10.1001/jama.1996.03540110043030>.

Ding, Peng, and Fan Li. 2018. “Causal Inference: A Missing Data
Perspective.” *Statistical Science* 33 (2): 214–37.
<https://doi.org/10.1214/18-STS645>.

Dorie, Vincent, Jennifer Hill, Uri Shalit, Marc Scott, and Dan Cervone.
2019. “Automated Versus Do-It-Yourself Methods for Causal Inference:
Lessons Learned from a Data Analysis Competition.” *Statistical Science*
34 (1): 43–68. <https://doi.org/10.1214/18-STS667>.

Hahn, P. Richard, Jared S. Murray, and Carlos M. Carvalho. 2020.
“Bayesian Regression Tree Models for Causal Inference: Regularization,
Confounding, and Heterogeneous Effects (with Discussion).” *Bayesian
Analysis* 15 (3): 965–1056. <https://doi.org/10.1214/19-BA1195>.

Hill, Jennifer L. 2011. “Bayesian Nonparametric Modeling for Causal
Inference.” *Journal of Computational and Graphical Statistics* 20 (1):
217–40. <https://doi.org/10.1198/jcgs.2010.08162>.

Li, Fan, Peng Ding, and Fabrizia Mealli. 2023. “Bayesian Causal
Inference: A Critical Review.” *Philosophical Transactions of the Royal
Society A: Mathematical, Physical and Engineering Sciences* 381 (2247):
20220153. <https://doi.org/10.1098/rsta.2022.0153>.

Linero, Antonio R. 2024. “In Nonparametric and High-Dimensional Models,
Bayesian Ignorability Is an Informative Prior.” *Journal of the American
Statistical Association* 119 (548): 2785–98.
<https://doi.org/10.1080/01621459.2023.2278202>.

McCulloch, Robert E., Rodney A. Sparapani, Brent R. Logan, and
Purushottam W. Laud. n.d. *Causal Inference with the Instrumental
Variable Approach and Bayesian Nonparametric Machine Learning*.
<https://doi.org/10.48550/arXiv.2102.01199>.

Snowden, Jonathan M., Sherri Rose, and Kathleen M. Mortimer. 2011.
“Implementation of g-Computation on a Simulated Data Set: Demonstration
of a Causal Inference Technique.” *American Journal of Epidemiology* 173
(7): 731–38. <https://doi.org/10.1093/aje/kwq472>.

Souto, Hugo Gobato, and Francisco Louzada. 2024. *Ablation Studies for
Novel Treatment Effect Estimation Models*.
<https://doi.org/10.48550/arXiv.2410.15560>.

[^1]: The `focal` argument is not needed because a 0/1 treatment settles
    which level is the treated one. When the treatment has more than two
    levels, `focal` is required and naming it is the whole of the
    choice, since `"ATT"` and `"ATC"` then mean the same thing: the
    effect among the units in the level named. When it has two levels
    whose labels say nothing about which is which, the level order is
    assumed and a message says so, which `focal` silences.
