# Frequently Asked Questions

## Introduction

These are the questions that come up most often about BART and about
*bartisan*. Each answer stands on its own, so they can be read in any
order, and each one points at the vignette or help page that works its
topic through in full.

## What is BART?

BART stands for “Bayesian Additive Regression Trees”. It’s a Bayesian
machine learning method used to model the relationship between
predictors and outcomes. It’s similar to GBM/xgboost, extraTrees, and
random forests in that the fit is built out of many small trees rather
than one big one. BART sums their contributions, as boosting does, where
random forests and extraTrees average theirs. This lets it approximate
complex and nonlinear functions without you needing to specify their
form.

## How/When should I use BART?

Use BART exactly the same way you would use any regression method or any
machine learning method.

## How can a machine learning method be Bayesian?

Each parameter of the model, including which variables are used to split
the trees, how deep each tree is, the predicted value in the tree, the
residual variance, etc., are assigned a prior, and through running the
model on the data, each model output gains a posterior. In traditional
machine learning models, these parameters are given an individual value,
or the values are tuned across a grid. The posterior enables Bayesian
inference on BART’s predictions.

## How can a Bayesian method be machine learning?

Instead of specifying a specific functional form as you would in a
Bayesian generalized linear model (GLM), you specify priors governing a
tree structure (e.g., the depth of trees, the variables to split on),
and that structure is used to model the regression function. It’s still
Bayesian because all components of the model have a prior and posterior,
but rather than putting priors on GLM coefficients, you put priors on
aspects of the machine learning model.

## How can BART be used for inference?

Because each prediction from the model has a posterior, any quantity
derived from the predictions, like average marginal effects, has one
too. These can be used as part of a modern regression workflow through
tools like *marginaleffects* ([Arel-Bundock et al.
2024](#ref-arel-bundockHowInterpretStatistical2024)), which produce
interpretable estimates from regression models in model-agnostic ways
([Rohrer and Arel-Bundock
2026](#ref-rohrerModelsPredictionMachines2026)). These tools make the
black-box nature of machine learning models no longer problematic for
arriving at interpretable model summaries, and the posteriors of these
quantities allow for Bayesian inference on them.

## How do I get estimates out of the model?

[`predict()`](https://rdrr.io/r/stats/predict.html) gives you the fitted
values or the full posterior of them.
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
turns those into a treatment effect: the ATE, the ATT, or one effect per
unit, with a credible interval and on the response scale, which is the
scale on which averaging unit-level contrasts gives the marginal effect.
[`partial_dependence()`](https://ngreifer.github.io/bartisan/reference/partial_dependence.md)
shows how the prediction moves with one or two predictors, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a fit draws
it.
[`variable_importance()`](https://ngreifer.github.io/bartisan/reference/variable_importance.md)
reports how often each predictor is split on, which is a rough guide to
what the model is using but not a test of anything.

For the quantities those do not cover,
[`marginaleffects::avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
and its relatives work on a fit directly, including arbitrary contrasts
between covariate values, hypotheses comparing one estimate to another,
and slopes, which want `x_transform = "range"` in
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md),
since the default quantile transform leaves the fit a step function of
the original predictor. Because every draw goes through the same
machinery either way, the result carries a posterior rather than a point
estimate with a delta-method standard error. See
[`vignette("effects")`](https://ngreifer.github.io/bartisan/articles/effects.md)
for the worked versions and
[`?estimate_effect`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
for the estimands.

## What if I have missing data?

Missing values in the predictors are fine and nothing is dropped for
them, though [`predict()`](https://rdrr.io/r/stats/predict.html) accepts
a missing value only in a column that had one at fitting time, since
only those columns’ rules carry an answer for it. Each splitting rule
carries its own answer for what to do with a missing value, drawn from
the prior alongside the variable and the cutpoint, so the missingness is
part of the model rather than something to impute first. This is
missingness incorporated in attributes, and it means a covariate can be
informative through whether it is observed as well as through its value.
A missing *outcome* is different: those rows are dropped, with a warning
saying how many.

## Can I supply my own likelihood?

Yes, through
[`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
which is the feature the package is named for. You give it a function of
the response and the additive predictors, `function(y, eta)`, returning
one log density per observation, and the sampler does the rest; you do
not have to derive a conjugate update or write any C++. The trade-off is
that a log density says nothing about how to draw from the distribution
it describes, so a custom family has no posterior predictive draws, and
everything built on them, including
[`pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)
and `posterior_predict()`, is unavailable. Fitting, prediction on the
link scale, and [`loo()`](https://mc-stan.org/loo/reference/loo.html)
all work, though `predict(type = "response")` hands back the additive
predictors rather than a fitted mean, since the package cannot know
where the mean of a supplied density sits.

## Why is my model slow, and what should I change first?

The three settings that matter most, in the order worth trying them.
`num_trees` is close to linear in cost, and a second forest can usually
be smaller than the first, so `num_trees = c(mean = 50, log_sd = 10)` on
a location-scale model is much cheaper than 50 for both. `gate = "hard"`
is several times faster than the soft default and costs nothing if the
relationship is well approximated by a step function. `bandwidth_every`
controls how often the soft-rule bandwidth is resampled, and raising it
keeps soft rules while recovering some of that speed; it is not free,
though, since the bandwidth update is what lets the rules sharpen toward
a step, so it costs some mixing and some accuracy on a mean function
with jumps. `update_bandwidth = FALSE` stops resampling it altogether,
which is faster again and is *more* accurate on smooth functions and
much worse on nonsmooth ones. Beyond those, run the chains in parallel
with a *future* plan, which costs nothing in draws.

## What if I’m a frequentist?

You don’t have to use the full posterior for inference; you can use the
predictions from a BART model (e.g., the posterior mean for each
prediction) in a frequentist estimator like double/debiased machine
learning (DML) ([Ahrens et al.,
n.d.](#ref-ahrensIntroductionDoubleDebiased2026); [Chernozhukov et al.
2018](#ref-chernozhukov2018)) or targeted maximum likelihood estimation
(TMLE) ([Schuler and Rose
2017](#ref-schulerTargetedMaximumLikelihood2017); [Laan and Rubin
2006](#ref-vanderlaanTargetedMaximumLikelihood2006)). These methods only
require predictions from models and allow you to perform valid
frequentist inference on specific estimands. That said, some also treat
the resulting credible intervals as confidence intervals, but the
frequentist operating characteristics of these intervals are not
guaranteed.

## Is BART good?

Yes! BART or a related method won the American Causal Inference
Conference Data Competition in both of the years whose results have been
written up ([Dorie et al.
2019](#ref-dorieAutomatedDoityourselfMethods2019); [Hahn et al.
2019](#ref-hahnAtlanticCausalInference2019)). This is a competition to
see which method can estimate treatment effects as accurately as
possible (including with accurate inference) across a wide variety of
data-generating processes. For general prediction, BART has been shown
to do as well or better than popular methods like GBM and random forests
([Chipman et al. 2010](#ref-chipmanBARTBayesianAdditive2010)).

## How can BART be used for causal inference?

The only step BART can be used for in causal inference is causal effect
estimation, i.e., estimating causal effects given causal assumptions on
the data. It works like any other regression method as part of a causal
effect estimator. The award-winning use of BART often involves
BART-based g-computation ([Hill 2011](#ref-hill2011); [Snowden et al.
2011](#ref-snowdenImplementationGComputationSimulated2011); [Dorie et
al. 2019](#ref-dorieAutomatedDoityourselfMethods2019); [Carnegie
2019](#ref-carnegieCommentContributionsModel2019)), which involves
modeling the relationship between the outcome, treatment, and
confounders, and using that model to predict the counterfactual outcomes
under each treatment for each unit. See
[`vignette("causal")`](https://ngreifer.github.io/bartisan/articles/causal.md)
for a walkthrough with *bartisan*. BART can also be used to estimate
propensity scores ([Hill et al.
2011](#ref-hillChallengesPropensityScore2011)), which, together with a
BART outcome model, can be used in doubly robust estimators like DML and
TMLE. BART can also be used in instrumental variables analysis
([McCulloch et al.,
n.d.](#ref-mccullochCausalInferenceInstrumental2021)), regression
discontinuity ([Alcantara et al.,
n.d.](#ref-alcantaraModifiedBARTLearning2024)), and
difference-in-differences ([Souto and Neto,
n.d.](#ref-soutoForestsDifferencesRobust2025)).

## What are Bayesian Causal Forests?

Bayesian Causal Forests (BCF) are a modification of BART used to
estimate treatment effects ([Hahn et al. 2020](#ref-hahn2020)). In
addition to flexibly modeling the relationship between the outcome and
the covariates, it also flexibly models the relationship between the
*magnitude of the treatment effect* and the covariates. In this way, BCF
is a varying-coefficient model ([Deshpande et al.
2026](#ref-deshpandeVCBARTBayesianTrees2026)). Compared to traditional
BART, BCF tends to have better calibrated intervals and regularization
of heterogeneous treatment effects. BCF tends to outperform traditional
BART in most causal inference contexts. It also tends to outperform
generalized random forests (GRF), which serve a similar function.

## How do I get a treatment effect?

[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
reports it, with the two average potential outcomes printed beneath. A
fit from [`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md)
carries its own treatment, so nothing else has to be named:

``` r

fit <- bcf(y ~ x1 + x2, treat = ~ z, data = d)

estimate_effect(fit)
```

The estimand is an argument, so `estimand = "ATT"` averages over the
treated instead of over everyone, `estimand = "CATE"` returns one effect
per unit, `by = ~ g` gives subgroup effects, and `comparison = "ratio"`
reports a ratio rather than a difference. On a fit from
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
rather than
[`bcf()`](https://ngreifer.github.io/bartisan/reference/bcf.md) it works
the same way once `treat` names the column.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws whatever
was asked for, and
[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
on the result says whether the estimand itself has mixed, which is not
implied by the fit’s own diagnostics.

[`vignette("causal")`](https://ngreifer.github.io/bartisan/articles/causal.md)
is the worked version, including the assumptions that turn the estimate
into a causal effect, which no model supplies.

## What papers should I read to better understand BART?

BART has a growing literature spread across multiple fields. We
recommend reading the original BART paper by Chipman et al.
([2010](#ref-chipmanBARTBayesianAdditive2010)), the paper introducing
BART for causal inference by Hill ([2011](#ref-hill2011)), an accessible
paper on BART for social scientists Green and Kern
([2012](#ref-greenModelingHeterogeneousTreatment2012)), and the Annual
Reviews paper on BART by Hill et al.
([2020](#ref-hillBayesianAdditiveRegression2020)).

## What does the name “bartisan” mean?

*bartisan* is a portmanteau of BART and “artisan”. *bartisan* allows you
to supply your own likelihood to fit the BART model of your choice
without needing to derive or program the whole machinery to do so. In
this way, you can build an artisanal model for your specific scenario.
In addition, the package offers a high degree of flexibility and
customization for its built-in models, including hard and smooth gates
for the trees, random effects, varying coefficient models, many model
families, and control of sparsity.

## How does *bartisan* differ from other BART implementations in R?

There are many R packages that implement BART, each of which has its own
strengths. *bartisan* aims to be highly general but with intelligent
defaults that make it easy to simply plug in in place of a GLM or
generalized additive model (GAM) as part of a data analysis workflow.
*bartisan* offers many more model types and families than other
packages. Often, R packages will be built to implement a few specific
features of BART; for example,
[*SoftBart*](https://cran.r-project.org/package=SoftBart) implements
soft decision trees,
[*VCBART*](https://cran.r-project.org/package=VCBART) implements the
varying coefficient model,
[*dbarts*](https://cran.r-project.org/package=dbarts) implements Normal
and probit BART regression with optional random intercepts.
[`vignette("implementation")`](https://ngreifer.github.io/bartisan/articles/implementation.md)
has a feature-by-feature comparison against the most popular BART
packages, along with the relative cost of each family. *bartisan* was
designed to incorporate as many recent advancements in BART theory and
practice as possible, including models for different outcome types
(e.g., multinomial, ordinal, \\\left\[ 0, 1 \right\]\\-bounded,
zero-inflated), soft trees, sparsity priors, varying coefficients,
random effects, and more, all in a single package[^1]. There are ways in
which *bartisan* is inferior to these more specialized packages, but as
a general-purpose tool, I hope you’ll find it effective.

## What model family should I use?

This of course depends on the outcome and the assumptions you are
willing to make about it. In general, the defaults are a good starting
point, though for numeric outcomes, bounding (e.g., at 0 or 1, inclusive
or exclusive, or integer or continuous) might suggest different
families. Check out
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
and `?bartisan-families` for more details and help with this decision.

## How should I select the values of the hyperparameters?

One of the great advantages of BART is that performance tends to be
strong using the default hyperparameters, and tuning them with
cross-validation tends to buy very little for the increase in
computation. We recommend changing the default parameters only if you
understand very well what each one does.

That said, there are a few hyperparameters that might be considered.
`gate` controls whether hard or soft decision trees are used; the
default, soft decision trees, tend to perform better when modeling
smooth relationships, but they cost substantially more computation,
because every observation reaches every leaf with some weight rather
than taking one side of each split. `gate = "hard"` is the fast option
when the relationship is expected to be a step function anyway.
`sparsity` controls how much sparsity is induced in the covariates
chosen for splitting; the default is to have sparsity on, which is
especially useful with many predictors, but sometimes it can be valuable
to turn the sparsity prior off to ensure all variables are used in the
model. `split_prior` can also be used instead of `sparsity` to manually
decide which variables should be split on more often, though in general
the Bayesian updates can discover this from the data anyway. See
[`?bartisan_control`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
for a full list.

## How should I choose how many chains and burn-in and posterior draws?

For all of these parameters, more is always better but requires more
computation and more memory usage. It can be a good idea to build your
model on synthetic data with the default number of burn-in and posterior
draws, and then to use a much higher number for the final analysis on
the real data. More chains come almost for free in wall-clock terms if
you have parallel processing, as each chain can be run in parallel,
whereas burn-in and posterior draws must be run sequentially. They are
not free diagnostically, though, and it is worth knowing which way that
cuts. Effective sample size depends on the total number of draws you
keep and not on how they are split up, so four chains of 10,000 and
eight of 5,000 give you about the same precision. But R-hat compares
chains against each other, and the value it takes under a perfectly
behaved sampler rises with the number of chains being compared: it sits
near \\1 + m/S\\ for \\m\\ chains and effective sample size \\S\\.
Splitting a fixed budget into more chains therefore raises the bar R-hat
has to clear without buying any more information, and a fit that passes
at four chains can warn at twelve. Two to four chains is a good default,
and when a quantity fails
[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)’s
R-hat threshold on too few effective draws for the comparison to mean
anything, the report adds a line saying the number cannot be read yet
rather than letting it stand on its own. See
[`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md)
for more information on diagnosing and fixing poor mixing.

## References

Ahrens, Achim, Victor Chernozhukov, Christian Hansen, Damian Kozbur,
Mark Schaffer, and Thomas Wiemann. n.d. *An Introduction to
Double/Debiased Machine Learning*.
<https://doi.org/10.48550/arXiv.2504.08324>.

Alcantara, Rafael, Meijia Wang, P. Richard Hahn, and Hedibert Lopes.
n.d. *Modified BART for Learning Heterogeneous Effects in Regression
Discontinuity Designs*. <https://doi.org/10.48550/arXiv.2407.14365>.

Arel-Bundock, Vincent, Noah Greifer, and Andrew Heiss. 2024. “How to
Interpret Statistical Models Using Marginaleffects for r and Python.”
*Journal of Statistical Software* 111 (November): 1–32.
<https://doi.org/10.18637/jss.v111.i09>.

Carnegie, Nicole Bohme. 2019. “Comment: Contributions of Model Features
to BART Causal Inference Performance Using ACIC 2016 Competition Data.”
*Statistical Science* 34 (1): 90–93.
<https://doi.org/10.1214/18-STS682>.

Chernozhukov, Victor, Denis Chetverikov, Mert Demirer, et al. 2018.
“Double/Debiased Machine Learning for Treatment and Structural
Parameters.” *The Econometrics Journal* 21 (1): C1–68.
<https://doi.org/10.1111/ectj.12097>.

Chipman, Hugh A., Edward I. George, and Robert E. McCulloch. 2010.
“BART: Bayesian Additive Regression Trees.” *The Annals of Applied
Statistics* 4 (1): 266–98. <https://doi.org/10.1214/09-AOAS285>.

Deshpande, Sameer K., Ray Bai, Cecilia Balocchi, Jennifer E. Starling,
and Jordan Weiss. 2026. “VCBART: Bayesian Trees for Varying
Coefficients.” *Bayesian Analysis* 21 (1): 281–308.
<https://doi.org/10.1214/24-BA1470>.

Dorie, Vincent, Jennifer Hill, Uri Shalit, Marc Scott, and Dan Cervone.
2019. “Automated Versus Do-It-Yourself Methods for Causal Inference:
Lessons Learned from a Data Analysis Competition.” *Statistical Science*
34 (1): 43–68. <https://doi.org/10.1214/18-STS667>.

Green, Donald P., and Holger L. Kern. 2012. “Modeling Heterogeneous
Treatment Effects in Survey Experiments with Bayesian Additive
Regression Trees.” *Public Opinion Quarterly* 76 (3): 491–511.
<https://doi.org/10.1093/poq/nfs036>.

Hahn, P. Richard, Vincent Dorie, and Jared S. Murray. 2019. “Atlantic
Causal Inference Conference (ACIC) Data Analysis Challenge 2017.”
*arXiv:1905.09515 \[Stat\]*, May. <http://arxiv.org/abs/1905.09515>.

Hahn, P. Richard, Jared S. Murray, and Carlos M. Carvalho. 2020.
“Bayesian Regression Tree Models for Causal Inference: Regularization,
Confounding, and Heterogeneous Effects (with Discussion).” *Bayesian
Analysis* 15 (3): 965–1056. <https://doi.org/10.1214/19-BA1195>.

Hill, Jennifer L. 2011. “Bayesian Nonparametric Modeling for Causal
Inference.” *Journal of Computational and Graphical Statistics* 20 (1):
217–40. <https://doi.org/10.1198/jcgs.2010.08162>.

Hill, Jennifer, Antonio Linero, and Jared Murray. 2020. “Bayesian
Additive Regression Trees: A Review and Look Forward.” *Annual Review of
Statistics and Its Application* 7 (1): annurev-statistics-031219-041110.
<https://doi.org/10.1146/annurev-statistics-031219-041110>.

Hill, Jennifer, Christopher Weiss, and Fuhua Zhai. 2011. “Challenges
with Propensity Score Strategies in a High-Dimensional Setting and a
Potential Alternative.” *Multivariate Behavioral Research* 46 (3):
477–513. <https://doi.org/10.1080/00273171.2011.570161>.

Laan, Mark J. van der, and Daniel Rubin. 2006. “Targeted Maximum
Likelihood Learning.” *The International Journal of Biostatistics* 2
(1). <https://doi.org/10.2202/1557-4679.1043>.

McCulloch, Robert E., Rodney A. Sparapani, Brent R. Logan, and
Purushottam W. Laud. n.d. *Causal Inference with the Instrumental
Variable Approach and Bayesian Nonparametric Machine Learning*.
<https://doi.org/10.48550/arXiv.2102.01199>.

Rohrer, Julia M., and Vincent Arel-Bundock. 2026. “Models as Prediction
Machines: How to Convert Confusing Coefficients into Clear Quantities.”
*Advances in Methods and Practices in Psychological Science* 9 (2):
25152459261424825. <https://doi.org/10.1177/25152459261424825>.

Schuler, Megan S., and Sherri Rose. 2017. “Targeted Maximum Likelihood
Estimation for Causal Inference in Observational Studies.” *American
Journal of Epidemiology* 185 (1): 65–73.
<https://doi.org/10.1093/aje/kww165>.

Snowden, Jonathan M., Sherri Rose, and Kathleen M. Mortimer. 2011.
“Implementation of g-Computation on a Simulated Data Set: Demonstration
of a Causal Inference Technique.” *American Journal of Epidemiology* 173
(7): 731–38. <https://doi.org/10.1093/aje/kwq472>.

Souto, Hugo Gobato, and Francisco Louzada Neto. n.d. *Forests for
Differences: Robust Causal Inference Beyond Parametric DiD*.
<https://doi.org/10.48550/arXiv.2505.09706>.

[^1]: I am well aware of [this xkcd comic](https://xkcd.com/927/).
