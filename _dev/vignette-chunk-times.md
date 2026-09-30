# Vignette chunk times

Measured 2026-09-24 with `_dev/vignette-chunk-timing.R`: each vignette rendered in a fresh R process against the installed package, with `_R_CHECK_LIMIT_CORES_=TRUE`. Seconds are wall-clock per chunk; *overhead* is the render time outside the chunks (R startup, loading, pandoc).

**These times were taken on a loaded machine** (load average about 3.9, with Spotlight's knowledge updater at 120% CPU). Chunks whose code did not change since the quiet run earlier the same day ran a median 1.44 times slower, rising from 1.13 in the first vignette to about 1.6 in the later ones. Dividing each vignette by its own factor estimates 7.6 minutes in total on a quiet machine, against 9.3 before the cuts. Read the ranking, not the absolute numbers, and re-time on a quiet machine before quoting any of them.

## Totals

| Vignette | Render (s) | In chunks (s) | Overhead (s) | Chunks |
| --- | ---: | ---: | ---: | ---: |
| diagnostics | 217.4 | 215.9 | 1.5 | 21 |
| effects | 102.3 | 101.5 | 0.8 | 14 |
| comparison | 101.9 | 100.6 | 1.3 | 17 |
| varying | 77.4 | 76.7 | 0.7 | 20 |
| bartisan | 63.3 | 61.7 | 1.6 | 19 |
| causal | 61.4 | 60.1 | 1.3 | 19 |
| importance | 13.9 | 13.2 | 0.7 | 10 |
| families | 10.9 | 9.8 | 1.1 | 20 |
| survival | 5.8 | 4.8 | 1.0 | 23 |
| implementation | 3.9 | 3.1 | 0.8 | 7 |
| faq | 0.8 | 0.0 | 0.8 | 0 |
| **All** | **659.0** (11.0 min) | 647.4 | 11.6 | 170 |

## Chain-count variants

| Variant | Render (s) | Current file (s) |
| --- | ---: | ---: |

## Every chunk, in document order

### diagnostics (217.4 s)

| # | Chunk | Evaluated | Seconds | First line |
| ---: | --- | --- | ---: | --- |
| 1 | `setup` | yes | 0.5 | `knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6.5` |
| 2 | `libs` | yes | 9.2 | `library(bartisan)` |
| 3 | `diagnose` | yes | 2.3 | `diagnose(fit)` |
| 4 | `friedmanbad` | yes | 0.4 | `set.seed(3)` |
| 5 | `friedmangood` | yes | 2.7 | `long_enough <- bartisan(y ~ ., fr, family = gaussian(), chains = 4)` |
| 6 | `perobs` | yes | 1.6 | `per_observation_rhat <- function(fit) {` |
| 7 | `summarise` | yes | 0.1 | `library(posterior)` |
| 8 | `trace` | yes | 0.3 | `library(bayesplot)` |
| 9 | `diagnose_effect` | yes | 43.4 | `bcf_fit <- bcf(death ~ age + sex + race + edu + aps + meanbp + resp +` |
| 10 | `diagnose_bcf_fit` | yes | 3.9 | `diagnose(bcf_fit)$table` |
| 11 | `me_draws` | yes | 11.6 | `library(marginaleffects)` |
| 12 | `mixmodel` | yes | 11.3 | `model <- death ~ rhc + age + aps + meanbp + surv2m` |
| 13 | `mixfixed` | yes | 117.2 | `set.seed(2026)` |
| 14 | `priorsummary` | yes | 0.1 | `rstantools::prior_summary(fit)` |
| 15 | `prioronly` | yes | 4.5 | `set.seed(2026)` |
| 16 | `priorquant` | yes | 0.3 | `p <- c(.01, .1, .5, .9, .99)` |
| 17 | `ppc` | yes | 0.3 | `pp_check(fit)` |
| 18 | `ppcstat` | yes | 0.1 | `pp_check(fit, type = "stat", stat = "sd")` |
| 19 | `ppcloo` | yes | 2.7 | `pp_check(fit, type = "loo_pit_ecdf")` |
| 20 | `calibration` | yes | 3.3 | `pp_check(fit, type = "loo_calibration")` |
| 21 | `r2` | yes | 0.2 | `performance::r2(fit)` |

### effects (102.3 s)

| # | Chunk | Evaluated | Seconds | First line |
| ---: | --- | --- | ---: | --- |
| 1 | `setup` | yes | 0.1 | `knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6.5` |
| 2 | `libs` | yes | 3.6 | `library(bartisan)` |
| 3 | `avgcomp` | yes | 9.4 | `comp <- avg_comparisons(fit, variables = c("rhc", "age", "card"))` |
| 4 | `step` | yes | 3.3 | `avg_comparisons(fit, variables = list(aps = 10))` |
| 5 | `bycard` | yes | 3.3 | `avg_comparisons(fit, variables = "rhc", by = "card")` |
| 6 | `hypo` | yes | 3.2 | `avg_comparisons(fit, variables = "rhc", by = "card",` |
| 7 | `bypred` | yes | 1.7 | `avg_predictions(fit, by = "card")` |
| 8 | `unnamed-chunk-1` | yes | 4.2 | `avg_predictions(fit, variables = "card")` |
| 9 | `pdp` | yes | 7.0 | `library(ggplot2)` |
| 10 | `pdp2` | yes | 20.3 | `pd2 <- partial_dependence(fit, ~ aps + rhc)` |
| 11 | `pdp3` | yes | 19.1 | `plot(fit, ~ meanbp + aps) +` |
| 12 | `pdp4` | yes | 19.4 | `plot(fit, ~ meanbp + aps,` |
| 13 | `pdp5` | yes | 6.9 | `at <- quantile(rhc$aps, c(.1, .5, .9))` |
| 14 | `pdp6` | yes | <0.1 | `plot_predictions(fit, condition = list(aps = at), draw = FALSE)` |

### comparison (101.9 s)

| # | Chunk | Evaluated | Seconds | First line |
| ---: | --- | --- | ---: | --- |
| 1 | `setup` | yes | 2.1 | `knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6.5` |
| 2 | `libs` | yes | 3.0 | `library(bartisan)` |
| 3 | `loo` | yes | 1.0 | `loo(full)` |
| 4 | `holdout` | yes | 2.6 | `set.seed(2026)` |
| 5 | `holdoutscale` | yes | 0.8 | `elpd_total <- loo(full)$estimates["elpd_loo", "Estimate"]` |
| 6 | `compare` | yes | 4.5 | `set.seed(2026)` |
| 7 | `holdoutcompare` | yes | 2.4 | `demographics_train <- bartisan(death ~ rhc + age + sex + race + edu,` |
| 8 | `kfold` | yes | 14.5 | `set.seed(2026)` |
| 9 | `kfoldcompare` | yes | 13.4 | `kfold_demographics <- kfold(demographics, folds = folds)` |
| 10 | `kfoldcheck` | yes | <0.1 | `c(kfold = kfold_full$estimates["elpd_kfold", "Estimate"] / nrow(rhc),` |
| 11 | `links` | yes | 4.4 | `set.seed(2026)` |
| 12 | `survfits` | yes | 22.1 | `library(survival)` |
| 13 | `stanglm` | yes | 3.9 | `logistic <- rstanarm::stan_glm(model, data = rhc, family = binomial(),` |
| 14 | `stancompare` | yes | 1.8 | `loo_compare(list(bart = loo(full), logistic = loo(logistic)))` |
| 15 | `tunegrid` | yes | 13.9 | `trees <- c(20, 50, 200)` |
| 16 | `tunedheld` | yes | 1.7 | `sapply(tuned, function(f) sum(predict(f, newdata = held, type = "densi` |
| 17 | `dropone` | yes | 8.3 | `drop_one <- function(v) {` |

### varying (77.4 s)

| # | Chunk | Evaluated | Seconds | First line |
| ---: | --- | --- | ---: | --- |
| 1 | `setup` | yes | 0.1 | `knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6.5` |
| 2 | `libs` | yes | <0.1 | `library(bartisan)` |
| 3 | `vc` | yes | 17.4 | `set.seed(2026)` |
| 4 | `unnamed-chunk-1` | yes | <0.1 | `summary(fit_vc)` |
| 5 | `vcmod` | yes | 16.3 | `set.seed(2026)` |
| 6 | `unnamed-chunk-2` | yes | <0.1 | `summary(fit_aps)` |
| 7 | `vcself` | no | <0.1 | `bartisan(death ~ age + sex + vc(aps, ~ . + aps),` |
| 8 | `vcconstant` | yes | 12.9 | `set.seed(2026)` |
| 9 | `vccoef` | yes | <0.1 | `b <- coef(constant, draws = TRUE)[["rhc"]][, 1]` |
| 10 | `coef` | yes | <0.1 | `head(coef(fit_vc))` |
| 11 | `coefdraws` | yes | <0.1 | `draws <- coef(fit_vc, draws = TRUE)[["rhc"]]` |
| 12 | `cate` | yes | 7.8 | `cate <- estimate_effect(fit_vc, treat = "rhc", estimand = "CATE")` |
| 13 | `trees` | no | <0.1 | `bartisan(death ~ age + sex + race + edu + aps +` |
| 14 | `vcmulti` | no | <0.1 | `bartisan(list(mean = y ~ x1 + x2 + vc(z),` |
| 15 | `vcmultitrees` | no | <0.1 | `bartisan(y ~ x1 + x2 + vc(z), data = d, family = gaussian_ls(),` |
| 16 | `bcf` | yes | 14.8 | `set.seed(2026)` |
| 17 | `bcfhand` | no | <0.1 | `bartisan(death ~ age + sex + race + edu + aps +` |
| 18 | `bcfplot` | yes | 5.7 | `plot(fit_bcf)` |
| 19 | `resim` | yes | 1.6 | `set.seed(2026)` |
| 20 | `ranef` | yes | <0.1 | `quantile(fit_re$tau[[1L]][, "g"], c(.025, .5, .975))` |

### bartisan (63.3 s)

| # | Chunk | Evaluated | Seconds | First line |
| ---: | --- | --- | ---: | --- |
| 1 | `setup` | yes | 0.4 | `knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6.5` |
| 2 | `libs` | yes | <0.1 | `library(bartisan)` |
| 3 | `unnamed-chunk-1` | no | <0.1 | `future::plan("multisession")` |
| 4 | `unnamed-chunk-2` | no | <0.1 | `progressr::handlers(global = TRUE)` |
| 5 | `unnamed-chunk-3` | no | <0.1 | `progressr::with_progress(` |
| 6 | `data` | yes | <0.1 | `data(rhc)` |
| 7 | `fit` | yes | 2.5 | `set.seed(2026)` |
| 8 | `diagnose` | yes | 1.0 | `diagnose(fit)` |
| 9 | `diagnose_effect` | yes | 2.3 | `diagnose(estimate_effect(fit, treat = "rhc"))` |
| 10 | `ppcheck` | yes | 1.0 | `bayesplot::pp_check(fit, type = "loo_calibration")` |
| 11 | `varimp` | yes | <0.1 | `variable_importance(fit)` |
| 12 | `comparisons` | yes | 34.4 | `marginaleffects::avg_comparisons(fit)` |
| 13 | `comparisons_iqr` | yes | 2.3 | `marginaleffects::avg_comparisons(fit, variables = list(surv2m = "iqr")` |
| 14 | `pdp` | yes | 10.7 | `plot(fit, ~ surv2m) +` |
| 15 | `effects` | yes | 2.8 | `eff <- estimate_effect(fit, treat = "rhc")` |
| 16 | `pred` | yes | <0.1 | `new_patient <- rhc[1, ]` |
| 17 | `predint` | yes | <0.1 | `marginaleffects::predictions(fit, newdata = new_patient)` |
| 18 | `loo` | yes | 0.6 | `library(loo)` |
| 19 | `loocompare` | yes | 3.7 | `set.seed(2026)` |

### causal (61.4 s)

| # | Chunk | Evaluated | Seconds | First line |
| ---: | --- | --- | ---: | --- |
| 1 | `setup` | yes | 0.4 | `knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6.5` |
| 2 | `libs` | yes | 0.1 | `library(bartisan)` |
| 3 | `libs2` | yes | <0.1 | `if (has_cobalt) library(cobalt)` |
| 4 | `baltab` | yes | 0.1 | `bal.tab(rhc ~ age + sex + race + edu + aps + meanbp + resp +` |
| 5 | `naive` | yes | <0.1 | `with(rhc, tapply(death, rhc, mean))` |
| 6 | `psfit` | yes | 2.8 | `ps_fit <- bartisan(` |
| 7 | `balplot` | yes | 0.2 | `bal.plot(rhc ~ prop_score, data = rhc, type = "hist", mirror = TRUE)` |
| 8 | `po` | yes | 2.5 | `fit <- bartisan(death ~ rhc + age + sex + race + edu + aps + meanbp + ` |
| 9 | `ate` | yes | 2.6 | `ate <- estimate_effect(fit, treat = "rhc")` |
| 10 | `ate_or` | yes | 2.4 | `estimate_effect(fit, treat = "rhc", comparison = "lnor")` |
| 11 | `me_equiv` | yes | 2.9 | `options(marginaleffects_posterior_center = mean)` |
| 12 | `bcffit` | yes | 12.4 | `fit_bcf <- bcf(death ~ age + sex + race + edu + aps + meanbp + resp + ` |
| 13 | `bcfeff` | yes | 4.5 | `estimate_effect(fit_bcf)` |
| 14 | `bcfcate` | yes | 4.9 | `cate <- estimate_effect(fit_bcf, estimand = "CATE", comparison = "or")` |
| 15 | `bcfplot` | yes | 4.9 | `plot(fit_bcf)` |
| 16 | `lalonde` | yes | 0.1 | `data("lalonde", package = "cobalt")` |
| 17 | `lalondefit` | yes | 16.9 | `fit_earn_bcf <- bcf(` |
| 18 | `lalondecate` | yes | 2.5 | `att <- estimate_effect(fit_earn_bcf, estimand = "ATT")` |
| 19 | `lalondeplot` | yes | 0.1 | `plot(cate_att)` |

### importance (13.9 s)

| # | Chunk | Evaluated | Seconds | First line |
| ---: | --- | --- | ---: | --- |
| 1 | `setup` | yes | 0.1 | `knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6.5` |
| 2 | `libs` | yes | <0.1 | `library(bartisan)` |
| 3 | `vi` | yes | 4.0 | `set.seed(2026)` |
| 4 | `viplot` | yes | 0.5 | `plot(imp)` |
| 5 | `vidraws` | yes | <0.1 | `counts <- variable_importance(fit, draws = TRUE)` |
| 6 | `noise` | yes | 3.8 | `set.seed(11)` |
| 7 | `friedman` | yes | 1.4 | `friedman <- function(n) {` |
| 8 | `corr` | yes | 0.9 | `set.seed(5)` |
| 9 | `correff` | yes | 1.8 | `library(marginaleffects)` |
| 10 | `corrjoint` | yes | 0.7 | `lo <- transform(cd, x1 = 0.25, x1_copy = 0.25)` |

### families (10.9 s)

| # | Chunk | Evaluated | Seconds | First line |
| ---: | --- | --- | ---: | --- |
| 1 | `setup` | yes | <0.1 | `knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6.5` |
| 2 | `libs` | yes | 0.1 | `library(bartisan)` |
| 3 | `simdata` | yes | <0.1 | `n <- 300` |
| 4 | `loo` | yes | 1.2 | `loo::loo_compare(` |
| 5 | `default` | yes | 0.1 | `d$binary <- rbinom(n, 1, 0.4)` |
| 6 | `perforest` | yes | 0.3 | `d$het <- rnorm(n, 2 * d$x1, exp(-1 + 1.5 * d$x2))` |
| 7 | `sharesparsity` | yes | 1.0 | `set.seed(9)` |
| 8 | `emptyforest` | yes | 0.3 | `set.seed(8)` |
| 9 | `dpm` | yes | 0.1 | `d$heavy <- 2 * sin(pi * d$x1) + d$x2 + rt(n, df = 3)` |
| 10 | `dpmplot` | yes | 0.4 | `plot(error_density(fit_dpm))` |
| 11 | `zerosim` | yes | <0.1 | `set.seed(11)` |
| 12 | `zerohist` | yes | 0.2 | `library(ggplot2)` |
| 13 | `scalegap` | no | <0.1 | `bartisan(score ~ ., data = d, family = ordinal())` |
| 14 | `ordcont` | yes | 0.4 | `nbins <- 25` |
| 15 | `beta` | yes | 1.0 | `d$rate <- rbeta(n, plogis(1.5 * sin(pi * d$x1)) * 12,` |
| 16 | `tweedie` | yes | 0.3 | `mu <- exp(1.5 + sin(pi * d$x1))` |
| 17 | `aft` | yes | 0.3 | `d$time <- rexp(n, exp(-(1 + d$x1)))` |
| 18 | `phsurv` | yes | 0.3 | `fit_ph <- bartisan(survival::Surv(time, event) ~ x1 + x2,` |
| 19 | `custom` | yes | 2.8 | `pois_by_hand <- custom_family(` |
| 20 | `customaux` | yes | 1.0 | `gauss_by_hand <- custom_family(` |

### survival (5.8 s)

| # | Chunk | Evaluated | Seconds | First line |
| ---: | --- | --- | ---: | --- |
| 1 | `setup` | yes | 0.8 | `knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,` |
| 2 | `libs` | yes | 0.1 | `library(bartisan)` |
| 3 | `data` | yes | <0.1 | `n <- 400` |
| 4 | `aftfits` | yes | 0.6 | `fit_w <- bartisan(Surv(time, status) ~ x1 + x2 + x3 + trt, data = d,` |
| 5 | `dpmfit` | yes | 0.7 | `fit_dpm <- bartisan(Surv(time, status) ~ x1 + x2 + x3 + trt, data = d,` |
| 6 | `errdens` | yes | 0.3 | `ed <- error_density(fit_dpm)` |
| 7 | `phfit` | yes | 0.4 | `fit_ph <- bartisan(Surv(time, status) ~ x1 + x2 + x3 + trt, data = d,` |
| 8 | `survtype` | yes | <0.1 | `head(predict(fit_ph, type = "survival", times = c(1, 2, 5)), 3)` |
| 9 | `survest` | yes | 0.1 | `library(marginaleffects)` |
| 10 | `hazshapes` | yes | 0.2 | `tt <- seq(0.05, 3, length.out = 200)` |
| 11 | `rmsefig` | yes | 0.2 | `ggplot(agg, aes(s_rmse, family, color = family)) +` |
| 12 | `rmsetable` | yes | <0.1 | `knitr::kable(pivot("s_rmse"), digits = 3, row.names = FALSE,` |
| 13 | `curvefig` | yes | 0.3 | `keep <- c("truth", "weibull_aft()", "dpm_aft()", "ph()",` |
| 14 | `rankfig` | yes | 0.2 | `ggplot(agg, aes(rank, family, color = family)) +` |
| 15 | `densfig` | yes | 0.2 | `dn <- res$dens[res$dens$truth %in% c("bimodal errors", "log-normal err` |
| 16 | `lstable` | yes | <0.1 | `knitr::kable(pivot("logscore", setdiff(fam_levels, "discrete-time prob` |
| 17 | `censfig` | yes | 0.2 | `sw <- res$sweep_agg` |
| 18 | `dtcode` | yes | 0.4 | `edges <- quantile(d$time[d$status == 1], seq(0.05, 0.95, length.out = ` |
| 19 | `dtsurv` | yes | <0.1 | `grid_dat <- d[rep(1:3, each = length(edges)), c("x1", "x2", "x3", "trt` |
| 20 | `dtcompare` | yes | <0.1 | `sel <- agg[agg$truth %in% c("hazard turns over", "crossing hazards") &` |
| 21 | `jacobian` | yes | 0.1 | `log_score_T <- function(fit, newdata) {` |
| 22 | `binsweep` | yes | <0.1 | `knitr::kable(res$bins, digits = 3, row.names = FALSE,` |
| 23 | `timefig` | yes | 0.1 | `tm <- res$timing` |

### implementation (3.9 s)

| # | Chunk | Evaluated | Seconds | First line |
| ---: | --- | --- | ---: | --- |
| 1 | `setup` | yes | <0.1 | `knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6.5` |
| 2 | `libs` | yes | 0.1 | `library(bartisan)` |
| 3 | `control` | yes | <0.1 | `ctrl <- bartisan_control(num_trees = 20, num_burn = 300, num_draws = 3` |
| 4 | `softhard` | yes | 1.1 | `friedman <- function(n) {` |
| 5 | `sparsitytable` | yes | <0.1 | `knitr::kable(data.frame(` |
| 6 | `missing` | yes | 0.4 | `d_miss <- train` |
| 7 | `chains` | yes | 1.5 | `fit_chains <- bartisan(y ~ . - eta, data = train, family = gaussian(),` |

