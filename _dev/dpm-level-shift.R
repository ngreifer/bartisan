# Is the level of a bcf(family = dpm()) fit the sampler's raw forest level?
#
# What is tested. `_dev/TASKS.md` records that the lalonde `bcf(family = dpm())`
# fit gave treated-group levels of 3780 to 5650 across builds while the ATT
# barely moved. The response-scale mean of a `dpm()` fit is the forest plus the
# mixture's mean, and only the sum is identified: the forest's level and the
# mixture's center trade off along a direction the likelihood does not
# constrain, with a correlation of -0.987 across draws measured in TASKS.md.
# `DPMFamily::report_shift()` in `src/family.cpp` exists to report in the chart
# where the mixture has mean zero, so the recorded predictor is the sum.
# `VaryingCoefficientFamily`, which `bcf()` fits through, does not override
# `report_shift()`, so the base class's zero shift applies and the recorded
# predictor is the raw forest. The claim: a `bcf()` fit reports the
# unidentified coordinate and a plain `bartisan()` fit does not.
#
# What is measured. One simulated data set, skewed errors with mean zero so the
# raw mixture center is well away from zero, two chains of 200 + 400 sweeps.
# Two fits of the same mean structure: `bartisan()` with the treatment as a
# covariate, and `bcf()` with `propensity = FALSE`. Per draw, `level` is the
# mean over observations of the reported predictor and `center` is the raw
# mixture mean from `fit$aux`. The columns to read, per fit:
#
#   cor         correlation of level with center across draws
#   sd_level    sd of level across draws
#   sd_sum      sd of level + center, the sum the likelihood identifies
#   chain_gap   difference between the two chains' mean levels
#   dens_gap    largest gap between the per-draw sum of
#               predict(type = "density", log = TRUE) and fit$loglik
#   rhat_level, ess_level, rhat_center, ess_center  from diagnose()
#
# What each outcome means. If the bcf fit has cor near -1, sd_level many times
# sd_sum, a chain_gap of the size of sd_level and a dens_gap far above rounding
# while the plain fit has cor near 0, sd_level near sd_sum and dens_gap at
# rounding, then the lalonde level problem is the shift missing from the
# wrapper, and forwarding `report_shift()` through it is the fix. If the two
# fits look alike, the wrapper is not the cause and the cause is in the ridge
# dynamics or in the error shape varying with x on lalonde.
#
# Run with: DPM_LEVEL_TAG=before Rscript _dev/dpm-level-shift.R
# Writes:   _dev/dpm-level-shift-<tag>.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

TAG <- Sys.getenv("DPM_LEVEL_TAG", "run")
OUT <- sprintf("_dev/dpm-level-shift-%s.rds", TAG)

CHAINS <- 2L
BURN <- 200L
DRAWS <- 400L

set.seed(2026)
n <- 1000L
d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n),
                z = rbinom(n, 1L, 0.5))
truth <- 2 + 3 * sin(pi * d$x1) + 2 * d$x2 + 1 * d$z
# Skewed, mean zero, sd about 2.45.
d$y <- truth + 3 * (rgamma(n, 1.5, 1.5) - 1)

diagnosis_table <- function(fit) {
  dg <- diagnose(fit)
  tab <- if (is.data.frame(dg)) dg else dg[["table"]]
  as.data.frame(tab)
}

row_of <- function(tab, pattern) {
  at <- grep(pattern, tab[["quantity"]])
  if (!length(at)) return(c(rhat = NA_real_, ess = NA_real_))
  c(rhat = tab[["rhat"]][at[1L]], ess = tab[["ess_bulk"]][at[1L]])
}

measure <- function(fit, label) {
  link <- predict(fit, type = "link", draws = TRUE)
  level <- rowMeans(link)
  center <- fit[["aux"]][, "center"]
  chain <- rep(seq_len(CHAINS), each = DRAWS)
  chain_means <- tapply(level, chain, mean)

  dens <- predict(fit, type = "density", draws = TRUE, log = TRUE)
  dens_gap <- max(abs(rowSums(dens) - fit[["loglik"]]))

  tab <- diagnosis_table(fit)
  lv <- row_of(tab, "^eta.*average|average.*eta|^eta \\(average|eta.*mean")
  ce <- row_of(tab, "aux\\.center")

  data.frame(fit = label,
             cor = cor(level, center),
             sd_level = sd(level),
             sd_sum = sd(level + center),
             mean_level = mean(level),
             mean_y = mean(d$y),
             chain_gap = unname(diff(chain_means)),
             dens_gap = dens_gap,
             rhat_level = lv[["rhat"]], ess_level = lv[["ess"]],
             rhat_center = ce[["rhat"]], ess_center = ce[["ess"]],
             tag = TAG)
}

pr <- prog_init(total = 2L, title = sprintf("dpm level: plain vs bcf (%s)", TAG),
                unit = "fit", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
traces <- list()

plain <- bartisan(y ~ x1 + x2 + x3 + z, data = d, family = dpm(),
                  chains = CHAINS, num_burn = BURN, num_draws = DRAWS)
rows[["plain"]] <- measure(plain, "plain")
traces[["plain"]] <- cbind(level = rowMeans(predict(plain, type = "link", draws = TRUE)),
                           center = plain[["aux"]][, "center"])
prog_tick(pr, label = "plain")
saveRDS(list(rows = do.call(rbind, rows), traces = traces, complete = FALSE), OUT)

vc <- bcf(y ~ x1 + x2 + x3, treat = ~ z, data = d, family = dpm(),
          propensity = FALSE, chains = CHAINS, num_burn = BURN, num_draws = DRAWS)
rows[["bcf"]] <- measure(vc, "bcf")
traces[["bcf"]] <- cbind(level = rowMeans(predict(vc, type = "link", draws = TRUE)),
                         center = vc[["aux"]][, "center"])
prog_tick(pr, label = "bcf")

res <- do.call(rbind, rows)
rownames(res) <- NULL
saveRDS(list(rows = res, traces = traces, complete = TRUE), OUT)

on.exit()
prog_end(pr, "done")

print(res, digits = 4)
