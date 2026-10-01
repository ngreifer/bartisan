# The lalonde bcf(family = dpm()) levels, before and after the wrapper shift.
#
# What is tested. `_dev/TASKS.md` records the vignette("causal") lalonde fit,
# `bcf(re78 ~ ..., treat = ~ treat, family = dpm(), chains = 4)`, putting Y[0]
# among the treated at 5420 and 5650 in two four-chain builds, 4780 with two
# chains and 3780 with one, against observed treated earnings of 6349, with the
# ATT intervals all overlapping. `_dev/dpm-level-shift.R` found the cause:
# `VaryingCoefficientFamily` did not forward `report_shift()`, so a `bcf()` fit
# recorded the raw forest and its level was the coordinate the likelihood does
# not identify. The claim here: with the shift forwarded, the treated
# potential-outcome levels agree across seeds and chain counts to within their
# Monte Carlo error, and the ATT does not move.
#
# What is measured. The vignette chunk as it was before db7125b, at one and four
# chains, two seeds each. Per fit, from `estimate_effect(estimand = "ATT")` and
# `diagnose()` on its result:
#
#   att, att_lo, att_hi   the ATT and its 95% interval
#   y0, y1                the average potential outcomes among the treated
#   observed              the observed treated mean, 6349, for reference
#   cor_y1_center         correlation across draws of the Y[1] level with
#                         aux.center, the raw mixture mean
#   rhat_*, ess_*         split R-hat and bulk ESS for Y[1], Y[0] and the ATT
#   mcse_y1               sd(Y[1] draws) / sqrt(ess), the Monte Carlo error
#   *_raw                 the same for the Y[1] level in the raw chart, which
#                         is the recorded level minus aux.center per draw: the
#                         chart change does not touch the sampler, so this is
#                         exactly what the fit reported before the fix
#
# Read the spread of y1 across the four fits against mcse_y1, and the cor column.
#
# What each outcome means. Before the fix: y1 spreads over hundreds to thousands
# of dollars, cor near -1, ESS in single digits, which is the record in
# TASKS.md restated with seeds. After it: cor near 0, ESS in the hundreds, the
# spread of y1 within about two Monte Carlo errors, and the ATT column
# unchanged from before. If the spread survives the fix, something other than
# the chart moves the level. A y1 that is stable but sits away from 6349 is the
# model rather than the sampler: one error shape for every x, where the share
# of zero earnings depends on x.
#
# Run with: DPM_LEVEL_TAG=before Rscript _dev/dpm-lalonde-level.R
# Writes:   _dev/dpm-lalonde-level-<tag>.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

TAG <- Sys.getenv("DPM_LEVEL_TAG", "run")
OUT <- sprintf("_dev/dpm-lalonde-level-%s.rds", TAG)

data("lalonde", package = "cobalt")
observed <- mean(lalonde$re78[lalonde$treat == 1])

grid <- expand.grid(chains = c(1L, 4L), seed = c(11L, 12L))

# Exact names: the contrast row is "Y[1] - Y[0]", which a substring match on
# "Y[1]" would pick up first.
pick <- function(tab, name, col) {
  at <- which(tab[["quantity"]] == name)
  if (!length(at)) return(NA_real_)
  tab[[col]][at[1L]]
}

pr <- prog_init(total = nrow(grid), title = sprintf("lalonde dpm levels (%s)", TAG),
                unit = "fit", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
draws <- list()

for (k in seq_len(nrow(grid))) {
  chains <- grid$chains[k]
  seed <- grid$seed[k]

  set.seed(seed)
  t0 <- Sys.time()
  fit <- bcf(re78 ~ age + educ + race + married + nodegree + re74 + re75,
             treat = ~ treat, data = lalonde, family = dpm(), chains = chains)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  eff <- estimate_effect(fit, estimand = "ATT")
  po <- attr(eff, "po_draws")
  att <- attr(eff, "draws")[[1L]]
  center <- fit[["aux"]][, "center"]
  tab <- as.data.frame(diagnose(eff)[["table"]])
  contrast <- tab[grepl(" - ", tab[["quantity"]], fixed = TRUE), , drop = FALSE]

  ess_y1 <- pick(tab, "Y[1]", "ess_bulk")

  by_chain <- function(x) matrix(x, ncol = chains)
  raw_y1 <- po[["Y[1]"]] - center
  ess_y1_raw <- posterior::ess_bulk(by_chain(raw_y1))
  draws[[k]] <- list(po = po, att = att, center = center)

  rows[[k]] <- data.frame(
    tag = TAG, chains = chains, seed = seed, secs = secs,
    att = mean(att),
    att_lo = unname(quantile(att, 0.025)),
    att_hi = unname(quantile(att, 0.975)),
    y0 = mean(po[["Y[0]"]]),
    y1 = mean(po[["Y[1]"]]),
    observed = observed,
    sd_y1 = sd(po[["Y[1]"]]),
    cor_y1_center = cor(po[["Y[1]"]], center),
    rhat_y1 = pick(tab, "Y[1]", "rhat"),
    ess_y1 = ess_y1,
    mcse_y1 = sd(po[["Y[1]"]]) / sqrt(ess_y1),
    rhat_y0 = pick(tab, "Y[0]", "rhat"),
    ess_y0 = pick(tab, "Y[0]", "ess_bulk"),
    rhat_att = if (nrow(contrast)) contrast[["rhat"]][1L] else NA_real_,
    ess_att = if (nrow(contrast)) contrast[["ess_bulk"]][1L] else NA_real_,
    sd_y1_raw = sd(raw_y1),
    rhat_y1_raw = posterior::rhat(by_chain(raw_y1)),
    ess_y1_raw = ess_y1_raw,
    mcse_y1_raw = sd(raw_y1) / sqrt(ess_y1_raw))

  prog_tick(pr, label = sprintf("chains %d seed %d", chains, seed))
  saveRDS(list(rows = do.call(rbind, rows), draws = draws,
               complete = k == nrow(grid), done = k, total = nrow(grid)), OUT)
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
rownames(res) <- NULL
print(res[, c("chains", "seed", "secs", "att", "att_lo", "att_hi", "y0", "y1",
              "observed")], digits = 4)
print(res[, c("chains", "seed", "sd_y1", "cor_y1_center", "rhat_y1", "ess_y1",
              "mcse_y1", "rhat_y0", "ess_y0", "rhat_att", "ess_att")], digits = 3)
print(res[, c("chains", "seed", "sd_y1_raw", "rhat_y1_raw", "ess_y1_raw",
              "mcse_y1_raw")], digits = 3)
