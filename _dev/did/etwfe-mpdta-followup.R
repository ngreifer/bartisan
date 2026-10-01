# Two follow-ups to `_dev/did/etwfe-mpdta.R`: whether the per-period time
# structure is converged, and a pre-trend check from the ETWFE imputation.
#
# What is being tested (1). In `etwfe-mpdta.R` the ATT's bulk ESS for the
# additive specification, `vc(one, ~ period + lpop)`, fell from 137 to 114 as
# the chains went from 1000 to 2500 kept draws, which DID.md phase 13c reads as
# a chain still drifting, while the per-period specification, one coefficient
# forest per year over lpop (`vc(periodf, ~ lpop)`, the direct analogue of
# Wooldridge's fs_t and fs_t * x terms), reached 870 at 1000. The claim is that
# the per-period form is converged: its ESS rises at 2500. Read: the ATT's bulk
# ESS and R-hat at 2500 against 870 and 1.00 at 1000. If it rises, the
# per-period form is the time structure to recommend; if it falls too, neither
# form's ATT ESS can be read at one length and the notes say to run two.
#
# What is being tested (2). Wooldridge (2025, sec. 6 and 6.1) gets the
# event-study ("leads and lags") estimator from the same imputation by giving
# each cohort's periods s <= g - 2 indicators of their own in the first stage,
# which makes g - 1 the reference. In bartisan free indicators for those rows
# are the same as leaving them out of the fit and imputing them, so their gaps
# Y - Y0 are the placebo estimates. The claim is that these cover zero on mpdta
# and sit near the linear leads estimates. Unlike the Eq. 2 placebo of DID.md
# phase 8, which returned 0.67-0.99 where the truth was zero and so could not
# fail, this one can: a cohort whose trend departs from the controls' before
# treatment leaves gaps its level alone cannot absorb. Read: placebo ATTs at
# event times -4, -3 and -2 (cohort 2007 has leads at 2003-2005, cohort 2006 at
# 2003-2004, cohort 2004 none), posterior mean and 95% interval, against the
# linear imputation fitted to the same rows and did's universal-base placebos
# recorded in DID.md (+0.0069, +0.0276, +0.0235), under both time structures.
# If they cover zero and sit near the linear ones, the pre-trend check is
# available from the imputation, filling the gap DID.md records for the
# level-form models; if they are far from the linear ones or exclude zero where
# those do not, leaning on one reference period makes the forests' cohort level
# unstable and the check should use the linear leads regression instead.
#
# Both parts: 4 chains, hard gates, 500 warmup, the specifications of
# `etwfe-mpdta.R` with a county random intercept.

source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")
suppressMessages({
  library(bartisan)
  library(did)
})

out_file <- Sys.getenv("OUT", "_dev/did/etwfe-mpdta-followup.rds")
future::plan(future::multisession, workers = 4L)

data(mpdta, package = "did")
d <- mpdta
d$cohort <- factor(ifelse(d$first.treat == 0, "never", d$first.treat),
                   levels = c("never", "2004", "2006", "2007"))
d$w <- as.integer(d$first.treat > 0 & d$year >= d$first.treat)
d$lead <- as.integer(d$first.treat > 0 & d$year <= d$first.treat - 2)
d$period <- d$year
d$periodf <- factor(d$year)
d$one <- 1
d$e <- ifelse(d$first.treat > 0, d$year - d$first.treat, NA)

specs <- list(
  additive_re = list(f = lemp ~ cohort + lpop + vc(one, ~ period + lpop) +
                       (1 | countyreal), trees = c(100L, 50L)),
  per_period_re = list(f = lemp ~ cohort + lpop + vc(periodf, ~ lpop) +
                         (1 | countyreal), trees = NULL))

fit_one <- function(spec, data, draws) {
  args <- list(formula = spec$f, data = data, family = gaussian(), chains = 4L,
               gate = "hard", num_burn = 500L, num_draws = draws,
               verbose = FALSE)
  if (!is.null(spec$trees)) {
    args$num_trees <- spec$trees
  }
  set.seed(2026)
  do.call(bartisan, args)
}

pr <- prog_init(total = 3L, title = "ETWFE follow-ups on mpdta", unit = "fit",
                kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

# (1) The per-period form at a second chain length.
ctl <- d[d$w == 0, ]
tr <- d[d$w == 1, ]
fit <- fit_one(specs$per_period_re, ctl, 2500L)
y0 <- predict(fit, newdata = tr, draws = TRUE)
att <- rowMeans(sweep(-y0, 2L, tr$lemp, "+"))
m <- matrix(att, length(att) / 4L, 4L)
length_check <- data.frame(
  spec = "per_period_re", draws = 2500L, estimate = mean(att),
  lower = unname(stats::quantile(att, .025)),
  upper = unname(stats::quantile(att, .975)),
  rhat = posterior::rhat(m), ess_bulk = posterior::ess_bulk(m),
  ess_tail = posterior::ess_tail(m))
saveRDS(list(length_check = length_check, complete = FALSE), out_file)
prog_tick(pr, label = "per_period_re, 2500 draws")

# (2) Leads, under both time structures.
fit_rows <- d[d$w == 0 & d$lead == 0, ]
target <- d[d$w == 1 | d$lead == 1, ]
lin <- lm(lemp ~ cohort * lpop + periodf * lpop, data = fit_rows)
target$gap_lin <- target$lemp - predict(lin, newdata = target)

leads <- list()
for (nm in names(specs)) {
  fit <- fit_one(specs[[nm]], fit_rows, 1000L)
  y0 <- predict(fit, newdata = target, draws = TRUE)
  gap <- sweep(-y0, 2L, target$lemp, "+")
  for (e in sort(unique(target$e))) {
    on <- target$e == e
    v <- rowMeans(gap[, on, drop = FALSE])
    leads[[length(leads) + 1L]] <- data.frame(
      spec = nm, e = e, n = sum(on), bart = mean(v),
      lower = unname(stats::quantile(v, .025)),
      upper = unname(stats::quantile(v, .975)),
      linear = mean(target$gap_lin[on]))
  }
  saveRDS(list(length_check = length_check, leads = do.call(rbind, leads),
               complete = FALSE), out_file)
  prog_tick(pr, label = paste(nm, "leads"))
}

res <- list(length_check = length_check, leads = do.call(rbind, leads),
            complete = TRUE)
saveRDS(res, out_file)
on.exit()
prog_end(pr, "done")

print(res$length_check, row.names = FALSE, digits = 3)
cat("\n")
print(res$leads, row.names = FALSE, digits = 3)
