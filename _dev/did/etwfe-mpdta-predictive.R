# The mpdta ETWFE intervals with the treated rows' own noise put back in.
#
# What is being tested. `_dev/did/etwfe-sim.R` found that imputing the
# conditional mean of Y(inf) undercovers where the model is right (the linear
# cell: 90% for the overall ATT, 80% at e = 1) and that imputing posterior
# predictive draws of it restores 95%. The conditional mean leaves out the
# treated rows' own noise, whose average over n rows has sd sigma / sqrt(n). On
# mpdta the conditional-mean imputation gave an overall interval 0.9 times the
# width of the linear ETWFE's, and a placebo interval at e = -3 that excluded
# zero by 0.001 where did's covered it (DID.md, phase 15). The claim is that
# both come from the omitted noise: with predictive draws the overall interval
# widens to near the linear width and the e = -3 placebo interval covers zero.
#
# The measurement. The four fits of `etwfe-mpdta.R` and
# `etwfe-mpdta-followup.R` (additive and per-period time structures, each fit
# to all untreated rows and to the untreated rows less the leads), at the same
# seed and settings, so the conditional-mean numbers reproduce the recorded
# ones and check that the fits are the same. From each fit, both imputations:
# `predict(draws = TRUE)` and `rstantools::posterior_predict()`. Read the
# overall ATT's interval width against 0.049 (linear ETWFE, county-clustered)
# and 0.045 (the conditional mean), then the lower bound of the e = -3 placebo
# interval against zero, with did's placebo widths 0.093, 0.069 and 0.058 at
# e = -4, -3 and -2 for context.
#
# What each outcome means. Widths near the linear ones and e = -3 covering zero:
# the predictive imputation is the one to recommend, the notes withdraw "0.9
# times the linear width" for the width the predictive draws give, and the
# e = -3 exclusion was the omitted noise. Widths barely changed: the omitted
# noise is small on mpdta, the narrower intervals are the model's own, and the
# e = -3 exclusion stands as a borderline pre-trend signal that did's wider
# interval does not flag.

source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")
suppressMessages({
  library(bartisan)
  library(did)
})

out_file <- Sys.getenv("OUT", "_dev/did/etwfe-mpdta-predictive.rds")
# SMOKE=1 shrinks every chain to check the script end to end in seconds.
smoke <- nzchar(Sys.getenv("SMOKE"))
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
designs <- list(
  lags_only = list(fit = d$w == 0, target = d$w == 1),
  leads = list(fit = d$w == 0 & d$lead == 0, target = d$w == 1 | d$lead == 1))

fit_one <- function(spec, data) {
  args <- list(formula = spec$f, data = data, family = gaussian(), chains = 4L,
               gate = "hard", num_burn = if (smoke) 20L else 500L,
               num_draws = if (smoke) 20L else 1000L,
               verbose = FALSE)
  if (!is.null(spec$trees)) {
    args$num_trees <- spec$trees
  }
  set.seed(2026)
  do.call(bartisan, args)
}

summ <- function(v) {
  c(estimate = mean(v), lower = unname(stats::quantile(v, .025)),
    upper = unname(stats::quantile(v, .975)))
}

pr <- prog_init(total = 4L, title = "ETWFE on mpdta, predictive imputation",
                unit = "fit", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)
rows <- list()

for (dn in names(designs)) {
  fit_rows <- d[designs[[dn]]$fit, ]
  target <- d[designs[[dn]]$target, ]
  es <- sort(unique(target$e))
  groups <- c(list(overall = target$w == 1),
              setNames(lapply(es, function(e) target$e == e), paste0("e", es)))
  for (sn in names(specs)) {
    fit <- fit_one(specs[[sn]], fit_rows)
    y0 <- list(mean = predict(fit, newdata = target, draws = TRUE),
               predictive = rstantools::posterior_predict(fit, newdata = target))
    for (im in names(y0)) {
      gap <- sweep(-y0[[im]], 2L, target$lemp, "+")
      for (q in names(groups)) {
        rows[[length(rows) + 1L]] <- data.frame(
          design = dn, spec = sn, imputation = im, quantity = q,
          n = sum(groups[[q]]),
          t(summ(rowMeans(gap[, groups[[q]], drop = FALSE]))))
      }
    }
    saveRDS(list(res = do.call(rbind, rows), complete = FALSE), out_file)
    prog_tick(pr, label = paste(dn, sn))
  }
}

res <- do.call(rbind, rows)
res$width <- res$upper - res$lower
saveRDS(list(res = res, complete = TRUE), out_file)
on.exit()
prog_end(pr, "done")

print(res, row.names = FALSE, digits = 3)
