# Turns the raw simulation output into the compact object the survival vignette
# reads, so that the vignette builds without refitting anything.
#
#   Rscript _dev/survival-sim.R      # ~40 min
#   Rscript _dev/survival-bins.R     # ~5 min
#   Rscript _dev/survival-timing.R   # ~2 min, wants a quiet machine
#   Rscript _dev/survival-results.R  # writes vignettes/survival-results.rds

sim <- readRDS("_dev/survival-sim.rds")
bins <- readRDS("_dev/survival-bins.rds")
timings <- readRDS("_dev/survival-timing.rds")

fam_levels <- c("weibull_aft()", "loglogistic_aft()", "lognormal_aft()",
                "dpm_aft()", "ph()", "discrete-time probit")
truth_order <- c("Weibull PH", "hazard turns over", "log-normal errors",
                 "bimodal errors", "heavy-tailed errors", "crossing hazards")

se <- function(z) if (length(z) < 2) 0 else sd(z) / sqrt(length(z))

# Mean over replicates, with a standard error on the primary metric.
agg <- do.call(rbind, lapply(
  split(sim$main, list(sim$main$truth, sim$main$family), drop = TRUE),
  function(z) data.frame(
    truth = z$truth[1], family = z$family[1],
    s_rmse = mean(z$s_rmse), s_rmse_se = se(z$s_rmse),
    s_worst = mean(z$s_worst), rank = mean(z$rank),
    logscore = mean(z$logscore), secs = mean(z$secs))))
agg <- agg[order(match(agg$truth, truth_order),
                 match(agg$family, fam_levels)), ]
rownames(agg) <- NULL

sweep_agg <- do.call(rbind, lapply(
  split(sim$sweep, list(sim$sweep$censor, sim$sweep$family), drop = TRUE),
  function(z) data.frame(censor = z$censor[1], family = z$family[1],
                         s_rmse = mean(z$s_rmse), rank = mean(z$rank))))
rownames(sweep_agg) <- NULL

# From the dedicated timing run, not from sim$main, whose `secs` column is
# measured under whatever else the machine was doing.
timing <- do.call(rbind, lapply(split(timings, timings$family), function(z)
  data.frame(family = z$family[1], secs = mean(z$secs))))
timing <- timing[order(match(timing$family, fam_levels)), ]
rownames(timing) <- NULL

# The bin sweep, one column per truth, formatted as text so the two truths sit
# side by side in one narrow table.
bagg <- do.call(rbind, lapply(
  split(bins, list(bins$truth, bins$bins), drop = TRUE),
  function(z) data.frame(truth = z$truth[1], bins = z$bins[1],
                         s_rmse = mean(z$s_rmse), r_rmse = mean(z$r_rmse),
                         p_loo = mean(z$p_loo), bad_k = mean(z$bad_k))))

pair <- function(col, fmt) {
  vapply(sort(unique(bagg$bins)), function(b) {
    a <- bagg[bagg$bins == b & bagg$truth == "hazard turns over", col]
    w <- bagg[bagg$bins == b & bagg$truth == "Weibull PH", col]
    sprintf(paste0(fmt, " / ", fmt), a, w)
  }, character(1))
}

bagg$bad_pct <- 100 * bagg$bad_k

bins_tab <- data.frame(
  Bins = sort(unique(bagg$bins)),
  `S(t | x) RMSE` = pair("s_rmse", "%.3f"),
  `r(x) RMSE` = pair("r_rmse", "%.3f"),
  `effective parameters` = pair("p_loo", "%.0f"),
  `Pareto k above 0.7` = pair("bad_pct", "%.1f%%"),
  check.names = FALSE)

meta <- list(n_train = 700L, n_test = 700L, num_trees = 50L, num_draws = 500L,
             num_burn = 500L, n_rep = 5L, truth_order = truth_order,
             fam_levels = fam_levels, when = sim$when)

saveRDS(list(meta = meta, agg = agg, curves = sim$curves, dens = sim$dens,
             sweep_agg = sweep_agg, timing = timing, bins = bins_tab),
        "vignettes/survival-results.rds")

cat("wrote vignettes/survival-results.rds\n")
print(agg[, c("truth", "family", "s_rmse", "rank", "logscore")], digits = 3)
print(bins_tab)
