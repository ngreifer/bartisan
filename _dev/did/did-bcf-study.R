# DiD-BCF (Souto & Louzada Neto 2025, arXiv:2505.09706) fitted with bartisan.
#
# What is tested: that bartisan() + vc() fits their Eq. 3,
#   Y_it = mu(D_i, t, X_it) + tau(X_it, k_it) * D_it,
# and recovers ATT and CATT. Two ways it could fail, both from the
# specification rather than from BART: (i) mu sees the ever-treated indicator
# and calendar time, whose interaction IS D_it, so the model is not
# likelihood-identified and the mu/tau split rests on the priors; (ii) tau may
# vary with event time, so it could invent dynamics where the truth is flat.
#
# Outcome measured: bias, RMSE and 95% interval coverage of the ATT over 20
# replicates per cell; RMSE of the per-observation CATT against the known truth
# in cell B; the event-study profile over treated rows, which should be flat.
# Read against static TWFE with unit and period fixed effects.
#
# Cell A = their DGP 2 / Setting 1 (staggered, homogeneous tau = 3, linear;
# TWFE correctly specified). Cell B = their DGP 5 / Setting 3 (staggered +
# selection on statics, heterogeneous tau, all covariates non-linear,
# quadratic trend; TWFE badly misspecified).
#
# Reading: small ATT bias in both cells and a CATT RMSE well under the spread
# of the true CATT means the mapping works and the vignette is worth writing.
# Cell A clean but cell B biased toward zero means mu is absorbing the effect
# and the mapping needs the identification repaired first.

suppressPackageStartupMessages(library(bartisan))
SP <- "_dev"
source(file.path(SP, "did-bcf-dgp.R"))
source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")

REPS <- 20L
OUT <- file.path(SP, "did-study.rds")

fit_did <- function(d, covs) {
  rhs <- paste(c("cohort", "t", covs), collapse = " + ")
  mods <- paste(c(covs, "t", "cohort"), collapse = " + ")
  f <- stats::as.formula(sprintf("y ~ %s + vc(Dit, ~ %s)", rhs, mods))
  suppressMessages(bartisan(f, data = d, family = gaussian(),
    control = bartisan_control(num_trees = c(50L, 25L), num_burn = 400L,
                               num_draws = 600L, verbose = FALSE)))
}

twfe <- function(d, covs) {
  f <- stats::as.formula(sprintf("y ~ Dit + %s | id + t", paste(covs, collapse = " + ")))
  m <- fixest::feols(f, data = d)
  ci <- stats::confint(m)
  c(est = unname(stats::coef(m)["Dit"]),
    lo = unname(ci["Dit", 1]), hi = unname(ci["Dit", 2]))
}

one_rep <- function(cell, r) {
  gen <- if (cell == "A") dgp_A else dgp_B
  covs <- if (cell == "A") covs_A else covs_B
  d <- gen(seed = 1000L + r)
  tr <- d$Dit == 1

  fit <- fit_did(d, covs)
  cf <- coef(fit, draws = TRUE)[[1L]]              # draws x n, tau(X_it, k_it)
  att_draws <- rowMeans(cf[, tr, drop = FALSE])
  att_true <- mean(d$tau_true[tr])
  tau_hat <- colMeans(cf)

  tw <- twfe(d, covs)

  es <- tapply(tau_hat[tr], d$k[tr], mean)
  gatt <- tapply(tau_hat[tr], droplevels(d$cohort[tr]), mean)

  list(
    summary = data.frame(
      cell = cell, rep = r, att_true = att_true,
      bart_att = mean(att_draws),
      bart_lo = unname(quantile(att_draws, .025)),
      bart_hi = unname(quantile(att_draws, .975)),
      catt_rmse = sqrt(mean((tau_hat[tr] - d$tau_true[tr])^2)),
      catt_cor = if (stats::sd(d$tau_true[tr]) > 0)
        stats::cor(tau_hat[tr], d$tau_true[tr]) else NA_real_,
      catt_sd_true = stats::sd(d$tau_true[tr]),
      twfe_att = tw[["est"]], twfe_lo = tw[["lo"]], twfe_hi = tw[["hi"]],
      stringsAsFactors = FALSE),
    es = data.frame(cell = cell, rep = r, k = as.numeric(names(es)),
                    tau = as.vector(es), stringsAsFactors = FALSE),
    gatt = data.frame(cell = cell, rep = r, cohort = names(gatt),
                      tau = as.vector(gatt), stringsAsFactors = FALSE))
}

cells <- c("A", "B")
pr <- prog_init(total = length(cells) * REPS, title = "DiD-BCF via bartisan",
                unit = "replicate", kind = "simulation",
                root = "/Users/NoahGreifer/Dropbox/Research/R/bartisan/.progress-runs")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list(); ess <- list(); gatts <- list()
for (cell in cells) {
  for (r in seq_len(REPS)) {
    t0 <- proc.time()[[3]]
    out <- one_rep(cell, r)
    rows[[length(rows) + 1L]] <- out$summary
    ess[[length(ess) + 1L]] <- out$es
    gatts[[length(gatts) + 1L]] <- out$gatt
    prog_tick(pr, label = sprintf("cell %s rep %d (%.0fs)", cell, r,
                                  proc.time()[[3]] - t0))
    # Checkpoint every replicate: a killed run keeps everything finished.
    saveRDS(list(summary = do.call(rbind, rows), es = do.call(rbind, ess),
                 gatt = do.call(rbind, gatts), complete = FALSE,
                 done = length(rows), total = length(cells) * REPS), OUT)
  }
}
on.exit(); prog_end(pr, "done")
saveRDS(list(summary = do.call(rbind, rows), es = do.call(rbind, ess),
             gatt = do.call(rbind, gatts), complete = TRUE,
             done = length(rows), total = length(rows)), OUT)
cat("DONE\n")
