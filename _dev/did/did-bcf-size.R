# Size and power of the DiD-BCF model (Souto & Louzada Neto 2025,
# arXiv:2505.09706) fitted with bartisan, against the standard DiD estimators.
#
# What is tested. Whether the false-alarm rate of a bartisan DiD fit is at its
# nominal level. The concern is specific and came from a measurement: fitting
# their Eq. 2 form -- tau multiplied by the ever-treated indicator, so that tau
# is free before treatment -- returned tau-hat of 0.67 to 0.99 at every
# pre-treatment event time where the truth is exactly zero. If the forest
# manufactures an effect there, it may manufacture one under a true null too,
# in which case the Eq. 3 form's ATT interval would exclude zero more often
# than 5% of the time and the method would not be usable for testing.
#
# What is measured. Four cells, 100 Monte Carlo replicates each, N = 200 units
# by T = 8 periods, staggered adoption. A1/B1 carry a real effect and A0/B0 set
# it to exactly zero; A is their DGP 2 / Setting 1 (linear, homogeneous) and B
# their DGP 5 / Setting 3 (non-linear, selection into cohort, heterogeneous).
# The decisive endpoint on the null cells is the rejection rate of the ATT at
# two thresholds: the paper's own p_Bayes = min(P(tau>0), P(tau<0)) < 0.05,
# which corresponds to a 90% interval excluding zero, and the 95% credible
# interval excluding zero. Alongside, the share of individual CATT intervals
# excluding zero, since the paper claims per-observation directional evidence
# and that claim carries a multiplicity cost. On the effect cells the endpoints
# are ATT bias, RMSE, coverage and power. Read against TWFE (fixest),
# Callaway-Sant'Anna (did) and Gardner's two-stage (did2s).
#
# What each outcome means. A null rejection rate near 0.05 at the 95% interval
# means the fit can be used for testing and the pre-trend failure was a
# property of the Eq. 2 parameterization rather than of the method. Materially
# above 0.05 means it cannot, and the vignette must say so and recommend the
# conventional estimators for inference. A CATT false-alarm share far above
# 0.05 means the per-observation claim needs a multiplicity warning whatever
# the ATT does.
#
# Usage: Rscript _dev/did-bcf-size.R <cell: A1|A0|B1|B0> <reps> <out.rds>

suppressPackageStartupMessages(library(bartisan))
source("_dev/did/did-bcf-dgp.R")

args <- commandArgs(trailingOnly = TRUE)
cell <- args[1L]; REPS <- as.integer(args[2L]); OUT <- args[3L]
gen <- if (substr(cell, 1L, 1L) == "A") dgp_A else dgp_B
covs <- if (substr(cell, 1L, 1L) == "A") covs_A else covs_B
tau_mult <- if (substr(cell, 2L, 2L) == "0") 0 else 1

fit_did <- function(d) {
  rhs <- paste(c("cohort", "t", covs), collapse = " + ")
  mods <- paste(c(covs, "t", "cohort"), collapse = " + ")
  f <- stats::as.formula(sprintf("y ~ %s + vc(Dit, ~ %s)", rhs, mods))
  suppressMessages(bartisan(f, data = d, family = gaussian(),
    control = bartisan_control(num_trees = c(50L, 25L), num_burn = 400L,
                               num_draws = 600L, chains = 2L, verbose = FALSE)))
}

quiet <- function(e) suppressWarnings(suppressMessages(tryCatch(e, error = function(x) NULL)))

one <- function(r) {
  d <- gen(seed = 5000L + r, tau_mult = tau_mult)
  d$gname <- ifelse(is.finite(d$G), d$G, 0)
  tr <- d$Dit == 1
  att_true <- mean(d$tau_true[tr])

  fit <- fit_did(d)
  cf <- coef(fit, draws = TRUE)[[1L]]
  att <- rowMeans(cf[, tr, drop = FALSE])
  p_bayes <- min(mean(att > 0), mean(att < 0))
  lo <- unname(quantile(att, .025)); hi <- unname(quantile(att, .975))

  # Per-observation CATT intervals among treated rows: the multiplicity check.
  q <- apply(cf[, tr, drop = FALSE], 2L, quantile, probs = c(.025, .975))
  catt_excl <- mean(q[1L, ] > 0 | q[2L, ] < 0)

  out <- data.frame(cell = cell, rep = r, att_true = att_true,
    bart_att = mean(att), bart_lo = lo, bart_hi = hi, bart_p = p_bayes,
    bart_rej95 = as.integer(lo > 0 | hi < 0),
    bart_rej_pb = as.integer(p_bayes < 0.05),
    catt_excl = catt_excl,
    catt_rmse = sqrt(mean((colMeans(cf)[tr] - d$tau_true[tr])^2)),
    stringsAsFactors = FALSE)

  fx <- quiet(fixest::feols(stats::as.formula(sprintf(
    "y ~ Dit + %s | id + t", paste(covs, collapse = " + "))), data = d))
  if (!is.null(fx)) {
    ci <- stats::confint(fx)
    out$twfe_att <- unname(stats::coef(fx)["Dit"])
    out$twfe_rej <- as.integer(ci["Dit", 1] > 0 | ci["Dit", 2] < 0)
  } else { out$twfe_att <- NA_real_; out$twfe_rej <- NA_integer_ }

  # `did` adjusts for covariates at the base period, so handing it covariates
  # that are redrawn independently every period is a misspecification of its
  # own making: on one null replicate that returned -0.63 where the truth is 0.
  # It gets the covariates that are actually fixed per unit -- none in cell A,
  # where assignment is random, and the two statics that drive assignment in
  # cell B.
  cs_x <- if (substr(cell, 1L, 1L) == "A") ~ 1 else ~ x1 + x8
  cs <- quiet({
    a <- did::att_gt(yname = "y", tname = "t", idname = "id", gname = "gname",
                     xformla = cs_x, data = d,
                     control_group = "nevertreated", bstrap = FALSE)
    did::aggte(a, type = "simple", na.rm = TRUE)
  })
  if (!is.null(cs)) {
    out$cs_att <- cs$overall.att
    out$cs_rej <- as.integer(abs(cs$overall.att / cs$overall.se) > 1.96)
  } else { out$cs_att <- NA_real_; out$cs_rej <- NA_integer_ }

  g2 <- quiet(did2s::did2s(data = d, yname = "y",
    first_stage = stats::as.formula(paste("~", paste(covs, collapse = " + "), "| id + t")),
    second_stage = ~ Dit, treatment = "Dit", cluster_var = "id", verbose = FALSE))
  if (!is.null(g2)) {
    ci <- stats::confint(g2)
    out$g2_att <- unname(stats::coef(g2)["Dit"])
    out$g2_rej <- as.integer(ci["Dit", 1] > 0 | ci["Dit", 2] < 0)
  } else { out$g2_att <- NA_real_; out$g2_rej <- NA_integer_ }

  out
}

rows <- list()
for (r in seq_len(REPS)) {
  rows[[length(rows) + 1L]] <- one(r)
  saveRDS(list(res = do.call(rbind, rows), complete = FALSE,
               done = length(rows), total = REPS, cell = cell), OUT)
  cat(sprintf("[%s] %d/%d\n", cell, r, REPS)); utils::flush.console()
}
saveRDS(list(res = do.call(rbind, rows), complete = TRUE,
             done = length(rows), total = REPS, cell = cell), OUT)
