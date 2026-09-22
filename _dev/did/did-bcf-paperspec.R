# Does the paper's own specification, or the one the vignette uses, do better --
# and on the paper's own metric, can either reproduce the numbers they report?
#
# The paper's Eq. 3 is mu(D_i, t, X_it) with D_i the *ever-treated indicator*.
# The vignette instead gives mu the adoption cohort. That is a deviation, and it
# was found by comparing specifications on `did::mpdta` rather than derived, so
# it needs a test that could refute it.
#
# The theory that would justify it: the paper's own DGPs give every treated unit
# the same baseline shift (0.75 * D_i) and no cohort-specific level, so mu(D_i,
# t, X) is correctly specified *for their simulation*. Real staggered panels
# generally do have cohort-specific baselines -- the group-time framework of
# Callaway and Sant'Anna assumes each cohort its own level and parallel trend --
# and a mu that cannot tell cohorts apart leaves that difference to tau, which
# is the estimand.
#
# That predicts: the two specifications should perform about the same on the
# paper's DGPs, where no cohort-specific baseline exists, and the cohort
# version should win only where one does. If instead the cohort version also
# wins here, the justification is wrong and the choice was a specification
# search after all.
#
# Metrics are the paper's, computed over treated observations within a
# replicate and reported as mean and sd across replicates. Their MAE and MAPE
# formulas carry a stray square root that their own numbers do not reflect, so
# the standard definitions are used.
#
# Usage: Rscript _dev/did-bcf-paperspec.R <cell A1|B1> <spec ever|cohort> <reps> <out.rds>
suppressPackageStartupMessages(library(bartisan))
setwd("/Users/NoahGreifer/Dropbox/Research/R/bartisan"); source("_dev/did/did-bcf-dgp.R")

a <- commandArgs(trailingOnly = TRUE)
cell <- a[1L]; spec <- a[2L]; REPS <- as.integer(a[3L]); OUT <- a[4L]
gen  <- if (substr(cell, 1L, 1L) == "A") dgp_A else dgp_B
covs <- if (substr(cell, 1L, 1L) == "A") covs_A else covs_B

rows <- list()
for (r in seq_len(REPS)) {
  d <- gen(seed = 7000L + r)
  tr <- d$Dit == 1
  group <- if (spec == "ever") "ever" else "cohort"
  # The paper moderates tau by the covariates and the event time. Event time is
  # a modifier only; the control function has the group and the period.
  f <- stats::as.formula(sprintf("y ~ %s + t + %s + vc(Dit, ~ %s + k_fin)",
        group, paste(covs, collapse = " + "), paste(covs, collapse = " + ")))
  d$k_fin <- ifelse(is.finite(d$k), d$k, -99)
  fit <- suppressMessages(suppressWarnings(bartisan(f, data = d,
    family = gaussian(),
    control = bartisan_control(num_trees = c(50L, 25L), num_burn = 400L,
                               num_draws = 600L, chains = 2L, verbose = FALSE))))
  tau_hat <- colMeans(coef(fit, draws = TRUE)[[1L]])[tr]
  err <- tau_hat - d$tau_true[tr]
  rows[[length(rows) + 1L]] <- data.frame(cell = cell, spec = spec, rep = r,
    rmse = sqrt(mean(err^2)), mae = mean(abs(err)),
    mape = mean(abs(err) / abs(d$tau_true[tr])),
    att_err = mean(tau_hat) - mean(d$tau_true[tr]), stringsAsFactors = FALSE)
  saveRDS(list(res = do.call(rbind, rows), complete = FALSE), OUT)
  cat(sprintf("[%s/%s] %d/%d\n", cell, spec, r, REPS)); utils::flush.console()
}
saveRDS(list(res = do.call(rbind, rows), complete = TRUE), OUT)
