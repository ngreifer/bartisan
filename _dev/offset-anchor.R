# Does an anchor that ignores the offset cost a Poisson rate model anything?
#
# Run with: Rscript _dev/offset-anchor.R <mean log exposure> [replicates]
# Writes:   _dev/offset-anchor-mu<mean log exposure>.rds
#
# ---- what is tested ----------------------------------------------------------
#
# `prepare_response()` anchors a Poisson fit at `log(mean(y))` and
# `build_offset()` adds the user's offset on top, so with a log-exposure offset
# the anchor sits about `log(mean(exposure))` above the null fit, and the forest
# has to carry that gap by itself against a leaf prior centered at zero
# (`_dev/TASKS.md`, 2026-11-04, found through the Poisson SBC). The claim is
# that this costs a rate model with realistic exposures something a user would
# see: a worse estimate of the log rate, or intervals that cover worse. The
# suspected route is the drawn leaf scale. A fit with a 5-unit gap drew
# `sigma_mu` at a median of 0.33 against a prior median of 0.21, so the forest
# may loosen its own shrinkage to carry the level and fit the shape with less
# of it.
#
# ---- what is measured --------------------------------------------------------
#
# log rate(x) = b + f(x), with f a fixed function of five uniform predictors
# (an interaction, a linear term and a step; standard deviation about 0.4 on
# the log scale), and exposure lognormal with mean log `MU` and sd 1,
# independent of x. The baseline b is set so the expected count is 4 at every
# `MU`, so the arms differ only in how large the offset is and therefore how far
# the package's anchor is from the null fit: about MU + 0.5, so 0.5, 3 and 5.5
# at MU = 0, 2.5 and 5. n = 500, package defaults otherwise (50 trees, 200 +
# 800 draws, smoothstep gate, `sigma_mu` drawn).
#
# Each data set is fit twice with the same seed, so the chains share a random
# stream and the pair differs only in the anchor: once with the offset as given
# (the package as it is), and once with the anchor moved to the intercept-only
# fit that includes the offset, `log(sum(y) / sum(exposure))`. The second is
# done by adding the difference of the two anchors to the offset, which is what
# the fix in `prepare_response()` would do, without touching the package.
#
# Per fit, against the true log rate b + f(x) at the 500 training points: the
# RMSE of the posterior mean (the decisive column), the same after removing
# each side's mean (shape alone), the mean error (level alone), pointwise 95%
# coverage and mean interval width, and the posterior median of `sigma_mu`.
# Read the paired differences, anchor as is minus corrected, with their
# standard errors over replicates.
#
# ---- what each outcome would mean --------------------------------------------
#
# If the paired differences in RMSE and coverage are within about two standard
# errors of zero at MU = 5, the forest absorbs the gap at no cost to the
# estimate. The fix is then a matter of the documented claim ("the exact
# null-model estimate") and a tidier prior, a small change worth making but
# with no NEWS item about results. If RMSE rises or coverage moves with the
# anchor as is, and more at MU = 5 than at MU = 0, the bug changes results for
# anyone fitting a rate model, the fix goes in with a NEWS item, and the Gaussian
# branch, which ignores the offset the same way, needs the same check. Whichever
# way it comes out, a higher `sigma_mu` with the anchor as is, with no change in
# the estimate, would say the leaf scale is what absorbs the gap.

source(path.expand("~/.claude/skills/live-progress/assets/progress.R"))

suppressPackageStartupMessages(library(bartisan))

args <- commandArgs(trailingOnly = TRUE)
MU <- as.numeric(args[1L])
REPS <- if (length(args) > 1L) as.integer(args[2L]) else 100L
N <- 500L
SD_EXP <- 1
OUT <- sprintf("_dev/offset-anchor-mu%s%s.rds", format(MU),
               if (REPS != 100L) sprintf("-reps%d", REPS) else "")

f_true <- function(x) {
  0.8 * sin(pi * x[, 1L] * x[, 2L]) + 0.6 * (x[, 3L] - 0.5) +
    0.4 * (x[, 4L] > 0.5) - 0.2
}

# The baseline that puts the expected count at 4 for every MU.
B <- log(4) - MU - SD_EXP^2 / 2

score <- function(fit, offset, truth) {
  e <- fit[["eta"]][[1L]] - matrix(offset, nrow(fit[["eta"]][[1L]]), N,
                                   byrow = TRUE)
  est <- colMeans(e)
  lo <- apply(e, 2L, stats::quantile, 0.025, names = FALSE)
  hi <- apply(e, 2L, stats::quantile, 0.975, names = FALSE)

  c(rmse = sqrt(mean((est - truth)^2)),
    rmse_shape = sqrt(mean(((est - mean(est)) - (truth - mean(truth)))^2)),
    level = mean(est - truth),
    coverage = mean(lo <= truth & truth <= hi),
    width = mean(hi - lo),
    sigma_mu = stats::median(fit[["sigma_mu"]]))
}

res <- NULL
START <- 1L

if (file.exists(OUT)) {
  prev <- readRDS(OUT)

  if (!isTRUE(prev$complete) && isTRUE(prev$total == REPS)) {
    res <- prev$res
    START <- max(res$rep) + 1L
    cat(sprintf("resuming at replicate %d of %d\n", START, REPS))
  }
}

todo <- if (START > REPS) integer(0) else seq(START, REPS)

pr <- prog_init(total = max(length(todo), 1L),
                title = sprintf("Offset anchor: mean log exposure %s", format(MU)),
                unit = "replicate", kind = "simulation")
on.exit(prog_end(pr, "failed", "aborted before the last replicate"), add = TRUE)

for (r in todo) {
  t0 <- Sys.time()
  set.seed(7000L + r)

  x <- matrix(stats::runif(N * 5L), N, 5L)
  d <- as.data.frame(x)
  log_exp <- stats::rnorm(N, MU, SD_EXP)
  truth <- B + f_true(x)
  d$y <- stats::rpois(N, exp(log_exp + truth))

  # What `prepare_response()` adds, and what an intercept-only fit with the
  # offset would give.
  as_is <- log(max(mean(d$y), 1 / (2 * N)))
  correct <- log(sum(d$y) / sum(exp(log_exp)))

  fit_seed <- 9000L + r
  set.seed(fit_seed)
  fit_a <- bartisan(y ~ ., data = d, family = stats::poisson(), offset = log_exp)

  set.seed(fit_seed)
  fit_b <- bartisan(y ~ ., data = d, family = stats::poisson(),
                    offset = log_exp + correct - as_is)

  stopifnot(abs(fit_a[["intercept"]] - as_is) < 1e-12,
            abs(fit_b[["intercept"]] - as_is) < 1e-12)

  # Both are scored against the user's offset, since the shift added to the
  # second is part of its anchor, not of the data.
  row <- rbind(data.frame(rep = r, anchor = "as is", gap = as_is - correct,
                          t(score(fit_a, log_exp, truth))),
               data.frame(rep = r, anchor = "corrected", gap = 0,
                          t(score(fit_b, log_exp, truth))))
  res <- rbind(res, row)

  saveRDS(list(res = res, complete = FALSE, total = REPS, mu = MU), OUT)
  prog_tick(pr, label = sprintf("replicate %d", r),
            secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

saveRDS(list(res = res, complete = TRUE, total = REPS, mu = MU), OUT)
on.exit()
prog_end(pr, "done", sprintf("%d replicates", REPS))

# ---- the report --------------------------------------------------------------

a <- res[res$anchor == "as is", ]
b <- res[res$anchor == "corrected", ]
stopifnot(identical(a$rep, b$rep))

cat(sprintf("\nmean log exposure %s, %d replicates, anchor gap %.2f on average\n\n",
            format(MU), nrow(a), mean(a$gap)))
cat(sprintf("%-11s %9s %9s %9s %9s\n", "", "as is", "corrected", "diff", "SE"))

for (m in c("rmse", "rmse_shape", "level", "coverage", "width", "sigma_mu")) {
  d <- a[[m]] - b[[m]]
  cat(sprintf("%-11s %9.4f %9.4f %+9.4f %9.4f\n", m, mean(a[[m]]),
              mean(b[[m]]), mean(d), stats::sd(d) / sqrt(length(d))))
}
