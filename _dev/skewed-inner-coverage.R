# The skewed survival truth, scored away from the cells its error's lower edge
# makes uncoverable.
#
# What is tested. The skewed truth of `_dev/survival-skew-coverage.R` has
# log T = m(x) + 0.7 W with W a standardized gamma(1.5), whose support ends
# below at -0.7 sqrt(1.5). On the evaluation grid, 25% of the (subject, time)
# cells lie below that edge (measured over twenty replicates, 2026-11-03),
# where the true S(t | x) is exactly one; a model built from normal or
# log-normal components gives a value below one in every draw there, so those
# cells are missed by every family alike. Every skewed-truth coverage reported
# so far (about .67 for the Dirichlet process families) counted them. The
# claim: over the remaining cells, the cells with a true S(t | x) below one,
# `dpm_aft()` and `dpm_aft_ls()` cover near nominal once k_0 and lambda are
# drawn (`_dev/dpm-k0-oracle.R` found their error density covering .89 to .91
# there with the mean known), and the default arm covers less. It fails if the
# interior coverage is low in every arm, which would make the skewed shortfall
# a real calibration problem rather than a property of the truth's bounded
# support.
#
# The measurement. Two parts on the skewed truth's own data and seeds. Part
# `replicates`: n = 700, 30% censoring, replicates 6 to 25 of
# `_dev/dpm-k0-reps.R` (the same data sets), `dpm_aft()` and `dpm_aft_ls()`
# under `fixed`, `both` (k_0 and lambda drawn, gamma priors of shape 1.5
# centered on the fixed values) and `both_nu3` (the same, nu = 3). Part
# `oracle`: replicates 1 to 20 of `_dev/dpm-k0-oracle.R`, the true mean as an
# offset with a forest on x5, `dpm_aft()` under the same three arms. Per fit:
# `cover_inner`, the coverage of 90% intervals for S(t | x) over the cells
# whose true value is below one (READ THIS); `cover`, over all cells, which
# should reproduce the earlier studies' level; `edge_share`, the share of cells
# at one; `width`; `s_rmse`.
#
# What each outcome would mean.
#   - Interior coverage near nominal (about .85 or more) under the drawn arms:
#     the skewed truth's shortfall was its bounded support plus, at the
#     defaults, the variance prior's floor. The documentation can then say
#     precisely that intervals for S(t | x) cannot reach one before an error's
#     hard lower edge, and the comparison with `lognormal_aft()` uses
#     `cover_inner` on this truth.
#   - Interior coverage low in every arm: a calibration shortfall on skewed
#     errors that no prior here fixes, written up as a limitation.
#
# Run with: BARTISAN_LIB=<library holding the working tree> \
#             Rscript _dev/skewed-inner-coverage.R [workers]
# Writes:   _dev/skewed-inner-coverage.rds, and one file per data set in
#           _dev/skewed-inner-coverage/ as it finishes.

ROOT <- normalizePath(".")

# As in `_dev/dpm-k0-reps.R`; each worker runs it too.
setup <- function() {
  lib <- Sys.getenv("BARTISAN_LIB")

  if (nzchar(lib)) {
    .libPaths(c(lib, .libPaths()))
  }

  lines <- readLines(file.path(ROOT, "_dev/survival-sim.R"))
  eval(parse(text = lines[seq_len(grep("^# ---- one replicate", lines) - 1L)]),
       envir = globalenv())

  truth_skewed <- function(d, sigma = 0.7, shape = 1.5) {
    m <- signal(d)
    w <- (rgamma(nrow(d), shape) - shape) / sqrt(shape)
    list(time = exp(m + sigma * w),
         surv = function(tt) {
           1 - outer(m, log(tt), function(mm, lt) {
             pgamma(shape + sqrt(shape) * (lt - mm) / sigma, shape)
           })
         })
  }

  truths <- get("TRUTHS", envir = globalenv())
  truths[["skewed errors"]] <- truth_skewed
  assign("TRUTHS", truths, envir = globalenv())
}

setup()

args <- commandArgs(TRUE)
WORKERS <- if (length(args) >= 1L) as.integer(args[1L]) else 4L
OUT <- "_dev/skewed-inner-coverage.rds"
UNITS_DIR <- "_dev/skewed-inner-coverage"
LEVEL <- 0.90

# The arms name `update_baseline = FALSE` explicitly where they mean the fixed
# baseline, which was the default when this was run and is not since
# 2026-11-03.
ARMS <- list(fixed = list(update_baseline = FALSE),
             both = list(k0_shape = 1.5, lambda_shape = 1.5),
             both_nu3 = list(k0_shape = 1.5, lambda_shape = 1.5, nu = 3))

one_dataset <- function(unit) {
  if (!exists("TRUTHS", envir = globalenv(), inherits = FALSE)) {
    setup()
  }

  TRUTHS <- get("TRUTHS", envir = globalenv())
  truth <- "skewed errors"

  # The seed and the order of draws are those of `_dev/dpm-k0-reps.R` at
  # n = 700 and of `_dev/dpm-k0-oracle.R`, which agree.
  set.seed(20260828L + 1000L * match(truth, names(TRUTHS)) + unit$rep)

  gen <- TRUTHS[[truth]]
  train <- make_x(N_TRAIN)
  test  <- make_x(N_TEST)
  tr <- gen(train)
  te <- gen(test)

  ctr <- apply_censoring(tr$time, CENSOR)
  train$time <- ctr$time; train$status <- ctr$status
  cte <- apply_censoring(te$time, CENSOR)
  test$time <- cte$time; test$status <- cte$status

  edges <- unique(quantile(tr$time, seq(0.05, 0.95, length.out = N_GRID)))
  strue <- te$surv(edges)
  inner <- strue < 1 - 1e-9

  families <- if (unit$part == "oracle") "dpm_aft()" else c("dpm_aft()", "dpm_aft_ls()")
  out <- list()

  for (arm in names(ARMS)) {
    for (nm in families) {
      fam <- if (nm == "dpm_aft()") dpm_aft() else dpm_aft_ls()
      fam[names(ARMS[[arm]])] <- ARMS[[arm]]

      tick <- proc.time()[["elapsed"]]

      if (unit$part == "oracle") {
        fit <- bartisan(Surv(time, status) ~ x5, data = train, family = fam,
                        offset = signal(train), control = ctrl, verbose = FALSE)
        drawn <- predict(fit, newdata = test, type = "survival", times = edges,
                         draws = TRUE, offset = signal(test))
      }
      else {
        fit <- bartisan(FORM, data = train, family = fam, control = ctrl,
                        verbose = FALSE)
        drawn <- predict(fit, newdata = test, type = "survival", times = edges,
                         draws = TRUE)
      }

      secs <- proc.time()[["elapsed"]] - tick
      shat <- apply(drawn, c(2L, 3L), mean)
      lo <- apply(drawn, c(2L, 3L), stats::quantile, (1 - LEVEL) / 2)
      hi <- apply(drawn, c(2L, 3L), stats::quantile, 1 - (1 - LEVEL) / 2)
      hit <- strue >= lo & strue <= hi

      out[[length(out) + 1L]] <- data.frame(
        part = unit$part, rep = unit$rep, family = nm, arm = arm,
        cover = mean(hit), cover_inner = mean(hit[inner]),
        edge_share = mean(!inner), width = mean(hi - lo),
        s_rmse = sqrt(mean((shat - strue)^2)), secs = secs)
    }
  }

  out <- do.call(rbind, out)
  saveRDS(out, file.path(UNITS_DIR, sprintf("%s-%02d.rds", unit$part,
                                            unit$rep)))
  out
}

units <- rbind(data.frame(part = "replicates", rep = 6:25),
               data.frame(part = "oracle", rep = 1:20))
units <- split(units, seq_len(nrow(units)))

dir.create(UNITS_DIR, showWarnings = FALSE)

pr <- prog_init(total = length(units), title = "Skewed truth, interior cells",
                unit = "data set", workers = WORKERS, kind = "simulation")
on.exit(prog_end(pr, "failed", "aborted before the last data set"), add = TRUE)

future::plan(future::multisession, workers = WORKERS)

results <- prog_future_lapply(pr, units, one_dataset,
                              label = function(u) sprintf("%s rep %d", u$part, u$rep))

results <- do.call(rbind, Filter(Negate(is.null), results))
saveRDS(list(main = results, complete = TRUE, level = LEVEL,
             when = Sys.time()), OUT)
on.exit()
prog_end(pr, "done", sprintf("%d data sets", length(units)))
cat("wrote", OUT, "\n")
