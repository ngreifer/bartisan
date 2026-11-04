# Does `dpm_aft()` cover S(t | x) at least as well as `lognormal_aft()` on every
# survival truth?
#
# What is tested. `dpm_aft()` is the default family for a survival response.
# The maintainer's rule (2026-11-03): it stays the default as long as its
# coverage is at least that of its parametric counterpart, `lognormal_aft()`,
# everywhere, and where both fall short the shortfall is documented as a
# limitation with alternatives. The same is asked of `dpm_aft_ls()` against
# `lognormal_aft_ls()`. On five replicates (`_dev/survival-skew-coverage.R`,
# rerun with the drawn mixture weights) `dpm_aft()` covered more on five of
# the eight truths and less on three: crossing hazards (.57 against .66), the
# log-normal truth, where `lognormal_aft()` is the true model (.93 against .95),
# and the spread that varies (.50 against .52). The claim: on twenty-five
# replicates `dpm_aft()`'s coverage is at least `lognormal_aft()`'s, within two
# SE, on every truth but the log-normal one, where the true model may cover
# slightly better. It fails on any other truth where `lognormal_aft()` covers
# better by more than two SE.
#
# The measurement. The eight truths of `_dev/survival-skew-coverage.R`, n =
# 700, 30% censoring, replicates 6 to 25 of `_dev/dpm-k0-reps.R` (the same
# seeds and the same draws, so these fits pair with that study's by
# replicate), `lognormal_aft()` and `lognormal_aft_ls()` at their defaults,
# 50 trees, 500 + 500. Replicates 1 to 5 come from
# `_dev/survival-skew-coverage.rds`, fitted on a build whose two log-normal
# families are unchanged since. Per fit: `cover`, the share of (subject, time)
# cells whose 90% interval holds the true S(t | x) (READ THIS), `width`,
# `s_rmse`, and `cover_inner`, the coverage over the cells whose true S(t | x)
# is below one. That column was added before this run, on finding that 25% of
# the skewed truth's cells lie below its gamma error's lower edge, where the
# true S(t | x) is exactly one and no family with normal or log-normal
# components can cover it; on that truth read `cover_inner`. The comparison is
# the paired difference in `cover` against each Dirichlet process family's
# default arm, with its SE over replicates.
#
# What each outcome would mean.
#   - The claim holds: the defaults stay, and the truths where both fall short
#     (bimodal and skewed errors, crossing hazards, and for the single-scale
#     families a spread that varies) are written up as limitations, with the
#     location-scale families and the discrete-time route as the alternatives
#     where they cover better.
#   - `lognormal_aft()` covers better on some truth other than the log-normal
#     one: the default is reconsidered for that case, and what the drawn
#     baseline (`_dev/dpm-k0-reps.R`) does there decides whether a default
#     change fixes it.
#
# Run with: BARTISAN_LIB=<library holding the working tree> \
#             Rscript _dev/aft-lognormal-reps.R [workers]
# Writes:   _dev/aft-lognormal-reps.rds, and one file per data set in
#           _dev/aft-lognormal-reps/ as it finishes.

ROOT <- normalizePath(".")

# The setup of `_dev/survival-sim.R` plus the skewed truth, as in
# `_dev/dpm-k0-reps.R`; each worker runs it too.
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
OUT <- "_dev/aft-lognormal-reps.rds"
UNITS_DIR <- "_dev/aft-lognormal-reps"
LEVEL <- 0.90

PARAMETRIC <- list("lognormal_aft()" = lognormal_aft,
                   "lognormal_aft_ls()" = lognormal_aft_ls)

one_dataset <- function(unit) {
  if (!exists("TRUTHS", envir = globalenv(), inherits = FALSE)) {
    setup()
  }

  TRUTHS <- get("TRUTHS", envir = globalenv())

  # The seed and the order of draws are `_dev/dpm-k0-reps.R`'s at n = 700.
  set.seed(20260828L + 1000L * match(unit$truth, names(TRUTHS)) + unit$rep)

  gen <- TRUTHS[[unit$truth]]
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

  out <- lapply(names(PARAMETRIC), function(nm) {
    tick <- proc.time()[["elapsed"]]
    fit <- bartisan(FORM, data = train, family = PARAMETRIC[[nm]](),
                    control = ctrl, verbose = FALSE)
    secs <- proc.time()[["elapsed"]] - tick

    drawn <- predict(fit, newdata = test, type = "survival", times = edges,
                     draws = TRUE)
    shat <- apply(drawn, c(2L, 3L), mean)
    lo <- apply(drawn, c(2L, 3L), stats::quantile, (1 - LEVEL) / 2)
    hi <- apply(drawn, c(2L, 3L), stats::quantile, 1 - (1 - LEVEL) / 2)

    hit <- strue >= lo & strue <= hi
    inner <- strue < 1 - 1e-9

    data.frame(truth = unit$truth, rep = unit$rep, family = nm,
               s_rmse = sqrt(mean((shat - strue)^2)),
               cover = mean(hit), cover_inner = mean(hit[inner]),
               width = mean(hi - lo), secs = secs)
  })

  out <- do.call(rbind, out)
  saveRDS(out, file.path(UNITS_DIR, sprintf("%s-%02d.rds",
                                            gsub("[^a-z]", "", unit$truth),
                                            unit$rep)))
  out
}

units <- expand.grid(rep = 6:25, truth = names(TRUTHS), stringsAsFactors = FALSE)
units <- split(units, seq_len(nrow(units)))

dir.create(UNITS_DIR, showWarnings = FALSE)

pr <- prog_init(total = length(units), title = "Log-normal AFT, replicates",
                unit = "data set", workers = WORKERS, kind = "simulation")
on.exit(prog_end(pr, "failed", "aborted before the last data set"), add = TRUE)

future::plan(future::multisession, workers = WORKERS)

results <- prog_future_lapply(pr, units, one_dataset,
                              label = function(u) sprintf("%s rep %d", u$truth, u$rep))

results <- do.call(rbind, Filter(Negate(is.null), results))
saveRDS(list(main = results, complete = TRUE, level = LEVEL,
             when = Sys.time()), OUT)
on.exit()
prog_end(pr, "done", sprintf("%d data sets", length(units)))
cat("wrote", OUT, "\n")
