# `dpm_aft_ls()` against nftbart on censored survival data.
#
# What is tested. `dpm_aft_ls()` is documented as the NFT BART model of
# Sparapani et al. (2023), and nftbart (2.3, CRAN) is those authors'
# implementation of it, but the two have only been compared on uncensored
# data (`_dev/nftbart-compare.R`, where `dpm_ls()` was no worse on score or
# RMSE). They fit the same model and differ in the following: bartisan draws
# the mixture weights each iteration (since 2026-11-03), nftbart sets them to
# the components' count shares; bartisan draws k_0 and lambda under gamma
# priors centered on George et al.'s (2019) calibration, nftbart draws k0 and
# b0 under the LIO hyperpriors of Shi et al. (2019); nftbart uses k = 5, ten
# trees on the scale, and recenters the atoms inside the chain. Two claims,
# either of which could be false.
#   (1) `dpm_aft_ls()`'s 90% intervals for S(t | x) cover at least as well as
#       nftbart's on every truth (paired difference above -2 SE), and better
#       on the bimodal error, where fixed weights leave out uncertainty that
#       `_dev/dpm-k0-draw.R` and the weights study found matters.
#   (2) `dpm_aft_ls()`'s S(t | x) RMSE is within two SE of nftbart's or better
#       on every truth.
#
# The measurement. The eight truths of `_dev/dpm-k0-reps.R` (the seven of
# `_dev/survival-sim.R` plus skewed errors), n = 700, 30% censoring,
# replicates 6 to 25 with that script's seeds, so the same data sets. Both
# methods run 2000 burn-in iterations and keep 2000 draws (nftbart's defaults:
# 1000 adaptation and 1000 skipped), 50 trees on the mean, and otherwise each
# package's defaults; nftbart runs single-threaded. Per fit, on the first 200
# test subjects by the 20 times of that study: `cover`, the share of
# (subject, time) cells whose equal-tailed 90% interval from the draws holds
# the true S(t | x) (READ THIS); `cover_inner`, the same over cells whose true
# value is below one (READ THIS on the skewed truth, where a quarter of cells
# sit at one and no normal mixture can reach it); `width`; `s_rmse`, of the
# posterior mean against the truth; `secs`, fit plus prediction. Both
# methods' intervals are computed here from their draws in the same way. The
# comparison is paired within data set, mean and SE over the 20 replicates.
#
# What each outcome would mean.
#   - (1) and (2) hold: `dpm_aft_ls()` matches or beats the reference
#     implementation of its model, and nothing in nftbart needs carrying over.
#     Where both fall short of nominal (expected: spread varies, the skewed
#     edge), the shortfall is the model's, not bartisan's, and the documented
#     limitations stand as written.
#   - nftbart covers better, or is more accurate, beyond 2 SE on some truth:
#     one of its differences is helping there, and a follow-up isolates which
#     (k, the scale forest's size and prior, the LIO priors, recentering)
#     before anything in the package changes.
#
# Run with: BARTISAN_LIB=<library holding bartisan> \
#             Rscript _dev/nftbart-survival.R [workers] [reps]
# Writes:   _dev/nftbart-survival.rds, and one file per data set in
#           _dev/nftbart-survival/ as it finishes.

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

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

  suppressPackageStartupMessages(library(nftbart))
}

setup()

args <- commandArgs(TRUE)
WORKERS <- if (length(args) >= 1L) as.integer(args[1L]) else 6L
REPS <- if (length(args) >= 2L) eval(parse(text = args[2L])) else 6:25
OUT <- "_dev/nftbart-survival.rds"
UNITS_DIR <- "_dev/nftbart-survival"
LEVEL <- 0.90
N_EVAL <- 200L

one_dataset <- function(unit) {
  if (!exists("TRUTHS", envir = globalenv(), inherits = FALSE)) {
    setup()
  }

  TRUTHS <- get("TRUTHS", envir = globalenv())

  # `_dev/dpm-k0-reps.R`'s seed at n = 700, and its draws in its order, so
  # the data sets are that study's.
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
  keep <- seq_len(N_EVAL)
  strue <- te$surv(edges)[keep, , drop = FALSE]
  test <- test[keep, , drop = FALSE]
  inner <- strue < 1

  # Each method's draws as an array, draws x subjects x times.
  fitters <- list(
    "dpm_aft_ls()" = function() {
      ctrl <- bartisan::bartisan_control(num_trees = 50, num_burn = 2000,
                                         num_draws = 2000)
      fit <- bartisan::bartisan(FORM, data = train,
                                family = bartisan::dpm_aft_ls(),
                                control = ctrl, verbose = FALSE)
      predict(fit, newdata = test, type = "survival", times = edges,
              draws = TRUE)
    },
    "nftbart" = function() {
      x <- as.matrix(train[XNAMES])
      utils::capture.output({
        fit <- nftbart::nft(x.train = x, times = train$time,
                            delta = train$status, tc = 1, K = 0,
                            printevery = 100000L)
        # `na.rm = TRUE` only affects nftbart's own summaries: its hazard is
        # 0/0 where S underflows, and the quantile of that stops the call.
        # The draws of S used here have no missing values.
        p <- predict(fit, x.test = as.matrix(test[XNAMES]), tc = 1,
                     events = edges, na.rm = TRUE)
      })
      # Columns are subject-major with time fastest.
      aperm(array(p$surv.test, c(nrow(p$surv.test), length(edges), nrow(test))),
            c(1L, 3L, 2L))
    })

  out <- lapply(names(fitters), function(nm) {
    set.seed(unit$rep)
    tick <- proc.time()[["elapsed"]]
    drawn <- tryCatch(fitters[[nm]](), error = conditionMessage)
    secs <- proc.time()[["elapsed"]] - tick

    if (is.character(drawn) || anyNA(drawn)) {
      msg <- if (is.character(drawn)) drawn else "missing values in the draws"
      return(data.frame(truth = unit$truth, rep = unit$rep, method = nm,
                        cover = NA_real_, cover_inner = NA_real_,
                        width = NA_real_, s_rmse = NA_real_,
                        ndraws = NA_integer_, secs = secs, error = msg))
    }

    shat <- apply(drawn, c(2L, 3L), mean)
    lo <- apply(drawn, c(2L, 3L), stats::quantile, (1 - LEVEL) / 2)
    hi <- apply(drawn, c(2L, 3L), stats::quantile, 1 - (1 - LEVEL) / 2)
    hit <- strue >= lo & strue <= hi

    data.frame(truth = unit$truth, rep = unit$rep, method = nm,
               cover = mean(hit), cover_inner = mean(hit[inner]),
               width = mean(hi - lo), s_rmse = sqrt(mean((shat - strue)^2)),
               ndraws = dim(drawn)[1L], secs = secs, error = NA_character_)
  })

  out <- do.call(rbind, out)
  saveRDS(out, file.path(UNITS_DIR, sprintf("%s-%02d.rds",
                                            gsub("[^a-z]", "", unit$truth),
                                            unit$rep)))
  out
}

units <- expand.grid(rep = REPS, truth = names(TRUTHS), stringsAsFactors = FALSE)
units <- split(units, seq_len(nrow(units)))

dir.create(UNITS_DIR, showWarnings = FALSE)

pr <- prog_init(total = length(units), title = "dpm_aft_ls() against nftbart",
                unit = "data set", workers = WORKERS, kind = "simulation")
on.exit(prog_end(pr, "failed", "aborted before the last data set"), add = TRUE)

if (WORKERS > 1L) {
  future::plan(future::multisession, workers = WORKERS)
} else {
  future::plan(future::sequential)
}

results <- prog_future_lapply(pr, units, one_dataset,
                              label = function(u) sprintf("%s rep %d", u$truth, u$rep))

results <- do.call(rbind, Filter(Negate(is.null), results))
saveRDS(list(main = results, complete = TRUE, level = LEVEL,
             when = Sys.time()), OUT)
on.exit()
prog_end(pr, "done", sprintf("%d data sets", length(units)))
cat("wrote", OUT, "\n")
