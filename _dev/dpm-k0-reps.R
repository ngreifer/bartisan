# Drawing the Dirichlet process baseline's k_0, and its lambda, on fresh
# replicates and at three sample sizes.
#
# What is tested. Two claims, either of which could be false.
#   (1) `_dev/dpm-k0-draw.R` found, on five replicates per truth, that drawing
#       k_0 under a gamma prior centered on the fixed value raises the coverage
#       of 90% intervals for S(t | x) on a bimodal log-time error (.70 to .81
#       for `dpm_aft()`, .79 to .85 for `dpm_aft_ls()`) at one cost past two
#       SE: `dpm_aft_ls()` on a skewed error, .68 to .65. `_dev/dpm-k0-oracle.R`
#       then found, in a smoke test, that the components' variance prior (lambda,
#       calibrated from the overall residual spread and weighted as nu = 10
#       observations) biases the components' sds upward, and an independent
#       review of the code added that with nu = 10 a small component cannot be
#       narrow, which may be what misses the skewed error's sharp lower edge
#       (one replicate of the oracle study, nu = 3 with both drawn: error
#       density coverage .94). The claim: on twenty fresh replicates, drawing
#       both k_0 and lambda brings bimodal coverage to about nominal for both
#       families, nu = 3 with them brings the skewed truth up too, and no truth
#       loses coverage or accuracy beyond two SE; the skewed cost seen on five
#       replicates under `centered` is chance.
#   (2) Under the fixed priors the bimodal shortfall shrinks as n grows, since
#       a fixed prior's pull on each component is worth a fixed number of
#       observations, while on the skewed error it does not, since a mixture of
#       normals smoothing a sharp edge is a bias that does not shrink as the
#       intervals do.
#
# The measurement. Part 1: the eight truths of `_dev/survival-skew-coverage.R`,
# n = 700, 30% censoring, replicates 6 to 25 (the seeds continue that script's,
# so none of these data sets has been seen), `dpm_aft()` and `dpm_aft_ls()`,
# four arms: `fixed` (the defaults), `centered` (k_0 drawn, gamma prior with
# shape 1.5 and mean at the fixed value), `both` (k_0 and lambda drawn, each
# so) and `both_nu3` (the same with nu = 3). Part 2: the bimodal and skewed
# truths at n = 350 and 1400 (n = 700 is part 1), replicates 6 to 15,
# `dpm_aft()`, `fixed`, `both` and `both_nu3`. Per fit:
# `cover`, the share of (subject, time) cells whose 90% interval holds the true
# S(t | x) (READ THIS); `width`; `s_rmse`; the posterior means of k_0, lambda,
# and the number of clusters. Arms are compared within data set, and the
# change's mean and SE are taken over replicates.
#
# What each outcome would mean.
#   - (1) holds: `both` or `both_nu3`, whichever has no cost, is the candidate
#     default, to be checked on uncensored `dpm()` and `dpm_ls()`
#     (`_dev/dpm-k0-numeric.R`) before anything ships.
#   - `both` lifts bimodal coverage but a truth pays beyond two SE: a trade,
#     and the cost is located before any default changes.
#   - Bimodal coverage stays clearly short under `both` (under .85): the
#     baseline priors are not the whole story, and the code review and the
#     concentration prior come next.
#   - (2) holds: the fixed priors' cost is a small-sample one, which tells how
#     much it matters at the sample sizes people fit; the skewed shortfall is
#     the model's approximation and is documented, not fixed by priors. If the
#     fixed arm's bimodal coverage does not improve with n, the pull is not
#     fading the way the arithmetic says, which would point at the code.
#
# Run with: BARTISAN_LIB=<library holding the working tree> \
#             Rscript _dev/dpm-k0-reps.R [workers]
# Writes:   _dev/dpm-k0-reps.rds, and one file per data set in
#           _dev/dpm-k0-reps/ as it finishes.

ROOT <- normalizePath(".")

# The setup of `_dev/survival-sim.R` plus the skewed truth, in the global
# environment. Each worker runs it too, since the truths are closures over that
# script's helpers, which a worker would not otherwise have.
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
WORKERS <- if (length(args) >= 1L) as.integer(args[1L]) else 6L
OUT <- "_dev/dpm-k0-reps.rds"
UNITS_DIR <- "_dev/dpm-k0-reps"
LEVEL <- 0.90

# The arms name `update_baseline = FALSE` explicitly where they mean the fixed
# baseline, which was the default when this was run and is not since
# 2026-11-03.
ARMS <- list(fixed = list(update_baseline = FALSE),
             centered = list(update_baseline = FALSE, k0_shape = 1.5),
             both = list(k0_shape = 1.5, lambda_shape = 1.5),
             both_nu3 = list(k0_shape = 1.5, lambda_shape = 1.5, nu = 3))

DPM_FAMILIES <- list("dpm_aft()" = dpm_aft, "dpm_aft_ls()" = dpm_aft_ls)

one_dataset <- function(unit) {
  if (!exists("TRUTHS", envir = globalenv(), inherits = FALSE)) {
    setup()
  }

  TRUTHS <- get("TRUTHS", envir = globalenv())

  set.seed(20260828L + 1000L * match(unit$truth, names(TRUTHS)) + unit$rep +
             100L * (unit$n != N_TRAIN) * (1L + (unit$n > N_TRAIN)))

  gen <- TRUTHS[[unit$truth]]
  train <- make_x(unit$n)
  test  <- make_x(N_TEST)
  tr <- gen(train)
  te <- gen(test)

  ctr <- apply_censoring(tr$time, CENSOR)
  train$time <- ctr$time; train$status <- ctr$status
  cte <- apply_censoring(te$time, CENSOR)
  test$time <- cte$time; test$status <- cte$status

  edges <- unique(quantile(tr$time, seq(0.05, 0.95, length.out = N_GRID)))
  strue <- te$surv(edges)

  out <- list()

  for (arm in unit$arms[[1L]]) {
    for (nm in unit$families[[1L]]) {
      fam <- DPM_FAMILIES[[nm]]()
      fam[names(ARMS[[arm]])] <- ARMS[[arm]]

      tick <- proc.time()[["elapsed"]]
      fit <- bartisan(FORM, data = train, family = fam, control = ctrl,
                      verbose = FALSE)
      secs <- proc.time()[["elapsed"]] - tick

      drawn <- predict(fit, newdata = test, type = "survival", times = edges,
                       draws = TRUE)
      shat <- apply(drawn, c(2L, 3L), mean)
      lo <- apply(drawn, c(2L, 3L), stats::quantile, (1 - LEVEL) / 2)
      hi <- apply(drawn, c(2L, 3L), stats::quantile, 1 - (1 - LEVEL) / 2)
      draws <- seq_len(nrow(fit[["aux"]]))

      out[[length(out) + 1L]] <- data.frame(
        part = unit$part, truth = unit$truth, n = unit$n, rep = unit$rep,
        family = nm, arm = arm,
        s_rmse = sqrt(mean((shat - strue)^2)),
        cover = mean(strue >= lo & strue <= hi), width = mean(hi - lo),
        k_0 = mean(bartisan:::dpm_hyper(fit, draws, "k_0")),
        lambda = mean(bartisan:::dpm_hyper(fit, draws, "lambda")),
        clusters = mean(fit[["aux"]][, "clusters"]), secs = secs)
    }
  }

  out <- do.call(rbind, out)
  saveRDS(out, file.path(UNITS_DIR, sprintf("%s-%s-%d-%02d.rds", unit$part,
                                            gsub("[^a-z]", "", unit$truth),
                                            unit$n, unit$rep)))
  out
}

part1 <- expand.grid(rep = 6:25, truth = names(TRUTHS), stringsAsFactors = FALSE)
part1$n <- N_TRAIN
part1$part <- "replicates"
part1$arms <- list(names(ARMS))
part1$families <- list(names(DPM_FAMILIES))

part2 <- expand.grid(rep = 6:15, truth = c("bimodal errors", "skewed errors"),
                     n = c(350L, 1400L), stringsAsFactors = FALSE)
part2$part <- "sample size"
part2$arms <- list(c("fixed", "both", "both_nu3"))
part2$families <- list("dpm_aft()")

units <- rbind(part1, part2[names(part1)])
units <- split(units, seq_len(nrow(units)))

dir.create(UNITS_DIR, showWarnings = FALSE)

pr <- prog_init(total = length(units), title = "Drawing k_0 and lambda, replicates",
                unit = "data set", workers = WORKERS, kind = "simulation")
on.exit(prog_end(pr, "failed", "aborted before the last data set"), add = TRUE)

future::plan(future::multisession, workers = WORKERS)

results <- prog_future_lapply(pr, units, one_dataset,
                              label = function(u) {
                                sprintf("%s n=%d rep %d", u$truth, u$n, u$rep)
                              })

results <- do.call(rbind, Filter(Negate(is.null), results))
saveRDS(list(main = results, complete = TRUE, level = LEVEL,
             when = Sys.time()), OUT)
on.exit()
prog_end(pr, "done", sprintf("%d data sets", length(units)))
cat("wrote", OUT, "\n")
