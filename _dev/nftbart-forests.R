# Which of nftbart's forest settings give it its coverage under crossing
# hazards.
#
# What is tested. In `_dev/nftbart-crossing.R`, on the crossing-hazards truth
# (replicates 26 to 45), `dpm_aft_ls()` with hard rules covered .869 and with
# nftbart's forest settings added .909, a step of +.040 (SE .009) that changed
# four settings at once: k = 5 on the mean forest (`k = c(5, 2)`), the leaf
# scales fixed (`update_sigma_mu = FALSE`), ten trees on the scale
# (`num_trees = c(50, 10)`), and quantile cutpoints. The claim: fixed leaf
# scales carry most of the step, k = 5 matters only once the leaf scale is
# fixed (with the scale drawn, k only centers its prior), and the scale trees
# and cutpoints do nothing. Separately, whether the forest settings need hard
# rules: they are applied under the default soft rules too.
#
# The measurement. The same 20 data sets, seeds and build as
# `_dev/nftbart-crossing.R`, whose fits reproduce exactly (replicate 26 gave
# identical numbers in its smoke test and its full run), so the arms here are
# paired with its stored `hard` and `hard_forests`, which are not refitted.
# Arms, all 2000 + 2000: `add_k5`, `add_fixed`, `add_scale10`, `add_quantile`
# (one setting added to `hard`); `drop_k5`, `drop_fixed`, `drop_scale10`,
# `drop_quantile` (one setting dropped from `hard_forests`); `soft_forests`
# (all four under the default gate). Per fit: `cover` of 90% intervals for
# S(t | x) on 200 test subjects by 20 times (READ THIS), `width`, `s_rmse`.
# Read each setting's step when added and its cost when dropped.
#
# What each outcome would mean.
#   - One setting's added step is most of +.040 and dropping it loses most: it
#     is the cause. If `soft_forests` also reaches about .90, the fix does not
#     need hard rules and that setting is a default candidate, to be checked
#     on the other seven truths of `_dev/nftbart-survival.R`.
#   - Added steps small, dropped costs small, the set's step large: the
#     settings act jointly, and the limitation names the set as the
#     alternative.
#   - `soft_forests` near the default's .82: the settings work only with hard
#     rules.
#
# Run with: BARTISAN_LIB=<library holding bartisan> \
#             Rscript _dev/nftbart-forests.R [workers] [reps]
# Writes:   _dev/nftbart-forests.rds, and one file per data set in
#           _dev/nftbart-forests/ as it finishes.

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

# `setup()` from `_dev/nftbart-survival.R`: the truths, helpers and libraries.
lines <- readLines("_dev/nftbart-survival.R")
eval(parse(text = lines[grep("^ROOT <-", lines):grep("^setup\\(\\)$", lines)]))

args <- commandArgs(TRUE)
WORKERS <- if (length(args) >= 1L) as.integer(args[1L]) else 6L
REPS <- if (length(args) >= 2L) eval(parse(text = args[2L])) else 26:45
OUT <- "_dev/nftbart-forests.rds"
UNITS_DIR <- "_dev/nftbart-forests"
TRUTH <- "crossing hazards"
LEVEL <- 0.90
N_EVAL <- 200L

one_dataset <- function(unit) {
  if (!exists("TRUTHS", envir = globalenv(), inherits = FALSE)) {
    setup()
  }

  TRUTHS <- get("TRUTHS", envir = globalenv())

  set.seed(20260828L + 1000L * match(TRUTH, names(TRUTHS)) + unit$rep)

  gen <- TRUTHS[[TRUTH]]
  train <- make_x(N_TRAIN)
  test  <- make_x(N_TEST)
  tr <- gen(train)
  te <- gen(test)

  ctr <- apply_censoring(tr$time, CENSOR)
  train$time <- ctr$time; train$status <- ctr$status

  edges <- unique(quantile(tr$time, seq(0.05, 0.95, length.out = N_GRID)))
  keep <- seq_len(N_EVAL)
  strue <- te$surv(edges)[keep, , drop = FALSE]
  test <- test[keep, , drop = FALSE]

  hard <- list(num_trees = 50, num_burn = 2000, num_draws = 2000,
               gate = "hard")
  settings <- list(k5 = list(k = c(5, 2)),
                   fixed = list(update_sigma_mu = FALSE),
                   scale10 = list(num_trees = c(50, 10)),
                   quantile = list(x_transform = "quantile"))
  all_four <- do.call(c, unname(settings))

  arms <- c(
    stats::setNames(lapply(settings, function(s) {
      list(control = utils::modifyList(hard, s), family = list())
    }), paste0("add_", names(settings))),
    stats::setNames(lapply(names(settings), function(nm) {
      rest <- do.call(c, unname(settings[names(settings) != nm]))
      list(control = utils::modifyList(hard, rest), family = list())
    }), paste0("drop_", names(settings))),
    list(soft_forests = list(
      control = utils::modifyList(within(hard, rm(gate)), all_four),
      family = list())))

  bart_fit <- function(arm) {
    fam <- bartisan::dpm_aft_ls()
    fam[names(arm$family)] <- arm$family
    ctrl <- do.call(bartisan::bartisan_control, arm$control)
    fit <- bartisan::bartisan(FORM, data = train, family = fam,
                              control = ctrl, verbose = FALSE)
    drawn <- predict(fit, newdata = test, type = "survival", times = edges,
                     draws = TRUE)
    aux <- fit[["aux"]]
    list(drawn = drawn,
         k_0 = if ("k_0" %in% colnames(aux)) mean(aux[, "k_0"]) else NA_real_,
         lambda = if ("lambda" %in% colnames(aux)) mean(aux[, "lambda"]) else NA_real_,
         clusters = mean(aux[, "clusters"]))
  }

  fitters <- lapply(arms, function(arm) function() bart_fit(arm))
  out <- lapply(names(fitters), function(nm) {
    set.seed(unit$rep)
    tick <- proc.time()[["elapsed"]]
    res <- tryCatch(fitters[[nm]](), error = conditionMessage)
    secs <- proc.time()[["elapsed"]] - tick

    if (is.character(res) || anyNA(res$drawn)) {
      msg <- if (is.character(res)) res else "missing values in the draws"
      return(data.frame(rep = unit$rep, arm = nm, cover = NA_real_,
                        width = NA_real_, s_rmse = NA_real_, k_0 = NA_real_,
                        lambda = NA_real_, clusters = NA_real_, secs = secs,
                        error = msg))
    }

    drawn <- res$drawn
    shat <- apply(drawn, c(2L, 3L), mean)
    lo <- apply(drawn, c(2L, 3L), stats::quantile, (1 - LEVEL) / 2)
    hi <- apply(drawn, c(2L, 3L), stats::quantile, 1 - (1 - LEVEL) / 2)

    data.frame(rep = unit$rep, arm = nm,
               cover = mean(strue >= lo & strue <= hi),
               width = mean(hi - lo), s_rmse = sqrt(mean((shat - strue)^2)),
               k_0 = res$k_0, lambda = res$lambda, clusters = res$clusters,
               secs = secs, error = NA_character_)
  })

  out <- do.call(rbind, out)
  saveRDS(out, file.path(UNITS_DIR, sprintf("forests-%02d.rds", unit$rep)))
  out
}

units <- lapply(REPS, function(r) list(rep = r))

dir.create(UNITS_DIR, showWarnings = FALSE)

pr <- prog_init(total = length(units), title = "Crossing hazards, which forest settings",
                unit = "data set", workers = WORKERS, kind = "simulation")
on.exit(prog_end(pr, "failed", "aborted before the last data set"), add = TRUE)

if (WORKERS > 1L) {
  future::plan(future::multisession, workers = WORKERS)
} else {
  future::plan(future::sequential)
}

results <- prog_future_lapply(pr, units, one_dataset,
                              label = function(u) sprintf("rep %d", u$rep))

results <- do.call(rbind, Filter(Negate(is.null), results))
saveRDS(list(main = results, complete = TRUE, level = LEVEL,
             when = Sys.time()), OUT)
on.exit()
prog_end(pr, "done", sprintf("%d data sets", length(units)))
cat("wrote", OUT, "\n")
