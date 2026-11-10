# What gives nftbart its coverage under crossing hazards.
#
# What is tested. In `_dev/nftbart-survival.R`, on the crossing-hazards truth,
# `dpm_aft_ls()`'s 90% intervals for S(t | x) covered .838 against nftbart's
# .896 (paired difference -.058, SE .007, 20 data sets), with intervals 47%
# narrower and a smaller RMSE. The claim: the gap comes from settings that
# differ between the two rather than from anything else in their samplers, and
# moving `dpm_aft_ls()` onto nftbart's settings one group at a time closes it.
# It could be false in two ways: the gap may not reproduce on fresh data sets,
# or it may survive with every setting that can be matched matched.
#
# nftbart's settings, as they map onto bartisan's (read from nftbart 2.3's
# `nft2()` and `DPMLIOneal8.h`):
#   - hard rules, `gate = "hard"`;
#   - 50 trees on the mean and 10 on the scale, `num_trees = c(50, 10)`;
#   - k = 5 on the mean forest, where both packages follow BART's convention
#     (leaf sd of the sum half the response's range over k), `k = c(5, 2)`,
#     with the leaf scales fixed, `update_sigma_mu = FALSE`;
#   - cutpoints at the predictors' sorted unique values thinned to 100, closest
#     to `x_transform = "quantile"`;
#   - DPM precision tau ~ Gamma(a0 = 3, b0), b0 ~ Gamma(2, 1), and
#     mu | tau ~ N(0, 1 / (k0 tau)), k0 ~ Gamma(1.5, 7.5), on residuals
#     standardized by a scale centered on s, the scale of a Weibull `survreg()`
#     fit. bartisan's mixture is on log time, sigma^2 ~ nu lambda / chisq_nu,
#     so nu = 2 a0 = 6 and lambda = s^2 b0 / 3 ~ Gamma(2, rate 3 / s^2); k0 is
#     scale-free and carries over.
# Not matched, because bartisan has no setting for them: nftbart's Gamma(1, 0.1)
# prior on the concentration (bartisan's is Rossi's taper), its product-of-trees
# prior on the scale, its recentering of the atoms inside the chain, and its
# count-share mixture weights.
#
# The measurement. The crossing-hazards truth of `_dev/survival-sim.R`, n =
# 700, 30% censoring, replicates 26 to 45 with `_dev/dpm-k0-reps.R`'s seed
# formula, which no study has fitted, so whether the earlier gap is
# seed-dependent shows directly. Arms, each adding to the one before:
# `default` (`dpm_aft_ls()` as shipped), `hard` (hard rules), `hard_forests`
# (nftbart's forests: trees, k, fixed leaf scales, cutpoints), `nft` (its DPM
# hyperpriors as well), and nftbart itself at its defaults. All run 2000
# burn-in iterations and keep 2000 draws. Per fit, on 200 test subjects by 20
# times: `cover`, the share of cells whose equal-tailed 90% interval holds the
# true S(t | x) (READ THIS); `width`; `s_rmse`. Paired differences against
# nftbart within data set, mean and SE over the 20.
#
# What each outcome would mean.
#   - `default` within two SE of nftbart: the earlier -.058 was particular to
#     those 20 data sets, and there is nothing to fix.
#   - The gap closes along the arms: the arm at which it closes names the
#     setting that matters, and that setting is checked for cost on the other
#     seven truths before any default changes.
#   - The gap survives `nft`: it comes from what cannot be matched (recentering,
#     the concentration prior, the weights, the scale prior), and the shortfall
#     is documented as a limitation of `dpm_aft_ls()` on non-location-scale
#     truths.
#
# Run with: BARTISAN_LIB=<library holding bartisan> \
#             Rscript _dev/nftbart-crossing.R [workers] [reps]
# Writes:   _dev/nftbart-crossing.rds, and one file per data set in
#           _dev/nftbart-crossing/ as it finishes.

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

# `setup()` from `_dev/nftbart-survival.R`: the truths, helpers and libraries.
lines <- readLines("_dev/nftbart-survival.R")
eval(parse(text = lines[grep("^ROOT <-", lines):grep("^setup\\(\\)$", lines)]))

args <- commandArgs(TRUE)
WORKERS <- if (length(args) >= 1L) as.integer(args[1L]) else 6L
REPS <- if (length(args) >= 2L) eval(parse(text = args[2L])) else 26:45
OUT <- "_dev/nftbart-crossing.rds"
UNITS_DIR <- "_dev/nftbart-crossing"
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

  # nftbart's standardizing scale, from the same Weibull fit it makes.
  s2 <- survival::survreg(FORM, data = train, dist = "weibull")$scale^2

  base <- list(num_trees = 50, num_burn = 2000, num_draws = 2000)
  forests <- list(gate = "hard", num_trees = c(50, 10), k = c(5, 2),
                  update_sigma_mu = FALSE, x_transform = "quantile")
  hyper <- list(nu = 6, k0_shape = 1.5, k0_rate = 7.5,
                lambda_shape = 2, lambda_rate = 3 / s2)

  arms <- list(
    default = list(control = base, family = list()),
    hard = list(control = utils::modifyList(base, list(gate = "hard")),
                family = list()),
    hard_forests = list(control = utils::modifyList(base, forests),
                        family = list()),
    nft = list(control = utils::modifyList(base, forests), family = hyper))

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
  fitters[["nftbart"]] <- function() {
    utils::capture.output({
      fit <- nftbart::nft(x.train = as.matrix(train[XNAMES]),
                          times = train$time, delta = train$status, tc = 1,
                          K = 0, printevery = 100000L)
      p <- predict(fit, x.test = as.matrix(test[XNAMES]), tc = 1,
                   events = edges, na.rm = TRUE)
    })
    list(drawn = aperm(array(p$surv.test, c(nrow(p$surv.test), length(edges),
                                            nrow(test))), c(1L, 3L, 2L)),
         k_0 = if (is.null(fit$k0)) NA_real_ else mean(fit$k0),
         lambda = NA_real_,
         clusters = mean(apply(fit$dpC, 1L, function(v) length(unique(v)))))
  }

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
  saveRDS(out, file.path(UNITS_DIR, sprintf("crossing-%02d.rds", unit$rep)))
  out
}

units <- lapply(REPS, function(r) list(rep = r))

dir.create(UNITS_DIR, showWarnings = FALSE)

pr <- prog_init(total = length(units), title = "Crossing hazards, nftbart's settings",
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
