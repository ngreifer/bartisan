# With the mean known, is the mixture's own posterior calibrated once k_0, or
# k_0 and lambda, are drawn?
#
# What is tested. `_dev/dpm-aft-coverage.R` gave `dpm_aft()` the true mean as an
# offset, so that only the error distribution was estimated, and found the 90%
# intervals of `error_density()` holding the true density at 19% of grid points
# on a bimodal error that is itself a two-component normal mixture, which the
# model nests exactly. `_dev/dpm-k0-draw.R` then found that drawing k_0 lifts
# the full model's S(t | x) coverage on that truth from .70 to .81 but not to
# nominal. The claim: with k_0 drawn and the mean known, the mixture's
# posterior is close to calibrated on the bimodal error, so whatever remains in
# the full model comes from the mean forest rather than the mixture. It could
# fail: if `dens_cover` stays well below .90 with k_0 drawn, the mixture is still
# too sure of itself for another reason (the component variances' prior, which
# is calibrated from the overall residual spread with nu = 10 prior
# observations; the concentration's prior; or a fault in the code).
#
# The measurement. Three truths of `_dev/dpm-aft-coverage.R` (log-normal,
# bimodal, skewed errors, the same data and seeds), twenty replicates (the
# first five are that script's), 30% censoring, `dpm_aft()` with the true mean as an
# offset and a forest on x5 alone, which the mean does not use. Four arms:
# `fixed` (the defaults), `centered` (k_0 drawn, gamma prior with shape 1.5
# and mean at the fixed value), `both` (k_0 and lambda drawn, each so), and
# `both_nu3` (the same with nu = 3, BART's own value, in place of 10). The
# fourth was added on an independent review of the code, before this run: with
# nu = 10 a small component carries ten prior observations, so a narrow one at
# the gamma's sharp lower edge cannot form, and in a one-off `dpm()` run nu = 3
# gave the edge a component of sd .19 where nu = 10 gave .37. The third was
# added after a one-replicate smoke test of the first two, before
# this run: on the bimodal truth `centered` raised `dens_cover` only from .30 to
# .40, with the two main components' sd at .385 against a true .35. lambda is
# calibrated from the overall residual spread (1.35 here) and enters each
# component's variance as nu = 10 prior observations at that scale, which
# biases every sd upward; nftbart draws its counterpart (b0). Per fit: `dens_cover`, the share of 101
# points over the middle 96% of the true error distribution whose 90% interval
# from `error_density()` holds the true density (READ THIS FIRST); `cover`, the
# same for S(t | x) on the test set's grid; `z`, the RMS of standardized
# errors of S(t | x) (one when calibrated); `s_rmse`; and the posterior means
# of k_0, the number of clusters, and on the bimodal truth the standard
# deviation of the two largest components (true .35). The points of one
# replicate's grid miss together (in the smoke test's replicate the sample put
# 54% of its errors in the upper mode, a 2.1-SE fluctuation, and centering moved
# both modes by .13), so coverage rests on about twenty independent outcomes per
# truth, and its standard error is computed over replicates.
#
# What each outcome would mean.
#   - Bimodal `dens_cover` near .90 (about .80 or more) under `both`, and
#     `top_sd` near .35: the fixed baseline priors were the mixture's whole
#     problem on this truth, and the full model's remaining gap is the mean
#     forest's. Under `centered` alone the smoke test predicts little gain.
#   - Still well below under `both` (under .6): the mixture is overconfident for
#     a further reason; the code review and the concentration's prior come
#     next.
#   - Log-normal should cover near nominal in every arm; if a drawn arm drops
#     it, the looser baseline costs on the case the default must get right.
#   - Skewed: if `both_nu3` raises its `dens_cover` clearly above the other
#     arms, the skewed shortfall is the variance prior's floor, as the review
#     argues, and not only the approximation of a sharp edge by normals; if it
#     stays low in every arm, it is the approximation.
#   - `both_nu3` against `both` on the bimodal truth says whether nu matters
#     once lambda is drawn; its large clusters should make little of it.
#
# Run with: BARTISAN_LIB=<library holding the working tree> \
#             Rscript _dev/dpm-k0-oracle.R
# Writes:   _dev/dpm-k0-oracle.rds

lib <- Sys.getenv("BARTISAN_LIB")

if (nzchar(lib)) {
  .libPaths(c(lib, .libPaths()))
}

lines <- readLines("_dev/survival-sim.R")
eval(parse(text = lines[seq_len(grep("^# ---- one replicate", lines) - 1L)]))

OUT <- "_dev/dpm-k0-oracle.rds"
N_REP <- as.integer(Sys.getenv("SIM_REPS", 20L))
LEVEL <- 0.90

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

TRUTHS[["skewed errors"]] <- truth_skewed

error_truth <- list(
  "log-normal errors" = list(
    dens = function(w) dnorm(w, 0, 0.7),
    range = qnorm(c(0.02, 0.98), 0, 0.7)),
  "bimodal errors" = list(
    dens = function(w) 0.5 * dnorm(w, 1.3, 0.35) + 0.5 * dnorm(w, -1.3, 0.35),
    range = c(-1.3 - 2.05 * 0.35, 1.3 + 2.05 * 0.35)),
  "skewed errors" = list(
    dens = function(w) {
      dgamma(1.5 + sqrt(1.5) * w / 0.7, 1.5) * sqrt(1.5) / 0.7
    },
    range = 0.7 * (qgamma(c(0.02, 0.98), 1.5) - 1.5) / sqrt(1.5))
)

STUDY <- c("log-normal errors", "bimodal errors", "skewed errors")
# The arms name `update_baseline = FALSE` explicitly where they mean the fixed
# baseline, which was the default when this was run and is not since
# 2026-11-03.
ARMS <- list(fixed = list(update_baseline = FALSE),
             centered = list(update_baseline = FALSE, k0_shape = 1.5),
             both = list(k0_shape = 1.5, lambda_shape = 1.5),
             both_nu3 = list(k0_shape = 1.5, lambda_shape = 1.5, nu = 3))

# The standard deviations of the two largest components at each draw,
# averaged over draws: the bimodal truth's are .35.
top_sd <- function(fit) {
  draws <- seq_len(nrow(fit[["aux"]]))
  sds <- vapply(draws, function(s) {
    m <- bartisan:::mixture_at(fit, s)
    m <- m[order(-m[, "weight"]), , drop = FALSE]
    mean(m[seq_len(min(2L, nrow(m))), "sd"])
  }, numeric(1))
  mean(sds)
}

one_rep <- function(truth_name, rep) {
  out <- list()

  for (arm in names(ARMS)) {
    set.seed(20260828L + 1000L * match(truth_name, names(TRUTHS)) + rep)

    gen <- TRUTHS[[truth_name]]
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

    fam <- dpm_aft()
    fam[names(ARMS[[arm]])] <- ARMS[[arm]]
    label <- sprintf("%s / %s rep %d", truth_name, arm, rep)

    tick <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      bartisan(Surv(time, status) ~ x5, data = train, family = fam,
               offset = signal(train), control = ctrl, verbose = FALSE),
      error = function(e) e)
    secs <- proc.time()[["elapsed"]] - tick

    if (inherits(fit, "error")) {
      prog_tick(pr, ok = FALSE, secs = secs, msg = conditionMessage(fit),
                label = label)
      next
    }

    drawn <- predict(fit, newdata = test, type = "survival", times = edges,
                     draws = TRUE, offset = signal(test))
    shat <- apply(drawn, c(2L, 3L), mean)
    sdev <- apply(drawn, c(2L, 3L), stats::sd)
    lo <- apply(drawn, c(2L, 3L), stats::quantile, (1 - LEVEL) / 2)
    hi <- apply(drawn, c(2L, 3L), stats::quantile, 1 - (1 - LEVEL) / 2)

    et <- error_truth[[truth_name]]
    at <- seq(et$range[1L], et$range[2L], length.out = 101L)
    ed <- error_density(fit, at = at, level = LEVEL)
    aux <- fit[["aux"]]

    out[[arm]] <- data.frame(
      truth = truth_name, arm = arm, rep = rep,
      dens_cover = mean(et$dens(at) >= ed$lower & et$dens(at) <= ed$upper),
      cover = mean(strue >= lo & strue <= hi),
      z = sqrt(mean(((shat - strue) / pmax(sdev, 1e-8))^2)),
      s_rmse = sqrt(mean((shat - strue)^2)),
      k_0 = mean(bartisan:::dpm_hyper(fit, seq_len(nrow(aux)), "k_0")),
      lambda = mean(bartisan:::dpm_hyper(fit, seq_len(nrow(aux)), "lambda")),
      clusters = mean(aux[, "clusters"]), top_sd = top_sd(fit), secs = secs)

    prog_tick(pr, secs = secs, label = label)
  }

  do.call(rbind, out)
}

N_FITS <- length(STUDY) * N_REP * length(ARMS)

pr <- prog_init(total = N_FITS, title = "Drawing k_0, mean known",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed", "aborted before the last replicate"),
        add = TRUE)

main <- list()

checkpoint <- function(complete) {
  saveRDS(list(main = do.call(rbind, main), complete = complete,
               level = LEVEL, when = Sys.time()), OUT)
}

for (nm in STUDY) {
  for (r in seq_len(N_REP)) {
    main[[length(main) + 1L]] <- one_rep(nm, r)
    checkpoint(FALSE)
  }
}

checkpoint(TRUE)
on.exit()
prog_end(pr, "done", sprintf("%d replicates", length(main)))
cat("wrote", OUT, "\n")
