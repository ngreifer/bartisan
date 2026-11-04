# Drawing the Dirichlet process baseline's k_0 rather than fixing it, and the
# coverage of S(t | x) intervals.
#
# What is tested. With the mixture's weights drawn (2026-11-03), 90% intervals
# for S(t | x) from `dpm_aft()` and `dpm_aft_ls()` cover .70 and .79 on the
# bimodal log-time error of `_dev/survival-sim.R`, against .88 to .94 on the
# accelerated failure time truths whose error they also nest.
# `_dev/dpm-aft-coverage.R` traced part of the shortfall to the baseline's k_0:
# George et al. (2019) fix it so that the largest linear-fit residual sits ten
# prior scale units out, which pulls the component means toward zero and, through
# the conjugate baseline, widens the components that sit far out (fitted sds .41
# to .48 against a true .35). nftbart draws k_0 instead, under a gamma prior
# (Shi et al. 2019). The claim: drawing k_0 lets the data loosen that pull,
# which raises bimodal coverage toward nominal without costing coverage or point
# accuracy on the other truths. It could fail either way: the bimodal shortfall
# may not be the pull, or a looser baseline may cost elsewhere, for instance by
# letting the mixture's location trade against the mean forest's level.
#
# The measurement. The eight truths of `_dev/survival-skew-coverage.R` (its
# data and seeds), five replicates, `dpm_aft()` and `dpm_aft_ls()` with k_0
# drawn under two priors: `centered`, gamma with shape 1.5 and mean at the
# paper's fixed value, which keeps that calibration as the prior's center; and
# `lio`, nftbart's gamma(1.5, rate 7.5), mean .2, far looser than the fixed
# value (4.6 on a normal-error check data set of 200). The fixed arm is
# `_dev/survival-skew-coverage.rds`, the same data on the same build (default
# fits are bit-identical with the option off), paired by replicate. Per fit:
# `cover`, the share of (subject, time) cells whose 90% interval holds the true
# S(t | x); `width`; `s_rmse`; and the posterior mean of k_0 and of the number
# of clusters. Read `cover` on the bimodal truth first, then the paired change
# in `s_rmse` and `cover` on every other truth.
#
# What each outcome would mean.
#   - Bimodal coverage reaches about .85 or more under a drawn arm, with no
#     other truth losing coverage or accuracy beyond two SE: the pull was what
#     remained, and drawing k_0 is the candidate default. Adopting it means
#     rerunning the families and survival simulations the vignettes quote and a
#     NEWS item, since every Dirichlet process fit would change.
#   - Bimodal coverage rises but another truth pays: a trade, reported as such;
#     the fixed default stays and drawing k_0 is at most an option.
#   - Bimodal coverage does not rise, though the posterior of k_0 moves well
#     below the fixed value: the pull is not what remains, and the shortfall is
#     elsewhere (the concentration's prior, or the forest).
#   - If k_0 does not move off the fixed value under `centered`, the data do not
#     inform it at five or so components and only the prior's choice matters.
#   The skewed truth should not move under either arm: its misses come from a
#     mixture of normals smoothing a sharp edge, not from the baseline.
#
# Run with: BARTISAN_LIB=<library holding the working tree> \
#             Rscript _dev/dpm-k0-draw.R
# Writes:   _dev/dpm-k0-draw.rds

lib <- Sys.getenv("BARTISAN_LIB")

if (nzchar(lib)) {
  .libPaths(c(lib, .libPaths()))
}

lines <- readLines("_dev/survival-sim.R")
eval(parse(text = lines[seq_len(grep("^# ---- one replicate", lines) - 1L)]))

OUT <- "_dev/dpm-k0-draw.rds"
N_REP <- as.integer(Sys.getenv("SIM_REPS", N_REP))
LEVEL <- 0.90

# As in `_dev/survival-skew-coverage.R`, where it is the eighth truth; its
# seeds follow from that position.
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

# `k0_rate = NULL` centers the prior on the paper's fixed k_0.
ARMS <- list(centered = list(k0_shape = 1.5, k0_rate = NULL),
             lio = list(k0_shape = 1.5, k0_rate = 7.5))

DPM_FAMILIES <- list("dpm_aft()" = dpm_aft, "dpm_aft_ls()" = dpm_aft_ls)

one_rep <- function(truth_name, rep) {
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

  out <- list()

  for (arm in names(ARMS)) {
    for (nm in names(DPM_FAMILIES)) {
      fam <- DPM_FAMILIES[[nm]]()
      # The fixed baseline was the default when this was run; since
      # 2026-11-03 lambda is drawn unless this says otherwise.
      fam$update_baseline <- FALSE
      fam$k0_shape <- ARMS[[arm]]$k0_shape
      fam$k0_rate <- ARMS[[arm]]$k0_rate
      label <- sprintf("%s / %s %s rep %d", truth_name, nm, arm, rep)

      tick <- proc.time()[["elapsed"]]
      fit <- tryCatch(bartisan(FORM, data = train, family = fam,
                               control = ctrl, verbose = FALSE),
                      error = function(e) e)
      secs <- proc.time()[["elapsed"]] - tick

      if (inherits(fit, "error")) {
        prog_tick(pr, ok = FALSE, secs = secs, msg = conditionMessage(fit),
                  label = label)
        next
      }

      drawn <- predict(fit, newdata = test, type = "survival", times = edges,
                       draws = TRUE)
      shat <- apply(drawn, c(2L, 3L), mean)
      lo <- apply(drawn, c(2L, 3L), stats::quantile, (1 - LEVEL) / 2)
      hi <- apply(drawn, c(2L, 3L), stats::quantile, 1 - (1 - LEVEL) / 2)
      err <- surv_error(shat, strue)
      aux <- fit[["aux"]]

      out[[length(out) + 1L]] <- data.frame(
        truth = truth_name, family = nm, arm = arm, rep = rep,
        s_rmse = err[["rmse"]], rank = err[["rank"]],
        cover = mean(strue >= lo & strue <= hi), width = mean(hi - lo),
        k_0 = mean(aux[, "k_0"]), k_0_fixed = fit[["family_opts"]][["k_0"]],
        clusters = mean(aux[, "clusters"]), secs = secs)

      prog_tick(pr, secs = secs, label = label)
    }
  }

  do.call(rbind, out)
}

N_FITS <- length(TRUTHS) * N_REP * length(ARMS) * length(DPM_FAMILIES)

pr <- prog_init(total = N_FITS, title = "Drawing k_0", unit = "fit",
                kind = "simulation")
on.exit(prog_end(pr, "failed", "aborted before the last replicate"),
        add = TRUE)

main <- list()

checkpoint <- function(complete) {
  saveRDS(list(main = do.call(rbind, main), complete = complete,
               level = LEVEL, when = Sys.time()), OUT)
}

for (nm in names(TRUTHS)) {
  for (r in seq_len(N_REP)) {
    main[[length(main) + 1L]] <- one_rep(nm, r)
    checkpoint(FALSE)
  }
}

checkpoint(TRUE)
on.exit()
prog_end(pr, "done", sprintf("%d replicates", length(main)))
cat("wrote", OUT, "\n")
