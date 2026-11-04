# A skewed log-time error, and the coverage of every survival family's
# intervals for S(t | x).
#
# What is tested. Two claims, either of which could be false.
#   (1) On a skewed log-time error of constant spread, which none of the seven
#       truths of `_dev/survival-sim.R` has, `dpm_aft_ls()` costs little
#       against `dpm_aft()`. `dpm_ls()` pays exactly this cost without
#       censoring (`_dev/dpm-ls-forms.R`, `_dev/dpm-ls-homoskedastic.R`):
#       chance variation in its scale forest moves the estimated mean of a
#       skewed error, because the distance between its mean and its bulk grows
#       with its scale.
#   (2) Every family's 90% intervals for S(t | x) cover near nominal wherever
#       its model is right: `dpm_aft()` and `dpm_aft_ls()` on every
#       accelerated failure time truth they nest, `ph()` on the proportional
#       hazards truths, and each parametric family on its own truth. It could
#       fail by under-coverage (a prior that attenuates, a chain too short) or
#       by over-coverage, as the hard-rule mean intervals did elsewhere.
#
# The measurement. The eight truths (the seven, plus log T = m(x) + 0.7 W
# with W a standardized gamma(1.5) error) by the seven families, five
# replicates, with `_dev/survival-sim.R`'s data and settings, read from its
# definitions. Per fit: `s_rmse`, the RMSE of S(t | x) over the evaluation
# grid, which on the seven existing truths should reproduce that script's
# rerun exactly (same build, same data, same fitting order); `cover`, the
# share of (subject, time) cells whose 90% posterior interval holds the true
# S(t | x); and `width`, the mean width of those intervals. The columns to
# read are `dpm_aft_ls()` minus `dpm_aft()` on the skewed truth, paired by
# replicate, and the coverage table.
#
# What each outcome would mean.
#   - (1) within two SE: no skewed-error cost at this sample size, and the
#     vignettes stand. Behind by more: `vignette("survival")` gets one
#     sentence, as `vignette("families")` has for `dpm_ls()`.
#   - (2) coverage between about .85 and .97 where the model is right: the
#     intervals can be used as stated. Clear under-coverage there: a
#     calibration problem, recorded and investigated before any interval
#     claim is made. Under-coverage where a family is misspecified is expected
#     and is reported as such.
#
# Rerun 2026-11-03 after `DPMFamily` began drawing the mixture's reported
# weights from Dirichlet(n_k, alpha) given the counts, in place of the count
# shares (`_dev/dpm-aft-coverage.R` found the fixed weights understating the
# uncertainty of S(t | x) on bimodal errors). Prediction: the two Dirichlet
# process families' coverage rises on the bimodal truth, modestly at the
# default `k_s`, stays about where it was on the skewed one (approximation, not
# uncertainty) and elsewhere, and their point accuracy hardly moves; the four
# families fitted before `dpm_aft()` in each replicate reproduce exactly, and
# the two after it move by chance, since the draws shift the random number
# stream. If the bimodal coverage does not rise, the engine's draw is not doing
# what the recomputation in R did, and that comes first.
#
# Run with: BARTISAN_LIB=<library holding the working tree> \
#             Rscript _dev/survival-skew-coverage.R
# Writes:   _dev/survival-skew-coverage.rds

lib <- Sys.getenv("BARTISAN_LIB")

if (nzchar(lib)) {
  .libPaths(c(lib, .libPaths()))
}

lines <- readLines("_dev/survival-sim.R")
eval(parse(text = lines[seq_len(grep("^# ---- one replicate", lines) - 1L)]))

OUT <- "_dev/survival-skew-coverage.rds"
N_REP <- as.integer(Sys.getenv("SIM_REPS", N_REP))
LEVEL <- 0.90

# A standardized gamma(1.5) error on log T, mean zero and variance one, so
# that the predictor is still the conditional mean of log T.
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

  for (nm in names(FAMILIES)) {
    tick <- proc.time()[["elapsed"]]
    fit <- tryCatch(bartisan(FORM, data = train, family = FAMILIES[[nm]](),
                             control = ctrl, verbose = FALSE),
                    error = function(e) e)
    secs <- proc.time()[["elapsed"]] - tick

    if (inherits(fit, "error")) {
      prog_tick(pr, ok = FALSE, secs = secs, msg = conditionMessage(fit),
                label = sprintf("%s / %s rep %d", truth_name, nm, rep))
      next
    }

    drawn <- predict(fit, newdata = test, type = "survival", times = edges,
                     draws = TRUE)
    shat <- apply(drawn, c(2L, 3L), mean)
    lo <- apply(drawn, c(2L, 3L), stats::quantile, (1 - LEVEL) / 2)
    hi <- apply(drawn, c(2L, 3L), stats::quantile, 1 - (1 - LEVEL) / 2)
    err <- surv_error(shat, strue)

    out[[nm]] <- data.frame(
      truth = truth_name, family = nm, rep = rep,
      s_rmse = err[["rmse"]], rank = err[["rank"]],
      cover = mean(strue >= lo & strue <= hi), width = mean(hi - lo),
      secs = secs)

    prog_tick(pr, secs = secs,
              label = sprintf("%s / %s rep %d", truth_name, nm, rep))
  }

  do.call(rbind, out)
}

N_FITS <- length(TRUTHS) * N_REP * length(FAMILIES)

pr <- prog_init(total = N_FITS, title = "Survival, skewed error and coverage",
                unit = "fit", kind = "simulation")
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
