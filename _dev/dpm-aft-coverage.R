# Why `dpm_aft()` under-covers S(t | x) on bimodal and skewed log-time errors.
#
# What is tested. `_dev/survival-skew-coverage.R` found 90% intervals for
# S(t | x) covering .55 on the bimodal truth and .67 on the skewed one, with
# the narrowest intervals and the best point accuracy of any family, against
# .88 to .93 on truths whose error is normal in shape. Four explanations, each
# with an arm that could rule it out:
#   H1 the chain is too short for the mixture's labels to mix, so the
#      posterior is too narrow. Arm `long`: four chains of 2000 + 1000.
#   H2 the Dirichlet process posterior settles on one error shape more firmly
#      than the data justify. Arm `oracle`: the true mean as an offset, with a
#      forest on x5 alone, which the mean does not use (a forest needs a
#      predictor), so only the error is estimated in effect; it also scores
#      the coverage of `error_density()` against the true error density.
#   H3 the mean forest's shrinkage biases f(x), and a steep S(t | x) turns
#      small location errors into misses. If `oracle` restores coverage and
#      `long` does not, the cause is in f.
#   H4 the imputation of censored times. Arm `uncensored`: no censoring.
#
# The measurement. Three truths of `_dev/survival-sim.R` and
# `_dev/survival-skew-coverage.R` (bimodal and skewed, and log-normal as the
# control that already covers), five replicates, their data and seeds. Arm
# `default` is the earlier fit on the same data, with a different chain, since
# there `dpm_aft()` was fitted after four other families. Per
# arm: `cover`, the share of (subject, time) cells whose 90% interval holds
# the true S(t | x); `bias`, the mean signed error of the posterior mean;
# `sd`, the mean posterior standard deviation; `z`, the RMSE of the
# standardized errors (one for a calibrated posterior, above one when the
# intervals are too narrow); `s_rmse`; and for `oracle`, `dens_cover`, the
# share of points on a grid over the middle 96% of the true error
# distribution whose 90% interval holds the true density. Read `cover` and `z`
# first.
#
# What each outcome would mean.
#   - `long` restores coverage (about .85 or more): mixing; interval claims for
#     `dpm_aft()` need more draws, and `diagnose()` or the defaults should say
#     so.
#   - `oracle` restores it and `long` does not: bias in f, which BART's
#     shrinkage causes in general; document, since the mixture is not at
#     fault.
#   - neither does, and `dens_cover` is low: the mixture's shape posterior is
#     overconfident; look at the baseline and concentration priors.
#   - only `uncensored` restores it: the censoring imputation.
#
# Run with: BARTISAN_LIB=<library holding the working tree> \
#             Rscript _dev/dpm-aft-coverage.R
# Writes:   _dev/dpm-aft-coverage.rds

lib <- Sys.getenv("BARTISAN_LIB")

if (nzchar(lib)) {
  .libPaths(c(lib, .libPaths()))
}

lines <- readLines("_dev/survival-sim.R")
eval(parse(text = lines[seq_len(grep("^# ---- one replicate", lines) - 1L)]))

OUT <- "_dev/dpm-aft-coverage.rds"
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

# The true density of the error W in log T = m(x) + W, and the range its
# middle 96% covers, for scoring `error_density()`.
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
ARMS <- c("default", "long", "oracle", "uncensored")

ctrl_long <- bartisan_control(num_trees = 50, num_burn = 2000, num_draws = 1000,
                              chains = 4)

score <- function(drawn, strue) {
  shat <- apply(drawn, c(2L, 3L), mean)
  sdev <- apply(drawn, c(2L, 3L), stats::sd)
  lo <- apply(drawn, c(2L, 3L), stats::quantile, (1 - LEVEL) / 2)
  hi <- apply(drawn, c(2L, 3L), stats::quantile, 1 - (1 - LEVEL) / 2)
  z <- (shat - strue) / pmax(sdev, 1e-8)

  c(cover = mean(strue >= lo & strue <= hi), bias = mean(shat - strue),
    sd = mean(sdev), z = sqrt(mean(z^2)),
    s_rmse = sqrt(mean((shat - strue)^2)))
}

one_rep <- function(truth_name, rep) {
  out <- list()

  for (arm in ARMS) {
    set.seed(20260828L + 1000L * match(truth_name, names(TRUTHS)) + rep)

    gen <- TRUTHS[[truth_name]]
    train <- make_x(N_TRAIN)
    test  <- make_x(N_TEST)
    tr <- gen(train)
    te <- gen(test)

    censor <- if (arm == "uncensored") 0 else CENSOR
    ctr <- apply_censoring(tr$time, censor)
    train$time <- ctr$time; train$status <- ctr$status
    cte <- apply_censoring(te$time, censor)
    test$time <- cte$time; test$status <- cte$status

    edges <- unique(quantile(tr$time, seq(0.05, 0.95, length.out = N_GRID)))
    strue <- te$surv(edges)

    tick <- proc.time()[["elapsed"]]
    fit <- tryCatch({
      if (arm == "oracle") {
        # The main forest needs a predictor to split on, so it gets x5, which
        # the true mean does not use: it can fit only noise around the truth.
        bartisan(Surv(time, status) ~ x5, data = train,
                 family = dpm_aft(update_baseline = FALSE),
                 offset = signal(train), control = ctrl, verbose = FALSE)
      }
      else {
        # The fixed baseline was the default when this was run.
        bartisan(FORM, data = train, family = dpm_aft(update_baseline = FALSE),
                 control = if (arm == "long") ctrl_long else ctrl,
                 verbose = FALSE)
      }
    }, error = function(e) e)
    secs <- proc.time()[["elapsed"]] - tick

    if (inherits(fit, "error")) {
      prog_tick(pr, ok = FALSE, secs = secs, msg = conditionMessage(fit),
                label = sprintf("%s / %s rep %d", truth_name, arm, rep))
      next
    }

    drawn <- predict(fit, newdata = test, type = "survival", times = edges,
                     draws = TRUE,
                     offset = if (arm == "oracle") signal(test))
    sc <- score(drawn, strue)

    dens_cover <- NA_real_

    if (arm == "oracle") {
      et <- error_truth[[truth_name]]
      at <- seq(et$range[1L], et$range[2L], length.out = 101L)
      ed <- error_density(fit, at = at, level = LEVEL)
      dens_cover <- mean(et$dens(at) >= ed$lower & et$dens(at) <= ed$upper)
    }

    aux <- fit[["aux"]]
    out[[arm]] <- data.frame(
      truth = truth_name, arm = arm, rep = rep, t(sc),
      dens_cover = dens_cover, alpha = mean(aux[, "alpha"]),
      clusters = mean(aux[, "clusters"]), secs = secs)

    prog_tick(pr, secs = secs,
              label = sprintf("%s / %s rep %d", truth_name, arm, rep))
  }

  do.call(rbind, out)
}

N_FITS <- length(STUDY) * N_REP * length(ARMS)

pr <- prog_init(total = N_FITS, title = "dpm_aft() coverage", unit = "fit",
                kind = "simulation")
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
