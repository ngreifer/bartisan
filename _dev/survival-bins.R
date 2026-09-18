# The bin-count sensitivity check behind vignettes/survival.Rmd.
#
# Sweeps ph(num_bins = ) over a sixty-fold range on the two proportional-hazards
# truths from _dev/survival-sim.R, recording the error in S(t | x), the error in
# the log hazard ratio r(x) -- which is exactly signal() for these two truths --
# and the loo diagnostics. Writes _dev/survival-bins.rds.

suppressPackageStartupMessages({
  library(bartisan)
  library(survival)
  library(loo)
})

src <- readLines("_dev/survival-sim.R")
cut <- grep("^# ---- run", src)
eval(parse(text = paste(src[1:(cut - 1)], collapse = "\n")), envir = globalenv())

BINS <- c(4, 9, 20, 50, 100, 250)
TRUTH_PH <- c("hazard turns over", "Weibull PH")
N_REP_B <- 3L

# The head of survival-sim.R, evaluated above, brings progress.R with it along
# with its own `OUT`.
OUT <- "_dev/survival-bins.rds"
N_FITS <- length(TRUTH_PH) * N_REP_B * length(BINS)

out <- list()

pr <- prog_init(total = N_FITS, title = "Survival: how many bins ph() needs",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed", "aborted before the last replicate"),
        add = TRUE)

# `survival-results.R` reads this file as a data frame, so how far the run got
# is carried in attributes rather than by wrapping it in a list. Written after
# every replicate, so a killed run still leaves its finished ones readable.
checkpoint <- function(complete) {
  res <- do.call(rbind, out)
  attr(res, "complete") <- complete
  attr(res, "done") <- length(out)
  attr(res, "total") <- N_FITS
  saveRDS(res, OUT)
}

for (tn in TRUTH_PH) {
  for (r in seq_len(N_REP_B)) {
    set.seed(770011L + 100L * match(tn, TRUTH_PH) + r)
    gen <- TRUTHS[[tn]]
    train <- make_x(N_TRAIN); test <- make_x(N_TEST)
    tr <- gen(train); te <- gen(test)
    ctr <- apply_censoring(tr$time, CENSOR)
    train$time <- ctr$time; train$status <- ctr$status
    cte <- apply_censoring(te$time, CENSOR)
    test$time <- cte$time; test$status <- cte$status

    edges <- unique(quantile(tr$time, seq(0.05, 0.95, length.out = N_GRID)))
    strue <- te$surv(edges)
    # The true log hazard ratio, centered the way ph() reports its predictor.
    r_true <- signal(test) - mean(signal(test))

    for (b in BINS) {
      t0 <- Sys.time()
      fit <- bartisan(FORM, data = train, family = ph(num_bins = b),
                      control = ctrl, verbose = FALSE)
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      shat <- predict(fit, newdata = test, type = "survival", times = edges)
      rhat <- drop(predict(fit, newdata = test, type = "link"))
      rhat <- rhat - mean(rhat)
      lo <- suppressWarnings(loo(fit))
      out[[length(out) + 1L]] <- data.frame(
        truth = tn, bins = b, rep = r,
        s_rmse = sqrt(mean((shat - strue)^2)),
        r_rmse = sqrt(mean((rhat - r_true)^2)),
        p_loo = lo$estimates["p_loo", "Estimate"],
        bad_k = mean(pareto_k_values(lo) > 0.7))
      cat(sprintf("[bins] %-18s rep %d bins %3d done\n", tn, r, b))
      flush(stdout())
      prog_tick(pr, i = length(out), secs = secs,
                label = sprintf("%s rep %d, %d bins", tn, r, b))
    }

    checkpoint(FALSE)
  }
}

checkpoint(TRUE)
on.exit()
prog_end(pr, "done", sprintf("%d fits", length(out)))
cat("wrote", OUT, "\n")
