# Clean timing run for the survival vignette's speed figure. Kept separate from
# _dev/survival-sim.R because the `secs` column recorded there is measured while
# other work may be running, and the speed comparison is only meaningful when
# nothing is competing for the machine.

suppressPackageStartupMessages({
  library(bartisan)
  library(survival)
})

src <- readLines("_dev/survival-sim.R")
cut <- grep("^# ---- run", src)
eval(parse(text = paste(src[1:(cut - 1)], collapse = "\n")), envir = globalenv())

N_TIME_REP <- 3L

# The head of survival-sim.R, evaluated above, brings progress.R with it along
# with its own `OUT`.
OUT <- "_dev/survival-timing.rds"
N_FITS <- N_TIME_REP * (length(FAMILIES) + 1L)

out <- list()

pr <- prog_init(total = N_FITS, title = "Survival families: clean timing",
                unit = "fit", kind = "benchmark")
on.exit(prog_end(pr, "failed", "aborted before the last replicate"),
        add = TRUE)

# `survival-results.R` reads this file as a data frame, so how far the run got
# is carried in attributes rather than by wrapping it in a list. Written after
# every replicate, so a killed run still leaves its finished ones readable.
checkpoint <- function(complete, done) {
  res <- do.call(rbind, out)
  attr(res, "complete") <- complete
  attr(res, "done") <- done
  attr(res, "total") <- N_FITS
  saveRDS(res, OUT)
}

for (r in seq_len(N_TIME_REP)) {
  set.seed(4242L + r)
  train <- make_x(N_TRAIN)
  tr <- TRUTHS[["hazard turns over"]](train)
  ctr <- apply_censoring(tr$time, CENSOR)
  train$time <- ctr$time; train$status <- ctr$status

  for (nm in names(FAMILIES)) {
    tick <- proc.time()[["elapsed"]]
    bartisan(FORM, data = train, family = FAMILIES[[nm]](), control = ctrl,
             verbose = FALSE)
    out[[length(out) + 1L]] <- data.frame(
      family = nm, rep = r, secs = proc.time()[["elapsed"]] - tick)
    cat(sprintf("[time] rep %d %-18s %.2fs\n", r, nm,
                out[[length(out)]]$secs)); flush(stdout())
    prog_tick(pr, i = length(out), secs = out[[length(out)]]$secs,
              label = sprintf("%s rep %d", nm, r))
  }

  # The discrete-time route, whose cost is the expansion rather than the family.
  edges <- unique(quantile(tr$time, seq(0.05, 0.95, length.out = N_GRID)))
  tick <- proc.time()[["elapsed"]]
  long <- expand_dt(train, edges, XNAMES)
  bartisan(ev ~ ., data = long, family = binomial("probit"), control = ctrl,
           verbose = FALSE)
  out[[length(out) + 1L]] <- data.frame(
    family = "discrete-time probit", rep = r,
    secs = proc.time()[["elapsed"]] - tick)
  cat(sprintf("[time] rep %d %-18s %.2fs (%d rows)\n", r,
              "discrete-time probit", out[[length(out)]]$secs, nrow(long)))
  flush(stdout())
  prog_tick(pr, i = length(out), secs = out[[length(out)]]$secs,
            label = sprintf("discrete time rep %d", r))

  checkpoint(FALSE, length(out))
}

checkpoint(TRUE, length(out))
on.exit()
prog_end(pr, "done", sprintf("%d fits", length(out)))
cat("wrote", OUT, "\n")
