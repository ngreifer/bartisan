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

out <- list()
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
}

saveRDS(do.call(rbind, out), "_dev/survival-timing.rds")
cat("wrote _dev/survival-timing.rds\n")
