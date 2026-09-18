# Reproducing the un-augmented ordinal probit stall.
#
# `_dev/augment-benchmark.R` measured one `ordinal("probit")` fit at 3882
# seconds where the other fourteen replicates of the same cell took 30 to 43.
# This script replays that one cell. The replay is exact rather than
# approximate, which is worth stating because it is not obvious: a fit restores
# the session's RNG state when it returns, so the state entering the second arm
# of a replicate is the state `make()` left, and nothing earlier in the
# benchmark reaches it. Replaying one replicate therefore needs only its own
# seed, its own first arm, and the same control.
#
# `STALL_REPS` selects replicates (default the slow one), `STALL_GATE` the gate.
A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))
Sys.setenv(PROGRESS_ROOT = "/Users/NoahGreifer/Dropbox/Research/R/bartisan/.progress-runs")

suppressMessages(library(bartisan))

# So a profiler can be attached to the fit while it is still running.
writeLines(as.character(Sys.getpid()), "_dev/ordinal-stall.pid")

N <- 400L
P <- 8L

# Verbatim from `_dev/augment-benchmark.R`. It has to be verbatim: the draws it
# takes are what leave the stream where the fit picks it up, so a shortened
# version that generates only `y_ord` would fit a different chain.
make <- function(seed) {
  set.seed(seed)
  X <- matrix(runif(N * P), N, P)
  colnames(X) <- paste0("x", seq_len(P))
  d <- as.data.frame(X)
  mu <- 2 * sin(pi * d$x1) + 1.5 * d$x2
  lin <- mu - mean(mu)

  d$y_bin <- rbinom(N, 1L, stats::plogis(lin))
  d$y_ord <- factor(cut(mu + rnorm(N, 0, 0.5), 3L), ordered = TRUE)
  d$y_cnt <- rpois(N, exp(0.4 + 0.5 * mu))
  d$y_od <- rnbinom(N, size = 2, mu = exp(0.4 + 0.5 * mu))
  d$y_zi <- ifelse(runif(N) < 0.3, 0L, rpois(N, exp(0.4 + 0.5 * mu)))
  d$y_zo <- ifelse(runif(N) < 0.3, 0L, rnbinom(N, size = 2, mu = exp(0.4 + 0.5 * mu)))
  d$y_mn <- factor(sample(3L, N, TRUE, prob = c(0.4, 0.35, 0.25)))
  tm <- rexp(N, exp(-0.5 - mu / 2))
  ev <- rbinom(N, 1L, 0.7)
  d$y_surv <- survival::Surv(tm, ev)
  d
}

rhs <- paste(paste0("x", seq_len(P)), collapse = " + ")
form <- stats::as.formula(paste("y_ord ~", rhs))

quiet <- function(expr) {
  invisible(utils::capture.output(out <- suppressMessages(suppressWarnings(expr))))
  out
}

ess_pair <- function(fit) {
  tab <- diagnose(fit)[["table"]]
  tab <- tab[!grepl("average over observations", tab[["quantity"]], fixed = TRUE), ]
  c(worst = min(tab[["ess_bulk"]], na.rm = TRUE),
    median = stats::median(tab[["ess_bulk"]], na.rm = TRUE))
}

reps <- as.integer(strsplit(Sys.getenv("STALL_REPS", "13"), ",")[[1L]])
gate <- Sys.getenv("STALL_GATE", "hard")
OUT <- sprintf("_dev/ordinal-stall-%s.rds", gate)

fam <- ordinal("probit")

pr <- prog_init(total = length(reps) * 2L,
                title = sprintf("Ordinal probit stall, %s rules", gate),
                unit = "fit", kind = "benchmark")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
k <- 0L

for (r in reps) {
  d <- make(600L + r)

  for (aug in c(TRUE, FALSE)) {
    smoke <- nzchar(Sys.getenv("STALL_SMOKE"))
    ctrl <- bartisan_control(num_trees = if (smoke) 5L else 50L,
                             num_burn = if (smoke) 20L else 500L,
                             num_draws = if (smoke) 20L else 1000L,
                             chains = 2L, gate = gate, augment = aug)
    t0 <- Sys.time()
    fit <- quiet(bartisan(form, d, family = fam, control = ctrl))
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ess <- quiet(ess_pair(fit))

    k <- k + 1L
    rows[[k]] <- data.frame(rep = r, gate = gate, augment = aug, secs = secs,
                            ess_worst = ess[["worst"]],
                            ess_median = ess[["median"]],
                            stringsAsFactors = FALSE)
    prog_tick(pr, i = k, secs = secs,
              label = sprintf("rep %d augment=%s", r, aug))
    saveRDS(list(res = do.call(rbind, rows), complete = FALSE, done = k,
                 total = length(reps) * 2L), OUT)
    cat(sprintf("rep %2d augment=%-5s %10.2f s\n", r, aug, secs))
  }
}

saveRDS(list(res = do.call(rbind, rows), complete = TRUE, done = k,
             total = k), OUT)
on.exit()
prog_end(pr, "done")
