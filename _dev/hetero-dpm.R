# Does a treatment-arm-specific residual variance earn a place in bcf(), and
# does it make sense with a DPM error?
#
# In this package a per-arm variance needs no new feature: `gaussian_ls()` with
# `log_sd = ~ z` puts the treatment in the scale forest's formula and yields one
# variance per arm. What is NOT available is a DPM error whose scale varies,
# because dpm() has a single additive predictor and its scale is one drawn
# scalar. So the question is which half-right model costs less when the truth is
# both heteroscedastic AND non-normal:
#
#   dpm()                    right error shape, one shared scale
#   gaussian_ls(log_sd ~ z)  right scale per arm, normal shape
#   dpm() fitted per arm     both right -- the stand-in for the family that does
#                            not exist, bounding what adding one could buy
#
# A plain gaussian() arm was dropped from an earlier version of this script: it
# answered nothing the other three did not, and it cost a third of the runtime.
#
# The prediction worth testing is that the error model shows up in interval
# coverage rather than in the point estimate, because both mistakes misprice
# uncertainty more than they bias the mean. Read the coverage column first.
#
# Results are written after every (cell, replicate) rather than once at the end.
# An earlier run of this script was killed at 14 of 320 fits and every accuracy
# number was lost, because the only saveRDS() was after the loop. The partial
# file is the point: an interrupted simulation should still be worth reading.
A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))
Sys.setenv(PROGRESS_ROOT = "/Users/NoahGreifer/Dropbox/Research/R/bartisan/.progress-runs")

suppressMessages(library(bartisan))

REPS <- 10L
N <- 600L
OUT <- file.path("_dev", "hetero-dpm.rds")

draw_err <- function(n, shape) {
  if (identical(shape, "t3")) stats::rt(n, df = 3) / sqrt(3) else stats::rnorm(n)
}

make <- function(n, seed, shape, hetero) {
  set.seed(seed)
  d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
  d$z <- rbinom(n, 1L, stats::plogis(-0.4 + 1.6 * d$x1 - 0.8 * d$x2))

  prog <- 2 * sin(pi * d$x1) + 1.5 * d$x2
  tau <- 1 + 1.5 * d$x1
  sd_z <- if (hetero) exp(-0.7 + 1.1 * d$z) else exp(-0.15)

  d$y <- prog + d$z * tau + sd_z * draw_err(n, shape)
  list(data = d, tau = tau, ate = mean(tau))
}

# The ATE by g-computation over the posterior: predict every unit under both
# arms, difference, average, summarize the draws.
ate_draws <- function(fit, d) {
  p1 <- stats::predict(fit, newdata = transform(d, z = 1L), type = "response",
                       draws = TRUE)
  p0 <- stats::predict(fit, newdata = transform(d, z = 0L), type = "response",
                       draws = TRUE)
  list(ate = rowMeans(p1 - p0), cate = colMeans(p1 - p0))
}

# One DPM per arm, then g-computation using each arm's own model for its own
# potential outcome.
ate_draws_split <- function(d, ctrl) {
  p <- lapply(0:1, function(a) {
    f <- bartisan(y ~ x1 + x2 + x3, d[d$z == a, ], family = dpm(),
                  control = ctrl)
    stats::predict(f, newdata = d, type = "response", draws = TRUE)
  })
  k <- min(nrow(p[[1L]]), nrow(p[[2L]]))
  diff <- p[[2L]][seq_len(k), , drop = FALSE] - p[[1L]][seq_len(k), , drop = FALSE]
  list(ate = rowMeans(diff), cate = colMeans(diff))
}

score <- function(dr, dat, label, shape, hetero, rep, secs) {
  q <- stats::quantile(dr$ate, c(0.025, 0.975))
  data.frame(arm = label, shape = shape, hetero = hetero, rep = rep,
             secs = secs, ate_est = mean(dr$ate), ate_truth = dat$ate,
             ate_bias = mean(dr$ate) - dat$ate,
             ate_covered = q[[1L]] <= dat$ate && dat$ate <= q[[2L]],
             ate_width = q[[2L]] - q[[1L]],
             cate_rmse = sqrt(mean((dr$cate - dat$tau)^2)),
             stringsAsFactors = FALSE)
}

quiet <- function(expr) {
  invisible(utils::capture.output(out <- suppressMessages(suppressWarnings(expr))))
  out
}

ctrl <- bartisan_control(num_trees = 50L, num_burn = 500L, num_draws = 1000L,
                         chains = 2L)
# The scale forest is given far fewer trees than the mean: a spread needs less
# resolution than a mean surface, and at 50 trees apiece this arm took three
# times as long as any other.
ctrl_ls <- bartisan_control(num_trees = c(mean = 50L, log_sd = 10L),
                            num_burn = 500L, num_draws = 1000L, chains = 2L)

cells <- expand.grid(shape = "t3", hetero = c(TRUE, FALSE),
                     stringsAsFactors = FALSE)
arms <- c("dpm", "gaussian_ls(z)", "dpm per arm")

pr <- prog_init(total = nrow(cells) * REPS * length(arms),
                title = "Heteroscedastic non-normal errors under BCF",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
k <- 0L

for (i in seq_len(nrow(cells))) {
  shape <- cells$shape[i]
  hetero <- cells$hetero[i]

  for (r in seq_len(REPS)) {
    dat <- make(N, seed = 9000L + 137L * r + 11L * i, shape, hetero)
    d <- dat$data

    for (arm in arms) {
      t0 <- Sys.time()
      dr <- quiet(switch(
        arm,
        "dpm" = ate_draws(
          bcf(y ~ x1 + x2 + x3, treatment = ~ z, data = d, family = dpm(),
              control = ctrl), d),
        "gaussian_ls(z)" = ate_draws(
          bartisan(list(mean = y ~ x1 + x2 + x3 + vc(z), log_sd = ~ z),
                   data = d, family = gaussian_ls(), control = ctrl_ls), d),
        "dpm per arm" = ate_draws_split(d, ctrl)))

      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      k <- k + 1L
      rows[[k]] <- score(dr, dat, arm, shape, hetero, r, secs)
      prog_tick(pr, i = k, secs = secs,
                label = sprintf("%s / hetero=%s rep %d", arm, hetero, r))
    }

    # Written every replicate, not once at the end, so a killed run still
    # leaves every completed replicate on disk.
    saveRDS(list(res = do.call(rbind, rows), reps = REPS, n = N,
                 complete = FALSE, done = k,
                 total = nrow(cells) * REPS * length(arms)), OUT)
  }
}

res <- do.call(rbind, rows)
saveRDS(list(res = res, reps = REPS, n = N, complete = TRUE, done = k,
             total = k), OUT)
on.exit()
prog_end(pr, "done")

agg <- aggregate(cbind(ate_bias, ate_width, cate_rmse, secs) ~ arm + hetero,
                 data = res, FUN = mean)
cov <- aggregate(ate_covered ~ arm + hetero, data = res, FUN = mean)
out <- merge(agg, cov)
print(out[order(out$hetero, out$arm),
          c("hetero", "arm", "ate_bias", "ate_covered", "ate_width",
            "cate_rmse", "secs")], row.names = FALSE, digits = 3)
cat("\nwrote", OUT, "\n")
