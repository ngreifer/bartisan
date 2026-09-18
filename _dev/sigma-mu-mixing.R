# `sigma_mu` was the worst-mixing quantity in all three fits in
# `_dev/bcf-mixing.R`, at rhat 1.36 to 1.54 and ess_bulk under 10 out of 3000
# draws, including for the compiled `gaussian()`. That is not a BCF question. Is
# it this design, or is it general?

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

OUT <- "_dev/sigma-mu-mixing.rds"

pr <- prog_init(total = 5L, title = "Whether sigma_mu's mixing is general",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed", "aborted before the last design"), add = TRUE)

rows <- list()

# Written after every design rather than at the end, so a run that is killed
# leaves its finished designs readable. `complete` is what tells a reader which
# of the two it is looking at.
checkpoint <- function(complete) {
  saveRDS(list(res = do.call(rbind, rows), complete = complete,
               done = length(rows), total = 5L), OUT)
}

# `fit` is a promise, so it is the first line here that runs the sampler, which
# is what makes the elapsed time measurable at all.
worst <- function(label, fit) {
  t0 <- Sys.time()
  r <- fit$rhat
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  prog_tick(pr, secs = secs, label = label)

  sm <- r[startsWith(r$quantity, "sigma_mu"), ]
  et <- r[startsWith(r$quantity, "eta"), ]
  cat(sprintf("%-34s sigma_mu rhat %5.3f ess %6.1f | eta rhat %5.3f ess %6.1f\n",
              label, max(sm$rhat), min(sm$ess_bulk), max(et$rhat), min(et$ess_bulk)))

  rows[[label]] <<- data.frame(
    design = label, secs = secs, sigma_mu_rhat = max(sm$rhat),
    sigma_mu_ess = min(sm$ess_bulk), eta_rhat = max(et$rhat),
    eta_ess = min(et$ess_bulk))
  checkpoint(FALSE)
  invisible(NULL)
}

ctrl <- bartisan_control(num_burn = 750, num_draws = 750)
n <- 1000

set.seed(1)
# Friedman, the package's own benchmark design.
x <- matrix(runif(n * 10), n, 10, dimnames = list(NULL, paste0("x", 1:10)))
fr <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
  10 * x[, 4] + 5 * x[, 5]
d1 <- data.frame(y = fr + rnorm(n), x)
set.seed(1); worst("friedman, signal/noise high", bartisan(y ~ ., data = d1, family = gaussian(), control = ctrl, chains = 4))

# Same design, noise dominating.
d2 <- data.frame(y = fr / 10 + rnorm(n), x)
set.seed(1); worst("friedman, signal/noise low", bartisan(y ~ ., data = d2, family = gaussian(), control = ctrl, chains = 4))

# Pure noise: nothing for the forest to find.
d3 <- data.frame(y = rnorm(n), x)
set.seed(1); worst("pure noise", bartisan(y ~ ., data = d3, family = gaussian(), control = ctrl, chains = 4))

# A single strong linear predictor.
d4 <- data.frame(y = 3 * x[, 1] + rnorm(n), x)
set.seed(1); worst("one linear predictor", bartisan(y ~ ., data = d4, family = gaussian(), control = ctrl, chains = 4))

# The fixed-scale alternative the warning already recommends.
set.seed(1); worst("friedman, update_sigma_mu = FALSE",
                   bartisan(y ~ ., data = d1, family = gaussian(), chains = 4,
                            control = bartisan_control(num_burn = 750, num_draws = 750,
                                                       update_sigma_mu = FALSE)))

checkpoint(TRUE)
on.exit()
prog_end(pr, "done", sprintf("%d designs", length(rows)))
cat("\nwrote", OUT, "\n")
