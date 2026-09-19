# How long the burn-in transient is, as a function of n.
#
# `_dev/TASKS.md` records the transient as 36 to 61 sweeps, which is what struck
# both warm-start items: a warm start can remove at most that much of a
# 1000-sweep run. Those fits are at package scale, where `rhc` ships at n = 1500
# and `_dev/benchmark.Rmd` runs n = 5000. He and Hahn's warm start moved coverage
# from 0.74 to 0.96 at n = 10000, p = 30, and Krantsevich, He and Hahn's at
# n = 5000, p = 50, so the measurement and the claim are in different regimes.
#
# The transient should not be constant in n. `P_BIRTH_DEATH = 0.7` gives 0.35
# birth attempts per tree per sweep and `src/node.h` puts birth acceptance near
# a third, so a tree gains about 0.12 internal nodes per sweep before deaths.
# The transient is then proportional to how many internal nodes the data
# supports, which grows with n. This measures whether it does.
#
# Two statistics per fit, both read off draws retained from sweep 1 with
# `num_burn = 0`:
#
#   loglik    the first sweep whose log likelihood is within two standard
#             deviations of its eventual level, matching the recorded statistic.
#   share     the first sweep at which the share of splits falling on the five
#             predictors that matter is within two standard deviations of its
#             own plateau. This is the recorded sparse-case statistic, and it is
#             the one that took longest before.
#
# Both plateaus are the mean and standard deviation over the last half of the
# draws.

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

SMOKE <- nzchar(Sys.getenv("TRANSIENT_SMOKE"))

OUT <- "_dev/transient-scaling.rds"

if (SMOKE) {
  NS <- c(1500L)
  PS <- c(10L)
  REPS <- 1L
  DRAWS <- 60L
  OUT <- "_dev/transient-scaling-smoke.rds"
} else {
  NS <- c(1500L, 10000L, 50000L)
  PS <- c(10L, 30L)
  REPS <- 3L
  DRAWS <- 800L
}

GATES <- c(soft = "smoothstep", hard = "hard")
RELEVANT <- paste0("x", 1:5)

# The Friedman function, with `p - 5` irrelevant predictors after the five it
# uses. The record's sparse arm is 30 predictors with 25 irrelevant, which is
# p = 30 here.
sim <- function(n, p, rep) {
  set.seed(7000 + 97 * rep + p)
  X <- matrix(stats::runif(n * p), n, p,
              dimnames = list(NULL, paste0("x", seq_len(p))))
  truth <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 +
    10 * X[, 4] + 5 * X[, 5]

  cbind(data.frame(y = truth + stats::rnorm(n)), as.data.frame(X))
}

# The first index at which a rising series reaches within two standard
# deviations of the level it settles at. NA when it never does, which would mean
# the run is too short to have a plateau to compare against.
reach <- function(x) {
  half <- seq(floor(length(x) / 2) + 1L, length(x))
  level <- mean(x[half])
  spread <- stats::sd(x[half])
  hit <- which(x >= level - 2 * spread)

  if (!length(hit)) return(NA_integer_)

  as.integer(hit[1L])
}

run <- function(n, p, rule, rep) {
  dat <- sim(n, p, rep)

  elapsed <- system.time(
    fit <- bartisan(y ~ ., data = dat, family = gaussian(),
                    num_burn = 0L, num_draws = DRAWS, num_thin = 1L,
                    chains = 1L, gate = GATES[[rule]])
  )[["elapsed"]]

  counts <- fit[["counts"]][[1L]]
  share <- rowSums(counts[, RELEVANT, drop = FALSE]) / pmax(rowSums(counts), 1)

  data.frame(n = n, p = p, rule = rule, rep = rep,
             draws = DRAWS, seconds = elapsed,
             loglik_sweep = reach(fit[["loglik"]]),
             share_sweep = reach(share),
             loglik_plateau = mean(utils::tail(fit[["loglik"]], DRAWS %/% 2)),
             share_plateau = mean(utils::tail(share, DRAWS %/% 2)),
             share_first = share[1L],
             leaves_plateau = mean(utils::tail(
               rowSums(counts) / fit[["num_trees"]][1L], DRAWS %/% 2)))
}

cells <- expand.grid(rep = seq_len(REPS), rule = names(GATES), p = PS, n = NS,
                     stringsAsFactors = FALSE)
cells <- cells[order(cells$n, cells$p, cells$rule, cells$rep), ]

pr <- prog_init(total = nrow(cells), title = "Burn-in transient against n",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()

for (i in seq_len(nrow(cells))) {
  cell <- cells[i, ]
  label <- sprintf("n=%d p=%d %s rep %d", cell$n, cell$p, cell$rule, cell$rep)

  rows[[length(rows) + 1L]] <- prog_do(pr, i, function(.) {
    run(cell$n, cell$p, cell$rule, cell$rep)
  }, label = label)

  # After every fit, so a killed run leaves every finished one behind.
  saveRDS(list(res = do.call(rbind, rows), complete = FALSE,
               done = length(rows), total = nrow(cells)), OUT)
}

res <- do.call(rbind, rows)
saveRDS(list(res = res, complete = TRUE, done = nrow(res),
             total = nrow(cells)), OUT)

on.exit()
prog_end(pr, "done")

print(res)
