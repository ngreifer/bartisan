# The shape of the burn-in transient, and whether the recorded statistic is
# scale-free.
#
# `_dev/transient-scaling.R` reports the first sweep whose log likelihood is
# within two standard deviations of its eventual level, which is the statistic
# `_dev/TASKS.md` records as 36 to 61 sweeps. It came back at 152 to 370. Before
# reading that as a real difference, it has to be separated from the statistic
# itself: the "eventual level" is the mean over the last half of the draws, so a
# longer run has a higher level and a smaller spread and the threshold is harder
# to reach. The statistic is therefore not scale-free, and the recorded numbers
# may simply come from a shorter run.
#
# This saves the whole log likelihood trace for a few cells and reports the same
# statistic recomputed on the first 200, 400 and 800 sweeps, so the run-length
# sensitivity is visible rather than inferred. It also reports how far the trace
# has risen by a given sweep as a fraction of its total rise, which does not
# depend on where the run was stopped.

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

DRAWS <- 800L
REPS <- 2L
OUT <- "_dev/transient-trace.rds"

CELLS <- expand.grid(rep = seq_len(REPS),
                     rule = c("soft", "hard"),
                     n = c(1500L, 10000L, 50000L),
                     stringsAsFactors = FALSE)

GATES <- c(soft = "smoothstep", hard = "hard")

sim <- function(n, rep) {
  set.seed(7000 + 97 * rep + 10L)
  p <- 10L
  X <- matrix(stats::runif(n * p), n, p,
              dimnames = list(NULL, paste0("x", seq_len(p))))
  truth <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 +
    10 * X[, 4] + 5 * X[, 5]

  cbind(data.frame(y = truth + stats::rnorm(n)), as.data.frame(X))
}

reach <- function(x) {
  half <- seq(floor(length(x) / 2) + 1L, length(x))
  level <- mean(x[half])
  spread <- stats::sd(x[half])
  hit <- which(x >= level - 2 * spread)

  if (!length(hit)) return(NA_integer_)

  as.integer(hit[1L])
}

# Where the trace stands at `at`, as a fraction of the whole rise from its first
# value to its final plateau. Unlike `reach()` this does not move when the run
# is truncated, provided the plateau itself is reached.
risen <- function(x, at) {
  plateau <- mean(x[seq(floor(length(x) / 2) + 1L, length(x))])
  (x[at] - x[1L]) / (plateau - x[1L])
}

pr <- prog_init(total = nrow(CELLS), title = "Transient trace shape",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
traces <- list()

for (i in seq_len(nrow(CELLS))) {
  cell <- CELLS[i, ]
  label <- sprintf("n=%d %s rep %d", cell$n, cell$rule, cell$rep)

  out <- prog_do(pr, i, function(.) {
    dat <- sim(cell$n, cell$rep)
    fit <- bartisan(y ~ ., data = dat, family = gaussian(),
                    num_burn = 0L, num_draws = DRAWS, num_thin = 1L,
                    chains = 1L, gate = GATES[[cell$rule]])
    as.numeric(fit[["loglik"]])
  }, label = label)

  if (is.null(out)) next

  traces[[label]] <- out

  rows[[length(rows) + 1L]] <- data.frame(
    n = cell$n, rule = cell$rule, rep = cell$rep,
    reach_200 = reach(out[1:200]),
    reach_400 = reach(out[1:400]),
    reach_800 = reach(out),
    risen_36 = round(risen(out, 36L), 3),
    risen_61 = round(risen(out, 61L), 3),
    risen_200 = round(risen(out, 200L), 3),
    risen_400 = round(risen(out, 400L), 3))

  saveRDS(list(res = do.call(rbind, rows), traces = traces,
               complete = FALSE, done = length(rows),
               total = nrow(CELLS)), OUT)
}

res <- do.call(rbind, rows)
saveRDS(list(res = res, traces = traces, complete = TRUE,
             done = nrow(res), total = nrow(CELLS)), OUT)

on.exit()
prog_end(pr, "done")

print(res, row.names = FALSE)
