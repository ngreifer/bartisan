# What sharing a topology costs in time, measured on a quiet machine.
#
# Separate from `_dev/shared-topology-sim.R` for the reason `_dev/survival-
# timing.R` is separate from the rest of the survival scripts: a timing taken
# while eight workers are fitting other models is not a timing. This runs one
# fit at a time, in one process, and nothing else should be running.
#
# The prediction going in was that sharing should be roughly H times faster,
# since one topology is grown where H were. It is not, and the accounting is
# worth stating because it says where the saving actually is. Per sweep, a
# forest of T trees costs
#
#   T tree moves        one node's support each, plus a Laplace fit per
#                       component, plus the rule and the gates
#   T leaf refreshes    every leaf of every tree, so a full pass over the data
#                       per tree per component
#   T/every bandwidth   a full rebuild of the tree's supports and a full
#                       likelihood pass, per tree per component
#
# Sharing removes H - 1 of the H rule draws and H - 1 of the H gate
# evaluations, and it collapses the bandwidth move from H rebuilds to one. It
# does *not* touch the leaf refresh: the model has exactly as many leaf values
# as it had before, and each still has to be drawn against its own component.
# So the saving is bounded by the share of the run that is not leaf refresh,
# and the grid below is chosen to show that boundary rather than a single
# number: `num_trees` and `n` move the balance between the two, and
# `bandwidth_every` and `gate` move the cost of the part that does shrink.

library(bartisan)

# Sharing the topology is not an argument: it is engaged by the package's
# internal flag, for the reasons recorded beside it in `R/control.R`. These
# scripts are the record of how it was measured, so they reach it the same way
# the tests do.
share_forests <- function(on) {
  flags <- bartisan:::the
  flags$share_forests <- isTRUE(on)
}

# Live progress, per the `/live-progress` skill. One tick per fit: the grid's
# rows differ by more than an order of magnitude in cost, so a per-row estimate
# would be wrong for most of the run.
source(path.expand(Sys.getenv(
  "LIVE_PROGRESS", "~/.claude/skills/live-progress/assets/progress.R")))

REPS <- 3L
OUT <- file.path("_dev", "shared-topology-timing.rds")

# Fits: every grid row runs two arms of REPS fits, and so does every width of
# the varying-coefficient series below.
N_FITS <- (11L + 3L) * 2L * REPS

sim_ls <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  d <- as.data.frame(X)
  d$y <- 3 * sin(pi * d$x1) + 2 * d$x2 +
    rnorm(n, 0, exp(-0.9 + 0.9 * d$x1 + 0.7 * d$x2))
  d
}

sim_zi <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  d <- as.data.frame(X)
  lam <- exp(0.4 + 1.2 * d$x1)
  pz <- plogis(-0.6 + 1.6 * (0.9 * d$x1 + 0.7 * d$x2))
  d$y <- ifelse(runif(n) < pz, 0L, rpois(n, lam))
  d
}

sim_vc <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  d <- as.data.frame(X)
  d$z <- rnorm(n)
  d$y <- 3 * sin(pi * d$x1) + d$z * (1 + 1.5 * d$x2) + rnorm(n)
  d
}

# One row of the grid, timed both ways with everything else held equal.
time_cell <- function(pr, model, n, p, num_trees, num_burn, num_draws,
                      bandwidth_every, gate, seed) {
  dat <- switch(model, gaussian_ls = sim_ls(n, p, seed),
                zi_poisson = sim_zi(n, p, seed),
                vc = sim_vc(n, p, seed))

  rhs <- paste(paste0("x", seq_len(p)), collapse = " + ")
  form <- stats::as.formula(if (identical(model, "vc")) {
    paste("y ~", rhs, "+ vc(z)")
  } else {
    paste("y ~", rhs)
  })

  fam <- switch(model, gaussian_ls = gaussian_ls(), zi_poisson = zi_poisson(),
                vc = gaussian())

  out <- numeric(0)

  for (shared in c(FALSE, TRUE)) {
    # Before the control is built: the flag is read there, not at fit time.
    share_forests(shared)
    ctrl <- bartisan_control(num_trees = num_trees, num_burn = num_burn,
                             num_draws = num_draws, gate = gate,
                             bandwidth_every = bandwidth_every)
    reps <- vapply(seq_len(REPS), function(r) {
      secs <- system.time(bartisan(form, data = dat, family = fam,
                                   control = ctrl))[["elapsed"]]
      prog_tick(pr, secs = secs,
                label = sprintf("%s n=%d p=%d trees=%d %s %s rep %d", model, n,
                                p, num_trees, gate,
                                if (shared) "shared" else "separate", r))
      secs
    }, numeric(1L))
    out <- c(out, min(reps))
  }

  data.frame(model = model, n = n, p = p, num_trees = num_trees,
             num_draws = num_draws, bandwidth_every = bandwidth_every,
             gate = gate, separate = out[1], shared = out[2],
             speedup = out[1] / out[2], stringsAsFactors = FALSE)
}

# The grid. `n` and `num_trees` move the balance between the tree moves, which
# sharing shrinks, and the leaf refresh, which it does not; `bandwidth_every`
# and `gate` move the cost of the part that shrinks.
grid <- rbind(
  data.frame(model = "gaussian_ls", n = 400,  p = 20, num_trees = 50, bandwidth_every = 1, gate = "smoothstep"),
  data.frame(model = "gaussian_ls", n = 2000, p = 20, num_trees = 50, bandwidth_every = 1, gate = "smoothstep"),
  data.frame(model = "gaussian_ls", n = 400,  p = 20, num_trees = 200, bandwidth_every = 1, gate = "smoothstep"),
  data.frame(model = "gaussian_ls", n = 400,  p = 20, num_trees = 50, bandwidth_every = 10, gate = "smoothstep"),
  data.frame(model = "gaussian_ls", n = 400,  p = 20, num_trees = 50, bandwidth_every = 1, gate = "logistic"),
  data.frame(model = "gaussian_ls", n = 400,  p = 20, num_trees = 50, bandwidth_every = 1, gate = "hard"),
  data.frame(model = "gaussian_ls", n = 400,  p = 100, num_trees = 50, bandwidth_every = 1, gate = "smoothstep"),
  data.frame(model = "zi_poisson",  n = 400,  p = 20, num_trees = 50, bandwidth_every = 1, gate = "smoothstep"),
  data.frame(model = "zi_poisson",  n = 2000, p = 20, num_trees = 50, bandwidth_every = 1, gate = "smoothstep"),
  data.frame(model = "vc",          n = 400,  p = 20, num_trees = 50, bandwidth_every = 1, gate = "smoothstep"),
  data.frame(model = "vc",          n = 2000, p = 20, num_trees = 50, bandwidth_every = 1, gate = "smoothstep"),
  stringsAsFactors = FALSE
)

pr <- prog_init(total = N_FITS, title = "Shared forests: timing",
                unit = "fit", kind = "benchmark")
on.exit(prog_end(pr, "failed"), add = TRUE)

# Written after every grid row and every width rather than once at the end, so
# a run that is killed still leaves its finished rows on disk. `complete` is
# what tells a reader which of the two it has.
grid_rows <- list()
wide_rows <- list()

checkpoint <- function(complete) {
  saveRDS(list(grid = do.call(rbind, grid_rows),
               vc_wide = if (length(wide_rows)) do.call(rbind, wide_rows),
               reps = REPS, complete = complete,
               done = length(grid_rows) + length(wide_rows),
               total = nrow(grid) + 3L), OUT)
}

for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  r <- time_cell(pr, g$model, g$n, g$p, g$num_trees, num_burn = 200L,
                 num_draws = 400L, g$bandwidth_every, g$gate,
                 seed = 4000L + i)
  cat(sprintf("%2d/%d  %s n=%d p=%d trees=%d every=%d %s: separate %.2fs  shared %.2fs  ->  %.2fx\n",
              i, nrow(grid), g$model, g$n, g$p, g$num_trees,
              g$bandwidth_every, g$gate, r$separate, r$shared, r$speedup))
  grid_rows[[i]] <- r
  checkpoint(FALSE)
}

res <- do.call(rbind, grid_rows)

# A varying-coefficient model with several coefficients, where the number of
# forests is larger than two and the saving should be correspondingly bigger.
many_vc <- function(n, p, n_slope, seed) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  d <- as.data.frame(X)
  eta <- 3 * sin(pi * d$x1)

  for (j in seq_len(n_slope)) {
    d[[paste0("z", j)]] <- rnorm(n)
    eta <- eta + d[[paste0("z", j)]] * (1 + d[[paste0("x", j + 1)]])
  }

  d$y <- eta + rnorm(n)
  d
}

for (n_slope in c(1L, 3L, 7L)) {
  d <- many_vc(400L, 20L, n_slope, seed = 8000L + n_slope)
  rhs <- paste(paste0("x", 1:20), collapse = " + ")
  vc <- paste(sprintf("vc(z%d)", seq_len(n_slope)), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs, "+", vc))

  out <- vapply(c(FALSE, TRUE), function(shared) {
    share_forests(shared)
    ctrl <- bartisan_control(num_trees = 50L, num_burn = 200L,
                             num_draws = 400L)
    min(vapply(seq_len(REPS), function(r) {
      secs <- system.time(bartisan(form, data = d, control = ctrl,
                                   family = gaussian()))[["elapsed"]]
      prog_tick(pr, secs = secs,
                label = sprintf("vc %d forests %s rep %d", n_slope + 1L,
                                if (shared) "shared" else "separate", r))
      secs
    }, numeric(1L)))
  }, numeric(1L))

  cat(sprintf("vc with %d coefficients (%d forests): separate %.2fs  shared %.2fs  ->  %.2fx\n",
              n_slope, n_slope + 1L, out[1], out[2], out[1] / out[2]))

  wide_rows[[length(wide_rows) + 1L]] <-
    data.frame(forests = n_slope + 1L, separate = out[1], shared = out[2],
               speedup = out[1] / out[2])
  checkpoint(FALSE)
}

vc_wide <- do.call(rbind, wide_rows)
checkpoint(TRUE)

on.exit()
prog_end(pr, "done", sprintf("%d cells", nrow(grid) + length(wide_rows)))

cat("\nwrote", OUT, "\n")
print(res)
print(vc_wide)
