# Is the soft-rule Poisson SBC deviation in Poisson regression, or in bartisan?
#
# What is tested. `_dev/sbc.R`'s soft-rule Poisson arm gives posteriors slightly
# too narrow whenever the bandwidth is drawn: a dispersion contrast near +2 on
# 600 replicates, across bandwidth priors, warmup, chain length, response level
# and thinning, where the arm with the bandwidth fixed and the soft-rule logit
# arms are uniform (`_dev/TASKS.md`, 2026-10-01 and 2026-10-02). Every part of
# bartisan's sampler read since is exact as written, and mixing was excluded on
# 2026-10-01 (`_dev/sbc-mixing.R`). The claim here is that the deviation belongs
# to bartisan's sampler and not to the problem: that a general-purpose sampler
# known to be correct, given the same data, the same priors and the same
# harness, produces uniform ranks.
#
# What is measured. The generator of `_dev/sbc.R`, its functions evaluated from
# that file, with its default seed block, n = 400 on its fixed design, and the
# soft-rule Poisson settings, so each replicate here is the same data set and
# the same truth as the replicate of that number there. The fit is Stan's NUTS
# (`_dev/sbc-stan.stan`, through *rstan*) on the same model with the tree
# structures given as data: leaf values N(0, 0.35^2) and one bandwidth per tree
# from an exponential prior with mean 0.1, the smoothstep gate of `node.h`.
# Conditioning on part of the truth as data keeps the ranks uniform for a
# correct sampler, so the test is valid; what it leaves out is the tree moves.
# Per replicate: the rank of the true contrast eta[A] - eta[B] among 100 draws
# thinned from the kept ones, as `_dev/sbc.R` takes it; the contrast's bulk
# effective sample size, and the count of divergent transitions, which say
# whether Stan's own draws can be trusted. NUTS takes minutes a replicate on
# this model, so the run is replicates 1 to 300 with 500 warmup and 500 kept
# draws. Read the dispersion contrast and the corrected chi-square over them,
# against bartisan's ranks on the same 300 data sets, matched by the true
# contrast, which the two generators produce identically.
#
# What each outcome would mean. Uniform ranks (a dispersion contrast near zero)
# mean Poisson regression with drawn bandwidths calibrates under a correct
# sampler, so the deviation is bartisan's, and with the leaf and bandwidth
# updates verified the tree moves are where it lives. A dispersion contrast
# near bartisan's means the deviation is a property of the problem or of the
# harness, and not a defect of the package.
#
# Run with: Rscript _dev/sbc-stan.R <first replicate> <last replicate> [offset]
#   STAN_WARMUP and STAN_DRAWS change the chain, default 500 and 500.
# Writes:   _dev/sbc-stan-off<offset>-w<warmup>-d<draws>-r<first>-<last>.rds

A_ <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A_, "progress.R"))

suppressPackageStartupMessages(library(rstan))
rstan_options(auto_write = TRUE)

args <- commandArgs(trailingOnly = TRUE)
FIRST <- as.integer(args[1L])
LAST <- as.integer(args[2L])
OFFSET <- if (length(args) > 2L) as.numeric(args[3L]) else 0
SEED_BASE <- 5000L
L <- 100L
# Stan's draws of the contrast came out nearly independent (an effective sample
# size of about 1000 per 1000) on the first replicates, so 500 kept draws still
# leave five times what 100 ranks need. Both counts are in the file name.
WARMUP <- as.integer(Sys.getenv("STAN_WARMUP", "500"))
DRAWS <- as.integer(Sys.getenv("STAN_DRAWS", "500"))
OUT <- sprintf("_dev/sbc-stan-off%s-w%d-d%d-r%d-%d.rds", format(OFFSET), WARMUP,
               DRAWS, FIRST, LAST)

# The generator, as `_dev/sbc.R` has it for this arm.
P <- 2L
TREES <- 20L
SIGMA_MU <- 0.35
GAMMA <- 0.95
BETA <- 2
BANDWIDTH <- 0.1
SOFT <- TRUE
FIX_BANDWIDTH <- FALSE
N <- 400L
LAST_BW <- numeric(0)

wanted <- c("HALF_WIDTH", "draw_tree", "left_prob", "eval_tree", "draw_forest")
for (e in parse("_dev/sbc.R")) {
  if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.name(e[[2L]]) &&
      as.character(e[[2L]]) %in% wanted) {
    eval(e, globalenv())
  }
}
stopifnot(all(vapply(wanted, exists, logical(1L))))

set.seed(20260905)
u <- matrix(stats::runif(N * P), N, P)
colnames(u) <- paste0("x", seq_len(P))
A <- which.min(u[, 1L])
B <- which.max(u[, 1L])

# `draw_forest()`, keeping the trees: the same draws in the same order, a
# bandwidth and then a tree for each of the twenty.
draw_trees <- function() {
  lapply(seq_len(TREES), function(i) {
    b <- stats::rexp(1L, rate = 1 / BANDWIDTH)
    list(bandwidth = b, tree = draw_tree())
  })
}

# Every leaf as the path that reaches it.
leaf_paths <- function(node, path = list()) {
  if (isTRUE(node$leaf)) {
    return(list(list(mu = node$mu, path = path)))
  }
  c(leaf_paths(node$left, c(path, list(c(node$var, node$val, 1)))),
    leaf_paths(node$right, c(path, list(c(node$var, node$val, 0)))))
}

stan_data <- function(forest, y) {
  leaves <- do.call(c, lapply(seq_along(forest), function(t) {
    lapply(leaf_paths(forest[[t]]$tree), function(l) c(l, tree = t))
  }))
  K <- length(leaves)
  depth <- vapply(leaves, function(l) length(l$path), integer(1L))
  D <- max(1L, depth)
  var <- matrix(1L, K, D)
  val <- matrix(0, K, D)
  left <- matrix(1L, K, D)
  for (k in seq_len(K)) {
    for (s in seq_len(depth[k])) {
      step <- leaves[[k]]$path[[s]]
      var[k, s] <- as.integer(step[1L])
      val[k, s] <- step[2L]
      left[k, s] <- as.integer(step[3L])
    }
  }
  list(N = N, P = P, X = u, y = as.integer(y), log_offset = rep(OFFSET, N),
       T = TREES, K = K, D = D,
       leaf_tree = vapply(leaves, function(l) as.integer(l$tree), integer(1L)),
       depth = depth, split_var = var, split_val = val, goes_left = left,
       sigma_mu = SIGMA_MU, bandwidth_mean = BANDWIDTH,
       half_width = HALF_WIDTH, A = A, B = B)
}

model <- stan_model("_dev/sbc-stan.stan")

rows <- list()
start <- FIRST
if (file.exists(OUT)) {
  prev <- readRDS(OUT)
  if (!isTRUE(prev$complete) && NROW(prev$rows)) {
    rows <- split(prev$rows, seq_len(nrow(prev$rows)))
    start <- max(prev$rows$rep) + 1L
    cat(sprintf("resuming at replicate %d\n", start))
  }
}

pr <- prog_init(total = LAST - FIRST + 1L,
                title = sprintf("SBC with Stan: soft Poisson, offset %s, %d-%d",
                                format(OFFSET), FIRST, LAST),
                unit = "replicate", kind = "simulation")
if (start > FIRST) prog_tick(pr, i = start - FIRST, label = "resumed")
on.exit(prog_end(pr, "failed"), add = TRUE)

for (r in seq(start, LAST)) {
  set.seed(SEED_BASE + r)
  forest <- draw_trees()
  eta <- Reduce(`+`, lapply(forest, function(f) eval_tree(f$tree, u, f$bandwidth)),
                numeric(N))

  # The same draws as `draw_forest()` under the same seed, or this is not the
  # replicate of that number in `_dev/sbc.R`. Both consume the same random
  # numbers, so the stream is where it would have been either way.
  if (r == start) {
    set.seed(SEED_BASE + r)
    stopifnot(identical(eta, draw_forest(u)))
  }

  truth <- eta[A] - eta[B]
  y <- stats::rpois(N, exp(eta + OFFSET))

  if (length(unique(y)) < 2L) {
    prog_tick(pr, ok = FALSE, label = sprintf("replicate %d", r))
    next
  }

  sd <- stan_data(forest, y)
  t0 <- Sys.time()
  fit <- sampling(model, data = sd, chains = 1L,
                  warmup = WARMUP, iter = WARMUP + DRAWS, refresh = 0L,
                  seed = SEED_BASE + r, show_messages = FALSE)

  contrast <- as.vector(rstan::extract(fit, "contrast")[[1L]])
  thin <- contrast[seq(1L, length(contrast), length.out = L)]
  ci <- stats::quantile(contrast, c(0.025, 0.975), names = FALSE)

  rows[[length(rows) + 1L]] <- data.frame(
    rep = r, rank = sum(thin < truth), truth = truth,
    covered = as.numeric(ci[1L] <= truth && truth <= ci[2L]),
    ess = posterior::ess_bulk(matrix(contrast, ncol = 1L)),
    divergent = rstan::get_num_divergent(fit),
    leaves = sd$K,
    steps = mean(rstan::get_num_leapfrog_per_iteration(fit)),
    secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))

  prog_tick(pr, label = sprintf("replicate %d, rank %d", r, rows[[length(rows)]]$rank))
  saveRDS(list(rows = do.call(rbind, rows), complete = r == LAST), OUT)
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
cat(sprintf("\n%d replicates; median ESS %.0f of 1000; %d with divergences\n",
            nrow(res), stats::median(res$ess), sum(res$divergent > 0)))
