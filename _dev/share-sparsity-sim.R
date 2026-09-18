# Does pooling the splitting proportions across a family's forests help, and
# what does it cost when the assumption behind it is false?
#
# `share_sparsity = TRUE` says the same predictors are relevant to every
# additive predictor of the family. That is an assumption, and it has two sides
# worth measuring separately:
#
#   agree      the components really do depend on the same predictors, which is
#              where sharing should win, and win more as P grows
#   disagree   they depend on disjoint sets, which is where sharing is wrong and
#              the question is what it costs
#
# Scored on what each family is for: the mean and the log standard deviation for
# `gaussian_ls()`, and the zero probability and the count mean for
# `zi_poisson()`. Root mean squared error against the truth on a held-out set,
# so the numbers are comparable across arms.
#
# Usage: Rscript _dev/share-sparsity-sim.R <reps> <family>
options(parallelly.availableCores.fallback = 4, parallelly.maxWorkers.localhost = Inf)
A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

suppressMessages({
  library(bartisan); library(future)
})
plan(multisession, workers = 4)

ROOT <- "/Users/NoahGreifer/Dropbox/Research/R/bartisan"
args <- commandArgs(TRUE)
REPS <- if (length(args) > 0) as.integer(args[[1]]) else 20L
FAM  <- if (length(args) > 1) args[[2]] else "gaussian_ls"

# Two smooth functions of five predictors each, so that "the same predictors"
# and "disjoint predictors" are the only thing that changes between designs.
f1 <- function(x, j) {
  1.4 * sin(pi * x[, j[1L]] * x[, j[2L]]) + 1.6 * (x[, j[3L]] - 0.5)^2 +
    0.9 * x[, j[4L]] + 0.5 * x[, j[5L]]
}

design <- function(n, p, agree, seed) {
  set.seed(seed)
  x <- matrix(stats::runif(n * p), n, p)
  colnames(x) <- paste0("x", seq_len(p))
  ja <- 1:5
  jb <- if (agree) 1:5 else 6:10
  list(x = x, a = f1(x, ja), b = f1(x, jb))
}

run_gaussian_ls <- function(train, test, share) {
  # Mean from one set, log standard deviation from the other.
  mu <- 2 + train[["a"]]
  lsd <- -0.7 + 0.5 * train[["b"]]
  d <- data.frame(y = stats::rnorm(nrow(train[["x"]]), mu, exp(lsd)), train[["x"]])

  fit <- bartisan(y ~ ., d, family = gaussian_ls(), share_sparsity = share,
                  num_trees = 50L, num_burn = 500L, num_draws = 1000L,
                  verbose = FALSE)

  nd <- data.frame(test[["x"]])
  eta <- predict(fit, newdata = nd, type = "link")

  c(mean = sqrt(mean((eta[, 1L] - (2 + test[["a"]]))^2)),
    log_sd = sqrt(mean((eta[, 2L] - (-0.7 + 0.5 * test[["b"]]))^2)))
}

run_zi_poisson <- function(train, test, share) {
  # The zero-inflation probability from one set, the count mean from the other.
  # `zi_poisson()`'s first predictor is the count log mean and the second the
  # logit of the zero-inflation probability.
  lam <- 1.1 + train[["a"]]
  zi <- -0.6 + train[["b"]]
  n <- nrow(train[["x"]])
  y <- ifelse(stats::rbinom(n, 1L, stats::plogis(zi)) == 1L, 0L,
              stats::rpois(n, exp(lam)))
  d <- data.frame(y = y, train[["x"]])

  fit <- bartisan(y ~ ., d, family = zi_poisson(), share_sparsity = share,
                  num_trees = 50L, num_burn = 500L, num_draws = 1000L,
                  verbose = FALSE)

  nd <- data.frame(test[["x"]])
  eta <- predict(fit, newdata = nd, type = "link")

  c(count = sqrt(mean((eta[, 1L] - (1.1 + test[["a"]]))^2)),
    zero = sqrt(mean((eta[, 2L] - (-0.6 + test[["b"]]))^2)))
}

runner <- switch(FAM, gaussian_ls = run_gaussian_ls, zi_poisson = run_zi_poisson,
                 stop("family?"))

cells <- expand.grid(p = c(5L, 25L, 100L), agree = c(TRUE, FALSE),
                     KEEP.OUT.ATTRS = FALSE)
# With disjoint sets the truth needs ten predictors, so P = 5 cannot hold it.
cells <- cells[!(cells[["p"]] == 5L & !cells[["agree"]]), ]

OUT <- file.path(ROOT, sprintf("_dev/share-sparsity-%s.rds", FAM))
n_total <- nrow(cells) * REPS * 2L

pr <- prog_init(total = n_total, unit = "fit", kind = "simulation",
                title = sprintf("Shared sparsity: %s", FAM))
on.exit(prog_end(pr, "failed", "aborted before the last replicate"),
        add = TRUE)

out <- list()
done <- 0L

# Written after every replicate rather than at the end, so a run that is killed
# leaves its finished replicates readable. `complete` is what tells a reader
# which of the two it is looking at.
checkpoint <- function(complete) {
  saveRDS(list(res = do.call(rbind, out), complete = complete, done = done,
               total = n_total, reps = REPS, family = FAM, cells = cells), OUT)
}

for (i in seq_len(nrow(cells))) {
  p <- cells[["p"]][i]
  agree <- cells[["agree"]][i]

  for (r in seq_len(REPS)) {
    seed <- 7000L * i + r
    train <- design(400L, p, agree, seed)
    test <- design(1000L, p, agree, seed + 400000L)

    for (share in c(FALSE, TRUE)) {
      set.seed(seed + 900000L)
      t0 <- Sys.time()
      got <- runner(train, test, share)
      out[[length(out) + 1L]] <- data.frame(
        p = p, agree = agree, rep = r, share = share,
        component = names(got), rmse = as.vector(got))

      done <- done + 1L
      prog_tick(pr, i = done,
                secs = as.numeric(difftime(Sys.time(), t0, units = "secs")),
                label = sprintf("p=%d agree=%s %s rep %d", p, agree,
                                if (share) "shared" else "separate", r))
    }

    checkpoint(FALSE)
  }
}

checkpoint(TRUE)
res <- do.call(rbind, out)

on.exit()
prog_end(pr, "done", sprintf("%d fits", done))

cat(sprintf("\n%s: RMSE on the held-out truth, mean over %d reps, n = 400\n\n",
            FAM, REPS))
agg <- aggregate(rmse ~ p + agree + component + share, data = res, FUN = mean)
w <- reshape(agg, idvar = c("p", "agree", "component"), timevar = "share",
             direction = "wide")
names(w) <- sub("^rmse\\.", "share_", names(w))
w[["ratio"]] <- w[["share_FALSE"]] / w[["share_TRUE"]]
w <- w[order(w[["agree"]], w[["component"]], w[["p"]]), ]
print(w, row.names = FALSE, digits = 3)
cat("\nratio > 1 means sharing helped; < 1 means it hurt.\n")
cat("DONE\n")
