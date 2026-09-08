# Shared tree topology: does it cost anything?
#
# Sharing the topology gives every additive predictor of a multi-forest family
# the same trees, so the components partition the covariate space identically
# and differ only in the value each leaf carries. That is a real restriction,
# and the question this script answers is what it costs when it is wrong and
# what it buys when it is right.
#
# Three models, three shapes of forest:
#
#   gaussian_ls()  two parameters, one strong (the mean) and one weak (the log
#                  standard deviation)
#   zi_poisson()   two parameters, likewise -- the count mean and the
#                  zero-inflation probability
#   vc()           two forests feeding *one* parameter: a control function and a
#                  varying coefficient, which is the case sharing was not
#                  designed for and so the one worth checking
#
# and two data-generating shapes: the components driven by the same three
# predictors, and by disjoint sets of three. The disjoint case is the adversary.
# Sharing says the components split on the same predictors, and there the truth
# says they do not, so anything sharing costs shows up there at full size.
#
# Errors are decomposed rather than only averaged. The test set is held fixed
# across replicates, so at each test point the spread of the fitted value over
# replicates is variance and its distance from the truth is bias, and a
# restriction that hurts can be told apart from one that helps. Timings here are
# taken under contention and are not comparable; `_dev/shared-topology-timing.R`
# measures those on a quiet machine.

library(bartisan)

# Sharing the topology is not an argument: it is engaged by the package's
# internal flag, for the reasons recorded beside it in `R/control.R`. These
# scripts are the record of how it was measured, so they reach it the same way
# the tests do.
share_forests <- function(on) {
  flags <- bartisan:::the
  flags$share_forests <- isTRUE(on)
}
library(future.apply)

# Live progress, per the `/live-progress` skill. `prog_future_lapply()` ticks
# from inside each worker, which is what makes a static-chunked `future` run
# visible at all: the workers otherwise return nothing until their whole chunk
# is finished, and an earlier version of this script printed a single line and
# then went silent for forty-four minutes.
source(path.expand(Sys.getenv(
  "LIVE_PROGRESS", "~/.claude/skills/live-progress/assets/progress.R")))

N_TRAIN <- 400L
N_TEST <- 1000L
REPS <- 20L
P_GRID <- c(10L, 50L)

NUM_TREES <- 50L
NUM_BURN <- 500L
NUM_DRAWS <- 1000L

OUT <- file.path("_dev", "shared-topology-sim.rds")

# The strong component: a Friedman-flavored surface over three predictors, with
# enough curvature that a forest has something to find.
f_strong <- function(X, a) {
  3 * sin(pi * X[, a[1]]) + 2 * X[, a[2]] + 3 * (X[, a[3]] - 0.5)^2
}

# The weak one: additive and shallow. The point of it is that a forest looking
# only at its own residual signal has a hard time telling which predictors it
# involves, which is the situation sharing is supposed to rescue.
f_weak <- function(X, b) {
  0.9 * X[, b[1]] + 0.7 * X[, b[2]] - 0.6 * X[, b[3]]
}

active <- function(overlap) {
  list(a = 1:3, b = if (identical(overlap, "same")) 1:3 else 4:6)
}

design <- function(n, p) {
  X <- matrix(runif(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  as.data.frame(X)
}

# One cell's data: a fixed test set, drawn once per cell so that the bias and
# variance below are taken at the same points every replicate, and a training
# set that changes with the replicate.
make_data <- function(model, overlap, p, n, seed) {
  set.seed(seed)
  who <- active(overlap)
  d <- design(n, p)
  X <- as.matrix(d)
  strong <- f_strong(X, who$a)
  weak <- f_weak(X, who$b)

  if (identical(model, "gaussian_ls")) {
    eta <- cbind(mean = strong, log_sd = -0.9 + weak)
    d$y <- eta[, 1] + rnorm(n, 0, exp(eta[, 2]))
    return(list(data = d, eta = eta,
                response = eta[, 1]))
  }

  if (identical(model, "zi_poisson")) {
    eta <- cbind(count = 0.4 + 0.55 * strong, zero = -0.6 + 1.6 * weak)
    zero <- runif(n) < plogis(eta[, 2])
    d$y <- ifelse(zero, 0L, rpois(n, exp(eta[, 1])))
    return(list(data = d, eta = eta,
                response = (1 - plogis(eta[, 2])) * exp(eta[, 1])))
  }

  # The varying-coefficient model. `z` is the covariate whose coefficient
  # varies; both forests are functions of the same x, and "disjoint" means the
  # control function and the coefficient involve different ones.
  eta <- cbind(control = strong, coef = 1 + 1.5 * weak)
  d$z <- rnorm(n)
  d$y <- eta[, 1] + d$z * eta[, 2] + rnorm(n)
  list(data = d, eta = eta, response = eta[, 1] + d$z * eta[, 2])
}

formula_of <- function(model, p) {
  rhs <- paste(paste0("x", seq_len(p)), collapse = " + ")

  if (identical(model, "vc")) {
    # `z` is deliberately left out of the fixed part: with it in, the control
    # function and the coefficient are not separately identified.
    return(stats::as.formula(paste("y ~", rhs, "+ vc(z)")))
  }

  stats::as.formula(paste("y ~", rhs))
}

family_of <- function(model) {
  switch(model,
         gaussian_ls = gaussian_ls(),
         zi_poisson = zi_poisson(),
         vc = gaussian())
}

# The two component surfaces a fit implies, on the scales `make_data()` built
# them on. For a varying-coefficient model the two forests are not separate
# columns of the predictor -- the model is one parameter -- so they are
# recovered by evaluating it at two values of the varying covariate, which the
# identity link makes exact.
components <- function(fit, model, newdata) {
  if (!identical(model, "vc")) {
    return(unname(as.matrix(predict(fit, newdata = newdata, type = "link"))))
  }

  at0 <- transform(newdata, z = 0)
  at1 <- transform(newdata, z = 1)
  e0 <- as.numeric(predict(fit, newdata = at0, type = "link"))
  e1 <- as.numeric(predict(fit, newdata = at1, type = "link"))
  cbind(e0, e1 - e0)
}

# Average splits per draw in each forest, which says whether sharing changes how
# complicated the trees get as well as whether they agree.
splits_of <- function(fit) {
  vapply(fit[["counts"]], function(m) mean(rowSums(m)), numeric(1L))
}

one_fit <- function(model, p, dat, shared) {
  share_forests(shared)
  ctrl <- bartisan_control(num_trees = NUM_TREES, num_burn = NUM_BURN,
                           num_draws = NUM_DRAWS)
  t <- system.time(fit <- bartisan(formula_of(model, p), data = dat[["data"]],
                                   family = family_of(model), control = ctrl))
  list(fit = fit, elapsed = t[["elapsed"]])
}

# One replicate of one cell, fitting both ways on the same data so the
# comparison is paired.
#
# `test` is passed in rather than regenerated here, and that is not a
# convenience. `set.seed()` keeps whatever RNG kind is current, and
# `future.seed = TRUE` puts every worker on L'Ecuyer-CMRG while the main
# process stays on Mersenne-Twister -- so the same seed builds one dataset in a
# worker and a different one outside. Regenerating the test set to score
# against scored each fit on a different draw of the design, which shows up as
# fitted values with the right marginal distribution and no correlation with
# the truth at all. The test sets are built once, in the main process, and
# handed over.
one_rep <- function(cell, rep, test) {
  model <- cell[["model"]]
  overlap <- cell[["overlap"]]
  p <- cell[["p"]]

  train <- make_data(model, overlap, p, N_TRAIN, seed = 500000L + rep * 977L +
                       17L * p + 100L * match(overlap, c("same", "disjoint")) +
                       1000L * match(model, c("gaussian_ls", "zi_poisson", "vc")))

  out <- list()

  for (shared in c(FALSE, TRUE)) {
    got <- one_fit(model, p, train, shared)
    eta <- components(got[["fit"]], model, test[["data"]])
    resp <- as.numeric(predict(got[["fit"]], newdata = test[["data"]],
                               type = "response"))

    out[[length(out) + 1L]] <- list(
      model = model, overlap = overlap, p = p, rep = rep, shared = shared,
      elapsed = got[["elapsed"]],
      splits = splits_of(got[["fit"]]),
      eta = eta,
      resp = resp
    )
  }

  out
}

cells <- expand.grid(model = c("gaussian_ls", "zi_poisson", "vc"),
                     overlap = c("same", "disjoint"),
                     p = P_GRID,
                     stringsAsFactors = FALSE)

jobs <- do.call(rbind, lapply(seq_len(REPS), function(r) {
  cbind(cells, rep = r, stringsAsFactors = FALSE)
}))

cat("cells:", nrow(cells), " jobs:", nrow(jobs), "\n")

# The test sets, built once here so that the fits and the truth they are scored
# against are the same draw of the design. Keyed the way the scoring reads them.
cell_key <- function(model, overlap, p) paste(model, overlap, p)

truth <- list()

for (i in seq_len(nrow(cells))) {
  key <- cell_key(cells[["model"]][i], cells[["overlap"]][i], cells[["p"]][i])
  truth[[key]] <- make_data(cells[["model"]][i], cells[["overlap"]][i],
                            cells[["p"]][i], N_TEST,
                            seed = 10000L + 17L * cells[["p"]][i] +
                              100L * match(cells[["overlap"]][i],
                                           c("same", "disjoint")) +
                              1000L * match(cells[["model"]][i],
                                            c("gaussian_ls", "zi_poisson",
                                              "vc")))
}

# Workers from the environment, so the same script serves a dry run in one
# process and the real thing on however many cores are free.
workers <- as.integer(Sys.getenv("SHARED_WORKERS", "8"))

if (is.na(workers) || workers < 2L) {
  future::plan(future::sequential)
} else {
  future::plan(future::multisession, workers = workers)
}

pr <- prog_init(total = nrow(jobs), title = "Shared forests: accuracy",
                unit = "job", workers = max(workers, 1L), kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

res <- prog_future_lapply(pr, seq_len(nrow(jobs)), function(j) {
  cell <- as.list(jobs[j, c("model", "overlap", "p")])
  one_rep(cell, jobs[["rep"]][j],
          truth[[cell_key(cell[["model"]], cell[["overlap"]], cell[["p"]])]])
}, label = function(j) sprintf("%s %s p=%d rep %d", jobs[["model"]][j],
                               jobs[["overlap"]][j], jobs[["p"]][j],
                               jobs[["rep"]][j]))

# `prog_do()` turns a failed job into a NULL rather than killing the run, so a
# short result list means jobs died; the widget's failures panel says which.
bad <- sum(vapply(res, is.null, logical(1L)))

if (bad > 0) {
  prog_log(pr, sprintf("%d of %d jobs failed", bad, nrow(jobs)), level = "warn")
}

saveRDS(list(res = unlist(res, recursive = FALSE), truth = truth, cells = cells,
             settings = list(n_train = N_TRAIN, n_test = N_TEST, reps = REPS,
                             num_trees = NUM_TREES, num_burn = NUM_BURN,
                             num_draws = NUM_DRAWS)),
        OUT)

on.exit()
prog_end(pr, if (bad > 0) "failed" else "done")

cat("wrote", OUT, "\n")
