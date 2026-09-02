# What the sampler costs, with the settings actually matched.
#
# `_dev/benchmark.Rmd` compares bartisan against dbarts at the same tree count
# and the same number of iterations, and with `gate = "hard"` so both use hard
# decision rules. Four things were still not matched, and they pull in both
# directions:
#
#   1. Test predictions. `bart(xtr, ytr, xte, ...)` computes 1000 test points at
#      every draw inside the timed call; bartisan's `predict()` ran after the
#      timer stopped. dbarts was charged for work bartisan was not.
#   2. Retaining trees. dbarts ran with `keeptrees = FALSE`, so it discarded what
#      bartisan serializes on every draw to make `predict()` on new data work.
#      Again dbarts was doing less.
#   3. The sparsity prior. bartisan defaults to DART (Linero 2018), which draws
#      splitting probabilities every sweep. dbarts has no such prior, so bartisan
#      was doing more.
#   4. The leaf scale. bartisan draws `sigma_mu`; dbarts fixes it from `k`. Again
#      bartisan was doing more.
#
# So the honest question is not one number but a decomposition. Every arm below
# fits on training data only and predicts nothing, which isolates the sampler.

library(bartisan)
suppressMessages(library(dbarts))

REPS <- 6L
N <- 1000L
P <- 10L
TREES <- 50L
DRAWS <- 1000L

friedman <- function(n, p, seed) {
  set.seed(seed)
  x <- matrix(stats::runif(n * p), n, p, dimnames = list(NULL, paste0("x", seq_len(p))))
  eta <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
    10 * x[, 4] + 5 * x[, 5]
  list(x = x, eta = eta, y = eta + stats::rnorm(n, sd = 0.5))
}

d <- friedman(N, P, seed = 12)
frame <- data.frame(d[["x"]], y = d[["y"]])

# Best of REPS, which is what a timing claim should rest on: the minimum is the
# least contaminated by whatever else the machine was doing.
#
# `eval.parent(substitute(expr))` rather than `force(expr)`: a promise is
# evaluated once and cached, so forcing it in a loop times the first iteration
# and returns instantly for the rest, and the minimum over that is zero.
best_of <- function(label, expr) {
  e <- substitute(expr)
  seconds <- Inf

  for (i in seq_len(REPS)) {
    t <- system.time(eval(e, parent.frame()))[["elapsed"]]
    seconds <- min(seconds, t)
  }

  cat(sprintf("  %-52s %7.3f s\n", label, seconds))
  invisible(seconds)
}

ctrl <- function(...) {
  bartisan_control(num_trees = TREES, num_burn = DRAWS, num_draws = DRAWS,
                   gate = "hard", verbose = FALSE, ...)
}

# Named explicitly, and that is the point. A numeric response with no `family`
# reaches `default_family()`, which returns `dpm()` -- a Dirichlet process
# mixture for the error distribution, not a Gaussian. `_dev/benchmark.Rmd`
# omitted it in the Gaussian section, so every row labelled "Gaussian" there was
# timing a DPM against three packages fitting a Gaussian.
fit_bartisan <- function(...) {
  bartisan(y ~ ., data = frame, family = gaussian(), ...)
}

cat(sprintf("\nFriedman, n = %d, p = %d, %d trees, %d warmup + %d draws, best of %d\n\n",
            N, P, TREES, DRAWS, DRAWS, REPS))

cat("dbarts\n")
dbarts_notrees <- best_of("bart(), keeptrees = FALSE, no test data",
  dbarts::bart(d[["x"]], d[["y"]], ntree = TREES, nskip = DRAWS,
               ndpost = DRAWS, verbose = FALSE, keeptrees = FALSE))
dbarts_trees <- best_of("bart(), keeptrees = TRUE, no test data",
  dbarts::bart(d[["x"]], d[["y"]], ntree = TREES, nskip = DRAWS,
               ndpost = DRAWS, verbose = FALSE, keeptrees = TRUE))
best_of("bart(), keeptrees = FALSE, WITH test data (as benchmarked)",
  dbarts::bart(d[["x"]], d[["y"]], d[["x"]], ntree = TREES, nskip = DRAWS,
               ndpost = DRAWS, verbose = FALSE, keeptrees = FALSE))
best_of("bart(), usequants = TRUE, keeptrees = TRUE",
  dbarts::bart(d[["x"]], d[["y"]], ntree = TREES, nskip = DRAWS,
               ndpost = DRAWS, verbose = FALSE, keeptrees = TRUE,
               usequants = TRUE))

cat("\nbartisan, hard rules\n")
bart_default <- best_of("gaussian, defaults (sparsity + drawn leaf scale)",
  fit_bartisan(control = ctrl()))
bart_nosparse <- best_of("gaussian, sparsity = FALSE",
  fit_bartisan(control = ctrl(sparsity = FALSE)))
bart_matched <- best_of("gaussian, sparsity = FALSE, update_sigma_mu = FALSE",
  fit_bartisan(control = ctrl(sparsity = FALSE, update_sigma_mu = FALSE)))
bart_soft <- best_of("gaussian, soft rules (the package default gate)",
  bartisan(y ~ ., data = frame, family = gaussian(),
           control = bartisan_control(num_trees = TREES, num_burn = DRAWS,
                                      num_draws = DRAWS, verbose = FALSE)))
bart_dpm <- best_of("dpm, hard rules  <- what benchmark.Rmd timed",
  bartisan(y ~ ., data = frame, control = ctrl()))

cat("\nRatios against dbarts\n")
cat(sprintf("  as benchmarked  (bartisan default / dbarts keeptrees=FALSE)  %.2f\n",
            bart_default / dbarts_notrees))
cat(sprintf("  matched models  (bartisan matched / dbarts keeptrees=TRUE)   %.2f\n",
            bart_matched / dbarts_trees))
cat(sprintf("  what DART costs                                             %.2f x\n",
            bart_default / bart_nosparse))
cat(sprintf("  what drawing the leaf scale costs                           %.2f x\n",
            bart_nosparse / bart_matched))
cat(sprintf("  what retaining trees costs dbarts                           %.2f x\n",
            dbarts_trees / dbarts_notrees))
cat(sprintf("  what the missing `family = gaussian()` cost (dpm / gaussian) %.2f x\n",
            bart_dpm / bart_default))
cat(sprintf("  what soft rules cost over hard                              %.2f x\n",
            bart_soft / bart_default))

if (requireNamespace("BART", quietly = TRUE)) {
  cat("\nBART\n")
  best_of("wbart(), no test data",
    invisible(utils::capture.output(
      BART::wbart(d[["x"]], d[["y"]], ntree = TREES, nskip = DRAWS,
                  ndpost = DRAWS, printevery = DRAWS + 1L))))
}

if (requireNamespace("stochtree", quietly = TRUE)) {
  cat("\nstochtree\n")
  best_of("bart(), num_gfr = 0, no test data",
    stochtree::bart(X_train = d[["x"]], y_train = d[["y"]], num_gfr = 0L,
                    num_burnin = DRAWS, num_mcmc = DRAWS,
                    mean_forest_params = list(num_trees = TREES)))
}
