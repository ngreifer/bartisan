# Is a formal variable-selection test missing, and would the permutation one buy
# anything?
#
# `vignette("implementation")` lists "formal variable-selection test" as a gap,
# filled by a helper package, and names two things the other packages do:
# *SoftBart* reports posterior inclusion probabilities and *bartMachine*
# permutes the response. Reading the source settles the first half without any
# fitting -- `SoftBart::posterior_probs()` is
#
#     post_probs <- colMeans(fit$var_counts > 0)
#     median_probability_model <- which(post_probs > 0.5)
#
# which is the `prop_used` column of `variable_importance()` and a cut at .5.
# So this script is about the second half: does the permutation null of Bleich
# et al. (2014) select a different set than the median probability model, and
# does either hold its size?
#
# The comparison is deliberately unfair to the permutation test in one way and
# to the median probability model in another. `bartMachine` has no sparsity
# prior, so its inclusion proportions are what a uniform splitting prior
# produces and its threshold has to come from somewhere; the Dirichlet prior
# does the shrinking inside the model instead. Both arms are therefore run both
# ways, so the question "is the test doing the work, or is the prior" has an
# answer rather than a confound.
#
# Four selectors, all reading the same fitted counts:
#
#   mpm      prop_used > .5, the median probability model
#   local    prop_splits above its own predictor's 1-alpha permutation quantile
#   max      prop_splits above the 1-alpha quantile of the permutation row maxima
#   se       prop_splits above perm mean + c * perm SD, c bisected for
#            1-alpha simultaneous coverage -- bartMachine's third procedure
#
# Two data-generating processes: three real predictors among p, and a null where
# none of them matter, which is where a selector's size is visible at all.
A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))
Sys.setenv(PROGRESS_ROOT = "/Users/NoahGreifer/Dropbox/Research/R/bartisan/.progress-runs")

suppressMessages(library(bartisan))

REPS <- 20L
N <- 300L
P <- 25L
NPERM <- 40L        # bartMachine's default is 100; 40 keeps 1600 fits affordable
TREES <- 20L        # bartMachine's `num_trees_for_permute`, and its reason
ALPHA <- 0.05
TRUE_VARS <- c("x1", "x2", "x3")
OUT <- file.path("_dev", "varsel-check.rds")

ctrl <- function(sparsity) {
  bartisan_control(num_trees = TREES, num_burn = 150L, num_draws = 300L,
                   sparsity = sparsity)
}

make <- function(seed, signal) {
  set.seed(seed)
  X <- matrix(runif(N * P), N, P)
  colnames(X) <- paste0("x", seq_len(P))
  d <- as.data.frame(X)
  mu <- if (signal) {
    2 * sin(pi * d$x1) + 1.5 * d$x2 + 2 * (d$x3 - 0.5)^2
  } else {
    rep(0, N)
  }
  d$y <- mu + rnorm(N, 0, 0.7)
  d
}

rhs <- paste(paste0("x", seq_len(P)), collapse = " + ")
form <- stats::as.formula(paste("y ~", rhs))

quiet <- function(expr) {
  invisible(utils::capture.output(out <- suppressMessages(suppressWarnings(expr))))
  out
}

# The two functionals of the split counts the two literatures use. `used` is
# SoftBart's inclusion probability; `share` is bartMachine's inclusion
# proportion, which is a share of the forest's rules rather than a probability.
counts_of <- function(fit) {
  m <- variable_importance(fit, draws = TRUE)
  tot <- rowSums(m)
  list(used = colMeans(m > 0),
       share = colMeans(m / ifelse(tot > 0, tot, NA_real_), na.rm = TRUE))
}

# bartMachine's `bisectK`: the multiplier that makes mean + c * SD cover
# 1 - alpha of the permutation draws simultaneously across predictors.
bisect_k <- function(perm, coverage, tol = 0.01, lo = 1, hi = 20, limit = 100) {
  m <- colMeans(perm)
  s <- apply(perm, 2L, stats::sd)
  cover <- function(k) mean(apply(perm, 1L, function(r) all(r <= m + k * s)))
  for (i in seq_len(limit)) {
    mid <- (lo + hi) / 2
    cv <- cover(mid)
    if (abs(cv - coverage) < tol) return(mid)
    if (cv < coverage) lo <- mid else hi <- mid
  }
  mid
}

select_all <- function(obs, perm, alpha) {
  nm <- names(obs$share)
  q <- apply(perm, 2L, stats::quantile, probs = 1 - alpha)
  max_cut <- stats::quantile(apply(perm, 1L, max), probs = 1 - alpha)
  k <- bisect_k(perm, coverage = 1 - alpha)
  thr_se <- colMeans(perm) + k * apply(perm, 2L, stats::sd)
  list(mpm = nm[obs$used > 0.5],
       local = nm[obs$share > q & obs$share > 0],
       max = nm[obs$share >= max_cut & obs$share > 0],
       se = nm[obs$share >= thr_se & obs$share > 0])
}

cells <- expand.grid(sparsity = c(TRUE, FALSE), signal = c(TRUE, FALSE),
                     KEEP.OUT.ATTRS = FALSE)

pr <- prog_init(total = nrow(cells) * REPS,
                title = "Permutation null vs the median probability model",
                unit = "replicate", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
k <- 0L

for (i in seq_len(nrow(cells))) {
  sp <- cells$sparsity[i]
  sig <- cells$signal[i]

  for (r in seq_len(REPS)) {
    t0 <- Sys.time()
    d <- make(50000L + 137L * r + (if (sig) 0L else 7L), signal = sig)

    fit <- quiet(bartisan(form, d, family = stats::gaussian(),
                          control = ctrl(sp)))
    obs <- counts_of(fit)

    # The null: the same fit on a response with its link to X broken. Only the
    # response is permuted, so the predictors keep their joint distribution and
    # their correlations, which is what makes the null the right one.
    perm <- matrix(NA_real_, NPERM, P, dimnames = list(NULL, paste0("x", 1:P)))
    for (b in seq_len(NPERM)) {
      db <- d
      db$y <- sample(d$y)
      pf <- quiet(bartisan(form, db, family = stats::gaussian(),
                           control = ctrl(sp)))
      perm[b, ] <- counts_of(pf)$share
    }

    sel <- select_all(obs, perm, ALPHA)

    k <- k + 1L
    rows[[k]] <- do.call(rbind, lapply(names(sel), function(nm) {
      s <- sel[[nm]]
      data.frame(
        sparsity = sp, signal = sig, rep = r, method = nm,
        n_selected = length(s),
        n_true = if (sig) length(intersect(s, TRUE_VARS)) else NA_integer_,
        n_false = length(setdiff(s, if (sig) TRUE_VARS else character())),
        any_false = length(setdiff(s, if (sig) TRUE_VARS else character())) > 0,
        stringsAsFactors = FALSE)
    }))

    prog_tick(pr, i = k, label = sprintf(
      "sparsity=%s signal=%s r%d (%.0fs)", sp, sig, r,
      as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }

  saveRDS(list(res = do.call(rbind, rows), reps = REPS, n = N, p = P,
               nperm = NPERM, trees = TREES, alpha = ALPHA,
               true_vars = TRUE_VARS, complete = FALSE,
               done = i, total = nrow(cells)), OUT)
}

res <- do.call(rbind, rows)
saveRDS(list(res = res, reps = REPS, n = N, p = P, nperm = NPERM,
             trees = TREES, alpha = ALPHA, true_vars = TRUE_VARS,
             complete = TRUE, done = nrow(cells), total = nrow(cells)), OUT)
on.exit()
prog_end(pr, "done")

key <- paste(res$sparsity, res$signal, res$method)
out <- do.call(rbind, lapply(split(res, key), function(z) {
  data.frame(sparsity = z$sparsity[1], signal = z$signal[1],
             method = z$method[1],
             power = if (isTRUE(z$signal[1])) mean(z$n_true) / 3 else NA_real_,
             false_pos = mean(z$n_false),
             any_false = mean(z$any_false),
             selected = mean(z$n_selected),
             stringsAsFactors = FALSE)
}))
out <- out[order(out$signal, out$sparsity, out$method), ]
print(out, row.names = FALSE, digits = 3)
cat("\n", REPS, " reps, n = ", N, ", p = ", P, ", ", TREES, " trees, ",
    NPERM, " permutations, alpha = ", ALPHA,
    "\n`any_false` on the null rows is the family-wise error rate.\n", sep = "")
cat("wrote", OUT, "\n")
