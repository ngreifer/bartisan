# Should `sparsity = TRUE` stay the default?
#
# The atom-at-zero finding says the sparsity prior costs something for a
# contrast estimand. It does not say what the prior buys elsewhere, and the
# default should turn on the trade rather than on either half of it. Two
# questions, measured separately because they have different answers:
#
#   (1) Prediction. Does concentrating splits improve accuracy, and does the
#       answer depend on how many irrelevant predictors there are?
#   (2) A contrast. What does each setting do to the ATE: its bias, its interval,
#       and the mass the prior puts at exactly zero?
#
# Everything below uses default `num_burn` and `num_draws` so that this script
# does not depend on their names.

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)
library(marginaleffects)

LEVELS <- c("none", "weak", "moderate", "strong")
REPS <- 3
P_GRID <- c(10L, 50L)
OUT <- "_dev/sparsity-default.rds"

# Four sparsity settings per replicate, at each of two widths, in each of the
# two sections.
N_FITS <- 2L * length(P_GRID) * REPS * length(LEVELS)

pr <- prog_init(total = N_FITS, title = "The sparsity default", unit = "fit",
                kind = "simulation")
on.exit(prog_end(pr, "failed", "aborted before the last replicate"),
        add = TRUE)

results <- list()
done <- 0L

# Written after every replicate rather than at the end, so a run that is killed
# leaves its finished replicates readable. `complete` is what tells a reader
# which of the two it is looking at.
checkpoint <- function(complete) {
  saveRDS(list(res = results, complete = complete, done = done,
               total = N_FITS, reps = REPS, levels = LEVELS), OUT)
}

friedman <- function(x) {
  10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
    10 * x[, 4] + 5 * x[, 5]
}

# ---- (1) prediction ------------------------------------------------------

predict_rep <- function(p, rep) {
  set.seed(1000 * p + rep)
  n <- 500
  xtr <- matrix(runif(n * p), n, p, dimnames = list(NULL, paste0("x", seq_len(p))))
  xte <- matrix(runif(1000 * p), 1000, p, dimnames = list(NULL, colnames(xtr)))
  ftr <- friedman(xtr)
  fte <- friedman(xte)
  dtr <- data.frame(y = ftr + rnorm(n), xtr)
  dte <- as.data.frame(xte)

  vapply(LEVELS, function(lv) {
    set.seed(7)
    t0 <- Sys.time()
    fit <- bartisan(y ~ ., data = dtr, family = gaussian(), chains = 2,
                    control = bartisan_control(sparsity = lv, verbose = FALSE))
    done <<- done + 1L
    prog_tick(pr, i = done,
              secs = as.numeric(difftime(Sys.time(), t0, units = "secs")),
              label = sprintf("prediction p=%d rep %d / %s", p, rep, lv))
    # Scored against the true regression function, not the noisy outcome, so
    # this is error in the fit rather than irreducible noise.
    sqrt(mean((predict(fit, newdata = dte) - fte)^2))
  }, numeric(1))
}

cat("== Prediction: RMSE against the true function, mean of", REPS, "reps ==\n")
cat(sprintf("%-8s %s\n", "p", paste(sprintf("%8s", LEVELS), collapse = "")))
for (p in P_GRID) {
  key <- sprintf("prediction, p = %d", p)
  each <- matrix(NA_real_, length(LEVELS), REPS,
                 dimnames = list(LEVELS, NULL))

  for (r in seq_len(REPS)) {
    each[, r] <- predict_rep(p, r)
    results[[key]] <- each
    checkpoint(FALSE)
  }

  cat(sprintf("%-8d %s\n", p,
              paste(sprintf("%8.3f", rowMeans(each)), collapse = "")))
}

# ---- (2) a contrast ------------------------------------------------------

TAU <- 0.5

contrast_rep <- function(p, rep) {
  set.seed(2000 * p + rep)
  n <- 800
  x <- matrix(runif(n * p), n, p, dimnames = list(NULL, paste0("x", seq_len(p))))
  # Confounded assignment, and a constant treatment effect so the truth is one
  # number the interval either covers or does not.
  z <- rbinom(n, 1, plogis(2 * (x[, 1] - 0.5)))
  d <- data.frame(y = friedman(x) / 5 + TAU * z + rnorm(n), z = z, x)

  t(vapply(LEVELS, function(lv) {
    set.seed(7)
    t0 <- Sys.time()
    fit <- bartisan(y ~ ., data = d, family = gaussian(), chains = 2,
                    control = bartisan_control(sparsity = lv, verbose = FALSE))
    done <<- done + 1L
    prog_tick(pr, i = done,
              secs = as.numeric(difftime(Sys.time(), t0, units = "secs")),
              label = sprintf("contrast p=%d rep %d / %s", p, rep, lv))
    a <- avg_comparisons(fit, variables = "z")
    post <- attr(a, "posterior_draws")
    vi <- variable_importance(fit)

    c(est = a$estimate, lo = a$conf.low, hi = a$conf.high,
      width = a$conf.high - a$conf.low,
      covers = as.numeric(a$conf.low <= TAU && a$conf.high >= TAU),
      atom = if (is.null(post)) NA_real_ else mean(post == 0),
      prop_used = vi[["prop_used"]][vi[["variable"]] == "z"])
  }, numeric(7)))
}

for (p in P_GRID) {
  key <- sprintf("contrast, p = %d", p)
  each <- list()

  for (r in seq_len(REPS)) {
    each[[r]] <- contrast_rep(p, r)
    results[[key]] <- each
    checkpoint(FALSE)
  }

  acc <- Reduce(`+`, each) / REPS
  cat(sprintf("\n== A contrast, p = %d, truth %.2f, mean of %d reps ==\n",
              p, TAU, REPS))
  cat(sprintf("%-10s %7s %7s %7s %7s %7s %9s\n",
              "sparsity", "est", "lo", "hi", "width", "atom", "prop_used"))
  for (i in seq_along(LEVELS)) {
    cat(sprintf("%-10s %7.3f %7.3f %7.3f %7.3f %7.3f %9.3f\n",
                LEVELS[i], acc[i, "est"], acc[i, "lo"], acc[i, "hi"],
                acc[i, "width"], acc[i, "atom"], acc[i, "prop_used"]))
  }
  cat(sprintf("coverage of the truth: %s\n",
              paste(sprintf("%s %.2f", LEVELS, acc[, "covers"]), collapse = "  ")))
}

checkpoint(TRUE)
on.exit()
prog_end(pr, "done", sprintf("%d fits", done))
cat("\nwrote", OUT, "\n")
