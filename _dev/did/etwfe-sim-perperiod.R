# Does one forest per period remove the bias the additive time forest leaves
# where linearity fails?
#
# What is being tested. In `_dev/did/etwfe-sim.R`'s nonlinear cell the
# additive forest ETWFE, `y ~ cohort + x1 + x2 + vc(one, ~ period + x1 + x2) +
# (1 | id)`, cut linear ETWFE's bias in the overall ATT from +0.233 to +0.037
# but did not remove it, and what remained grew with exposure (+0.010, +0.055,
# +0.096 at e = 0, 1, 2) and was not in the intervals (coverage 55% for the
# conditional-mean imputation, 75% for the predictive one). The suspected
# mechanism: the trend is m(x) (t - 1) / 5, an interaction of period and the
# covariates, which a single forest over (period, x) has to build from splits
# on both and which its regularization attenuates most in the late periods,
# where only the never-treated inform it. The per-period form,
# `vc(periodf, ~ x1 + x2)`, gives every period its own forest over x (the
# analogue of Wooldridge's fs_t and fs_t * x terms), so the interaction needs no
# period splits. The claim is that it removes most of the remaining bias. It
# could be false if the bias comes instead from thin support: cohort 4 is
# selected on the extremes of x1, where m is largest and the never-treated are
# fewest, and a forest shrinks toward the mean there whatever its structure.
#
# The measurement. The nonlinear cell's 20 datasets, regenerated from the same
# seeds by the same code (checked against the recorded truths), so the
# comparison with the additive form is paired. Per-period + random intercept,
# 4 chains, hard gates, 500 + 1000 draws, 100 trees on the cohort forest and 50
# on each period forest (the additive run had 100 and 50). Both imputations.
# Read the bias of the overall ATT and at e = 2, paired against the additive
# form's +0.037 and +0.096 on the same datasets; bias is the decisive endpoint
# because at 20 replicates coverage is known only to about 0.1. Then coverage
# and width.
#
# What each outcome means. Bias at e = 2 down to a third of the additive form's
# or less, with coverage near 95%: the per-period form is the time structure to
# recommend, and the additive form's bias was the attenuated interaction. Bias
# about where it was: the problem is the covariate support, not the forest's
# structure, and the notes say forest ETWFE reduces but does not remove the
# bias from a nonlinear trend, with intervals that do not cover it where the
# treated cohorts sit at the edges of the controls' covariate distribution.

source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")
suppressMessages(library(bartisan))

rep_from <- as.integer(Sys.getenv("REP_FROM", "1"))
rep_to <- as.integer(Sys.getenv("REP_TO", "20"))
smoke <- nzchar(Sys.getenv("SMOKE"))
out_file <- Sys.getenv("OUT", "_dev/did/etwfe-sim-perperiod.rds")
future::plan(future::multisession, workers = 4L)

# ---- Verbatim from etwfe-sim.R, so the datasets are the same -------------------
n <- 500L
periods <- 1:6
cohorts <- c(4L, 5L, 6L)

m_fun <- list(
  nonlinear = function(x1, x2) 6 * (x1 - 0.5)^2 + 2 * sin(2 * pi * x2),
  linear = function(x1, x2) 1.5 * x1 + x2)

simulate <- function(cell, seed) {
  set.seed(seed)
  x1 <- stats::runif(n)
  x2 <- stats::runif(n)
  score <- cbind(0, 12 * (x1 - 0.5)^2 - 1.5, 1.5 * sin(2 * pi * x2) - 0.5,
                 x1 - 0.5)
  prob <- exp(score) / rowSums(exp(score))
  g <- c(Inf, cohorts)[apply(prob, 1L, function(p) sample.int(4L, 1L, prob = p))]
  level <- c(0, 0.8, -0.5, 0.3)[match(g, c(Inf, cohorts))]
  alpha <- stats::rnorm(n)
  m <- m_fun[[cell]](x1, x2)

  d <- expand.grid(id = seq_len(n), period = periods)
  i <- d$id
  d$x1 <- x1[i]
  d$x2 <- x2[i]
  d$g <- g[i]
  d$w <- as.integer(is.finite(d$g) & d$period >= d$g)
  d$tau <- ifelse(d$w == 1L, 1 + 0.25 * (d$period - d$g) + 0.5 * (d$x2 - 0.5), 0)
  d$y <- alpha[i] + level[i] + 0.5 * d$x2 + 0.3 * d$period +
    m[i] * (d$period - 1) / 5 + d$tau + stats::rnorm(nrow(d), sd = 0.3)
  d$cohort <- factor(ifelse(is.finite(d$g), d$g, 0), levels = c(0, cohorts))
  d$periodf <- factor(d$period)
  d$one <- 1
  d$cell <- ifelse(d$w == 1L, paste(d$g, d$period, sep = ":"), "none")
  d$e <- ifelse(d$w == 1L, d$period - d$g, NA)
  d
}

# Weights over the treated cells for each reported quantity.
aggregations <- function(tr) {
  cells <- sort(unique(tr$cell))
  n_gt <- as.numeric(table(tr$cell)[cells])
  e_of <- tapply(tr$e, tr$cell, `[`, 1L)[cells]
  out <- list(overall = n_gt)
  for (e in 0:2) {
    out[[paste0("e", e)]] <- ifelse(e_of == e, n_gt, 0)
  }
  lapply(out, function(w) w / sum(w))
}

truth <- function(d) {
  tr <- d[d$w == 1L, ]
  cells <- sort(unique(tr$cell))
  cell_truth <- tapply(tr$tau, tr$cell, mean)[cells]
  vapply(aggregations(tr), function(w) sum(w * cell_truth), numeric(1L))
}
# ---- End of the verbatim block --------------------------------------------------

bart_per_period <- function(d) {
  ctl <- d[d$w == 0L, ]
  tr <- d[d$w == 1L, ]
  fit <- bartisan(y ~ cohort + x1 + x2 + vc(periodf, ~ x1 + x2) + (1 | id),
                  data = ctl, family = gaussian(), chains = 4L, gate = "hard",
                  num_trees = c(100L, rep(50L, length(periods))),
                  num_burn = if (smoke) 20L else 500L,
                  num_draws = if (smoke) 20L else 1000L, verbose = FALSE)
  cells <- sort(unique(tr$cell))
  w <- aggregations(tr)
  summarize <- function(y0) {
    te <- sweep(-y0, 2L, tr$y, "+")
    cell_draws <- sapply(cells, function(k) rowMeans(te[, tr$cell == k, drop = FALSE]))
    lapply(w, function(wt) {
      v <- drop(cell_draws %*% wt)
      c(estimate = mean(v), lower = unname(stats::quantile(v, .025)),
        upper = unname(stats::quantile(v, .975)))
    })
  }
  list(mean = summarize(predict(fit, newdata = tr, draws = TRUE)),
       predictive = summarize(rstantools::posterior_predict(fit, newdata = tr)))
}

reps <- seq(rep_from, rep_to)
pr <- prog_init(total = length(reps),
                title = sprintf("Reps %d-%d, ETWFE per-period forests", rep_from, rep_to),
                unit = "replicate", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)
rows <- list()

for (r in reps) {
  # etwfe-sim.R seeds cell k of replicate r with 1000 r + k, and "nonlinear" is
  # its first cell.
  d <- simulate("nonlinear", seed = 1000L * r + 1L)
  tv <- truth(d)
  t0 <- Sys.time()
  est <- setNames(bart_per_period(d), c("per_period_mean", "per_period_predictive"))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  for (method in names(est)) {
    for (q in names(tv)) {
      e <- est[[method]][[q]]
      rows[[length(rows) + 1L]] <- data.frame(
        rep = r, cell = "nonlinear", method = method, quantity = q,
        truth = tv[[q]], estimate = e[["estimate"]], lower = e[["lower"]],
        upper = e[["upper"]], seconds = round(secs, 1))
    }
  }
  saveRDS(list(res = do.call(rbind, rows), complete = FALSE,
               done = length(rows) / 8L, total = length(reps)), out_file)
  prog_tick(pr, label = sprintf("rep %d", r))
}

res <- do.call(rbind, rows)
saveRDS(list(res = res, complete = TRUE, done = length(reps),
             total = length(reps)), out_file)
on.exit()
prog_end(pr, "done")
