# Timing and bit-exactness benchmark for the C++ changes.
#
# Usage: Rscript _dev/review-bench.R <lib.loc> <out.rds> <label> [config]
#
# Install each build into its own library first (R CMD INSTALL --library=...)
# and pass that library; the two builds are never loaded into one session.
#
# What is tested: that the sampler changes reduce wall time per fit and leave the
# draws bit-identical, on a fixed set of (family, gate) configurations. Outcome:
# seconds per fit (three replicates, distinct seeds, identical across builds) and
# a hash of the first predictor's stored draws plus the nuisance parameters.
# Reading: compare `seconds` between builds per configuration (median), and
# `hash_eta`/`hash_aux` must be identical between builds for every row.
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1L]; out <- args[2L]; label <- args[3L]
# An optional fourth argument names one configuration, for a driver that
# interleaves builds per configuration; progress is then the driver's job.
only <- if (length(args) >= 4L) args[4L] else NA_character_

source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")
suppressPackageStartupMessages(library(bartisan, lib.loc = lib))
stopifnot(bartisan:::.bartisan_optimized())
cat("bartisan loaded from:", dirname(system.file(package = "bartisan")), "\n")

sim <- function(n = 1000L, seed = 1L) {
  set.seed(seed)
  d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n), x4 = runif(n),
                  x5 = runif(n), x6 = runif(n),
                  f1 = factor(sample(letters[1:4], n, TRUE)),
                  f2 = factor(sample(c("u", "v", "w"), n, TRUE)))
  d$mu <- 2 * sin(pi * d$x1 * d$x2) + 4 * (d$x3 - .5)^2 + d$x4 +
    0.3 * as.integer(d$f1)
  d
}

configs <- list(
  gaussian_soft   = list(family = gaussian(), gate = "smoothstep", augment = TRUE,
                         y = function(d) d$mu + rnorm(nrow(d))),
  gaussian_hard   = list(family = gaussian(), gate = "hard", augment = TRUE,
                         y = function(d) d$mu + rnorm(nrow(d))),
  binom_aug_soft  = list(family = binomial(), gate = "smoothstep", augment = TRUE,
                         y = function(d) rbinom(nrow(d), 1, plogis(d$mu - 2))),
  binom_noaug_soft = list(family = binomial(), gate = "smoothstep", augment = FALSE,
                         y = function(d) rbinom(nrow(d), 1, plogis(d$mu - 2))),
  negbin_soft     = list(family = negbin(), gate = "smoothstep", augment = TRUE,
                         y = function(d) rnbinom(nrow(d), mu = exp(d$mu / 2), size = 2)),
  poisson_soft    = list(family = poisson(), gate = "smoothstep", augment = TRUE,
                         y = function(d) rpois(nrow(d), exp(d$mu / 2))),
  gamma_soft      = list(family = Gamma("log"), gate = "smoothstep", augment = TRUE,
                         y = function(d) rgamma(nrow(d), shape = 2, rate = 2 / exp(d$mu / 2))),
  beta_soft       = list(family = Beta(), gate = "smoothstep", augment = TRUE,
                         y = function(d) { m <- plogis(d$mu - 2); rbeta(nrow(d), m * 10, (1 - m) * 10) }),
  ordinal_noaug_soft = list(family = ordinal(), gate = "smoothstep", augment = FALSE,
                         y = function(d) cut(d$mu + rlogis(nrow(d)), c(-Inf, 1, 2.5, 4, Inf), ordered_result = TRUE)),
  multinom_noaug_soft = list(family = multinomial(), gate = "smoothstep", augment = FALSE,
                         y = function(d) factor(apply(cbind(0, d$mu - 2, 0.5 * d$mu - 1) + matrix(rlogis(3 * nrow(d)), ncol = 3), 1, which.max))),
  dpm_soft        = list(family = dpm(), gate = "smoothstep", augment = TRUE,
                         y = function(d) d$mu + rt(nrow(d), 3))
)

if (!is.na(only)) configs <- configs[only]
reps <- 3L
rows <- list()
pr <- if (is.na(only)) prog_init(total = length(configs) * reps, title = paste("bench", label),
                unit = "fit", kind = "benchmark",
                root = "/Users/NoahGreifer/Dropbox/Research/R/bartisan/.progress-runs")
if (!is.null(pr)) on.exit(prog_end(pr, "failed"), add = TRUE)

for (nm in names(configs)) {
  cf <- configs[[nm]]
  for (r in seq_len(reps)) {
    d <- sim(seed = 100L + r)
    set.seed(1000L + r)
    d$y <- cf$y(d)
    ctrl <- bartisan_control(num_trees = 50L, num_burn = 300L, num_draws = 300L,
                             gate = cf$gate, augment = cf$augment, verbose = FALSE)
    set.seed(2000L + r)
    t0 <- proc.time()[["elapsed"]]
    fit <- suppressMessages(suppressWarnings(
      bartisan(y ~ x1 + x2 + x3 + x4 + x5 + x6 + f1 + f2, data = d,
               family = cf$family, control = ctrl)))
    secs <- proc.time()[["elapsed"]] - t0
    rows[[length(rows) + 1L]] <- data.frame(
      config = nm, rep = r, seconds = secs,
      hash_eta = rlang::hash(lapply(fit$eta, unname)),
      hash_aux = rlang::hash(list(fit$aux, fit$sigma_mu, fit$bandwidth, fit$loglik)),
      stringsAsFactors = FALSE)
    if (!is.null(pr)) prog_tick(pr, label = sprintf("%s rep %d (%.1fs)", nm, r, secs))
    saveRDS(list(res = do.call(rbind, rows), complete = FALSE, label = label), out)
  }
}
on.exit()
if (!is.null(pr)) prog_end(pr, "done")
saveRDS(list(res = do.call(rbind, rows), complete = TRUE, label = label), out)
print(do.call(rbind, rows))
