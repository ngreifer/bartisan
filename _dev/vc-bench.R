# Timing and bit-exactness for the varying-coefficient family, across the
# configurations that route through `VaryingCoefficientFamily`. The sibling of
# `_dev/review-bench.R`, which has no `vc()` configuration.
#
# Usage: Rscript _dev/vc-bench.R <lib.loc> <out.rds> <label>, once per build,
# each build installed into its own library; then merge on config and rep and
# require the `hash` columns to agree before reading any timing.
args <- commandArgs(TRUE); lib <- args[1L]; out <- args[2L]; label <- args[3L]
suppressPackageStartupMessages(library(bartisan, lib.loc = lib))
stopifnot(bartisan:::.bartisan_optimized())

sim <- function(n, seed) {
  set.seed(seed)
  d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n), x4 = runif(n),
                  f1 = factor(sample(letters[1:4], n, TRUE)),
                  z = rbinom(n, 1L, 0.5), zc = rnorm(n))
  d$mu <- 2 * sin(pi * d$x1 * d$x2) + 0.5 * as.integer(d$f1)
  d
}
cfg <- list(
  vc_gaussian   = list(f = y ~ x1+x2+x3+x4+f1 + vc(z), fam = gaussian(),
                       y = function(d) d$mu + (1 + d$x1) * d$z + rnorm(nrow(d))),
  vc_continuous = list(f = y ~ x1+x2+x3+x4+f1 + vc(zc), fam = gaussian(),
                       y = function(d) d$mu + (1 + d$x1) * d$zc + rnorm(nrow(d))),
  vc_binomial   = list(f = y ~ x1+x2+x3+x4+f1 + vc(z), fam = binomial(),
                       y = function(d) rbinom(nrow(d), 1, plogis(d$mu - 1 + d$z))),
  vc_poisson    = list(f = y ~ x1+x2+x3+x4+f1 + vc(z), fam = poisson(),
                       y = function(d) rpois(nrow(d), exp(d$mu/2 + 0.3*d$z))),
  vc_lsd        = list(f = y ~ x1+x2+x3+x4+f1 + vc(z), fam = gaussian_ls(),
                       y = function(d) d$mu + d$z + rnorm(nrow(d))),
  bcf_default   = list(bcf = TRUE, fam = gaussian(),
                       y = function(d) d$mu + (1 + d$x1) * d$z + rnorm(nrow(d))))

rows <- list()
for (nm in names(cfg)) {
  cf <- cfg[[nm]]
  for (r in 1:3) {
    d <- sim(1000L, 100L + r); set.seed(1000L + r); d$y <- cf$y(d)
    ctl <- bartisan_control(num_trees = 50L, num_burn = 300L,
                            num_draws = 300L, verbose = FALSE)
    set.seed(2000L + r); t0 <- proc.time()[[3]]
    fit <- suppressMessages(suppressWarnings(
      if (isTRUE(cf$bcf))
        bcf(y ~ x1+x2+x3+x4+f1, treat = ~ z, data = d, family = cf$fam,
            propensity = FALSE, control = ctl)
      else bartisan(cf$f, data = d, family = cf$fam, control = ctl)))
    secs <- proc.time()[[3]] - t0
    rows[[length(rows)+1L]] <- data.frame(config = nm, rep = r, seconds = secs,
      hash = rlang::hash(list(lapply(fit$eta, unname), fit$aux, fit$sigma_mu,
                              fit$bandwidth, fit$loglik)),
      stringsAsFactors = FALSE)
    cat(sprintf("%-14s rep %d  %.2fs\n", nm, r, secs)); utils::flush.console()
  }
}
saveRDS(list(res = do.call(rbind, rows), label = label), out)
