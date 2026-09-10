# What the data augmentation buys, family by family.
#
# This is the benchmark behind the "What Data Augmentation Buys" table in
# `vignette("implementation")`. Each row is one family and gate, fitted twice on
# the same data -- `augment = TRUE` and `augment = FALSE` -- and the table
# reports three ratios of the second to the first:
#
#   Speed    seconds(FALSE) / seconds(TRUE)          higher = augmentation faster
#   ESS      ESS(TRUE) / ESS(FALSE)                  lower  = augmentation mixes worse
#   ESS/sec  their product, which is what to judge on
#
# The augmentation is a reparameterization, not an approximation, so both arms
# target the same posterior and the only question is how efficiently. ESS is
# taken on the worst-mixing quantity, dropping the averages over observations,
# for the reason the other scripts here give: an average mixes about as well as
# the draw count whatever the sampler is doing.
#
# Paired: both arms see the same data in each replicate.
A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))
Sys.setenv(PROGRESS_ROOT = "/Users/NoahGreifer/Dropbox/Research/R/bartisan/.progress-runs")

suppressMessages(library(bartisan))

REPS <- 15L
N <- 400L
P <- 8L
OUT <- file.path("_dev", "augment-benchmark.rds")

make <- function(seed) {
  set.seed(seed)
  X <- matrix(runif(N * P), N, P)
  colnames(X) <- paste0("x", seq_len(P))
  d <- as.data.frame(X)
  mu <- 2 * sin(pi * d$x1) + 1.5 * d$x2
  lin <- mu - mean(mu)

  d$y_bin <- rbinom(N, 1L, stats::plogis(lin))
  d$y_ord <- factor(cut(mu + rnorm(N, 0, 0.5), 3L), ordered = TRUE)
  d$y_cnt <- rpois(N, exp(0.4 + 0.5 * mu))
  d$y_od <- rnbinom(N, size = 2, mu = exp(0.4 + 0.5 * mu))
  d$y_zi <- ifelse(runif(N) < 0.3, 0L, rpois(N, exp(0.4 + 0.5 * mu)))
  d$y_zo <- ifelse(runif(N) < 0.3, 0L, rnbinom(N, size = 2, mu = exp(0.4 + 0.5 * mu)))
  d$y_mn <- factor(sample(3L, N, TRUE, prob = c(0.4, 0.35, 0.25)))
  tm <- rexp(N, exp(-0.5 - mu / 2))
  ev <- rbinom(N, 1L, 0.7)
  d$y_surv <- survival::Surv(tm, ev)
  d
}

rhs <- paste(paste0("x", seq_len(P)), collapse = " + ")
form <- function(y) stats::as.formula(paste(y, "~", rhs))

# One row per (family, gate) the vignette's table reports.
cells <- list(
  list(row = "lognormal_aft(), hard rules",   y = "y_surv", gate = "hard",       fam = quote(lognormal_aft())),
  list(row = "lognormal_aft(), soft rules",   y = "y_surv", gate = "smoothstep", fam = quote(lognormal_aft())),
  list(row = "loglogistic_aft(), hard rules", y = "y_surv", gate = "hard",       fam = quote(loglogistic_aft())),
  list(row = "loglogistic_aft(), soft rules", y = "y_surv", gate = "smoothstep", fam = quote(loglogistic_aft())),
  list(row = "ordinal(\"probit\"), hard rules", y = "y_ord", gate = "hard",       fam = quote(ordinal("probit"))),
  list(row = "ordinal(\"probit\"), soft rules", y = "y_ord", gate = "smoothstep", fam = quote(ordinal("probit"))),
  list(row = "ordinal(\"logit\"), hard rules",  y = "y_ord", gate = "hard",       fam = quote(ordinal("logit"))),
  list(row = "ordinal(\"logit\"), soft rules",  y = "y_ord", gate = "smoothstep", fam = quote(ordinal("logit"))),
  list(row = "binomial(\"probit\")",            y = "y_bin", gate = "smoothstep", fam = quote(stats::binomial("probit"))),
  list(row = "binomial(\"logit\")",             y = "y_bin", gate = "smoothstep", fam = quote(stats::binomial("logit"))),
  list(row = "multinomial(), hard rules",       y = "y_mn",  gate = "hard",       fam = quote(multinomial())),
  list(row = "multinomial(), soft rules",       y = "y_mn",  gate = "smoothstep", fam = quote(multinomial())),
  list(row = "zi_poisson(), hard rules",        y = "y_zi",  gate = "hard",       fam = quote(zi_poisson())),
  list(row = "zi_poisson(), soft rules",        y = "y_zi",  gate = "smoothstep", fam = quote(zi_poisson())),
  list(row = "zi_negbin(), hard rules",         y = "y_zo",  gate = "hard",       fam = quote(zi_negbin())),
  list(row = "negbin()",                        y = "y_od",  gate = "hard",       fam = quote(negbin()))
)

quiet <- function(expr) {
  invisible(utils::capture.output(out <- suppressMessages(suppressWarnings(expr))))
  out
}

# Both readings, because they have very different variances. The worst-quantity
# ESS is the conservative number and the one a user is bounded by, but it is a
# minimum over many quantities and so is a high-variance statistic: at three
# replicates its ratio between the augmented and unaugmented arms swung by a
# median factor of 4.5 within a single cell, which is why this script now runs
# fifteen. The median quantity is far stabler and says whether the whole chain
# moved or only its worst corner.
ess_pair <- function(fit) {
  tab <- diagnose(fit)[["table"]]
  tab <- tab[!grepl("average over observations", tab[["quantity"]], fixed = TRUE), ]
  c(worst = min(tab[["ess_bulk"]], na.rm = TRUE),
    median = stats::median(tab[["ess_bulk"]], na.rm = TRUE))
}

# Warm-up, so the first cell does not absorb process start-up.
invisible(quiet(bartisan(y_bin ~ x1 + x2, make(1L), family = stats::binomial(),
                         control = bartisan_control(num_trees = 5L,
                                                    num_burn = 20L,
                                                    num_draws = 20L))))

pr <- prog_init(total = length(cells) * REPS * 2L,
                title = "What data augmentation buys", unit = "fit",
                kind = "benchmark")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
k <- 0L

for (i in seq_along(cells)) {
  cl <- cells[[i]]
  fam <- eval(cl$fam)

  for (r in seq_len(REPS)) {
    d <- make(600L + r)

    for (aug in c(TRUE, FALSE)) {
      ctrl <- bartisan_control(num_trees = 50L, num_burn = 500L,
                               num_draws = 1000L, chains = 2L,
                               gate = cl$gate, augment = aug)
      t0 <- Sys.time()
      fit <- tryCatch(quiet(bartisan(form(cl$y), d, family = fam,
                                     control = ctrl)),
                      error = function(e) NULL)
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      ess <- if (is.null(fit)) c(worst = NA_real_, median = NA_real_) else quiet(ess_pair(fit))

      k <- k + 1L
      rows[[k]] <- data.frame(row = cl$row, gate = cl$gate, augment = aug,
                              rep = r, secs = secs,
                              ess_worst = ess[["worst"]],
                              ess_median = ess[["median"]],
                              ok = !is.null(fit), stringsAsFactors = FALSE)
      prog_tick(pr, i = k, secs = secs, ok = !is.null(fit),
                label = sprintf("%s augment=%s r%d", cl$row, aug, r))
    }
  }

  saveRDS(list(res = do.call(rbind, rows), reps = REPS, n = N, p = P,
               complete = FALSE, done = i, total = length(cells)), OUT)
}

res <- do.call(rbind, rows)
saveRDS(list(res = res, reps = REPS, n = N, p = P, complete = TRUE,
             done = length(cells), total = length(cells)), OUT)
on.exit()
prog_end(pr, "done")

# Ratios of means rather than means of ratios: a per-replicate ratio of two
# noisy ESS estimates has a heavy right tail, and averaging those tails is what
# made the three-replicate reading unusable.
out <- do.call(rbind, lapply(split(res, res$row), function(z) {
  a <- z[z$augment, ]
  b <- z[!z$augment, ]
  rat <- function(col) mean(a[[col]], na.rm = TRUE) / mean(b[[col]], na.rm = TRUE)
  # Bootstrap over replicates, to say how far the ratio can be trusted.
  boot <- replicate(2000L, {
    i <- sample(nrow(a), replace = TRUE)
    mean(a[["ess_worst"]][i], na.rm = TRUE) /
      mean(b[["ess_worst"]][i], na.rm = TRUE)
  })
  data.frame(row = z$row[1],
             speed = mean(b$secs) / mean(a$secs),
             ess = rat("ess_worst"),
             ess_lo = stats::quantile(boot, 0.1, names = FALSE),
             ess_hi = stats::quantile(boot, 0.9, names = FALSE),
             ess_med_q = rat("ess_median"),
             secs_aug = mean(a$secs), secs_raw = mean(b$secs),
             stringsAsFactors = FALSE)
}))
out$ess_per_sec <- out$speed * out$ess
print(out[order(-out$ess_per_sec), ], row.names = FALSE, digits = 3)
cat("\n`ess` is the worst-quantity ESS ratio with an 80% bootstrap interval",
    "over replicates;\n`ess_med_q` is the same ratio for the median quantity.",
    "\n", REPS, "replicates.\n")
cat("\nwrote", OUT, "\n")
