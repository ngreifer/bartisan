# Do 95% credible intervals for the additive predictor cover it 95% of the time?
#
# This is the benchmark behind the coverage sentence in the Correctness section
# of `vignette("implementation")`, which currently reports .95 (Gaussian), .91
# (binomial), .96 (Poisson) and .96 (gamma), and attributes the binomial
# shortfall to a bias-to-posterior-SD ratio near .8.
#
# Two things are measured per family, over replicate datasets at the package
# defaults:
#
#   coverage   the share of observations whose true eta lies in its own 95%
#              interval, averaged over replicates
#   bias/sd    mean |posterior mean - truth| over mean posterior SD, which is
#              the quantity the vignette's explanation rests on
#
# The claim is about the predictor, not the response, so everything is on the
# link scale and the truth is the linear predictor the data were generated from.
# Both are pointwise: this is not a statement about an interval for an average.
#
# The vignette also says running more or longer chains does not help, which
# rules out mixing as the cause. That is a checkable claim rather than an aside,
# so a second arm refits the binomial at four chains and five times the draws.
A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))
Sys.setenv(PROGRESS_ROOT = "/Users/NoahGreifer/Dropbox/Research/R/bartisan/.progress-runs")

suppressMessages(library(bartisan))

REPS <- 40L
N <- 500L
P <- 10L
OUT <- file.path("_dev", "coverage-calibration.rds")

make <- function(seed) {
  set.seed(seed)
  X <- matrix(runif(N * P), N, P)
  colnames(X) <- paste0("x", seq_len(P))
  d <- as.data.frame(X)
  # One mean function for every family, centered, then mapped to each scale.
  eta <- 2 * sin(pi * d$x1) + 1.5 * d$x2 + 1.5 * (d$x3 - 0.5)^2
  eta <- eta - mean(eta)

  d$eta <- eta
  d$y_num <- eta + rnorm(N, 0, 0.7)
  d$y_bin <- rbinom(N, 1L, stats::plogis(eta))
  d$y_cnt <- rpois(N, exp(eta + 1))
  d$y_pos <- rgamma(N, shape = 3, rate = 3 / exp(eta + 1))
  d
}

rhs <- paste(paste0("x", seq_len(P)), collapse = " + ")
form <- function(y) stats::as.formula(paste(y, "~", rhs))

# The truth on each family's own link scale, matching how the response was made.
truth_of <- function(d, fam) {
  switch(fam,
         gaussian = d$eta,
         binomial = d$eta,
         poisson = d$eta + 1,
         gamma = d$eta + 1)
}

cells <- list(
  list(fam = "gaussian", y = "y_num", f = quote(stats::gaussian()),      long = FALSE),
  list(fam = "binomial", y = "y_bin", f = quote(stats::binomial()),      long = FALSE),
  list(fam = "poisson",  y = "y_cnt", f = quote(stats::poisson()),       long = FALSE),
  list(fam = "gamma",    y = "y_pos", f = quote(stats::Gamma("log")),    long = FALSE),
  # The mixing control: same family, four chains and five times the draws.
  list(fam = "binomial (long run)", y = "y_bin", f = quote(stats::binomial()), long = TRUE)
)

quiet <- function(expr) {
  invisible(utils::capture.output(out <- suppressMessages(suppressWarnings(expr))))
  out
}

pr <- prog_init(total = length(cells) * REPS, title = "Interval coverage for eta",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
k <- 0L

for (i in seq_along(cells)) {
  cl <- cells[[i]]
  fam <- eval(cl$f)
  base <- sub(" .*", "", cl$fam)

  ctrl <- if (cl$long) {
    bartisan_control(num_trees = 50L, chains = 4L, num_burn = 1000L,
                     num_draws = 4000L)
  } else {
    bartisan_control(num_trees = 50L)          # package defaults otherwise
  }

  for (r in seq_len(REPS)) {
    d <- make(30000L + 71L * r)
    tru <- truth_of(d, base)

    fit <- quiet(bartisan(form(cl$y), d, family = fam, control = ctrl))
    dr <- stats::predict(fit, type = "link", draws = TRUE)
    if (is.list(dr)) dr <- dr[[1L]]

    post <- colMeans(dr)
    sdv <- apply(dr, 2L, stats::sd)
    lo <- apply(dr, 2L, stats::quantile, 0.025)
    hi <- apply(dr, 2L, stats::quantile, 0.975)

    k <- k + 1L
    rows[[k]] <- data.frame(
      family = cl$fam, rep = r,
      coverage = mean(lo <= tru & tru <= hi),
      bias_sd = mean(abs(post - tru)) / mean(sdv),
      mean_sd = mean(sdv), mean_abs_bias = mean(abs(post - tru)),
      stringsAsFactors = FALSE)

    prog_tick(pr, i = k, label = sprintf("%s rep %d", cl$fam, r))
  }

  saveRDS(list(res = do.call(rbind, rows), reps = REPS, n = N, p = P,
               complete = FALSE, done = i, total = length(cells)), OUT)
}

res <- do.call(rbind, rows)
saveRDS(list(res = res, reps = REPS, n = N, p = P, complete = TRUE,
             done = length(cells), total = length(cells)), OUT)
on.exit()
prog_end(pr, "done")

out <- do.call(rbind, lapply(split(res, res$family), function(z) {
  data.frame(family = z$family[1],
             coverage = mean(z$coverage),
             mcse = stats::sd(z$coverage) / sqrt(nrow(z)),
             bias_sd = mean(z$bias_sd),
             mean_sd = mean(z$mean_sd),
             stringsAsFactors = FALSE)
}))
print(out, row.names = FALSE, digits = 3)
cat("\n", REPS, " replicates, n = ", N, ", p = ", P,
    "; non-long arms at package defaults (200 warmup + 800 draws, 1 chain)\n",
    sep = "")
cat("wrote", OUT, "\n")
