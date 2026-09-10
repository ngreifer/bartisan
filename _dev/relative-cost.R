# Cost of each family relative to a Gaussian fit of the same size.
#
# This is the benchmark behind the "Relative cost" table in
# `vignette("implementation")`. That table and the prose beneath it currently
# disagree: the table says `Gamma("log")` and `tweedie()` cost about 2x a
# Gaussian, while the paragraph below says tweedie costs 5.5x where
# `Gamma("log")` costs 5.2x. Both cannot be right, so both are re-measured.
# `dpm()` was missing from the table and is included here.
#
# Everything is held equal except the family: same n, same predictors, same
# tree count, same sweeps, same gate. Only the response changes, because each
# family needs one it is defined on. The reported figure is the MINIMUM over
# replicates rather than the mean: for a timing the minimum is the estimate
# least contaminated by other load, and the vignette already warns that
# absolute times move by about 20% run to run.
A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))
Sys.setenv(PROGRESS_ROOT = "/Users/NoahGreifer/Dropbox/Research/R/bartisan/.progress-runs")

suppressMessages(library(bartisan))

REPS <- 5L
N <- 500L
P <- 10L
OUT <- file.path("_dev", "relative-cost.rds")

set.seed(20260910)
X <- matrix(runif(N * P), N, P)
colnames(X) <- paste0("x", seq_len(P))
d <- as.data.frame(X)
mu <- 2 * sin(pi * d$x1) + 1.5 * d$x2

# One response per family, on the scale that family is defined on.
d$y_num <- mu + rnorm(N, 0, 0.7)
d$y_t <- mu + rt(N, df = 3) / sqrt(3)                       # dpm has a reason to exist
d$y_pos <- rgamma(N, shape = 2, rate = 2 / exp(mu / 2))
d$y_cnt <- rpois(N, exp(0.4 + 0.5 * mu))
d$y_od <- rnbinom(N, size = 2, mu = exp(0.4 + 0.5 * mu))
d$y_bin <- rbinom(N, 1L, plogis(mu - mean(mu)))
d$y_ord <- factor(cut(mu + rnorm(N, 0, 0.5), 3L), ordered = TRUE)
d$y_zi <- ifelse(runif(N) < 0.3, 0L, rpois(N, exp(0.4 + 0.5 * mu)))
d$y_01 <- pmin(pmax(plogis(mu - mean(mu) + rnorm(N, 0, 0.4)), 0), 1)
d$y_01[d$y_01 < 0.15] <- 0
d$y_01[d$y_01 > 0.9] <- 1
d$y_beta <- plogis(mu - mean(mu) + rnorm(N, 0, 0.4))
d$y_tw <- ifelse(runif(N) < 0.25, 0, rgamma(N, shape = 2, rate = 2 / exp(mu)))

rhs <- paste(paste0("x", seq_len(P)), collapse = " + ")
form <- function(y) stats::as.formula(paste(y, "~", rhs))

cells <- list(
  list(label = "gaussian()",             y = "y_num",  fam = quote(stats::gaussian())),
  list(label = "dpm()",                  y = "y_t",    fam = quote(dpm())),
  list(label = "gaussian_ls()",          y = "y_num",  fam = quote(gaussian_ls())),
  list(label = "Gamma(\"log\")",         y = "y_pos",  fam = quote(stats::Gamma("log"))),
  list(label = "poisson()",              y = "y_cnt",  fam = quote(stats::poisson())),
  list(label = "negbin()",               y = "y_od",   fam = quote(negbin())),
  list(label = "binomial(\"logit\")",    y = "y_bin",  fam = quote(stats::binomial())),
  list(label = "binomial(\"probit\")",   y = "y_bin",  fam = quote(stats::binomial("probit"))),
  list(label = "ordinal(\"probit\")",    y = "y_ord",  fam = quote(ordinal("probit"))),
  list(label = "ordinal(\"logit\")",     y = "y_ord",  fam = quote(ordinal("logit"))),
  list(label = "zi_poisson()",           y = "y_zi",   fam = quote(zi_poisson())),
  list(label = "Beta()",                 y = "y_beta", fam = quote(Beta())),
  list(label = "ordbeta()",              y = "y_01",   fam = quote(ordbeta())),
  list(label = "tweedie()",              y = "y_tw",   fam = quote(tweedie())),
  list(label = "tweedie(power drawn)",   y = "y_tw",
       fam = quote(tweedie(power = NULL)))
)

ctrl <- bartisan_control(num_trees = 50L, num_burn = 200L, num_draws = 400L)

# One throwaway fit before any timing. The first fit in a fresh process pays
# for loading and for first-touch allocation, and whichever family is measured
# first would otherwise absorb it -- which is `gaussian()`, the baseline every
# other number is divided by. A dry run at tiny settings showed the baseline at
# 0.30s against 0.01s for every other family, entirely from this.
invisible(suppressMessages(suppressWarnings(
  bartisan(form("y_num"), d, family = stats::gaussian(),
           control = bartisan_control(num_trees = 5L, num_burn = 20L,
                                      num_draws = 20L)))))

pr <- prog_init(total = length(cells) * REPS, title = "Cost relative to gaussian",
                unit = "fit", kind = "benchmark")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
k <- 0L

for (i in seq_along(cells)) {
  cl <- cells[[i]]
  fam <- eval(cl$fam)

  for (r in seq_len(REPS)) {
    t <- tryCatch(
      system.time(suppressMessages(suppressWarnings(
        bartisan(form(cl$y), d, family = fam, control = ctrl))))[["elapsed"]],
      error = function(e) NA_real_)

    k <- k + 1L
    rows[[k]] <- data.frame(label = cl$label, rep = r, secs = t,
                            stringsAsFactors = FALSE)
    prog_tick(pr, i = k, secs = if (is.na(t)) NULL else t,
              ok = !is.na(t), label = sprintf("%s rep %d", cl$label, r))
  }

  # Checkpointed per family, so an interrupted run still reports the families
  # it finished.
  saveRDS(list(res = do.call(rbind, rows), reps = REPS, n = N, p = P,
               complete = FALSE, done = i, total = length(cells)), OUT)
}

res <- do.call(rbind, rows)
on.exit()
prog_end(pr, "done")

# The vignette's absolute anchor, re-measured at the REAL defaults. It
# currently says "the default 1000 warmup plus 1000 saved draws"; the defaults
# are 200 and 800, so the anchor was quoting settings no caller gets by default.
anchor <- do.call(rbind, lapply(c(500L, 5000L), function(n) {
  set.seed(11L)
  Xa <- matrix(runif(n * P), n, P)
  colnames(Xa) <- paste0("x", seq_len(P))
  da <- as.data.frame(Xa)
  da$y <- 2 * sin(pi * da$x1) + 1.5 * da$x2 + rnorm(n, 0, 0.7)

  do.call(rbind, lapply(c("smoothstep", "hard"), function(g) {
    tt <- min(vapply(seq_len(3L), function(r) {
      system.time(suppressMessages(bartisan(
        form("y"), da, family = stats::gaussian(),
        control = bartisan_control(num_trees = 50L, gate = g))))[["elapsed"]]
    }, numeric(1L)))
    data.frame(n = n, gate = g, secs = tt, stringsAsFactors = FALSE)
  }))
}))

cat("
== absolute anchor, 50 trees, package defaults (200 warmup + 800 draws) ==
")
print(anchor, row.names = FALSE, digits = 3)

best <- tapply(res$secs, res$label, function(v) min(v, na.rm = TRUE))
base <- best[["gaussian()"]]
out <- data.frame(family = names(best), secs = round(as.numeric(best), 2),
                  relative = round(as.numeric(best) / base, 1),
                  stringsAsFactors = FALSE)
print(out[order(out$relative), ], row.names = FALSE)
saveRDS(list(res = res, anchor = anchor, reps = REPS, n = N, p = P,
             complete = TRUE, done = length(cells), total = length(cells)), OUT)
cat("\nGaussian baseline:", round(base, 2), "s at n =", N, ", p =", P,
    ", 50 trees, 200+400 draws\n")
cat("wrote", OUT, "\n")
