# Which family should model a continuous treatment's propensity score in bcf()?
#
# What is being tested. For a continuous treatment, bcf() adds the fitted
# conditional mean E[A | X] to the control function, fit with gaussian(). The
# question (2026-09-25) is whether dpm(), which bartisan() picks for a numeric
# response, should be the default instead. The claim at stake: dpm() recovers
# E[A | X] at least as well as gaussian() under normal errors and better under
# heavy-tailed or skewed ones. It could be false because dpm() assumes an error
# distribution that does not depend on x, so a spread that varies with x or a
# point mass could bias its fitted mean.
#
# The measurement. A nonlinear E[A | X] in five uniform covariates, n = 500, ten
# replications per design: normal errors, t3 errors, centered log-normal
# (skewed) errors, errors whose sd grows with x1, and a semi-continuous
# treatment with a covariate-dependent share of zeros and a gamma positive part,
# where tweedie() is fit as well. Read the Spearman correlation between the
# fitted and the true E[A | X], since a tree uses only the ordering of a
# covariate; also the RMSE after centering each at its mean, and seconds.
#
# What each outcome means. dpm() at least as good everywhere and better under
# heavy tails or skew: it becomes the default. Worse anywhere, particularly
# under the varying spread or the point mass: gaussian() stays, and the docs
# point to the override. The semi-continuous design says whether tweedie() is
# the family to recommend for such treatments.

source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")
suppressMessages(library(bartisan))

reps <- as.integer(Sys.getenv("REPS", "10"))
n <- 500L
designs <- c("normal", "t3", "skewed", "heteroskedastic", "semicontinuous")
families <- function(design) {
  if (design == "semicontinuous") c("gaussian", "dpm", "tweedie") else c("gaussian", "dpm")
}
total <- reps * sum(vapply(designs, function(d) length(families(d)), 1L))

m_fun <- function(x) 2 * sin(pi * x$x1 * x$x2) + 4 * (x$x3 - 0.5)^2 + 2 * x$x4 + x$x5

simulate <- function(design) {
  x <- as.data.frame(matrix(stats::runif(n * 5), n, 5, dimnames = list(NULL, paste0("x", 1:5))))
  m <- m_fun(x)
  a <- switch(design,
    normal = m + stats::rnorm(n),
    t3 = m + stats::rt(n, 3) / sqrt(3),
    skewed = {
      e <- exp(stats::rnorm(n, 0, 0.75)); m + (e - exp(0.75^2 / 2)) / sd(e)
    },
    heteroskedastic = m + stats::rnorm(n, 0, 0.3 + 1.5 * x$x1),
    semicontinuous = {
      p <- stats::plogis(-1 + 2 * x$x2 + x$x3)
      mu <- 0.5 + m / 2
      pos <- stats::rbinom(n, 1L, p)
      a <- pos * stats::rgamma(n, shape = 2, rate = 2 / mu)
      m <- p * mu          # the true conditional mean
      a
    })
  x$a <- a
  list(data = x, truth = m)
}

fam <- list(gaussian = gaussian(), dpm = dpm(), tweedie = tweedie())

pr <- prog_init(total = total, title = "Propensity family for a continuous treatment",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)
out_file <- Sys.getenv("OUT", "_dev/propensity-family.rds")
rows <- list()

for (design in designs) {
  for (r in seq_len(reps)) {
    set.seed(100 * match(design, designs) + r)
    s <- simulate(design)
    for (f in families(design)) {
      t0 <- Sys.time()
      fit <- suppressMessages(bartisan(a ~ x1 + x2 + x3 + x4 + x5, data = s$data,
                                       family = fam[[f]], verbose = FALSE))
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      est <- as.numeric(stats::fitted(fit))
      rows[[length(rows) + 1L]] <- data.frame(
        design = design, rep = r, family = f, seconds = round(secs, 2),
        spearman = stats::cor(est, s$truth, method = "spearman"),
        rmse_centered = sqrt(mean(((est - mean(est)) - (s$truth - mean(s$truth)))^2)))
      prog_tick(pr, label = sprintf("%s %d %s", design, r, f))
    }
    saveRDS(list(res = do.call(rbind, rows), complete = FALSE), out_file)
  }
}

res <- do.call(rbind, rows)
saveRDS(list(res = res, complete = TRUE), out_file)
on.exit()
prog_end(pr, "done")

a <- aggregate(cbind(spearman, rmse_centered, seconds) ~ design + family, res, mean)
a$design <- factor(a$design, designs)
print(a[order(a$design, a$family), ], row.names = FALSE, digits = 3)
