# Confirmation for raising the `bandwidth_every` default from 1 to 10.
#
# `_dev/bandwidth-every.R` found the candidate: at every = 10 a fit is about
# 1.6x faster with no detectable change in held-out RMSE or interval coverage,
# and better effective sample size per second in both data shapes it tried.
# That was two shapes at one n, one p and one family, which is not a range a
# default should rest on. This varies n, the predictor count, the data shape and
# the family, and compares only the incumbent against the candidate so the
# replicate budget goes into breadth instead of into settings already ruled out.
#
# Paired throughout: both settings see the same training data in each replicate,
# so the difference is a within-replicate contrast and the test is a paired one.
A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))
Sys.setenv(PROGRESS_ROOT = "/Users/NoahGreifer/Dropbox/Research/R/bartisan/.progress-runs")

suppressMessages(library(bartisan))

REPS <- 10L
N_TEST <- 500L
OUT <- file.path("_dev", "bandwidth-confirm.rds")

mu_of <- function(X, shape) {
  if (identical(shape, "smooth")) {
    3 * sin(pi * X[, 1L]) + 2 * X[, 2L] + 3 * (X[, 3L] - 0.5)^2
  }
  else {
    3 * (X[, 1L] > 0.5) + 2 * (X[, 2L] > 0.3) - 1.5 * (X[, 3L] > 0.7)
  }
}

make <- function(n, p, shape, family, seed) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  d <- as.data.frame(X)
  mu <- mu_of(X, shape)

  if (identical(family, "binomial")) {
    d$y <- rbinom(n, 1L, stats::plogis(mu - mean(mu)))
    return(list(data = d, mu = stats::plogis(mu - mean(mu)), type = "response"))
  }

  if (identical(family, "gaussian_ls")) {
    d$y <- mu + rnorm(n, 0, exp(-0.9 + 0.8 * X[, 1L]))
    return(list(data = d, mu = mu, type = "link"))
  }

  d$y <- mu + rnorm(n, 0, 0.7)
  list(data = d, mu = mu, type = "link")
}

fam_of <- function(family) {
  switch(family, gaussian = stats::gaussian(), binomial = stats::binomial(),
         gaussian_ls = gaussian_ls())
}

quiet <- function(expr) {
  invisible(utils::capture.output(out <- suppressMessages(suppressWarnings(expr))))
  out
}

worst_ess <- function(fit) {
  tab <- diagnose(fit)[["table"]]
  tab <- tab[!grepl("average over observations", tab[["quantity"]], fixed = TRUE), ]
  min(tab[["ess_bulk"]], na.rm = TRUE)
}

cells <- rbind(
  expand.grid(n = c(250L, 1000L), p = c(10L, 50L),
              shape = c("smooth", "step"), family = "gaussian",
              stringsAsFactors = FALSE),
  expand.grid(n = 500L, p = 10L, shape = c("smooth", "step"),
              family = c("gaussian_ls", "binomial"), stringsAsFactors = FALSE)
)

n_total <- nrow(cells) * REPS * 2L

pr <- prog_init(total = n_total,
                title = "bandwidth_every = 10 confirmation", unit = "fit",
                kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
k <- 0L

for (i in seq_len(nrow(cells))) {
  cl <- cells[i, ]
  test <- make(N_TEST, cl$p, cl$shape, cl$family, seed = 400L + i)

  for (r in seq_len(REPS)) {
    train <- make(cl$n, cl$p, cl$shape, cl$family,
                  seed = 20000L + 211L * r + 13L * i)

    for (every in c(1L, 10L)) {
      ctrl <- bartisan_control(num_trees = 50L, num_burn = 1000L,
                               num_draws = 1000L, chains = 4L,
                               bandwidth_every = every)
      t0 <- Sys.time()
      fit <- quiet(bartisan(y ~ ., train$data, family = fam_of(cl$family),
                            control = ctrl))
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

      dr <- stats::predict(fit, newdata = test$data, type = test$type,
                           draws = TRUE)
      # A location-scale fit returns a list of predictors; the mean is first.
      if (is.list(dr)) dr <- dr[[1L]]

      post <- colMeans(dr)
      lo <- apply(dr, 2L, stats::quantile, 0.025)
      hi <- apply(dr, 2L, stats::quantile, 0.975)
      ess <- quiet(worst_ess(fit))
      k <- k + 1L

      rows[[k]] <- data.frame(
        n = cl$n, p = cl$p, shape = cl$shape, family = cl$family,
        every = every, rep = r, secs = secs,
        rmse = sqrt(mean((post - test$mu)^2)),
        coverage = mean(lo <= test$mu & test$mu <= hi),
        ess_min = ess, ess_per_sec = ess / secs,
        stringsAsFactors = FALSE)

      prog_tick(pr, i = k, secs = secs,
                label = sprintf("%s %s n=%d p=%d every=%d r%d", cl$family,
                                cl$shape, cl$n, cl$p, every, r))
    }

    # After every cell rather than after the last one, so a killed run leaves
    # its finished cells readable. `complete` is what tells a reader which it
    # is looking at.
    saveRDS(list(res = do.call(rbind, rows), reps = REPS,
                 complete = FALSE, done = k, total = n_total), OUT)
  }
}

res <- do.call(rbind, rows)
saveRDS(list(res = res, reps = REPS, complete = TRUE, done = k,
             total = n_total), OUT)
on.exit()
prog_end(pr, "done")

# Paired contrasts, per cell.
key <- with(res, paste(family, shape, n, p))
out <- do.call(rbind, lapply(split(res, key), function(z) {
  w <- reshape(z[, c("rep", "every", "rmse", "coverage", "secs", "ess_per_sec")],
               idvar = "rep", timevar = "every", direction = "wide")
  tt <- stats::t.test(w[["rmse.10"]] - w[["rmse.1"]])
  data.frame(cell = z$family[1], shape = z$shape[1], n = z$n[1], p = z$p[1],
             speedup = mean(w[["secs.1"]]) / mean(w[["secs.10"]]),
             rmse_ratio = mean(w[["rmse.10"]]) / mean(w[["rmse.1"]]),
             rmse_p = tt$p.value,
             cover_1 = mean(w[["coverage.1"]]),
             cover_10 = mean(w[["coverage.10"]]),
             esss_ratio = mean(w[["ess_per_sec.10"]]) / mean(w[["ess_per_sec.1"]]),
             stringsAsFactors = FALSE)
}))
print(out[order(out$cell, out$shape, out$n, out$p), ], row.names = FALSE, digits = 3)
cat("\nwrote", OUT, "\n")
