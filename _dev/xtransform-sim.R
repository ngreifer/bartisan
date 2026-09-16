# Is `x_transform = "range"` an acceptable default?
#
# The transform maps each numeric predictor to [0, 1]; cutpoints are uniform
# there and, with soft rules, the gate bandwidth is measured there too. So the
# two options fail for opposite reasons, and the scenarios below are built to
# make each failure happen.
#
#   `range` spends the cutpoint prior on the RANGE of x, so a heavy tail or a
#   handful of outliers leaves the bulk of the data inside a sliver of [0, 1]:
#   most candidate cutpoints fall in empty space, and a bandwidth that is a
#   fixed fraction of [0, 1] is enormous relative to the bulk.
#
#   `quantile` spends it on the RANKS, so the resolution in x is the spacing of
#   the order statistics and nothing finer can be expressed. Where the data are
#   sparse the model is coarse, and a truth that is smooth in x becomes sharp in
#   u wherever the density changes quickly.
source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")
suppressMessages(library(bartisan))

SCRATCH <- "/private/tmp/claude-501/-Users-NoahGreifer-Dropbox-Research-R-bartisan/42a89c8f-bd1a-4a13-9717-48692000a78c/scratchpad/xt"
OUT <- file.path(SCRATCH, "results.rds")

N_REP <- 20L
N_TEST <- 1000L
SIZES <- c(500L, 2000L)

# Each generator returns the predictors and the true regression function.
scenarios <- list(
  uniform = list(
    why = "neutral: the ECDF is nearly the identity, so the two should tie",
    gen = function(n) {
      x <- data.frame(x1 = runif(n), x2 = runif(n))
      list(x = x, f = sin(2 * pi * x$x1) + 0.5 * x$x2)
    }),

  lognormal = list(
    why = "mild skew, the ordinary case quantile is meant for",
    gen = function(n) {
      x <- data.frame(x1 = rlnorm(n, 0, 1), x2 = runif(n))
      list(x = x, f = 2 * plogis(log(x$x1)) + 0.5 * x$x2)
    }),

  outliers = list(
    why = "3% of x is 25x wider; range's cutpoints and bandwidth land in empty space",
    gen = function(n) {
      big <- rbinom(n, 1L, 0.03)
      x <- data.frame(x1 = ifelse(big == 1L, rnorm(n, 0, 25), rnorm(n, 0, 1)),
                      x2 = runif(n))
      list(x = x, f = 2 * plogis(x$x1) + 0.5 * x$x2)
    }),

  pareto = list(
    why = "infinite-variance tail, the extreme version of the same problem",
    gen = function(n) {
      x <- data.frame(x1 = (1 - runif(n))^(-1 / 1.5), x2 = runif(n))
      list(x = x, f = 2 * plogis(log(x$x1) - 0.5) + 0.5 * x$x2)
    }),

  bimodal_linear = list(
    why = "two tight clusters far apart with a LINEAR truth; the ECDF jumps .5 across the gap, so linear in x is a step in u",
    gen = function(n) {
      hi <- rbinom(n, 1L, 0.5)
      x <- data.frame(x1 = ifelse(hi == 1L, rnorm(n, 10, 0.1), rnorm(n, 0, 0.1)),
                      x2 = runif(n))
      list(x = x, f = 0.3 * x$x1 + 0.5 * x$x2)
    }),

  sparse_tail = list(
    why = "the interesting feature sits in the sparse upper tail, where quantile has few order statistics to resolve it",
    gen = function(n) {
      x <- data.frame(x1 = rlnorm(n, 0, 1), x2 = runif(n))
      list(x = x, f = 2 * plogis(5 * (x$x1 - 3.6)) + 0.5 * x$x2)
    }),

  outliers_fine = list(
    why = "1% of x is 300x wider AND the truth oscillates inside the bulk, so range must resolve fine structure in 1% of its coordinate",
    gen = function(n) {
      big <- rbinom(n, 1L, 0.01)
      x <- data.frame(x1 = ifelse(big == 1L, rnorm(n, 0, 300), rnorm(n, 0, 1)),
                      x2 = runif(n))
      list(x = x, f = sin(3 * x$x1) + 0.5 * x$x2)
    }),

  mixed = list(
    why = "five predictors of different shapes at once, which is what real data looks like",
    gen = function(n) {
      big <- rbinom(n, 1L, 0.03)
      x <- data.frame(x1 = rlnorm(n, 0, 1),
                      x2 = ifelse(big == 1L, rnorm(n, 0, 25), rnorm(n, 0, 1)),
                      x3 = runif(n), x4 = runif(n), x5 = rnorm(n))
      list(x = x, f = 1.5 * plogis(log(x$x1)) + sin(2 * pi * x$x3) + 0.1 * x$x2)
    })
)

ctl <- function(tr) bartisan_control(num_trees = 50L, num_burn = 200L,
                                     num_draws = 500L, x_transform = tr,
                                     verbose = FALSE)

one_cell <- function(scn, n, rep) {
  set.seed(1000L * rep + n)
  tr <- scenarios[[scn]]$gen(n)
  te <- scenarios[[scn]]$gen(N_TEST)
  sdf <- stats::sd(tr$f)
  y <- tr$f + rnorm(n, sd = 0.3 * sdf)
  train <- cbind(tr$x, y = y)

  out <- NULL
  for (tf in c("quantile", "range")) {
    fit <- bartisan(y ~ ., train, family = gaussian(), control = ctl(tf))
    dr <- predict(fit, newdata = te$x, type = "response", draws = TRUE)
    fhat <- colMeans(dr)
    lo <- apply(dr, 2L, quantile, 0.025)
    hi <- apply(dr, 2L, quantile, 0.975)
    out <- rbind(out, data.frame(
      scenario = scn, n = n, rep = rep, transform = tf,
      rmse = sqrt(mean((fhat - te$f)^2)) / sdf,     # relative to the signal sd
      cover = mean(te$f >= lo & te$f <= hi),
      width = mean(hi - lo) / sdf))
  }
  out
}

cells <- expand.grid(scn = names(scenarios), n = SIZES,
                     stringsAsFactors = FALSE)
total <- nrow(cells) * N_REP

pr <- prog_init(total = total, title = "x_transform simulation",
                unit = "fit pair", workers = 1L, kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
for (k in seq_len(nrow(cells))) {
  for (r in seq_len(N_REP)) {
    res <- tryCatch(one_cell(cells$scn[k], cells$n[k], r),
                    error = function(e) {
                      data.frame(scenario = cells$scn[k], n = cells$n[k],
                                 rep = r, transform = NA_character_,
                                 rmse = NA_real_, cover = NA_real_,
                                 width = NA_real_)
                    })
    rows[[length(rows) + 1L]] <- res
    prog_tick(pr, label = sprintf("%s n=%d rep %d", cells$scn[k], cells$n[k], r))
  }
  # After every cell, so a killed run still leaves every finished one behind.
  saveRDS(list(res = do.call(rbind, rows), complete = FALSE,
               scenarios = lapply(scenarios, `[[`, "why")),
          OUT)
}

saveRDS(list(res = do.call(rbind, rows), complete = TRUE,
             scenarios = lapply(scenarios, `[[`, "why")), OUT)
on.exit()
prog_end(pr, "done")
cat("done:", nrow(do.call(rbind, rows)), "rows\n")
