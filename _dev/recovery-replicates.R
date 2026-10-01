# Replicates for the two cells of _dev/recovery-matrix.R that read off their row.
#
# What is tested. In the single-replicate matrix, ordbeta() recovered an effect
# of 1.00 in both random-intercept cells against 0.79 and 0.86 without them,
# and zi_negbin() recovered 0.71 and 0.67 under vc() against 0.82 and 0.84
# without it. The truth is 0.8 in every cell. The claim: both are Monte Carlo
# noise from one replicate of 500 draws, as the zi_negbin cell flagged in
# _dev/AUDITING.md once was, and not a defect specific to those compositions.
#
# What is measured. The two families under plain, vc(), random intercepts and
# both, three data-and-fit seeds each, 50 trees, hard rules, 300 + 800 sweeps.
# Per fit, the posterior mean and posterior sd of the effect, the mean over
# units of link(z = 1) - link(z = 0) by draw. Read the cell mean over seeds
# against the family's other cells, and the posterior sd against the gap.
#
# What each outcome means. A cell mean within about 0.1 of its row with 0.8
# inside about 1.5 posterior sds was noise, and the matrix row stands. A gap
# that persists across the three seeds with the posterior sd well below it is
# a defect in that composition, to be investigated before shipping.
#
# Run with: Rscript _dev/recovery-replicates.R
# Writes:   _dev/recovery-replicates.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

OUT <- "_dev/recovery-replicates.rds"
SEEDS <- c(41L, 42L, 43L)
control <- bartisan_control(num_trees = 50L, num_burn = 300L, num_draws = 800L,
                            gate = "hard", verbose = FALSE)

simulate <- function(seed, n = 1000L) {
  set.seed(seed)
  d <- as.data.frame(matrix(runif(n * 3L), n, 3L))
  names(d) <- paste0("x", 1:3)
  d$z <- rbinom(n, 1L, 0.5)
  d$g <- factor(sample(letters[1:8], n, TRUE))
  u <- setNames(rnorm(8L, 0, 0.6), letters[1:8])
  f <- 1.2 * d$x1 - 1.2 * d$x2
  lp <- f - mean(f) + 0.8 * d$z + u[as.character(d$g)]

  mu <- plogis(lp)
  d$yob <- pmin(pmax(rbeta(n, mu * 10, (1 - mu) * 10), 1e-4), 1 - 1e-4)
  d$yob[lp < quantile(lp, 0.08)] <- 0
  d$yob[lp > quantile(lp, 0.92)] <- 1
  d$yzinb <- rnbinom(n, mu = exp(lp - 0.5), size = 3) * rbinom(n, 1L, 0.75)
  d
}

cases <- list(
  list(name = "ordbeta", family = ordbeta(), response = "yob"),
  list(name = "zi_negbin", family = zi_negbin(), response = "yzinb")
)
structures <- c("plain", "vc", "ranef", "vc + ranef")

formula_for <- function(response, structure) {
  terms <- switch(structure,
                  plain = c("z", "x1", "x2", "x3"),
                  vc = c("x1", "x2", "x3", "vc(z)"),
                  ranef = c("z", "x1", "x2", "x3", "(1 | g)"),
                  `vc + ranef` = c("x1", "x2", "x3", "vc(z)", "(1 | g)"))
  as.formula(paste(response, "~", paste(terms, collapse = " + ")),
             env = globalenv())
}

first_predictor <- function(p) {
  if (is.list(p)) p <- p[[1L]]
  p
}

grid <- expand.grid(case = seq_along(cases), structure = structures,
                    seed = SEEDS, stringsAsFactors = FALSE)

pr <- prog_init(total = nrow(grid), title = "Recovery replicates", unit = "fit",
                kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()

for (k in seq_len(nrow(grid))) {
  case <- cases[[grid$case[k]]]
  structure <- grid$structure[k]
  seed <- grid$seed[k]
  d <- simulate(seed)
  d1 <- d
  d1$z <- 1L
  d0 <- d
  d0$z <- 0L

  t0 <- Sys.time()
  set.seed(seed)
  fit <- bartisan(formula_for(case$response, structure), data = d,
                  family = case$family, control = control)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  e1 <- first_predictor(predict(fit, newdata = d1, type = "link", draws = TRUE))
  e0 <- first_predictor(predict(fit, newdata = d0, type = "link", draws = TRUE))
  if (length(dim(e1)) == 3L) {
    e1 <- e1[, , 1L]
    e0 <- e0[, , 1L]
  }
  by_draw <- rowMeans(e1 - e0)

  rows[[k]] <- data.frame(family = case$name, structure = structure, seed = seed,
                          effect = mean(by_draw), post_sd = sd(by_draw),
                          secs = secs)
  prog_tick(pr, label = sprintf("%s / %s / seed %d", case$name, structure, seed))
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid),
               done = k, total = nrow(grid)), OUT)
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
summary_tab <- aggregate(cbind(effect, post_sd) ~ family + structure, res, mean)
spread <- aggregate(effect ~ family + structure, res,
                    function(x) diff(range(x)))
names(spread)[3L] <- "seed_range"
summary_tab <- merge(summary_tab, spread)
print(summary_tab[order(summary_tab$family, summary_tab$structure), ],
      digits = 3, row.names = FALSE)
