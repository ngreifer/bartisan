# Follow-up to _dev/recovery-replicates.R: the two cells that stayed off.
#
# What is tested. With three seeds and 800 draws, ordbeta() overshot an effect
# of 0.8 in every structure (0.86 and 0.93 without random intercepts, 1.04 and
# 1.06 with them, posterior sd 0.055) and its point estimates swung across
# seeds by six to eight posterior sds; zi_negbin() read about 0.1 lower under
# vc() than without (0.65 and 0.62 against 0.76 and 0.72, posterior sd 0.12).
# Two claims. First, the ordbeta() inflation is the generator's: the recovery
# matrix set y to 0 or 1 by a deterministic cut on the predictor, a step the
# ordered beta model represents with logistic boundary probabilities, so the
# fit inflates the latent scale to approximate it. Second, the zi_negbin() dip
# is weak identification between the count and zero processes, which a
# dedicated coefficient forest resolves differently from a shared covariate
# and which improves with n.
#
# What is measured. ordbeta() under plain, vc(), random intercepts and both,
# three seeds, with the response drawn from the ordered beta model itself:
# P(y = 0) = 1 - plogis(lp - c1), P(y = 1) = plogis(lp - c2), c1 = -1.5,
# c2 = 1.5, and Beta(mu phi, (1 - mu) phi) with mu = plogis(lp), phi = 10
# otherwise. zi_negbin() under plain and vc() at n = 4000, two seeds. 50
# trees, hard rules, 300 + 800 draws. Per fit the effect's posterior mean and
# posterior sd, as in the replicates.
#
# What each outcome means. ordbeta() effects within about 0.1 of 0.8 in all
# four cells, with seed ranges of the order of the posterior sd, clear the row
# and fault the generator; inflation that persists is a defect in the family,
# and the random-intercept composition is where to look first. A zi_negbin()
# gap between vc() and plain that halves or vanishes at n = 4000 is
# identification and the row stands; one that holds at 0.1 is a defect in the
# composition.
#
# Run with: Rscript _dev/recovery-followup.R
# Writes:   _dev/recovery-followup.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

OUT <- "_dev/recovery-followup.rds"
control <- bartisan_control(num_trees = 50L, num_burn = 300L, num_draws = 800L,
                            gate = "hard", verbose = FALSE)

simulate <- function(seed, n) {
  set.seed(seed)
  d <- as.data.frame(matrix(runif(n * 3L), n, 3L))
  names(d) <- paste0("x", 1:3)
  d$z <- rbinom(n, 1L, 0.5)
  d$g <- factor(sample(letters[1:8], n, TRUE))
  u <- setNames(rnorm(8L, 0, 0.6), letters[1:8])
  f <- 1.2 * d$x1 - 1.2 * d$x2
  lp <- f - mean(f) + 0.8 * d$z + u[as.character(d$g)]

  # The ordered beta model's own boundary mechanism.
  mu <- plogis(lp)
  p0 <- 1 - plogis(lp + 1.5)
  p1 <- plogis(lp - 1.5)
  which <- runif(n)
  d$yob <- pmin(pmax(rbeta(n, mu * 10, (1 - mu) * 10), 1e-4), 1 - 1e-4)
  d$yob[which < p0] <- 0
  d$yob[which > 1 - p1] <- 1

  d$yzinb <- rnbinom(n, mu = exp(lp - 0.5), size = 3) * rbinom(n, 1L, 0.75)
  d
}

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
  if (length(dim(p)) == 3L) p <- p[, , 1L]
  p
}

grid <- rbind(
  expand.grid(family = "ordbeta", response = "yob", n = 1000L,
              structure = c("plain", "vc", "ranef", "vc + ranef"),
              seed = c(41L, 42L, 43L), stringsAsFactors = FALSE),
  expand.grid(family = "zi_negbin", response = "yzinb", n = 4000L,
              structure = c("plain", "vc"), seed = c(41L, 42L),
              stringsAsFactors = FALSE))

families <- list(ordbeta = ordbeta(), zi_negbin = zi_negbin())

pr <- prog_init(total = nrow(grid), title = "Recovery follow-up", unit = "fit",
                kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()

for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]
  d <- simulate(g$seed, g$n)
  d1 <- d
  d1$z <- 1L
  d0 <- d
  d0$z <- 0L

  t0 <- Sys.time()
  set.seed(g$seed)
  fit <- bartisan(formula_for(g$response, g$structure), data = d,
                  family = families[[g$family]], control = control)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  e1 <- first_predictor(predict(fit, newdata = d1, type = "link", draws = TRUE))
  e0 <- first_predictor(predict(fit, newdata = d0, type = "link", draws = TRUE))
  by_draw <- rowMeans(e1 - e0)

  rows[[k]] <- data.frame(family = g$family, n = g$n, structure = g$structure,
                          seed = g$seed, effect = mean(by_draw),
                          post_sd = sd(by_draw), secs = secs)
  prog_tick(pr, label = sprintf("%s / n %d / %s / seed %d", g$family, g$n,
                                g$structure, g$seed))
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid),
               done = k, total = nrow(grid)), OUT)
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
summary_tab <- aggregate(cbind(effect, post_sd) ~ family + n + structure, res, mean)
spread <- aggregate(effect ~ family + n + structure, res, function(x) diff(range(x)))
names(spread)[4L] <- "seed_range"
summary_tab <- merge(summary_tab, spread)
print(summary_tab[order(summary_tab$family, summary_tab$structure), ],
      digits = 3, row.names = FALSE)
