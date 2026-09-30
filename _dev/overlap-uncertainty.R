# Does the interval for tau(x) widen where the treatment groups do not overlap?
#
# What is being tested. Li, Ding and Mealli (2023, Example 4.1, reproduced from
# Papadogeorgou and Li's discussion of Hahn, Murray and Carvalho 2020) show a
# BART interval for the CATE staying the same width whatever the overlap, which
# makes it overconfident where only one group has data, and say soft decision
# trees (Linero and Yang 2018) can mitigate this. bartisan's rules are soft by
# default. The claim that could be false: bartisan's soft-rule interval widens
# in the one-group tails and the hard-rule interval does not.
#
# The measurement. Their design: 250 treated with X ~ Gamma(mean 35, sd 8), 250
# controls with X ~ Gamma(mean 60, sd 8), Y(z) = 10 + 5z - 0.3X + N(0, 1), so
# tau(x) = 5 everywhere. Ten replications. Fits: bartisan() with z as a
# predictor under hard and under soft rules, bcf(), and the paper's own design,
# a separate fit per arm (a T-learner), under hard and soft rules. Read the posterior sd
# of tau(x) and the coverage of 5 by the 95% interval in three regions: the
# treated-only tail (x 20-25), the overlap (45-50) and the control-only tail
# (75-80).
#
# What each outcome means. Soft sd rising in the tails and hard staying flat:
# the mitigation holds, causal.Rmd's "the interval will not tell you" is too
# strong for the default rules, and an overlap diagnostic built on posterior
# uncertainty is viable. Both flat: the sentence stands and an overlap check
# has to come from the propensity score.
#
# Result (2026-09-25): the paper's finding reproduces for its own design, a
# separate fit per arm, and soft rules do not fix it there: posterior sd 0.58
# in the overlap against 0.68 and 0.86 in the tails (hard), 0.39 against 0.52
# and 0.68 (soft), with bias of 5 to 6 and coverage 0 in both tails. The
# single-forest designs bartisan uses do widen: 0.34 against 1.43 and 1.20
# (hard), 0.26 against 1.48 and 1.08 (soft), 0.31 against 2.86 and 2.09
# (bcf), covering in the tails, with bias 0.5 to 0.6 (hard), under 0.1
# (soft) and 1.4 to 1.9 (bcf). One simple design (linear in x, constant
# effect); the single-forest result may flatter a problem with less shared
# structure across the arms.

source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")
suppressMessages(library(bartisan))

reps <- 10L
grid <- data.frame(x = c(20, 22.5, 25, 45, 47.5, 50, 75, 77.5, 80),
                   region = rep(c("treated only", "overlap", "control only"), each = 3))
pr <- prog_init(total = reps, title = "Overlap and CATE uncertainty", unit = "replicate",
                kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)
out_file <- "_dev/overlap-uncertainty.rds"
rows <- list()

tau_draws <- function(fit, extra = NULL) {
  g1 <- cbind(grid["x"], z = 1); g0 <- cbind(grid["x"], z = 0)
  predict(fit, newdata = g1, draws = TRUE) - predict(fit, newdata = g0, draws = TRUE)
}

for (r in seq_len(reps)) {
  set.seed(1000 + r)
  z <- rep(1:0, each = 250)
  x <- ifelse(z == 1, rgamma(500, shape = (35 / 8)^2, rate = 35 / 64),
              rgamma(500, shape = (60 / 8)^2, rate = 60 / 64))
  d <- data.frame(x = x, z = z, y = 10 + 5 * z - 0.3 * x + rnorm(500))

  fits <- list(
    hard = bartisan(y ~ z + x, data = d, family = gaussian(), gate = "hard", verbose = FALSE),
    soft = bartisan(y ~ z + x, data = d, family = gaussian(), verbose = FALSE),
    bcf  = suppressMessages(bcf(y ~ x, treat = ~ z, data = d, family = gaussian(),
                                verbose = FALSE))
  )

  # The paper's design: one model per arm, so each extrapolates on its own.
  arm <- function(k, gate) {
    ctrl <- if (gate == "hard") bartisan_control(gate = "hard", verbose = FALSE) else
      bartisan_control(verbose = FALSE)
    bartisan(y ~ x, data = d[d$z == k, ], family = gaussian(), control = ctrl)
  }
  t_fits <- list(t_hard = list(arm(1, "hard"), arm(0, "hard")),
                 t_soft = list(arm(1, "soft"), arm(0, "soft")))

  for (m in names(fits)) {
    td <- tau_draws(fits[[m]])
    q <- apply(td, 2, quantile, c(.025, .975))
    rows[[length(rows) + 1L]] <- data.frame(rep = r, model = m, grid,
                                            mean = colMeans(td), sd = apply(td, 2, sd),
                                            covered = q[1, ] <= 5 & 5 <= q[2, ])
  }
  for (m in names(t_fits)) {
    td <- predict(t_fits[[m]][[1]], newdata = grid["x"], draws = TRUE) -
      predict(t_fits[[m]][[2]], newdata = grid["x"], draws = TRUE)
    q <- apply(td, 2, quantile, c(.025, .975))
    rows[[length(rows) + 1L]] <- data.frame(rep = r, model = m, grid,
                                            mean = colMeans(td), sd = apply(td, 2, sd),
                                            covered = q[1, ] <= 5 & 5 <= q[2, ])
  }
  prog_tick(pr, label = sprintf("replicate %d", r))
  saveRDS(list(res = do.call(rbind, rows), complete = FALSE, done = r, total = reps), out_file)
}

res <- do.call(rbind, rows)
saveRDS(list(res = res, complete = TRUE, done = reps, total = reps), out_file)
on.exit()
prog_end(pr, "done")

res$region <- factor(res$region, c("treated only", "overlap", "control only"))
a <- aggregate(cbind(sd, covered, bias = mean - 5) ~ model + region, res, mean)
print(a[order(a$model, a$region), ], row.names = FALSE, digits = 3)
