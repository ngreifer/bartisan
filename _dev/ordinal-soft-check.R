# Is the soft-rule attenuation of the ordinal effect mixing or structural?
#
# What is tested. The soft-rule arm of _dev/recovery-matrix.R (2026-10-01) put
# the three ordinal rows about 30% below their hard-rule values in every
# structure alike (logit 0.54 to 0.67 against 0.82 to 0.89 for a latent effect
# of 0.8; probit 0.34 to 0.37 against 0.50 to 0.53), while binomial() under the
# same latent samplers did not move. Two readings: the cutpoints and the
# predictor are still settling at 800 draws under soft rules, or the soft-rule
# ordinal fit is biased toward zero.
#
# What is measured. ordinal() and ordinal("probit"), plain structure, hard
# against smoothstep, 300 + 800 against 300 + 3200 draws, on the matrix's own
# data set and on two fresh ones. Per fit the effect's posterior mean and sd, the mean over units
# of link(z = 1) - link(z = 0) by draw.
#
# What each outcome means. A soft-rule effect that moves toward the hard-rule
# value at 3200 draws is mixing, and the default length is what to record. One
# that stays near 70% of it is a bias of the soft-rule ordinal fit, to be
# investigated on its own before shipping.
#
# Run with: Rscript _dev/ordinal-soft-check.R
# Writes:   _dev/ordinal-soft-check.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

simulate <- function(seed, n = 1000L) {
  set.seed(seed)
  d <- as.data.frame(matrix(runif(n * 3L), n, 3L))
  names(d) <- paste0("x", 1:3)
  d$z <- rbinom(n, 1L, 0.5)
  d$g <- factor(sample(letters[1:8], n, TRUE))
  u <- setNames(rnorm(8L, 0, 0.6), letters[1:8])
  f <- 1.2 * d$x1 - 1.2 * d$x2
  lp <- f - mean(f) + 0.8 * d$z + u[as.character(d$g)]
  d$yord <- cut(lp + rlogis(n), c(-Inf, -0.5, 0.5, Inf), labels = 1:3,
                ordered_result = TRUE)
  d
}

# The matrix's own data set, which draws the ordinal response after several
# other responses in its generator, so that a seed alone does not reproduce it.
# Reproduced here step for step so that the gates can be compared on the data
# that showed the gap; seed 0 names it below.
simulate_matrix <- function(n = 1000L) {
  set.seed(41)
  d <- as.data.frame(matrix(runif(n * 3L), n, 3L))
  names(d) <- paste0("x", 1:3)
  d$z <- rbinom(n, 1L, 0.5)
  d$g <- factor(sample(letters[1:8], n, TRUE))
  u <- setNames(rnorm(8L, 0, 0.6), letters[1:8])
  f <- 1.2 * d$x1 - 1.2 * d$x2
  lp <- f - mean(f) + 0.8 * d$z + u[as.character(d$g)]
  d$ynorm <- lp + rnorm(n, 0, 0.5)
  d$ybin <- rbinom(n, 1L, plogis(lp))
  d$yprob <- rbinom(n, 1L, pnorm(lp))
  d$ycll <- rbinom(n, 1L, 1 - exp(-exp(lp)))
  d$ycnt <- rpois(n, exp(lp - 0.5))
  d$ynb <- rnbinom(n, mu = exp(lp - 0.5), size = 3)
  d$yzi <- d$ycnt * rbinom(n, 1L, 0.75)
  d$yzinb <- d$ynb * rbinom(n, 1L, 0.75)
  d$ygam <- rgamma(n, shape = 4, rate = 4 / exp(lp - 0.5))
  mu <- plogis(lp)
  d$ybeta <- pmin(pmax(rbeta(n, mu * 10, (1 - mu) * 10), 1e-4), 1 - 1e-4)
  p0 <- 1 - plogis(lp + 1.5)
  p1 <- plogis(lp - 1.5)
  which_bound <- runif(n)
  d$yob <- d$ybeta
  d$yob[which_bound < p0] <- 0
  d$yob[which_bound > 1 - p1] <- 1
  events <- rpois(n, exp(lp - 0.5))
  d$ytw <- ifelse(events > 0, rgamma(n, shape = 2 * pmax(events, 1), rate = 2), 0)
  d$yord <- cut(lp + rlogis(n), c(-Inf, -0.5, 0.5, Inf), labels = 1:3,
                ordered_result = TRUE)
  d
}

grid <- expand.grid(link = c("logit", "probit"), gate = c("hard", "smoothstep"),
                    draws = c(800L, 3200L), seed = c(0L, 41L, 42L),
                    stringsAsFactors = FALSE)

pr <- prog_init(total = nrow(grid), title = "Ordinal under soft rules",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()

for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]
  d <- if (g$seed == 0L) simulate_matrix() else simulate(g$seed)
  d1 <- d
  d1$z <- 1L
  d0 <- d
  d0$z <- 0L

  t0 <- Sys.time()
  set.seed(g$seed)
  fit <- bartisan(yord ~ z + x1 + x2 + x3, data = d, family = ordinal(g$link),
                  control = bartisan_control(num_trees = 50L, num_burn = 300L,
                                             num_draws = g$draws, gate = g$gate,
                                             verbose = FALSE))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  by_draw <- rowMeans(predict(fit, newdata = d1, type = "link", draws = TRUE) -
                        predict(fit, newdata = d0, type = "link", draws = TRUE))

  rows[[k]] <- data.frame(link = g$link, gate = g$gate, draws = g$draws,
                          seed = g$seed, effect = mean(by_draw),
                          post_sd = sd(by_draw), secs = secs)
  prog_tick(pr, label = sprintf("%s / %s / %d / seed %d", g$link, g$gate,
                                g$draws, g$seed))
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid)),
          "_dev/ordinal-soft-check.rds")
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
cat("the matrix's data (seed 0):\n")
print(res[res$seed == 0L, c("link", "gate", "draws", "effect", "post_sd")],
      digits = 3, row.names = FALSE)
cat("\nfresh data, mean of seeds 41 and 42:\n")
out <- aggregate(cbind(effect, post_sd, secs) ~ link + gate + draws,
                 res[res$seed != 0L, ], mean)
print(out[order(out$link, out$gate, out$draws), ], digits = 3, row.names = FALSE)
