# Does the augmented binomial mix worse under hard rules than under soft?
#
# What is tested. The weights arm of the reseeded recovery matrix
# (`_dev/recovery-matrix.R`, 2026-10-01) read 1.53 for binomial logit under
# hard rules against 1.13 under soft, and `_dev/weights-check.R` traced that to
# chain-to-chain variation rather than to weights: the weighted and replicated
# fits agreed under soft rules to 0.03 while under hard rules their gap ran
# -0.25 to +0.18 across four seeds with no consistent sign, against posterior
# standard deviations of 0.1 to 0.17. The claim: the hard-rule fits of the
# augmented binomial families mix badly enough on this design that the
# posterior mean of a treatment effect moves by about a posterior standard
# deviation from seed to seed, and the soft-rule fits do not.
#
# What is measured. The recovery matrix's data and its plain structure,
# `y ~ z + x1 + x2 + x3` at n = 1000, 50 trees, 300 warmup and 500 kept draws,
# one chain. Binomial logit and complementary log-log, with the augmentation on
# (the default, a Polya-Gamma draw for the logit and an exponential waiting
# time for the complementary log-log) and off, under hard and smoothstep rules,
# four seeds each. Per fit, over the 500 retained draws of the effect, the mean
# over units of link(z = 1) - link(z = 0):
#
#   effect      its posterior mean
#   ess         bulk effective sample size of the effect
#   post_sd     its posterior standard deviation
#
# and per cell the spread of `effect` across the four seeds, read against the
# median `post_sd`.
#
# What each outcome would mean. Effective sample sizes in the tens for the
# hard-rule augmented cells where the soft-rule ones are in the hundreds makes
# the seed spread a mixing property of that combination, which belongs in
# `vignette("diagnostics")` as a caution, since hard rules are what the speed
# benchmarks recommend. Comparable effective sample sizes across the four cells
# would mean the spread is something else and the augmentation is not what to
# look at.
#
# Run with: Rscript _dev/augmented-mixing.R
# Writes:   _dev/augmented-mixing.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

drawn <- function(k, expr) {
  set.seed(41L + k)
  expr
}
n <- 1000L
d <- drawn(0L, {
  d <- as.data.frame(matrix(runif(n * 3L), n, 3L))
  names(d) <- paste0("x", 1:3)
  d$z <- rbinom(n, 1L, 0.5)
  d$g <- factor(sample(letters[1:8], n, TRUE))
  d$w <- sample(1:3, n, TRUE)
  d
})
u <- drawn(1L, setNames(rnorm(8L, 0, 0.6), letters[1:8]))
f <- 1.2 * d$x1 - 1.2 * d$x2
lp <- f - mean(f) + 0.8 * d$z + u[as.character(d$g)]
d$ybin <- drawn(3L, rbinom(n, 1L, plogis(lp)))
d$ycll <- drawn(5L, rbinom(n, 1L, 1 - exp(-exp(lp))))

d1 <- d
d1$z <- 1L
d0 <- d
d0$z <- 0L

cases <- list(
  list(name = "logit", family = binomial(), response = "ybin"),
  list(name = "cloglog", family = binomial("cloglog"), response = "ycll")
)
grid <- expand.grid(case = seq_along(cases), gate = c("hard", "smoothstep"),
                    augment = c(TRUE, FALSE), seed = 1:4,
                    stringsAsFactors = FALSE)

pr <- prog_init(total = nrow(grid), title = "Augmented binomial: mixing by gate",
                unit = "fit", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]
  case <- cases[[g$case]]
  form <- as.formula(paste(case$response, "~ z + x1 + x2 + x3"), env = globalenv())
  control <- bartisan_control(num_trees = 50L, num_burn = 300L, num_draws = 500L,
                              gate = g$gate, augment = g$augment, verbose = FALSE)
  t0 <- Sys.time()
  set.seed(100L + g$seed)
  fit <- bartisan(form, data = d, family = case$family, control = control)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  by_draw <- rowMeans(predict(fit, newdata = d1, type = "link", draws = TRUE) -
                        predict(fit, newdata = d0, type = "link", draws = TRUE))

  rows[[k]] <- data.frame(family = case$name, gate = g$gate, augment = g$augment,
                          seed = g$seed, effect = mean(by_draw),
                          ess = posterior::ess_bulk(matrix(by_draw, ncol = 1L)),
                          post_sd = sd(by_draw), secs = secs)
  prog_tick(pr, label = sprintf("%s / %s / augment %s / seed %d", case$name,
                                g$gate, g$augment, g$seed), secs = secs)
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid)),
          "_dev/augmented-mixing.rds")
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
out <- do.call(rbind, lapply(split(res, list(res$family, res$gate, res$augment)),
  function(x) {
    if (!nrow(x)) return(NULL)
    data.frame(family = x$family[1L], gate = x$gate[1L], augment = x$augment[1L],
               effect_mean = mean(x$effect), seed_sd = sd(x$effect),
               post_sd = median(x$post_sd),
               seed_sd_over_post_sd = sd(x$effect) / median(x$post_sd),
               ess_median = median(x$ess), secs = mean(x$secs))
  }))
print(out[order(out$family, out$augment, out$gate), ], digits = 3, row.names = FALSE)
