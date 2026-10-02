# How fast does a treatment contrast stop mixing as the sample grows?
#
# What is tested. `_dev/weights-mixing.R` found that the bulk effective sample
# size of a g-computation contrast on a binary covariate, over 500 retained
# draws of a binomial logit fit, is 5 under hard rules and 338 under soft once
# the information is doubled by frequency weights, where the same fit without
# weights gives 362 and 457. The record's entry "Mixing against sample size,
# against Tan et al. (2026)" measures a different statistic, the held-out root
# mean squared error of each draw, and finds soft rules ahead of hard by only
# about 1.5 times. The claim: the contrast is a far more sensitive quantity
# than that one, because it moves only through trees that split on the
# treatment, which is a discrete feature of the forest under hard rules and a
# smooth one under soft; so its effective sample size collapses with the sample
# size under hard rules while quantities averaged over the whole fit do not.
#
# What is measured. One generator, `lp = f(x) + 0.8 z` with f centered and
# three uniform predictors, a binomial logit response, the plain structure
# `y ~ z + x1 + x2 + x3`, 50 trees, 300 warmup and 500 kept draws, one chain,
# the default augmentation. Sample sizes 500, 1000, 2000 and 4000, both gates,
# three data-and-fit seeds. Per fit, the bulk effective sample size over the
# 500 retained draws of:
#
#   ess_contrast   the mean over units of link(z = 1) - link(z = 0)
#   ess_eta        the additive predictor averaged over observations
#   ess_loglik     the reported log likelihood
#
# with the share of splitting rules that fall on the treatment alongside.
#
# What each outcome would mean. A contrast whose effective sample size
# collapses with n under hard rules while the averaged predictor's and the log
# likelihood's hold up says the slow quantity is the estimand and not the fit,
# which is a caution for `vignette("diagnostics")` and for `?estimate_effect`,
# and an argument for the soft default beyond accuracy. All three falling
# together would make this the general scaling the record already carries, and
# nothing new.
#
# Run with: Rscript _dev/contrast-mixing.R
# Writes:   _dev/contrast-mixing.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

simulate <- function(n, seed) {
  set.seed(seed)
  d <- as.data.frame(matrix(runif(n * 3L), n, 3L))
  names(d) <- paste0("x", 1:3)
  d$z <- rbinom(n, 1L, 0.5)
  f <- 1.2 * d$x1 - 1.2 * d$x2
  lp <- f - mean(f) + 0.8 * d$z
  d$y <- rbinom(n, 1L, plogis(lp))
  d
}

grid <- expand.grid(n = c(500L, 1000L, 2000L, 4000L),
                    gate = c("hard", "smoothstep"), seed = 1:3,
                    stringsAsFactors = FALSE)

pr <- prog_init(total = nrow(grid), title = "Contrast mixing against sample size",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

ess1 <- function(x) posterior::ess_bulk(matrix(x, ncol = 1L))

rows <- list()
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]
  d <- simulate(g$n, 300L + g$seed)
  d1 <- d
  d1$z <- 1L
  d0 <- d
  d0$z <- 0L

  t0 <- Sys.time()
  set.seed(400L + g$seed)
  fit <- bartisan(y ~ z + x1 + x2 + x3, data = d, family = binomial(),
                  control = bartisan_control(num_trees = 50L, num_burn = 300L,
                                             num_draws = 500L, gate = g$gate,
                                             verbose = FALSE))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  contrast <- rowMeans(predict(fit, newdata = d1, type = "link", draws = TRUE) -
                         predict(fit, newdata = d0, type = "link", draws = TRUE))
  counts <- fit[["counts"]][[1L]]

  rows[[k]] <- data.frame(
    n = g$n, gate = g$gate, seed = g$seed,
    effect = mean(contrast), post_sd = sd(contrast),
    ess_contrast = ess1(contrast),
    ess_eta = ess1(rowMeans(fit[["eta"]][[1L]])),
    ess_loglik = ess1(fit[["loglik"]]),
    z_share = mean(counts[, "z"] / pmax(rowSums(counts), 1)),
    secs = secs)
  prog_tick(pr, label = sprintf("n %d / %s / seed %d", g$n, g$gate, g$seed),
            secs = secs)
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid)),
          "_dev/contrast-mixing.rds")
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
out <- do.call(rbind, lapply(split(res, list(res$n, res$gate)), function(x) {
  if (!nrow(x)) return(NULL)
  data.frame(n = x$n[1L], gate = x$gate[1L],
             ess_contrast = round(median(x$ess_contrast)),
             ess_eta = round(median(x$ess_eta)),
             ess_loglik = round(median(x$ess_loglik)),
             z_share = median(x$z_share), effect = mean(x$effect),
             secs = round(mean(x$secs), 1))
}))
print(out[order(out$gate, out$n), ], digits = 3, row.names = FALSE)
