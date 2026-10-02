# Weights or misspecification: which stops the contrast mixing under hard rules?
#
# What is tested. The bulk effective sample size of a g-computation contrast on
# a binary covariate, over 500 retained draws of a binomial logit fit under
# hard rules, was 5 in `_dev/weights-mixing.R` and about 400 at every sample
# size from 500 to 4000 in `_dev/contrast-mixing.R`. Two things differ between
# those designs: the first fits data whose generator carries a grouping factor
# the model does not include (eight levels, standard deviation 0.6), and the
# first carries frequency weights of 1 to 3. The claim being tested is which of
# the two is responsible, with no prediction between them; the third
# possibility is that only the pair does it.
#
# What is measured. Binomial logit, n = 1000, the plain structure
# `y ~ z + x1 + x2 + x3`, 50 trees, 300 warmup and 500 kept draws, one chain,
# the default augmentation. An omitted group effect present or absent crossed
# with frequency weights present or absent, under hard and smoothstep rules,
# three data-and-fit seeds each. Per fit the bulk effective sample size of the
# contrast, of the predictor averaged over observations and of the log
# likelihood, with the effect and its posterior standard deviation.
#
# What each outcome would mean. A collapse confined to the cells with weights
# puts the fault in how a weight enters the hard-rule leaf target, which is a
# defect to read in the engine. One confined to the cells with the omitted
# group effect makes it a caution about a misspecified mean under hard rules,
# which belongs in `vignette("diagnostics")`. One confined to the cell with
# both makes it an interaction, and the caution has to name both halves.
#
# Run with: Rscript _dev/contrast-mixing-cause.R
# Writes:   _dev/contrast-mixing-cause.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

simulate <- function(seed, group) {
  n <- 1000L
  set.seed(seed)
  d <- as.data.frame(matrix(runif(n * 3L), n, 3L))
  names(d) <- paste0("x", 1:3)
  d$z <- rbinom(n, 1L, 0.5)
  d$g <- factor(sample(letters[1:8], n, TRUE))
  d$w <- sample(1:3, n, TRUE)
  u <- setNames(rnorm(8L, 0, 0.6), letters[1:8])
  f <- 1.2 * d$x1 - 1.2 * d$x2
  lp <- f - mean(f) + 0.8 * d$z + if (group) u[as.character(d$g)] else 0
  d$y <- rbinom(n, 1L, plogis(lp))
  d
}

grid <- expand.grid(group = c(TRUE, FALSE), weights = c(TRUE, FALSE),
                    gate = c("hard", "smoothstep"), seed = 1:3,
                    stringsAsFactors = FALSE)

pr <- prog_init(total = nrow(grid), title = "Contrast mixing: weights or misspecification",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

ess1 <- function(x) posterior::ess_bulk(matrix(x, ncol = 1L))

rows <- list()
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]
  d <- simulate(300L + g$seed, g$group)
  d1 <- d
  d1$z <- 1L
  d0 <- d
  d0$z <- 0L

  args <- list(y ~ z + x1 + x2 + x3, data = d, family = binomial(),
               control = bartisan_control(num_trees = 50L, num_burn = 300L,
                                          num_draws = 500L, gate = g$gate,
                                          verbose = FALSE))
  if (g$weights) args[["weights"]] <- d$w

  t0 <- Sys.time()
  set.seed(400L + g$seed)
  fit <- do.call(bartisan, args)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  contrast <- rowMeans(predict(fit, newdata = d1, type = "link", draws = TRUE) -
                         predict(fit, newdata = d0, type = "link", draws = TRUE))

  rows[[k]] <- data.frame(
    group = g$group, weights = g$weights, gate = g$gate, seed = g$seed,
    effect = mean(contrast), post_sd = sd(contrast),
    ess_contrast = ess1(contrast),
    ess_eta = ess1(rowMeans(fit[["eta"]][[1L]])),
    ess_loglik = ess1(fit[["loglik"]]), secs = secs)
  prog_tick(pr, label = sprintf("group %s / weights %s / %s / seed %d",
                                g$group, g$weights, g$gate, g$seed), secs = secs)
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid)),
          "_dev/contrast-mixing-cause.rds")
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
out <- do.call(rbind, lapply(split(res, list(res$group, res$weights, res$gate)),
  function(x) {
    if (!nrow(x)) return(NULL)
    data.frame(omitted_group = x$group[1L], weights = x$weights[1L],
               gate = x$gate[1L],
               ess_contrast = round(median(x$ess_contrast)),
               ess_eta = round(median(x$ess_eta)),
               ess_loglik = round(median(x$ess_loglik)),
               effect = round(mean(x$effect), 2),
               post_sd = round(median(x$post_sd), 3))
  }))
print(out[order(out$gate, out$omitted_group, out$weights), ], row.names = FALSE)
