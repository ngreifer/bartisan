# Is the weighted binomial fit's stall a runaway leaf scale?
#
# What is tested. `_dev/weights-warmup.R` found that a binomial logit fit with
# frequency weights of 1 to 3 on n = 1000 leaves the reported log likelihood
# with a bulk effective sample size of 1 to 5 and a drift of about one standard
# deviation across the retained draws, and that a tenfold longer warmup barely
# moves either; the treatment contrast is inflated to 1.15 and 1.50 against a
# truth of 0.8, worse at the longer warmup. Those runs also printed
# `warn_runaway_scale()`, the package's own check that the leaf scale is
# growing without bound, and the explanation it offers fits: weight on a binary
# response pushes each leaf toward separation, where the likelihood is
# maximized at an infinite predictor, so the drawn leaf scale keeps climbing
# and the predictor follows it. The claim: holding the leaf scale fixed removes
# the drift, restores the effective sample sizes and returns the contrast to
# about 0.8.
#
# What is measured. Binomial logit, n = 1000, `y ~ z + x1 + x2 + x3`, weights
# of 1 to 3, 50 trees, 300 warmup and 1000 kept draws, one chain, both gates,
# two data-and-fit seeds, with `update_sigma_mu` true (the default) and false.
# The unweighted fit is included as the reference. Per fit:
#
#   sigma_mu_first, sigma_mu_last   the leaf scale averaged over the first and
#                                   last tenth of the retained draws
#   ess_loglik, drift               as in `_dev/weights-warmup.R`
#   ess_contrast, effect            the treatment contrast, truth 0.8
#   warned                          whether the runaway-scale check fired
#
# What each outcome would mean. If the fixed scale restores the effective
# sample sizes and the effect while the drawn scale climbs through the retained
# draws, the weights arm of `_dev/recovery-matrix.R` was measuring this
# pathology and not how weights enter the target, the entry for it has to say
# so, and the caution belongs with `weights` in the documentation. If the fixed
# scale changes nothing, the stall is something else and the leaf scale is a
# symptom.
#
# Run with: Rscript _dev/weights-runaway.R
# Writes:   _dev/weights-runaway.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

simulate <- function(seed) {
  n <- 1000L
  set.seed(seed)
  d <- as.data.frame(matrix(runif(n * 3L), n, 3L))
  names(d) <- paste0("x", 1:3)
  d$z <- rbinom(n, 1L, 0.5)
  d$w <- sample(1:3, n, TRUE)
  f <- 1.2 * d$x1 - 1.2 * d$x2
  lp <- f - mean(f) + 0.8 * d$z
  d$y <- rbinom(n, 1L, plogis(lp))
  d
}

grid <- expand.grid(weighted = c(TRUE, FALSE), update_scale = c(TRUE, FALSE),
                    gate = c("hard", "smoothstep"), seed = 1:2,
                    stringsAsFactors = FALSE)

pr <- prog_init(total = nrow(grid), title = "Weighted binomial: the leaf scale",
                unit = "fit", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

ess1 <- function(x) posterior::ess_bulk(matrix(x, ncol = 1L))

rows <- list()
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]
  d <- simulate(300L + g$seed)
  d1 <- d
  d1$z <- 1L
  d0 <- d
  d0$z <- 0L

  args <- list(y ~ z + x1 + x2 + x3, data = d, family = binomial(),
               control = bartisan_control(num_trees = 50L, num_burn = 300L,
                                          num_draws = 1000L, gate = g$gate,
                                          update_sigma_mu = g$update_scale,
                                          verbose = FALSE))
  if (g$weighted) args[["weights"]] <- d$w

  warned <- FALSE
  set.seed(400L + g$seed)
  fit <- withCallingHandlers(do.call(bartisan, args),
                             warning = function(w) {
                               if (grepl("leaf scale", conditionMessage(w),
                                         fixed = TRUE)) {
                                 warned <<- TRUE
                                 invokeRestart("muffleWarning")
                               }
                             })

  ll <- fit[["loglik"]]
  half <- seq_len(length(ll) / 2L)
  sm <- fit[["sigma_mu"]][, 1L]
  tenth <- seq_len(length(sm) %/% 10L)
  contrast <- rowMeans(predict(fit, newdata = d1, type = "link", draws = TRUE) -
                         predict(fit, newdata = d0, type = "link", draws = TRUE))

  rows[[k]] <- data.frame(
    weighted = g$weighted, update_scale = g$update_scale, gate = g$gate,
    seed = g$seed,
    sigma_mu_first = mean(sm[tenth]),
    sigma_mu_last = mean(sm[length(sm) - rev(tenth) + 1L]),
    ess_loglik = ess1(ll), drift = (mean(ll[-half]) - mean(ll[half])) / sd(ll),
    ess_contrast = ess1(contrast), effect = mean(contrast), warned = warned)
  prog_tick(pr, label = sprintf("weighted %s / scale %s / %s / seed %d",
                                g$weighted, g$update_scale, g$gate, g$seed))
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid)),
          "_dev/weights-runaway.rds")
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
out <- do.call(rbind, lapply(split(res, list(res$weighted, res$update_scale, res$gate)),
  function(x) {
    if (!nrow(x)) return(NULL)
    data.frame(wts = x$weighted[1L], drawn_scale = x$update_scale[1L],
               gate = x$gate[1L],
               sm_first = round(median(x$sigma_mu_first), 3),
               sm_last = round(median(x$sigma_mu_last), 3),
               ess_ll = round(median(x$ess_loglik), 1),
               drift = round(median(x$drift), 2),
               ess_contrast = round(median(x$ess_contrast)),
               effect = round(mean(x$effect), 2), warned = any(x$warned))
  }))
print(out[order(out$gate, out$wts, out$drawn_scale), ], row.names = FALSE)
