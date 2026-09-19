# Does Tan, Ronen, Saarinen and Yu's mixing trend hold for this sampler, and
# how far apart do the reported and internal rows of `diagnose()` run?
#
# WHAT IS BEING TESTED
#
# Two claims, one borrowed and one local.
#
# Borrowed. Tan, Ronen, Saarinen and Yu (2026, arXiv:2406.19958) prove that the
# BART sampler's hitting time for a high posterior density set grows with the
# sample size, because the posterior over tree structures is multimodal, and
# measure the consequence empirically as the Gelman-Rubin statistic at a fixed
# iteration budget: R-hat rises with n "for the in-practice BART sampler and
# across all datasets tested". Ronen, Saarinen, Tan, Duncan and Yu (2022,
# arXiv:2210.09352) give the matching exponential mixing-time lower bound and
# conclude by recommending more chains as n grows.
#
# Their analysis is of the Chipman sampler: hard decision rules, conjugate
# Gaussian leaf updates, one response family. **This package is none of those.**
# Rules are soft by default, the kernel is Linero's Laplace-approximation
# reversible jump, and the families include augmented ones whose leaf target is
# quadratic for reasons the original sampler had nothing to say about. Their
# paper reports the trend as robust to the feature-selection prior, to
# initialization and to burn-in length, but it says nothing about soft rules,
# and `_dev/TASKS.md` measures soft rules as reaching their fit three to four
# times faster than hard ones. So whether the trend transfers here is open.
#
# Local. The grading added to `diagnosis_checks()` reads R-hat and effective
# sample size on the quantities a fit reports and gives the splitting rules a
# check of their own, on the argument that a sum of trees reaches one function
# through many partitions and the reported quantities are integrals over that
# structure. That argument is about the parameterization, not a measurement.
# What is unmeasured is how far apart the two actually run, whether the gap is
# stable across families and sample sizes, and -- the case that would send this
# back to the code -- whether keying the separation on the `splits.` prefix
# misses other rows that behave the same way.
#
# `ordinal()` is in the design for exactly that reason. `_dev/TASKS.md` records
# its cutpoints as giving "a median effective sample size of only 7 to 21 per 100
# draws" while "the regression function itself mixes fine", so it is a family
# where a *reported* parameter mixes badly. If the gap is about reported against
# internal, ordinal's cutpoints should sit with the reported rows and mix badly
# anyway, and the advice needs to name them.
#
# WHAT THE OUTCOME MEASUREMENT IS
#
# Three families: `gaussian()`, `binomial("probit")`, and `ordinal()` on a
# 20-category response, which is fast per sweep here because the cutpoint block
# is O(n) but has 19 tightly coupled thresholds to move. Hard and soft rules.
# n of 500, 2000 and 8000. Eight chains of 1000 draws after 200 warmup, which is
# the shape Tan et al. use, scaled down. Eight replicates. `num_trees` is left at
# the default throughout and is not an arm: it is a component of the model rather
# than a sampler setting, and the record already establishes the default as
# sufficient or better.
#
# The primary column is **R-hat on a scalar functional of the fit**, computed the
# way they compute it: the held-out root mean squared error of each draw against
# the truth, one number per draw, reshaped into chains and passed through this
# package's own `diagnosis_stats()` so the statistic is the one `diagnose()`
# reports. Reading it against n is the borrowed claim.
#
# The secondary columns are R-hat and bulk effective sample size for every row of
# `diagnose()$table`, from which the reported-against-internal gap follows. The
# gap, not its two halves, is the quantity of interest, and it is paired within a
# fit so it carries no between-replicate noise.
#
# Coverage and mean interval width come along for interpretation. For `ordinal()`
# the additive predictor is identified only up to a shift, so its coverage is
# computed after centering the draws and the truth alike; the other two are
# compared directly.
#
# WHAT EACH OUTCOME WOULD MEAN
#
#   Functional R-hat rises with n under both rule types. Their finding transfers
#   and the sampler here is no exception. `?bartisan_control` should say that
#   larger n wants more chains, with the citation, and the chains recommendation
#   gains support that does not depend on this package's own measurements.
#
#   It rises for hard rules and not for soft. Soft rules are a mitigation their
#   paper does not cover, which is a finding about this package's default worth
#   stating plainly, and it argues for the default on mixing grounds and not only
#   on the accuracy grounds already recorded.
#
#   It does not rise for either. Something about this kernel or these families
#   breaks the mechanism, and the next question is which -- the Laplace
#   reversible jump, or the augmentations that make the leaf target quadratic.
#   That would be a claim against a published result on a neighboring sampler,
#   so it would need a second design before being written down as one.
#
#   The reported-against-internal gap is large and stable across families and n.
#   The grading is justified and `?diagnose` can quantify it instead of arguing
#   from the parameterization alone.
#
#   The gap is erratic, or ordinal's cutpoints mix as badly as its splits do.
#   Then "reported against internal" is the wrong cut. Keying on the `splits.`
#   prefix would be too narrow, the advice has to name the cutpoints, and the
#   change just made to `diagnosis_checks()` needs revisiting rather than
#   documenting.

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

SMOKE <- nzchar(Sys.getenv("MIXING_SMOKE"))

N_TEST <- 500L
P <- 10L
K_ORD <- 20L
LEVEL <- 0.95

CHAINS <- 8L
BURN <- 200L
DRAWS <- 1000L
NS <- c(500L, 2000L, 8000L)
REPS <- 8L
OUT <- "_dev/mixing-against-n.rds"

if (SMOKE) {
  CHAINS <- 4L
  BURN <- 20L
  DRAWS <- 60L
  NS <- c(500L)
  REPS <- 1L
  OUT <- "_dev/mixing-against-n-smoke.rds"
}

FAMILIES <- c("gaussian", "probit", "ordinal")
GATES <- c(hard = "hard", soft = "smoothstep")

# One latent function for all three families, so that the only thing changing
# across them is the likelihood. Friedman on the first five of ten predictors,
# scaled to a standard deviation of one so the probit and ordinal signal-to-noise
# ratios are comparable to the Gaussian one.
sim <- function(n, rep, family) {
  set.seed(52000 + 211 * rep + match(family, FAMILIES))
  total <- n + N_TEST
  X <- matrix(stats::runif(total * P), total, P,
              dimnames = list(NULL, paste0("x", seq_len(P))))
  f <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 +
    10 * X[, 4] + 5 * X[, 5]
  eta <- (f - mean(f)) / stats::sd(f)

  y <- switch(family,
              gaussian = eta + stats::rnorm(total, sd = 0.5),
              probit = stats::rbinom(total, 1L,
                                     stats::pnorm(eta)),
              ordinal = {
                latent <- eta + stats::rlogis(total)
                cuts <- stats::quantile(latent,
                                        seq_len(K_ORD - 1L) / K_ORD)
                factor(findInterval(latent, cuts), levels = 0:(K_ORD - 1L),
                       ordered = TRUE)
              })

  dat <- cbind(data.frame(y = y), as.data.frame(X))
  train <- seq_len(n)

  list(train = dat[train, ], test = dat[-train, ], truth = eta[-train])
}

fit_family <- function(family) {
  switch(family,
         gaussian = stats::gaussian(),
         probit = stats::binomial("probit"),
         ordinal = ordinal())
}

# Tan et al.'s statistic: one scalar per draw, then R-hat across chains. Draws
# are stacked chain by chain, which is the layout `diagnose()`'s own
# `as_chains()` assumes, so the same reshape applies here.
functional_stats <- function(eta_draws, truth, chains) {
  per_draw <- sqrt(rowMeans((eta_draws - rep(truth, each = nrow(eta_draws)))^2))
  per <- length(per_draw) %/% chains
  wide <- matrix(per_draw[seq_len(per * chains)], nrow = per, ncol = chains)

  bartisan:::diagnosis_stats(wide)
}

run <- function(family, rule, n, rep) {
  d <- sim(n, rep, family)

  elapsed <- system.time(
    fit <- bartisan(y ~ ., data = d$train, family = fit_family(family),
                    chains = CHAINS, num_burn = BURN, num_draws = DRAWS,
                    num_thin = 1L, gate = GATES[[rule]])
  )[["elapsed"]]

  eta <- predict(fit, newdata = d$test, type = "link", draws = TRUE)

  if (is.list(eta)) {
    eta <- eta[[1L]]
  }

  truth <- d$truth

  # The ordinal additive predictor is identified up to a shift, so it and the
  # truth are both centered before being compared. The other two are on the
  # scale they were generated on.
  if (identical(family, "ordinal")) {
    eta <- eta - rowMeans(eta)
    truth <- truth - mean(truth)
  }

  fs <- functional_stats(eta, truth, CHAINS)

  lo <- apply(eta, 2L, stats::quantile, probs = (1 - LEVEL) / 2)
  hi <- apply(eta, 2L, stats::quantile, probs = 1 - (1 - LEVEL) / 2)
  mid <- colMeans(eta)

  tb <- diagnose(fit)[["table"]]

  # Every row, long, so the reported-against-internal gap can be computed any
  # way later rather than being fixed here.
  rows <- data.frame(
    family = family, rule = rule, n = n, rep = rep,
    quantity = tb[["quantity"]], rhat = tb[["rhat"]],
    ess_bulk = tb[["ess_bulk"]], ess_tail = tb[["ess_tail"]])

  summary <- data.frame(
    family = family, rule = rule, n = n, rep = rep,
    rhat_functional = as.numeric(fs[["rhat"]]),
    ess_functional = as.numeric(fs[["ess_bulk"]]),
    coverage = mean(truth >= lo & truth <= hi),
    width = mean(hi - lo),
    abs_bias = mean(abs(mid - truth)),
    rmse = sqrt(mean((mid - truth)^2)),
    seconds = elapsed)

  list(summary = summary, rows = rows)
}

cells <- expand.grid(rep = seq_len(REPS), n = NS, rule = names(GATES),
                     family = FAMILIES, stringsAsFactors = FALSE)
cells <- cells[order(cells$rep, cells$family, cells$rule, cells$n), ]

pr <- prog_init(total = nrow(cells), title = "Mixing against sample size",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

summaries <- list()
all_rows <- list()

for (i in seq_len(nrow(cells))) {
  cell <- cells[i, ]
  label <- sprintf("rep %d: %s %s n=%d", cell$rep, cell$family, cell$rule,
                   cell$n)

  out <- prog_do(pr, i, function(.) {
    run(cell$family, cell$rule, cell$n, cell$rep)
  }, label = label)

  if (!is.null(out)) {
    summaries[[length(summaries) + 1L]] <- out$summary
    all_rows[[length(all_rows) + 1L]] <- out$rows
  }

  # After every fit, so a killed run leaves every finished one behind.
  saveRDS(list(res = do.call(rbind, summaries),
               rows = do.call(rbind, all_rows), complete = FALSE,
               done = length(summaries), total = nrow(cells)), OUT)
}

saveRDS(list(res = do.call(rbind, summaries),
             rows = do.call(rbind, all_rows), complete = TRUE,
             done = length(summaries), total = nrow(cells)), OUT)

on.exit()
prog_end(pr, "done")

cat("done:", length(summaries), "of", nrow(cells), "\n")
