# Does `sparsity = "gibbs"` lose real predictors at moderate p, and is zeta = 1
# the wrong default here?
#
# WHAT IS BEING TESTED
#
# `_dev/gibbs-prior-smoke.R` found the Gibbs prior recovering 4.5 of 5 real
# predictors at p = 20 where DART recovered 5.0, on six replicates, while
# cutting false positives from 1.00 to 0.67. Six replicates cannot tell that
# from noise, and the two readings imply different actions.
#
#   Noise. Three of six replicates happened to drop a predictor and the true
#   recovery is at or near 5.0. Then zeta = 1 stands as the default, the
#   prototype's only cost is the 18-28% the paper reports, and nothing in the
#   documentation has to warn about anything.
#
#   Real. The Gibbs prior over-penalizes as p grows, trading recall for the
#   precision it buys. Linero and Du's own Table 1 is consistent with this: their
#   recall is 1.00 for both priors at P = 7 and slips to 0.94 against DART's 0.96
#   by P = 50, sigma = 5. If it reproduces here it is a property rather than a
#   bug, but it is one that has to be either priced into a lower default zeta or
#   stated plainly in `?bartisan_control`.
#
# The second reading names a mechanism that is directly testable rather than
# only observable: zeta *is* the penalty on model size, so if over-penalization
# is the cause then lowering it should recover the lost predictor. If recovery
# is flat in zeta, the loss is coming from somewhere else and the diagnosis is
# wrong.
#
# WHAT THE OUTCOME MEASUREMENT IS
#
# Friedman at n = 400 with 5 real predictors, crossed over p in {7, 20, 50} and
# over five arms: DART, and the Gibbs prior at zeta in {0, 0.5, 1, 2}. Forty
# replicates a cell, 20 trees, 200 warmup and 400 draws. Selection is the median
# probability model, `prop_used >= 0.5`, which is the rule `variable_importance()`
# reports and the rule `_dev/varsel-check.R` measured.
#
# Two columns, and neither is decisive alone because the whole question is a
# trade:
#
#   recovered  how many of the 5 real predictors are selected. The column the
#              question is about.
#   noise      how many predictors with no signal are selected. The column that
#              says whether any recall the prior gives back was paid for.
#
# Read them together. A zeta that restores recovery by also restoring DART's
# false-positive rate has not found anything; the arm worth having is one whose
# recovery matches DART at noise below DART's.
#
# WHAT EACH OUTCOME WOULD MEAN
#
# Recovery at zeta = 1 within a replicate or two of DART at every p: the smoke
# test's 4.5 was noise, the default stands, nothing to document.
#
# Recovery at zeta = 1 below DART and rising as zeta falls: over-penalization
# confirmed and the mechanism confirmed with it. Then the question is whether
# some lower zeta holds the precision gain, and the noise column answers it. A
# default change would be argued from this table.
#
# Recovery at zeta = 1 below DART and *flat* in zeta: the diagnosis is wrong,
# the loss is not the model-size penalty, and the prototype has a defect to find
# before any of it is reported as a property of the prior.
#
# Recovery falling with p under every arm including DART: the design is too hard
# at p = 50 and that column says nothing about either prior. The p = 7 and p = 20
# columns still stand.

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

SMOKE <- nzchar(Sys.getenv("GIBBS_RECALL_SMOKE"))

REPS <- if (SMOKE) 2L else 40L
PS <- if (SMOKE) c(20L) else c(7L, 20L, 50L)
N <- 400L
REAL <- 5L
OUT <- "_dev/gibbs-recall.rds"

ARMS <- list(
  list(name = "dart", sparsity = TRUE, zeta = 1),
  list(name = "gibbs z=0", sparsity = "gibbs", zeta = 0),
  list(name = "gibbs z=0.5", sparsity = "gibbs", zeta = 0.5),
  list(name = "gibbs z=1", sparsity = "gibbs", zeta = 1),
  list(name = "gibbs z=2", sparsity = "gibbs", zeta = 2)
)

sim <- function(rep, p) {
  set.seed(9000L + 131L * rep + p)
  X <- matrix(stats::runif(N * p), N, p,
              dimnames = list(NULL, paste0("x", seq_len(p))))
  truth <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 +
    10 * X[, 4] + 5 * X[, 5]
  data.frame(y = truth + stats::rnorm(N, 0, 3), X)
}

score <- function(fit) {
  imp <- variable_importance(fit)
  used <- imp[["prop_used"]]
  names(used) <- imp[["variable"]]
  picked <- names(used)[used >= 0.5]
  real <- paste0("x", seq_len(REAL))
  c(recovered = sum(real %in% picked), noise = sum(!picked %in% real))
}

cells <- expand.grid(arm = seq_along(ARMS), p = PS, KEEP.OUT.ATTRS = FALSE)
total <- nrow(cells) * REPS

pr <- prog_init(total = total,
                title = "Gibbs prior: recall against zeta and p",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()

for (i in seq_len(nrow(cells))) {
  arm <- ARMS[[cells[["arm"]][i]]]
  p <- cells[["p"]][i]

  ctl <- bartisan_control(num_trees = 20L, num_burn = 200L, num_draws = 400L,
                          sparsity = arm[["sparsity"]], zeta = arm[["zeta"]])

  for (rep in seq_len(REPS)) {
    d <- sim(rep, p)

    out <- prog_do(pr, i = rep,
                   f = function(r) {
                     set.seed(50L * r + p)
                     score(bartisan(y ~ ., d, family = gaussian(), control = ctl))
                   },
                   label = sprintf("%s, p=%d, rep %d", arm[["name"]], p, rep))

    if (!is.null(out)) {
      rows[[length(rows) + 1L]] <- data.frame(arm = arm[["name"]], p = p,
                                              rep = rep, t(out))
    }
  }

  # After every cell, so a killed run still leaves every finished cell behind.
  saveRDS(list(res = do.call(rbind, rows), complete = FALSE,
               done = length(rows), total = total), OUT)
}

res <- do.call(rbind, rows)
saveRDS(list(res = res, complete = TRUE, done = nrow(res), total = total), OUT)

on.exit()
prog_end(pr, "done")

cat("\n--- mean over", REPS, "replicates ---\n")
agg <- aggregate(cbind(recovered, noise) ~ arm + p, res, mean)
print(agg[order(agg[["p"]], agg[["arm"]]), ], row.names = FALSE)

cat("\n--- recovery against DART, paired within replicate ---\n")
for (p in PS) {
  base <- res[res[["arm"]] == "dart" & res[["p"]] == p, ]
  for (a in vapply(ARMS, `[[`, character(1L), "name")) {
    if (a == "dart") next
    w <- res[res[["arm"]] == a & res[["p"]] == p, ]
    m <- merge(base, w, by = "rep", suffixes = c(".d", ".g"))
    dr <- m[["recovered.g"]] - m[["recovered.d"]]
    dn <- m[["noise.g"]] - m[["noise.d"]]
    se <- stats::sd(dr) / sqrt(length(dr))
    cat(sprintf("p=%-3d %-13s recovered %+5.2f (se %.2f, t %+5.2f)   noise %+5.2f\n",
                p, a, mean(dr), se,
                if (se > 0) mean(dr) / se else 0, mean(dn)))
  }
}
