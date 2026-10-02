# Does recording a ph() draw in the centered chart fix what diagnose() reports,
# and change nothing else?
#
# What is tested. `ph()` models lambda(t | x) = lambda_0(t) exp(eta(x)), and
# multiplying every bin hazard by e^c while subtracting c from eta leaves the
# likelihood unchanged, so the level of the predictor and the level of the
# baseline are identified only jointly. `PHFamily` in `src/family.cpp` records
# each draw in the sampler's own chart, where that level drifts. On the default
# `ph()` fit to `rhc` (the fit in `vignettes/survival.Rmd`), `diagnose()`
# reports R-hat about 1.6 and bulk ESS 3 for every `aux.lambda*` row and for the
# predictor's average, and advises more chains and draws, while
# log(lambda_b) + mean(eta) has bulk ESS 450 to 790 and the survival
# probabilities 660 to 810 (`_dev/TASKS.md`, To Do). The claim: recording each
# draw where the predictor has mean zero over the fitted sample, with the bin
# hazards scaled by exp(mean) and their prior's rate by exp(-mean), through the
# `report_shift()` and `aux_values_shifted()` hooks that `ordinal()` already
# uses, is a change of chart and not of model. The sampler's path is the same
# draw for draw, so every identified quantity is unchanged, and the hazard rows
# become identified quantities that mix like the survival probabilities.
#
# What is measured. One fit per arm, the same call and seed: `before` loads the
# installed package, which was built from the unpatched source, and `after`
# loads the patched source with `pkgload::load_all()`. The columns to read, in
# the comparison printed by the `after` run:
#
#   s_gap        largest difference between the arms in any draw of S(t | x)
#                at t = 30 and 180 days, over 50 patients
#   avg_gap      the same for the per-draw average of S(180 | x) over all 1500
#   loglik_gap   largest difference in the recorded log likelihood per draw
#   level        largest absolute per-draw mean of the recorded predictor
#   rhat, ess_bulk for aux.lambda1, aux.lambda6, aux.lambda12,
#                aux.lambda_rate, the predictor's average and worst 5%, and
#                loglik, from diagnose(), per arm
#   ess_ident    in the before arm, bulk ESS of log(lambda_b) + mean(eta), the
#                identified quantity the after arm's lambda rows should match
#
# What each outcome means. If s_gap, avg_gap and loglik_gap are at rounding
# (below 1e-8), the level is zero to rounding, and the after arm's lambda rows
# have R-hat near 1 with a bulk ESS of the order of ess_ident while the before
# arm's sit near 3, then the misgrading was the chart alone, the patch removes
# it, and nothing a user computes from the fit has moved. If any gap is above
# rounding, the shift is not a pure change of chart (the trees, the predictor
# and the hazards were not moved together), and the patch is wrong however the
# table looks. If the lambda rows still mix poorly after, something besides the
# level drives them and a chart change is not the fix. The loglik row does not
# depend on the chart, so it must read the same in both arms; if it is poor in
# both, that is a separate question for its own measurement, not this one.
#
# Run with: PH_LEVEL_TAG=before Rscript _dev/ph-level-shift.R
#           PH_LEVEL_TAG=after  Rscript _dev/ph-level-shift.R
# Writes:   _dev/ph-level-shift-<tag>.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

TAG <- Sys.getenv("PH_LEVEL_TAG", "run")
OUT <- sprintf("_dev/ph-level-shift-%s.rds", TAG)

if (identical(TAG, "after")) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(bartisan)
}

suppressPackageStartupMessages(library(survival))
data("rhc", package = "bartisan")

pr <- prog_init(total = 4L, title = sprintf("ph() level shift: %s", TAG),
                unit = "step", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

# The vignette's fit: default settings, one chain.
set.seed(2026)
fit <- bartisan(Surv(days, death) ~ rhc + age + sex + race + edu + aps +
                  meanbp + resp + hema + pafi + paco2 + crea + surv2m + card,
                data = rhc, family = ph())
prog_tick(pr, label = "fit")

table <- diagnose(fit)[["table"]]
prog_tick(pr, label = "diagnose")

eta <- predict(fit, type = "link", draws = TRUE)
level <- rowMeans(eta)
lambda <- fit[["aux"]][, grep("^lambda[0-9]+$", colnames(fit[["aux"]])),
                       drop = FALSE]

# The survival draws for the first 50 patients, and the average over all of
# them at 180 days, which is the quantity the vignette reports.
s_some <- predict(fit, newdata = rhc[1:50, ], type = "survival",
                  times = c(30, 180), draws = TRUE)
s_avg <- rowMeans(predict(fit, type = "survival", times = 180,
                          draws = TRUE)[, , 1L])
prog_tick(pr, label = "predict")

ess_ident <- vapply(c(1L, 6L, ncol(lambda)), function(b) {
  posterior::ess_bulk(log(lambda[, b]) + level)
}, numeric(1L))
names(ess_ident) <- paste0("lambda", c(1L, 6L, ncol(lambda)))

res <- list(table = table, loglik = fit[["loglik"]], level = level,
            lambda = lambda, s_some = s_some, s_avg = s_avg,
            ess_ident = ess_ident)
saveRDS(res, OUT)
prog_tick(pr, label = "save")

rows <- c("aux.lambda1", "aux.lambda6", sprintf("aux.lambda%d", ncol(lambda)),
          "aux.lambda_rate", "eta.eta (average over observations)",
          "eta.eta (worst 5% of observations)", "loglik")

cat(sprintf("\n%s: %s\n", TAG, OUT))
print(table[table$quantity %in% rows, c("quantity", "rhat", "ess_bulk")],
      row.names = FALSE)
cat(sprintf("largest |level|: %.3g\n", max(abs(level))))
cat("bulk ESS of log(lambda_b) + level:\n")
print(round(ess_ident))

before_file <- "_dev/ph-level-shift-before.rds"

if (identical(TAG, "after") && file.exists(before_file)) {
  before <- readRDS(before_file)

  cat("\nAgainst the before arm:\n")
  cat(sprintf("s_gap      %.3g\n", max(abs(res$s_some - before$s_some))))
  cat(sprintf("avg_gap    %.3g\n", max(abs(res$s_avg - before$s_avg))))
  cat(sprintf("loglik_gap %.3g\n", max(abs(res$loglik - before$loglik))))

  both <- merge(before$table[before$table$quantity %in% rows,
                             c("quantity", "rhat", "ess_bulk")],
                res$table[res$table$quantity %in% rows,
                          c("quantity", "rhat", "ess_bulk")],
                by = "quantity", suffixes = c("_before", "_after"),
                sort = FALSE)
  print(both, row.names = FALSE, digits = 3)
}

on.exit()
prog_end(pr, "done")
