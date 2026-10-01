# Known-truth recovery across the family-by-structure matrix.
#
# What is tested. `_dev/AUDITING.md` sec. 3 fills the binomial row of the family
# by structure matrix with recovery of a known effect and leaves every other
# row at the two identities. The bug of 2026-10-01, a `bcf()` fit under `dpm()`
# whose recorded level was the unidentified ridge coordinate, would have failed
# a recovery check on the level and passed every identity then in the suite.
# The claim: every family recovers a link-scale effect of 0.8 and, where its
# predictor's level is identified, a known level, under each of five
# structures, to the same degree across structures. A hook dropped by a wrapper
# shows up as one structure's cell being off where the family's other cells are
# not; shrinkage that is the same in every cell of a row is the prior.
#
# What is measured. n = 1000, three uniform predictors, a binary z, a grouping
# factor with eight levels and random intercepts of sd 0.6, and
# lp = f(x) + 0.8 z + u_g with f centered. Each family's response is drawn from
# lp on its own scale with a family-specific offset c (-0.5 for the count and
# positive families so that means stay moderate, 0 otherwise). The survival
# responses use each family's own error, so the predictor is the location of
# log T in every case. Structures: plain (z a covariate), vc(z) with the fixed
# coding, vc(z, center = "estimate") (the drawn coding bcf() uses, where the
# family's leaf target allows it), random intercepts, and vc(z) with random
# intercepts. 50 trees, hard rules, 300 + 500 draws, one chain. Per cell:
#
#   effect   mean over units of link(z = 1) - link(z = 0); the truth is 0.8
#   level    mean over units of the fitted link minus mean(lp + c); NA for the
#            ordinal families, whose chart is centered, and for ph(), whose
#            baseline hazard absorbs it
#   secs     the fit's wall time
#
# What each outcome means. A row whose cells agree within Monte Carlo error
# (effects within about 0.15 of one another, levels within about 0.1) passes:
# the wrappers compose for that family. One cell off from the rest of its row
# by more than that is a dropped hook or a wrong chart in that composition and
# is what to investigate first, as the dpm cell under the drawn coding would
# have read before the fix. A whole row shrunk alike (an effect near 0.7 in
# every cell) is the prior, which the binomial row already shows at 88 to 102%.
#
# Run with: Rscript _dev/recovery-matrix.R          RECOVERY_SMOKE=1 for 4 cells
#           RECOVERY_GATE=smoothstep Rscript _dev/recovery-matrix.R   soft rules
# Writes:   _dev/recovery-matrix.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

SMOKE <- nzchar(Sys.getenv("RECOVERY_SMOKE"))
# The gate: "hard" for the 2026-10-01 run, "smoothstep" for the soft-rule arm,
# which is the package default and a different code path through the leaves.
GATE <- Sys.getenv("RECOVERY_GATE", "hard")
OUT <- {
  if (SMOKE) "_dev/recovery-matrix-smoke.rds"
  else if (identical(GATE, "hard")) "_dev/recovery-matrix.rds"
  else sprintf("_dev/recovery-matrix-%s.rds", GATE)
}

N <- if (SMOKE) 300L else 1000L
control <- bartisan_control(num_trees = if (SMOKE) 10L else 50L,
                            num_burn = if (SMOKE) 50L else 300L,
                            num_draws = if (SMOKE) 50L else 500L,
                            gate = GATE, verbose = FALSE)

set.seed(41)
n <- N
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
# The ordered beta model's own boundary mechanism. The 2026-10-01 hard-rule run
# set the boundaries by a deterministic cut on lp instead, a step the model
# represents with logistic probabilities, and the ordbeta row overshot 0.8 in
# every cell and swung across seeds; `_dev/recovery-followup.R` reran that row
# this way. The change added the `runif()` call below, which moved the random
# stream for every response drawn after it, so the soft-rule arm of the same
# day fit the tweedie, ordinal and survival rows to different data from the
# hard arm; `_dev/ordinal-soft-check.R` has what that looked like. A seed per
# response is the fix to make before the next cross-arm comparison.
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
d$ydpm <- lp + 2 * (rgamma(n, 1.5, 1.5) - 1)

cens <- lp + rnorm(n, 0.6, 0.7)
for (nm in c("weib", "llog", "lnorm", "dpma")) {
  e <- switch(nm,
              weib = log(rexp(n)),
              llog = rlogis(n),
              lnorm = rnorm(n),
              dpma = 1.2 * (rgamma(n, 1.5, 1.5) - 1))
  lt <- lp + 0.5 * e
  d[[paste0("t_", nm)]] <- exp(pmin(lt, cens))
  d[[paste0("e_", nm)]] <- as.integer(lt <= cens)
}
lt <- log((-log(runif(n)) / exp(lp))^(1 / 1.5))
d$t_ph <- exp(pmin(lt, cens))
d$e_ph <- as.integer(lt <= cens)

# `offset` is c, the shift between lp and the family's link-scale location;
# `level = NA` marks a family whose level is not identified.
cases <- list(
  list(name = "gaussian", family = gaussian(), response = "ynorm", offset = 0),
  list(name = "gaussian_ls", family = gaussian_ls(), response = "ynorm", offset = 0),
  list(name = "binomial logit", family = binomial(), response = "ybin", offset = 0),
  list(name = "binomial probit", family = binomial("probit"), response = "yprob", offset = 0),
  list(name = "binomial cloglog", family = binomial("cloglog"), response = "ycll", offset = 0),
  list(name = "poisson", family = poisson(), response = "ycnt", offset = -0.5),
  list(name = "negbin", family = negbin(), response = "ynb", offset = -0.5),
  list(name = "zi_poisson", family = zi_poisson(), response = "yzi", offset = -0.5),
  list(name = "zi_negbin", family = zi_negbin(), response = "yzinb", offset = -0.5),
  list(name = "Gamma", family = Gamma("log"), response = "ygam", offset = -0.5),
  list(name = "Gamma_ls", family = Gamma_ls(), response = "ygam", offset = -0.5),
  list(name = "Beta", family = Beta(), response = "ybeta", offset = 0),
  list(name = "ordbeta", family = ordbeta(), response = "yob", offset = 0),
  list(name = "tweedie", family = tweedie(), response = "ytw", offset = -0.5),
  list(name = "ordinal logit", family = ordinal(), response = "yord", level = NA),
  list(name = "ordinal probit", family = ordinal("probit"), response = "yord", level = NA),
  list(name = "ordinal cloglog", family = ordinal("cloglog"), response = "yord", level = NA),
  list(name = "dpm", family = dpm(), response = "ydpm", offset = 0),
  list(name = "weibull_aft", family = weibull_aft(), response = "cbind(t_weib, e_weib)", offset = 0),
  list(name = "loglogistic_aft", family = loglogistic_aft(), response = "cbind(t_llog, e_llog)", offset = 0),
  list(name = "lognormal_aft", family = lognormal_aft(), response = "cbind(t_lnorm, e_lnorm)", offset = 0),
  list(name = "dpm_aft", family = dpm_aft(), response = "cbind(t_dpma, e_dpma)", offset = 0),
  list(name = "ph", family = ph(), response = "cbind(t_ph, e_ph)", level = NA)
)

structures <- c("plain", "vc", "drawn", "ranef", "vc + ranef")

if (SMOKE) {
  cases <- cases[c(1L, 18L)]
  structures <- c("plain", "drawn")
}

formula_for <- function(response, structure) {
  terms <- switch(structure,
                  plain = c("z", "x1", "x2", "x3"),
                  vc = c("x1", "x2", "x3", "vc(z)"),
                  drawn = c("x1", "x2", "x3", 'vc(z, center = "estimate")'),
                  ranef = c("z", "x1", "x2", "x3", "(1 | g)"),
                  `vc + ranef` = c("x1", "x2", "x3", "vc(z)", "(1 | g)"))

  as.formula(paste(response, "~", paste(terms, collapse = " + ")),
             env = globalenv())
}

first_predictor <- function(p) {
  if (is.list(p)) p <- p[[1L]]
  if (is.matrix(p)) p[, 1L] else p
}

d1 <- d
d1$z <- 1L
d0 <- d
d0$z <- 0L

grid <- expand.grid(case = seq_along(cases), structure = structures,
                    stringsAsFactors = FALSE)

pr <- prog_init(total = nrow(grid),
                title = sprintf("Recovery matrix%s", if (SMOKE) " (smoke)" else ""),
                unit = "cell", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()

for (k in seq_len(nrow(grid))) {
  case <- cases[[grid$case[k]]]
  structure <- grid$structure[k]
  label <- sprintf("%s / %s", case$name, structure)

  form <- formula_for(case$response, structure)
  t0 <- Sys.time()

  fit <- tryCatch(
    bartisan(form, data = d, family = case$family, control = control),
    error = function(e) {
      if (identical(structure, "drawn") &&
          grepl("leaf target is\\s+quadratic", conditionMessage(e))) {
        return(NULL)
      }
      structure(list(message = conditionMessage(e)), class = "cell_error")
    })

  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (is.null(fit)) {
    rows[[k]] <- data.frame(family = case$name, structure = structure,
                            effect = NA_real_, level = NA_real_, secs = secs,
                            note = "drawn coding refused: target not quadratic")
  }
  else if (inherits(fit, "cell_error")) {
    rows[[k]] <- data.frame(family = case$name, structure = structure,
                            effect = NA_real_, level = NA_real_, secs = secs,
                            note = paste("ERROR:", fit$message))
  }
  else {
    l1 <- first_predictor(predict(fit, newdata = d1, type = "link"))
    l0 <- first_predictor(predict(fit, newdata = d0, type = "link"))
    lfit <- first_predictor(predict(fit, type = "link"))

    level <- {
      if (isTRUE(is.na(case$level))) NA_real_
      else mean(lfit) - mean(lp + case$offset)
    }

    rows[[k]] <- data.frame(family = case$name, structure = structure,
                            effect = mean(l1 - l0), level = level, secs = secs,
                            note = "")
  }

  prog_tick(pr, label = label)
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid),
               done = k, total = nrow(grid), truth = 0.8), OUT)
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
rownames(res) <- NULL

wide <- function(col) {
  out <- reshape(res[, c("family", "structure", col)], idvar = "family",
                 timevar = "structure", direction = "wide")
  names(out) <- sub(paste0("^", col, "\\."), "", names(out))
  out
}

cat("\nEffect (truth 0.8):\n")
print(wide("effect"), digits = 2, row.names = FALSE)
cat("\nLevel error (fitted link minus truth):\n")
print(wide("level"), digits = 2, row.names = FALSE)
cat("\nSeconds:\n")
print(wide("secs"), digits = 2, row.names = FALSE)
if (any(nzchar(res$note))) {
  cat("\nNotes:\n")
  print(res[nzchar(res$note), c("family", "structure", "note")], row.names = FALSE)
}
