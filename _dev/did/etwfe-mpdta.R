# Wooldridge's extended TWFE, fit in bartisan by cohort imputation, on mpdta.
#
# What is being tested. Wooldridge (2025, eq. 4.10; 2023, eq. 3.4) writes the
# never-treated mean as a cohort part plus a time part, each with its own
# covariate slopes, and shows (Proposition 5.2) that fitting that on the
# untreated rows and imputing the treated ones (Procedure 4.1) is numerically
# ETWFE. The claim is that the linear covariate terms can be replaced by two
# additive forests, a(cohort, x) and b(period, x), written
# `y ~ cohort + x + vc(one, ~ period + x)` and fit on the untreated rows, and
# that on mpdta, where a linear specification in log population is plausible,
# this reproduces the linear ETWFE's ATT(g, t) within their intervals. It could
# be false if regularization shrinks the cohort level shifts DiD rests on, or if
# the two forests' shared covariate makes the split mix too badly to use. Two
# side claims: on a panel this persistent (level sd 1.51 against innovation sd
# 0.16) a unit random intercept is what brings the width down near the linear
# estimator's clustered width; and a single forest over cohort, period and x,
# which drops the additivity that encodes parallel trends, extrapolates freely
# into the treated cells and should not be trusted.
#
# The measurement. Four bartisan variants (additive + county random intercept;
# additive without it; per-period coefficient forests; one joint forest),
# 4 chains, hard gates. Overall ATT, the 7 ATT(g, t) cells and 4 event-time
# averages, posterior mean and 95% interval, against the linear imputation
# estimate (checked identical to the POLS/ETWFE regression) with county-
# clustered SEs, and did::att_gt(est_method = "reg") for context. The
# additive + random-intercept variant is run at 1000 and at 2500 kept draws per
# chain, and R-hat and bulk ESS are read on the ATT draws at both, since an ESS
# on these designs means nothing at one chain length (DID.md, phase 13c).
#
# What each outcome means. Additive + RE matching linear ETWFE at comparable
# width with ESS rising with length: ETWFE by imputation is straightforward in
# bartisan and the notes prefer it to Category B. Much wider, or ESS falling:
# the forest split is a ridge and imputation with forests needs more work. No-RE
# much wider than RE: the random intercept is load-bearing on persistent panels.
# Joint forest agreeing: additivity is not binding on a 5-period, 3-cohort
# panel; departing: it shows why the additive form is the identifying one.

source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")
suppressMessages({
  library(bartisan)
  library(did)
  library(fixest)
})

out_file <- Sys.getenv("OUT", "_dev/did/etwfe-mpdta.rds")
# SMOKE=1 shrinks every chain to check the script end to end in seconds.
smoke <- nzchar(Sys.getenv("SMOKE"))
future::plan(future::multisession, workers = 4L)

data(mpdta, package = "did")
d <- mpdta
d$cohort <- factor(ifelse(d$first.treat == 0, "never", d$first.treat),
                   levels = c("never", "2004", "2006", "2007"))
d$w <- as.integer(d$first.treat > 0 & d$year >= d$first.treat)
d$period <- d$year
d$periodf <- factor(d$year)
d$one <- 1
d$e <- ifelse(d$w == 1, d$year - d$first.treat, NA)
d$cell <- ifelse(d$w == 1, paste(d$first.treat, d$year, sep = ":"), "none")

ctl <- d[d$w == 0, ]
tr <- d[d$w == 1, ]
cells <- sort(unique(tr$cell))

# ---- The linear reference: Procedure 4.1 and Procedure 5.1 --------------------
lin <- lm(lemp ~ cohort * lpop + periodf * lpop, data = ctl)
tr$te_lin <- tr$lemp - predict(lin, newdata = tr)
imp_gt <- tapply(tr$te_lin, tr$cell, mean)[cells]

# POLS on cohort dummies (5.3): cell dummies and their cohort-centered covariate
# interactions on the treated rows, and the controls on all rows.
xbar <- tapply(d$lpop, d$cohort, mean)
d$xdot <- d$lpop - xbar[as.character(d$cohort)]
for (k in cells) {
  nm <- paste0("c", gsub(":", "_", k))
  d[[nm]] <- as.integer(d$cell == k)
  d[[paste0(nm, "_x")]] <- d[[nm]] * d$xdot
}
cell_terms <- paste0("c", gsub(":", "_", cells))
pols_f <- stats::as.formula(paste(
  "lemp ~", paste(c(cell_terms, paste0(cell_terms, "_x")), collapse = " + "),
  "+ cohort * lpop + periodf * lpop"))
pols <- feols(pols_f, data = d, cluster = ~ countyreal)
pols_gt <- coef(pols)[cell_terms]
names(pols_gt) <- cells
V <- vcov(pols)[cell_terms, cell_terms]
n_gt <- table(tr$cell)[cells]

lin_summary <- function(wts) {
  wts <- wts / sum(wts)
  est <- sum(wts * pols_gt)
  se <- sqrt(drop(t(wts) %*% V %*% wts))
  c(estimate = est, lower = est - 1.96 * se, upper = est + 1.96 * se)
}

cell_g <- as.integer(sub(":.*", "", cells))
cell_t <- as.integer(sub(".*:", "", cells))
cell_e <- cell_t - cell_g
aggs <- c(list(overall = rep(TRUE, length(cells))),
          setNames(lapply(sort(unique(cell_g)), function(g) cell_g == g),
                   paste0("g", sort(unique(cell_g)))),
          setNames(lapply(sort(unique(cell_e)), function(e) cell_e == e),
                   paste0("e", sort(unique(cell_e)))))
linear <- do.call(rbind, lapply(names(aggs), function(a) {
  wts <- ifelse(aggs[[a]], as.numeric(n_gt), 0)
  data.frame(method = "linear ETWFE", quantity = a, t(lin_summary(wts)))
}))
cells_lin <- data.frame(method = "linear ETWFE", quantity = cells,
                        estimate = pols_gt,
                        lower = pols_gt - 1.96 * sqrt(diag(V)),
                        upper = pols_gt + 1.96 * sqrt(diag(V)))

equiv <- max(abs(pols_gt - imp_gt))

cs <- att_gt(yname = "lemp", tname = "year", idname = "countyreal",
             gname = "first.treat", xformla = ~ lpop, data = mpdta,
             control_group = "nevertreated", est_method = "reg")
cs_simple <- aggte(cs, type = "simple")

# ---- bartisan by imputation ----------------------------------------------------
variants <- list(
  additive_re = list(f = lemp ~ cohort + lpop + vc(one, ~ period + lpop) +
                       (1 | countyreal), trees = c(100L, 50L)),
  additive = list(f = lemp ~ cohort + lpop + vc(one, ~ period + lpop),
                  trees = c(100L, 50L)),
  per_period_re = list(f = lemp ~ cohort + lpop + vc(periodf, ~ lpop) +
                         (1 | countyreal), trees = NULL),
  joint_re = list(f = lemp ~ cohort + period + lpop + (1 | countyreal),
                  trees = 150L))
runs <- list(list(v = "additive_re", draws = 1000L),
             list(v = "additive_re", draws = 2500L),
             list(v = "additive", draws = 1000L),
             list(v = "per_period_re", draws = 1000L),
             list(v = "joint_re", draws = 1000L))

pr <- prog_init(total = length(runs), title = "ETWFE by imputation on mpdta",
                unit = "fit", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list(linear, cells_lin)
diag_rows <- list()

for (r in runs) {
  spec <- variants[[r$v]]
  set.seed(2026)
  t0 <- Sys.time()
  args <- list(formula = spec$f, data = ctl, family = gaussian(), chains = 4L,
               gate = "hard", num_burn = if (smoke) 20L else 500L,
               num_draws = if (smoke) 20L else r$draws,
               verbose = FALSE)
  if (!is.null(spec$trees)) {
    args$num_trees <- spec$trees
  }
  fit <- do.call(bartisan, args)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  p0 <- predict(fit, newdata = tr, draws = TRUE)
  te <- sweep(-p0, 2L, tr$lemp, "+")
  cell_draws <- sapply(cells, function(k) rowMeans(te[, tr$cell == k, drop = FALSE]))

  label <- sprintf("%s, %d draws", r$v, r$draws)
  summ <- function(v) c(estimate = mean(v), lower = unname(quantile(v, .025)),
                        upper = unname(quantile(v, .975)))
  for (a in names(aggs)) {
    wts <- ifelse(aggs[[a]], as.numeric(n_gt), 0)
    v <- drop(cell_draws %*% (wts / sum(wts)))
    rows[[length(rows) + 1L]] <- data.frame(method = label, quantity = a,
                                            t(summ(v)))
  }
  for (k in cells) {
    rows[[length(rows) + 1L]] <- data.frame(method = label, quantity = k,
                                            t(summ(cell_draws[, k])))
  }

  att <- drop(cell_draws %*% (as.numeric(n_gt) / sum(n_gt)))
  per <- length(att) / 4L
  m <- matrix(att, per, 4L)
  diag_rows[[length(diag_rows) + 1L]] <- data.frame(
    method = label, seconds = round(secs),
    rhat = posterior::rhat(m), ess_bulk = posterior::ess_bulk(m),
    ess_tail = posterior::ess_tail(m),
    sigma_unit = if (!is.null(fit[["tau"]])) mean(fit[["tau"]][[1L]][, 1L]) else NA)

  saveRDS(list(att = do.call(rbind, rows), diag = do.call(rbind, diag_rows),
               equiv = equiv, cs = cs_simple, complete = FALSE), out_file)
  prog_tick(pr, label = label)
}

res <- list(att = do.call(rbind, rows), diag = do.call(rbind, diag_rows),
            equiv = equiv, cs = cs_simple, complete = TRUE)
saveRDS(res, out_file)
on.exit()
prog_end(pr, "done")

cat(sprintf("Imputation against POLS/ETWFE, max |difference| over cells: %.2e\n",
            equiv))
cat(sprintf("did att_gt reg, never-treated, simple ATT: %.4f (se %.4f)\n\n",
            cs_simple$overall.att, cs_simple$overall.se))
print(res$diag, row.names = FALSE, digits = 3)
cat("\n")
a <- res$att
a$width <- a$upper - a$lower
print(a[a$quantity %in% names(aggs), ], row.names = FALSE, digits = 3)
cat("\n")
print(a[!a$quantity %in% names(aggs), ], row.names = FALSE, digits = 3)
