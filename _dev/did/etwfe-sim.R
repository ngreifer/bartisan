# Does the forest version of Wooldridge's ETWFE remove the bias his linear
# version has when its linearity assumption fails?
#
# What is being tested. Wooldridge (2025) identifies ATT(g, t) under
# conditional parallel trends plus Assumption LIN, linearity of the
# never-treated mean in the covariates, and his concluding remarks name LIN as
# the one thing the POLS/ETWFE approach does not relax. The bartisan version
# (`_dev/did/etwfe-mpdta.R`) replaces the linear terms with two additive
# forests, a(cohort, x) + b(period, x), fit on the untreated rows, and imputes
# the treated ones. The claim is that it is close to unbiased where LIN fails
# and linear ETWFE is not, and costs little where LIN holds. It could be false
# if the period forest, learned only from untreated rows, cannot carry a
# covariate-dependent trend into the treated cells, whose covariate mix differs
# from the controls'.
#
# The measurement. 500 units by 6 periods, cohorts first treated at 4, 5 and 6
# plus never-treated, unit effects of sd 1 and noise of sd 0.3, cohort
# membership depending on the covariates. Conditional parallel trends holds
# exactly: the never-treated mean moves by m(x) (t - 1) / 5. Two cells:
# "nonlinear", m(x) = 6 (x1 - 1/2)^2 + 2 sin(2 pi x2), where LIN is false and
# selection is on the same features (cohort 4 U-shaped in x1, cohort 5 on
# sin(2 pi x2)), which is what makes a linear-in-x trend biased rather than
# merely noisy; "linear", m(x) = 1.5 x1 + x2, where LIN holds. A first design
# selected on features m did not share (cohort 4 U-shaped in x1 against
# sin(2 pi x1), which averages to zero over it) and left linear ETWFE unbiased
# and only wider, so it could not test the claim. The effect is 1 + 0.25 (t - g) + 0.5 (x2 - 1/2). 20 replicates per
# cell. Estimators: linear ETWFE by imputation, with unit-clustered intervals
# from the equivalent POLS regression; bartisan imputing the conditional mean
# of Y(inf); and bartisan imputing posterior predictive draws of it. Read the
# bias and coverage of the overall ATT first (the truth is the average effect
# over the treated rows), then the event-time ATTs at e = 0, 1, 2, then width.
#
# What each outcome means. bartisan near-unbiased in the nonlinear cell where
# linear ETWFE is not, and close to it in the linear cell: the forest version
# earns its place and the notes recommend it where trends may depend on
# covariates nonlinearly. bartisan biased too: the period forest cannot carry
# the trend into the treated cells and imputation with forests needs more (for
# example, more trees on the period forest). bartisan much wider in the linear
# cell: flexibility has a real price there and the notes say to prefer the
# linear form when a linear-in-x trend is defensible. Of the two imputations,
# the conditional mean leaves out the treated rows' own noise and may
# undercover, the predictive draws treat that noise as unrelated to the
# observed outcome and may overcover; whichever is near 95% is the one to
# recommend.

source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")
suppressMessages({
  library(bartisan)
  library(fixest)
})

reps <- as.integer(Sys.getenv("REPS", "20"))
smoke <- nzchar(Sys.getenv("SMOKE"))
out_file <- Sys.getenv("OUT", "_dev/did/etwfe-sim.rds")
future::plan(future::multisession, workers = 4L)

n <- 500L
periods <- 1:6
cohorts <- c(4L, 5L, 6L)

m_fun <- list(
  nonlinear = function(x1, x2) 6 * (x1 - 0.5)^2 + 2 * sin(2 * pi * x2),
  linear = function(x1, x2) 1.5 * x1 + x2)

simulate <- function(cell, seed) {
  set.seed(seed)
  x1 <- stats::runif(n)
  x2 <- stats::runif(n)
  score <- cbind(0, 12 * (x1 - 0.5)^2 - 1.5, 1.5 * sin(2 * pi * x2) - 0.5,
                 x1 - 0.5)
  prob <- exp(score) / rowSums(exp(score))
  g <- c(Inf, cohorts)[apply(prob, 1L, function(p) sample.int(4L, 1L, prob = p))]
  level <- c(0, 0.8, -0.5, 0.3)[match(g, c(Inf, cohorts))]
  alpha <- stats::rnorm(n)
  m <- m_fun[[cell]](x1, x2)

  d <- expand.grid(id = seq_len(n), period = periods)
  i <- d$id
  d$x1 <- x1[i]
  d$x2 <- x2[i]
  d$g <- g[i]
  d$w <- as.integer(is.finite(d$g) & d$period >= d$g)
  d$tau <- ifelse(d$w == 1L, 1 + 0.25 * (d$period - d$g) + 0.5 * (d$x2 - 0.5), 0)
  d$y <- alpha[i] + level[i] + 0.5 * d$x2 + 0.3 * d$period +
    m[i] * (d$period - 1) / 5 + d$tau + stats::rnorm(nrow(d), sd = 0.3)
  d$cohort <- factor(ifelse(is.finite(d$g), d$g, 0), levels = c(0, cohorts))
  d$periodf <- factor(d$period)
  d$one <- 1
  d$cell <- ifelse(d$w == 1L, paste(d$g, d$period, sep = ":"), "none")
  d$e <- ifelse(d$w == 1L, d$period - d$g, NA)
  d
}

# Weights over the treated cells for each reported quantity.
aggregations <- function(tr) {
  cells <- sort(unique(tr$cell))
  n_gt <- as.numeric(table(tr$cell)[cells])
  e_of <- tapply(tr$e, tr$cell, `[`, 1L)[cells]
  out <- list(overall = n_gt)
  for (e in 0:2) {
    out[[paste0("e", e)]] <- ifelse(e_of == e, n_gt, 0)
  }
  lapply(out, function(w) w / sum(w))
}

linear_etwfe <- function(d) {
  tr <- d[d$w == 1L, ]
  cells <- sort(unique(tr$cell))
  xbar1 <- tapply(d$x1, d$cohort, mean)
  xbar2 <- tapply(d$x2, d$cohort, mean)
  d$xd1 <- d$x1 - xbar1[as.character(d$cohort)]
  d$xd2 <- d$x2 - xbar2[as.character(d$cohort)]
  terms <- character()
  for (k in cells) {
    nm <- paste0("c", gsub(":", "_", k))
    d[[nm]] <- as.integer(d$cell == k)
    d[[paste0(nm, "_x1")]] <- d[[nm]] * d$xd1
    d[[paste0(nm, "_x2")]] <- d[[nm]] * d$xd2
    terms <- c(terms, nm)
  }
  f <- stats::as.formula(paste(
    "y ~", paste(c(terms, paste0(terms, "_x1"), paste0(terms, "_x2")),
                 collapse = " + "),
    "+ cohort * (x1 + x2) + periodf * (x1 + x2)"))
  fit <- feols(f, data = d, cluster = ~ id)
  b <- coef(fit)[terms]
  V <- vcov(fit)[terms, terms]
  lapply(aggregations(tr), function(w) {
    est <- sum(w * b)
    se <- sqrt(drop(t(w) %*% V %*% w))
    c(estimate = est, lower = est - 1.96 * se, upper = est + 1.96 * se)
  })
}

bart_etwfe <- function(d) {
  ctl <- d[d$w == 0L, ]
  tr <- d[d$w == 1L, ]
  fit <- bartisan(y ~ cohort + x1 + x2 + vc(one, ~ period + x1 + x2) + (1 | id),
                  data = ctl, family = gaussian(), chains = 4L, gate = "hard",
                  num_trees = c(100L, 50L),
                  num_burn = if (smoke) 20L else 500L,
                  num_draws = if (smoke) 20L else 1000L, verbose = FALSE)
  cells <- sort(unique(tr$cell))
  w <- aggregations(tr)
  summarize <- function(y0) {
    te <- sweep(-y0, 2L, tr$y, "+")
    cell_draws <- sapply(cells, function(k) rowMeans(te[, tr$cell == k, drop = FALSE]))
    lapply(w, function(wt) {
      v <- drop(cell_draws %*% wt)
      c(estimate = mean(v), lower = unname(stats::quantile(v, .025)),
        upper = unname(stats::quantile(v, .975)))
    })
  }
  list(mean = summarize(predict(fit, newdata = tr, draws = TRUE)),
       predictive = summarize(rstantools::posterior_predict(fit, newdata = tr)))
}

truth <- function(d) {
  tr <- d[d$w == 1L, ]
  cells <- sort(unique(tr$cell))
  cell_truth <- tapply(tr$tau, tr$cell, mean)[cells]
  vapply(aggregations(tr), function(w) sum(w * cell_truth), numeric(1L))
}

cells_run <- names(m_fun)
pr <- prog_init(total = reps * length(cells_run),
                title = "ETWFE, linear against forests", unit = "replicate",
                kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)
rows <- list()

for (r in seq_len(reps)) {
  for (cell in cells_run) {
    d <- simulate(cell, seed = 1000L * r + match(cell, cells_run))
    tv <- truth(d)
    t0 <- Sys.time()
    est <- c(list(linear = linear_etwfe(d)),
             setNames(bart_etwfe(d), c("bart_mean", "bart_predictive")))
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    for (method in names(est)) {
      for (q in names(tv)) {
        e <- est[[method]][[q]]
        rows[[length(rows) + 1L]] <- data.frame(
          rep = r, cell = cell, method = method, quantity = q,
          truth = tv[[q]], estimate = e[["estimate"]], lower = e[["lower"]],
          upper = e[["upper"]], seconds = round(secs, 1))
      }
    }
    prog_tick(pr, label = sprintf("rep %d %s", r, cell))
  }
  saveRDS(list(res = do.call(rbind, rows), complete = FALSE,
               done = r, total = reps), out_file)
}

res <- do.call(rbind, rows)
saveRDS(list(res = res, complete = TRUE, done = reps, total = reps), out_file)
on.exit()
prog_end(pr, "done")

res$bias <- res$estimate - res$truth
res$covers <- res$lower <= res$truth & res$truth <= res$upper
res$width <- res$upper - res$lower
a <- aggregate(cbind(bias, covers, width) ~ cell + quantity + method, res, mean)
a$rmse <- aggregate(bias ~ cell + quantity + method, res,
                    function(b) sqrt(mean(b^2)))$bias
a <- a[order(a$cell, a$quantity, a$method), ]
print(a, row.names = FALSE, digits = 3)
