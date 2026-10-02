# Which gate handles frequency weights correctly?
#
# What is tested. In the reseeded recovery matrix of 2026-10-01 the two gates
# agreed within 0.1 in every cell but the weights column, where on identical
# data binomial logit read 1.53 under hard rules against 1.13 under soft,
# cloglog 1.15 against 0.65, ordinal logit 1.15 against 0.97, ph() 0.93 against
# 0.73 and ordbeta 0.85 against 0.71. Gates that agree everywhere else cannot
# both be right with weights. The package's own claim is that a weighted fit
# and a fit on row-replicated data are the same fit, which is the discriminator.
#
# What is measured. binomial(), binomial("cloglog"), ordbeta() and ph() on the
# matrix's data, four fits each: weights w under hard rules, the rows
# replicated w times under hard rules, and the same two under soft rules. Per
# fit the effect (mean over units of link(z = 1) - link(z = 0)) and its
# posterior sd, 50 trees, 300 + 500 draws.
#
# What each outcome means. A gate whose weighted fit matches its replicated
# fit within Monte Carlo error handles weights as documented. One whose
# weighted fit sits away from its replicated fit has a defect in how weights
# enter that family's leaf target under that gate; the identity test's weights
# arm passing then says the density route carries the same weight the sampler
# did, right or wrong.
#
# Run with: Rscript _dev/weights-check.R
# Writes:   _dev/weights-check.rds

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
mu <- plogis(lp)
d$ybeta <- drawn(11L, pmin(pmax(rbeta(n, mu * 10, (1 - mu) * 10), 1e-4), 1 - 1e-4))
d$yob <- drawn(12L, {
  p0 <- 1 - plogis(lp + 1.5)
  p1 <- plogis(lp - 1.5)
  which_bound <- runif(n)
  yob <- d$ybeta
  yob[which_bound < p0] <- 0
  yob[which_bound > 1 - p1] <- 1
  yob
})
cens <- drawn(16L, lp + rnorm(n, 0.6, 0.7))
lt <- drawn(21L, log((-log(runif(n)) / exp(lp))^(1 / 1.5)))
d$t_ph <- exp(pmin(lt, cens))
d$e_ph <- as.integer(lt <= cens)

expanded <- d[rep(seq_len(n), d$w), ]

cases <- list(
  list(name = "binomial logit", family = binomial(), response = "ybin"),
  list(name = "binomial cloglog", family = binomial("cloglog"), response = "ycll"),
  list(name = "ordbeta", family = ordbeta(), response = "yob"),
  list(name = "ph", family = ph(), response = "cbind(t_ph, e_ph)")
)
# WEIGHTS_SEEDS names fit seeds beyond the first, for the two augmented binomial
# families only, to tell a repeatable gap from Monte Carlo error.
extra <- as.integer(strsplit(Sys.getenv("WEIGHTS_SEEDS", ""), ",")[[1L]])
grid <- expand.grid(case = seq_along(cases), gate = c("hard", "smoothstep"),
                    route = c("weighted", "replicated"), seed = 7L,
                    stringsAsFactors = FALSE)
if (length(extra)) {
  grid <- rbind(grid, expand.grid(case = 1:2, gate = "hard",
                                  route = c("weighted", "replicated"),
                                  seed = extra, stringsAsFactors = FALSE))
}
if (nzchar(Sys.getenv("WEIGHTS_EXTRA_ONLY"))) grid <- grid[grid$seed != 7L, ]

pr <- prog_init(total = nrow(grid), title = "Weights against replication",
                unit = "fit", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]
  case <- cases[[g$case]]
  form <- as.formula(paste(case$response, "~ z + x1 + x2 + x3"), env = globalenv())
  control <- bartisan_control(num_trees = 50L, num_burn = 300L, num_draws = 500L,
                              gate = g$gate, verbose = FALSE)
  set.seed(g$seed)
  fit <- {
    if (identical(g$route, "weighted"))
      bartisan(form, data = d, family = case$family, weights = d$w, control = control)
    else
      bartisan(form, data = expanded, family = case$family, control = control)
  }
  d1 <- d
  d1$z <- 1L
  d0 <- d
  d0$z <- 0L
  by_draw <- rowMeans(predict(fit, newdata = d1, type = "link", draws = TRUE) -
                        predict(fit, newdata = d0, type = "link", draws = TRUE))
  rows[[k]] <- data.frame(family = case$name, gate = g$gate, route = g$route,
                          seed = g$seed, effect = mean(by_draw),
                          post_sd = sd(by_draw))
  prog_tick(pr, label = sprintf("%s / %s / %s / seed %d", case$name, g$gate,
                                g$route, g$seed))
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
out <- if (nzchar(Sys.getenv("WEIGHTS_EXTRA_ONLY"))) "_dev/weights-check-seeds.rds" else "_dev/weights-check.rds"
saveRDS(list(rows = res, complete = TRUE), out)
w <- reshape(res[, c("family", "gate", "seed", "route", "effect")],
             idvar = c("family", "gate", "seed"), timevar = "route", direction = "wide")
names(w) <- sub("effect.", "", names(w), fixed = TRUE)
w$gap <- w$weighted - w$replicated
print(w[order(w$family, w$gate, w$seed), ], digits = 3, row.names = FALSE)
cat("\nposterior sd of the effect, by fit:\n")
print(res[, c("family", "gate", "route", "seed", "post_sd")], digits = 2, row.names = FALSE)
