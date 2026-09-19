# Is the bad end of the cutpoint profile the pinned index, or the low end of the
# response?
#
# WHAT IS BEING TESTED
#
# `_dev/ordinal-cut-chart.R` establishes that the reported chart is not what
# produces the degradation: the sampler's own cutpoints mix worse than the
# reported ones, uniformly, and the reported chart cancels a slow component that
# all of them carry. What is left to settle is which feature of the model that
# slow component belongs to, and two mechanisms predict the same profile.
#
#   The pin. `update_ordinal_cuts()` never writes `cuts(0)`, so threshold 1 is
#   the one threshold that receives no update of its own. Its position relative
#   to the data can move only when the forest's level moves, which happens
#   through individual tree updates rather than through a draw of its own.
#
#   The response. The lowest category is an end of the latent scale, where the
#   fewest observations lie near the boundary, so its threshold is the least
#   informed whatever the sampler does with it.
#
# The claim under test: **the profile is anchored to the pinned index and not to
# the low end of the response.** It could be false, and the measurement is built
# so that it can come out either way rather than confirming the first reading.
#
# WHAT THE OUTCOME MEASUREMENT IS
#
# The same data as the cheapest cell of the main run (n = 500, 20 categories,
# the Friedman latent function, hard rules, eight chains of 1000 draws after 200
# warmup), fitted twice in each of three replicates: once as generated, and once
# with the category order reversed. Reversing carries the original boundary
# between categories 1 and 2 onto new threshold 19, while the pin stays where it
# always is, on new threshold 1. The two codings describe the same data and the
# same partition of the latent scale; only which end the sampler pins changes.
#
# **The column to read is `ess_bulk` by cutpoint index in each coding**, and the
# question is where its minimum sits. `cor(aux.cutk, aux.cut1)` is recorded
# beside it, since the leak into the neighbors is what the gradient is made of.
#
# WHAT EACH OUTCOME WOULD MEAN
#
#   The reversed fit is still worst at k = 1. The pin causes it. The slow
#   quantity is the position of whichever threshold the sampler holds fixed, the
#   low categories have nothing to do with it, and the remedy is a move for the
#   location rather than anything about the cutpoint block. The tridiagonal
#   entry addresses the leak into the neighboring cutpoints and not the scalar
#   itself, and the record should say so rather than leaving it as the natural
#   fix.
#
#   The reversed fit is worst at k = 19. The response causes it, the pinned
#   index is incidental, and what the profile records is that a threshold in the
#   tail of the latent distribution is poorly determined. That is a property of
#   the model rather than of the sampler, nothing needs fixing, and the
#   documentation should say that the extreme cutpoints are the least informed.
#
#   Both ends are bad in both codings. The two mechanisms are separate and the
#   inverse-U in the main run is their sum, which would need the size of each
#   before either could be acted on.

library(bartisan)

N <- 500L
P <- 10L
K_ORD <- 20L
CHAINS <- 8L
REPS <- 3L

FAMILIES <- c("gaussian", "probit", "ordinal")

sim <- function(rep) {
  set.seed(52000 + 211 * rep + match("ordinal", FAMILIES))
  X <- matrix(stats::runif(N * P), N, P,
              dimnames = list(NULL, paste0("x", seq_len(P))))
  f <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 +
    10 * X[, 4] + 5 * X[, 5]
  eta <- (f - mean(f)) / stats::sd(f)

  latent <- eta + stats::rlogis(N)
  cuts <- stats::quantile(latent, seq_len(K_ORD - 1L) / K_ORD)
  y <- factor(findInterval(latent, cuts), levels = 0:(K_ORD - 1L),
              ordered = TRUE)

  cbind(data.frame(y = y), as.data.frame(X))
}

stats_of <- function(aux, k) {
  s <- bartisan:::diagnosis_stats(bartisan:::as_chains(aux[, k], CHAINS))

  data.frame(k = k, rhat = s[["rhat"]], ess_bulk = s[["ess_bulk"]],
             sd = stats::sd(aux[, k]),
             cor_cut1 = stats::cor(aux[, k], aux[, 1L]))
}

rows <- list()

for (rep in seq_len(REPS)) {
  d <- sim(rep)

  for (coding in c("as generated", "reversed")) {
    dd <- d

    if (identical(coding, "reversed")) {
      # The same partition of the latent scale, read from the other end, so the
      # original boundary between categories 1 and 2 becomes threshold 19.
      dd$y <- factor(dd$y, levels = rev(levels(dd$y)), ordered = TRUE)
    }

    set.seed(3300 + rep)
    fit <- bartisan(y ~ ., data = dd, family = ordinal(), chains = CHAINS,
                    num_burn = 200L, num_draws = 1000L, num_thin = 1L,
                    gate = "hard")

    aux <- fit[["aux"]]
    part <- do.call(rbind, lapply(seq_len(ncol(aux)), function(k) {
      stats_of(aux, k)
    }))

    rows[[length(rows) + 1L]] <- cbind(data.frame(coding = coding, rep = rep),
                                       part)
    saveRDS(do.call(rbind, rows), "_dev/ordinal-cut-reversed.rds")
  }
}

res <- do.call(rbind, rows)
agg <- aggregate(cbind(ess_bulk, rhat, sd, cor_cut1) ~ coding + k, res, mean)

cat("bulk effective sample size by cutpoint index, mean of", REPS,
    "replicates\n\n")
cat(" k | as generated   reversed |  cor(.,cut1) gen / rev\n")

for (k in sort(unique(agg$k))) {
  g <- agg[agg$coding == "as generated" & agg$k == k, ]
  r <- agg[agg$coding == "reversed" & agg$k == k, ]
  cat(sprintf("%2d | %12.0f %10.0f | %7.3f %7.3f\n", k, g$ess_bulk, r$ess_bulk,
              g$cor_cut1, r$cor_cut1))
}

for (coding in c("as generated", "reversed")) {
  a <- agg[agg$coding == coding, ]
  cat(sprintf("\n%-13s worst index k = %d (ess %.0f), best k = %d (ess %.0f)\n",
              coding, a$k[which.min(a$ess_bulk)], min(a$ess_bulk),
              a$k[which.max(a$ess_bulk)], max(a$ess_bulk)))
}
