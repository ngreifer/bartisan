# Summarize _dev/bcf-effect-sparsity.rds: cell means for each setting, and the
# paired (on minus off) difference per replicate with its standard error, since
# both settings see the same data and propensity score in every replicate.
x <- readRDS("_dev/bcf-effect-sparsity.rds")
res <- x$res
cat(sprintf("complete: %s, rows: %d, replicates: %d\n\n", x$complete, nrow(res),
            length(unique(res$rep))))

cols <- c("cate_rmse", "cate_cover", "cate_width", "ate_bias", "ate_cover",
          "ate_width", "share_true", "ess_ate", "atom", "seconds")
m <- aggregate(res[cols], res[c("p", "tau", "modset", "sparsity")], mean,
               na.rm = TRUE)
m <- m[order(m$p, m$tau, m$modset, m$sparsity), ]
cat("Cell means\n")
print(format(m, digits = 3), row.names = FALSE)

# Paired differences: on minus off, per replicate and cell.
key <- c("rep", "p", "tau", "modset")
on <- res[res$sparsity == "on", ]
off <- res[res$sparsity == "off", ]
mm <- merge(on, off, by = key, suffixes = c(".on", ".off"))
diffs <- data.frame(mm[c("p", "tau", "modset")],
                    d_cate_rmse = mm$cate_rmse.on - mm$cate_rmse.off,
                    d_abs_bias = abs(mm$ate_bias.on) - abs(mm$ate_bias.off),
                    d_ate_width = mm$ate_width.on - mm$ate_width.off,
                    d_cate_width = mm$cate_width.on - mm$cate_width.off,
                    d_ess = mm$ess_ate.on - mm$ess_ate.off,
                    on_better = mm$cate_rmse.on < mm$cate_rmse.off)
se <- function(v) sd(v) / sqrt(length(v))
agg <- aggregate(cbind(d_cate_rmse, d_abs_bias, d_ate_width, d_cate_width,
                       d_ess, on_better) ~ p + tau + modset, diffs, mean)
agg_se <- aggregate(cbind(d_cate_rmse, d_ess) ~ p + tau + modset, diffs, se)
names(agg_se)[4:5] <- c("se_d_cate_rmse", "se_d_ess")
out <- merge(agg, agg_se)
out <- out[order(out$p, out$tau, out$modset), ]
cat("\nPaired differences, on minus off (negative d_cate_rmse favors on)\n")
print(format(out, digits = 3), row.names = FALSE)
