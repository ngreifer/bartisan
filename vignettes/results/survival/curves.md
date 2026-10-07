
``` r
curves_rhc <- avg_predictions(fit_ph, variables = "rhc", type = "survival",
                              times = seq(5, 365, by = 10))

curves_rhc$time <- as.numeric(curves_rhc$group)

ggplot(curves_rhc, aes(x = time, y = estimate, ymin = conf.low,
                       ymax = conf.high, color = factor(rhc),
                       fill = factor(rhc))) +
  geom_ribbon(alpha = .3, color = NA) +
  geom_line() +
  labs(x = "Days since admission", y = "Average survival probability",
       color = "rhc", fill = "rhc")
```

<img src="results/survival/curves-1.png" alt="" style="display: block; margin: auto;" />