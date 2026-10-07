
``` r
library(ggplot2)

times <- seq(5, 730, by = 25)

# A fit's predicted survival curve, averaged over the patients
avg_curve <- function(fit, family) {
  data.frame(family = family, time = times,
             surv = colMeans(predict(fit, type = "survival",
                                     times = times)))
}

curves <- rbind(avg_curve(fit_dpm, "dpm_aft()"),
                avg_curve(fit_ln, "lognormal_aft()"),
                avg_curve(fit_ph, "ph()"))

ggplot(curves, aes(x = time, y = surv)) +
  geom_step(data = data.frame(time = c(0, km$time), surv = c(1, km$surv))) +
  geom_line(aes(color = family)) +
  coord_cartesian(xlim = c(0, 730), ylim = c(0, 1)) +
  labs(x = "Days since admission", y = "Survival probability",
       color = "Family")
```

<img src="results/survival/kmplot-1.png" alt="" style="display: block; margin: auto;" />