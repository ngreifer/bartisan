
``` r
plot(fit, ~ meanbp + aps,
     values = list(aps = function(x) quantile(x, c(.05, .5, .95)))) +
  labs(x = "Mean arterial blood pressure", y = "Fitted probability of death",
       color = "APACHE III", fill = "APACHE III")
```

<img src="results/effects/pdp4-1.png" alt="" style="display: block; margin: auto;" />