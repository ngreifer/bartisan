
``` r
plot(fit, ~ meanbp + aps) +
  labs(x = "Mean arterial blood pressure", y = "Fitted probability of death",
       color = "APACHE III", fill = "APACHE III")
#> ℹ Grouping by `aps` at 41, 54, and 68, three of its values near its quartiles.
#> ℹ Set `values` to choose them yourself.
```

<img src="results/effects/pdp3-1.png" alt="" style="display: block; margin: auto;" />