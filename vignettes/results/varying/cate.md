
``` r
cate <- estimate_effect(fit_vc, treat = "rhc", estimand = "CATE")

plot(cate)
```

<img src="results/varying/cate-1.png" alt="" style="display: block; margin: auto;" />