
``` r
plot_predictions(fit, condition = list(aps = at), draw = FALSE)
#>   rowid estimate conf.low conf.high  df   age card  crea   edu  hema meanbp
#> 1     1   0.6299   0.4969    0.7594 Inf 61.42   no 2.132 11.64 31.68     78
#> 2     2   0.6817   0.5604    0.7867 Inf 61.42   no 2.132 11.64 31.68     78
#> 3     3   0.7386   0.6018    0.8416 Inf 61.42   no 2.132 11.64 31.68     78
#>   paco2  pafi  race resp rhc  sex surv2m  aps
#> 1 38.85 217.4 white   28   0 male  0.587 29.9
#> 2 38.85 217.4 white   28   0 male  0.587 54.0
#> 3 38.85 217.4 white   28   0 male  0.587 83.0
```