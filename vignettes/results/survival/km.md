
``` r
km <- survfit(Surv(days, death) ~ 1, data = rhc)

summary(km, times = c(30, 180, 365))
#> Call: survfit(formula = Surv(days, death) ~ 1, data = rhc)
#> 
#>  time n.risk n.event survival std.err lower 95% CI upper 95% CI
#>    30    995     513    0.658  0.0123        0.634        0.682
#>   180    705     245    0.493  0.0129        0.468        0.519
#>   365    163      91    0.342  0.0166        0.311        0.376
```