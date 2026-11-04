
``` r
predict(bart_fit, rhc_test[1:5, ], type = "prob")
#> # A tibble: 5 × 2
#>   .pred_yes .pred_no
#>       <dbl>    <dbl>
#> 1     0.792   0.208 
#> 2     0.202   0.798 
#> 3     0.518   0.482 
#> 4     0.924   0.0760
#> 5     0.459   0.541

predict(bart_fit, rhc_test[1:5, ], type = "conf_int", level = .9)
#> # A tibble: 5 × 4
#>   .pred_lower_yes .pred_lower_no .pred_upper_yes .pred_upper_no
#>             <dbl>          <dbl>           <dbl>          <dbl>
#> 1           0.684         0.117            0.883          0.316
#> 2           0.108         0.681            0.319          0.892
#> 3           0.389         0.360            0.640          0.611
#> 4           0.864         0.0351           0.965          0.136
#> 5           0.308         0.377            0.623          0.692
```