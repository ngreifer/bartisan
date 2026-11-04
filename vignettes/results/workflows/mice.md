
``` r
library(mice)

set.seed(1234)
imp <- mice(rhc_mis, method = "bartisan", m = 10, maxit = 10,
            blots = list(pafi = list(family = Gamma("log")),
                         crea = list(family = Gamma("log"))),
            printFlag = FALSE)

imp$method
#>      death        age        sex        aps     meanbp       pafi       crea 
#>         ""         ""         ""         ""         "" "bartisan" "bartisan" 
#>       card 
#>         ""
```