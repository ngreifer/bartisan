
``` r
library(ggplot2)

ggplot(dz, aes(count)) +
  geom_histogram(binwidth = 1, fill = "grey70", color = "white") +
  labs(x = "Count", y = "Observations",
       subtitle = "Poisson, no zero inflation, marginal mean 1.9") +
  theme_bw()
```

![](results/families/zerohist-1.png)<!-- -->