
``` r
# 0-10 Likert where nobody picked 3 or 7: the gaps are ignored
bartisan(score ~ ., data = d, family = ordinal())

# the same scale, with every point modeled
d$score <- ordered(d$score, levels = 0:10)
bartisan(score ~ ., data = d, family = ordinal())
```