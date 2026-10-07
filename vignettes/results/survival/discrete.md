
``` r
# Interval edges at quantiles of the observed death times
edges <- quantile(rhc$days[rhc$death == 1], seq(.05, .95, length.out = 15))

# How many intervals each patient's time reaches
reached <- pmax(findInterval(rhc$days, edges), 1)

# One row per patient per interval reached, with the event, if
# there was one, in the last of them
long <- rhc[rep(seq_len(nrow(rhc)), reached), ]
long$interval <- sequence(reached)
long$event <- 0
long$event[cumsum(reached)] <- rhc$death

fit_dt <- bartisan(event ~ interval + rhc + age + sex + race + edu + aps +
                     meanbp + resp + hema + pafi + paco2 + crea + surv2m +
                     card,
                   data = long, family = binomial("probit"))
```