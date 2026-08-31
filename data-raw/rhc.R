# Builds data/rhc.rda from the original SUPPORT right heart catheterization file.
#
#   source: https://hbiostat.org/data/repo/rhc.csv
#   study:  Connors et al. (1996), JAMA 276(11), 889-897
#
# The covariate set is the thirteen variables used in
# https://iqss.github.io/dss-ps/example.html, all measured before catheterization.
# The outcome is death during follow-up, and `days` is the time to it, so the same
# event supports a binary analysis and a right-censored survival one.

raw <- read.csv("https://hbiostat.org/data/repo/rhc.csv", stringsAsFactors = FALSE)

rhc <- data.frame(
  # Treatment, and the outcome in its two forms. Both are 0/1 so that a contrast
  # between them is a single number rather than one per level of a factor.
  rhc    = as.integer(raw$swang1 == "RHC"),
  death  = as.integer(raw$death == "Yes"),
  days   = ifelse(raw$death == "Yes",
                  raw$dthdte - raw$sadmdte,
                  raw$lstctdte - raw$sadmdte),

  # Demographics.
  age    = raw$age,
  sex    = factor(raw$sex, levels = c("Female", "Male"), labels = c("female", "male")),
  race   = factor(raw$race, levels = c("white", "black", "other")),
  edu    = raw$edu,

  # Physiology on day 1.
  aps    = raw$aps1,
  meanbp = raw$meanbp1,
  resp   = raw$resp1,
  hema   = raw$hema1,
  pafi   = raw$pafi1,
  paco2  = raw$paco21,
  crea   = raw$crea1,

  # Prognosis and comorbidity.
  surv2m = raw$surv2md1,
  card   = factor(tolower(raw$card), levels = c("no", "yes"))
)

stopifnot(nrow(rhc) == 5735, !anyNA(rhc))

# A random 1500 of the 5735, which is the sample the package ships. The seed is
# fixed so that the shipped data can be reproduced from this script.
set.seed(2026)
rhc <- rhc[sort(sample(nrow(rhc), 1500L)), ]
rownames(rhc) <- NULL

stopifnot(nrow(rhc) == 1500L, !anyNA(rhc))

usethis::use_data(rhc, overwrite = TRUE)
