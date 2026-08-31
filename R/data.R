#' Right heart catheterization in critically ill patients
#'
#' Data from the SUPPORT study on whether right heart catheterization within 24
#' hours of admission to an intensive care unit affects survival (Connors et al.,
#' 1996). Catheterization was not randomized, so the comparison is
#' confounded by how sick each patient was on admission, which is what the
#' physiological covariates are for.
#'
#' The outcome appears in two forms. `death` is whether the patient died during
#' follow-up, and `days` is how long that took, so the same event supports a
#' binary analysis that ignores timing and a right-censored survival analysis
#' that does not. A patient who did not die is censored at their last contact.
#'
#' @format A data frame with 1500 rows and 16 columns.
#' \describe{
#'   \item{rhc}{whether the patient received right heart catheterization, 1 or 0.
#'     This is the treatment.}
#'   \item{death}{whether the patient died during follow-up, 1 or 0.}
#'   \item{days}{days from admission to death, or to last contact for a patient
#'     who did not die. Together with `death` this is the survival outcome.}
#'   \item{age}{age in years.}
#'   \item{sex}{`"female"` or `"male"`.}
#'   \item{race}{`"white"`, `"black"` or `"other"`.}
#'   \item{edu}{years of education.}
#'   \item{aps}{APACHE III score on day 1, ignoring coma. Higher is sicker.}
#'   \item{meanbp}{mean blood pressure on day 1.}
#'   \item{resp}{respiratory rate on day 1.}
#'   \item{hema}{hematocrit on day 1.}
#'   \item{pafi}{ratio of arterial oxygen to inspired oxygen on day 1.}
#'   \item{paco2}{arterial carbon dioxide on day 1.}
#'   \item{crea}{serum creatinine on day 1.}
#'   \item{surv2m}{the study's own model-based estimate, made on day 1, of the
#'     probability of surviving two months.}
#'   \item{card}{whether cardiovascular disease was a diagnosis, `"no"` or
#'     `"yes"`.}
#' }
#'
#' @details
#' These are a random 1500 of the 5735 patients in the original file. The
#' covariates are the thirteen used in the worked example at
#' <https://iqss.github.io/dss-ps/example.html>, all recorded before
#' catheterization; the full study collected many more.
#'
#' The treatment and the outcome are coded 0 and 1 rather than as factors, so
#' that a contrast between them is a single number rather than one per level.
#'
#' Nothing here makes the causal assumptions hold. Whether the effect of `rhc`
#' on `death` can be read causally depends on whether these covariates account
#' for how patients were selected for catheterization, which is a question about
#' the study rather than about any model. See `vignette("causal")`.
#'
#' @source Assembled from <https://hbiostat.org/data/repo/rhc.csv> by
#'   `data-raw/rhc.R`.
#'
#' @references
#' Connors, A. F., Speroff, T., Dawson, N. V., et al. (1996). The effectiveness
#' of right heart catheterization in the initial care of critically ill
#' patients. *JAMA*, 276(11), 889--897. \doi{10.1001/jama.1996.03540110043030}
#'
#' @examples
#' data(rhc)
#'
#' # The binary outcome.
#' table(rhc$rhc, rhc$death)
#'
#' # The same event as a survival outcome.
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   with(rhc, summary(survival::Surv(days, death)))
#' }
"rhc"
