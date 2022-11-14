#' A toy data set used for modeling the FHC model
#'
#' This toy data set has 500 patients' survival and baseline information and
#' contains a cure fraction.
#'
#'@format A data frame with 500 observations on the following 5 variables.
#' \describe{
#'     \item{id}{patient id}
#'     \item{age}{a continuous variable, which is centered to the mean}
#'     \item{sex}{a binary variable with 0 for male and 1 for female}
#'     \item{time}{observed survival time}
#'     \item{event}{an event indicator with 1 for dead and 0 for censored}
#'     }
#'
#'@source{Generated from a FHC model to serve as an example.}
#'
"fhc_dat"
