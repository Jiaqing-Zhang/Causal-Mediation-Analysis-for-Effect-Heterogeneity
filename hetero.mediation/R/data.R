#' Example Data: Add Health Wave I and Wave II Data
#'
#' A data set containing adolescent smoke use, neighborhood poverty, neighborhood urbanity, friends' substance use, sample weight and other covariates. The total sample size is 2,579.
#'
#' @format A data frame with 2,579 rows and 36 variables:
#'  \describe{
#'   \item{AID}{subject ID}
#'   \item{Urbanity}{living in urban areas, 0=no, 1=yes}
#'   \item{neighborhood_poverty}{reside in disadvantaged neighborhood, 0=no, 1=yes}
#'   \item{wave_II_friend_smoke}{number of close friends smoke}
#'   \item{wave_II_friend_drug}{number of close friends drink}
#'   \item{wave_II_exercise}{engagement in after-school sports or clubs, 0=no, 1=yes}
#'   \item{wave_II_school_safe}{feel secure at school, 0=no, 1=yes}
#'   \item{ysmoke}{whether the subject smoke, yes or no}
#'   \item{weight}{Inverse Propensity Score Weighting}
#' }
"smoke"

