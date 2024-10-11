#' @title Length Weight and Age data
#'
#' @description Length Weight and Age data from summer ground surveys conducted by the NOAA NMFS Alaska Fisheries Science Center
#'
#' @format A dataframe with the following columns
#' \describe{
#'   \item{REGION}{Alaska basing (GOA, AI, or BS)}
#'   \item{YEAR}{Year of the survey}
#'   \item{CRUISE_TYPE}{Specifying the survey cruise}
#'   \item{HAULJOIN}{Internally joining unique ID}
#'   \item{STATIONID}{Survey station}
#'   \item{STRATUM}{Survey stratum}
#'   \item{age}{Age of the fish in years}
#'   \item{weight}{Weight of the fish in grams}
#'   \item{Temp}{Bottom temperature at the survey location}
#'   \item{SPECIES_CODE}{Internal code to match species info}
#'   \item{CN}{Common name}
#'   \item{SN}{Species name}
#'   \item{Sp}{Secondary common name}
#'
#' }
#' @source K. Holsman et al. 2024
"LWA"
