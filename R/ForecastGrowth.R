#' Forecast growth
#'
#' Fits various weight-at-age models given input data
#' @param form model formula for machine learning based methods and selection of covariates
#' @param data data.frame with columns "year", "age", and "weight" and any covariates. Weight should be in kg so that uncertainty calculations do not cause an issue.
#' @param weights likelihood weights for model.
#' @param n_proj_years the max number of years to project forward (default = 2)
#'
#' @return
#' @export
#'
#' @examples
ForestGrowth <- function(form = formula(weight~age+year), data = NULL, n_proj_years = 2, peels = 3){

  # Data checks ----
  if(sum(c("weight", "age", "year") %in% colnames(data)) != 3){
    stop(print("Columns do not include `weight`, `age`, or `year`"))
  }

  if(peels >= unique(data$year)){
    stop(print("Number of retrospective peels is greater than the number of years"))
  }

  # Fit models ----
  # * VBGF ----
  vbgf <- FitVBGF(
    data = data,
    weights=NULL,
    n_proj_years = n_proj_years
  )

  # * GMRF ----
  gmrf <- FitGMRF(
    data = data,
    weights=NULL,
    n_proj_years = n_proj_years
  )

  # * NN ----
  #FIXME: no likelihood weights
  nn <- fit_nn(
    data,
    n_proj_years = n_proj_years)

}
