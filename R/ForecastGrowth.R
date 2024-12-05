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

  if(peels >= length(unique(data$year))){
    stop(print("Number of retrospective peels is greater than the number of years"))
  }


  # Run peels
  peel_list <- list()
  projection_list <- list()
  for(i in 1:peels){

    # Peel data ----
    last_year <- max(data$year, na.rm = TRUE) - i
    train <- data %>%
      dplyr::filter(year <= last_year)

    test <- data %>%
      dplyr::filter(year <= (last_year + n_proj_years) & year > last_year) # Only projection years


    # Fit models ----
    # * VBGF ----
    vbgf <- FitVBGF(
      data = train,
      weights=NULL,
      n_proj_years = n_proj_years,
      last_year = last_year
    )

    # * GMRF ----
    gmrf <- FitGMRF(
      data = train,
      weights=NULL,
      n_proj_years = n_proj_years,
      last_year = last_year
    )

    # * NN ----
    #FIXME: no likelihood weights
    nn <- fit_nn(
      data = train,
      n_proj_years = n_proj_years,
      last_year = last_year)


    # Combine ----
    peel_list[[i]] <- list(vbgf = vbgf,
                           gmrf1 = gmrf[[1]],
                           gmrf2 = gmrf[[2]],
                           gmrf3 = gmrf[[3]],
                           gmrf4 = gmrf[[4]],
                           nn = nn)

    projection_list[[i]] <- do.call("rbind",
                                    lapply(peel_list[[i]],
                                           function(x) x$prediction)
    ) %>%
      tidyr::pivot_wider(names_from = c(model), values_from = pred_weight)

  }
  names(peel_list) <- 1:peels
  names(projection_list) <- 1:peels


  # Performance metrics ----

}
