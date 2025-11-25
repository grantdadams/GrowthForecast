

#' Function to fit XGBoost (eXtreme Gradient Boosting) to weight-at-age data
#'
#' @param form model formula
#' @param data data.frame with columns "year", "age", and "weight" and other potential covariates. Weight should be in kg so that uncertainty calculations do not cause an issue.
#' @param weights likelihood weights for model.
#' @param n_proj_years the number of years to project forward
#' @param max.depth maximum depth of a tree. Default: 6
#' @param nthread description
#' @param nrounds max number of boosting iterations. Default: 2
#'
#' @return an XGBoost model object
#' @description
#' This fits a basic XGBoost model to weight-at-age data using squared log error. The model specification can be optimized for the data-set. See ?xgboost for more information
#' https://xgboost.readthedocs.io/en/stable/R-package/xgboostPresentation.html
#'
#' @export
#'
FitXGBoost <- function(
    data = NULL,
    n_proj_years = 2,
    last_year = NA,
    max.depth = 10,
    nthread = 1,
    nrounds = 20){

  # - Fit model
  bst_model <- xgboost::xgboost(data = data %>%
                         dplyr::select(-weight) %>%
                         as.matrix(),
                       label = data$weight,
                       max.depth = max.depth,
                       nthread = nthread,
                       nrounds = nrounds,
                       objective = "reg:squaredlogerror")

  # - Predict
  # -- Build data set for X year projection
  years <- do.call(seq, as.list(range(data$year)))
  proj_years <- (max(years) + 1):(max(years) + n_proj_years)
  ages <- do.call(seq, as.list(range(data$age)))

  # - Projections and fill missing years
  years_ages <- expand.grid(c(years, proj_years), ages)
  colnames(years_ages) <- c("year", "age")

  # - Combine
  proj_data <- data %>%
    dplyr::select(-weight) %>%
    dplyr::full_join(years_ages) %>%
    dplyr::arrange(year, age)

  # -- Fill projection data columns with NA with column mean
  proj_data[] <- lapply(proj_data, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))
  proj_data <- proj_data %>%
    dplyr::filter(year >  max(years))

  # -- Predict
  proj_data$pred_weight <-  xgboost::predict(bst_model, as.matrix(proj_data))

  # Predicted weight-at-age matrix ----
  predicted_weight = proj_data %>%
    dplyr::mutate(projection = year %in% proj_years,
                  last_year = last_year,
                  model = "XGBOOST") %>%
    dplyr::select(model, year, age, pred_weight, last_year, projection) %>%
    arrange(year, age) %>%
    as.data.frame()

  # Predicted observations ----
  data$predicted_weight <-  xgboost::predict(bst_model)
  data <- data %>%
    arrange(year, age)

  return(list(data = data, model = bst_model, predicted_weight = predicted_weight))
}
