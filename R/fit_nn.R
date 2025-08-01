# https://www.tidymodels.org/learn/models/parsnip-nnet/

#' Fit deep neural net
#'
#' @param form model formula to use for neural net
#' @param data data.frame with columns "year", "age", and "weight". Weight should be in kg so that uncertainty calculations do not cause an issue.
#' @param n_proj_years the number of years to project forward
#' @param startweights starting parameters
#' @param last_year final year of data
#'
#' @return
#' @export
#'
#' @examples
fit_nn <- function(
    form = formula(log(weight) ~ age + year),
    data,
    n_proj_years = 2,
    startweights = NULL,
    last_year = NA){
  require(neuralnet)

  # Fit deep NN ----
  nn <- neuralnet(formula = form,
                  data = data,
                  hidden = c(5,5,5,5,5),
                  startweights = startweights,
                  linear.output = TRUE,
                  stepmax = 1e6,
                  lifesign = 'minimal',
                  rep=1)

  # Predict ----
  years <- do.call(seq, as.list(range(data$year)))
  proj_years <- (max(years) + 1):(max(years) + n_proj_years)
  ages <- do.call(seq, as.list(range(data$age)))

  # - Projections and fill missing years
  years_ages <- expand.grid(c(years, proj_years), ages)
  colnames(years_ages) <- c("year", "age")

  # - Prediction for each obs
  # data$pred_weight = predict(nn, newdata = data)

  # - Predicted for hind/forecast
  years_ages$pred_weight <- exp(predict(nn, newdata = years_ages))

  years_ages <- years_ages %>%
    # dplyr::filter(year %in% proj_years) %>%
    dplyr::select(year, age, pred_weight) %>%
    merge(., data, by = c('year','age'), all = TRUE) %>%
    dplyr::mutate(model = "nn",
                  projection = year %in% proj_years,
                  last_year = last_year) %>%
    as.data.frame()%>%
    dplyr::arrange(year, age) %>%
    dplyr::select(model, last_year, year, age, obs_weight = weight, pred_weight, projection)

  # Return ----
  return(list(obj = nn, data = data, prediction = years_ages))
}

