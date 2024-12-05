

#' Fit deep neural net
#'
#' @param data
#' @param n_proj_years
#'
#' @return
#' @export
#'
#' @examples
fit_nn <- function(data,
                   n_proj_years = 2){
  require(neuralnet)

  # Fit deep NN ----
  nn <- neuralnet(log(weight) ~ age+year,
                  data = data,
                  hidden = c(5,5,5,5),
                  linear.output = TRUE,
                  stepmax = 1e6,
                  lifesign = 'minimal',
                  rep=1)

  # Predict ----
  years <- do.call(seq, as.list(range(data$year)))
  proj_years <- (max(years) + 1):(max(years) + n_proj_years)
  ages <- do.call(seq, as.list(range(data$age)))
  years_ages <- expand.grid(proj_years, ages)
  colnames(years_ages) <- c("year", "age")

  years_ages$pred_weight <- exp(predict(nn, newdata = years_ages))
  years_ages <- years_ages %>%
    dplyr::select(year, age, pred_weight) %>%
    as.data.frame()

  # Return ----
  return(list(obj = nn, data = data, prediction = years_ages))
}

