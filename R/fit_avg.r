#' Calculate average X years of weight-at-age
#'
#' Uses spline interpolation to fill in missing ages
#'
#' @param data data.frame with columns "year", "age", and "weight". Weight should be in kg so that uncertainty calculations do not cause an issue.
#' @param nyrs number of previous years to calculate average
#' @param n_proj_years the number of years to project forward
#' @param last_year final year of data
average_weight <- function(data,
                           nyrs = 5,
                           n_proj_years = 2,
                           last_year = NA){

  # Matrix of all years/ages  ----
  years <- do.call(seq, as.list(range(data$year)))
  proj_years <- (max(years) + 1):(max(years) + n_proj_years)
  ages <- do.call(seq, as.list(range(data$age)))

  years_ages <- expand.grid(c(years, proj_years), ages)
  colnames(years_ages) <- c("year", "age")

  # Take mean of X-prev years
  mean_weight <- data %>%
    dplyr::filter(year %in% ((max(data$year)-nyrs):max(data$year))) %>%
    dplyr::summarise(pred_weight = mean(weight),.by = age) %>%
    dplyr::arrange(age) %>%
    dplyr::mutate(last_year = last_year,
                  model = paste0("Mean-", nyrs))

  # Predicted weight-at-age matrix ----
  predicted_weight = mean_weight %>%
    dplyr::select(-last_year) %>% # remove last year in case missing year/age combo (same name problem)
    dplyr::full_join(years_ages, by = join_by(age)) %>%
    dplyr::mutate(projection = year %in% proj_years,
                  last_year = last_year,
                  model = paste0("Mean-", nyrs)) %>%
    group_by(year) %>%
    arrange(year, age) %>%
    mutate(ip.value = spline(age, pred_weight, n=n())$y,                           # Spline interpolation of missing ages
           pred_weight = ifelse(is.na(pred_weight), ip.value, pred_weight)) %>%
    dplyr::select(model, year, age, pred_weight, last_year, projection) %>%
    dplyr::arrange(year, age) %>%
    as.data.frame()

  # Predicted observations ----
  data <-  predicted_weight %>%
    dplyr::full_join(data, by = join_by(year, age)) %>%
    arrange(year, age)

  return(list(data = data, model = mean_weight, predicted_weight = predicted_weight))
}
