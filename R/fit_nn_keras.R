#' Fit deep neural net in Keras
#'
#' @param form model formula of features to use for neural net
#' @param data data.frame with columns "year", "age", and "weight". Weight should be in kg so that uncertainty calculations do not cause an issue.
#' @param n_proj_years the number of years to project forward
#' @param last_year final year of data
#'
#' @return
#' @export
#'
#' @examples
fit_nn_keras <- function(
    form = formula(~ age * year + as.factor(age)),
    data,
    n_proj_years = 2,
    last_year = NA){
  require(neuralnet)

  # Data ----
  # Matrix of all years/ages  ----
  years <- do.call(seq, as.list(range(data$year)))
  proj_years <- (max(years) + 1):(max(years) + n_proj_years)
  ages <- do.call(seq, as.list(range(data$age)))

  years_ages <- expand.grid(c(years, proj_years), ages)
  colnames(years_ages) <- c("year", "age")

  # - Data
  x_train <- model.matrix(form, data)
  x_test <- model.matrix(form, years_ages) %>%
    as.matrix()
  y_train <- data$weight


  # Model ----
  model <- keras3::keras_model_sequential() %>%
    keras3::layer_dense(units = 256, activation = 'relu') %>%
    keras3::layer_dropout(rate = 0.2) %>%
    keras3::layer_dense(units = 128, activation = 'relu') %>%
    keras3::layer_dropout(rate = 0.2) %>%
    keras3::layer_dense(units = 1, activation = 'linear')

  model %>% compile(
    loss = 'mean_squared_error',
    optimizer = optimizer_adam(learning_rate = 0.001),
    metrics = c('mse')
  )

  model %>% fit(
    x_train, y_train,
    epochs = 50,
    batch_size = 512,
    validation_split = 0.2,
    verbose=0
  )

  model %>% evaluate(x_train, y_train)


  # Prediction ----
  years_ages$pred_weight <- model |> predict(x_test) # Just year x age combo

  # Predicted weight-at-age matrix ----
  predicted_weight <- years_ages %>%
    dplyr::mutate(model = "NN",
                  projection = year %in% proj_years,
                  last_year = last_year) %>%
    dplyr::select(model, year, age, pred_weight, last_year, projection) %>%
    dplyr::arrange(year, age) %>%
    as.data.frame()

  # Predicted observations ----
  data <-  predicted_weight %>%
    dplyr::full_join(data, by = join_by(year, age)) %>%
    dplyr::arrange(year, age)

  # Return ----
  return(list(obj = model, data = data, predicted_weight = predicted_weight))
}

