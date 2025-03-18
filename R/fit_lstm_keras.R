#' Function to fit LSTM recurrent neural net in RTMB
#'
#' @param data data.frame with weight, year, and age
#' @param nhidden_layer number of hidden layers (this is on a per obs basis)
#' @param hidden_dim dimension of hidden layers (nhidden_layer + 2 for in/out)
#'
#' @return
#' @export
#'
#' @examples
fit_lstm_keras <- function(data,
                          hidden_dim = 5,
                          n_proj_years = 2,
                          last_year = NA){

  # - Settings
  nnform = formula(~age+year)

  # - Projection data
  years <- do.call(seq, as.list(range(data$year)))
  proj_years <- (max(years) + 1):(max(years) + n_proj_years)
  ages <- do.call(seq, as.list(range(data$age)))

  # - Projections and fill missing years
  years_ages <- expand.grid(c(years, proj_years), ages)
  colnames(years_ages) <- c("year", "age")
  years_ages$weight = NA

  # - Data
  X_train <- model.matrix(nnform, data)
  X_test <- model.matrix(nnform, years_ages)
  y_train <- data$weight


  # Define the LSTM model ----
  model <- keras::keras_model_sequential()
  model %>%
    keras::layer_lstm(units = hidden_dim, input_shape = c(nrow(X_test)),
                      activation = "tanh",
                      recurrent_activation = "sigmoid") %>%  # units = number of LSTM cells %>%
    keras::layer_lstm(units = hidden_dim, input_shape = c(nrow(X_test)),
                      activation = "tanh",
                      recurrent_activation = "sigmoid") %>%  # units = number of LSTM cells
    layer_dense(units = 1, activation = 'linear')

  # Compile the model
  model %>% compile(
    loss = 'mse',  # Mean squared error
    optimizer = optimizer_adam(),  # Optimizer like Adam
  )
  summary(model)

  # Train the model
  history <- model %>% fit(
    X_train, y_train,
    epochs = 10,  # Number of training epochs
    batch_size = 32,  # Batch size
    validation_split = 0.2  # Percentage of data to use for validation
  )


  # Prediction ----
  # - Prediction for each obs
  data$pred_weight = model %>% predict(X_train)
  years_ages$pred_weight = model %>% predict(X_test)

  # - Predicted for projection
  pred_weight <- rbind(data, years_ages) %>%
    # dplyr::filter(year %in% proj_years) %>%
    # dplyr::select(-weight) %>%
    mutate(model = "lstm keras",
           projection = year %in% proj_years,
           last_year = last_year) %>%
    as.data.frame()%>%
    arrange(year, age) %>%
    dplyr::select(model, last_year, year, age, obs_weight = weight, pred_weight, projection)


  # Return ----
  return(list(obj = model, data = data, report = report, prediction = pred_weight))
}



