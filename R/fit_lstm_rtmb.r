#' LSTM Cell Function
#'
#' @param x
#' @param h_prev
#' @param c_prev
#' @param input_weights
#' @param hidden_weights
#' @param bias_vector
#'
#' @export
lstm_neuron <- function(x, h_prev, c_prev, input_weights, hidden_weights, bias_vector) {

  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  "diag<-" <- ADoverload("diag<-")

  act_neuron <- x %*% input_weights + h_prev %*% hidden_weights + bias_vector # Compute the combined input
  act_neuron_sigmoid <- 1 / (1 + exp(-act_neuron[, 1:(ncol(input_weights)/4)]))  # sig transform entire neuron

  forget_gate <- 1 / (1 + exp(-act_neuron[, ((ncol(input_weights)/4) + 1):(2 * (ncol(input_weights)/4))]))  # Compute the forget gate
  output_gate <- 1 / (1 + exp(-act_neuron[, ((2 * (ncol(input_weights)/4)) + 1):(3 * (ncol(input_weights)/4))]))  # Compute the output gate
  candidate <- tanh(act_neuron[, ((3 * (ncol(input_weights)/4)) + 1):(4 * (ncol(input_weights)/4))])  # Compute the candidate cell state

  input_gate <- forget_gate * c_prev + act_neuron_sigmoid * candidate  # Update the cell state
  h_update <- output_gate * tanh(input_gate)  # Compute the updated hidden state
  return(list(hidden = h_update, cell = input_gate))  # Return the hidden and cell states
}


#' Function to fit LSTM in RTMB
#'
#' @param pars
#' @param data_list
#'
#' @export
lstm_fun_rtmb <- function(pars, data_list) {

  require(RTMB)
  RTMB::getAll(pars, data_list)


  # "[<-" <- ADoverload("[<-") # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  "diag<-" <- ADoverload("diag<-")

  # Data transform ----
  nll = 0
  tiny = 1.0e-6 # Parameter values
  output <- matrix(1e-4, nrow =  dim(hidden)[1], ncol = dim(last_layer)[2]) ## initialize output matrix

  # Process data sequentially by timestep
  for (y in 1:dim(hidden)[1]) {

    if (y == 1) {
      h_prev = matrix(0.1, nrow = 1, ncol =  dim(hidden)[2])  # Hidden states
      c_prev = matrix(0.1, nrow = 1, ncol = dim(hidden)[2]) # Cell states
    }

    x_t <- data_list[[y]] ## pull out the matrix for the current timestep
    # x_t[,'age'] <- (x_t[,'age'] - mean(x_t[,'age'])) / sd(x_t[,'age']) ## normalize
    # Call the LSTM neuron function
    lstm_out <- lstm_neuron(
      x = t(x_t[,'age']), # Input features for the current observation
      h_prev = h_prev,  # Previous hidden state
      c_prev = c_prev,  # Previous cell state
      input_weights = input_weights,          # Weight matrix for input
      hidden_weights = hidden_weights,          # Weight matrix for hidden state
      bias_vector = bias_vector           # Bias vector
    )

    # Store the outputs for the current observation
    h_prev <- lstm_out$hidden  # Updated hidden state
    c_prev <- lstm_out$cell  # Updated cell state

    # Store the hidden and cell states for the current timestep
    hidden[y, ] <- h_prev
    cell[y, ] <- c_prev

    ## dense layer for output, unique age-year combos
    output[y, ] <- hidden[y, ] %*% last_layer

    ## update likelihood
    # - Likelihood

    for (i in 1:length(output[y,])) {
      logvalue <- log(x_t[i,'value'])
      if(is.na(logvalue)) next()
      nll <- nll + (logvalue - log(output[y,i]))^2
    }
    nll <- nll + tiny * sum(input_weights^2)  # Add regularization for input_weights
    nll <- nll + tiny * sum(hidden_weights^2)  # Add regularization for hidden_weights

  } ## end timesteps

  # Report
  RTMB::REPORT(hidden)  # Report hidden states
  RTMB::REPORT(output)  # Report output
  RTMB::REPORT(nll)  # Report negative log-likelihood

  return(nll)  # Return the negative log-likelihood

}

#' Function to fit LSTM recurrent neural net in RTMB
#'
#' @param data data.frame with value, year, and age
#' @param nhidden_layer number of hidden layers (this is on a per obs basis)
#' @param hidden_dim dimension of hidden layers (nhidden_layer + 2 for in/out)
#'
#' @export
fit_lstm_rtmb <- function(data,
                          nhidden_layer = 3,
                          hidden_dim = 5,
                          input_par = NULL,
                          n_proj_years = 2,
                          last_year = NA){

  # - Projection data
  years <- do.call(seq, as.list(range(data$year)))
  proj_years <- (max(years) + 1):(max(years) + n_proj_years)
  ages <- do.call(seq, as.list(range(data$age)))

  # - Projections and fill missing years
  years_ages <- expand.grid(c(years, proj_years), ages)
  colnames(years_ages) <- c("year", "age")
  years_ages$weight = NA
  timesteps <- length(unique(c(years, proj_years)))

  # - Expand for missing years & ages
  data <- data %>% dplyr::full_join(years_ages)

  # - Create list of matrices with imputed data (mean)
  unique_years <- sort(unique(data$year))  # Get unique years

  data_list <- lapply(unique_years, function(yr) {
    # Filter data for the year
    year_data <- data %>%
      dplyr::filter(year == yr) %>%
      dplyr::group_by(age) %>%  # Group by age
      dplyr::summarize(weight = mean(weight, na.rm = TRUE), .groups = "drop")  # Calculate mean weight

    # Ensure missing ages are included with NA weights
    all_ages <- data.frame(age = ages)  # Create a complete list of ages
    year_data <- dplyr::full_join(all_ages, year_data, by = "age")  # Merge with year_data

    # Create the matrix for the current year
    matrix(
      c(year_data$age, year_data$weight),  # Combine age and weight columns
      ncol = 2,  # Two columns: age and weight
      dimnames = list(NULL, c("age", "value"))  # Add column names
    )
  })
  names(data_list) <- unique_years


  # Parameters ----
  if(is.null(input_par)){
    ## note: the dims are multiplied by 4 to account for the memory "gates", this will not change
    par_list <- list(
      input_weights = matrix(0.1, nrow =  length(ages), ncol = 4 * hidden_dim) , # Initialize values for input
      hidden_weights = matrix(0.1, nrow = hidden_dim ,ncol = 4 * hidden_dim),  # Initialize values for hidden state
      hidden = matrix(0.1, nrow = timesteps, ncol = hidden_dim),  # Hidden states
      cell = matrix(0.1, nrow = timesteps, ncol = hidden_dim), # Cell states
      bias_vector = matrix(0.1, nrow = 1, ncol = 4 * hidden_dim),  # Bias vector
      last_layer =  matrix(0.1, nrow = hidden_dim, ncol = length(ages))
    )
  } else{
    par_list <- input_par
  }



  # Build and fit ----
  cmb <- function(f, d) function(p) f(p, d) ## Helper to make closure
  obj <- RTMB::MakeADFun(cmb(lstm_fun_rtmb, data_list), par_list, silent = TRUE)
  gc()

  fit <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e7))

  report <- obj$report(obj$env$last.par.best)
  dimnames(report$output) <- list(year = c(years, proj_years), age = ages)
  par_list <- obj$env$parList(obj$env$last.par.best)

  # Prediction ----
  pred_value <- reshape2::melt(report$output)  %>%
    dplyr::rename(pred_value = value) %>%
    merge(., data, by = c("year", "age")) %>%
    dplyr::mutate(model = "lstm",
           projection = year %in% proj_years,
           last_year = last_year) %>%
    as.data.frame()%>%
    arrange(year, age) %>%
    dplyr::select(model, last_year, year, age, obs_weight = weight, pred_weight = pred_value, projection)

  # plot(pred_value$obs_value, pred_value$pred_value); abline(0,1,col = 'red')
  #
  # tt <- pred_value %>% summarise(meanpred =mean(pred_value, na.rm = TRUE),
  #                                meanobs = mean(obs_value, na.rm = TRUE), .by = c(year, age))
  # plot(tt$meanobs, tt$meanpred); abline(0,1,col = 'red')
  #
  # ggplot(pred_value, aes(x = age,)) +
  #   geom_point(aes(y = obs_value)) +
  #   geom_line(aes(y = pred_value), col = 'blue') +
  #   facet_grid(~year)


  # Return ----
  return(list(obj = obj, data = data, fit = fit, report = report, prediction = pred_value))
}

