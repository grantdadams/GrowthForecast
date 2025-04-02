# LSTM Cell Function ----
lstm_neuron <- function(x, h_prev, c_prev, W, U, b) {
  act_neuron <- x %*% W + h_prev %*% U + b # Compute the combined input
  act_neuron_sigmoid <- 1 / (1 + exp(-act_neuron[, 1:(ncol(W)/4)]))  # sig transform entire neuron

  forget_gate <- 1 / (1 + exp(-act_neuron[, ((ncol(W)/4) + 1):(2 * (ncol(W)/4))]))  # Compute the forget gate
  output_gate <- 1 / (1 + exp(-act_neuron[, ((2 * (ncol(W)/4)) + 1):(3 * (ncol(W)/4))]))  # Compute the output gate
  candidate <- tanh(act_neuron[, ((3 * (ncol(W)/4)) + 1):(4 * (ncol(W)/4))])  # Compute the candidate cell state

  input_gate <- forget_gate * c_prev + act_neuron_sigmoid * candidate  # Update the cell state
  h_update <- output_gate * tanh(input_gate)  # Compute the updated hidden state
  return(list(h = h_update, c = input_gate))  # Return the hidden and cell states
}

# Function to fit LSTM in RTMB ----
lstm_fun_rtmb <- function(pars, data_list) {

  require(RTMB)
  RTMB::getAll(pars, data_list)

  # Data transform ----
  nll = 0
  tiny = 1.0e-6 # Parameter values
  output <- matrix(1e-4, nrow =  dim(h)[1], ncol = dim(last_layer)[2]) ## initialize output matrix

  # Process data sequentially by timestep
  for (y in 1:dim(h)[1]) {

    if (y == 1) {
      h_prev = matrix(0.1, nrow = 1, ncol =  dim(h)[2])  # Hidden states
      c_prev = matrix(0.1, nrow = 1, ncol = dim(h)[2]) # Cell states
    }

    x_t <- data_list[[y]] ## pull out the matrix for the current timestep
    # x_t[,'age'] <- (x_t[,'age'] - mean(x_t[,'age'])) / sd(x_t[,'age']) ## normalize
    # Call the LSTM neuron function
    lstm_out <- lstm_neuron(
      x = t(x_t[,'age']), # Input features for the current observation
      h_prev = h_prev,  # Previous hidden state
      c_prev = c_prev,  # Previous cell state
      W = W,          # Weight matrix for input
      U = U,          # Weight matrix for hidden state
      b = b           # Bias vector
    )

    # Store the outputs for the current observation
    h_prev <- lstm_out$h  # Updated hidden state
    c_prev <- lstm_out$c  # Updated cell state

    # Store the hidden and cell states for the current timestep
    h[y, ] <- h_prev
    c[y, ] <- c_prev

    ## dense layer for output, unique age-year combos
    output[y, ] <- h[y, ] %*% last_layer

    ## update likelihood
    # - Likelihood

    for (i in 1:length(output[y,])) {
      logvalue <- log(x_t[i,'value'])
      if(is.na(logvalue)) next()
      nll <- nll + (logvalue - log(output[y,i]))^2
    }
    nll <- nll + tiny * sum(W^2)  # Add regularization for W
    nll <- nll + tiny * sum(U^2)  # Add regularization for U

  } ## end timesteps

  # Report
  RTMB::REPORT(h)  # Report hidden states
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
#' @return
#' @export
#'
#' @examples
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
      W = matrix(0.1, nrow =  length(ages), ncol = 4 * hidden_dim) , # Initialize values for input
      U = matrix(0.1, nrow = hidden_dim ,ncol = 4 * hidden_dim),  # Initialize values for hidden state
      h = matrix(0.1, nrow = timesteps, ncol = hidden_dim),  # Hidden states
      c = matrix(0.1, nrow = timesteps, ncol = hidden_dim), # Cell states
      b = matrix(0.1, nrow = 1, ncol = 4 * hidden_dim),  # Bias vector
      last_layer =  matrix(0.1, nrow = hidden_dim, ncol = length(ages))
    )
  } else{
    par_list <- input_par
  }



  # Build and fit ----
  cmb <- function(f, d) function(p) f(p, d) ## Helper to make closure
  obj <- RTMB::MakeADFun(cmb(lstm_fun_rtmb, data_list), par_list, silent = FALSE)
  gc()

  fit <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e7))

  report <- obj$report(obj$env$last.par.best)
  par_list <- obj$env$parList(obj$env$last.par.best)

  # Prediction ----
  pred_value <- reshape2::melt(report$output)  %>%
    select(year = Var1, age = Var2, pred_value = value) %>%
    merge(., data, by = c("year", "age")) %>%
    mutate(model = "lstm",
           projection = year %in% proj_years,
           last_year = last_year) %>%
    as.data.frame()%>%
    arrange(year, age) %>%
    dplyr::select(model, last_year, year, age, obs_value = weight, pred_value, projection)

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

