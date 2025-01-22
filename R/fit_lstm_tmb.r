# LSTM Cell Function ----
ltsm_neuron <- function(x, h_prev, c_prev, W, U, b) {
  z <- x %*% W + h_prev %*% U +  matrix(b, nrow = 1,  byrow = TRUE)  # Compute the combined input  # Compute the combined input
  input_gate <- 1 / (1 + exp(-z[, 1:dim(W)[2]]))  # Compute the input gate
  forget_gate <- 1 / (1 + exp(-z[, ((ncol(W)/4) + 1):(2 * (ncol(W)/4))]))  # Compute the forget gate
  output_gate <-1 / (1 + exp(-z[, ((2 * (ncol(W)/4)) + 1):(3 * (ncol(W)/4))])) # Compute the output gate
  candidate <- tanh(z[, ((3 * (ncol(W)/4)) + 1):(4 * (ncol(W)/4))])  # Compute the candidate cell state
  c <- forget_gate * c_prev + input_gate * candidate  # Update the cell state
  h <- output_gate * tanh(c)  # Compute the hidden state
  return(list(h = h, c = c))  # Return the hidden and cell states
}

# Function to fit LSTM in RTMB ----
lstm_fun_rtmb <- function(data, nhidden_layer = 2, hidden_dim = 5, input_par = NULL) {

  require(RTMB)
  # RTMB::getAll(pars, data_list)

  # Data transform ----
  logweight = log(weight)
  tiny = 1.0e-6 # Parameter weights


  # Model ----
  # Initialize LSTM parameters
  W <- matrix(rnorm(ncol(data_list$mat) * 4 * hidden_dim), ncol = 4 * hidden_dim)  # Initialize weights for input
  U <- matrix(rnorm(hidden_dim * 4 * hidden_dim), ncol = 4 * hidden_dim)  # Initialize weights for hidden state
  b <- rnorm(4 * hidden_dim)  # Initialize biases

  # Forward pass
  h <- matrix(0, nrow = nrow(data_list$mat), ncol = hidden_dim)  # Initialize hidden states
  c <- matrix(0, nrow = nrow(data_list$mat), ncol = hidden_dim)  # Initialize cell states


  for (t in 1:nrow(data_list$mat)) {  # Loop through each time step
    lstm_out <- ltsm_neuron(data_list$mat[t, ],
                          h[t, ],
                          c[t, ], W, U, b)  # Compute LSTM cell output

    # Update cell & hidden states and ensure matrix structure
    h  <- matrix(lstm_out$h, nrow = nrow(data_list$mat), ncol = hidden_dim)
    c  <- matrix(lstm_out$c, nrow = nrow(data_list$mat), ncol = hidden_dim)
  }

  # Final layer
  logoutput <- h %*% matrix(rnorm(hidden_dim), ncol = 1)  # Compute the final output
  output <- exp(logoutput)  # Apply exponential to the output

  # Likelihood
  nll <- sum((log(data_list$weight) - logoutput)^2)  # Compute the negative log-likelihood
  nll <- nll + tiny * sum(W^2)  # Add regularization for W
  nll <- nll + tiny * sum(U^2)  # Add regularization for U
  nll <- nll + tiny * sum(b^2)  # Add regularization for b

  # Report
  RTMB::REPORT(h)  # Report hidden states
  RTMB::REPORT(output)  # Report output
  RTMB::REPORT(nll)  # Report negative log-likelihood

  return(nll)  # Return the negative log-likelihood

}

#' Function to fit LSTM recurrent neural net in RTMB
#'
#' @param data data.frame with weight, year, and age
#' @param nhidden_layer number of hidden layers
#' @param hidden_dim dimension of hidden layers
#'
#' @return
#' @export
#'
#' @examples
fit_lstm_rtmb <- function(data, nhidden_layer = 3, hidden_dim = 5, input_par = NULL){

  # - Rearrange data
  nlayer = nhidden_layer + 2
  nnform = formula(~age+year)
  data_list <- list(
    weight = data$weight,
    mat = model.matrix(nnform, data),
    nlayer = nlayer
  )


  # Parameters ----
  if(is.null(input_par)){
    par_list <- list(
      layer1   = matrix(0, ncol(data_list$mat), hidden_dim),
      hidden  = array(0, dim = c(hidden_dim+1, hidden_dim, nhidden_layer)),
      last_layer = matrix(0, hidden_dim+1, 1)
    )
  } else{
    par_list <- input_par
  }



  # Build and fit ----
  cmb <- function(f, d) function(p) f(p, d) ## Helper to make closure
  obj <- RTMB::MakeADFun(cmb(lstm_fun_rtmb, data_list), par_list, silent = FALSE)

  if(is.null(input_par)){
    fit <- nlminb(obj$par, obj$fn, obj$gr,
                  control=list(eval.max=200000, iter.max=100000, trace=0))
  }else{
    fit <- NULL
  }
  report <- obj$report(obj$env$last.par.best)
  par_list <- obj$env$parList(obj$env$last.par.best)

  # Return ----
  return(list(obj = obj, data = data, fit = fit, report = report,
              parList = par_list, input_par = input_par))
}


test <- fit_lstm_rtmb(data=dat, input_par = NULL)

dat$predicted <- test$report$output

ggplot(dat, aes(x = age, y = weight, colour = year)) +
  geom_point(size = 2) +
  scale_color_discrete() +
  geom_line


