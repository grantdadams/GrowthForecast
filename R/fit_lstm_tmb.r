# LSTM Cell Function ----
lstm_neuron <- function(x, h_prev, c_prev, W, U, b) {
  act_neuron <- x %*% W + h_prev %*% U # Compute the combined input
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
  logweight = log(weight)
  tiny = 1.0e-6 # Parameter weights


  # Model ----
  for (i in 1:nrow(data_list$mat)) {
    lstm_out <- lstm_neuron(x= data_list$mat[i, ],
                          h_prev=h[i, ], ## output layer
                          c_prev=c[i, ], ## long term memory
                          W,
                          U,
                          b)  # Compute LSTM cell output
    h[i,] <- lstm_out$h ## update stm
    c[i,] <- lstm_out$c ## update cell
  }

  # Final layer
  logoutput <- h %*% last_layer  # Compute the final output
  output <- exp(logoutput)  # Apply exponential to the output

  # Likelihood
  nll <- sum((logweight - logoutput)^2)  # Compute the negative log-likelihood
  nll <- nll + tiny * sum(W^2)  # Add regularization for W
  nll <- nll + tiny * sum(U^2)  # Add regularization for U

  # Report
  RTMB::REPORT(h)  # Report hidden states
  RTMB::REPORT(output)  # Report output
  RTMB::REPORT(nll)  # Report negative log-likelihood

  return(nll)  # Return the negative log-likelihood

}

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
fit_lstm_rtmb <- function(data, nhidden_layer = 2, hidden_dim = 4, input_par = NULL){

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
    ## note: the dims are multiplied by 4 to account for the memory "gates", this will not change
    par_list <- list(
      W = matrix(0.1,nrow =  ncol(data_list$mat), ncol = 4 * hidden_dim) , # Initialize weights for input
      U = matrix(0.1,nrow = hidden_dim ,ncol = 4 * hidden_dim),  # Initialize weights for hidden state
      h = matrix(0.1, nrow = nrow(data_list$mat), ncol = hidden_dim),  # Initialize hidden states
      c = matrix(0.1, nrow = nrow(data_list$mat), ncol = hidden_dim),  # Initialize cell states
      last_layer =  matrix(0.1, hidden_dim, ncol = 1)
    )
  } else{
    par_list <- input_par
  }



  # Build and fit ----
  cmb <- function(f, d) function(p) f(p, d) ## Helper to make closure
  obj <- RTMB::MakeADFun(cmb(lstm_fun_rtmb, data_list), par_list, silent = FALSE)

  if(is.null(input_par)){
    fit <- nlminb(obj$par, obj$fn, obj$gr)
    # fit <- nlminb(obj$par, obj$fn, obj$gr,control=list(eval.max=200000, iter.max=100000, trace=0))
  }else{
    fit <- NULL
  }
  report <- obj$report(obj$env$last.par.best)
  par_list <- obj$env$parList(obj$env$last.par.best)

  # Return ----
  return(list(obj = obj, data = data, fit = fit, report = report,
              parList = par_list, input_par = input_par))
}


# test <- fit_lstm_rtmb(data, input_par = NULL)
# data$predicted <- test$report$output
# ggplot(data, aes(x = age, y = weight, colour = year)) +
#   geom_point(size = 2) +
#   scale_color_discrete() +
#   geom_line(aes(y = predicted))


