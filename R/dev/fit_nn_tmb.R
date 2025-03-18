nn_fun_rtmb <- function(pars, data_list){
  require(RTMB)
  RTMB::getAll(pars, data_list)

  # Parameter transform ----
  # sigmaObs = exp(log_sigma_obs)


  # Data transform ----
  logweight = log(weight)
  neurons <- list()
  tiny = 1.0e-6 # Parameter weights


  # Model ----
  # - Layer 1
  neurons[[1]] <- mat
  nueron <- neurons[[1]] %*% layer1
  act.nueron <- 1/(1+exp(-nueron))

  # - Hidden layers
  for(i in 2:(nlayer-1)){
    neurons[[i]] = cbind(1, act.nueron)
    nueron <- neurons[[i]] %*% hidden[,,i-1]
    act.nueron = 1/(1 + exp(-nueron))
  }


  # - Final layer
  neurons[[nlayer]] = cbind(1, act.nueron)
  logoutput = neurons[[nlayer]] %*% last_layer # Linear output activation function
  output = exp(logoutput)


  # Likelihood ----
  nll = sum((logweight-logoutput)^2) # Sum of squares, dnorm may be worth testing


  nll = nll + tiny * sum(layer1^2)
  nll = nll + tiny * sum(last_layer^2)
  nll = nll + tiny * sum(hidden^2)


  # Report ----
  # RTMB::REPORT(sigmaObs)
  RTMB::REPORT(neurons)
  RTMB::REPORT(output)
  RTMB::REPORT(nll)

  return(nll)
}



#' Function to fit neural net in RTMB
#'
#' @param data data.frame with weight, year, and age
#' @param nhidden_layer number of hidden layers
#' @param hidden_dim dimension of hidden layers
#'
#' @return
#' @export
#'
#' @examples
fit_nn_rtmb <- function(data, nhidden_layer = 3, hidden_dim = 5, input_par = NULL){

  # - Rearrange data
  nlayer = nhidden_layer + 2
  nnform = formula(~age+year) #TODO adjust to use
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
      last_layer = matrix(0, hidden_dim+1, 1) # ,
      # log_sigma_obs   = -.3
    )
  } else{
    par_list <- input_par
  }


  # library(abind)
  # nn <- neuralnet(log(weight) ~ age+year,
  #                 data = data %>%
  #                   filter(year <= ngroup_hind),
  #                 hidden = c(5,5,5,5),
  #                 linear.output = TRUE,
  #                 stepmax = 1e6,
  #                 lifesign = 'minimal',
  #                 rep=1)
  # par_list2 <- list(
  #   layer1   = nn$weights[[1]][[1]],
  #   hidden  = abind(nn$weights[[1]][2:4], along = 3),
  #   last_layer = nn$weights[[1]][[5]],
  #   log_sigma_obs   = -.3
  # )


  # Build and fit ----
  cmb <- function(f, d) function(p) f(p, d) ## Helper to make closure
  obj <- RTMB::MakeADFun(cmb(nn_fun_rtmb, data_list), par_list, silent = FALSE)
  if(is.null(input_par)){
    fit <- nlminb(obj$par, obj$fn, obj$gr,
                  control=list(eval.max=200000, iter.max=100000, trace=0))
  }else{
    fit <- NULL
  }
  report <- obj$report(obj$env$last.par.best)
  par_list <- obj$env$parList(obj$env$last.par.best)

  # Return ----
  return(list(obj = obj, data = data, fit = fit, report = report, parList = par_list, input_par = input_par))
}

