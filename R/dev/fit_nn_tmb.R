nn_fun <- function(pars){
  require(RTMB)
  RTMB::getAll(pars, data_list)

  # Parameter transform ----
  sigmaObs = exp(log_sigma_obs)


  # Data transform ----
  logweight = log(weight)
  neurons <- list()


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
  nll = -sum((logweight-logoutput)^2)

  # Report ----
  RTMB::REPORT(sigmaObs)
  RTMB::REPORT(neurons)
  RTMB::REPORT(output)
  RTMB::REPORT(nll)

  return(nll)
}



fit_nn <- function(dat, nhidden_layer = 3, hidden_dim = 5){

  # - Rearrange data
  nlayer = nhidden_layer + 2
  nnform = formula(~age+year)
  data_list <- list(
    weight = data$weight,
    mat = model.matrix(nnform, data),
    nlayer = nlayer
  )


  # Parameters ----
  par_list <- list(
    layer1   = matrix(0, ncol(data_list$mat), hdim),
    hidden  = array(0, dim = c(hdim+1, hdim, nlayer)),
    last_layer = matrix(0, hdim+1, 1),
    log_sigma_obs   = -.3
  )

  library(abind)
  # nn <- neuralnet(log(weight) ~ age+year,
  #                 data = data %>%
  #                   filter(year <= ngroup_hind),
  #                 hidden = c(5,5,5,5),
  #                 linear.output = TRUE,
  #                 stepmax = 1e6,
  #                 lifesign = 'minimal',
  #                 rep=1)
  par_list2 <- list(
    layer1   = nn$weights[[1]][[1]],
    hidden  = abind(nn$weights[[1]][2:4], along = 3),
    last_layer = nn$weights[[1]][[5]],
    log_sigma_obs   = -.3
  )


  # Build and fit ----
  obj <- MakeADFun(nn_fun, par_list, silent = FALSE)
  fit <- optim(par = obj$par,
               fn = obj$fn,
               gr = obj$gr,
               control = list(maxit = 1e6))
  report <- obj$report(obj$env$last.par.best)


  parListFit <- obj$env$parList(obj$env$last.par.best)
  plot(data$age, report$output)
  plot(data$weight, report$output)

  plot(data$weight, exp(predict(nn, newdata = data)))

  # Return ----
  return(list(obj = obj, data = dat, fit = fit, report = report))
}

