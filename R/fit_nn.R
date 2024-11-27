nn_fun <- function(pars){
  require(RTMB)
  RTMB::getAll(pars, data_list)

  # Parameter transform ----
  sigmaObs = exp(log_sigma_obs)

  # Data transform ----
  logweight = (weight)

  # Model ----
  hidden = 1/(1 + exp(-mat %*% layer1))
  for(i in 1:nlayer){
    hidden = 1/(1 + exp(-(hidden %*% hidden_layer[,,1])))
  }

  logoutput = hidden %*% last_layer
  output = (logoutput)

  # Likelihood ----
  nll = -sum(dnorm(logweight, logoutput, sigmaObs, TRUE))

  # Report ----
  RTMB::REPORT(sigmaObs)
  RTMB::REPORT(output)

  return(nll)
}



fit_nn <- function(dat){

  # - Rearrange data
  hdim = 5
  nlayer = 5
  nnform = formula(~-1+age+year)
  data_list <- list(
    weight = data$weight,
    mat = model.matrix(nnform, data),
    nlayer = nlayer
  )


  # Parameters ----
  par_list <- list(
    layer1   = matrix(0, ncol(data_list$mat), hdim),
    hidden_layer  = array(0, dim = c(hdim, hdim, nlayer)),
    last_layer = matrix(0, hdim, 1),
    log_sigma_obs   = -.3
  )

  # Build and fit ----
  obj <- MakeADFun(nn_fun, par_list, silent = FALSE)
  fit <- optim(par = obj$par,
               fn = obj$fn,
               gr = obj$gr,
               control = list(maxit = 1e6))
  report <- obj$report(obj$env$last.par.best)
  plot(data$age, report$output)
  plot(data$weight, report$output)

  # Return ----
  return(list(obj = obj, data = dat, fit = fit, report = report))
}

