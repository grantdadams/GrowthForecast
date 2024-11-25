nn_fun <- function(pars){
  require(RTMB)
  RTMB::getAll(pars, data_list)

  # Parameter transform ----
  sigmaObs = exp(log_sigma_obs)

  # Data transform ----
  logweight = log(weight)

  # Model ----
  hidden = mat %*% layer1 + beta1
  for(i in 1:nlayer){
    hidden = (hidden %*% hidden_layer[,,1]) + hidden_beta[i]
  }

  logoutput = hidden %*% last_layer + last_beta
  output = exp(logoutput)

  # Likelihood ----
  nll = -sum(dnorm(logweight, logoutput, sigmaObs, TRUE))

  # Report ----
  RTMB::REPORT(sigmaObs)
  RTMB::REPORT(output)

  return(nll)
}



fit_nn <- function(dat){

  # - Rearrange data
  hdim = 10
  nlayer = 3
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
    beta1 = 0,
    last_beta = 0,
    hidden_beta = rep(0, nlayer),
    last_layer = matrix(0, hdim, 1),
    log_sigma_obs   = -.3
  )

  # Build and fit ----
  obj <- MakeADFun(nn_fun, par_list, silent = TRUE)
  fit <- optim(par = obj$par,
               fn = obj$fn,
               gr = obj$gr,
               control = list(maxit = 1e6))
  report <- obj$report(obj$env$last.par.best)
  plot(data$age, report$output)

  # Return ----
  return(list(obj = obj, data = dat, fit = fit, report = report))
}

