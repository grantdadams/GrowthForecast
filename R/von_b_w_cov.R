vonBEE <- function(pars){
  require(RTMB)
  RTMB::getAll(pars, data)

  # Parameter transform ----
  sigmaProc = exp(logSigmaProc)
  sigmaObs = exp(logSigmaObs)
  H = exp(logH)
  K = exp(logK)

  # Data transform ----
  logWeightObs = log(weight)

  # Model ----
  x = mu_d + u_y(YrIndex) + beta_0(cohortYr) + covars
  d          = ( 1 - d_offset )/( 1 + exp(-x) )
  Winf       = pow((H/K),1.0/(1.0 - d) )
  logWhat    = log(Winf) + (1.0/(1.0 - d))*log(1.0 - exp(-K * (1.0 - d) * (age - t0)))
  What       = exp(logWhat)


  # Likelihood ----
  nll_obs = -sum(dnorm(logWobs, logWhat, sigma_obs, TRUE))
  nllProc = -sum(dnorm(u_y, m, sigmaProc, TRUE))
  nll = nll_obs + nll_proc;
  return(nll)
}
