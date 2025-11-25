#' von Bert in RTMB with random effects on Winf and Kappa
#'
#' @param pars list of von bertalannfy parameters
#' @param data_list data list
#'
#' @return
#' @export
#'
#' @examples
vbgf_re <- function(pars, data_list){
  require(RTMB)
  RTMB::getAll(pars, data_list)

  # Parameter transform
  # - Mean
  Winf =  exp(logWinf)  # mean asymptoptic weight
  Kappa = exp(logKappa) # mean growth rate
  sigma = exp(logSigma)
  sigmaWinf = exp(logSigmaWinf)   # SD of annual Winf
  sigmaKappa = exp(logSigmaKappa) # SD of annual Kappa

  # - Annual
  Winf_y =  exp(logWinf + logWinf_dev)   # annual asymptoptic weight
  Kappa_y = exp(logKappa + logKappa_dev) # annual growth rate

  # Output
  nll = 0
  WeightPred = rep(0, length(weight))

  # Random effects likelihood
  for(yr in 1:nyrs){
    nll = nll - dnorm(logWinf_dev[yr], 0, sigmaWinf, TRUE)
    nll = nll - dnorm(logKappa_dev[yr], 0, sigmaKappa, TRUE)
  }

  # Model
  for(i in 1:length(weight)){
    tmp = Winf_y[year[i]] * (1.0 - exp(- Kappa_y[year[i]] *(age[i] - t0)))
    WeightPred[i] = tmp

    # Likelihood
    if(!is.na(weight[i])){
      nll = nll - dnorm(log(weight[i]), log(tmp), sigma, TRUE) * weights[i]
    }
  }

  # Report ----
  RTMB::REPORT(sigma)
  RTMB::REPORT(Winf)
  RTMB::REPORT(Kappa)
  RTMB::REPORT(Winf_y)
  RTMB::REPORT(Kappa_y)
  RTMB::REPORT(WeightPred)

  return(nll)
}



#' Fit 3-parameter von-bertalanfy growth function with random effects on Winf and Kappa
#'
#' @param data
#' @param n_proj_years number of years to project forward (default = 2)
#'
#' @return
#' @export
#'
#' @examples
FitVBGF_RE <- function(data,
                    weights=NULL,
                    n_proj_years = 2,
                    last_year = NA){

  # Data ranges ----
  years <- do.call(seq, as.list(range(data$year)))
  proj_years <- (max(years) + 1):(max(years) + n_proj_years)
  ages <- do.call(seq, as.list(range(data$age)))

  # - Projections and fill missing years
  years_ages <- expand.grid(c(years, proj_years), ages)
  colnames(years_ages) <- c("year", "age")

  # - Deal with weights
  if(is.null(weights)){
    weights = rep(1, nrow(data))
  }
  data$weights <- weights


  # Reformat data ----
  data <- data %>%
    dplyr::select(year, age, weight, weights) %>%
    dplyr::full_join(years_ages) %>%
    dplyr::arrange(year, age)

  data_list <- list(
    age        = data$age,
    weight     = data$weight,
    weights    = data$weights,
    year = factor(data$year),
    nyrs = length(unique(data$year))
  )

  # Parameters ----
  par_list <- list(
    logWinf   = log(max(data$weight, na.rm = TRUE)),
    logKappa   = log(0.3),
    logWinf_dev   = rep(0, data_list$nyrs),
    logKappa_dev   = rep(0, data_list$nyrs),
    t0  = -0.25,
    logSigma = -0.5,
    logSigmaWinf = -0.5,
    logSigmaKappa = -0.5
  )

  # Build and fit ----
  cmb <- function(f, d) function(p) f(p, d) ## Helper to make closure
  obj <- RTMB::MakeADFun(cmb(vbgf_re, data_list), par_list, silent = TRUE)
  fit <- nlminb(obj$par, obj$fn, obj$gr,
                control=list(eval.max=200000, iter.max=100000, trace=0))
  report <- obj$report(obj$env$last.par.best)

  # Predicted observations ----
  data <- data %>%
    dplyr::mutate(pred_weight = report$WeightPred,
                  model = "VBGF-RE",
                  projection = year %in% proj_years,
                  last_year = last_year) %>%
    as.data.frame()

  # Predicted weight-at-age matrix ----
  predicted_weight <- data %>%
    dplyr::group_by(year, age) %>%
    dplyr::slice(1) %>%
    dplyr::select(model, year, age, pred_weight, last_year, projection) %>%
    dplyr::arrange(year, age) %>%
    as.data.frame()

  # Return ----
  return(list(obj = obj, data = data, fit = fit, report = report, predicted_weight = predicted_weight))
}

