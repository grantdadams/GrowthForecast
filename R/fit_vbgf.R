vbgf <- function(pars, data_list){
  require(RTMB)
  RTMB::getAll(pars, data_list)

  # Parameter transform
  Winf =  exp(logWinf) # asymptoptic weight
  Kappa = exp(logKappa) # growth rate
  sigma = exp(logSigma)

  # Output
  nll = 0
  WeightPred = rep(0, length(data_list$weight))

  for(i in 1:length(weight)){
    # Model
    tmp = Winf * (1.0 - exp(- Kappa *(age[i] - t0)))
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
  RTMB::REPORT(WeightPred)

  return(nll)
}



#' Fit 3-parameter von-bertalanfy growth function
#'
#' @param data
#' @param n_proj_years number of years to project forward (default = 2)
#'
#' @return
#' @export
#'
#' @examples
FitVBGF <- function(data,
                    weights=NULL,
                    n_proj_years = 2,
                    last_year = NA){

  # Projection ----
  years <- do.call(seq, as.list(range(data$year)))
  proj_years <- (max(years) + 1):(max(years) + n_proj_years)
  ages <- do.call(seq, as.list(range(data$age)))
  pred_df <- expand.grid(proj_years, ages)
  colnames(pred_df) <- c("year", "age")
  pred_df$weight = NA
  pred_df$weights = 1

  # Deal with weights
  if(is.null(weights)){
    weights = rep(1, nrow(data))
  }
  data$weights <- weights


  # Data ----
  data <- data %>%
    dplyr::select(year, age, weight, weights) %>%
    rbind(pred_df)

  data_list <- list(
    age        = data$age,
    weight     = data$weight,
    weights    = data$weights
  )

  # Parameters ----
  par_list <- list(
    logWinf   = log(5),
    logKappa   = log(0.15),
    t0  = -0.25,
    logSigma = -0.5
  )

  # Build and fit ----
  cmb <- function(f, d) function(p) f(p, d) ## Helper to make closure
  obj <- RTMB::MakeADFun(cmb(vbgf, data_list), par_list, silent = TRUE)
  fit <- nlminb(obj$par, obj$fn, obj$gr,
                control=list(eval.max=200000, iter.max=100000, trace=0))
  report <- obj$report(obj$env$last.par.best)

  # Prediction ----
  # - Prediction for each obs
  data$pred_weight = report$WeightPred

  # - Predicted for projection
  pred_weight <- data %>%
    dplyr::filter(year %in% proj_years) %>%
    dplyr::select(-weight, -weights) %>%
    mutate(model = "vbgf",
           last_year = last_year) %>%
    as.data.frame()

  # Return ----
  return(list(obj = obj, data = data, fit = fit, report = report, prediction = pred_weight))
}

