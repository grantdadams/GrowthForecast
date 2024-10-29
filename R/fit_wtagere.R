WtAgeRE <- function(pars){
  require(RTMB)
  RTMB::getAll(pars, data_list)

  # Parameter transform ----
  L1 = exp(log_L1)
  L2 = exp(log_L2)
  sigma_coh = exp(log_sd_coh)
  sigma_yr = exp(log_sd_yr )
  K = exp(log_K)
  alpha = exp(log_alpha)


  # Data transform ----
  nages = length(ages)
  nyrs = length(years)


  # Model ----
  wt_hat = X_at
  mnwt    = alpha * (L1 + L2 * (1.0 - (K^(ages-ages[1]))) / (1.0-(K^(ages[nages]-1))))^3
  wt_inc = mnwt[2:nages] - mnwt[1:(nages-1)]

  # -  Initialize first year
  wt_hat[,1] = mnwt

  # - Subsequent years
  for(i in 2:nyrs){
    wt_hat[1, i] = mnwt[1]*exp(c_eta[i])
    wt_hat[2:nages, i] = (wt_hat[1:(nages-1), i-1] + wt_inc * exp(y_eta[i]))
  }


  # Likelihood ----
  # - Observation error
  nll = 0
  for (i in 1:nyrs){ # Ages
    for (j in 1:nages){ # Years
      if(!is.na(X_at[j,i])){
        if(Xsd_at[j,i] > 0){
          nll = nll - dnorm(log(X_at[j,i]), log(mnwt[j]), Xsd_at[j,i], TRUE)
          nll = nll - dnorm(log(X_at[j,i]), log(wt_hat[j,i]), Xsd_at[j,i], TRUE)
        }
      }
    }
  }

  # - Random effects
  nll = nll - sum(dnorm(c_eta, 0, sigma_coh, TRUE))
  nll = nll - sum(dnorm(y_eta, 0, sigma_yr, TRUE))


  # Report ----
  RTMB::REPORT(wt_hat)
  RTMB::REPORT(mnwt)
  RTMB::REPORT(sigma_coh)
  RTMB::REPORT(sigma_yr)
  RTMB::REPORT(K)
  RTMB::REPORT(alpha)
  RTMB::REPORT(L1)
  RTMB::REPORT(L2)

  return(nll)
}



#' Function to fit weight-at-age model from Jim Ianelli
#'
#' @param data data.frame with columns "year", "age", and "weight". Weight should be in kg so that uncertainty calculations do not cause an issue.
#' @param weights likelihood weights for model.
#' @param n_proj_years the number of years to project forward
#'
#' @return a list models
#' @description
#' https://github.com/afsc-assessments/WtAgeRe/tree/main
#'
#' @export
#'
FitWtAgeRE <- function(
    data = NULL,
    weights=NULL,
    # - Number of projection years
    n_proj_years = 2
){

  # Reformat data for GMRF ----
  colnames(data) <- tolower(colnames(data))

  years <- do.call(seq, as.list(range(data$year)))
  proj_years <- (max(years) + 1):(max(years) + n_proj_years)
  ages <- do.call(seq, as.list(range(data$age)))
  years_ages <- expand.grid(years, ages)
  colnames(years_ages) <- c("year", "age")

  gmrf_data <- data %>%
    mutate(weight = weight/1000) %>%
    dplyr::group_by(year, age) %>%
    dplyr::summarise(mn_weight = #ifelse(!is.null(weights), weighted.mean(weight, weights, na.rm = TRUE),
                       mean(weight, na.rm = TRUE), #),
                     sd = # ifelse(!is.null(weights), sqrt(sum(weights * (weight - mn_weight)^2, na.rm = TRUE)),
                       sd(weight, na.rm = TRUE), #),
                     n = n()
    )

  gmrf_mn <- gmrf_data %>%
    dplyr::filter(!is.na(sd)) %>%
    dplyr::select(year, age, mn_weight) %>%
    dplyr::full_join(years_ages) %>%
    dplyr::arrange(year, age) %>%
    tidyr::pivot_wider(names_from = age, values_from = mn_weight)

  gmrf_sd <- gmrf_data %>%
    dplyr::filter(!is.na(sd)) %>%
    dplyr::select(year, age, sd) %>%
    dplyr::full_join(years_ages) %>%
    dplyr::arrange(year, age) %>%
    tidyr::pivot_wider(names_from = age, values_from = sd)


  # Set up data for TMB ----
  # - Read in data weight at age matrix
  X_at <- t(as.matrix(gmrf_mn[,-1])) # removing first col (year column)

  # - Read in standard deviations for weight at age matrix
  Xse_at <- t(as.matrix(gmrf_sd[,c(-1)])) # removing first col (year column)

  # - Convert to CV and sd in lognormal space
  Xcv_at <- sqrt( (exp(Xse_at^2) - 1) )
  Xsd_at <- sqrt((log((Xcv_at)^2 + 1))/(log(10)^2))

  # - Create projection columns
  proj_cols <- matrix(NA, nrow = length(ages), ncol = n_proj_years)

  # - Append NA for projection year
  X_at <- cbind(X_at, proj_cols)
  Xsd_at <- cbind(Xsd_at, proj_cols)


  # Set up TMB inputs  ----
  data_list <- list( years = c(years, proj_years),
                     ages = ages,
                     X_at = X_at,
                     Xsd_at = Xsd_at
  )

  # # Set up TMB inputs  ----
  # data_list <- list( years = meancod$Year,
  #                    ages = 0:10,
  #                    X_at = t(as.matrix(meancod[,-1])),
  #                    Xsd_at = t(as.matrix(stdcod[,-1]))
  # )

  par_list <- list(
    y_eta = rep(0, length(data_list$years)),
    c_eta = rep(0, length(data_list$years)),
    log_K = log(0.3),
    log_alpha = log(1),
    log_L1 = log(10),
    log_L2 = log(80),
    log_sd_coh = 0,
    log_sd_yr = 0
  )


  # Build and fit ----
  obj <- MakeADFun(WtAgeRE, par_list, silent = TRUE, random = c("y_eta", "c_eta"))
  fit <- optim(par = obj$par,
               fn = obj$fn,
               gr = obj$gr,
               control = list(maxit = 1e6))
  report <- obj$report(obj$env$last.par.best)


  # Return ----
  return(list(obj = obj, data = dat, fit = fit, report = report))
}
