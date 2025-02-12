WtAgeRE <- function(pars, data_list){
  require(RTMB)
  RTMB::getAll(pars)

  # Parameter transform ----
  L1 = exp(log_L1)
  L2 = exp(log_L2)
  sigma_coh = exp(log_sd_coh)
  sigma_yr = exp(log_sd_yr )
  K = exp(log_K)
  alpha = exp(log_alpha)


  # Data transform ----
  years = data_list$years
  ages = data_list$ages
  X_at = data_list$X_at
  Xsd_at = data_list$Xsd_at
  nages = length(ages)
  nyrs = length(years)

  # Model ----
  wt_hat = X_at
  mnwt    = alpha * (L1 + (L2 - L1) * (1.0 - (K^(ages-ages[1]))) / (1.0-(K^(ages[nages]-1))))^3
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
          nll = nll - dnorm((X_at[j,i]), (mnwt[j]), Xsd_at[j,i], TRUE) # Mean
          nll = nll - dnorm((X_at[j,i]), (wt_hat[j,i]), Xsd_at[j,i], TRUE) # Annual deviation
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
  years <- do.call(seq, as.list(range(data$year)))
  proj_years <- (max(years) + 1):(max(years) + n_proj_years)
  ages <- do.call(seq, as.list(range(data$age)))
  years_ages <- expand.grid(proj_years, ages)
  colnames(years_ages) <- c("year", "age")

  # Deal with weights
  if(is.null(weights)){
    weights = rep(1, nrow(data))
  }
  data$weights <- weights

  gmrf_data <- data %>%
    dplyr::group_by(year, age) %>%
    dplyr::summarise(mn_weight = weighted.mean(weight, weights, na.rm = TRUE),
                     sd = sqrt(sum(weights * (weight - mn_weight)^2, na.rm = TRUE)),
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

  par_list <- list(
    y_eta = rep(0, length(data_list$years)),
    c_eta = rep(0, length(data_list$years)),
    log_K = log(0.3),
    log_alpha = -11, # Fixed
    log_L1 = log(10),
    log_L2 = log(80),
    log_sd_coh = 0,
    log_sd_yr = 0
  )



  # Build and fit ----
  cmb <- function(f, d) function(p) f(p, d) ## Helper to make closure
  map <- list(log_alpha = factor(NA))
  obj <- RTMB::MakeADFun(cmb(WtAgeRE, data_list), par_list, silent = TRUE, map = map, random = c("y_eta", "c_eta"))
  fit <- optim(par = obj$par,
               fn = obj$fn,
               gr = obj$gr,
               control = list(maxit = 1e6))
  report <- obj$report(obj$env$last.par.best)
  colnames(report$wt_hat) <- data_list$years
  rownames(report$wt_hat) <- data_list$age

  # Prediction ----
  # - Prediction for each obs
  pred_weight <- reshape2::melt(report$wt_hat)
  colnames(pred_weight) <- c("age", "year", "pred_weight")
  data <- merge(data, pred_weight, all.x = TRUE)

  # - Predicted for forecast
  pred_weight <- data %>%
    dplyr::filter(year %in% proj_years) %>%
    dplyr::select(-weight, -weights) %>%
    mutate(model = "WtAgeRe",
           last_year = last_year) %>%
    as.data.frame()

  # Return ----
  return(list(obj = obj, data = data, fit = fit, report = report, prediction = pred_weight))
}

# # Set up TMB inputs  ----
# # - GOA pollock example
# data <- read.csv("data/wtagere_pk.csv")
# data_list <- list( years = data$Cohort,
#                    ages = 1:10,
#                    X_at = data %>% dplyr::select(contains("Age")) %>% t(),
#                    Xsd_at = data %>% dplyr::select(contains("Std")) %>% t()
# )
#
# par_list <- list(
#   y_eta = data$y_eta * exp(-0.00172266199104) + 0.5 * exp(-0.00172266199104)^2,
#   c_eta = data$c_eta * exp(0.0113343421692) + 0.5 * exp(0.0113343421692)^2,
#   log_K = (-0.245120653667),
#   log_alpha = (-11.0000000000),
#   log_L1 = log(18.0855652966),
#   log_L2 = log(44.2873288952),
#   log_sd_coh = 0.0113343421692,
#   log_sd_yr = -0.00172266199104
# )

# for (int j=age_st;j<=age_end;j++)
# {
#   mnwt(j)    = alpha * pow(L1 + (L2-L1)*(1.-pow(K,double(j-age_st))) / (1.-pow(K,double(nages-1))) ,3);
# }
# wt_inc       = --mnwt(age_st+1,age_end) - mnwt(age_st,age_end-1);
#
# // Initialize first year
# wt_pre(styr)    = mnwt;
#
# // subsequent years
# for (int i=styr+1;i<=endyr;i++)
# {
#   wt_pre(i,age_st) = mnwt(age_st)*mfexp(square(sigma_coh)/2.+sigma_coh*coh_eff(i));
#   if (last_phase())
#     wt_pre(i)(age_st+1,age_end) = ++(wt_pre(i-1)(age_st,age_end-1) + wt_inc*mfexp(square(sigma_yr)/2. + sigma_yr*yr_eff(i)));
#   else
#     wt_pre(i)(age_st+1,age_end) = ++(wt_pre(i-1)(age_st,age_end-1) + wt_inc*mfexp(sigma_yr*yr_eff(i)));
# }
# int iyr;
#
# // Fit global mean to all years...
# for (int h = 1;h<=ndat;h++) // Loop over number of data sets (assuming 1st is target)
# {
#   for (int i=1;i<=nyrs_data(h);i++) // Loop over how many observations w/in a data set
#   {
#     iyr = yrs_data(h,i);
#     if (h>1) // First data set is correct scale, latter sets need a multiplier (d_scale)
#     wt_hat(h,i) = elem_prod(d_scale(h-1) , wt_pre(iyr) );
#     else
#       wt_hat(h,i) = wt_pre(iyr);
#
#     for (int j=age_st;j<=age_end;j++)
#     {
#       // cout<<nll<<endl;
#       nll += square(wt_obs(h,i,j) - mnwt(j))      /(2.*square(sd_obs(h,i,j)));
#       nll += square(wt_obs(h,i,j) - wt_hat(h,i,j))/(2.*square(sd_obs(h,i,j)));
#     }
#   }
# }
# // Random effects are on N(0,1) scale
# nll += 0.5*norm2(coh_eff);
# nll += 0.5*norm2( yr_eff);
#
# if (sd_phase())
# {
#   wt_cur  = wt_pre(cur_yr);
#   wt_next = wt_pre(endyr-1);
#   wt_yraf = wt_pre(endyr);
#   wt_hist = wt_pre;
# }

