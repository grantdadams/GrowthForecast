# Cheng, M. L., Thorson, J. T., Ianelli, J. N., & Cunningham, C. J. (2023). Unlocking the triad of age, year, and cohort effects for stock assessment: demonstration of a computationally efficient and reproducible framework using weight-at-age. Fisheries Research, 266, 106755.
# https://www.sciencedirect.com/science/article/abs/pii/S0165783623001480
# https://github.com/chengmatt/GMRF_WAA/blob/master/R_scripts/model%20runs/RTMB_GrowthModel.R


#' Title Constructor algorithin for correlations within ages, years, and cohort
#'
#' @param n_ages Number of ages
#' @param n_yrs Number of years
#' @param pcorr_age correlations for age
#' @param pcorr_year correaltions for year
#' @param pcorr_cohort correlaitons for cohort
#' @param ln_var_value log space variance
#' @param Var_Param variance type == 0, marginal (stationary and slower run time), == 1 conditional (non-statationary, faster run time)
#'
#' @returns Sparse precision matrix dimensioned by n_ages * n_years, n_ages * n_years
#' @export
#'
Get_3d_precision <- function(n_ages, n_yrs, pcorr_age, pcorr_year, pcorr_cohort, ln_var_value, Var_Param){

  require(Matrix)
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  "diag<-" <- ADoverload("diag<-")

  index = expand.grid(seq_len(n_ages), seq_len(n_yrs)) # create index combinations to loop through
  i = j = x = numeric(0) # initialize posiiton to fill in precision matrix
  var_value = exp(ln_var_value) # transform to normal space

  for(n in 1:nrow(index)){
    age = index[n,1] # get age index out of all index combinations
    year = index[n,2] # get year index out of all index combinations
    if(age > 1 ){
      i = c(i, n)
      j = c(j, which(index[,1] == (age-1) & index[,2] == year))
      x = c(x, pcorr_year) # year correaltion indexing
    }
    if(year > 1){
      i = c(i, n)
      j = c(j, which(index[,1]==age & index[,2]==(year-1)) )
      x = c(x, pcorr_age) # age correlation indexing
    }
    if( age>1 & year>1 ){
      i = c(i, n)
      j = c(j, which(index[,1]==(age-1) & index[,2] == (year-1)) )
      x = c(x, pcorr_cohort) # cohort correlation indexing
    }
  } # end n loop

  # create B path matrix
  B = matrix(0, nrow = n_ages * n_yrs, ncol = n_ages * n_yrs)
  B[cbind(i, j)] = x
  B = as(B, "sparseMatrix")

  # identity matrix
  I = as(diag(1, n_ages * n_yrs, n_ages * n_yrs), "sparseMatrix")

  if(Var_Param == 0) d = var_value # conditional variance (non-stationary variance)

  # Solve Omega recursively for stationary variance (accumulator function)
  if(Var_Param == 1) {
    L = solve(I-B) # solve to get accumulator function for stationary variance
    d = rep(0, nrow(index))
    for(n in 1:nrow(index) ){
      if(n==1){
        d[n] = var_value
      }else{
        cumvar = sum(L[n,seq_len(n-1)] * d[seq_len(n-1)] * L[n,seq_len(n-1)])
        d[n] = (var_value-cumvar) / L[n,n]^2
      }
    } # end n loop
  } # end marginal variance (stationary variance)

  # omega matrix
  Omega_inv = diag(1/d, n_ages * n_yrs, n_ages * n_yrs)
  Q = as((I-t(B)) %*% Omega_inv %*% (I-B), "sparseMatrix") # solve for precision

  return(Q)
}


#' GMRF (Cheng et al 2023)
#'
#' @param pars
#'
#' @return
#' @export
#'
#' @examples
growth_3d = function(pars, data_list) {

  RTMB::getAll(pars, data_list) # load in starting values and data

  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  "diag<-" <- ADoverload("diag<-")

  # make precision
  Q_sparse = Get_3d_precision(nrow(ln_Y_at), ncol(ln_Y_at), rho_a, rho_y, rho_c, log_sigma2, Var_Param)

  jnLL = 0 # set up jnLL

  # set up containers
  mu_at = array(0, dim = c(nrow(ln_Y_at), ncol(ln_Y_at)))
  eps_at = array(0, dim = c(nrow(ln_Y_at), ncol(ln_Y_at)))

  # transform pars
  L0 = exp(ln_L0);
  Linf = exp(ln_Linf);
  k = exp(ln_k);
  alpha = exp(ln_alpha);
  beta = exp(ln_beta);

  # get parametric form
  for(a in 1:nrow(X_at)) {
    for(t in 1:ncol(X_at)) {
      mu_at[a,t] = Linf - (Linf - L0) * exp(-k * a)
      mu_at[a,t] = alpha * mu_at[a,t]^beta
    }
  }

  # observation likelihood
  for(a in 1:nrow(X_at)) {
    for(t in 1:ncol(X_at)) {
      if(!is.na(X_at[a,t])) {
        if(Xsd_at[a,t] > 0) jnLL = jnLL + -dnorm(log(X_at[a,t]), ln_Y_at[a,t], Xsd_at[a,t], TRUE)
      }
    }
  }

  # process error liklelihood
  eps_at = ln_Y_at - log(mu_at)
  jnLL = jnLL - dgmrf(x = as.vector(eps_at), mu = 0, Q = Q_sparse, TRUE)

  RTMB::REPORT(jnLL);
  RTMB::REPORT(Q_sparse);
  RTMB::REPORT(mu_at);
  RTMB::REPORT(ln_Y_at);

  return(jnLL)
}


#' Function to fit GMRF weight-at-age model from Cheng et al 2023
#' https://github.com/chengmatt/GMRF_WAA/blob/master/R_scripts/model%20runs/RTMB_GrowthModel.R
#'
#' @param data data.frame with columns "year", "age", and "weight". Weight should be in kg so that uncertainty calculations do not cause an issue.
#' @param weights likelihood weights for model.
#' @param n_proj_years the number of years to project forward
#' @param last_year
#' @param n.newton number of extra newton steps we want to take
#'
#' @return a list of GMRF models
#' @description
#' https://doi.org/10.1016/j.fishres.2023.106755
#'
#' @export
#'
FitGMRF_RTMB <- function(
    data = NULL,
    weights = NULL,
    n_proj_years = 2,
    last_year = NA,
    n.newton = 3
){

  # Reformat data for GMRF ----
  years <- do.call(seq, as.list(range(data$year)))
  proj_years <- (max(years) + 1):(max(years) + n_proj_years)
  ages <- do.call(seq, as.list(range(data$age)))

  # - Projections and fill missing years
  years_ages <- expand.grid(c(years, proj_years), ages)
  colnames(years_ages) <- c("year", "age")

  # Weight data
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
  Xsd_at[Xsd_at == 1] <- max(Xsd_at, na.rm = TRUE) # Set 0 sd to max

  # - Create an index for ages and years to feed into TMB, which helps construct the precision matrix
  ay_Index <- as.matrix(expand.grid("age" = seq_len(length(ages)),
                                    "year" = seq_len(length(years) + n_proj_years) ))


  # Set up TMB inputs  ----
  data_list <- list( years = years,
                     ages = ages,
                     X_at = X_at,
                     Xsd_at = Xsd_at,
                     ay_Index = ay_Index,
                     n_proj_years = n_proj_years,
                     Var_Param = 0) # Var_Param == 0 Conditional, == 1 Marginal

  parameters <- list( rho_a = 0,
                      rho_y = 0,
                      rho_c = 0,
                      log_sigma2 = log(0.1),
                      ln_L0 = log(min(data$weight, na.rm = TRUE)),
                      ln_Linf = log(max(data$weight, na.rm = TRUE)),  # Fixed at arbitrary value
                      ln_k = log(0.15),
                      ln_alpha = log(3.5e-7), # Start alpha at a reasonable space
                      # Starting value for alpha derived from a run where none of the rhos were estimated.
                      ln_beta = log(3), # Fix at isometric
                      ln_Y_at = array(0,dim=dim(X_at)) )


  # Run factorial models ----
  map_factorial <- tidyr::crossing(rho_y = 1, rho_c = 0:1, rho_a = 0:1) %>%
    # dplyr::filter(rowSums(.) > 1) %>%
    dplyr::slice(n()) %>%
    data.frame()

  # - Empty list to store model objects
  models <- list()

  for(n_fact in 1:nrow(map_factorial)) {

    # Create empty map list object
    map <- list()

    # Extract combinations of parameters estimated here
    map_fact <- map_factorial[n_fact,]

    # Create our mapping list to feed into MakeADFun
    for(m in 1:length(names(map_fact))) {

      if(map_fact[1,names(map_fact)[m]] == 0) { # if factorial = 0, turn estimation off
        map[[m]] <- factor(NA)
        names(map)[m] <- names(map_fact)[m] # Name list object
      } else{ # factorial == 1
        map[[m]] <- factor(1)
        names(map)[m] <- names(map_fact)[m] # Name list object
      } # ifelse statement for constructing map list object
    } # end m loop

    map <- c(map, list("ln_Linf" = factor(NA),
                       "ln_beta" = factor(NA)))

    # Build AD model function
    cmb <- function(f, d) function(p) f(p, d) ## Helper to make closure
    obj <- RTMB::MakeADFun(cmb(growth_3d, data_list), parameters = parameters, map = map, random = 'ln_Y_at', silent = TRUE)

    # Fit model
    kill_mod <- tryCatch({
      R.utils::withTimeout({
        fit <- stats::nlminb(obj$par, obj$fn, obj$gr,
                             control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))
        return(kill_mod = FALSE)
      },
      timeout = 60*5)
    },
    error = function(e){
      return(kill_mod = TRUE)
    },
    TimeoutException = function(e){
      return(kill_mod = TRUE)
    })

    # Save if optimized
    if(!kill_mod){

      # Take some additional newton steps to make sure we reach a minimum
      tryCatch(expr = for(i in 1:n.newton) {
        g = as.numeric(obj$gr(fit$par))
        h = optimHess(fit$par, fn = obj$fn, gr = obj$gr)
        fit$par = fit$par - solve(h,g)
        fit$objective = obj$fn(fit$par)
      }, error = function(e){e})

      # Save optimized model results
      # Get report
      report <- obj$report(obj$env$last.par.best)
      report$Y_at <- exp(report$ln_Y_at)
      colnames(report$Y_at) <- c(years, proj_years)
      rownames(report$Y_at) <- ages

      # Predicted weight-at-age matrix ----
      predicted_weight <- report$Y_at %>%
        as.data.frame() %>%
        mutate(age = ages) %>%
        tidyr::pivot_longer(!age) %>%
        dplyr::rename(pred_weight = value, year = name) %>%
        dplyr::mutate(
          year = as.numeric(year),
          model = paste0("GMRF", n_fact),
          projection = year %in% proj_years,
          last_year = last_year) %>%
        dplyr::select(model, year, age, pred_weight, last_year, projection) %>%
        dplyr::arrange(year, age) %>%
        as.data.frame()

      # Predicted observations ----
      data_tmp <- predicted_weight %>%
        merge(., data, by = c('year','age'), all = TRUE) %>%
        as.data.frame()

      # Save
      models[[n_fact]] <- list(obj = obj, map = map_factorial[n_fact,], opt = fit, report = report,
                               data = data_tmp,
                               predicted_weight = predicted_weight)
    }

    # print(n_fact)
  } # loop through to run multiple models

  return(models)
}
