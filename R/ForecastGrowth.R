

FitGMRF <- function(
    data = NULL,
    weights=NULL){


  # Reformat data for GMRF ----
  colnames(data) <- tolower(colnames(data))
  years <- do.call(seq, as.list(range(data$year)))
  ages <- do.call(seq, as.list(range(data$age)))
  years_ages <- expand.grid(years, ages)
  colnames(years_ages) <- c("year", "age")

  gmrf_data <- data %>%
    dplyr::group_by(year, age) %>%
    dplyr::summarise(mn_weight = ifelse(!is.null(weights), weighted.mean(weight, weights, na.rm = TRUE),
                              mean(weight, na.rm = TRUE)),
                     sd = ifelse(!is.null(weights), sqrt(sum(weights * (weight - mn_weight)^2, na.rm = TRUE)),
                                     sd(weight, na.rm = TRUE)),
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
  # - Number of projection years
  n_proj_years <- 3

  # - Years in data
  years <- gmrf_mn$year

  # - Ages
  ages <- as.numeric(colnames(gmrf_mn)[-1])

  # - Read in data weight at age matrix
  X_at <- t(as.matrix(gmrf_mn[,-1])) # removing first col (year column)

  # - Create projection columns
  proj_cols <- matrix(NA, nrow = length(ages), ncol = n_proj_years)

  # - Append NA for projection year
  X_at <- cbind(X_at, proj_cols)

  # - Read in standard deviations for weight at age matrix
  Xse_at <- t(as.matrix(gmrf_sd[,c(-1)])) # removing first col (year column)

  # - Convert to CV and sd in lognormal space
  Xcv_at <- sqrt( (exp(Xse_at^2) - 1) )
  Xsd_at <- sqrt((log((Xcv_at)^2 + 1))/(log(10)^2))

  # - Create an index for ages and years to feed into TMB, which helps construct the precision matrix
  ay_Index <- as.matrix(expand.grid("age" = seq_len(length(ages)),
                                    "year" = seq_len(length(years) + n_proj_years) ))


  # Set up TMB inputs  ----
  data <- list( years = years,
                ages = ages,
                X_at = X_at,
                Xsd_at = Xsd_at,
                ay_Index = ay_Index,
                n_proj_years = n_proj_years,
                Var_Param = 1) # Var_Param == 0 Conditional, == 1 Marginal

  parameters <- list( rho_y = 0,
                      rho_a = 0,
                      rho_c = 0,
                      log_sigma2 = log(0.1),
                      ln_L0 = log(45),
                      ln_Linf = log(80),  # Fixed at arbitrary value
                      ln_k = log(0.15),
                      ln_alpha = log(3.5e-7), # Start alpha at a reasonable space
                      # Starting value for alpha derived from a run where none of the rhos were estimated.
                      ln_beta = log(3), # Fix at isometric
                      ln_Y_at = array(0,dim=dim(X_at)) )


  # Run factorial models ----
  map_factorial <- tidyr::crossing(rho_y = 0:1, rho_c = 0:1, rho_a = 0:1) %>%
    data.frame()

  # - Define number of extra newton steps we want to take
  n.newton <- 3

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

    map <- c(map, list(map, "ln_Linf" = factor(NA),
                              "ln_beta" = factor(NA)))

    # Now, make AD model function
    waa_model <- TMB::MakeADFun(data = data, parameters = parameters,
                           map = map, random = "ln_Y_at",
                           DLL = "GMRF_WAA", silent = TRUE)

    # Now, optimize the function
    waa_optim <- stats::nlminb(waa_model$par, waa_model$fn, waa_model$gr,
                               control = list(iter.max = 1e5, eval.max = 1e5))

    # Take some additional newton steps to make sure we reach a minimum
    tryCatch(expr = for(i in 1:n.newton) {
      g = as.numeric(waa_model$gr(waa_optim$par))
      h = optimHess(waa_optim$par, fn = waa_model$fn, gr = waa_model$gr)
      waa_optim$par = waa_optim$par - solve(h,g)
      waa_optim$objective = waa_model$fn(waa_optim$par)
    }, error = function(e){e})

    # Save optimized model results
    waa_model$optim <- waa_optim

    # Get report
    waa_model$rep <- waa_model$report(waa_model$env$last.par.best)

    # Get sd report
    waa_model$sd_rep <- TMB::sdreport(waa_model)

    models[[n_fact]] <- waa_model

    print(n_fact)

  } # loop through to run multiple models

}
