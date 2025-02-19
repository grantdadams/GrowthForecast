#' Forecast growth
#'
#' Fits various weight-at-age models given input data
#' @param form model formula for machine learning based methods and selection of covariates
#' @param data data.frame with columns "year", "age", and "weight" and any covariates. Weight should be in kg so that uncertainty calculations do not cause an issue.
#' @param weights likelihood weights for model.
#' @param n_proj_years the max number of years to project forward (default = 2)
#'
#' @return
#' @export
#'
#' @examples
ForestGrowth <- function(form = formula(weight~age+year), data = NULL, n_proj_years = 2, peels = 4){

  # Data checks ----
  if(sum(c("weight", "age", "year") %in% colnames(data)) != 3){
    stop(print("Columns do not include `weight`, `age`, or `year`"))
  }

  if(peels >= length(unique(data$year))){
    stop(print("Number of retrospective peels is greater than the number of years"))
  }

  if(sum(data$age %% 1 != 0)){
    stop(print("'age' contains non-integer ages"))
  }


  # Run peels
  peel_list <- list()
  projection_list <- list()
  test_list <- list()

  for(i in 1:peels){

    # Peel data ----
    last_year <- max(data$year, na.rm = TRUE) - (i-1)
    # last_year <- max(data$year, na.rm = TRUE) - (n_proj_years+(i-1))

    train <- data %>%
      dplyr::filter(year <= last_year)

    test <- data %>%
      dplyr::filter(year <= (last_year + n_proj_years) & year > last_year) # Only projection years


    # Fit models ----
    # * VBGF ----
    vbgf <- FitVBGF(
      data = train,
      weights=NULL,
      n_proj_years = n_proj_years,
      last_year = last_year
    )

    # * WtAgeRe ----
    wtagere <- FitWtAgeRE(
      data = train,
      weights=NULL,
      # - Number of projection years
      n_proj_years = n_proj_years
    )

    # * GMRF ----
    gmrf <- FitGMRF(
      data = train,
      weights=NULL,
      n_proj_years = n_proj_years,
      last_year = last_year
    )

    # * NN ----
    #FIXME: no likelihood weights
    nn_init <- NULL
    if(i != 1){nn_init = nn$obj$weights}
    nn <- fit_nn(
      data = train,
      startweights = nn_init,
      n_proj_years = n_proj_years,
      last_year = last_year)

    # * LSTM ----
    lstm <- fit_lstm_rtmb(data = train,
                          nhidden_layer = 2,
                          hidden_dim = 5,
                          input_par = NULL,
                          last_year = last_year)


    # Combine ----
    peel_list[[i]] <- list(vbgf = vbgf,
                           wtagere = wtagere,
                           gmrf1 = gmrf[[1]],
                           gmrf2 = gmrf[[2]],
                           gmrf3 = gmrf[[3]],
                           gmrf4 = gmrf[[4]],
                           nn = nn,
                           lstm = lstm)

    ## retrospective forecast performance (train data, from models)
    projection_list[[i]]   <-   do.call("rbind",
                                        lapply(peel_list[[i]],
                                               function(x){x$prediction %>%
                                                   dplyr::select(year, age, obs_weight,
                                                                 pred_weight, model, projection)}
                                        )) %>%
      # tidyr::pivot_wider(names_from = c(model), values_from = pred_weight) %>%
      mutate(terminal_train_year = last_year,
             peel_id = i)

    ## test performance (terminal train year + 2)

    test_list[[i]]   <-   test %>%
      select(year, age, obs_weight= weight   )   %>%
      merge(.,  projection_list[[i]]  %>% filter(projection) %>% select(-obs_weight) %>%
              mutate(year = as.numeric(year)),
            by = c('year','age'))


  }

  names(peel_list) <- 1:peels
  names(projection_list) <- 1:peels
  names(test_list) <- 1:peels


  # Performance metrics ----

  # Calculate overall RSE for each model and peel
  rse_table   <- do.call("rbind", lapply(1:length(test_list), function(i) {
    test_list[[i]]  %>%
      dplyr::mutate(YID = paste0('y+',year - terminal_train_year)) %>%
      dplyr::summarise(RSE = sqrt(sum((obs_weight - pred_weight)^2) / n()), .by = c(model,YID,peel_id))
  })) %>%
    summarise(mean_RSE = mean(RSE),.by = c(model, YID)) %>% ## average across peels
    arrange(mean_RSE)

  # Calculate RSE by model and age for each peel
  rse_table_by_age <- do.call("rbind", lapply(1:length(test_list), function(i) {
    test_list[[i]]  %>%
      dplyr::mutate(YID = paste0('y+',year - terminal_train_year)) %>%
      dplyr::summarise(RSE = sqrt(sum((obs_weight - pred_weight)^2) / n()), .by = c(model,age,YID,peel_id))
  })) %>%
    summarise(mean_RSE = mean(RSE),.by = c(model,age, YID)) %>% ## average across peels
    arrange(age,mean_RSE)


  # Pick best model based on overall RSE
  best_mods <- rse_table %>%
    group_by(YID) %>%
    filter(mean_RSE == min(mean_RSE)) %>%
    ungroup() %>%
    select(YID, model, mean_RSE)

  # - Output (lookup projection in first peel)
  test_list[[1]] %>%
    ## only good performers
    filter(model %in% best_mods$model) %>%
    ## get all future (test) years
    filter(year %in% max(data$year, na.rm = TRUE):(max(data$year, na.rm = TRUE)+n_proj_years)) %>%
    filter(!duplicated(.$pred_weight)) %>%
    select(year, age, pred_weight) %>%
    tidyr::pivot_wider(., names_from = age, values_from = pred_weight)



}


