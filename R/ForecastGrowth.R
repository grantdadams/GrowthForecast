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
ForecastGrowth <- function(form = formula(weight~age+year), data = NULL, n_proj_years = 2, peels = 4){

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
                          n_proj_years = n_proj_years,
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
      dplyr::mutate(terminal_train_year = last_year,
             peel_id = i)

    ## fill in test data (what the model has not seen, for RSE calc)
    if(i != 1){
      test_list[[i]]   <-   test %>%
        dplyr::select(year, age, obs_weight= weight   )   %>%
        merge(.,  projection_list[[i]]  %>%
                filter(projection) %>%
                dplyr::select(-obs_weight) %>%
                mutate(year = as.numeric(year)),
              by = c('year','age'))

    } else{
      ## first peel is pure forecast, obs_weight is NA
      test_list[[i]]   <-     projection_list[[i]]  %>%
        dplyr::filter(projection) %>%
        dplyr::mutate(year = as.numeric(year))

    } ## end test_list
  } ## end peels

  names(peel_list) <- 1:peels
  names(projection_list) <- 1:peels
  names(test_list) <- 1:peels


  # Performance metrics ----

  # Calculate overall RSE for each model and peel
  ## create master dataframes
  ## TODO consider if/where to include 5-year mean calc, should it be moving with peel?
  peel_list_summary <-do.call("rbind", lapply(1:length(peel_list), function(i) { peel_list[[i]]}))
  test_list_summary <-do.call("rbind", lapply(1:length(test_list), function(i) { test_list[[i]]}))
  projection_list_summary <-do.call("rbind", lapply(1:length(projection_list), function(i) { projection_list[[i]]}))

  ## drop first peel as it is raw projection (no observations)
  rse_table   <-  test_list_summary %>%
    filter(peel_id > 1) %>%
    dplyr::mutate(YID = paste0('y+',year - terminal_train_year)) %>%
    dplyr::summarise(RSE = sqrt(sum((obs_weight - pred_weight)^2) / n()), .by = c(model,YID,peel_id)) %>%
    summarise(mean_RSE = mean(RSE),.by = c(model, YID)) %>% ## average across peels
    arrange(mean_RSE)

  # Calculate RSE by model and age for each peel
  rse_table_by_age <- test_list_summary %>%
    dplyr::filter(peel_id > 1) %>%
    dplyr::mutate(YID = paste0('y+',year - terminal_train_year)) %>%
    dplyr::summarise(RSE = sqrt(sum((obs_weight - pred_weight)^2) / n()), .by = c(model,age,YID,peel_id)) %>%
    dplyr::summarise(mean_RSE = mean(RSE),.by = c(model,age, YID)) %>% ## average across peels
    arrange(age,mean_RSE)


  # Pick best model based on overall RSE
  best_mods <- rse_table %>%
    dplyr::group_by(YID) %>%
    dplyr::filter(mean_RSE == min(mean_RSE)) %>%
    ungroup() %>%
    dplyr::select(YID, model, mean_RSE)

  # - Output (lookup projection in first peel)
  projected_waa <- test_list[[1]] %>%
    dplyr::mutate(YID = paste0('y+',year - terminal_train_year)) %>%
    ## pull best model-forecast combination
    inner_join(best_mods, by = c("YID", "model")) %>%
    ## get all future (test) years
    filter(year %in% max(data$year, na.rm = TRUE):(max(data$year, na.rm = TRUE)+n_proj_years)) %>%
    # filter(!duplicated(.$pred_weight)) %>%
    dplyr::select(year, age, pred_weight) %>%
    tidyr::pivot_wider(., names_from = age, values_from = pred_weight) %>%
    arrange(year) %>%
    mutate(model = paste0("#",best_mods$model))

  # summary figures ----
  projection_list_summary$year <- as.numeric(projection_list_summary$year)
  plotdf0 <- bind_rows(projection_list_summary,test_list_summary) %>%
    distinct()

  plot_list <- list()

  ## plot historical data and fits
  for(i in unique(plotdf0$model)){
    plotdf <- plotdf0 %>%
      ## TODO customize the year bounds or change how this is displayed
      filter(model == i & year > (max(data$year)-(n_proj_years*peels)) ) ## truncate historical years
    plot_list[[i]] <- ggplot(data= NULL, aes(x = age)) +
      geom_point(data = plotdf,
                 alpha = 0.05,
                 size = 2,aes(y = obs_weight, color = factor(projection))) +
      geom_line(data = plotdf,
                aes(y = pred_weight,
                    color = factor(projection),
                    group = interaction(year, peel_id)))+
      scale_color_manual(values = c('grey22','blue'))+
      facet_grid(peel_id ~ year, scales = 'free_y') +
      theme_minimal() +
      theme(legend.position = 'none') +
      labs(x = 'age', y = 'weight', title = i)
  }

  ## plot projected waa
  plot_list[[length(plot_list)+1]] <- reshape2::melt(projected_waa, id = c('year','model')) %>%
    ungroup() %>%
    mutate(age = as.numeric(variable)) %>%
    ggplot(., aes(x = age, y = value, color = interaction(year, model))) +
    scale_color_manual(values = grey.colors(n = n_proj_years, start = 0.3, end = 0.5)) +
    geom_line(lwd = 1)+
    theme_minimal() +
    theme(legend.position = 'bottom')+
    labs(x = 'age', y = 'weight', title = 'projected weight at age',
         color = '')


  return(list(peel_list_summary = peel_list_summary,
              projection_list_summary=projection_list_summary,
              test_list_summary= test_list_summary,
              rse_table = rse_table,
              rse_table_by_age =rse_table_by_age,
              best_mods = best_mods,
              projected_waa = projected_waa,
              plot_list = plot_list))



} ## end function


