#' Forecast growth
#'
#' Fits various weight-at-age models given input data and evaluates forecast skill using retrospective forecasting. Outputs the forecasted weight-at-age from the "best" performing model.
#'
#' @param form model formula for machine learning based methods and selection of covariates
#' @param data data.frame with columns "year", "age", and "weight" and any covariates. Weight should be in kg so that uncertainty calculations do not cause an issue.
#' @param n_proj_years the max number of years to project forward (default = 2)
#' @param peels the number of years to peel back for retrospective forecast skill evaluation (default = 3)
#' @param maturity_vec a vector of length 1:max(age) indicating proportion mature
#' @description
#' The function takes a \code{data.frame} of weight-at-age data and compares the retrospective forecast skill of a variety of weight-at-age models. This is done by peeling the data \code{peels} times and comparing the observed vs forecasted data using root mean squared error (RMSE) across peels and all ages (or by each age) for each of the projection years set by \code{n_proj_years}. The function then outputs the forecasted weight-at-age based on the model with the lowest RMSE for each projection year.
#'
#'
#' @return
#' A list with:
#' * all_predictions = data.frame of observed and model predicted (hindcasts & forecasts) weight-at-age for each model and retrospective peel
#' * all_forecasts = data.frame of forecasted weight-at-age for each model and retrospective peel
#' * rmse_table = data.frame of root mean squared error across all ages/peels from retrospective forecasts
#' * rmse_table_by_age = data.frame of root mean squared error by age across peels from retrospective forecasts
#' * best_mods = data.frame with of the best model from\code{rmse_table}
#' * best_mods_by_age = data.frame with of the best model from\code{rmse_table_by_age} for each age
#' * best_forecast = data.frame with forecasted weight-at-age from\code{best_mods}
#' * plots = list of plots from each model
#'
#' @export
#'
#'
#' @examples
ForecastGrowth <- function(form = formula(weight~age+year), data = NULL, n_proj_years = 2, peels = 3,
                           maturity_vec = NULL){

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

  # Data characteristics ----
  years <- sort(unique(data$year))
  ages <- seq(range(data$age)[1], range(data$age)[2], by = 1 )

  # Run peels ----
  peels = peels + 1
  non_converged = c(NA)
  peel_list <- list() # List of models for each peel
  pred_list <- list() # List of hindcast/forecast predictions for each peel
  test_list <- list() # List of forecast predictions for each peel

  for(i in 1:peels){

    # Peel data ----
    last_year <- rev(years)[i]

    train <- data %>%
      dplyr::filter(year <= last_year)

    test <- data %>%
      dplyr::filter(year <= (last_year + n_proj_years) & year > last_year) # Only projection years


    # Fit models ----
    # * VBGF ----
    vbgf <-  tryCatch(
      suppressMessages(
        FitVBGF(
          data = train,
          weights=NULL,
          n_proj_years = n_proj_years,
          last_year = last_year
        )
      ),
      error = function(e) return(NULL)
    )

    # * WtAgeRe ----
    wtagere <-
      tryCatch(
        suppressMessages(
          FitWtAgeRE(
            data = train,
            weights=NULL,
            # - Number of projection years
            n_proj_years = n_proj_years,
            last_year = last_year
          )
        ),
        error = function(e) return(NULL)
      )

    # * GMRF ----
    gmrf <-  tryCatch(
      suppressMessages(
        FitGMRF(
          data = train,
          weights=NULL,
          n_proj_years = n_proj_years,
          last_year = last_year
        )
      ),
      error = function(e) return(NULL)
    )

    # * NN ----
    #FIXME: no likelihood weights
    #FIXME: determine whether this model is worth keeping
    # nn_init <- NULL
    # if(i != 1){nn_init = nn$obj$weights}
    # nn <-  tryCatch(
    #   suppressMessages(
    #     fit_nn(
    #       data = train,
    #       startweights = nn_init,
    #       n_proj_years = n_proj_years,
    #       last_year = last_year)
    #   ),
    #   error = function(e) return(NULL)
    # )

    # * LSTM ----
    lstm <- if(!"lstm" %in% non_converged){
      tryCatch(
        suppressMessages(
          fit_lstm_rtmb(data = train,
                        nhidden_layer = 2,
                        hidden_dim = 3,
                        n_proj_years = n_proj_years,
                        input_par = NULL,
                        last_year = last_year)
        ),
        error = function(e) return(NA)
      )
    }else{
      NA
    }

    # * Average of previous 5-years ----
    avg5 <- list(mean = train %>%
                   dplyr::filter(year %in% ((max(train$year)-5):max(train$year))) %>%
                   dplyr::summarise(pred_weight = mean(weight),.by = age) %>%
                   dplyr::arrange(age))

    avg5[['hindcast']] <-  avg5$mean %>%
      merge(., train) %>%
      dplyr::mutate(last_year = last_year,
                    projection = FALSE,
                    model = 'avg5') %>%
      dplyr::select(model, last_year, year, age, obs_weight = weight, pred_weight, projection) %>%
      arrange(year,age)

    avg5[['forecast']] <-  avg5$mean %>%
      merge(., data.frame(age = rep(ages, 2),
                          year = rep((last_year+1):(last_year + n_proj_years), each = length(ages)))) %>%
      dplyr::mutate(last_year = last_year,
                    weight = NA,
                    projection = TRUE,
                    model = 'avg5') %>%
      dplyr::select(model, last_year, year, age, obs_weight = weight, pred_weight, projection) %>%
      arrange(year,age)

    avg5[['prediction']] = rbind(avg5[['hindcast']],
                                 avg5[['forecast']])

    # Combine ----
    peel_list[[i]] <- list(vbgf = vbgf,
                           wtagere = wtagere,
                           gmrf1 = gmrf[[1]],
                           gmrf2 = gmrf[[2]],
                           gmrf3 = gmrf[[3]],
                           gmrf4 = gmrf[[4]],
                           # nn = nn,
                           lstm = lstm,
                           avg5 = avg5)

    # * Filter non-converged models ----
    non_converged <- c(non_converged,
                       names(peel_list[[i]])[sapply(peel_list[[i]], function(x) length(x) == 1)] # remove
    )
    peel_list[[i]] <- peel_list[[i]][sapply(peel_list[[i]], function(x) length(x) != 1)] # keep

    # * Pull predictions (hindcast and forecast) ----
    pred_list[[i]]   <-   do.call("rbind",
                                  lapply(peel_list[[i]],
                                         function(x){x$prediction %>%
                                             dplyr::select(year, age, obs_weight,
                                                           pred_weight, model, projection)}
                                  )) %>%
      dplyr::mutate(terminal_train_year = last_year,
                    peel_id = i)

    # * Fill in test data ----
    # - what the model has not seen, for rmse calc
    if(i != 1){
      test_list[[i]]   <-   test %>%
        dplyr::select(year, age, obs_weight= weight   )   %>%
        merge(.,  pred_list[[i]]  %>%
                filter(projection) %>%
                dplyr::select(-obs_weight) %>%
                mutate(year = as.numeric(year)),
              by = c('year','age'),
              all = TRUE)

    } else{
      ## first peel is pure forecast, obs_weight is NA
      test_list[[i]]   <-     pred_list[[i]]  %>%
        dplyr::filter(projection) %>%
        dplyr::mutate(year = as.numeric(year))

    } ## end test_list
  } ## end peels

  names(peel_list) <- 1:peels
  names(pred_list) <- 1:peels
  names(test_list) <- 1:peels
  non_converged <- unique(non_converged)

  # Performance metrics ----
  # * Create master dataframes ----
  test_list_summary <- do.call("rbind", lapply(1:length(test_list), function(i) { test_list[[i]]})) %>%
    dplyr::filter(!model %in% non_converged)
  pred_list_summary <- do.call("rbind", lapply(1:length(pred_list), function(i) { pred_list[[i]]})) %>%
    dplyr::filter(!model %in% non_converged)

  # * Calculate overall rmse for each model and peel ----
  ## drop first peel as it is raw projection (no observations)
  maa <- cbind('age' = 1:length(maturity_vec), 'maturity' = maturity_vec)

  rmse_table   <-  test_list_summary %>%
    filter(peel_id > 1) %>%
    merge(.,maa, by = 'age') %>%
    dplyr::mutate(YID = paste0('y+',year - terminal_train_year)) %>%
    dplyr::summarise(RMSE = sqrt(mean((obs_weight - pred_weight)^2, na.rm = TRUE)),
                     RMSE_mat =  sqrt(mean(maturity*(obs_weight - pred_weight)^2, na.rm = TRUE)),
                     .by = c(model, YID)) %>%
    arrange(YID, RMSE_mat) %>%
    dplyr::select(YID, model, RMSE, RMSE_mat)

  # * Calculate RSE by model and age for each peel ----
  rmse_table_by_age <- test_list_summary %>%
    dplyr::filter(peel_id > 1) %>%
    dplyr::mutate(YID = paste0('y+', year - terminal_train_year)) %>%
    dplyr::group_by(age, YID) %>%
    dplyr::slice_sample(n = 50 ) %>%
    ungroup() %>%
    dplyr::summarise(RMSE = sqrt(mean((obs_weight - pred_weight)^2, na.rm = TRUE)), .by = c(model, YID, age)) %>%
    dplyr::mutate(var_RMSE = var(RMSE, na.rm = TRUE), .by = c(YID, age)) %>%
    filter(!is.na(RMSE)) %>%
    arrange(YID, age, RMSE) %>%
    dplyr::select(YID, age, model, RMSE, var_RMSE)

  # Pick best model ----
  # * Based on overall RSE ----
  best_mods <- rmse_table %>%
    dplyr::group_by(YID) %>%
    dplyr::filter(RMSE == min(RMSE, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::select(YID, model, RMSE)

  # - Expand if biennial or triennial survey
  best_mods <- best_mods %>%
    dplyr::full_join(rmse_table %>%
                       dplyr::distinct(YID)) %>%
    dplyr::arrange(YID) %>%
    dplyr::mutate(model = ifelse(is.na(model), lead(model), model)) %>%
    dplyr::mutate(model = ifelse(is.na(model), lead(model), model))

  # * Based on age-specific RSE ----
  best_mods_by_age <- rmse_table_by_age %>%
    dplyr::group_by(YID, age) %>%
    dplyr::filter(RMSE == min(RMSE, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::select(YID, age, model, RMSE)

  # - Expand if biennial or triennial survey
  best_mods_by_age <- best_mods_by_age %>%
    dplyr::full_join(rmse_table_by_age %>%
                       dplyr::distinct(YID, age)) %>%
    dplyr::arrange(age, YID) %>%
    dplyr::mutate(model = ifelse(is.na(model), lead(model), model)) %>%
    dplyr::mutate(model = ifelse(is.na(model), lead(model), model)) %>%
    dplyr::select(YID, age, model, RMSE)


  # * Based on maturity-weighted RSE ----
  best_mods_by_mat <- rmse_table %>%
    dplyr::group_by(YID) %>%
    dplyr::filter(RMSE_mat == min(RMSE_mat, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::select(YID, model, RMSE_mat)

  # - Expand if biennial or triennial survey
  best_mods_by_mat <- best_mods_by_mat %>%
    dplyr::full_join(rmse_table %>%
                       dplyr::distinct(YID)) %>%
    dplyr::arrange(YID) %>%
    dplyr::mutate(model = ifelse(is.na(model), lead(model), model)) %>%
    dplyr::mutate(model = ifelse(is.na(model), lead(model), model))


  # - Output (lookup projection in first peel)
  projected_waa <- test_list[[1]] %>%
    dplyr::mutate(YID = paste0('y+',year - terminal_train_year)) %>%
    ## pull best model-forecast combination
    inner_join(best_mods, by = c("YID", "model")) %>%
    ## get all future (test) years
    filter(year %in% max(data$year, na.rm = TRUE):(max(data$year, na.rm = TRUE)+n_proj_years)) %>%
    dplyr::select(year, age, pred_weight) %>%
    tidyr::pivot_wider(., names_from = age, values_from = pred_weight) %>%
    arrange(year) %>%
    mutate(model = paste0("#",best_mods$model))

  # Summary figures ----
  pred_list_summary$year <- as.numeric(pred_list_summary$year)
  plotdf0 <- bind_rows(pred_list_summary,test_list_summary) %>%
    distinct()

  plot_list <- list()

  # * Plot historical data and fits ----
  for(i in unique(plotdf0$model)){
    plotdf <- plotdf0 %>%
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

  # * single plot with best projection
  plotdf <- plotdf0 %>%
    filter(model  %in% best_mods$model &
             year > (max(data$year)-(n_proj_years*peels)) ) ## truncate historical years

  if(length(unique(best_mods$model))>1) colVec <- c('blue','dodgerblue','#9A031E','#E36923')
  if(length(unique(best_mods$model))==1) colVec <- c('grey22','blue')

  plot_list[[length(plot_list)+1]] <- ggplot(data= NULL, aes(x = age)) +
    geom_point(data = plotdf,
               alpha = 0.05,
               color = 'grey22',
               size = 2,aes(y = obs_weight)) +
    geom_line(data = plotdf,
              lwd = 1,
              aes(y = pred_weight,
                  color = interaction(projection, model),
                  group = interaction(model, year, peel_id)))+
    scale_color_manual(values = colVec)+
    facet_grid(peel_id ~ year, scales = 'free_y') +
    theme_minimal() +
    theme(legend.position = 'bottom') +
    labs(x = 'age', y = 'weight', title = 'best model(s)', color = '')


  # * Plot projected waa ----
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

  names(plot_list) <- c(names(peel_list[[peels]]), "projection")


  # Return ----
  return(list(all_predictions = pred_list_summary,
              all_forecasts = test_list_summary,
              rmse_table = rmse_table,
              rmse_table_by_age = rmse_table_by_age,
              best_mods = best_mods,
              best_mods_by_age = best_mods_by_age,
              best_mods_by_mat = best_mods_by_mat,
              best_forecast = projected_waa,
              plots = plot_list))

} ## end function


