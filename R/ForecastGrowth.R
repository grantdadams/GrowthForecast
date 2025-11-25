#' Compare NN and process-based growth models
#'
#' Fits various weight-at-age models given input data and evaluates hindcast and forecast skill using retrospective forecasting.
#'
#' @param data data.frame with columns "year", "age", and "weight" and any covariates. Weight should be in kg so that uncertainty calculations do not cause an issue.
#' @param n_proj_years the max number of years to project forward (default = 2)
#' @param peels the number of years to peel back for retrospective forecast skill evaluation (default = 3) or number of k-folds to do for cross-validation
#' @param forecast evaluate forecast performance (TRUE) or hindcast performance (FALSE)
#' @param stop_if_non_converged stops if any of the models do not converge (default = FALSE)
#' @description
#' The function takes a \code{data.frame} of weight-at-age data and compares the hindcast and forecast forecast skill of a variety of neural net and parametric weight-at-age models. This is done by peeling the data \code{peels} times and comparing the observed vs forecasted data using root mean squared error (RMSE) across peels and all ages (or by each age) for each of the projection years set by \code{n_proj_years}.
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
#' * plots = list of plots from each model
#'
#' @export
#'
#'
#' @examples
ForecastGrowth <- function(data = NULL,
                           n_proj_years = 2,
                           peels = 3,
                           forecast = TRUE,
                           stop_if_non_converged = FALSE,
                           seed = 123){
  set.seed(seed)

  # Data checks ----
  if(any(!c("weight", "age", "year") %in% colnames(data))){
    stop(print("Columns do not include `weight`, `age`, or `year`"))
  }

  if(any(sapply(data, class) == "character")){
    warning(print("Converting character columns to numeric"))
    data <- data %>%
      dplyr::mutate_if(is.character, as.numeric)
  }

  if(peels >= length(unique(data$year))){
    stop(print("Number of retrospective peels is greater than the number of years"))
  }

  if(sum(data$age %% 1 != 0)){
    stop(print("'age' contains non-integer ages"))
  }

  years <- sort(unique(data$year))
  kfold_id = sample(1:peels, length(years), replace = TRUE)
  ages <- seq(range(data$age)[1], range(data$age)[2], by = 1)


  # Compare prediction skill ----
  if(forecast){peels = peels + 1}
  non_converged = c()
  peel_list <- list()        # List of models for each peel
  pred_weight_list <- list() # List of predicted weight-at-age for each age/yeal from each model and peel
  train_list <- list()       # List of data with predicted weight from each peel (replicated for each model)
  test_list <- list()        # List of forecast predictions for each peel

  # - Retrospective peels
  for(i in 1:peels){

    # * Subset data ----
    last_year <- rev(years)[i]
    # - if evaluating forecast skill
    if(forecast){
      train <- data %>%
        dplyr::filter(year <= last_year)

      test <- data %>%
        dplyr::filter(year <= (last_year + n_proj_years) & year > last_year) # Only projection years
    } else {
      # - K-fold cross validation
      train <- data %>%
        dplyr::filter(!year %in% years[which(kfold_id == i)])

      test <- data %>%
        dplyr::filter(year %in% years[which(kfold_id == i)])
    }


    # * Fit models ----
    # ** VBGF ----
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

    # ** VBGF-RE ----
    vbgf_re <-
      tryCatch(
        suppressMessages(
          FitVBGF_RE(
            data = train,
            weights=NULL,
            n_proj_years = n_proj_years,
            last_year = last_year
          )
        ),
        error = function(e) return(NULL)
      )

    # ** VBGF-RE (Ianelli) ----
    vbgf_ianelli <-
      tryCatch(
        suppressMessages(
          FitWtAgeRE(
            data = train,
            weights=NULL,
            n_proj_years = n_proj_years,
            last_year = last_year
          )
        ),
        error = function(e) return(NULL)
      )


    # ** GMRF ----
    gmrf <-  tryCatch(
      suppressWarnings(
        suppressMessages(
          FitGMRF_RTMB(
            data = train,
            weights=NULL,
            n_proj_years = n_proj_years,
            last_year = last_year
          )
        )
      ),
      error = function(e) return(NULL)
    )

    # ** NN ----
    nn <-
      tryCatch(
        suppressMessages(
          fit_nn_keras(
            form = formula(~age + year),
            data = train,
            n_proj_years = n_proj_years,
            last_year = last_year)
        ),
        error = function(e) return(NULL)
      )

    # ** LSTM ----
    lstm <-
      tryCatch(
        suppressWarnings(
          suppressMessages(
            fit_lstm_rtmb(data = train,
                          nhidden_layer = 2,
                          hidden_dim = 3,
                          n_proj_years = n_proj_years,
                          input_par = NULL,
                          last_year = last_year)
          )
        ),
        error = function(e) return(NULL)
      )


    # ** XGboost
    XGBoost <-  tryCatch(
      suppressWarnings(
        suppressMessages(
          FitXGBoost(
            data = train,
            weights=NULL,
            n_proj_years = n_proj_years,
            last_year = last_year
          )
        )
      ),
      error = function(e) return(NULL)
    )



    # ** Average of previous 5-years ----
    avg5 <- average_weight(
      data = train,
      nyrs=5,
      n_proj_years = n_proj_years,
      last_year = last_year)


    # Combine ----
    peel_list[[i]] <- list(vbgf = vbgf,
                           vbgf_re = vbgf_re,
                           vbgf_ianelli = vbgf_ianelli,
                           gmrf1 = gmrf[[1]],
                           # gmrf2 = gmrf[[2]],
                           # gmrf3 = gmrf[[3]],
                           # gmrf4 = gmrf[[4]],
                           nn = nn,
                           lstm = lstm,
                           XGBoost = XGBoost,
                           avg5 = avg5)


    # Stop procedure if any models dont converge
    if(stop_if_non_converged){
      if(length(non_converged) > 0){
        error("Some models did not converge")
      }
    }


    # * Filter non-converged models ----
    non_converged <- c(non_converged,
                       names(peel_list[[i]])[sapply(peel_list[[i]], function(x) length(x) == 1 | is.null(x))] # remove
    )
    peel_list[[i]] <- peel_list[[i]][sapply(peel_list[[i]], function(x) length(x) != 1 & !is.null(x))] # keep


    # * Pull predictions to training data ----
    train_list[[i]]   <-   do.call("rbind",
                                   lapply(peel_list[[i]],
                                          function(x){x$data %>%
                                              dplyr::select(year, age, weight,
                                                            pred_weight, model, projection)}
                                   )) %>%
      dplyr::mutate(terminal_train_year = last_year,
                    peel_id = i) %>%
      dplyr::filter(!is.na(weight))


    # * Pull predictions for each age/year ----
    pred_weight_list[[i]] <- do.call("rbind",
                                     lapply(peel_list[[i]], function(x){x$predicted_weight})
    ) %>%
      dplyr::mutate(terminal_train_year = last_year,
                    peel_id = i) %>%
      dplyr::mutate(year = as.numeric(year))


    # * Fill in test data ----
    # - what the model has not seen, for RMSE calc
    if(i == 1 & forecast){
      ## first peel is pure forecast, weight is NA
      test_list[[i]] <- pred_weight_list[[i]]  %>%
        dplyr::filter(projection) %>%
        dplyr::mutate(weight = NA) %>%
        dplyr::select(model, year, age, weight, pred_weight, terminal_train_year, projection, peel_id)
    } else{
      if(forecast){
        pred_tmp <- pred_weight_list[[i]]  %>%
          filter(projection)
      } else{
        pred_tmp <- pred_weight_list[[i]]
      }
      test_list[[i]] <- test %>%
        dplyr::select(year, age, weight) %>%
        dplyr::full_join(pred_tmp,
                         by = join_by(year, age), relationship = "many-to-many") %>%
        dplyr::select(model, year, age, weight, pred_weight, terminal_train_year, projection, peel_id)
    } ## end test_list

    test_list[[i]] <- test_list[[i]] %>%
      dplyr::mutate(YID = paste0('y+', year - terminal_train_year))
  } ## end peels

  # Rename ----
  names(peel_list) <- 1:peels
  names(pred_weight_list) <- 1:peels
  names(train_list) <- 1:peels
  names(test_list) <- 1:peels
  non_converged <- unique(non_converged)


  # Forecast performance metrics ----
  # * Create master dataframes ----
  test_data_combined <- do.call("rbind", lapply(1:length(test_list), function(i) { test_list[[i]]})) %>%
    dplyr::filter(!model %in% non_converged)

  if(!forecast){
    test_data_combined <- test_data_combined %>%
      dplyr::mutate(YID = "hindcast")
  }

  train_data_combined <- do.call("rbind", lapply(1:length(train_list), function(i) { train_list[[i]]})) %>%
    dplyr::filter(!model %in% non_converged)

  predicted_weight_combined <- do.call("rbind", lapply(1:length(pred_weight_list), function(i) { pred_weight_list[[i]]})) %>%
    dplyr::filter(!model %in% non_converged)

  # * Calculate average RMSE by model and forecast year ----
  rmse_table   <-  test_data_combined %>%
    filter(peel_id > 1) %>%
    dplyr::summarise(RMSE = sqrt(mean((weight - pred_weight)^2, na.rm = TRUE)),
                     .by = c(model, YID)) %>%
    arrange(YID, RMSE) %>%
    dplyr::select(YID, model, RMSE)

  # * Calculate average RMSE by age, model, and forecast year ----
  rmse_table_by_age <- test_data_combined %>%
    dplyr::filter(peel_id > 1) %>%
    dplyr::group_by(age, YID) %>%
    dplyr::slice_sample(n = 50 ) %>%
    ungroup() %>%
    dplyr::summarise(RMSE = sqrt(mean((weight - pred_weight)^2, na.rm = TRUE)),
                     .by = c(model, YID, age)) %>%
    filter(!is.na(RMSE)) %>%
    arrange(YID, age, RMSE) %>%
    dplyr::select(YID, age, model, RMSE)

  # Pick best model ----
  # * Based on overall RMSE ----
  best_mods <- rmse_table %>%
    dplyr::group_by(YID) %>%
    dplyr::filter(RMSE == min(RMSE, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::select(YID, model, RMSE) %>%
    suppressWarnings()

  # - Expand if biennial or triennial data
  best_mods <- best_mods %>%
    dplyr::full_join(rmse_table %>%
                       dplyr::distinct(YID)) %>%
    dplyr::arrange(YID) %>%
    dplyr::mutate(model = ifelse(is.na(model), lead(model), model)) %>%
    dplyr::mutate(model = ifelse(is.na(model), lead(model), model))

  # * Based on age-specific RMSE ----
  best_mods_by_age <- rmse_table_by_age %>%
    dplyr::group_by(YID, age) %>%
    dplyr::filter(RMSE == min(RMSE, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::select(YID, age, model, RMSE)

  # - Expand if biennial or triennial data
  best_mods_by_age <- best_mods_by_age %>%
    dplyr::full_join(rmse_table_by_age %>%
                       dplyr::distinct(YID, age)) %>%
    dplyr::arrange(age, YID) %>%
    dplyr::mutate(model = ifelse(is.na(model), lead(model), model)) %>%
    dplyr::mutate(model = ifelse(is.na(model), lead(model), model)) %>%
    dplyr::select(YID, age, model, RMSE)


  # Output ----
  # (lookup projection in first peel)
  plot_list <- list()
  return_list <- list(training_data = train_data_combined,
                      test_data = test_data_combined,
                      predictions = predicted_weight_combined,
                      rmse_table = rmse_table,
                      rmse_table_by_age = rmse_table_by_age,
                      best_mods = best_mods,
                      best_mods_by_age = best_mods_by_age)


  if(forecast){
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
    train_data_combined$year <- as.numeric(train_data_combined$year)
    plotdf0 <- bind_rows(train_data_combined, test_data_combined) %>%
      distinct()

    # * Plot historical data and fits ----
    for(mod in unique(plotdf0$model)){
      plotdf <- plotdf0 %>%
        filter(model == mod & year > (max(data$year)-(n_proj_years*peels)) ) ## truncate historical years

      plot_list[[mod]] <- ggplot2::ggplot(data= NULL, aes(x = age)) +
        ggplot2::geom_point(data = plotdf,
                            alpha = 0.05,
                            size = 2,aes(y = weight, color = factor(projection))) +
        ggplot2::geom_line(data = plotdf,
                           aes(y = pred_weight,
                               color = factor(projection),
                               group = interaction(year, peel_id)))+
        ggplot2::scale_color_manual(values = c('grey22','blue'))+
        ggplot2::facet_grid(peel_id ~ year, scales = 'free_y') +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = 'none') +
        ggplot2::labs(x = 'age', y = 'weight', title = mod) +
        theme_classic()
    }

    # * single plot with best projection
    plotdf <- plotdf0 %>%
      dplyr::filter(model  %in% best_mods$model &
                      year > (max(data$year)-(n_proj_years*peels)) ) ## truncate historical years

    if(length(unique(best_mods$model))>1) colVec <- c('grey30','dodgerblue','grey30','blue')
    if(length(unique(best_mods$model))==1) colVec <- c('grey22','blue')

    plot_list[[length(plot_list)+1]] <- ggplot2::ggplot(data= NULL, aes(x = age)) +
      ggplot2::geom_point(data = plotdf,
                          alpha = 0.05,
                          color = 'grey22',
                          size = 2,aes(y = weight)) +
      ggplot2::geom_line(data = plotdf,
                         lwd = 1,
                         aes(y = pred_weight,
                             color = interaction(projection, model),
                             group = interaction(model, year, peel_id)))+
      ggplot2::scale_color_manual(values = colVec)+
      ggplot2::facet_grid(peel_id ~ year, scales = 'free_y') +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = 'bottom') +
      ggplot2::labs(x = 'age', y = 'weight', title = 'best model(s)', color = '') +
      theme_classic()


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
           color = '') +
      theme_classic()

    names(plot_list) <- c(names(peel_list[[peels]]), "projection")

    return_list <- c(return_list, list(best_forecast = projected_waa,
                                       plots = plot_list)
    )
  }

  # Return ----
  return(return_list)
} ## end function


