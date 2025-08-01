
#' Obtain and format growth parameters (mean and covariance) from FishLife
#'
#' @param species_name one of "pop", "pcod", "flathead" OR latin genus-species
#' @description
#' @return
#' @export
#'


get_growthPars <- function(species_name){

  # Get bivariate info ----
  # devtools::install_github("James-Thorson-NOAA/FishLife")
  if(species_name == 'pcod'){
    predict0 <- FishLife::Plot_taxa( FishLife::Search_species(Genus="Gadus", Species = "macrocephalus")$match_taxonomy)
  } else if(species_name == 'pop'){
    predict0 <- FishLife::Plot_taxa( FishLife::Search_species(Genus="Sebastes", Species = "alutus")$match_taxonomy)
  } else if(species_name == 'flathead'){
    predict0 <- FishLife::Plot_taxa( FishLife::Search_species(Genus="Hippoglossoides", Species = "elassodon")$match_taxonomy)
  } else if(species_name == 'pollock'){
    predict0 <- FishLife::Plot_taxa( FishLife::Search_species(Genus="Gadus")$match_taxonomy)
  } else{
    predict0 <- FishLife::Plot_taxa( FishLife::Search_species(Genus= strsplit(species_name,' ')[[1]][1],
                                          Species =  strsplit(species_name,' ')[[1]][2])$match_taxonomy)
  }


  # Get mean ----
  mu.Linf <- GrowthForecast::reformat_FishLife(predict = predict0, trait = 'Loo')$mean
  mu.k <- GrowthForecast::reformat_FishLife(predict = predict0, trait = 'K')$mean
  mu.t0 <-  -0.1 # FishLife fixes it t0 at -0.1
  mu.Winf <- GrowthForecast::reformat_FishLife(predict = predict0, trait = 'Winfinity')$mean ## also known as a*linf^b
  mu.parms <- c(mu.Winf, mu.k, mu.t0, mu.Linf)
  npar <- length(mu.parms)

  # Covariance ----
  # * Winf, k, t0, Loo
  cov.group.mat <-  cbind(predict0[[1]]$Cov_pred[c(3,2,1),c(3,2,1)]) #FIXME: tm = t0?
  names(mu.parms) <- rownames(cov.group.mat)

  # Combine and save
  dataList <- list(mu = mu.parms, vcov = cov.group.mat)
  save(dataList, file = file.path('data', paste0(species_name, "_parameters.Rda")))

  return(dataList)
}


#' Bias correct mean and covariance values from FishLife (utility function)
#'
#' @param predict list object from of FishLife::PlotTaxa
#' @param trait specific trait eg K, M, Winfinity, Loo
#' @description
#' @return
#' @export


reformat_FishLife <- function(predict=predict0, trait='Loo'){
  log_mean<-predict[[1]]$Mean_pred[which(names(predict[[1]]$Mean_pred)==trait)]
  log_var <- diag(predict[[1]]$Cov_pred)[which(rownames(predict[[1]]$Cov_pred)==trait)]
  log_sd <- sqrt(log_var)
  mean_fl <- exp(log_mean + 0.5*log_var) ## get back the mean via bias correction
  sd_fl <- mean_fl*sqrt(exp(log_var)-1)
  cv_fl <- sd_fl/mean_fl ## the model
  return(list('mean'=mean_fl,  'sd'=sd_fl,'cv'= cv_fl))
}


## Not run:
# get_growthPars(species_name = 'flathead')


