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
