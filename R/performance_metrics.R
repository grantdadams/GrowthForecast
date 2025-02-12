calc_pm <- function(observed = NULL, prediction = NULL){

  unlist(x$data)

  n <- length(observed)
  RSE <- sqrt(sum(observed-prediction)^2/ n)
  return(RSE)
}
