
#' Obtain and format growth parameters (mean and covariance) from FishLife
#'
#' @param species_name one of "pop", "pcod", "flathead" OR latin genus-species
#' @description
#' @return
#' @export
#'


get_growthPars <- function(species_name){
  if(species_name == 'pcod'){
    predict0 <- Plot_taxa( Search_species(Genus="Gadus", Species = "macrocephalus")$match_taxonomy)
  } else if(species_name == 'pop'){
    predict0 <- Plot_taxa( Search_species(Genus="Sebastes", Species = "alutus")$match_taxonomy)
  } else if(species_name == 'flathead'){
    predict0 <- Plot_taxa( Search_species(Genus="Hippoglossoides", Species = "elassodon")$match_taxonomy)
  } else{
    predict0 <- Plot_taxa( Search_species(Genus= strsplit(species_name,' ')[[1]][1],
                                          Species =  strsplit(species_name,' ')[[1]][2])$match_taxonomy)
  }

  mu.Linf <- get_mean_sd_fishlife(predict = predict0, trait = 'Loo')$mean
  mu.k <- get_mean_sd_fishlife(predict = predict0, trait = 'K')$mean
  mut.t0 <- get_mean_sd_fishlife(predict = predict0, trait = 'tm')$mean
  mu.Winf <- get_mean_sd_fishlife(predict = predict0, trait = 'Winfinity')$mean ## also known as a*linf^b
  mu.parms <- c(mu.Linf, mu.k, mut.t0, mu.Winf)
  npar <- length(mu.parms)

  ## Loo, k, t0, Winf
  cov.group.mat <-  cbind(predict0[[1]]$Cov_pred[c(1,2,4,3),c(1,2,4,3)])
  names(mu.parms) <- rownames(cov.group.mat)

  save(mu.parms, file = file.path('data', paste0(species_name, "_muparms.Rda")))
  save(cov.group.mat, file = file.path('data', paste0(species_name, "_covgroupmat.rda")))

  return(list(mu.parms, cov.group.mat))

}

## Not run:
# get_growthPars(species_name = 'flathead')


