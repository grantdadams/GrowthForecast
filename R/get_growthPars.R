
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
  } else{
    predict0 <- Plot_taxa( Search_species(Genus="Sebastes", Species = "alutus")$match_taxonomy)
  } else{
    predict0 <- Plot_taxa( Search_species(Genus="Hippoglossoides", Species = "elassodon")$match_taxonomy)
  } else{
    predict0 <- Plot_taxa( Search_species(Genus= strsplit(species_name,' ')[[1]][1],
                                          Species =  strsplit(species_name,' ')[[1]][2])$match_taxonomy)
  }




}

