## Generate and save selectivity curves

logistic <- function(ages=1:20, a50=10, delta=5){
  1/(1+exp(-log(19)*(ages-a50)/delta))
}

## GOA POP: https://www.npfmc.org/wp-content/PDFdocuments/SAFE/2023/GOApop.pdf
pop_sel <- logistic(ages=1:29,a50=5.80564854327391, delta=6.32499372501583)
save(pop_sel, file = file.path('data', paste0(species_name, "_pop-selex.rda")))
## GOA FLATHEAD: https://apps-afsc.fisheries.noaa.gov/Plan_Team/2022/GOAflathead.pdf
## had to approximate from double normal
fhs_sel <- logistic(ages=1:29, a50=4.420456 , delta=6.212803)
save(fhs_sel, file = file.path('data', paste0(species_name, "_flathead-selex.rda")))

## GOA PCod: https://www.npfmc.org/wp-content/PDFdocuments/SAFE/2024/GOApcod.pdf
## TODO: the assessment has this length based, do we want to proceed?
