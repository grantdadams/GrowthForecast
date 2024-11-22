# Load data and packages
library(GrowthForecast)
data(LWA)

# Look at data
sort(unique(LWA$REGION))
sort(unique(LWA$CN))
sort(unique(LWA$YEAR))

# Select species and region
pollock <- LWA %>%
  dplyr::filter(CN == "walleye pollock" & REGION == "GOA") %>%
  dplyr::select(YEAR, age, weight, Temp) %>%
  dplyr::mutate(weight = weight/1000) # Grams to kg

colnames(pollock) <- tolower(colnames(pollock))

# Fit models ----
FitWtAgeRE(
  data = pollock,
  weights=NULL,
  # - Number of projection years
  n_proj_years = 2
)
