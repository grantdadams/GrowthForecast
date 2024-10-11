# Load data and packages
library(GrowthForecast)
data(LWA)

# Look at data
sort(unique(LWA$REGION))
sort(unique(LWA$CN))
sort(unique(LWA$YEAR))

# Select species and region
pollock <- LWA %>%
  dplyr::filter(CN == "walleye pollock" & REGION == "GOA")
