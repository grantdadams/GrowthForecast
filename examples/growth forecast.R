# Load data and packages
library(GrowthForecast)
library(dplyr)
library(ggplot2)
library(tidyr)

# Observed data ----
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

pollock %>%
  count(age, year) %>%
  arrange(year) %>%
  pivot_wider(names_from = year, values_from = n)

pollockforecast <- GrowthForecast::ForecastGrowth(form = formula(weight~age+year), data = pollock, n_proj_years = 2, peels = 3)
pollockforecast$rmse_table
pollockforecast$best_forecast

# Simulated data ----
# - Simulate data
set.seed(123)
data = simulate_weight(seed = 2,
                       nyrs = 10,
                       rho_ar1 = .99, # Time series rho
) %>%
  mutate(age = round(age))

# - Plot the data
ggplot(data, aes(x = age, y = weight, colour = year)) +
  geom_point(size = 2) +
  scale_color_continuous()

# - Test forecasting
simforecast <- ForecastGrowth(form = formula(weight~age+year), data = data, n_proj_years = 2, peels = 3)
simforecast$best_mods



# Simulated data w/ trend ----
# - Simulate data
set.seed(123)
data = simulate_weight(seed = 2,
                       nyrs = 10,
                       mu_trend = c(0.2, 0, 0), # 20% increase in Winf across years
                       vcov = diag(c(0, 0, 0)),
                       rho_ar1 = 0, # Time series rho
                       sigma_obs = 0.1
)

data$data <- data$data %>%
  mutate(age = round(age))

# - Plot the data
ggplot(data$data, aes(x = age, y = weight, colour = year)) +
  geom_point(size = 2) +
  scale_color_continuous()

ggplot(data$data %>%
         filter(round(age) == 20), aes(x = year, y = weight)) +
  geom_point(size = 2)

# - Test forecasting
simforecast <- ForecastGrowth(form = formula(weight~age+year), data = data, n_proj_years = 2, peels = 3)
simforecast$best_mods
