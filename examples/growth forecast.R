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

pollockforecast <- GrowthForecast::ForecastGrowth(form = formula(weight~age+year), data = pollock, n_proj_years = 4, peels = 3)
pollockforecast$best_mods[c(2,4),]
pollockforecast$rmse_table
pollockforecast$best_forecast
pollockforecast$plots

# Simulated populations ----
# - for final SSB comparison

# Initialize parameters
years <- 250
ages <- 20
M <- 0.2
init_pop <- 1e4
# - maturity curve, assume a50 = 8
maturity <- 1 / (1 + exp(-0.5 * (1:ages - 8)))

# Create a matrix to store the population projections
population_matrix <- matrix(0, nrow = years, ncol = ages)

# Set the initial population for age 1 at year 1
population_matrix[1, ] <- exp(-M*c(1:20))*(init_pop)

# Loop through each year to project the population forward
for (year in 2:years) {
  # Randomly sample a new value for age 1 using a normal distribution
  population_matrix[year, 1] <- 5000#rnorm(1, mean = init_pop, sd = 500)

  # Loop through each age to update the population
  for (age in 2:ages) {
    # Decrement the previous year's population at age-1 by exp(-0.2*age)
    population_matrix[year, age] <- population_matrix[year - 1, age - 1] * exp(-M * (age - 1))
  }
}


# - calculate true SSB and extract last 15 years
# - TODO make this longer and use earliest DAT waa to hindcast for completeness
population_matrix <- tail(population_matrix,10)

# - get putative WAA curve to make true OM/Assessment weights
sim_data <- simulate_weight(seed = 2,
                            nyrs = 10,
                            mu = c(5, 0.15, -0.5), # Winf, K, t0
                            sigma_obs = 0.15,
                            rho_ar1 = 0.9999  # Time series rho

)

waa_true <- sim_data$ewaa

# - calculate true SSB
true_ssb <- round((population_matrix *  waa_true) %*% maturity) %>%
  data.frame() %>%
  mutate(year = 1:nrow(waa_true)) %>%
  dplyr::select(year, ssb_true = ".")

plot(true_ssb$year, true_ssb$ssb_true, type = 'l')

# Simulated data ----
# - Simulate data
set.seed(123)

data <-  sim_data$data %>% dplyr::mutate(age = round(age))

# - Plot the data
ggplot(data, aes(x = age, y = weight, colour = year)) +
  geom_point(size = 2) +
  scale_color_continuous()

# - Test forecasting
simforecast <- ForecastGrowth(form = formula(weight~age+year),
                              data = data,
                              n_proj_years = 2,
                              peels = 3,
                              maturity_vec =  maturity)
simforecast$best_mods

simforecast$plots



# - plot ssb comparisons

predicted_values <- population_matrix %>%
  data.frame() %>%
  mutate(year = row_number(.)) %>%
  reshape2::melt(id = 'year', value.name = 'numbers') %>%
  mutate(age = as.numeric(gsub("X","", variable))) %>%
  select(year,age,numbers) %>%
  arrange(year) %>%
  merge(.,simforecast$all_predictions %>%
          filter(projection),
        by = c('year','age'), all.y = FALSE) %>%
  merge(.,cbind(age = 1:20, maturity), by = 'age') %>%
  dplyr::mutate(YID = paste0('y+',year - terminal_train_year)) %>%
  mutate(pred_sb = numbers*pred_weight *maturity) %>%
  dplyr::summarise(pred_ssb = sum(pred_sb),.by = c(projection,year, YID,model)) %>%
  merge(., simforecast$best_mods, by = c('YID', 'model'), all.x = TRUE) %>%
  mutate(best = ifelse(is.na(RMSE), FALSE, TRUE)) %>%
  select(-RMSE)

ggplot(data = NULL, aes(x = year)) +
  geom_line(data = true_ssb, aes(y = ssb_true), lwd = 1) +
  geom_point(data = predicted_values, aes(y = pred_ssb,   color = interaction(model, best)) ,size = 2, shape = 4) +
  # geom_point(data = subset(predicted_values,best), aes(y = pred_ssb,   color = interaction(model, best)) ,size = 2, shape = 4) +
  facet_wrap(~YID) +
  scale_color_manual(values = c('red',
                                grey.colors(6, start = 0.5, end = 0.7),
                                'blue','purple'))+
  theme_minimal() +
  labs(x = 'year', y= 'SSB', title = 'forecasted vs true ssb')


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
