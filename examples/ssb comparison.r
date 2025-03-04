## Simulated comparison of SSB using growthForecast method
library(GrowthForecast)
library(dplyr)
library(ggplot2)
library(tidyr)

## NAA matrix
# Initialize parameters
years <- 250
ages <- 20
M <- 0.2
init_pop <- 1e4

# Create a matrix to store the population projections
population_matrix <- matrix(0, nrow = years, ncol = ages)

# Set the initial population for age 1 at year 1
population_matrix[1, ] <- exp(-M*c(1:20))*(init_pop)

# Loop through each year to project the population forward
for (year in 2:years) {
  # Randomly sample a new value for age 1 using a normal distribution
  population_matrix[year, 1] <- rnorm(1, mean = init_pop, sd = 0.05*init_pop)

  # Loop through each age to update the population
  for (age in 2:ages) {
    # Decrement the previous year's population at age-1 by exp(-0.2*age)
    population_matrix[year, age] <- population_matrix[year - 1, age - 1] * exp(-M * (age - 1))
  }
}

## maturity curve, assume a50 = 8
maturity <- 1 / (1 + exp(-0.5 * (1:ages - 8)))

## calculate true SSB and extract last 15 years
## TODO make this longer and use earliest DAT waa to hindcast for completeness

population_matrix <- tail(population_matrix,10)

## simulated WAA data used for forecasting
simulated_data <- simulate_weight(nyrs = 10, nsamples = 2500)
dat <- simulated_data$data %>%
 arrange(year, age, weight)
dat$age <- round(dat$age)
dat %>% summarise(n(),.by = age)

## get putative WAA curve to make true OM/Assessment weights
waa_true <- simulated_data$ewaa

## calculate true SSB
true_ssb <- (population_matrix *  waa_true) %*% maturity %>%
  data.frame() %>%
  mutate(year = 1:nrow(waa_true)) %>%
  dplyr::select(year, ssb_true = ".")

## do forecasting exercise
# simulated_forecast <- ForecastGrowth(data = dat)


simulated_forecast$best_mods

predicted_values <- population_matrix %>%
  data.frame() %>%
  mutate(year = row_number(.)) %>%
  reshape2::melt(id = 'year', value.name = 'numbers') %>%
  mutate(age = as.numeric(gsub("X","", variable))) %>%
  select(year,age,numbers) %>%
  arrange(year) %>%
  merge(.,simulated_forecast$all_predictions %>%
          filter(projection),
        by = c('year','age'), all.y = FALSE) %>%
  merge(.,cbind(age = 1:20, maturity), by = 'age') %>%
  dplyr::mutate(YID = paste0('y+',year - terminal_train_year)) %>%
  mutate(pred_sb = numbers*pred_weight *maturity) %>%
  dplyr::summarise(pred_ssb = sum(pred_sb),.by = c(projection,year, YID,model)) %>%
  merge(., simulated_forecast$best_mods, by = c('YID', 'model'), all.x = TRUE) %>%
  mutate(best = ifelse(is.na(RMSE), FALSE, TRUE)) %>%
  select(-RMSE)


## plot comparisons

ggplot(data = NULL, aes(x = year)) +
  geom_line(data = true_ssb, aes(y = ssb_true), lwd = 1) +
  geom_point(data = predicted_values, aes(y = pred_ssb,
                                          color = interaction(model, best)) ,size = 2, shape = 4) +
  # geom_point(data = subset(predicted_values,best), aes(y = pred_ssb), color = 'red', size = 2, shape = 4) +
  facet_wrap(~YID) +
  # scale_color_manual(values = c(grey.colors(8, start = 0.2, end = 0.6),'red','blue'))+
  theme_minimal()


