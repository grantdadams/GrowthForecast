# Load data and packages
library(GrowthForecast)
library(tidymodels)
library(agua)
library(ggplot2)
data(LWA)

# Look at data
sort(unique(LWA$REGION))
sort(unique(LWA$CN))
yrs <- sort(unique(LWA$YEAR))

# Select species and region
pollock <- LWA %>%
  dplyr::filter(CN == "walleye pollock" & REGION == "GOA")



theme_set(theme_bw())
h2o_start()

# Example ----
data(concrete)
set.seed(4595)
concrete_split <- initial_split(concrete, strata = compressive_strength)
concrete_train <- training(concrete_split)
concrete_test <- testing(concrete_split)

# run for a maximum of 120 seconds
auto_spec <-
  auto_ml() %>%
  set_engine("h2o", max_runtime_secs = 120, seed = 1) %>%
  set_mode("regression")

normalized_rec <-
  recipe(compressive_strength ~ ., data = concrete_train) %>%
  step_normalize(all_predictors())

auto_wflow <-
  workflow() %>%
  add_model(auto_spec) %>%
  add_recipe(normalized_rec)

auto_fit <- fit(auto_wflow, data = concrete_train)
extract_fit_parsnip(auto_fit)


# WAA ----
# * 1-year forecast
pollock <- pollock %>%
  dplyr::mutate(Cohort = YEAR - age)

yrs <- sort(unique(pollock$YEAR))
pollock_train <- pollock %>%
  filter(YEAR < max(yrs))
pollock_test <-  pollock %>%
  filter(YEAR == max(yrs))

# run for a maximum of 120 seconds
auto_spec <-
  auto_ml() %>%
  set_engine("h2o", max_runtime_secs = 120, seed = 1) %>%
  set_mode("regression")

rec <-
  recipe(weight ~ age + Temp + YEAR + Cohort, data = pollock_train)

auto_wflow <-
  workflow() %>%
  add_model(auto_spec) %>%
  add_recipe(rec)

auto_fit_pk1 <- fit(auto_wflow, data = pollock_train)
extract_fit_parsnip(auto_fit_pk1)

# - Compare with mean
mean_wt <- pollock_train %>%
  group_by(age) %>%
  summarise(mean_weight = mean(weight, na.rm = TRUE))

pollock_test <- merge(pollock_test, mean_wt, by = "age")
sqrt(mean((pollock_test$weight - pollock_test$mean_weight)^2)) # RMSE
mean((pollock_test$weight - pollock_test$mean_weight)^2) # MSE
mean(abs(pollock_test$weight - pollock_test$mean_weight)) # MAE

# Should be weighted by numbers


# * 2-year forecast
pollock <- pollock %>%
  dplyr::mutate(Cohort = YEAR - age)

yrs <- sort(unique(pollock$YEAR))
pollock_train <- pollock %>%
  filter(YEAR < (max(yrs)-1))
pollock_test <-  pollock %>%
  filter(YEAR == max(yrs))

# run for a maximum of 120 seconds
auto_spec <-
  auto_ml() %>%
  set_engine("h2o", max_runtime_secs = 120, seed = 1) %>%
  set_mode("regression")

rec <-
  recipe(weight ~ age + Temp + YEAR + Cohort, data = pollock_train)

auto_wflow <-
  workflow() %>%
  add_model(auto_spec) %>%
  add_recipe(rec)

auto_fit_pk2 <- fit(auto_wflow, data = pollock_train)
extract_fit_parsnip(auto_fit_pk2)

# - Compare with mean
mean_wt <- pollock_train %>%
  group_by(age) %>%
  summarise(mean_weight = mean(weight, na.rm = TRUE))

pollock_test <- merge(pollock_test, mean_wt, by = "age")
sqrt(mean((pollock_test$weight - pollock_test$mean_weight)^2)) # RMSE
mean((pollock_test$weight - pollock_test$mean_weight)^2) # MSE
mean(abs(pollock_test$weight - pollock_test$mean_weight)) # MAE
