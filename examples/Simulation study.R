# Load data and packages
library(GrowthForecast)
library(dplyr)
library(ggplot2)
library(tidyr)
library(FishLife) # devtools::install_github("James-Thorson-NOAA/FishLife")
data("FishBase_and_Morphometrics")

# Simulation set-up
nsim = 300
edge_names = c( FishBase_and_Morphometrics$tree$tip.label,
                FishBase_and_Morphometrics$tree$node.label[-1] ) # Removing root



# Life history 1 (Pcod) ----
cod_pars <- get_growthPars(species_name = 'pcod')
cod_pars$mu[1] <- cod_pars$mu[1]/1000 # G to KG
cod_pars$mu <- cod_pars$mu[1:3]
cod_pars$vcov <- cod_pars$vcov[1:3,1:3]

which_g = match( "Gadus macrocephalus", edge_names )
lifehist_table = cbind(
  Mean = FishBase_and_Morphometrics$beta_gv[which_g,],
  SE = sqrt(diag(FishBase_and_Morphometrics$Cov_gvv[which_g,,]))
) %>%
  as.data.frame()

max_age <- exp(lifehist_table$Mean[1])
# get_growthPars(species_name = 'pop')
# get_growthPars(species_name = 'flathead')


# * Static ----
static_data = simulate_weight(
  nyrs = 40,
  nsamples = 10000,
  nages = ceiling(max_age),
  mu = cod_pars$mu, # Winf, K, t0
  trend_beta = c(0, 0, 0),
  vcov = diag(c(0, 0, 0)),
  sigma_obs = 0.1,
  rho_ar1 = 0,
  seed = 1234
)

sim_data <- static_data$data %>%
  mutate(age = round(age))

# - Plot the data
ggplot(sim_data, aes(x = age, y = weight, colour = year)) +
  geom_point(size = 2) +
  scale_color_continuous()

# - Forecast skill
simforecast <- ForecastGrowth(form = formula(weight~age+year), data = sim_data, n_proj_years = 2, peels = 10)
simforecast$best_mods


# * N(0, 1) noise ----
vbgm_n01 = simulate_weight(
  nyrs = 40,
  nsamples = 10000,
  nages = ceiling(max_age),
  mu = cod_pars$mu, # Winf, K, t0
  trend_beta = c(0, 0, 0),
  vcov = cod_pars$vcov,
  sigma_obs = 0.2,
  rho_ar1 = 0,
  seed = 1234
)

sim_data <- vbgm_n01$data %>%
  mutate(age = round(age))

# - Plot the data
ggplot(sim_data, aes(x = age, y = weight, colour = year)) +
  geom_point(size = 2) +
  scale_color_continuous()

# - Forecast skill
simforecast <- ForecastGrowth(form = formula(weight~age+year), data = sim_data, n_proj_years = 2, peels = 10)
simforecast$best_mods


# * AR1 ----
vbgm_ar1 = simulate_weight(
  nyrs = 40,
  nsamples = 10000,
  nages = ceiling(max_age),
  mu = cod_pars$mu, # Winf, K, t0
  trend_beta = c(0, 0, 0),
  vcov = cod_pars$vcov,
  sigma_obs = 0.2,
  rho_ar1 = 0.95,
  seed = 1234
)

sim_data <- vbgm_ar1$data %>%
  mutate(age = round(age))

# - Plot the data
ggplot(sim_data, aes(x = age, y = weight, colour = year)) +
  geom_point(size = 2) +
  scale_color_continuous()

# - Forecast skill
simforecast <- ForecastGrowth(form = formula(weight~age+year), data = sim_data, n_proj_years = 2, peels = 10)
simforecast$best_mods


# * Trend ----
vbgm_trend = simulate_weight(
  nyrs = 40,
  nsamples = 10000,
  nages = ceiling(max_age),
  mu = cod_pars$mu, # Winf, K, t0
  trend_beta = c(-0.2, 0.2, 0),
  vcov = cod_pars$vcov,
  sigma_obs = 0.2,
  rho_ar1 = 0.95,
  seed = 1234
)

sim_data <- vbgm_trend$data %>%
  mutate(age = round(age))

# - Plot the data
ggplot(sim_data, aes(x = age, y = weight, colour = year)) +
  geom_point(size = 2) +
  scale_color_continuous()

# - Forecast skill
simforecast <- ForecastGrowth(form = formula(weight~age+year), data = sim_data, n_proj_years = 2, peels = 10)
simforecast$best_mods
