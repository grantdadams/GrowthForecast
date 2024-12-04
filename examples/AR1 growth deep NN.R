# This code simulates weight-at-age data for individuals from a population where there are multiple years
# based on a von Bertalanfy growth curve with AR1 processes on the parameters


####################################################################################
# Simulate weight-at-age data ######################################################
####################################################################################
# Set Seed
library(GrowthForecast)
library(MASS)
library(ggplot2)
library(dplyr)
library(RTMB)
set.seed(1234)

## Specify data size ----
ngroup_hind = 5
ngroup_pred = 2
ngroup = ngroup_hind + ngroup_pred
nsamples = 1000  # Total number of samples


## Mu VBGM hyperparameters ----
mu.Winf = 500
mu.k = 0.3
mut.t0 = -0.5
mu.parms <- c(mu.Winf, mu.k, mut.t0)
sigma = 0.05 # Observation error


## Group level random effects ----
rho_ar1 <- 0.93 # Time series rho
sigma.group = c(0.1, 0.05, 0.2) # Variance in group level parameters
rho = 0.3 # Correlation between group level parameters
cor.group.mat = matrix(rho, 3, 3)
diag(cor.group.mat) <- 1
cov.group.mat <- diag(sigma.group) %*% cor.group.mat %*% diag(sigma.group) # Get covariance


## Simulate parameters for groups----
# - Empty matrix and vectors to fill with parameters and data, respectively
group.param.mat <- group.re.mat <- matrix(NA,ngroup,3,byrow = T)

# - Random effects
colnames(group.re.mat) <- c("log.Winf.group.re", "log.k.group.re", "t0.group.re")

# - On VBGF scale
colnames(group.param.mat) <- c("Winf.group", "k.group", "t0.group")


# - Simulate group level parameters
# -- Year 1
group.re.mat[1,] <- mvrnorm(1, rep(0,3), cov.group.mat)
group.param.mat[1,1:2] <- mu.parms[1:2] * exp(group.re.mat[1,1:2]) # Log to natural scale
group.param.mat[1,3] <- mu.parms[3] + group.re.mat[1,3]

# -- Year 2+
for(i in 2:ngroup){
  # - Sim from mvnorm
  group.re.mat[i,] <- mvrnorm(1, rep(0,3), cov.group.mat)

  # - Convert to parameter space
  group.param.mat[i,1:2] <- rho_ar1 * group.param.mat[i-1,1:2] * exp(group.re.mat[i,1:2]) # Log to natural scale
  group.param.mat[i,3] <- rho_ar1 * group.param.mat[i-1,3] + group.re.mat[i,3]
}


## Simulate weight-at-age data ----
nages = 20
group <- sample(1:ngroup, nsamples, replace = TRUE) # Group for individual X
age = runif(nsamples, 1, nages) # Sample random age from age range

weight = (group.param.mat[group, 1] * (1 - exp(-group.param.mat[group, 2]*(age-group.param.mat[group, 3])))) * exp(rnorm(nsamples, 0, sigma))/100


## Get expected curve -----
pred <- expand.grid(1:nages, 1:ngroup)
colnames(pred) <- c("age", "year")
pred$true_weight = (group.param.mat[pred$year,1] * (1 - exp(-group.param.mat[pred$year,2]*(pred$age-group.param.mat[pred$year,3]))))/100


# - Assign data to data frame
dat = data.frame(age = age, weight = weight, year = as.factor(group)) %>%
  arrange(year, age, weight)

# - Plot the data
ggplot(dat, aes(x = age, y = weight, colour = year)) +
  geom_point(size = 2) +
  scale_color_discrete()


## Fit deep NN models ----
library(neuralnet)
data = data.frame(weight = weight, age = round(age), year = group)

nn <- neuralnet(log(weight) ~ age+year,
               data = data %>%
                 filter(year <= ngroup_hind),
               hidden = c(5,5,5,5),
               linear.output = TRUE,
               stepmax = 1e6,
               lifesign = 'minimal',
               rep=1)

pred$pred_weight <- exp(predict(nn, newdata = pred))

# - TMB version
nn_rtmb <- fit_nn_rtmb(data %>%
                         filter(year <= 6) %>%
                         as.data.frame(),
                       nhidden_layer = 3, hidden_dim = 5)
pred$weight = 1
pred_nn_rtmb <- fit_nn_rtmb(pred, nhidden_layer = 3, hidden_dim = 5, input_par = nn_rtmb$parList)
pred$tmb_weight <- as.numeric(pred_nn_rtmb$report$output)
pred$weight = NULL




## Combine and plot ----
melted_pred <- merge(pred, data, by = c("age", "year"), all = TRUE) %>%
  mutate(year = as.factor(year))

ggplot(data = melted_pred, aes(x = age, y = true_weight, color = year, group = year)) +
  geom_line() +
  geom_line(aes(x = age, y = tmb_weight), lty = 2, size = 1.5) +
  geom_point(aes(x = age, y = weight)) +
  ylab("Weight")


gmrf <- FitGMRF(
  data = data %>%
    filter(year <= ngroup_hind),
  weights=NULL,
  # - Number of projection years
  n_proj_years = 2
)

# - Combine and plot
gmrf_pred <- as.data.frame(exp(gmrf[[4]]$report$ln_Y_at))
colnames(gmrf_pred) <- 1:7
gmrf_pred$age <- rownames(gmrf_pred)
gmrf_pred <- gmrf_pred %>%
  tidyr::pivot_longer(!age, values_to = "gmrf_wt", names_to = "year") %>%
  mutate(year = as.factor(year))

gmrf_pred <- merge(melted_pred, gmrf_pred, by = c("age", "year"), all = TRUE)

ggplot(data = gmrf_pred, aes(x = age, y = true_weight, color = year, group = year)) +
  geom_line() +
  geom_line(aes(x = age, y = gmrf_wt), lty = 2, size = 1.5) +
  geom_point(aes(x = age, y = weight)) +
  ylab("Weight")
