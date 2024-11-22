# This code simulates weight-at-age data for individuals from a population where there are multiple sub-groups
# based on a von Bertalanfy growth curve

# Weight.obs.i = Winf.i * (1-exp(-k.i * (age.obs.i - t0.i))) + eps.i
# Winf_i = exp(mu.winf + winf.group)
# K_i  = exp(mu.k + k.group)
# t0_i = mu.t0 + t0.group

# The group level parameters are assumed to be MVN with mean 0
# Eps.i is the normally distributed error term


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
G=10
nsamples = 300
group <- sample(1:6, nsamples, replace = TRUE) # Group for individual X
N <- nsamples # Total number of samples


## Mu VBGM hyperparameters ----
mu.Winf = 500
mu.k = 0.3
mut.t0 = -0.5
mu.parms <- c(mu.Winf, mu.k, mut.t0)
sigma = 10 # Observation error


## Group level random effects ----
sigma.group = c(0.1, 0.05, 0.2)
rho = 0.3 # Correlation between group level parameters
cor.group.mat = matrix(rho, 3, 3)
diag(cor.group.mat) <- 1
cov.group.mat <- diag(sigma.group) %*% cor.group.mat %*% diag(sigma.group) # Get covariance


## Simulate parameters for groups----
# - Empty matrix and vectors to fill with parameters and data, respectively
group.param.mat <- group.re.mat <- matrix(NA,G,3,byrow = T)

# - Random effects
colnames(group.re.mat) <- c("log.Winf.group.re", "log.k.group.re", "t0.group.re")

# - On VBGF scale
colnames(group.param.mat) <- c("Winf.group", "k.group", "t0.group")


# - Simulate group level parameters
for(i in 1:G){
  # - Sim from mvnorm
  group.re.mat[i,] <- mvrnorm(1, rep(0,3), cov.group.mat)

  # - Convert to parameter space
  group.param.mat[i,1:2] <- mu.parms[1:2] * exp(group.re.mat[i,1:2]) # Log to natural scale
  group.param.mat[i,3] <- mu.parms[3] + group.re.mat[i,3]
}


## Simulate weight-at-age data ----
ages = seq(from=1,to=20, by = .05)
age = c()
weight = c()
for(j in 1:N) {
  age[j] = sample(ages, 1) # Sample random age from age range
  weight[j] = (group.param.mat[group[j],1] * (1 - exp(-group.param.mat[group[j],2]*(age[j]-group.param.mat[group[j],3])))) + rnorm(1,0,sigma) # Add normally distributed random error
  #FIXME: may want lognormal error
}


# - Assign data to data frame
dat = data.frame(age = age, weight = weight/100, year = as.factor(group))
dat <- dat[which(dat$weight > 0),] # Make sure all weights are positive
dat <- dat %>% arrange(year, age, weight)

# - Plot the data
cols <- c("#86BBD8","#2F4858", "#F6AE2D", "#F26419", "#E86A92", "#57A773") # Colors for the VBGF lines (Females/Males)
ggplot(dat, aes(x = age, y = weight, colour = year)) +
  geom_point(size = 2) +
  scale_colour_manual(values=cols)


## Fit estimation models ----
data = data.frame(age = round(age), weight = weight/100, year = group)

wtagere <- FitWtAgeRE(
  data = data,
  weights=NULL,
  # - Number of projection years
  n_proj_years = 2
)

gmrf <- FitGMRF(
  data = data,
  weights=NULL,
  # - Number of projection years
  n_proj_years = 2
)
