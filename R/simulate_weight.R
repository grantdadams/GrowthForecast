

#' Simulate weight-at-age data
#'
#' @param nyrs Number of years
#' @param nsamples Number of samples across years and ages
#' @param nages Number of ages (starting at 1)
#' @param mu Vector of mean asymptoptic weight, growth coefficient, and t0 on natural scale
#' @param trend_beta Percent linear increase or decrease in mu: \code{mu * seq(1, 1 + trend_beta, length.out = nyrs)}. Vector of length 3.
#' @param vcov Variance-covariance matrix for von Bertalanfy parameters (Winf, K, t0)
#' @param sigma_obs Variance of lognormal observation error "sdlog"
#' @param rho Correlation for von Bertalanfy growth parameters (Winf, K, t0) MVN annual deviates
#' @param rho_ar1 Vector of length 3 of AR1 correlation for von Bertalanfy growth parameters (Winf, K, t0)
#' @param seed
#'
#' @description
#' This code simulates weight-at-age data for individuals from a population where there are multiple years
#' based on a von Bertalanfy growth curve with environmental processes on the parameters:
#'
#' \deqn{Weight_i = Winf_{yr[i]} * (1 - exp(k_{yr[i]} * (age_i - t0_{yr[i]})))}
#'
#' \deqn{VGBM_{yr} = AR1 * VGBM_{yr-1} + MVN(0, \Sigma)}
#'
#' @return
#' @export
#'
#' @examples
#'
#' dat = simulate_weight() %>%
#'  mutate(year = as.factor(year)) %>%
#'  arrange(year, age, weight)

#' # - Plot the data
#' ggplot(dat, aes(x = age, y = weight, colour = year)) +
#'  geom_point(size = 2) +
#'  scale_color_discrete()

simulate_weight <- function(
    nyrs = 10,
    nsamples = 1000,
    nages = 20,
    mu = c(5, 0.3, -0.5), # Winf, K, t0
    trend_beta = c(0, 0, 0),
    vcov = diag(c(0.05, 0.05, 0.2)),# Variance-covariance of vgbm parameters (mu)
    sigma_obs = 0.1,
    rho_ar1 = 0.95, # Time series rho
    seed = 1234){
  set.seed(seed)

  ## TODO decide whether to call get_growthPars here, pass species as arg
  ## similarly look up or read in nages and selex vector

  ## Year AR1 ----
  exponent <- abs(matrix(1:nyrs - 1, nrow = nyrs, ncol = nyrs, byrow = TRUE) -
                    (1:nyrs - 1))
  cor.ar1.mat <- rho_ar1^exponent
  cov.ar1.mat <- diag(1, nyrs) %*% cor.ar1.mat %*% diag(1, nyrs) # Get covariance


  # Trend ----
  # - Multiply trend by mu and convert to vector
  trend_df <- sapply(trend_beta, function(x) seq(1, 1 + x, length.out = nyrs))
  mut_df <- t(t(trend_df) * mu)
  mut_df[,1:2] <- log(mut_df[,1:2]) # Move winf and k to log scale
  mu_vec <- c(mut_df) # Convert to vector

  # - Create "indices"
  trend_df <- cbind(1:nyrs, trend_df)
  colnames(trend_df) <- c("year", "winf_trend", "k_trend", "t0_trend")


  # Simulate year specific parameters ----
  cov.mat <- vcov %x% cov.ar1.mat

  year.param.mat <- MASS::mvrnorm(1, mu_vec, cov.mat)
  year.param.mat <- matrix(year.param.mat, nrow = nyrs, ncol = 3, byrow = FALSE) # Rearrange
  year.param.mat[,1] <- exp(year.param.mat[,1]) # Move to natural scale
  year.param.mat[,2] <- exp(year.param.mat[,2])
  colnames(year.param.mat) <- c("Winf.year", "k.year", "t0.year")


  ## Simulate weight-at-age data ----
  year <- sample(1:nyrs, nsamples, replace = TRUE) # year for individual X
  age = runif(nsamples, 1, nages) # Sample random age from age range

  true_weight = (year.param.mat[year, 1] * (1 - exp(-year.param.mat[year, 2]*(age-year.param.mat[year, 3]))))
  obs_weight = true_weight * rlnorm(nsamples, meanlog = 0, sdlog = sigma_obs)

  # Expected weight at age ----
  year_age <- expand.grid(year = 1:nyrs, age = 1:nages)
  true_weight_mat <- (year.param.mat[year_age$year, 1] * (1 - exp(-year.param.mat[year_age$year, 2]*(year_age$age-year.param.mat[year_age$year, 3]))))
  true_weight_mat <- matrix(true_weight_mat, nrow = nyrs)


  # Return data object ----
  data <- data.frame(weight = obs_weight, age = age, year = as.numeric(year), true_weight = true_weight)
  data <- merge(data, trend_df, by = "year")
  return(list(data = data, ewaa = true_weight_mat))
}


