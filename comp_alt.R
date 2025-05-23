# This file runs the simulations shown in Section 3.2 of the paper
# "Improving the (approximate) sequential probability ratio test by avoiding overshoot"

# Load the boosting functions
source("boosting_functions.R")

# Set simulation parameters
alpha <- 0.05
mu_N <- 0
m <- 10000
n <- 10000
mus <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

# Initialize vectors for stopping times and power
mean_stop_sprt <- rep(0, length(mus))
mean_stop_boosted <- rep(0, length(mus))
power_sprt <- rep(0, length(mus))
power_boosted <- rep(0, length(mus))

# Set seed for reproducibility
set.seed(123)

count <- 1
for (mu_A in mus) {
  stop_sprt <- rep(n, m) # Stopping times for a specific mu_A
  stop_boosted <- rep(n, m) # Stopping times for a specific mu_A
  decision_sprt <- rep(0, m) # decisions for a specific mu_A
  decision_boosted <- rep(0, m) # decisions for a specific mu_A

  for (j in 1:m) {
    data <- rnorm(n, mean = mu_A, sd = 1) # Create normally distributed data with mean mu_A and variance 1
    m_t <- 1 / alpha # Used for boosting
    b <- rep(1, n) # boosting factors
    lr <- rep(1, n) # Likelihood ratio
    lr_boosted <- rep(1, n) # boosted likelihood ratio

    mu_A_pred <- mu_N # Estimate for mu

    for (i in 1:n) {
      # Calculate likelihood ratio
      lr[i] <- dnorm(data[i], mu_A_pred, 1) / dnorm(data[i], mu_N, 1)
      # Calculate boosting factor
      if (m_t > 1) {
        b[i] <- e_boosted(m_t, mu_A_pred)
      }
      # Calculate boosted likelihood ratio
      lr_boosted[i] <- b[i] * lr[i]

      # Calculate m_t (used for boosting)
      m_t <- 1 / (alpha * prod(lr_boosted[1:i]))
      # Calculate estimate for mu
      mu_A_pred <- max(mu_N, (sum(data[1:i]) + mu_N) / i)

      # Set stopping times and decisions
      if (prod(lr_boosted[1:i]) >= 1 / alpha & stop_boosted[j] == n) {
        stop_boosted[j] <- i
        decision_boosted[j] <- (prod(lr_boosted[1:i]) >= 1 / alpha)
      }
      if (prod(lr[1:i]) >= 1 / alpha & stop_sprt[j] == n) {
        stop_sprt[j] <- i
        decision_sprt[j] <- (prod(lr[1:i]) >= 1 / alpha)
        break
      }
    }
  }
  # Save the mean power and stops for each mu_A
  power_boosted[count] <- mean(decision_boosted)
  power_sprt[count] <- mean(decision_sprt)
  mean_stop_boosted[count] <- mean(stop_boosted)
  mean_stop_sprt[count] <- mean(stop_sprt)
  count <- count + 1
}

# Save the results
results_df <- data.frame(idx = mus, mean_stop_boosted, mean_stop_sprt, power_boosted, power_sprt)
save(results_df, file = "results/comp_alt.rda")
