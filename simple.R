# This file runs the simulations of the paper "Improving the (approximate) sequential probability ratio
# test by avoiding overshoot" that consider simple null and alternative hypotheses

# Load the boosting functions
source("boosting_functions.R")

# Set simulation parameters
alpha <- 0.05
mu_N <- 0
m <- 10000
n <- 10000
mus <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

# Initialize vectors for stopping times, power and type I error
mean_stop_sprt <- rep(0, length(mus))
mean_stop_sprt_ds <- rep(0, length(mus))
mean_stop_boosted <- rep(0, length(mus))
power_sprt <- rep(0, length(mus))
power_sprt_ds <- rep(0, length(mus))
power_boosted <- rep(0, length(mus))
mean_type_I_sprt <- rep(0, length(mus))
mean_type_I_sprt_ds <- rep(0, length(mus))
mean_type_I_boosted <- rep(0, length(mus))
mean_time_sprt <- rep(0, length(mus))
mean_time_sprt_ds <- rep(0, length(mus))
mean_time_boosted <- rep(0, length(mus))

# Set seed for reproducibility
set.seed(123)

count <- 1
for (mu_A in mus) {
  stop_sprt <- rep(n, m) # Stopping times for a specific mu_A
  stop_sprt_ds <- rep(n, m) # Stopping times for a specific mu_A
  stop_boosted <- rep(n, m) # Stopping times for a specific mu_A
  decision_sprt <- rep(0, m) # decisions for a specific mu_A
  decision_sprt_ds <- rep(0, m) # decisions for a specific mu_A
  decision_boosted <- rep(0, m) # decisions for a specific mu_A
  type_I_sprt <- rep(0, m) # Type I errors for a specific mu_A
  type_I_sprt_ds <- rep(0, m) # Type I errors for a specific mu_A
  type_I_boosted <- rep(0, m) # Type I errors for a specific mu_A
  time_sprt <- rep(0, m) # Stopping times for a specific mu_A
  time_sprt_ds <- rep(0, m) # Stopping times for a specific mu_A
  time_boosted <- rep(0, m) # Stopping times for a specific mu_A

  for (j in 1:m) {
    data <- rnorm(n, mean = mu_A, sd = 1) # Create normally distributed data with mean mu_A and variance 1
    m_t <- 1 / alpha # Used for boosting
    b <- rep(1, n) # boosting factors
    lr <- rep(1, n) # Likelihood ratio
    lr_boosted <- rep(1, n) # boosted likelihood ratio

    time_boosted[j] <- as.numeric(system.time({
      for (i in 1:n) {
        # Calculate likelihood ratio
        lr[i] <- dnorm(data[i], mu_A, 1) / dnorm(data[i], mu_N, 1)
        # Calculate boosting factor
        if (m_t > 1) {
          b[i] <- e_boosted(m_t, mu_A)
        }
        # Calculate boosted likelihood ratio
        lr_boosted[i] <- b[i] * lr[i]

        # Calculate m_t (used for boosting)
        m_t <- 1 / (alpha * prod(lr_boosted[1:i]))

        # Set stopping times and decisions
        if (prod(lr_boosted[1:i]) >= 1 / alpha & stop_boosted[j] == n) {
          stop_boosted[j] <- i
          decision_boosted[j] <- (prod(lr_boosted[1:i]) >= 1 / alpha)
          break
        }
      }
      })[3])
    
    time_sprt_ds[j] <- as.numeric(system.time({
      for (i in 1:n) {
        # Calculate likelihood ratio
        lr[i] <- dnorm(data[i], mu_A, 1) / dnorm(data[i], mu_N, 1)
        
        if (prod(lr[1:i]) >= 1 / (alpha * exp(mu_A * 0.583)) & stop_sprt_ds[j] == n) {
          stop_sprt_ds[j] <- i
          decision_sprt_ds[j] <- (prod(lr[1:i]) >= 1 / (alpha * exp(mu_A * 0.583)))
          break
        }
      }
    })[3])
    time_sprt[j] <- as.numeric(system.time({
      for (i in 1:n) {
        # Calculate likelihood ratio
        lr[i] <- dnorm(data[i], mu_A, 1) / dnorm(data[i], mu_N, 1)
        if (prod(lr[1:i]) >= 1 / alpha & stop_sprt[j] == n) {
          stop_sprt[j] <- i
          decision_sprt[j] <- (prod(lr[1:i]) >= 1 / alpha)
          break
        }
      }
    })[3])
    # Calculate the estimate for the type I error in each trial
    type_I_sprt[j] <- (1 / prod(lr[1:stop_sprt[j]])) * decision_sprt[j]
    type_I_sprt_ds[j] <- (1 / prod(lr[1:stop_sprt_ds[j]])) * decision_sprt_ds[j]
    type_I_boosted[j] <- (1 / prod(lr[1:stop_boosted[j]])) * decision_boosted[j]
  }
  # Save the mean power, stops and type I error for each mu_A
  power_boosted[count] <- mean(decision_boosted)
  power_sprt[count] <- mean(decision_sprt)
  power_sprt_ds[count] <- mean(decision_sprt_ds)
  mean_stop_boosted[count] <- mean(stop_boosted)
  mean_stop_sprt[count] <- mean(stop_sprt)
  mean_stop_sprt_ds[count] <- mean(stop_sprt_ds)
  mean_type_I_boosted[count] <- mean(type_I_boosted)
  mean_type_I_sprt[count] <- mean(type_I_sprt)
  mean_type_I_sprt_ds[count] <- mean(type_I_sprt_ds)
  mean_time_boosted[count] <- mean(time_boosted)
  mean_time_sprt[count] <- mean(time_sprt)
  mean_time_sprt_ds[count] <- mean(time_sprt_ds)
  count <- count + 1
}

# Save the results
results_df <- data.frame(
  idx = mus, mean_stop_boosted, mean_stop_sprt, mean_type_I_boosted, mean_type_I_sprt,
  power_boosted, power_sprt, mean_stop_sprt_ds, mean_type_I_sprt_ds, power_sprt_ds,
  mean_time_boosted, mean_time_sprt, mean_time_sprt_ds
)
save(results_df, file = "results/simple.rda")


############################### same simulations for alpha=0.01

alpha <- 0.01

# Initialize vectors for stopping times, power and type I error
mean_stop_sprt <- rep(0, length(mus))
mean_stop_sprt_ds <- rep(0, length(mus))
mean_stop_boosted <- rep(0, length(mus))
power_sprt <- rep(0, length(mus))
power_sprt_ds <- rep(0, length(mus))
power_boosted <- rep(0, length(mus))
mean_type_I_sprt <- rep(0, length(mus))
mean_type_I_sprt_ds <- rep(0, length(mus))
mean_type_I_boosted <- rep(0, length(mus))
mean_time_sprt <- rep(0, length(mus))
mean_time_sprt_ds <- rep(0, length(mus))
mean_time_boosted <- rep(0, length(mus))

# Set seed for reproducibility
set.seed(123)

count <- 1
for (mu_A in mus) {
  stop_sprt <- rep(n, m) # Stopping times for a specific mu_A
  stop_sprt_ds <- rep(n, m) # Stopping times for a specific mu_A
  stop_boosted <- rep(n, m) # Stopping times for a specific mu_A
  decision_sprt <- rep(0, m) # decisions for a specific mu_A
  decision_sprt_ds <- rep(0, m) # decisions for a specific mu_A
  decision_boosted <- rep(0, m) # decisions for a specific mu_A
  type_I_sprt <- rep(0, m) # Type I errors for a specific mu_A
  type_I_sprt_ds <- rep(0, m) # Type I errors for a specific mu_A
  type_I_boosted <- rep(0, m) # Type I errors for a specific mu_A
  time_sprt <- rep(0, m) # Stopping times for a specific mu_A
  time_sprt_ds <- rep(0, m) # Stopping times for a specific mu_A
  time_boosted <- rep(0, m) # Stopping times for a specific mu_A
  
  for (j in 1:m) {
    data <- rnorm(n, mean = mu_A, sd = 1) # Create normally distributed data with mean mu_A and variance 1
    m_t <- 1 / alpha # Used for boosting
    b <- rep(1, n) # boosting factors
    lr <- rep(1, n) # Likelihood ratio
    lr_boosted <- rep(1, n) # boosted likelihood ratio
    
    time_boosted[j] <- as.numeric(system.time({
      for (i in 1:n) {
        # Calculate likelihood ratio
        lr[i] <- dnorm(data[i], mu_A, 1) / dnorm(data[i], mu_N, 1)
        # Calculate boosting factor
        if (m_t > 1) {
          b[i] <- e_boosted(m_t, mu_A)
        }
        # Calculate boosted likelihood ratio
        lr_boosted[i] <- b[i] * lr[i]
        
        # Calculate m_t (used for boosting)
        m_t <- 1 / (alpha * prod(lr_boosted[1:i]))
        
        # Set stopping times and decisions
        if (prod(lr_boosted[1:i]) >= 1 / alpha & stop_boosted[j] == n) {
          stop_boosted[j] <- i
          decision_boosted[j] <- (prod(lr_boosted[1:i]) >= 1 / alpha)
          break
        }
      }
    })[3])
    
    time_sprt_ds[j] <- as.numeric(system.time({
      for (i in 1:n) {
        # Calculate likelihood ratio
        lr[i] <- dnorm(data[i], mu_A, 1) / dnorm(data[i], mu_N, 1)
        
        if (prod(lr[1:i]) >= 1 / (alpha * exp(mu_A * 0.583)) & stop_sprt_ds[j] == n) {
          stop_sprt_ds[j] <- i
          decision_sprt_ds[j] <- (prod(lr[1:i]) >= 1 / (alpha * exp(mu_A * 0.583)))
          break
        }
      }
    })[3])
    time_sprt[j] <- as.numeric(system.time({
      for (i in 1:n) {
        # Calculate likelihood ratio
        lr[i] <- dnorm(data[i], mu_A, 1) / dnorm(data[i], mu_N, 1)
        if (prod(lr[1:i]) >= 1 / alpha & stop_sprt[j] == n) {
          stop_sprt[j] <- i
          decision_sprt[j] <- (prod(lr[1:i]) >= 1 / alpha)
          break
        }
      }
    })[3])
    # Calculate the estimate for the type I error in each trial
    type_I_sprt[j] <- (1 / prod(lr[1:stop_sprt[j]])) * decision_sprt[j]
    type_I_sprt_ds[j] <- (1 / prod(lr[1:stop_sprt_ds[j]])) * decision_sprt_ds[j]
    type_I_boosted[j] <- (1 / prod(lr[1:stop_boosted[j]])) * decision_boosted[j]
  }
  # Save the mean power, stops and type I error for each mu_A
  power_boosted[count] <- mean(decision_boosted)
  power_sprt[count] <- mean(decision_sprt)
  power_sprt_ds[count] <- mean(decision_sprt_ds)
  mean_stop_boosted[count] <- mean(stop_boosted)
  mean_stop_sprt[count] <- mean(stop_sprt)
  mean_stop_sprt_ds[count] <- mean(stop_sprt_ds)
  mean_type_I_boosted[count] <- mean(type_I_boosted)
  mean_type_I_sprt[count] <- mean(type_I_sprt)
  mean_type_I_sprt_ds[count] <- mean(type_I_sprt_ds)
  mean_time_boosted[count] <- mean(time_boosted)
  mean_time_sprt[count] <- mean(time_sprt)
  mean_time_sprt_ds[count] <- mean(time_sprt_ds)
  count <- count + 1
}

# Save the results
results_df <- data.frame(
  idx = mus, mean_stop_boosted, mean_stop_sprt, mean_type_I_boosted, mean_type_I_sprt,
  power_boosted, power_sprt, mean_stop_sprt_ds, mean_type_I_sprt_ds, power_sprt_ds,
  mean_time_boosted, mean_time_sprt, mean_time_sprt_ds
)
save(results_df, file = "results/simple_alpha001.rda")
