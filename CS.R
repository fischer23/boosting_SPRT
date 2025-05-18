# This file runs the simulations shown in Section S.1 of the paper
# "Improving the (approximate) sequential probability ratio test by avoiding overshoot"

# Load the boosting functions
source("boosting_functions.R")
source("functions_CS.R")

# Set simulation parameters
mu_A <- 2
n <- 50
alpha <- 0.05
nsim <- 100
delta <- sqrt(8 * log(1 / alpha) / n)

# Initialize CS matrices
CSs <- matrix(, nrow = nsim, ncol = n)
boosted_CSs <- matrix(, nrow = nsim, ncol = n)

# Set seed for reproducibility
set.seed(123)

for (k in 1:nsim) {
  # Generate data
  data <- rnorm(n, mean = mu_A, sd = 1)
  # Calculate boosted CS
  boosted_CSs[k, ] <- boosted_CS(alpha, n, delta, data)
  # Calculate CS
  CSs[k, ] <- cummax(cumsum(data) / (1:n) - log(1 / alpha) / ((1:n) * delta) - delta / 2)
}

# Save results
results_df <- data.frame(idx = (1:n), CS = colMeans(CSs), boosted_CS = colMeans(boosted_CSs))
save(results_df, file = "results/CS.rda")
