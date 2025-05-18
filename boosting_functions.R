# This file contains the functions for boosting the martingale factors in the paper
# "Improving the (approximate) sequential probability ratio test by avoiding overshoot"

### Function to compute boosting factors for Gaussian testing problems
## Input:
# m_t:           cutoff value >=1 for the truncation function (m_t=1/(M*\alpha)).
# delta:         Real parameter used for boosting the martingale factors. In case of the simple null H_i:X_i~N(mu_N,1)
#               vs. alternative H_i^A:X_i~N(mu_A,1) delta should be set to delta=mu_A-mu_N.

### Output:      boosting factor.

e_boosted <- function(m_t, delta) {
  b_factor <- function(b) {
    return(b * (1 - pnorm(delta / 2 - log(m_t / b) / delta)) + m_t * (1 - pnorm((log(m_t / b) + delta^2 / 2) / delta)) - 1)
  }
  return(uniroot(b_factor, lower = 0, upper = 100000000000)$root)
}

### Function to compute boosting for Vovk's conformal martingales
## Input:
# m_t:           cutoff value for the truncation function (m_t=1/(M*alpha)).
# kappa:         Parameter in (0,1) used for Vovk's calibrator (f_{\kappa}(u)=\kappa u^{\kappa-1})

### Output:      boosting factor.

e_boosted_calib <- function(m_t, kappa) {
  b_factor <- function(b) {
    return(b * max((1 - (m_t / (b * kappa))^(kappa / (kappa - 1))), 0) + m_t * min((m_t / (b * kappa))^(1 / (kappa - 1)), 1) - 1)
  }
  return(uniroot(b_factor, lower = 0, upper = 100000000000)$root)
}


### Function to compute boosting factors for Gaussian testing problems with predefined stop for futility \nu_t.
## Input:
# m_t:           cutoff value >=1 for the truncation function (m_t=1/(M*\alpha)).
# delta:         Real parameter used for boosting the martingale factors. In case of the simple null H_i:X_i~N(mu_N,1)
#               vs. alternative H_i^A:X_i~N(mu_A,1) delta should be set to delta=mu_A-mu_N.
# l_t:           lower cutoff value <1 for the truncation function (l_t=\nu_t/M).

### Output:      boosting factor.


e_boosted_futility <- function(m_t, l_t, delta) {
  b_factor <- function(b) {
    return(b * (pnorm(delta / 2 - log(l_t / b) / delta) - pnorm(delta / 2 - log(m_t / b) / delta)) + m_t * (1 - pnorm((log(m_t / b) + delta^2 / 2) / delta)) - 1)
  }
  return(uniroot(b_factor, lower = 0, upper = 100000000000)$root)
}


### Function to compute boosting factors and inverse boosting factors for Gaussian testing problems with desired type II error control.
## Input:
# LR:            Current value of the boosted likelihood ratio process.
# LR_inv:        Current value of the boosted invers likelihood ratio process.
# b_past:        Product of all previous boosting factors
# b_inv_past:    Product of all previous inverse boosting factors
# alpha:         Desired type I error probability
# beta:          Desired type II error probability
# delta:         Real parameter used for boosting the martingale factors. In case of the simple null H_i:X_i~N(mu_N,1)
#               vs. alternative H_i^A:X_i~N(mu_A,1) delta should be set to delta=mu_A-mu_N.
# delta_inv:     Real parameter used for boosting the inverse martingale factors. In case of the simple null H_i:X_i~N(mu_N,1)
#               vs. alternative H_i^A:X_i~N(mu_A,1) delta_inv should be set to delta_inv=mu_N-mu_A.

### Output:      3-dimensional vector containing the boosting factor, the inverse boosting factor and the stop for futility \nu_t.

e_boosted_type2 <- function(LR, LR_inv, b_past, b_inv_past, alpha, beta, delta, delta_inv) {
  obj <- function(b_full) {
    b_inv <- b_full[2]
    b <- b_full[1]
    return(-b - b_inv)
  }

  ineq_func <- function(b_full) {
    b_inv <- b_full[2]
    b <- b_full[1]
    return(c(
      b * (pnorm(delta / 2 - log(min(beta * b_inv * b_inv_past * b_past * b, 1 / alpha) / (b * LR)) / delta) - pnorm(delta / 2 - log(1 / (alpha * LR * b)) / delta)) + 1 / (alpha * LR) * (1 - pnorm((log(1 / (alpha * LR * b)) + delta^2 / 2) / delta)) - 1,
      b_inv * (pnorm(delta_inv / 2 - log(1 / (beta * LR_inv * b_inv)) / delta_inv) - pnorm(delta_inv / 2 - log(min(alpha * b * b_past * b_inv_past * b_inv, 1 / beta) / (b_inv * LR_inv)) / delta_inv)) + 1 / (beta * LR_inv) * (pnorm((log(1 / (beta * LR_inv * b_inv)) + delta_inv^2 / 2) / delta_inv)) - 1
    ))
  }
  result <- nloptr::nloptr(x0 = c(1, 1), eval_f = obj, lb = c(1, 1), eval_g_ineq = ineq_func, opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1.0e-8, "constraints_ineq" = 1.0e-8))
  if (max(ineq_func(result$solution)) > 1e-7 | min(result$solution) < 1) {
    result$solution <- c(1, 1)
  }
  return(c(result$solution, min(1 / alpha, beta * b_past * b_inv_past * result$solution[1] * result$solution[2])))
}
