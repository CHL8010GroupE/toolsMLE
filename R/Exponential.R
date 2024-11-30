# MLE for Exponential Distribution
# Example of data: c(1.2, 2.4, 0.8, 1.5, 3.6), a vector that records observed data

# Exponential log-likelihood function
MLE.expll <- function(lambda, data) {
  if (lambda <= 0) {
    stop("Lambda must be greater than 0")
  }
  
  log_likelihood <- sum(log(lambda) - lambda * data)
  return(log_likelihood)
}

# Find the optimal MLE function
MLE.exponential <- function(data) {
  optimize(MLE.expll, c(1e-10, 1e10), data = data, maximum = TRUE)
}

# Likelihood Ratio Test
MLE.exponentialLRT <- function(data, lambda_test) {
  if (lambda_test <= 0) {
    stop("lambda_test must be greater than 0")
  }
  
  # Log-likelihood under null hypothesis (fixed lambda_test)
  log_test <- MLE.expll(lambda_test, data)
  
  # Log-likelihood under alternative hypothesis (MLE of lambda)
  estimated_parameter <- MLE.exponential(data)
  log_alter <- estimated_parameter$objective
  
  # Likelihood ratio statistic
  lR_statistic <- 2 * (log_alter - log_test)
  
  # Assume degree of freedom is 1
  p_value <- pchisq(lR_statistic, df = 1, lower.tail = FALSE)
  
  return(p_value)
}

# Goodness of Fit Test
MLE.exponentialGOF <- function(data) {
  # Estimate lambda using MLE
  loglikelihood <- MLE.exponential(data)
  lambda_hat <- loglikelihood$maximum
  
  # Divide data into bins for chi-square test
  observed_frequency <- table(cut(data, breaks = "Sturges", include.lowest = TRUE))
  observed_frequency <- as.numeric(observed_frequency)
  
  # Calculate expected frequencies
  breaks <- seq(min(data), max(data), length.out = length(observed_frequency) + 1)
  expected <- diff(pexp(breaks, rate = lambda_hat)) * length(data)
  
  # Chi-square statistic
  chi_square <- sum((observed_frequency - expected)^2 / expected)
  
  # Degrees of freedom = number of bins - 1 - number of estimated parameters (1 for lambda)
  df <- length(observed_frequency) - 1 - 1
  
  p_value <- pchisq(chi_square, df = df, lower.tail = FALSE)
  
  return(p_value)
}

