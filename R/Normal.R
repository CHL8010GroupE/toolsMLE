set.seed(123)
x <- rnorm(100, mean = 5, sd = 2)
n <- length(x)

#Likelihood
norml <- function(mu, sigma2, x) {
  if (sigma2 <= 0) {
    stop("Variance (sigma2) must be greater than 0.")
  }
  n <- length(x)
  prod((1 / sqrt(2 * pi * sigma2)) * exp(-((x - mu)^2) / (2 * sigma2)))
}

#Log-likelihood for normal distribution.
normll <- function(mu, sigma2, x) {
  n <- length(x)
  (-n/2) * log(sigma2) - (1/(2 * sigma2)) * sum((x - mu)^2)
}

#Maximum likelihood function to find mu and sigma for normal distribution.
mle_normal <- function(x) {
  n <- length(x)
  mu_hat <- mean(x)
  sigma_hat <- sqrt(sum((x - mu_hat)^2) / n)

  return(list(mean = mu_hat, sigma = sigma_hat, sigma_square = sigma_hat^2))
}


# Likelihood Ratio Function
likelihood_ratio <- function(mu, sigma2, mu_hat, sigma2_hat, x) {
  ll_mu_sigma2 <- normll(mu, sigma2, x)  # Log-likelihood for given mu and sigma2
  ll_mu_hat_sigma2_hat <- normll(mu_hat, sigma2_hat, x)  # Log-likelihood for MLE estimates

  ratio <- exp(ll_mu_sigma2 - ll_mu_hat_sigma2_hat)  # Exponentiate the difference to get the ratio
  return(ratio)
}

#Confidence interval
confidence_interval <- function(x, alpha = 0.05) {
  n <- length(x)
  mu_hat <- mean(x)
  sigma_hat <- sqrt(sum((x - mu_hat)^2) / n)
  z <- qnorm(1 - alpha / 2)

  mu_ci <- c(mu_hat - z * sigma_hat / sqrt(n), mu_hat + z * sigma_hat / sqrt(n))
  sigma_ci <- c(sigma_hat * sqrt((n - 1) / qchisq(1 - alpha / 2, df = n - 1)),
                sigma_hat * sqrt((n - 1) / qchisq(alpha / 2, df = n - 1)))

  return(list(mean_ci = mu_ci, sigma_ci = sigma_ci))
}



