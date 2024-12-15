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

# Likelihood Ratio Test Function

normD <- function(x, mu_null, sigma2_null, mu_hat, sigma_hat,mu_known, sigma_known, case = "mean") {
  n <- length(x)


  if (case == "mean") {
    if (is.null(mu_null) || is.null(sigma2_null)) {
      stop("Provide both mu_null and sigma2_null for testing mean.")
    }
    # Log-likelihood under H0: mu = mu_null, sigma^2 = sigma2_null
    ll_null <- normll(mu_null, sigma_known, x)
    # Log-likelihood under H1: mu = mu_hat, sigma^2 = sigma2_null
    ll_alt <- normll(mu_hat, sigma_known, x)
    df <- 1  # Degrees of freedom: testing one parameter (mu)
  } else if (case == "variance") {
    if (is.null(sigma2_null) || is.null(mu_null)) {
      stop("Provide both sigma2_null and mu_null for testing variance.")
    }
    # Log-likelihood under H0: sigma^2 = sigma2_null, mu = mu_null
    ll_null <- normll(mu_known, sigma2_null, x)
    # Log-likelihood under H1: sigma^2 = sigma2_hat, mu = mu_null
    ll_alt <- normll(mu_known, sigma2_hat, x)
    df <- 1  # Degrees of freedom: testing one parameter (sigma^2)
  } else if (case == "both") {
    if (is.null(sigma2_null) || is.null(mu_null)) {
      stop("Provide both sigma2_null and mu_null for testing variance.")
    }
    # Log-likelihood under H0: sigma^2 = sigma2_null, mu = mu_null
    ll_null <- normll(mu_null, sigma2_null, x)
    # Log-likelihood under H1: sigma^2 = sigma2_hat, mu = mu_hat
    ll_alt <- normll(mu_hat, sigma2_hat, x)
    df <- 2
  }
  else {
    stop("Invalid case. Use 'mean' or 'variance'.")
  }

  # Likelihood ratio test statistic
  test_stat <- -2 * (ll_null - ll_alt)

  # P-value and critical value
  p_value <- 1 - pchisq(test_stat, df)

  # Return results
  result <- list(
    test_statistic = test_stat,
    p_value = p_value
  )

  return(result)
}


