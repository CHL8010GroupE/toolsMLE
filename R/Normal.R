set.seed(123)
x <- rnorm(100, mean = 5, sd = 2)
n <- length(x)
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




#plot
library(ggplot2)
visualize_gaussian_mle <- function(data) {
  est <- mle_normal(data)
  mu_hat <- est$mean
  sigma_hat <- est$sd

  ggplot(data.frame(x = data), aes(x = x)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.5, fill = "lightblue", color = "black") +
    stat_function(fun = dnorm, args = list(mean = mu_hat, sd = sigma_hat), color = "red", size = 1) +
    labs(title = "Gaussian MLE Fit", x = "Data", y = "Density")
}



