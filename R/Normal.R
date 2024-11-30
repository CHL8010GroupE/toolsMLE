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

#optim find MLE
# Use optim to find MLE. It finds minimum by default so easiest is to use negative log-likelihood
normalllik <- function(mu , sigma, x) {
  n <- length(x)
  -1 * ((-n/2) * log(sigma) - (1/(2 * sigma)) * sum((x - mu)^2))
}

mle_normal2 <- function(mu, sigma, x)
  optim_results <- optim(c(50, 10), normalllik, x = x)
  return(print(print(optim_results)))



# Likelihood ratio surface calculation
RLL <- expand.grid(mu = seq(95, 110, length.out = 50), sigma2 = (seq(12, 20, length.out = 50))^2)
RLL$r <- sapply(1:nrow(RLL), function(i) normll(RLL$mu[i], RLL$sigma2[i], x) - normll(mhat, sig2hat, x))
RLL$R <- exp(RLL$r)

# Visualization of the likelihood ratio
library(lattice)
wireframe(R ~ mu + sigma2, data = RLL, drape = TRUE)

gg <- ggplot(RLL, aes(x = mu, y = sigma2, z = R)) +
  geom_contour(breaks = c(0.01, 0.1, 0.5)) +
  geom_point(data = data.frame(mu = mhat, sigma2 = sig2hat, R = exp(normll(mhat, sig2hat, x) - normll(mhat, sig2hat, x)))) +
  annotate("text", x = c(101.5, 101.5, 101.5), y = c(16, 17.4, 18.7)^2, label = c("R = 0.5", "R = 0.1", "R = 0.01")) +
  xlab(expression(mu)) + ylab(expression(sigma^2))
gg

# Function to compute contours directly
mle_normald <- function(phi, muhat, sigma2hat, x, p) {
  g <- function(d, phi, muhat, sigma2hat, x, p) {
    normll(muhat + d * cos(phi), sigma2hat + d * sin(phi), x) - normll(muhat, sigma2hat, x) - log(p)
  }
  uniroot(g, c(0, 50), phi = phi, muhat = muhat, sigma2hat = sigma2hat, x = x, p = p)$root
}

