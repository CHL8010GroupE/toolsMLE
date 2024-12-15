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


# Relative Log-likelihood Function
normrll <- function(mu, sigma2, mu_hat, sigma2_hat, x) {
  ll_mu_sigma2 <- normll(mu, sigma2, x)  # Log-likelihood for given mu and sigma2
  ll_mu_hat_sigma2_hat <- normll(mu_hat, sigma2_hat, x)  # Log-likelihood for MLE estimates

  return(ll_mu_sigma2 - ll_mu_hat_sigma2_hat)
}

# Likelihood Interval Function
normLI <- function(mu_hat, sigma2_hat, x, p, case = c("mean", "variance")){
  normLIc <- function(mu, sigma2, mu_hat, sigma2_hat, x, p){
    normrll(mu, sigma2, mu_hat, sigma2_hat, x)-log(p)
  }
  if(case=="mean"){
    lower <- uniroot(normLIc,c(10^(-8),mu_hat),sigma2=sigma2_hat,mu_hat=mu_hat,sigma2_hat=sigma2_hat,x=x,p=p)$root
    upper <- uniroot(normLIc,c(10^(-8),mu_hat),sigma2=sigma2_hat,mu_hat=mu_hat,sigma2_hat=sigma2_hat,x=x,p=p)$root
  }
  if(case=="variance"){
    lower <- uniroot(normLIc,c(10^(-8),sigma2_hat),mu=mu_hat,mu_hat=mu_hat,sigma2_hat=sigma2_hat,x=x,p=p)$root
    upper <- uniroot(normLIc,c(10^(-8),sigma2_hat),mu=mu_hat,mu_hat=mu_hat,sigma2_hat=sigma2_hat,x=x,p=p)$root
  }
  return(c(lower,upper))
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

normD <- function(x, mu_null=NULL, sigma2_null=NULL, mu_hat=NULL,
                  sigma2_hat=NULL,mu_known=NULL, sigma2_known=NULL,
                  case) {

  n <- length(x)


  if (case == "mean") {
    # Log-likelihood under H0: mu = mu_null, sigma^2 = sigma2_known
    ll_null <- normll(mu_null, sigma2_known, x)
    # Log-likelihood under H1: mu = mu_hat, sigma^2 = sigma2_known
    ll_alt <- normll(mu_hat, sigma2_known, x)
    #testing one parameter (mu)
    df <- 1

  } else if (case == "variance") {
    # Log-likelihood under H0: sigma^2 = sigma2_null, mu_known
    ll_null <- normll(mu_known, sigma2_null, x)
    # Log-likelihood under H1: sigma^2 = sigma2_hat, mu = mu_known
    ll_alt <- normll(mu_known, sigma2_hat, x)
    # testing one parameter (sigma^2)
    df <- 1

  } else if (case == "both") {
    # Log-likelihood under H0: sigma^2 = sigma2_null, mu = mu_null
    ll_null <- normll(mu_null, sigma2_null, x)
    # Log-likelihood under H1: sigma^2 = sigma2_hat, mu = mu_hat
    ll_alt <- normll(mu_hat, sigma2_hat, x)
    df <- 2
  }
  else {
    stop("Invalid case. Use 'mean' or 'variance'.")
  }

  # -2(l_0 - l_1)
  test_stat <- -2 * (ll_null - ll_alt)


  p_value <- pchisq(test_stat, df, lower.tail = F)


  result <- list(
    test_statistic = test_stat,
    p_value = p_value
  )

  return(result)
}


