# MLE functions for Exponential Distribution

# Likelihood function
expl <- function(lambda, x) {
  if (lambda <= 0) {
    stop("Lambda must be greater than 0.")
  }
  if (any(x < 0)) {
    stop("There is at least one negative value in the dataset.")
  }
  prod(lambda * exp(-lambda * x))
}

# Log likelihood function
expll <- function(lambda, x) {
  if (lambda <= 0) {
    stop("Lambda must be greater than 0.")
  }
  if (any(x < 0)) {
    stop("There is at least one negative value in the dataset.")
  }
  n <- length(x)
  sum_x <- sum(x)
  n * log(lambda) - lambda * sum_x
}

# Function to find MLE
expMLE <- function(x) {
  optimize(expll, c(1e-10, 1e10), x = x, maximum = TRUE)
}

# Score function
expS <- function(lambda, x) {
  if (lambda <= 0) {
    stop("Lambda must be greater than 0.")
  }
  if (any(x < 0)) {
    stop("There is at least one negative value in the dataset.")
  }
  sum_x <- sum(x)
  n <- length(x)
  n / lambda - sum_x
}

#MLE
mle_exp <- function(x){
  result <- uniroot(expS,c(0.000001,max(x)),x=x)$root
  return(list(rate = result))
}

# Information function
expI <- function(lambda, x) {
  if (lambda <= 0) {
    stop("Lambda must be greater than 0.")
  }
  if (any(x < 0)) {
    stop("There is at least one negative value in the dataset.")
  }
  n <- length(x)
  n / lambda^2
}

# Relative log likelihood function
exprll <- function(lambda, lambdahat, x) {
  if (lambda <= 0) {
    stop("Lambda must be greater than 0.")
  }
  if (lambdahat <= 0) {
    stop("Lambdahat must be greater than 0.")
  }
  if (any(x < 0)) {
    stop("There is at least one negative value in the dataset.")
  }
  expll(lambda, x) - expll(lambdahat, x)
}

# Function for likelihood intervals
expLI <- function(lambdahat, x, p) {
  if (lambdahat <= 0) {
    stop("Lambdahat must be greater than 0.")
  }
  if (any(x < 0)) {
    stop("There is at least one negative value in the dataset.")
  }
  expLIc <- function(lambda, lambdahat, x, p) {
    exprll(lambda, lambdahat, x) - log(p)
  }
  lower <- uniroot(expLIc, c(1e-10, lambdahat), lambdahat = lambdahat, x = x, p = p)$root
  upper <- uniroot(expLIc, c(lambdahat, lambdahat * 1e10), lambdahat = lambdahat, x = x, p = p)$root
  return(c(lower, upper))
}

# Likelihood ratio statistic
expD <- function(lambda, lambdahat, x) {
  if (lambda <= 0) {
    stop("Lambda must be greater than 0.")
  }
  if (lambdahat <= 0) {
    stop("Lambdahat must be greater than 0.")
  }
  if (any(x < 0)) {
    stop("There is at least one negative value in the dataset.")
  }
  -2 * exprll(lambda, lambdahat, x)
}

# Confidence interval based on likelihood ratio statistic
expDci <- function(lambda, lambdahat, x, conf) {
  if (lambda <= 0) {
    stop("Lambda must be greater than 0.")
  }
  if (lambdahat <= 0) {
    stop("Lambdahat must be greater than 0.")
  }
  if (any(x < 0)) {
    stop("There is at least one negative value in the dataset.")
  }
  -2 * exprll(lambda, lambdahat, x) - qchisq(conf, df = 1)
}

# Example usage
# x <- c(1.2, 2.4, 0.8, 1.5, 3.6)

# Find MLE for lambda
# mle_result <- expMLE(x)
# cat("MLE of lambda:", mle_result$maximum, "\n")

# Likelihood interval
# likelihood_interval <- expLI(mle_result$maximum, x, p = 0.95)
# cat("Likelihood interval for lambda:", likelihood_interval, "\n")
