# MLE for Binomial Distribution
# example of data: c(1,2,4,6,2,2,1,5,6,2,4,5,3), a list of vector that record success for each trail

#binomial likelihood function
MLE.binol <- function(p, data, trials){
  if (p < 0 || p > 1) {
    stop("p must be within (0, 1)")
  }
  prod(p^(data) * (1 - p)^(trials - data) * factorial(trials) / (factorial(data) * factorial(trials - data)))
}



# binomial log likelihood function
MLE.binoll <- function(p, data, trials){
  if (p < 0 || p > 1) {
    stop("p must be within (0, 1)")
  }

  binomial_log <- lfactorial(trials) - lfactorial(data) - lfactorial(trials - data)
  sum(data * log(p) + (trials - data) * log(1 - p) + binomial_log)
}


# find the optimal MLE function
MLE.binomial <- function(data, trials) {
  optimize(MLE.binoll, c(1e-10, 1 - 1e-10), data = data, trials = trials, maximum = TRUE)
}

# Score Function
MLE.binomialS <-function(p, data, trials){
  sum(data / p - (trials - data) / (1 - p))
}

# Information function
MLE.binomialI <- function(p, data, trials){
  -sum((trials * p  / p^2) + ((trials) * (1 - p) / (1 - p)^2))
}

# Relative Likelihood Function
MLE.binomialrll <- function(p, phat, data, trials){
  MLE.binoll(p, data, trials) - MLE.binoll(phat, data, trials)

}

# Function for likelihood intervals
MLE.binomialLI <- function(data, trials, alpha){
  estimated <- MLE.binomial(data, trials)$maximum
  MLE.binomialLIc <- function(data, trials, p, phat, alpha){
    MLE.binomialrll(p, phat, data, trials) - log(1 - alpha)
  }
  lower <- uniroot(MLE.binomialLIc, interval = c(1e-10, estimated), phat = estimated, data = data, trials = trials, alpha = alpha)$root

  upper <- uniroot(MLE.binomialLIc, interval = c(estimated, 1 - 1e-10), phat = estimated, data = data, trials = trials, alpha = alpha)$root

  c(lower = lower, upper = upper)
}

# Likelihood Ratio Test
MLE.binomialD <- function(data, trials, p_test){
  if (p_test < 0 || p_test > 1){
    stop("p-test must be within (0, 1)")
  }
  # generate log likelihood for null, p_test
  log_test <- MLE.binoll(p_test, data, trials)

  # generate log likelihood for p_alter
  estimated_parameter <- MLE.binomial(data, trials)


  log_alter <- estimated_parameter$objective

  # generate the likelihood ratio
  lR_statistic <- 2 * (log_alter - log_test)

  # assume degree of freedom is 1
  p_value <- pchisq(lR_statistic, df = 1, lower.tail = FALSE)

  p_value

}


# confidence interval based on likelihood ratio statistic
MLE.binomialDci <- function(ptest, data, trials, conf){
  # To do
}

# Goodness of fit test: test whether data follows a binomial distribution, print p-value
MLE.binomialGOF <- function(data, trials){
  # Check for valid input
  if (any(data > trials)) {
    stop("Data values cannot exceed the number of trials.")
  }

  # estimate p
  loglikelihood <- MLE.binomial(data, trials)
  p_hat <- loglikelihood$maximum

  possible_event <-  0:trials # this represent all possible event

  observed_frequency <- as.numeric(table(factor(data, levels = possible_event)))

  expected <- dbinom(possible_event, size = trials, prob = p_hat) * length(data)

  chi_square <- sum((observed_frequency - expected)^2 / expected)

  p_value <- pchisq(chi_square, df = length(possible_event) - 1, lower.tail = FALSE)

  p_value

}


