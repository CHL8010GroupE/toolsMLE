# MLE for Binomial Distribution
# example of data: c(1,2,4,6,2,2,1,5,6,2,4,5,3), a list of vector that record success for each trail

# binomial likelihood function
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


# Likelihood Ratio Test
MLE.binomialD <- function(data, trails, p_test){
  if (p_test < 0 || p_test > 1){
    stop("p-test must be within (0, 1)")
  }
  # generate log likelihood for null, p_test
  log_test <- MLE.binoll(p_test, data, trails)

  # generate log likelihood for p_alter
  estimated_parameter <- MLE.binomial(data, trails)


  log_alter <- estimated_parameter$objective

  # generate the likelihood ratio
  lR_statistic <- 2 * (log_alter - log_test)

  # assume degree of freedom is 1
  p_value <- pchisq(lR_statistic, df = 1, lower.tail = FALSE)

  p_value

}


# Goodness of fit test: test whether data follows a binomial distribution, print p-value
MLE.binomialGOF <- function(data, trails){
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


