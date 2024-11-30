library(testthat)

source("R/Binomial.R")

test_that("test MLE.binol", {

  # Test One, test if it is the same
  data <- c(3,3,3,3,3,3,3)
  trials <- c(5,5,5,5,5,5,5)

  p <- 0.5
  estimate <- prod(p^(data) * (1 - p)^(trials - data) * factorial(trials) / (factorial(data) * factorial(trials - data)))
  expect_equal(MLE.binol(p, data, trials), estimate)

  # Test Two, test extreme value
  data <- c(0,0,0,0,0,0,0)
  trials <- c(5,5,5,5,5,5,5)
  p <- 0
  estimate <- prod(p^(data) * (1 - p)^(trials - data) * factorial(trials) / (factorial(data) * factorial(trials - data)))
  expect_equal(MLE.binol(p, data, trials), estimate)

  # Test Three, test bug
  expect_error(MLE.binol(-0.1, data, trials), "p must be within \\(0, 1\\)")

})

test_that("test MLE.binomial", {

  # Test 1, with fixed trial
  data <- c(3,3,3,3,3,3,3)
  trials <- 5

  result <- MLE.binomial(data, trials)

     # Since it is approximation, result$maximum may not strictly equal to 0.6
  expect_gt(result$maximum, 0.599)
  expect_lt(result$maximum, 0.611)


  # Test 2, with list of trials
  data <- c(3,6,9,12,15)
  trials <- c(5,10,15,20,25)

  result <- MLE.binomial(data, trials)
  expect_gt(result$maximum, 0.599)
  expect_lt(result$maximum, 0.611)

  # Test 3, with a fixed trial and a random data
  data <- c(8,3,7,3,2,6,7,2,7,1,6)
  trials <- 10

  expected <- sum(data) / (length(data) * trials)


  result <- MLE.binomial(data, trials)
     # allow a small deviance
  expect_equal(result$maximum, expected, tolerance = 0.0001)

  # Test 4, test with extreme case, probability = 1
  data <- c(10,10,10,10,10,10,10,10,10)
  trials <- 10

  expected <- sum(data) / (length(data) * trials)


  result <- MLE.binomial(data, trials)
  expect_equal(result$maximum, expected, tolerance = 0.0001)

  # Test 5, test with extreme case, probability = 0
  data <- c(0,0,0,0,0,0,0,0,0,0)
  trials <- 10

  expected <- sum(data) / (length(data) * trials)


  result <- MLE.binomial(data, trials)
  expect_equal(result$maximum, expected, tolerance = 0.0001)


  # Test 5, bug, with trial smaller than data
  data <- c(5,5,5)
  trials <- c(6,3,10)
  expect_error(MLE.binomial(data, trials), "trials must be greater or equal than data")
})

test_that("test MLE.binomialLI", {

  # Test 1: test if interval includes the true probability at (1 - alpha)% confidence level
  p <- 0.8
  n <- 30
  size <- 20
  trials <- rep(size, n)

  alpha <- 0.05

  count <- 0
  iterations <- 1000
  for (i in 1:iterations) {
    data <- rbinom(n, size, p)
    interval <- MLE.binomialLI(data, trials, alpha)
    if (p >= interval["lower"] && p <= interval["upper"]) {
      count <- count + 1
    }
  }
  proportion <- count / iterations

  print(proportion)

  print(proportion)
  expect_gt(proportion, 0.93)
  expect_lt(proportion, 0.97)


})

test_that("test MLE.binomialGOF", {

  n <- 100
  size <- 100
  # Test 1: test with a true binomial distributed data
  p <- 0.8
  trials <- rep(size, n)
  data <- rbinom(n, size, p)
  pvalue <- MLE.binomialGOF(data, trials)
  expect_gte(pvalue, 0.5)
  expect_lte(pvalue, 1)


  # Test 2: test with a randomness distributed data
  random_data <- sample(0:size, n, replace = TRUE)
  pvalue_random <- MLE.binomialGOF(random_data, trials)
  expect_gte(pvalue_random, 0)
  expect_lte(pvalue_random, 0.0001)

  # Test 3: if success = trial in each observation
  data <- rep(size, n)
  trials <- rep(size, n)
  pvalue <- MLE.binomialGOF(data, trials)
  expect_equal(pvalue, 1, tolerance = 0.0001)

  # Test 4: if success very close to trial in each observation,
  # Ideally, this should not follow a binomial distribution, because even p = 0.99, it should have high chance of 100 and 98 in the trial
  data <- rep(99, n)
  trials <- rep(size, n)
  pvalue <- MLE.binomialGOF(data, trials)
  expect_gte(pvalue_random, 0)
  expect_lte(pvalue_random, 0.001)

  # Test 5: if a data follow a Poisson Distribution
  lambda <- 80
  poisson_data <- rpois(n, lambda)
  trials <- rep(size, n)
  pvalue_poisson <- MLE.binomialGOF(poisson_data, trials)
  expect_equal(pvalue_poisson, 0, tolerance = 0.0001)
})
