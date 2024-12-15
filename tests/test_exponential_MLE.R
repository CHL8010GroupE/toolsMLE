library(testthat)

source("R/Exponential.R")

test_that("test expl (Likelihood Function)", {
  # Test 1: Positive values
  lambda <- 2
  data <- c(1, 2, 3)
  expected <- prod(lambda * exp(-lambda * data))
  expect_equal(expl(lambda, data), expected)
  
  # Test 2: Negative values in data
  expect_error(expl(lambda, c(1, -2, 3)), "There is at least one negative value in the dataset.")
  
  # Test 3: Invalid lambda
  expect_error(expl(-1, data), "Lambda must be greater than 0.")
})

test_that("test expll (Log-Likelihood Function)", {
  # Test 1: Positive values
  lambda <- 2
  data <- c(1, 2, 3)
  n <- length(data)
  sum_x <- sum(data)
  expected <- n * log(lambda) - lambda * sum_x
  expect_equal(expll(lambda, data), expected)
  
  # Test 2: Negative values in data
  expect_error(expll(lambda, c(1, -2, 3)), "There is at least one negative value in the dataset.")
  
  # Test 3: Invalid lambda
  expect_error(expll(-1, data), "Lambda must be greater than 0.")
})

test_that("test expMLE (Find MLE)", {
  # Test 1: Regular data
  data <- c(1.2, 2.4, 0.8, 1.5, 3.6)
  result <- expMLE(data)
  expected_lambda <- 1 / mean(data)
  expect_equal(result$maximum, expected_lambda, tolerance = 0.0001)
  
  # Test 2: Negative values in data
  expect_error(expMLE(c(1.2, -0.8, 3.6)), "There is at least one negative value in the dataset.")
})

test_that("test expS (Score Function)", {
  # Test 1: Regular data
  lambda <- 2
  data <- c(1.2, 2.4, 0.8, 1.5, 3.6)
  n <- length(data)
  sum_x <- sum(data)
  expected <- n / lambda - sum_x
  expect_equal(expS(lambda, data), expected)
  
  # Test 2: Negative values in data
  expect_error(expS(lambda, c(1.2, -0.8, 3.6)), "There is at least one negative value in the dataset.")
  
  # Test 3: Invalid lambda
  expect_error(expS(-1, data), "Lambda must be greater than 0.")
})

test_that("test expLI (Likelihood Interval)", {
  # Test 1: Regular data
  data <- c(1.2, 2.4, 0.8, 1.5, 3.6)
  lambdahat <- 1 / mean(data)
  p <- 0.95
  interval <- expLI(lambdahat, data, p)
  expect_length(interval, 2)
  expect_gt(interval[1], 0)
  expect_gt(interval[2], interval[1])
  
  # Test 2: Negative values in data
  expect_error(expLI(lambdahat, c(1.2, -0.8, 3.6), p), "There is at least one negative value in the dataset.")
})

test_that("test expD (Likelihood Ratio Test)", {
  # Test 1: Regular data
  data <- c(1.2, 2.4, 0.8, 1.5, 3.6)
  lambdahat <- 1 / mean(data)
  lambda <- 2
  result <- expD(lambda, lambdahat, data)
  expect_named(result, c("D", "p_value"))
  expect_gt(result$D, 0)
  expect_gte(result$p_value, 0)
  expect_lte(result$p_value, 1)
  
  # Test 2: Invalid lambda
  expect_error(expD(-1, lambdahat, data), "Lambda must be greater than 0.")
})
