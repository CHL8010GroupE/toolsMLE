# toolsMLE

### Team Members (Alphabetize by Name)
- Bruce Liu
- Hanlong Chen
- Sze Wei Wang
- Yanyao Gu

## Overview
The toolsMLE package provides a suite of functions for performing Maximum Likelihood Estimation (MLE) and related statistical analyses for commonly used probability distributions. The package is designed to be user-friendly and flexible, offering interactive prompts for parameter inputs and streamlined workflows for parameter estimation and hypothesis testing.

### Course
CHL8010: Statistical Programming and Computation for Health Data    

### Instructor
Dr. Aya Mitani  

## Installation
You can install the development version from GitHub with:
```r
install.packages("devtools")
library(devtools)
devtools::install_github("CHL8010GroupE/toolsMLE")
```

## Features
- **`toolsMLE.l`**: Calculates the likelihood value for a given dataset and specified distribution.  
- **`toolsMLE.ll`**: Computes the log-likelihood value, offering enhanced numerical stability and interpretability.  
Supported distributions:
- Binomial
- Exponential
- Normal
- Poisson

### Maximum Likelihood Estimation (MLE)
- **`tools.MLE`**: Estimates the parameters of a specified distribution using MLE.  
Supports parameter estimation for:
  - Binomial: Probability of success
- Exponential: Rate parameter (λ)
- Normal: Mean (μ) and variance (σ²)
- Poisson: Rate parameter (λ)

### Likelihood Ratio Test (LRT)
- **`toolsMLE.lrt`**: Conducts a likelihood ratio test, comparing a null hypothesis parameter value with its MLE counterpart.  
Returns:
- Test statistic
- p-value  
Supports hypothesis testing for all included distributions.


## Example

data <- c(1.2, 2.4, 3.6)

Compute Likelihood
toolsMLE.l(data, dist = "exp") # For Exponential Distribution

Compute Log-Likelihood
toolsMLE.ll(data, dist = "norm") # For Normal Distribution

Find Maximum Likelihood Estimates
tools.MLE(data, dist = "pois") # For Poisson Distribution

Conduct a Likelihood Ratio Test
toolsMLE.lrt(data, dist = "bin") # For Binomial Distribution





## License
This project is licensed under the MIT License.
