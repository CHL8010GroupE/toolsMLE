# MLE functions for Poisson distribution

# likelihood function
poisl <- function(lambda, x){
  if(lambda<=0){stop("Lambda must be greater than 0.")}
  if(any(x<0)){stop("There is at least one negative value in the dataset.")}
  prod((lambda^x * exp(-lambda)) / factorial(x))
}

# log likelihood function
poisll <- function(lambda, x){
  if(lambda<=0){stop("Lambda must be greater than 0.")}
  if(any(x<0)){stop("There is at least one negative value in the dataset.")}
  n <- length(x)
  sum_x <- sum(x)
  sum_x * log(lambda) - n * lambda - sum(lfactorial(x))
}

# Score function
poisS <- function(lambda, x){
  if(lambda<=0){stop("Lambda must be greater than 0.")}
  if(any(x<0)){stop("There is at least one negative value in the dataset.")}
  sum_x <- sum(x)
  n <- length(x)
  sum_x / lambda - n
}

# function to find MLE
poisMLE <- function(x){
  lambda <- uniroot(poisS,c(0.000001,max(x)),x=x)$root
  return(lambda)
}

# Information function
poisI <- function(lambda, x){
  if(lambda<=0){stop("Lambda must be greater than 0.")}
  if(any(x<0)){stop("There is at least one negative value in the dataset.")}
  n <- length(x)
  n / lambda
}

# relative log likelihood function
poisrll <- function(lambda, lambdahat, x){
  if(lambda<=0){stop("Lambda must be greater than 0.")}
  if(lambdahat<=0){stop("Lambdahat must be greater than 0.")}
  if(any(x<0)){stop("There is at least one negative value in the dataset.")}
  poisll(lambda,x)-poisll(lambdahat,x)
}

# function for likelihood intervals
poisLI <- function(lambdahat, x, p){
  if(lambdahat<=0){stop("Lambdahat must be greater than 0.")}
  if(any(x<0)){stop("There is at least one negative value in the dataset.")}
  poisLIc <- function(lambda, lambdahat, x, p){
    poisrll(lambda,lambdahat,x)-log(p)
  }
  lower <- uniroot(poisLIc,c(10^(-8),lambdahat),lambdahat=lambdahat,x=x,p=p)$root
  upper <- uniroot(poisLIc,c(lambdahat,lambdahat*10^8),lambdahat=lambdahat,x=x,p=p)$root
  return(c(lower,upper))
}

# likelihood ratio statistic
poisD <- function(lambda, lambdahat, x){
  if(lambda<=0){stop("Lambda must be greater than 0.")}
  if(lambdahat<=0){stop("Lambdahat must be greater than 0.")}
  if(any(x<0)){stop("There is at least one negative value in the dataset.")}
  -2*poisrll(lambda,lambdahat,x)
}

# confidence interval based on likelihood ratio statistic
poisDci <- function(lambda, lambdahat, x, conf) {
  if(lambda<=0){stop("Lambda must be greater than 0.")}
  if(lambdahat<=0){stop("Lambdahat must be greater than 0.")}
  if(any(x<0)){stop("There is at least one negative value in the dataset.")}
  -2 * poisrll(lambda, lambdahat, x) -qchisq(conf, 1)
}
