# Main functions of toolsMLE

source("R/Binomial.R")
source("R/Exponential.R")
source("R/Normal.R")
source("R/Poisson.R")

#' @title getting likelihood value
#' @param data the dataset to be used
#' @param dist type of distribution of the dataset
#' @return likelihood value
#' @export
toolsMLE.l <- function(data, dist=c("bin","exp","norm","pois")){
  if(dist=="pois"){
    para <- as.numeric(readline("Please enter the lambda: "))
    print(poisl(para,data))
  }
  if(dist=="bin"){
    para1 <- as.numeric(readline("Please enter the probability: "))
    para2 <- as.numeric(readline("Please enter the number of trial: "))
    print(MLE.binol(para1,data,para2))
  }
  if(dist=="exp"){
    para <- as.numeric(readline("Please enter the lambda: "))
    print(expl(para,data))
  }
  if(dist=="norm"){
    para1 <- as.numeric(readline("Please enter the mean: "))
    para2 <- as.numeric(readline("Please enter the sigma squared: "))
    print(norml(para1,para2,data))
  }
}

#' @title getting log-likelihood value
#' @param data the dataset to be used
#' @param dist type of distribution of the dataset
#' @return log-likelihood value
#' @export
toolsMLE.ll <- function(data, dist=c("bin","exp","norm","pois")){
  if(dist=="pois"){
    para <- as.numeric(readline("Please enter the lambda: "))
    print(poisll(para,data))
  }
  if(dist=="bin"){
    para1 <- as.numeric(readline("Please enter the probability: "))
    para2 <- as.numeric(readline("Please enter the number of trial: "))
    print(MLE.binoll(para1,data,para2))
  }
  if(dist=="exp"){
    para <- as.numeric(readline("Please enter the lambda: "))
    print(expll(para,data))
  }
  if(dist=="norm"){
    para1 <- as.numeric(readline("Please enter the mean: "))
    para2 <- as.numeric(readline("Please enter the sigma squared: "))
    print(normll(para1,para2,data))
  }
}

#' @title getting Maximum Likelihood Estimator
#' @param data the dataset to be used
#' @param dist type of distribution of the dataset
#' @return Maximum Likelihood Estimator
#' @export
tools.MLE <- function(data, dist=c("bin","exp","norm","pois")){
  if(dist=="pois"){
    print(poisMLE(data))
  }
  if(dist=="bin"){
    para <- as.numeric(readline("Please enter the number of trial: "))
    print(MLE.binomial(data,para))
  }
  if(dist=="exp"){
    print(expMLE(data))
  }
  if(dist=="norm"){
    print(mle_normal(data))
  }
}

#' @title conducting likelihood ratio test
#' @param data the dataset to be used
#' @param dist type of distribution of the dataset
#' @return likelihood ratio test statistic and p-value
#' @export
toolsMLE.lrt <- function(data, dist=c("bin","exp","norm","pois")){
  if(dist=="pois"){
    lambda0 <- as.numeric(readline("Please enter the lambda0: "))
    lambdahat <- as.numeric(readline("Please enter the lambdahat: "))
    print(poisD(lambda0,lambdahat,data))
  }
  if(dist=="bin"){
    trials <- as.numeric(readline("Please enter the number of trial: "))
    p0 <- as.numeric(readline("Please enter the p0: "))
    phat <- as.numeric(readline("Please enter the phat: "))
    print(MLE.binomialD(data,trails,p0,phat))
  }
  if(dist=="exp"){
    lambda0 <- as.numeric(readline("Please enter the lambda0: "))
    lambdahat <- as.numeric(readline("Please enter the lambdahat: "))
    print(expD(lambda0,lambdahat,data))
  }
  if(dist=="norm"){
    case <- readline("mean for unknown mu but known sigma2, variance for unknown sigma2 but known mu, both for unknown mu and sigma2. Please enter: ")
    if(case=="mean"){
      mu0 <- as.numeric(readline("Please enter the mu0: "))
      muhat <- as.numeric(readline("Please enter the muhat: "))
      sigma2 <- as.numeric(readline("Please enter the variance: "))
      print(normD(data,mu_null=mu0,mu_hat=muhat,sigma2_known=sigma2,case=case))
    }
    if(case=="variance"){
      sigma20 <- as.numeric(readline("Please enter the sigma_square0: "))
      sigma2hat <- as.numeric(readline("Please enter the sigma_squarehat: "))
      mu <- as.numeric(readline("Please enter the mu: "))
      print(normD(data,sigma2_null=sigma20,sigma2_hat=sigma2hat,mu_known=mu,case=case))
    }
    if(case=="both"){
      mu0 <- as.numeric(readline("Please enter the mu0: "))
      muhat <- as.numeric(readline("Please enter the muhat: "))
      sigma20 <- as.numeric(readline("Please enter the sigma_square0: "))
      sigma2hat <- as.numeric(readline("Please enter the sigma_squarehat: "))
      print(normD(data,mu_null=mu0,mu_hat=muhat,sigma2_null=sigma20,sigma2_hat=sigma2hat,case=case))
    }
  }
}

toolsMLE.LI <- function(data,dist=c("bin","exp","norm","pois",p)){
  if(dist=="pois"){
    para <- as.numeric(readline("Please enter the lambdahat: "))
    print(poisLI(para,data,p))
  }
  if(dist=="bin"){
    para1 <- as.numeric(readline("Please enter the phat: "))
    para2 <- as.numeric(readline("Please enter the number of trial: "))
    print(MLE.binomialLI(para1,data,para2,p))
  }
  if(dist=="exp"){
    para <- as.numeric(readline("Please enter the lambdahat: "))
    print(expLI(para,data,p))
  }
  if(dist=="norm"){
    case <- readline("mean for unknown mu but known sigma2, variance for unknown sigma2 but known mu. Please enter: ")
    if(case=="mean"){
      muhat <- as.numeric(readline("Please enter the muhat: "))
      sigma2 <- as.numeric(readline("Please enter the variance: "))
      print(normLI(mu_hat = muhat, sigma2_hat = sigma2, data, p, case))
    }
    if(case=="variance"){
      mu <- as.numeric(readline("Please enter the mu: "))
      sigma2hat <- as.numeric(readline("Please enter the sigma2_hat: "))
      print(normLI(mu_hat = mu, sigma2_hat = sigma2hat, data, p, case))
    }
  }
}
