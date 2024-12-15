# Main functions of toolsMLE

source("R/Binomial.R")
source("R/Exponential.R")
source("R/Normal.R")
source("R/Poisson.R")

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
    mu0 <- as.numeric(readline("Please enter the mu0: "))
    muhat <- as.numeric(readline("Please enter the muhat: "))
    sigma20 <- as.numeric(readline("Please enter the sigma_square0: "))
    sigma2hat <- as.numeric(readline("Please enter the sigma_squarehat: "))
    print(normll(data,mu0,sigma20,muhat,sigma2hat))
  }
}
