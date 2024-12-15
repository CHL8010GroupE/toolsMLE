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
