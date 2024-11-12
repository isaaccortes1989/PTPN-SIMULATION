rm(list=ls(all=TRUE))
library(actuar)
library(pracma)
library(tpn)

#=======================#
#===Quantile Function===#
#=======================#

rPTPN <- function(n,sigma,lambda,gamma){
 U    <- runif(n)
 return(sigma*(qnorm(pnorm(lambda)*(U^(1/gamma)-1)+1)+lambda))
}

#=======================#
#===Sample size n=======#
#=======================#

rPTPNn <- function(n,sigma,lambda,gamma){
 x     <- numeric(n)
 for(i in 1:n){
 x[i] <- rPTPN(1,sigma,lambda,gamma)
 }
 return(x)
}

#=============================#
#===Log-likelihood function===#
#=============================#

log.likelihood = function(theta,x)
{
  sigma  <- theta[1]
  lambda <- theta[2]
  gamma  <- theta[3]
  logf   <- log(gamma)+dnorm(x/sigma-lambda,log=T)+(gamma-1)*log(pnorm(x/sigma-lambda)+pnorm(lambda)-1)-gamma*pnorm(lambda,log=T)-log(sigma)
  return(-sum(logf))
}

#=============================#
#===SIMULATION================#
#=============================#

simulation.PTPN<- function(n.iterations,sample.size,sigma,lambda){
  sigma.hat  <- rep(0, n.iterations)
  lambda.hat <- rep(0, n.iterations)
  mean.sigma <- 0 ; mean.lambda <- 0
  sd.sigma   <- 0 ; sd.lambda   <- 0
  iteration  <- 0
  while(iteration < n.iterations)
  {
    
    w       <- rPTPNn(sample.size,sigma,lambda)
    sample  <- w 
    result  <- optim(par = c(sigma,lambda), fn = log.likelihood, gr = NULL, method = c("L-BFGS-B"), 
                     lower = c(0.0001,0.001), upper = c(Inf, Inf), hessian = TRUE, x = sample)
    if (result$convergence == 0)
    {
      iteration             <- iteration + 1
      sigma.hat[iteration]  <- result$par[1]
      lambda.hat[iteration] <- result$par[2]
    }
  }
  mean.sigma = mean(sigma.hat) ; mean.lambda = mean(lambda.hat)
  sd.sigma   = sd(sigma.hat)   ; sd.lambda = sd(lambda.hat)
  cat("sigma.hat = ", mean.sigma, "sd.sigma = ", sd.sigma, fill = TRUE)
  cat("lambda.hat = ", mean.lambda, "sd.lambda = ", sd.lambda, fill = TRUE)
}


n.seq      <- c(150,300,600,1000)
sigma.seq  <- c(1,2,3)
lambda.seq <- c(1,2)
gamma.seq  <- c(0.75,2.5)
replicates <- 1000

for(a in 1:length(n.seq))
{
  for(b in 1:length(sigma.seq))
  {
    for(c in 1:length(lambda.seq))
    {
      for(d in 1:length(gamma.seq))
      {
        n       <- n.seq[a]
        sigma   <- sigma.seq[b]
        lambda  <- lambda.seq[c]
        gamma   <- gamma.seq[d]
        results <-c()
        m       <- 1
        while(m<=replicates)
        {
          y   <- rPTPNn(n,sigma,lambda,gamma)
          aux <- try(optim(par = c(sigma,lambda,gamma), fn = log.likelihood,x = y),silent=TRUE)
          
          if(!grepl("Error",aux)[1])
          {
            para <- c(aux$par[1],aux$par[2],aux$par[3])
            if(aux$conv==0)
            {
              hes <- try(solve(hessian(log.likelihood, x0=para, x=y)),silent=TRUE)
              if(!grepl("Error",hes)[1])
              {
                if(all(!is.na(hes)))
                {
                  if(min(diag(hes))>0)
                  {
                    
                    results <- rbind(results,c(para,sqrt(diag(hes))))
                    cc      <- paste(c("iteration",m,"sigma",sigma,"lambda",lambda,"gamma",gamma,"n",n),sep="")
                    print(cc)
                    m       <- m+1
                  }
                }
              }
            }
          }
        }
        name = paste("casessigma",sigma,"lambda",lambda,"gamma",gamma,"n",n,".txt",sep="")
        write(t(results), ncolumns = ncol(results),file=name)
      }
    }	
  }
}

