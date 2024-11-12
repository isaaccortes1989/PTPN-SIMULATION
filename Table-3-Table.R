#===============================================#
#===read the results of the simulation study====#
#===============================================#
rm(list=ls())

n.seq      <- c(150,300,600,1000)
sigma.seq  <- c(1,2,3)
lambda.seq <- c(1,2)
gamma.seq  <- c(0.75,2.5)
results    <- c()

for(b in 1:length(sigma.seq))
{
  for(c in 1:length(lambda.seq))
  {
    for (d in 1:length(gamma.seq)) 
      {
      res.n       <- c()
      for(a in 1:length(n.seq))
      {
        n         <- n.seq[a]
        sigma     <- sigma.seq[b]
        lambda    <- lambda.seq[c]
        gamma     <- gamma.seq[d]
        name      <- paste("casessigma",sigma,"lambda",lambda,"gamma",gamma,"n",n,".txt",sep="")
        base      <- read.table(name,h=F)
        real      <- c(sigma,lambda,gamma)
        sesgo.aux <- base[,1:3]-matrix(real,ncol=3,nrow=nrow(base),byrow=T)
        sesgo     <- apply(sesgo.aux,2, mean)
        se        <- apply(base[,4:6],2, mean)
        RMSE      <- sqrt(apply(sesgo.aux^2,2, mean))
        cp.aux    <- matrix(0,ncol=3,nrow=nrow(base))
        for(m in 1:nrow(base))
        {
          for(jj in 1:3)
          {
            lim.inf      <- sesgo.aux[m,jj]-1.96*base[m,jj+3]
            lim.sup      <- sesgo.aux[m,jj]+1.96*base[m,jj+3]
            cp.aux[m,jj] <- ifelse(lim.inf<0 & lim.sup>0,1,0)
          }
        }
        cp    <- apply(cp.aux, 2, mean)
        res.n <- cbind(res.n,cbind(sesgo,se,RMSE,cp))
      }
      results <- rbind(results, res.n)
    }
  }
}

colnames <- rep(c("bias","SE","RMSE","CP"),3)
rownames <- rep(c("sigma","lambda","gamma"),nrow(results)/3)
round(results,4)

