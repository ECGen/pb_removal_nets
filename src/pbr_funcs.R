library(bipartite)
library(methods)

###my module function
mm <- function(x){slot(computeModules(x),'likelihood')}
###sort network by degree
snort <- function(x){x[order(apply(x,1,sum),decreasing=TRUE),order(apply(x,2,sum),decreasing=TRUE)]}
###Mean networks
munet <- function(x,t){apply(x,2,function(x,t) tapply(x,t,mean),t=t)}
###SD networks
munet <- function(x,t){apply(x,2,function(x,t) tapply(x,t,sd),t=t)}
###module significance, only uses presence absence data
sigMods <- function(x,n=500){
  obs <- slot(computeModules(sign(x)),name='likelihood')
  sim <- list()
  for (i in 1:n){
    sim[[i]] <- commsimulator(sign(x),thin=50,method='r1')
  }
  sim <- unlist(lapply(sim,function(x) slot(computeModules(sign(x)),name='likelihood')))
  out <- c(obs=obs,mu.sim=mean(sim),sd.sim=sd(sim),
           z=((obs-mean(sim))/sd(sim)),
           p=((length(sim)[sim<=obs])/n)
           )
  return(out)
}
###contribution to modularity analysis
cMods <- function(x,n=100,thin=50,method='r1'){
  out <- 0
  obs <- slot(computeModules(x),name='likelihood')
  for (i in 1:nrow(x)){
    sim <- 0
    for (nit in 1:n){
      y <- x
      y[i,] <- commsimulator(x,method=method,thin=thin)[i,]
      sim[nit] <- slot(computeModules(y),name='likelihood')
    }
    out[i] <- (obs - mean(sim))/sd(sim)
  }
  return(out)
}
