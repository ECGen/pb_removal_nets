###Testing for the modularity of the network
##18 Feb 2014

library(bipartite)
source('../src/pbr_load_data.R')
m <- sign(pbr.09$c)[,apply(pbr.09$c,2,sum)!=0]
n <- 1000
                                        #generate randomized matrices
null <- list()
for (i in 1:n){null[[i]] <- commsimulator(m,method='r1',thin=100)}
mod.null <- lapply(null,function(x) slot(computeModules(x),name='likelihood'))
mod.null <- unlist(mod.null)
mod.obs <- computeModules(m)
z <- (mod.obs-mean(mod.null))/sd(mod.null)
p <- c(lower.p=length(mod.null[mod.null<=mod.obs]),upper.p=length(mod.null[mod.null>=mod.obs]))/n
out <- list(null=mod.null,results=c(obs=mod.obs,z=z,p))
dput(out,file='../results/pbr_modularity.R')

