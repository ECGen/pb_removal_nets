###Make a set of null communities using the 2009 control trees
##18 Feb 2014

library(bipartite)
library(methods)
source('../src/pbr_load_data.R')
m <- sign(pbr.09$c)[,apply(pbr.09$c,2,sum)!=0]
n <- 1000
print('Null Modules')
null <- list()
mod.null <- list()
for (i in 1:n){null[[i]] <- commsimulator(m,method='r1',thin=100);print(i)}
