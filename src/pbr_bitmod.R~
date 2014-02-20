###Testing for the modularity of the network
##19 Feb 2014

library(bipartite)
library(methods)
source('../src/pbr_load_data.R')
m <- sign(pbr.09$c)[,apply(pbr.09$c,2,sum)!=0]
print('Null Modules')
null <- commsimulator(m,method='r1',thin=100)
mod.null <- slot(computeModules(null),name='likelihood')
cat(mod.null,file='../results/bitmodout.txt',append=TRUE,sep=' ')
