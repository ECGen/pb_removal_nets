##mkNull09.R = making a set of null communities using the 2009 pbr control trees

#data input
library(vegan)
source('../src/pbr_load_data_hoth.R')
dir.create('../data/null09c')
com <- pbr.09$c
com[com!=0] <- 1
com <- com[,apply(com,2,sum)!=0]
method <- 'r1'
burnin <- 500
n.null <- 5000
                                        #burn
burn <- commsimulator(com,method=method)
null <- list()
for (i in 1:burnin){burn <- commsimulator(burn,method=method)}
for (i in 1:n.null){
  print(i)
  dir.create(paste('../data/null09c/',i,sep=''))
  write.csv(commsimulator(burn,method=method),paste(paste('../data/null09c',i,sep='')
                                 ,paste('null',i,sep=''),sep='/'),row.names=FALSE)
}
