##mkNull09.R = making a set of null communities using the 2009 pbr control trees

#data input
source('../src/pbr_load_data_hoth.R')
dir.create('../data/null09c')
com <- pbr.09$c
method <- 'r1'
burnin <- 500
n.null <- 5000
                                        #burn
burn <- commsimulator(com,method=method)
null <- list()
for (i in 1:burnin){burn <- commsimulator(burn,method=method)}
for (i in 1:n.null){
  dir.create(paste('../data/null09c/',i,sep=''))
  write.csv(commsimulator(burn,method=method),paste(paste('../data/null09c/',i,sep=''),'null.csv',sep='/'),row.names=FALSE)
}
