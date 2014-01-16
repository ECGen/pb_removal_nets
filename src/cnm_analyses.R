###Running co-occurrence analyses on Art's data
source('../src/pbr_load_data.R')
###Co-occurrence patterns
cnm.08c <- cnm.test(pbr.08[[1]],nits=1)
cnm.08x <- cnm.test(pbr.08[[2]])
cnm.09c <- cnm.test(pbr.09[[1]])
cnm.09x <- cnm.test(pbr.09[[2]])
write.csv(cbind(cnm.08c,cnm.08x,cnm.09c,cnm.09x),file='../results/cnm_results.csv',row.names=TRUE)
