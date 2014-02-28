library(bipartite)
source('../src/pbr_load_data_hoth.R')
prn.08 <- list(c=pbr.08$c,x=pbr.08$x,cnpb=pbr.08$c[,colnames(pbr.08$c)!='pb'],xnpb=pbr.08$x[,colnames(pbr.08$x)!='pb'])
                                        #run nestedness
nest08c <- oecosimu(prn.08$c,nestfun=nestedness,method='r1')
nest08x <- oecosimu(prn.08$x,nestfun=nestedness,method='r1')
nest08cnpb <- oecosimu(prn.08$cnpb,nestfun=nestedness,method='r1')
nest08xnpb <- oecosimu(prn.08$xnpb,nestfun=nestedness,method='r1')
                                        #write results
dput(nest08c,'../results/nest08c.results')
dput(nest08x,'../results/nest08x.results')
dput(nest08cnpb,'../results/nest08cnpb.results')
dput(nest08xnpb,'../results/nest08xnpb.results')
