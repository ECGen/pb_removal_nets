library(bipartite)
source('../src/pbr_load_data.R')
prn.09 <- list(c=pbr.09$c,x=pbr.09$x,cnpb=pbr.09$c[,colnames(pbr.09$c)!='pb'],xnpb=pbr.09$x[,colnames(pbr.09$x)!='pb'])
                                        #run nestedness
nest09c <- oecosimu(prn.09$c,nestfun=nestedness,method='r1')
nest09x <- oecosimu(prn.09$x,nestfun=nestedness,method='r1')
nest09cnpb <- oecosimu(prn.09$cnpb,nestfun=nestedness,method='r1')
nest09xnpb <- oecosimu(prn.09$xnpb,nestfun=nestedness,method='r1')
                                        #write results
dput(nest09c,'../results/nest09c.results')
dput(nest09x,'../results/nest09x.results')
dput(nest09cnpb,'../results/nest09cnpb.results')
dput(nest09xnpb,'../results/nest09xnpb.results')

