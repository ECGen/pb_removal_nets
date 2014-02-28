library(bipartite)
source('../src/pbr_load_data_hoth.R')
prn.09 <- list(c=pbr.09$c,x=pbr.09$x,cnpb=pbr.09$c[,colnames(pbr.09$c)!='pb'],xnpb=pbr.09$x[,colnames(pbr.09$x)!='pb'])
                                        #run nestedness
                                        #write results
print('2009 Control')
nest09c <- oecosimu(prn.09$c,nestfun=nested,method='r1')
dput(nest09c,'../results/nest09c.results')
print('2009 Exclusion')
nest09x <- oecosimu(prn.09$x,nestfun=nested,method='r1')
dput(nest09x,'../results/nest09x.results')
print('2009 Control no pb')
nest09cnpb <- oecosimu(prn.09$cnpb,nestfun=nested,method='r1')
dput(nest09cnpb,'../results/nest09cnpb.results')
print('2009 Exclusion no pb')
nest09xnpb <- oecosimu(prn.09$xnpb,nestfun=nested,method='r1')
dput(nest09xnpb,'../results/nest09xnpb.results')





