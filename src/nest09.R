library(bipartite)
source('../src/pbr_load_data_hoth.R')
prn.09 <- list(c=pbr.09$c,x=pbr.09$x,cnpb=pbr.09$c[,colnames(pbr.09$c)!='pb'],xnpb=pbr.09$x[,colnames(pbr.09$x)!='pb'])
                                        #run nestedness
nest09c <- oecosimu(prn.09$c,nestfun=nestedness,method='r1')
nest09x <- oecosimu(prn.09$x,nestfun=nestedness,method='r1')
nest09cnpb <- oecosimu(prn.09$cnpb,nestfun=nestedness,method='r1')
nest09xnpb <- oecosimu(prn.09$xnpb,nestfun=nestedness,method='r1')
                                        #write results
print('2009 Control')
dput(nest09c,'../results/nest09c.results')
print('2009 Exclusion')
dput(nest09x,'../results/nest09x.results')
print('2009 Control no pb')
dput(nest09cnpb,'../results/nest09cnpb.results')
print('2009 Exclusion no pb')
dput(nest09xnpb,'../results/nest09xnpb.results')

