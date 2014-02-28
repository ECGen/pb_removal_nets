source('../src/pbr_load_data.R')
prn.08 <- list(c=pbr.08$c,x=pbr.08$x,cnpb=pbr.08$c[,colnames(pbr.08$c)!='pb'],xnpb=pbr.08$x[,colnames(pbr.08$x)!='pb'])
prn.09 <- list(c=pbr.09$c,x=pbr.09$x,cnpb=pbr.09$c[,colnames(pbr.09$c)!='pb'],xnpb=pbr.09$x[,colnames(pbr.09$x)!='pb'])

#run nestedness
nest08c <- oecosimu(prn.08$c,nestfun=nestedness,method='r1')
nest08x <- oecosimu(prn.08$x,nestfun=nestedness,method='r1')
nest08cnpb <- oecosimu(prn.08$cnpb,nestfun=nestedness,method='r1')
nest08xnpb <- oecosimu(prn.08$xnpb,nestfun=nestedness,method='r1')
                                        #
nest09c <- oecosimu(prn.09$c,nestfun=nestedness,method='r1')
nest09x <- oecosimu(prn.09$x,nestfun=nestedness,method='r1')
nest09cnpb <- oecosimu(prn.09$cnpb,nestfun=nestedness,method='r1')
nest09xnpb <- oecosimu(prn.09$xnpb,nestfun=nestedness,method='r1')
