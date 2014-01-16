###Run the QAP test for Art's data

###Running co-occurrence analyses on Art's data
source('../src/pbr_load_data.R')
library(sna)
nits <- 9999

###Get networks
net.08 <- lapply(pbr.08,CoNetwork,threshold=1)
net.09 <- lapply(pbr.09,CoNetwork,threshold=1)

###QAP
reduce.net <- function(x){
  if (length(x)!=2){warning('Error: matrix list not of length 2!')}
  x1 <- x[[1]]
  x2 <- x[[2]]
  if ((dim(x1)[1]!=dim(x1)[2])|(dim(x2)[1]!=dim(x2)[2])){warning('Matrices are not square!')}
  if (any(dim(x1)!=dim(x2))){warning('Dimensions are not equal!')}
  if (any(rownames(x1)!=rownames(x2))){warning('Node names are mismatched!')}
  out <- list()
  out[[1]] <- x1[(apply(x1,1,sum)!=0|apply(x2,1,sum)!=0),(apply(x1,2,sum)!=0|apply(x2,2,sum)!=0)]
  out[[2]] <- x2[(apply(x1,1,sum)!=0|apply(x2,1,sum)!=0),(apply(x1,2,sum)!=0|apply(x2,2,sum)!=0)]
  names(out) <- names(x)
  return(x)
}
net.dif <- function(dat,dat2=NULL,g1=NULL,g2=NULL){
  sum((dat[,,g1]-dat[,,g2])^2)
}

net.08 <- reduce.net(net.08)
net.09 <- reduce.net(net.09)
qin.08 <- array(NA,dim=c(nrow(net.08[[1]]),ncol(net.08[[1]]),2))
qin.08[,,1] <- net.08[[1]]
qin.08[,,2] <- net.08[[2]]
qin.09 <- array(NA,dim=c(nrow(net.09[[1]]),ncol(net.09[[1]]),2))
qin.09[,,1] <- net.09[[1]]
qin.09[,,2] <- net.09[[2]]
rownames(qin.08) <- colnames(qin.08) <- rownames(net.08[[1]])
rownames(qin.09) <- colnames(qin.09) <- rownames(net.09[[1]])
qap08 <- qaptest(qin.08,net.dif,g1=1,g2=2,reps=nits)
qap09 <- qaptest(qin.09,net.dif,g1=1,g2=2,reps=nits)

###Output
write.csv(qap08$dist,file='../results/qap_dist_08.csv',row.names=FALSE)
write.csv(qap09$dist,file='../results/qap_dist_09.csv',row.names=FALSE)
writeLines(capture.output(qap08),con='../results/qap_results_08.txt')
writeLines(capture.output(qap09),con='../results/qap_results_09.txt')
