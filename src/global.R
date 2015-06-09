library(bipartite)

meanMat <- function(x,y){
    (equalMat(x,y)[[1]] + equalMat(x,y)[[2]]) / 2
}

equalMat <- function(x,y){
    x <- list(x,x)
    y <- list(y,y)
    x[[2]] <- matrix(0,nrow=nrow(x[[1]]),
                           ncol=length(
                               colnames(y[[1]])[colnames(y[[1]]) 
                                                %in% 
                                                colnames(x[[1]]) == FALSE])
                     )
    colnames(x[[2]]) <- colnames(y[[1]])[colnames(y[[1]]) 
                                         %in% 
                                         colnames(x[[1]]) == FALSE]
    x[[2]] <- cbind(x[[1]],x[[2]])

    y[[2]] <- matrix(0,nrow=nrow(y[[1]]),
                           ncol=length(
                               colnames(x[[1]])[colnames(x[[1]]) 
                                                %in% 
                                                colnames(y[[1]]) == FALSE])
                     )
    colnames(y[[2]]) <- colnames(x[[1]])[colnames(x[[1]]) 
                                         %in% 
                                         colnames(y[[1]]) == FALSE]
    y[[2]] <- cbind(y[[1]],y[[2]])
    return(list(x=x[[2]],y=y[[2]]))
}

mkNull <- function(x,n,null.dir){
    for (i in 1:n){
        write.table(nullmodel(x,method='r2d',N=1),
                    file=paste(null.dir,i,sep='/'),
                    row.names=FALSE,col.names=FALSE,
                    sep=',')
    }
    print("Done!")
}

## calculate modularity
getMods <- function(null){
    out <- lapply(null,computeModules)
    out <- unlist(lapply(out,function(x) 
        slot(x,'likelihood')))
    return(out)
}
    
## calculate contribution to modularity
getCM <- function(x,null){
    out <- list()
    for (i in 1:nrow(x)){
        m <- x
        mod <- numeric()
        for (j in 1:length(null)){
            m[,i] <- null[[j]][,i]
            mod[j] <- slot(computeModules(m),
                        "likelihood")
        }
        out[[i]] <- mod
    }
    return(out)
}
    
## calculate statistics
getStats <- function(obs,sim){
    mu.sim <- mean(sim)
    sd.sim <- sd(sim)
    z <- (obs - mu.sim) / sd.sim
    p <- length(sim[sim >= obs]) / length(sim)
    out <- c(obs,mu.sim,sd.sim,z,p)
    return(out)
}
