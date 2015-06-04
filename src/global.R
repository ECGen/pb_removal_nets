library(bipartite)

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
