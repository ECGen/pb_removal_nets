proc_bp <- function(x, reorder = TRUE){
    trt <- do.call(rbind, strsplit(x[, "Position"], " "))[, 2]
    bp <- split(x, trt)
    for (i in seq_along(bp)){
        bp[[i]] <- tree2geno(bp[[i]][, -1:-2], bp[[i]][, 1])
    }
    if (reorder){
        bp <- lapply(bp, msort)
    }
    return(bp)
}

tree2geno <- function(com, g, standardize = TRUE){
    if (standardize){com <- apply(com, 2, std.max)}
    apply(com, 2, function(x, g) tapply(x, g, mean), g = g)
}

msort <- function(x, decreasing = FALSE){
    x[order(apply(x, 1, sum)), order(apply(x, 2, sum))]
}

std.max <- function(x){
    x / max(x)
}

load_data <- function(dir, number.species = TRUE){
    files <- paste(dir, c("keith_pb_removal_2008.csv",
               "keith_pb_removal_2009.csv",
               "keith_pb_trt_effect_2008.csv",
               "keith_pb_trt_effect_2009.csv"), sep = "/")
    data <- lapply(files, read.csv)
    names(data) <- c("2008", "2009", "pb2008", "pb2009")
    ## Change empty/NA cells to zero
    for (i in seq_along(data)){data[[i]][is.na(data[[i]])] <- 0}
    ## Fix data entry errors
    data[["2008"]][, "Position"] <- fix_position(data[["2008"]][, "Position"])
    data[["2009"]][, "Position"] <- fix_position(data[["2009"]][, "Position"])
    ## Append PB data
    data[["2008"]] <- append_pb(data[["2008"]], data[["pb2008"]])
    data[["2009"]] <- append_pb(data[["2009"]], data[["pb2009"]])
    out <- data[c("2008", "2009")]
    if (number.species){
        for (i in seq_along(out)){
            colnames(out[[i]])[seq(3, ncol(out[[i]]))] <- seq(3, ncol(out[[i]]))
        }
    }
    return(out)
}

fix_position <- function(x){
    out <- gsub("x ", "x", x)
    out <- gsub("c ", "c", out)
    out <- gsub("x", " x", out)
    out <- gsub("c", " c", out)
    out <- gsub("  ", " ", out)
    return(out)
}

append_pb <- function(com, pb){
    com[, "Position"] <- as.character(com[, "Position"])
    colnames(pb) <- tolower(colnames(pb))
    if (all(c("x", "c") %in%  colnames(pb))){
        pbx <- pb[, c("tree", "x")]
        pbc <- pb[, c("tree", "c")]
        pbx[, "tree"] <- paste(pbx[, "tree"], "x")
        pbc[, "tree"] <- paste(pbc[, "tree"], "c")
    } else {
        pbx <- pb[pb[, "trt"] == 2, c("tree", "pb")]
        pbc <- pb[pb[, "trt"] == 1, c("tree", "pb")]
        pbx[, "tree"] <- paste(pbx[, "tree"], "x")
        pbc[, "tree"] <- paste(pbc[, "tree"], "c")
    }
    colnames(pbx) <- colnames(pbc) <- c("Position", "pb")
    pbxc <- rbind(pbx, pbc)
    pbxc <- pbxc[match(com[, "Position"], pbxc[, "Position"]), ]
    if(!all(com[, "Position"] == pbxc[, "Position"])){
        warning("Positions do not match")
    }
    com <- cbind(com, pb = pbxc[, "pb"])        
    return(com)
}
