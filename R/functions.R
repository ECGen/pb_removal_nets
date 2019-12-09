plot_bp <- function(net, mod, file = "results/bp_plot.pdf", 
                    width = 12, height = 4, renumber = TRUE){
    net <- net[apply(net, 1, sum) != 0, apply(net, 2, sum) != 0]
    mod.hi <- mod[match(colnames(net), names(mod))]
    mod.lo <- mod[match(rownames(net), names(mod))]
    col.pal <- grey(seq(0.9,0.1, length = max(mod)))
    col.hi <- col.pal[mod.hi]
    col.lo <- col.pal[mod.lo]
    if (renumber){colnames(net) <- seq(1, ncol(net))}
    pdf(file = file, width = width, height = height)
    plotweb(net,
            method = "normal", 
            text.rot = 90, 
            col.high = col.hi,
            col.low = col.lo,
            col.interaction = "grey20",
            bor.col.interaction = "grey20",
            bor.col.high = "grey90",
            bor.col.low = "grey90")
    dev.off()
}

get_bp <- function(x, reorder = TRUE){
    trt <- do.call(rbind, strsplit(x[, "Position"], " "))[, 2]
    bp <- split(x, trt)
    for (i in seq_along(bp)){
        bp[[i]] <- tree2geno(bp[[i]][, -1:-2], bp[[i]][, 1])
    }
    if (reorder){
        bp <- lapply(bp, msort)
    }
    if (any(is.na(bp))){warning("NA values produced.")}
    return(bp)
}

tree2geno <- function(com, g, standardize = FALSE){
    if (standardize){com <- apply(com, 2, std.max)}
    out <- apply(com, 2, function(x, g) tapply(x, g, mean), g = g)
    if (any(is.na(out))){warning("NA values produced.")}
    return(out)
}

msort <- function(x, decreasing = TRUE){
    x[order(apply(x, 1, sum), decreasing = decreasing), 
      order(apply(x, 2, sum), decreasing = decreasing)]
}

std.max <- function(x){
    if (sum(x) == 0){x}else{x / max(x)}
}

get_mods <- function(bp){
    mod <- lapply(bp, computeModules, deep = FALSE)
    out <- lapply(mod, get_modules)
    return(out)
}

get_modules <- function(mod){
    m <- listModuleInformation(mod)
    m <- m[[-1]]
    out <- m
    for (i in seq_along(m)){
        for (j in seq_along(m[[i]])){
            out[[i]][[j]] <- rep(i, length(m[[i]][[j]]))
            names(out[[i]][[j]]) <- m[[i]][[j]]
        }
    }
    out <- unlist(out)
    return(out)
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
