## Check for supporting packages
cran.pkgs <- c("magrittr", "devtools", "rmarkdown", "beepr",
               "readxl", "xtable", "reshape", "lubridate",
               "ggplot2", "enaR", "sna", "igraph", "Rgraphviz", 
               "gplots", "enaR", "grid", "bipartite", "RCurl", 
               "XML", "vegan", "RColorBrewer")
gh.pkgs <- c("ROpenSci/drake", "r-lib/styler")
gh.pkgs.names <- do.call(rbind, 
                         strsplit(gh.pkgs, 
                                  split = "/"))[, 2]
## install packages that are not installed
## CRAN
if (any(!(cran.pkgs %in% installed.packages()[, 1]))){
    sapply(cran.pkgs[which(!(cran.pkgs %in% 
                             installed.packages()[, 1]))], 
           install.packages, 
           dependencies = TRUE, 
           repos = 'http://cran.us.r-project.org')
}
## github 
if (any(!(gh.pkgs.names %in% installed.packages()[, 1]))){
    sapply(gh.pkgs[which(!(gh.pkgs.names %in% 
                           installed.packages()[, 1]))], 
           devtools::install_github,
           dependencies = TRUE)
}
## Load libraries
sapply(c(cran.pkgs, gh.pkgs.names), 
       library, quietly = TRUE, character.only = TRUE)


