plan <- drake_plan(
### Load data
    data = load_data("data"),
### Data wrangling
    bp = lapply(data, proc_bp),
### Modeling
### Analyses
### Analytical Checks
### Plots
    plotweb(bp[["2009"]][["c"]], method = "normal"),
    cgPlotweb(data[["2009"]][, -1:-2], data[["2009"]][, 1])
### Tables
### Tables and Figures for Manuscript
### Generate the manuscript
    ## update.manuscript = update_manuscript(
    ##     files = tables_figures, 
    ##     dir = "docs/lcn_manuscript", 
    ##     file.tex = "main.tex")
)
