plan <- drake_plan(
### Load data
    data.l = load_data("data"),
### Data wrangling
    bp.l = lapply(data.l, get_bp),
### Modeling
    mod.l = lapply(bp.l, get_mods), 
### Analyses
### Analytical Checks
### Plots
    plot_bp(bp.l[["2008"]][["c"]], 
            mod.l[["2008"]][["c"]],
            file = "results/bp_2008.pdf", 
            width = 12, 6)
### Tables
### Tables and Figures for Manuscript
### Generate the manuscript
    ## update.manuscript = update_manuscript(
    ##     files = tables_figures, 
    ##     dir = "docs/lcn_manuscript", 
    ##     file.tex = "main.tex")
)
