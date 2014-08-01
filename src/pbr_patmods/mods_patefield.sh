#! /usr/bin/env bash

#generate the null communities
Rscript null_gen.R

#run modularity
Rscript run_mods.R
