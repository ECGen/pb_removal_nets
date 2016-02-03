#! /bin/bash

#PBS -k o 
#PBS -l nodes=1:ppn=10,vmem=100gb,walltime=30:00:00
#PBS -M matthewklau@fas.harvard.edu
#PBS -m abe 
#PBS -N cmXavg
#PBS -j oe 
#PBS -d /N/u/mklau/Mason/mklau/pb_removal_nets/src/conMod/cmXavg

Rscript conMods.R
