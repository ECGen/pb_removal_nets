#! /bin/bash

#PBS -k o 
#PBS -l nodes=1:ppn=4,vmem=100gb,walltime=07:00:00
#PBS -M matthewklau@fas.harvard.edu
#PBS -m abe 
#PBS -N cm1
#PBS -j oe 
#PBS -d /N/u/mklau/Mason/mklau/pb_removal_nets/src

Rscript pb_cm.sh
