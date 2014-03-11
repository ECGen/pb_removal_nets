#!/usr/bin/env sh
#$ -N pbr09
#$ -pe mpi 1
#$ -M mkl48@nau.edu
#$ -m as
#Wall time:
#$ -l h_rt=0:01:00

###Submission loop:
###for case in $ql; do qsub -o ~/qmonitor -e ~/qmonitor -v case=$case $pbrsrc/qStats_pbr09c.sh; sleep 1; done
###Note: sleep can be used to delay submissions

echo "starting case number $case"
cd $case
cp /home/mlau/projects/cg_simulations/src/cgsMods.R ./
cp /home/mlau/projects/cg_simulations/src/cgsNest.R ./
input=$(ls)
echo $input
Rscript cgsMods.R $input ~/projects/pb_removal_nets/results/null_mods09c.txt
Rscript cgsNest.R $input ~/projects/pb_removal_nets/results/null_nest09c.txt
rm ./cgsMods.R
rm ./cgsNest.R
echo "completed case number $case"
