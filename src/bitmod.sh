#! /usr/bin/env bash

#clear the results location
rm ../results/bitmodout.txt

#loop the pbr_bitmod.R script
for i in {1..1000}
do
    echo $i
    Rscript pbr_bitmod.R
done

#email notification
echo 'Automated message' |mail -s 'hoth: bitmod.sh done.' mkl48@nau.edu