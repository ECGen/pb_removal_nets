#! /usr/bin/env bash

echo 2008
echo controls
echo '...' |mail -s 'hoth: running prn 2008 control' mkl48@nau.edu
Files=../data/null/r00/*
for f in $Files
do
    echo $f
    echo modularity
    Rscript getMods.R $f ../results/modsr00.txt
    echo nestedness
    Rscript getNest.R $f ../results/nestr00.txt
    echo SES
    Rscript getCscore.R $f ../results/csr00.txt    
done

echo r0
echo '...' |mail -s 'hoth: running r0' mkl48@nau.edu
Files=../data/null/r0/*
for f in $Files
do
    echo $f
    echo modularity
    Rscript getMods.R $f ../results/modsr0.txt
    echo nestedness
    Rscript getNest.R $f ../results/nestr0.txt
    echo SES
    Rscript getCscore.R $f ../results/csr0.txt    
done

echo c0
echo '...' |mail -s 'hoth: running c0' mkl48@nau.edu
Files=../data/null/c0/*
for f in $Files
do
    echo $f
    echo modularity
    Rscript getMods.R $f ../results/modsc0.txt
    echo nestedness
    Rscript getNest.R $f ../results/nestc0.txt
    echo SES
    Rscript getCscore.R $f ../results/csc0.txt    
done

echo r1
echo '...' |mail -s 'hoth: running r1' mkl48@nau.edu
Files=../data/null/r1/*
for f in $Files
do
    echo $f
    echo modularity
    Rscript getMods.R $f ../results/modsr1.txt
    echo nestedness
    Rscript getNest.R $f ../results/nestr1.txt
    echo SES
    Rscript getCscore.R $f ../results/csr1.txt    
done

#email notification
echo '...' |mail -s 'hoth: getStats.sh done.' mkl48@nau.edu