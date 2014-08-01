#! /usr/bin/env bash

for i  in $( ls ./data/null08xnpb/ ); do
Rscript get_mod.R ./data/null08xnpb/$i ./data/mods08xn.txt;
echo 'Now running modularity computation:'
echo $i;
done

