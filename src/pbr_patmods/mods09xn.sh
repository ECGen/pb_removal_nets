#! /usr/bin/env bash

for i  in $( ls ./data/null09xnpb/ ); do
Rscript get_mod.R ./data/null09xnpb/$i ./data/mods09xn.txt;
echo 'Now running modularity computation:'
echo $i;
done

