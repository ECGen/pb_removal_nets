#! /usr/bin/env bash

for i  in $( ls ./data/null09cnpb/ ); do
Rscript get_mod.R ./data/null09cnpb/$i ./data/mods09cn.txt;
echo 'Now running modularity computation:'
echo $i;
done

