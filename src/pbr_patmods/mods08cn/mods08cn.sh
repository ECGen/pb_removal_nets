#! /usr/bin/env bash

for i  in $( ls ../data/null08cnpb/ ); do
Rscript ../get_mod.R ../data/null08cnpb/$i ./mods08cn.txt;
echo 'Now running modularity computation:'
echo $i;
done

