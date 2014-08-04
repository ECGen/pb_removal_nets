#! /usr/bin/env bash

for i  in $( ls ../data/null08x/ ); do
Rscript ../get_mod.R ../data/null08x/$i ./mods08x.txt;
echo 'Now running modularity computation:'
echo $i;
done

