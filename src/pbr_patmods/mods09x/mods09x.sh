#! /usr/bin/env bash

for i  in $( ls ../data/null09x/ ); do
Rscript ../get_mod.R ../data/null09x/$i ./mods09x.txt;
echo 'Now running modularity computation:'
echo $i;
done

