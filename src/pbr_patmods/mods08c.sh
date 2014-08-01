#! /usr/bin/env bash

for i  in $( ls ./data/null08c/ ); do
Rscript get_mod.R ./data/null08c/$i ./data/mods08c.txt;
echo 'Now running modularity computation:'
echo $i;
done

