#! /usr/bin/env bash

for i  in $( ls ./data/null09c/ ); do
Rscript get_mod.R ./data/null09c/$i ./data/mods09c.txt;
echo 'Now running modularity computation:'
echo $i;
done

