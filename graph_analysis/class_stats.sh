#!/bin/bash
set -euo pipefail

# Compute the class stats (core, shell, cloud) for PanHordeum pangenome
# usage: bash class_stats.sh <PanHordeum.gfa>

graph=$1

# Compute node level stats
echo "Node level class stats"
echo "Core (c>=24), Shell (2<=c<24), Cloud (1<=c<2) nodes"
odgi paths -D "." -p 1 --coverage-levels=1,2,24 -i "$graph" | cut -f 3 | sort | uniq -c

# Compute base level stats 
echo "Base level class stats"
echo -e "Core\tShell\tCloud bases"

odgi stats -S -a".,0" -i "$graph" > "$graph.stats"
size=$(cut -f 1 "$graph.stats" | sed -n '2p')
tail -n +3 "$graph.stats" | awk -v size="$size" '{sum+=$3} END {print $2,size-$2-sum, sum}'

rm -f "$graph.stats"
