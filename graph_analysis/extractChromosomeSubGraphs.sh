#!/bin/bash

#SBATCH --mem=35g
#SBATCH -o slurm-%x_%A_%a.out 
#SBATCH -p short
#SBATCH --array=1-7

inputGraph=PanHordeum.pg

chrom="chr"$SLURM_ARRAY_TASK_ID"H"

path="Hvulga_MorexV3#1#"$chrom

source activate vg

echo "extract subgraph for chromosome $chrom"
date
vg chunk \
-x $inputGraph \
-c 20 \
-p $path \
-O pg \
> $chrom.pg


echo "done"
date
