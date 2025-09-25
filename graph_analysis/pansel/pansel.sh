#!/bin/bash

#SBATCH -o slurm-%x_%A_%a.out 
#SBATCH --mem=30g
#SBATCH -p long
#SBATCH --array=1-7


chrom="chr"$SLURM_ARRAY_TASK_ID"H"

gfaDir=./files/chunks/gfa
input=$gfaDir/$chrom.gfa
output=$chrom.pansel.txt
refPath="Hvulga_MorexV3#1#"$chrom

echo "run pansel on chromosome $chrom for path $refPath"
date

pansel \
-i $input \
-r $refPath \
> $output 2> $chrom.pansel.log

echo "done"
date

