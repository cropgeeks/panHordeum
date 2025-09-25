#!/bin/bash

#SBATCH --mem=600g
#SBATCH -o slurm-%x_%A.out  
#SBATCH --cpus-per-task=1
#SBATCH -p himem


numPermutations=100

ogFile=PanHordeum.og

source activate odgi

echo "run odgi heaps"
date
odgi heaps \
--threads $SLURM_CPUS_PER_TASK \
-i $ogFile \
-S \
-n $numPermutations \
> PanHordeum.heaps.txt


echo "done"
date
