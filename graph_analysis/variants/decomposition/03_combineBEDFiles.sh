#!/bin/bash

#SBATCH -o slurm-%x_%A.out   
#SBATCH -p short 
#SBATCH --mem=100m


echo "combine and format DEL BED files"
date
cat *.deletions.bed | sort -k1,1 -k2,2n | uniq > vg.DEL.bed

echo "combine and format INS BED files"
date
cat *.insertions.bed | sort -k1,1 -k2,2n | uniq > vg.INS.bed


echo "done"
date