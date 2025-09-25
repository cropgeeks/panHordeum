#!/bin/bash

#SBATCH -o slurm-%x_%A.out   
#SBATCH -p short 
#SBATCH --mem=100m

dataDir=./vg_decomposed/extractSVs

outFile_structuralVariants=allSVs.bed
outFile_smallVariants=allSmallVariants.bed

for i in {1..7}
do
	chrom="chr"$i"H"
	echo "process chromosome $chrom"
	#combine the SV BED files
	cat $dataDir/$chrom.deletions.bed >> tmp.SVs.bed
	cat $dataDir/$chrom.insertions.bed >> tmp.SVs.bed
	cat $dataDir/$chrom.small_variants.bed >> tmp.smallVariants.bed
done

echo "sorting and redundancy removal"
date
sort -k1,1 -k2,2n tmp.SVs.bed | uniq > $outFile_structuralVariants  
sort -k1,1 -k2,2n tmp.smallVariants.bed | uniq > $outFile_smallVariants

echo "done"
date
