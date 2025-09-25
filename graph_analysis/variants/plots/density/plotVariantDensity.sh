#!/bin/bash

#SBATCH -o slurm-%x_%A.out 
#SBATCH --mem=10g
#SBATCH -p short


source activate R_env

for variantType in SV SmallVariant
do

	echo "analyse variant type $variantType"

	echo "plot the data"
	Rscript \
	plotVariantDensity.r \
	"all"$variantType"s.bed" \
	$variantType \
	variantDensity.$variantType.png

done

echo "done"
