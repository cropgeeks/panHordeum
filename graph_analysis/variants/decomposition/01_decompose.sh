#!/bin/bash

#SBATCH -o slurm-%x_%A_%a.out 
#SBATCH --mem=10g
#SBATCH -p medium
#SBATCH --array=1-7

chrom="chr"$SLURM_ARRAY_TASK_ID"H"

vcfDir=./results/variants/chunks/vcf
inputVCF=$vcfDir/$chrom.vcf
prefix=$chrom


echo "count raw variants for chromosome $chrom"
date
numVars_inputVCF=`grep -v "^#" $inputVCF | wc -l `
echo "variant count in input = $numVars_inputVCF"

#remove variant sites with alleles larger than 100kb 
echo "run vcfbub"
date
vcfbub \
--input $inputVCF \
-l 0 \
-a 100000 \
> $prefix.vcfbub.vcf

conda deactivate


echo "count vcfbub variants for chromosome $chrom"
date
numVars_postVcfbub=`grep -v "^#" $prefix.vcfbub.vcf | wc -l `
echo "variant count post vcfbub = $numVars_postVcfbub"


#vcfwave realigns alternate alleles against the reference allele for each variant site using the bidirectional wavefront alignment (BiWFA) algorithm to decompose complex alleles into primitive ones
echo "run vcfwave"
date
source activate vcflib
vcfwave \
-I 1000 \
$prefix.vcfbub.vcf \
> $prefix.vcfbub.vcfwave.vcf


echo "count vcfwave variants for chromosome $chrom"
date
numVars_postvcfwave=`grep -v "^#" $prefix.vcfbub.vcfwave.vcf | wc -l `
echo "variant count post vcfwave = $numVars_postvcfwave"


echo "workflow complete"
date




