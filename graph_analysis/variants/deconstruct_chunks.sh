#!/bin/bash

#SBATCH -o slurm-%x_%A_%a.out 
#SBATCH --mem=70g
#SBATCH --partition=medium
#SBATCH --array=1-7

pgDir=./files/chunks

pathPrefix="Hvulga_MorexV3"

chrom="chr"$SLURM_ARRAY_TASK_ID"H"

echo "process chromosome $chrom"
echo "input file = $pgDir/$chrom.gfa"

echo "copy input data to $TMPDIR on host $HOSTNAME"
date
cp -v $pgDir/gfa/$chrom.gfa $TMPDIR
workingDir=`pwd`
cd $TMPDIR

source activate vg

echo "run vg deconstruct"
date
vg deconstruct \
--path-prefix $pathPrefix \
--path-traversals \
--all-snarls \
--verbose \
$chrom.gfa \
> $chrom.vcf


echo "compute VCF stats"
date
bcftools stats \
$chrom.vcf \
> $chrom.vcf.stats


echo "run vg paths"
date
vg paths \
-x $pgDir/pg/$chrom.pg \
-L \
> $chrom.paths.txt


echo "copy results back to main storage"
date
cp -v $chrom.vcf* $workingDir
cp -v $chrom.paths.txt $workingDir


echo "done"
date