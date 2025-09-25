#!/bin/bash

#SBATCH -o slurm-%x_%A_%a.out 
#SBATCH --mem=1g
#SBATCH -p short
#SBATCH --array=1-7

referenceFASTA=Barley_MorexV3_pseudomolecules.fasta

chrom="chr"$SLURM_ARRAY_TASK_ID"H"
prefix=$chrom

vcfFile=$prefix.normalised.vcf
header_file=$prefix.newVCFHeader.txt
output_vcf=$prefix.annotated_variants.vcf.gz


echo "run bcftools norm"
date
#need to change the chromosome IDs from those used in the graph to the standard designation
sed 's/Hvulga_MorexV3#1#//g' ../$chrom.vcfbub.vcfwave.vcf \
| bcftools norm \
--fasta-ref $referenceFASTA \
--check-ref w \
--multiallelics -any \
--output-type v \
- \
> $vcfFile


echo "Create a new header file with SVLEN definition"
date
echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">' > $header_file


echo "bgzip the VCF file"
date
bgzip \
-c \
$vcfFile \
> $vcfFile.bgz


echo "tabix the VCF file"
date
tabix \
--csi \
$vcfFile.bgz


echo "Add the new header to the VCF file"
date
bcftools annotate \
--header-lines $header_file \
-O z \
-o $prefix.variants_with_header.vcf.gz \
$vcfFile.bgz


echo "Calculate and add SVLEN values"
date
zcat $prefix.variants_with_header.vcf.gz | \
awk 'BEGIN {OFS="\t"} /^#/ {print} !/^#/ {svlen=length($5)-length($4); $8=$8";SVLEN="((svlen < 0) ? -svlen : svlen); print}' | \
bgzip -c > $output_vcf


echo "index the output VCF file"
date
bcftools index \
$output_vcf


echo "subset for structural variants"
date
bcftools view \
-i 'INFO/SVLEN>50' \
$output_vcf \
-O v \
-o $prefix.structural_variants.vcf


echo "extract SV insertions"
date
bcftools filter \
--include 'strlen(REF)<strlen(ALT)' \
--output-type v \
--output $prefix.insertions.vcf \
$prefix.structural_variants.vcf 


echo "convert insertions to BED format"
date
bcftools query \
-f '%CHROM\t%POS0\t%END\n' \
$prefix.insertions.vcf \
> $prefix.insertions.bed


echo "extract SV deletions"
date
bcftools filter \
--include 'strlen(REF)>strlen(ALT)' \
--output-type v \
--output $prefix.deletions.vcf \
$prefix.structural_variants.vcf


echo "convert deletions to BED format"
date
bcftools query \
-f '%CHROM\t%POS0\t%END\n' \
$prefix.deletions.vcf \
> $prefix.deletions.bed


echo "subset for small variants"
date
bcftools view \
-i 'INFO/SVLEN<=50' \
$output_vcf \
-O v \
-o $prefix.small_variants.vcf

echo "convert deletions to BED format"
date
bcftools query \
-f '%CHROM\t%POS0\t%END\n' \
$prefix.small_variants.vcf \
> $prefix.small_variants.bed


echo "workflow complete"
date




