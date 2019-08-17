#!/bin/sh
module load bcftools/1.8
module load parallel
#This script will perform two analyses on the vcfs generated for Solenopsis invicta populations in South
#and North America. The first analysis will use bcftools isec to determine which variants are in common between
#the two populations.
#The second analysis will use SNPeff on the South American VCF to check for the effect of the SNPs on the
#genes of S.invicta

#Comparison between populations
#Compress the files for bcftools to work
ls input/*vcf | cut -d '/' -f 2 | parallel 'bgzip -c input/{} > tmp/{}.gz' :::
#Index files
ls tmp/subset_*gz | parallel 'bcftools index {}' :::

#Output positions present in both files, write the common positions in the South American vcf
bcftools isec -n=2 -w1 tmp/subset_south_america.vcf.gz tmp/subset_north_america.vcf.gz -o tmp/common_snps.vcf
#Number of variants in common between the two populations
wc -l tmp/common_snps.vcf
#2877, ~ 92% of the SNPs in South American

#Use SNPeff to detect the effect of SNPs in genes (only common SNPs).
#The VCF is based on gnG (from RefSeq), but SNPeff uses annotations based on the Ensembl notation.
#Use 'flo' and CrossMap to liftover the RefSeq VCF notation to Ensembl
sh generate_ensembl_vcf.sh

#Run SNPeff on the South American vcf (made verbose by -v)
java -jar snpEff.jar Solenopsis_invicta -v\
          tmp/ensembl_sa.vcf > tmp/eff_supergene.vcf
