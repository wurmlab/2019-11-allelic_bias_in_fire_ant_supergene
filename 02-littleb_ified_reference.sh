#!/bin/sh
module load bcftools/1.8

#Commands to generate a vcf file with fix differences beteween B and b only
gunzip -c input/Si_gnG.fa.gz > tmp/Si_gnG.fa
ln -sfr input/subset.vcf tmp/
bgzip -f tmp/subset.vcf
tabix -fp vcf tmp/subset.vcf.gz

#The fixed difference between bb and lb have already been selected in the input vcf file.7
#No need therefore to select for fixed differences again in the input vcf.

#Build a consensus with the fix differences only (because the differences are fix any sample from the vcf file can be taken,
# it won't make a difference
#Select a random little b sample
sb_sample=$(bcftools query -l input/subset.vcf | grep 'littleb' | head -n 1)
bcftools consensus --fasta-ref tmp/Si_gnG.fa \
  -s ${sb_sample} tmp/subset.vcf.gz\
  > tmp/Si_gnG_littleb.fna

#When checked if correct
#mv tmp/Si_gnG_littleb.fna results/Si_gnG_littleb.fna
#gzip results/Si_gnG_littleb.fna
