#!/bin/sh
module load bcftools/1.8
module load R/3.4.3
module load parallel
module load/python/3.6.3
#Script to extract the indels with fixed differences between SB and Sb in the supergene region in common
#between North and South American populations of S.invicta

#Subset the vcfs to include only supergene regions
#Generate a list with the suepergene scaffolds
Rscript supergene_scaffolds.R

#Filter both files to include only scaffolds in the supergene
bcftools view input/unfiltered_fixed_different_indels_na.vcf.gz $(cat tmp/supergene_scaffolds.ids | tr '\n' ' ') > tmp/supergene_fixed_na.vcf
bcftools view input/unfiltered_fixed_different_indels_sa.vcf.gz $(cat tmp/supergene_scaffolds.ids | tr '\n' ' ') > tmp/supergene_fixed_sa.vcf

#Bgzip and index the new vcfs
ls tmp/*vcf | parallel 'bgzip {}' :::
ls tmp/*vcf.gz | parallel 'bcftools index' :::

#Merge the two vcfs
bcftools isec -n=2 -w1 tmp/supergene_fixed_na.vcf.gz tmp/supergene_fixed_sa.vcf.gz -o tmp/common_indels.vcf

#How many indels?
#176

#This takes into account unfiltered indels. It shouldn't matter because indels with fixed differences between SB and Sb
#(i.e. present in ALL SB individuals and ALL Sb individuals across populations) should be of good quality, but it is worth double-checking
#Filter out variants with quality under 20 and coverage under 60
bcftools query -i'QUAL>20 && DP>60' -f'%CHROM %POS %QUAL %DP\n' tmp/common_indels.vcf | wc -l
#176

## RUN SNPEff
#Generate an Ensembl vcf for all common indels using the chain file generated in ../2019-08-15-snps_effects/
gunzip -c input/Si_gnG_ensembl.fa.gz > tmp/Si_gnG_ensembl.fa
source ~/crossmaps/bin/activate
CrossMap.py vcf input/liftover.chn tmp/common_indels.vcf tmp/Si_gnG_ensembl.fa tmp/ensembl_common_indels.vcf
deactivate

#Run SNPEff
java -jar snpEff.jar Solenopsis_invicta -v\
          tmp/ensembl_common_indels.vcf > tmp/eff_common_indels.vcf
