#Script to generate a vcf with fixed differences between SB and Sb for the supergene region only.
module load R/3.4.3
module load bcftools/1.8
module load vcftools/0.1.15
#Generate the file 'supergene_scaffolds.ids' with only the names of the scaffolds within the supergene region
#in the gnG assembly of Solenopsis invicta.
Rscript supergene_scaffolds.R

mkdir tmp/fixed_diffs
#Index the original vcf
ln -sfr input/all.vcf.gz tmp/fixed_diffs/
tabix -f tmp/fixed_diffs/all.vcf.gz

#Filter the original vcf to include only scaffolds in the supergene
bcftools view tmp/fixed_diffs/all.vcf.gz $(cat tmp/supergene_scaffolds.ids | tr '\n' ' ') > tmp/fixed_diffs/supergene.vcf

#Use the R package VariantAnnotation to filter only for fixed differences
Rscript fixed_diffs.R

#When checked if results are correct:
#mv tmp/supergene_fixed_diffs_north_america.vcf results/supergene_fixed_diffs_north_america.vcf
#mv tmp/supergene_fixed_diffs_south_america.vcf results/supergene_fixed_diffs_south_america.vcf

#Number of fixed differences found:
#North America
#grep -vc '#' results/supergene_fixed_diffs_north_america.vcf
#25899

#South America
#grep -vc '#' results/supergene_fixed_diffs_south_america.vcf
#3129

#The vast majority of the reference variant belong to SB and the alternative, to Sb, but in a few cases this is not true.
#Find the SNPs that are in the 'wrong' category (alternative in all SBs or reference in all Sbs)
#In North America
vcftools --vcf results/supergene_fixed_diffs_north_america.vcf --indv f1-bigB --extract-FORMAT-info GT --out tmp/genotypes
grep '1$' tmp/genotypes.GT.FORMAT > tmp/misassigned_snps_na.tab
#How many are misassigned?
wc -l tmp/misassigned_snps_na.tab
#125

#In South America
vcftools --vcf results/supergene_fixed_diffs_south_america.vcf --indv AR102-1-bigB-p --extract-FORMAT-info GT --out tmp/genotypes
grep '1$' tmp/genotypes.GT.FORMAT > tmp/misassigned_snps_sa.tab
#How many are misassigned?
wc -l tmp/misassigned_snps_sa.tab
#33
