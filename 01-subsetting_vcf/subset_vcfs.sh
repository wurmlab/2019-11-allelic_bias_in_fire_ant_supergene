#Script to generate a vcf with fixed differences between SB and Sb for the supergene region only.
module load R/3.4.3
module load bcftools/1.8
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

#Number of fixed differences found:
#North America
#grep -vc '#' results/supergene_fixed_diffs_north_america.vcf
#25899

#South America
#grep -vc '#' results/supergene_fixed_diffs_south_america.vcf
#3129
