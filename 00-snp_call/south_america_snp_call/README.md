# Variant Calling

Rodrigo Pracana
August 2018

## Introduction

This work is part of the allele specific expression analysis in the social chromosomes of the fire ant being led by Carlos Martinez Ruiz. Here, I use freebayes to perform variant calling to identify SNP sites with fixed allelic differences between SB and Sb in Argentinian samples. These samples are listed in the file `argentina_samples`.

## Inputs

The alignments used for the variant calling have been previously performed, based on reads filtered by Anurag Priyam. The first step in the pipeline is to merge the bam files for each of the individuals.

```sh

ln -sfr ../env.sh .
source env.sh

mkdir -p input
mkdir -p tmp
mkdir -p results

# Get the alignments
ln -sfr ../2017-11-28-align/results input/alignments

# Get the reference
ln -sfr ../input_data/gnG_20161213.fa input/ref.fa

mkdir -p tmp/original_alignments
while read p; do
  ln -sfr input/alignments/$p.bam tmp/original_alignments/$p.bam
  ln -sfr input/alignments/$p.bam.bai tmp/original_alignments/$p.bam.bai
done < argentina_samples

# ls tmp/original_alignments/*.bam | wc -l
# 26
# ls tmp/original_alignments/*.bam.bai | wc -l
# 26

# Make files that look like:
# -b tmp/call/parse_bam/alignments/f1_B.bam       -s f1_B
# -b tmp/call/parse_bam/alignments/f1b.bam        -s f1b
ls tmp/original_alignments/*bam \
  | awk '{print "-b " $0}' \
  > tmp/bam

ls tmp/original_alignments/*bam \
  | cut -f 3 -d "/" \
  | cut -f 1 -d "." \
  | awk '{print "-s " $0}' \
  > tmp/samples

paste tmp/bam tmp/samples > tmp/command_part

module unload gcc
module load gcc/7.1.0
../2017-12-17-call//bin/bamaddrg/bamaddrg \
  $(cat tmp/command_part | tr "\n" " " | tr "\t" " ") \
  > tmp/joined.bam
samtools index tmp/joined.bam

ln -sfr tmp/joined.bam results


```

## Run freebayes

The variant calling pipeline is in the nextflow script `call.nf`, copied directly from `../2018-05-02-call`.

```sh

source env.sh
mkdir -p logs

nextflow run call.nf \
  -w 2018-08-30-run \
  -with-report logs/2018-08-30-00-run.html \
  -resume \
  > logs/2018-08-30-00-run.log

bgzip results/raw/1_1_0.2.vcf &
bgzip results/prefilter/q25_1_1_0.2.vcf

bcftools index results/raw/1_1_0.2.vcf.gz
bcftools index results/prefilter/q25_1_1_0.2.vcf.gz

```

## Choose filtering parameters

After running freebayes, we must filter the VCF, but with which parameters? This choice is guided by work in `../2018-05-03-choose_vcf_filter_parameters`.

```sh

source env.sh

# Per pop
mkdir -p tmp/choose_filter_samples
cp argentina_samples tmp/choose_filter_samples/all.samples
grep bigB tmp/choose_filter_samples/all.samples    > tmp/choose_filter_samples/bb.samples
grep littleb tmp/choose_filter_samples/all.samples > tmp/choose_filter_samples/lb.samples

nextflow run choose_vcf_filter_parameters.nf \
  -w 2018-08-31-choose_filter-run \
  -with-report logs/2018-08-31-01-choose_filter-run.html \
  -resume \
  > logs/2018-08-31-01-choose_filter-run.log

```

## Filtering the SNPs

The report previously made suggests that the filtering strategy use for all South American individuals (`../2018-05-03-stringent_vcf_filter`) is good for this dataset, so we can go ahead and filter the SNPs.

```sh

source env.sh

nextflow run stringent_filter.nf \
  -w 2018-08-31-00-filter_vcf-run \
  -with-report logs/2018-08-31-00-filter_vcf-run.html \
  -resume \
  > logs/2018-08-31-00-filter_vcf-run.log

```

FILTER | Number of Reads | Percentage of Prefilter
--- | --- | ---
pre-filter | 2231407 | 100
HIGH_DP | 27826 | 1.2
NO_COV | 1139585 | 51.1
AMBIGUOUS_SUPPORT | 172021 | 7.7
multiple_filters | 87479 | 3.9
total_filtered_out | 1229238 | 55.1
total_remaining | 1002169 | 44.9

After filtering, we should do a quality check, as in `../2018-05-10-stringent_vcf_filter-qual`.

```sh

nextflow run post_filter_quality.nf \
  -w 2018-08-31-00-post_filter_quality-run \
  -with-report logs/2018-08-31-00-post_filter_quality-run.html \
  -resume \
  > logs/2018-08-31-00-post_filter_quality-run.log


```

## Decompose complex variants

The software freebayes can call "complex" variants. I previously removed "complex" variants including indels, but I did not remove complex variants composed of multiple SNPs next to each other. I do this here, following the scripts in `../2018-05-25-decompose_vcf_into_snpsz`. For instance, the single variant with alleles ATA and CTG becomes two variants, the first with alleles A/C, the sceond with alleles A/G.


```sh

nextflow run decompose_vcf.nf \
  -w 2018-09-01-00-decompose_vcf-run \
  -with-report logs/2018-09-01-00-decompose_vcf-run.html \
  -resume \
  > logs/2018-09-01-00-decompose_vcf-run.log


ln -sfr \
  results/filtered_vcf_decomposed/all_q25_12_1_0.6_filtered.decomposed.vcf.gz \
  results/filtered_vcf_decomposed/all.vcf.gz

ln -sfr \
  results/filtered_vcf_decomposed/all_q25_12_1_0.6_filtered.decomposed.vcf.gz.csi \
  results/filtered_vcf_decomposed/all.vcf.gz.csi

```

## Find fixed differences

We finally find the SNPs with fixed allelic differences between SB and Sb, but only in scaffolds that map to the supergene.

```sh

ln -sfr ../2018-05-11-linkage-map/results/linkage_map_supergene.txt \
  input/genetic_map.txt

mkdir -p results/fixed_differences tmp/fixed_diffs

# Supergene region only
grep supergene input/genetic_map.txt | cut -f 6 -d ' ' > tmp/fixed_diffs/scaffold_list

# VCF of supergene only
ln -sfr results/filtered_vcf_decomposed/all.vcf.gz tmp/fixed_diffs
ln -sfr results/filtered_vcf_decomposed/all.vcf.gz.csi tmp/fixed_diffs

bcftools view tmp/fixed_diffs/all.vcf.gz $(cat tmp/fixed_diffs/scaffold_list | tr '\n' ' ') \
  > tmp/fixed_diffs/supergene.vcf

# Get VCF of the fixed differences between the B and the b individuals in the
#     supergene region
Rscript fixed_diff.R

# Number of lines
grep -vc "#" tmp/fixed_diffs/supergene.vcf
52812
grep -vc "#" results/fixed_differences/2018-09-01-AR_fixed_differences.vcf
3158
# lines

cp tmp/fixed_diffs/scaffold_list results/fixed_differences

bgzip results/fixed_differences/2018-09-01-AR_fixed_differences.vcf
bcftools index results/fixed_differences/2018-09-01-AR_fixed_differences.vcf.gz

ln -sfr results/fixed_differences/2018-09-01-AR_fixed_differences.vcf.gz \
  results/fixed_differences.vcf.gz
ln -sfr results/fixed_differences/2018-09-01-AR_fixed_differences.vcf.gz.csi \
  results/fixed_differences.vcf.gz.csi

```

## Find indels only

Note that we are doing this on the non-filtered dataset.

```sh

source env.sh
mkdir -p results/indels

bcftools view --types indels results/prefilter/q25_1_1_0.2.vcf.gz \
  | bgzip -c > results/indels/unfiltered_indels.vcf.gz

tabix -p vcf results/indels/unfiltered_indels.vcf.gz

Rscript fixed_diff_indels.R

bgzip results/indels/unfiltered_fixed_different_indels.vcf
tabix -p vcf results/indels/unfiltered_fixed_different_indels.vcf.gz

```
