# Variant Calling

Rodrigo Pracana
August 2018

## Introduction

This work is part of the allele specific expression analysis in the social chromosomes of the fire ant being led by Carlos Martinez Ruiz. Here, I use freebayes to perform variant calling to identify SNP sites with fixed allelic differences between SB and Sb in North American samples. These samples are listed in the file `north_american_samples`.

## Inputs

We have done two runs of sequencing, 1st_run and 2nd_run. Sequencing is paired-end, with an R1 and an R2 file for each sample. The R2 files in the first run are very short (which is why we did a second run), and therefore we use only the first reads of that sample, as though they were the result of single-end sequencing. For the first run of sequencing, we use the filtering done for the OBP project (Pracana et al 2017, Evolution Letters), which can be found [here](https://github.com/roddypr/obps/blob/master/bin/roddy/variants/read_filter.sh). For the second run, we use the reads filtered by Anurag Priyam using his pipeline [reads-qc](https://github.com/wurmlab/reads-qc).

In the following, I use links to transform the original file names so they have the file names listed in `north_american_samples`. I then copy the read files into `tmp`.

```sh

ln -sfr ../env.sh .
source env.sh

mkdir -p input
mkdir -p tmp
mkdir -p results


# Pair ended reads
#  But make their name similar to the sample index, rather than the original names
cat north_american_samples | cut -f 1 -d "-" | sort | uniq > tmp/sample_starts

ARCHIVE=/data/SBCS-WurmLab/archive/db/genomic/reads/S_invicta
ARCHIVE=${ARCHIVE}/2018-01-15-NorthAmerica_7B7b_males_HiSeq2000_filtered/merged_by_sample

mkdir -p input/pe_reads
mkdir -p tmp/pe_reads

while read p; do
  # BigB samples
  ln -sf ${ARCHIVE}/${p}B.R1.fastq.gz input/pe_reads/${p}-bigB-R1.fq.gz
  ln -sf ${ARCHIVE}/${p}B.R2.fastq.gz input/pe_reads/${p}-bigB-R2.fq.gz

  ln -sf ${ARCHIVE}/${p}b.R1.fastq.gz input/pe_reads/${p}-littleb-R1.fq.gz
  ln -sf ${ARCHIVE}/${p}b.R2.fastq.gz input/pe_reads/${p}-littleb-R2.fq.gz

done < tmp/sample_starts

mkdir -p tmp
cp -r input/pe_reads tmp

# Single end reads
ARCHIVE=/data/archive/archive-SBCS-WurmLab/rpracana/projects/obp/current/obp
ARCHIVE=${ARCHIVE}/analysis/01-gnh/analysis/00-read_filtering/results

mkdir -p input/se_reads
mkdir -p tmp/se_reads

while read p; do

  ln -sf ${ARCHIVE}/${p}_B.se.fq.gz input/se_reads/${p}-bigB.fq.gz
  ln -sf ${ARCHIVE}/${p}b.se.fq.gz input/se_reads/${p}-littleb.fq.gz

done < tmp/sample_starts

```

We also need the reference genome assembly. Carlos asked me to use Si_gnG.

```sh

mkdir -p tmp

ln -sfr \
  /data/SBCS-WurmLab/archive/db/genomic/S_invicta/Si_gnG/NCBI/GCF_000188075.1_Si_gnG_genomic.fna.gz \
  input/ref.fa.gz

gunzip -c input/ref.fa.gz > tmp/ref.fa

```


## Alignments

We align the reads of the reference genome assembly using bowtie. I performed the following in a nextflow script.

```sh

mkdir -p logs

nextflow run align.nf \
  -w 2018-09-05-00-align-run \
  -resume \
  -with-report logs/2018-09-05-01-align-run.html \
  > logs/2018-09-05-01-align-run.log

```

```sh
# ls tmp/original_alignments/*.bam | wc -l
# 26
# ls tmp/original_alignments/*.bam.bai | wc -l
# 26

mkdir -p tmp/joinbam
# Make files that look like:
# -b tmp/call/parse_bam/alignments/f1_B.bam       -s f1_B
# -b tmp/call/parse_bam/alignments/f1b.bam        -s f1b
ls results/alignments/*bam \
  | awk '{print "-b " $0}' \
  > tmp/joinbam/bam

ls results/alignments/*bam \
  | cut -f 3 -d "/" \
  | cut -f 1 -d "." \
  | awk '{print "-s " $0}' \
  > tmp/joinbam/samples

paste tmp/joinbam/bam tmp/joinbam/samples > tmp/joinbam/command_part

module unload gcc
module load gcc/7.1.0
../2017-12-17-call//bin/bamaddrg/bamaddrg \
  $(cat tmp/joinbam/command_part | tr "\n" " " | tr "\t" " ") \
  > tmp/joinbam/joined.bam
samtools index tmp/joinbam/joined.bam

ln -sfr tmp/joinbam/joined.bam results


```

## Run freebayes

The variant calling pipeline is in the nextflow script `call.nf`, copied directly from `../2018-05-02-call`.

```sh

source env.sh
mkdir -p logs

nextflow run call.nf \
  -w 2018-09-08-call-run \
  -with-report logs/2018-09-08-00-call-run.html \
  -resume \
  > logs/2018-09-08-00-call-run.log

bgzip results/raw/1_1_0.2.vcf
bgzip results/prefilter/q25_1_1_0.2.vcf

```

## Choose filtering parameters

After running freebayes, we must filter the VCF, but with which parameters? This choice is guided by work in `../2018-05-03-choose_vcf_filter_parameters`.

```sh

source env.sh

# Per pop
mkdir -p tmp/choose_filter_samples
cp north_american_samples tmp/choose_filter_samples/all.samples
grep bigB tmp/choose_filter_samples/all.samples    > tmp/choose_filter_samples/bb.samples
grep littleb tmp/choose_filter_samples/all.samples > tmp/choose_filter_samples/lb.samples

nextflow run choose_vcf_filter_parameters.nf \
  -w 2018-09-08-choose_filter-run \
  -with-report logs/2018-09-08-00-choose_filter-run.html \
  -resume \
  > logs/2018-09-08-00-choose_filter-run.log

## Repeat with changes to report
nextflow run choose_vcf_filter_parameters.nf \
  -w 2019-11-07-choose_filter-run \
  -with-report logs/2019-11-07-00-choose_filter-run.html \
  -resume \
  > logs/2018-11-07-00-choose_filter-run.log &

```

## Filtering the SNPs

The report previously made suggests that the filtering strategy use for all South American individuals (`../2018-05-03-stringent_vcf_filter`) is good for this dataset, so we can go ahead and filter the SNPs.

```sh

source env.sh

nextflow run stringent_filter.nf \
  -w 2018-09-08-00-filter_vcf-run \
  -with-report logs/2018-09-08-00-filter_vcf-run.html \
  -resume \
  > logs/2018-09-08-00-filter_vcf-run.log

```

FILTER | Number of Reads | Percentage of Prefilter
--- | --- | ---
pre-filter | 1011057 | 100
HIGH_DP | 22514 | 2.2
NO_COV | 51402 | 5.1
AMBIGUOUS_SUPPORT | 40409 | 4
multiple_filters | 6565 | 0.6
total_filtered_out | 90149 | 8.9
total_remaining | 920908 | 91.1

After filtering, we should do a quality check, as in `../2018-05-10-stringent_vcf_filter-qual`.

```sh

nextflow run post_filter_quality.nf \
  -w 2018-09-08-00-post_filter_quality-run \
  -with-report logs/2018-09-08-00-post_filter_quality-run.html \
  -resume \
  > logs/2018-09-08-00-post_filter_quality-run.log


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
  results/filtered_vcf_decomposed/all_q25_16_1_0.6_filtered.decomposed.vcf.gz \
  results/filtered_vcf_decomposed/all.vcf.gz

ln -sfr \
  results/filtered_vcf_decomposed/all_q25_16_1_0.6_filtered.decomposed.vcf.gz.csi \
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
