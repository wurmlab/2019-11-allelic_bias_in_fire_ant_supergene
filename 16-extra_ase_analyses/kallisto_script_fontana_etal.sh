#!/bin/sh
#Script to generate estimated read counts for all samples from the Fontana et al., 2019 RNAseq dataset for Solenopsis
#invicta. This dataset includes 4 replicates of queens from polygyne colonies with either SBSB or SBSb genotype (8 samples)
# in total.

gunzip -c input/GCF_000188075.1_Si_gnG_rna_from_genomic.fna.gz > tmp/reference_genes.fna

#Clean headers of reference to get read counts per transcript name only
cut -d ' ' -f 1 tmp/reference_genes.fna | sed 's/lcl|[A-Z]*_[0-9]*\.1_//' | sed 's/[a-z]*_//' | sed 's/_[0-9]*$//' > tmp/reference_genes_clean.fna
#Generate a Kallisto index:
kallisto index -i tmp/k_index_sinvicta tmp/reference_genes_clean.fna

#Run Kallisto in parallel in all samples
mkdir tmp/kallisto_counts
parallel -j 8 'kallisto quant -i tmp/k_index_sinvicta -o tmp/kallisto_counts/{} \
          --plaintext tmp/trimmed_reads_fontana_etal/trimmed_{}_1.fastq.gz tmp/trimmed_reads_fontana_etal/trimmed_{}_2.fastq.gz' \
          :::: tmp/samples.ids
