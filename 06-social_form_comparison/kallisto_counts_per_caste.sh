#!/bin/sh
#Script to generate estimated read counts for all samples from the Morandin et al., 2016 RNAseq dataset for Solenopsis
#invicta. This dataset includes 3 replicates of workers and queens from polygyne and monogyne colonies (12 samples)
# in total.

gunzip -c input/reference.fna.gz > tmp/reference.fna
#Generate a Kallisto index:
kallisto index -i tmp/k_index_ase_sinvicta tmp/reference.fna

#Generate a list with the names of the samples:
ls input/renamed_trimmed_reads/*gz | sed 's/_[12]\.t.*//' | cut -d '/' -f 3 | sort -u > tmp/samples.ids

mkdir tmp/kallisto_counts
parallel -j 12 'kallisto quant -i tmp/k_index_ase_sinvicta -o tmp/kallisto_counts/{} \
          --plaintext input/renamed_trimmed_reads/{}_1.t.fastq.gz input/renamed_trimmed_reads/{}_2.t.fastq.gz' \
          :::: tmp/samples.ids
