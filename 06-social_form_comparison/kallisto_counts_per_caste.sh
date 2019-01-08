#!/bin/sh
#Script to generate estimated read counts for all samples from the Morandin et al., 2016 RNAseq dataset for Solenopsis
#invicta. This dataset includes 3 replicates of workers and queens from polygyne and monogyne colonies (12 samples)
# in total.

gunzip -c input/GCF_000188075.1_Si_gnG_rna_from_genomic.fna.gz > tmp/reference.fna

#Clean headers of reference to get read counts per transcript name only
cut -d ' ' -f 1 tmp/reference.fna | sed 's/lcl|[A-Z]*_[0-9]*\.1_//' | sed 's/[a-z]*_//' | sed 's/_[0-9]*$//' > tmp/reference_clean.fna
#Generate a Kallisto index:
kallisto index -i tmp/k_index_sinvicta tmp/reference_clean.fna

#Generate a list with the names of the samples:
ls input/renamed_trimmed_reads/*gz | sed 's/_[12]\.t.*//' | cut -d '/' -f 3 | sort -u > tmp/samples.ids

mkdir tmp/kallisto_counts
parallel -j 12 'kallisto quant -i tmp/k_index_sinvicta -o tmp/kallisto_counts/{} \
          --plaintext input/renamed_trimmed_reads/{}_1.t.fastq.gz input/renamed_trimmed_reads/{}_2.t.fastq.gz' \
          :::: tmp/samples.ids

#Run the R markdown script and generate a pdf with the results
 R -e "rmarkdown::render('social_forms_comparison_from_rna.Rmd',output_file='social_forms_comparison_from_rna.pdf')"
