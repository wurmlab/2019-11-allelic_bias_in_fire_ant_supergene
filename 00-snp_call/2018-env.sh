#!/bin/sh

module load bcftools/1.4 htslib/1.6 samtools/1.6
module load vcftools/0.1.15

# Variant calling software
module load freebayes/1.1.0

# nextflow
module load nextflow/0.26.3

# R
module load R/3.4.3

# VCFLIB
PATH=$PATH:/data/home/btw749/software/vcflib/bin

# ruby
module load ruby/2.4.1

# bedtools
module load bedtools/2.26.0
