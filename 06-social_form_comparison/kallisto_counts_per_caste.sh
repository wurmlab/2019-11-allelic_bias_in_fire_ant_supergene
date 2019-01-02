#!/bin/sh
module load kallisto/0.42.5
#Generate a Kallisto index:

kallisto index -i input/k_index_ase_sinvicta input/GCF_000188075.1_Si_gnG_rna_from_genomic.fna

#Generate a list with the names of the samples:
ls input/trimmed_reads/*gz|cut -d 'R' -f 1 | cut -d '/' -f 3 | uniq > samples.txt

parallel -j 28 'kallisto quant -i input/k_index_ase_sinvicta -o ./{} \
          --plaintext input/trimmed_reads/{}R1_trimmed.fastq.gz input/trimmed_reads/{}R2_trimmed.fastq.gz' \
          :::: samples.txt
