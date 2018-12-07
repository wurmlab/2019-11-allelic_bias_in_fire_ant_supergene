#!/bin/bash
#Script to generate BAM files based on the refseq gnG assembly of the reference genome of Solenopsis invicta.
#This script generates an alignment for the SOUTH AMERICAN POPULATIONS of Solenopsis invcita
module load star/2.6.1a
#Generate a text file with the sample names
ls input/ar_trimmed_reads/*gz | cut -d "/"  -f3 | cut -d "R" -f1 | sort -u > tmp/samples.txt
mkdir tmp/star_reference
gunzip -c input/reference.fna.gz > tmp/reference.fna

#Generate index using refseq's gnG assembly and annotations
STAR --runThreadN 15 --runMode genomeGenerate --genomeDir tmp/star_reference \
     --genomeFastaFiles tmp/reference.fna \
     --sjdbGTFtagExonParentTranscript Parent \
     --sjdbGTFfile input/GCF_000188075.1_Si_gnG_genomic.gff --sjdbOverhang 149 \
     --limitGenomeGenerateRAM 60000000000

#Run actual alignment in an array job
qsub star_alignment_south_america.array.sh
