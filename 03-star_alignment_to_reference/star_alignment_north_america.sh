#!/bin/bash
#Script to generate BAM files aligning the RNAseq reads from Wurm et al 2011 to the gnG reference
#assembly of Solenopsis invcta.
module load star/2.6.1a
#Generate a text file with the sample names. Sample SRR619956 and SRR619959 seem to be the same
#by looking at the SRA metadata on NCBI. Remove one of them
ls input/wurm_2011_reads/*gz | cut -d "/" -f3 | cut -d '.' -f 1 | \
   sed '/SRR619959/d'  > tmp/samples.txt
gunzip -c input/reference.fna.gz > tmp/reference.fna

mkdir tmp/bam_files
#Run actual alignment in an array job
qsub star_alignment_north_america.array.sh
