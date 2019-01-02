#!/bin/bash
#Script to generate BAM files based on the refseq gnG assembly of the reference genome of Solenopsis invicta.
#This script generates an alignment for the SOUTH AMERICAN POPULATIONS of Solenopsis invcita
#Generate a text file with the sample names
ls input/ar_trimmed_reads/*gz | cut -d "/"  -f3 | cut -d "R" -f1 | sort -u > tmp/samples.txt
gunzip -c input/reference.fna.gz > tmp/reference.fna

mkdir tmp/bam_files
#Run actual alignment in an array job
qsub star_alignment_south_america.array.sh
