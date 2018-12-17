#!/bin/bash
#Script to generate a splice junction file (SJ) using all the available RNAseq samples from both
#North America and South America populations of Solenopsis invicta. The SJ file will be used for
#alignments downstream in the analysis.

#Align all South American samples to the gnG reference of Solenopsis invicta
qsub star_alignment_all_samples_south_america.sh

#Align all North American samples to the gnG reference of Solenopsis invicta
qsub star_alignment_all_samples_north_america.sh

#The results, including the SJ files, will be in tmp
