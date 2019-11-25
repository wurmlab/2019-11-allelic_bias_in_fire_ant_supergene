#!/bin/bash
#Script to perform the complete ASE analysis using a different dataset, in this case, the data from SBSB and SBSb virgin
#polygyne queens from Fontana et al 2019 Mol Ecol.

#FIRST STEP: READ CLEANING, QC AND ALIGNMENT--------------------------------------------------------------------------------

#Remove adapters from the raw RNAseq reads of the Solenopsis invicta queen samples from Fontana et al 2019
#Generate a list of sample names:
ls input/reads_fontana_etal | cut -d "_" -f 1 | sort -u > tmp/samples.ids

mkdir tmp/trimmed_reads_fontana_etal

#Run cutadapt
qsub cutadapt_script.array.sh
