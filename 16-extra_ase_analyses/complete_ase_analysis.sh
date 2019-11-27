#!/bin/bash
#Script to perform the complete ASE analysis using a different dataset, in this case, the data from SBSB and SBSb virgin
#polygyne queens from Fontana et al 2019 Mol Ecol.
module load fastqc
#FIRST STEP: READ CLEANING, QC AND ALIGNMENT--------------------------------------------------------------------------------

#Remove adapters from the raw RNAseq reads of the Solenopsis invicta queen samples from Fontana et al 2019
#Generate a list of sample names:
ls input/reads_fontana_etal | cut -d "_" -f 1 | sort -u > tmp/samples.ids

mkdir tmp/trimmed_reads_fontana_etal

#Run cutadapt
qsub cutadapt_script.array.sh

#QC the trimmed reads
mkdir tmp/fastqc_trimmed
ls tmp/trimmed_reads_fontana_etal/* | parallel 'fastqc -o tmp/fastqc_trimmed' :::
source ~/multiqc/bin/activate
multiqc tmp/fastqc_trimmed
deactivate
#All adapter removed
#Once the reads have been cleaned, align to the normal reference and to the Sb-ified reference.
#The samples used here were collected in Taiwan. We do not have SNPs between variants for this particular population,
#but the Taiwan populations of S.invicta are derived from the invasive North American ones (Ascunce et al 2011, Science).
#So here we will be using the Sb-ified reference generated using the information from North American populations:

#First, generate a splice junction file for each reference
qsub star_alignment_all_samples_bigb_fontana_etal.sh
qsub star_alignment_all_samples_littleb_fontana_etal.sh

#Once these files have been generated, run the alignment per sample
qsub star_alignment_fontana_etal_bigb.array.sh
qsub star_alignment_fontana_etal_littleb.array.sh
