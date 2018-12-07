#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=4:00:00
#$ -l h_vmem=2G   # Request 2GB RAM
#$ -t 1-6

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" path_to_bam.txt)
OUTPUT_FILE=$(sed -n "${SGE_TASK_ID}p" samples.ids)


java -jar GATK -R tmp/Si_gnG.fna -T ASEReadCounter -I ${INPUT_FILE} \
          -sites tmp/sorted_subset.vcf \
          -o results/${OUTPUT_FILE}.na_neutral.csv -U ALLOW_N_CIGAR_READS > results/${OUTPUT_FILE}.na_neutral.out
