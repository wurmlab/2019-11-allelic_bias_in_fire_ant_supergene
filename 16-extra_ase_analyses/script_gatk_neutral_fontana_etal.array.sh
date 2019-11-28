#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=4:00:00
#$ -l h_vmem=2G   # Request 2GB RAM
#$ -t 1-8

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" tmp/path_to_bam.ids)
OUTPUT_FILE=$(sed -n "${SGE_TASK_ID}p" tmp/samples.ids)


java -jar GATK -R tmp/reference.fna -T ASEReadCounter -I ${INPUT_FILE} \
          -sites tmp/sorted_subset.vcf \
          -o tmp/${OUTPUT_FILE}.sa_neutral.csv -U ALLOW_N_CIGAR_READS > tmp/${OUTPUT_FILE}.sa_neutral.out
