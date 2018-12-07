#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_rt=4:00:00
#$ -l h_vmem=8G
#$ -t 1-6
module load star/2.6.1a
module load samtools/1.9
INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" tmp/samples.txt)
#Run all samples array. First run
#For some reason the --readFilesCommand zcat does not work with these reads, so instead I directly use <(zcat sample_name) to get
#the same result
# by deault doesn't include MD (For mismatch), include here
STAR --runThreadN 4 --genomeDir tmp/star_reference \
     --readFilesIn <(zcat input/wurm_2011_reads/${INPUT_FILE}.fastq.gz) \
     --outSAMtype BAM SortedByCoordinate --outFileNamePrefix results/${INPUT_FILE} --outSAMattributes All

#Run all samples in parallel. Second run
STAR --runThreadN 4 --genomeDir tmp/star_reference \
     --readFilesIn <(zcat input/wurm_2011_reads/${INPUT_FILE}.fastq.gz) \
     --outSAMtype BAM SortedByCoordinate --outFileNamePrefix results/${INPUT_FILE} --outSAMattributes All \
     --sjdbFileChrStartEnd results/${INPUT_FILE}SJ.out.tab

samtools index results/${INPUT_FILE}Aligned.sortedByCoord.out.bam
