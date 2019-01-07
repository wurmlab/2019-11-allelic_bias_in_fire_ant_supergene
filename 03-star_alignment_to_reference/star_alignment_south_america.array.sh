#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_rt=12:00:00
#$ -l h_vmem=8G
#$ -t 1-28
module load star/2.6.1a
module load samtools/1.9
INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" tmp/samples.txt)

#Run all samples in parallel. Multiple mapping reads are not allowed. Only one random alignment of the multiply mapped reads
#will be considered#Run all samples in array. Second run

STAR --runThreadN 4 --genomeDir tmp/star_reference \
     --readFilesIn input/ar_trimmed_reads/${INPUT_FILE}R1_trimmed.fastq.gz input/ar_trimmed_reads/${INPUT_FILE}R2_trimmed.fastq.gz \
     --outSAMtype BAM SortedByCoordinate --outFileNamePrefix tmp/bam_files/${INPUT_FILE} --outSAMattributes All \
     --sjdbFileChrStartEnd input/all_samples_south_americaSJ.out.tab --outMultimapperOrder Random --outSAMmultNmax 1 \
     --readFilesCommand zcat

samtools index results/${INPUT_FILE}Aligned.sortedByCoord.out.bam
